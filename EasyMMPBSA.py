#!/usr/bin/env python
import os
import pygtk
import gtk
import gtk.glade
import glob
from AMBER import *
from MMPBSA import *

class App:
    def __init__(self):
        self.wTree = gtk.glade.XML(os.path.join(os.path.dirname(__file__), "MMPBSA.glade"))

        # Connect signals to handlers
        self.window = self.wTree.get_widget("main_window")
        self.wTree.signal_autoconnect(self)
        self.window.connect("destroy", gtk.main_quit)
        self.wTree.get_widget("quit1").connect("activate", gtk.main_quit)
        self.wTree.get_widget("about_dialog_ok").connect("clicked", self.hide_about_box)

        # Set defaults
        self.wTree.get_widget("decomp_type").set_active(1)
        self.wTree.get_widget("solvent_approx").set_active(0)
        self.wTree.get_widget("sasa_method").set_active(1)
        self.on_solvent_approx_changed()
        
        # Set window title to include current directory, with home directory
        # part replaced by a tilde.
        self.window.set_title("%s: %s" % (self.window.get_title(), \
            os.getcwd().replace(os.path.expanduser('~'), '~')))
        
        self.prmtop_names = {}
    
    def alert(self, msg):
        """Shows an alert dialog box."""
        dialog = gtk.Dialog(title="MM/PBSA Calculation Setup", parent=self.window, \
            flags = gtk.DIALOG_MODAL, buttons=("OK", 0))
        dialog.vbox.pack_start(gtk.Label(msg), True, True, 0)
        dialog.show_all()
        dialog.run()
        dialog.destroy()
        
    def prepare_residue_energy_calcuation(self):
        # Assemble list of frames by looking in the current directory
        # for extracted snapshots
        frames = {}
        ends_in_number = re.compile('[0-9]+$')
        for type in ["ligand", "receptor", "complex"]:
            frames[type] = glob.glob("snapshots/%s.crd.*" % type)
            frames[type] = [x for x in frames[type] if ends_in_number.search(x)]
            if len(frames[type]) == 0:
                del frames[type]
        
        # Determine the location of ligand and receptor within complex
        # The user should have specified this on the other tab
        params = {}
        prmtops = {}
        try:
            for type in frames:
                # Do we have any frames of that type?
                if len(frames[type]) > 0:
                    if type not in self.prmtop_names:
                        continue
                    prmtops[type] = AmberSystem(self.prmtop_names[type])
                    for blah in ["start", "end"]:
                        key = "%s_atom_%s" % (type, blah)
                        params[key] = int(self.wTree.get_widget(key).get_text())
        except ValueError:
            print "Did you specify all your prmtop files?!"
            return
            
        (complex_atoms, receptor_atoms, ligand_atoms) = self.get_wanted_atoms()
        
        try:
            max_cpus = int(self.wTree.get_widget("max_cpu").get_text())
        except ValueError:
            # User might have put garbage in max CPUs box
            print "How many CPUs do you want?"
            return

        runner = JobRunner(max_cpus)

        jobs_by_type = {}

        # Generate and run sander jobs
        for type in prmtops:
            assert(type in frames)
            jobs_by_type[type] = []
            for frame in frames[type]:
                num_residues = prmtops[type].num_residues()
                footer = '''Residues for decomposition
LRES %d %d
END
Residues to print
RES %d %d
END
END''' % (1, num_residues, 1, num_residues)        
                job = SanderJob.molecular_mechanics(unique_id=os.path.basename(frame), \
                    prmtop_file=prmtops[type].name, \
                    crd_file=frame, \
                    footer=footer)

                # Get parameters from the GUI
                solvent_approx = self.wTree.get_widget("solvent_approx").get_active()
                assert solvent_approx >= 0 and solvent_approx <= 2
                if solvent_approx == 0:
                    igb = 2
                elif solvent_approx == 1:
                    igb = 1
                else:
                    print "PB is not supported yet!"
                    return
                job.set_param('igb', igb)

                # SASA method: gbsa
                sasa_method = self.wTree.get_widget("sasa_method").get_active()
                assert sasa_method == 0 or sasa_method == 1
                job.set_param('gbsa', sasa_method + 1)

                # GB parameters
                try:
                    saltcon = float(self.wTree.get_widget("salt_concentration").get_text())
                    extdiel = float(self.wTree.get_widget("ext_dielectric").get_text())
                    intdiel = float(self.wTree.get_widget("int_dielectric").get_text())
                    surften = float(self.wTree.get_widget("surface_tension").get_text())
                    surfoff = float(self.wTree.get_widget("surfoff").get_text())
                except ValueError:
                    print "One of the GB parameters is messed up."
                    return

                job.set_param('saltcon', saltcon)
                job.set_param('extdiel', extdiel)
                job.set_param('intdiel', intdiel)
                job.set_param('surften', surften)
                # TODO: what do we do with surfoff?
                # idecomp=4 specifies pairwise residue decomposition
                job.set_param('idecomp', 4)
                runner.add(job)
                jobs_by_type[type].append(job)

        return (runner, prmtops, jobs_by_type)
        
    def calculate_residue_energies(self, widget=None):
        info("Calculating differences in residue interaction energies.")
        
        (runner, prmtops, jobs_by_prmtop) = self.prepare_residue_energy_calcuation()
        runner.go()
        
    def postprocess_residue_energies(self, widget=None):

        (runner, prmtops, jobs_by_prmtop) = self.prepare_residue_energy_calcuation()

        # Handle all energy terms
        for term in ['int', 'vdw', 'eel', 'pol', 'sa', 'total']:
            energies = {}
            for prmtop in prmtops:
                energies[prmtop] = PairwiseEnergies(prmtops[prmtop].num_residues())
                for job in jobs_by_prmtop[prmtop]:
                    info("Reading energies from %s..." % job.mdout_filename)
                    energies[prmtop].add_from_mdout_file(job.mdout_filename, \
                        job.frame_id, terms=[term])
                
                # Dump them all to a file.
                fp = open("%s.%s.out" % (prmtop, term), "w")
                for res_i in inclusive_range(1, prmtops[prmtop].num_residues()):
                    row = [str(energies[prmtop].mean_over_all_frames('TDC', term, res_i, res_j)) \
                        for res_j in inclusive_range(1, prmtops[prmtop].num_residues())]
                    print >>fp, " ".join(row)
                fp.close()
        
        info("Hooray!")
     
    # Convenience method to enable or disable text entry widgets
    def set_text_and_activation(self, widget, s):
        if s is None:
            s = ""
            widget.set_property("sensitive", False)
        else:
            s = str(s)
            widget.set_property("sensitive", True)

        widget.set_text(s)
    
    # Guesses the atom and residue ranges that map receptor and ligand into the complex    
    def guess_atom_ranges(self):
        cs, ce, rs, re, ls, le = None, None, None, None, None, None
        rrs, rre, rls, rle, tcs, tce = None, None, None, None, None, None
        complex_atom_names = None
        
        try:
            if "complex" in self.prmtop_names:
                # If the file is not valid, give up
                try:
                    complex_atom_names = get_atom_names_from_prmtop(self.prmtop_names["complex"])
                    complex_residue_names = get_residue_names_from_prmtop(self.prmtop_names["complex"])
                except:
                    return
                cs = 1
                ce = len(complex_atom_names) / 4
            else:
                # No point in doing the rest of it if there is no complex prmtop
                return
            
            # Aaah! Duplicated code! TODO: refactor
            
            if "fullsystem" in self.prmtop_names:
                fullsystem_atom_names = get_atom_names_from_prmtop(self.prmtop_names["fullsystem"])
                start = fullsystem_atom_names.find(complex_atom_names)
                tcs, tce = "", ""
                if start >= 0:
                    tcs = (start / 4) + 1
                    tce = tcs + len(complex_atom_names) / 4 - 1

            if "receptor" in self.prmtop_names:
                receptor_atom_names = get_atom_names_from_prmtop(self.prmtop_names["receptor"])
                receptor_residue_names = get_residue_names_from_prmtop(self.prmtop_names["receptor"])
                start = complex_atom_names.find(receptor_atom_names)
                rstart = complex_residue_names.find(receptor_residue_names)
                rs, re, rrs, rre = "", "", "", ""
                if start >= 0:
                    rs = (start / 4) + 1
                    re = rs + len(receptor_atom_names) / 4 - 1
                if rstart >= 0:
                    rrs = (rstart / 4) + 1
                    rre = rrs + len(receptor_residue_names) / 4 - 1

            if "ligand" in self.prmtop_names:
                ligand_atom_names = get_atom_names_from_prmtop(self.prmtop_names["ligand"])
                ligand_residue_names = get_residue_names_from_prmtop(self.prmtop_names["ligand"])
                start = complex_atom_names.find(ligand_atom_names)
                rstart = complex_residue_names.find(ligand_residue_names)
                ls, le, rls, rle = "", "", "", ""
                if start >= 0:
                    ls = (start / 4) + 1
                    le = ls + len(ligand_atom_names) / 4 - 1
                if rstart >= 0:
                    rls = (rstart / 4) + 1
                    rle = rls + len(ligand_residue_names) / 4 - 1
            
        finally:
            self.set_text_and_activation(self.wTree.get_widget("complex_atom_start"), cs)
            self.set_text_and_activation(self.wTree.get_widget("complex_atom_end"), ce)
            self.set_text_and_activation(self.wTree.get_widget("receptor_atom_start"), rs)
            self.set_text_and_activation(self.wTree.get_widget("receptor_atom_end"), re)
            self.set_text_and_activation(self.wTree.get_widget("ligand_atom_start"), ls)
            self.set_text_and_activation(self.wTree.get_widget("ligand_atom_end"), le)
            self.set_text_and_activation(self.wTree.get_widget("receptor_residue_start"), rrs)
            self.set_text_and_activation(self.wTree.get_widget("receptor_residue_end"), rre)
            self.set_text_and_activation(self.wTree.get_widget("ligand_residue_start"), rls)
            self.set_text_and_activation(self.wTree.get_widget("ligand_residue_end"), rle)
            self.set_text_and_activation(self.wTree.get_widget("traj_complex_atom_start"), tcs)
            self.set_text_and_activation(self.wTree.get_widget("traj_complex_atom_end"), tce)
                
    # Called when a prmtop file is chosen.
    # When a new prmtop file is chosen, it gives us the chance to re-guess the
    # atom and residue ranges of the receptor and ligand in the complex.    
    def on_choose_prmtop(self, widget, data=None):
        # Save the selected filenames into instance variables
        self.fullsystem_trajectory = self.wTree.get_widget("trajectory_chooser").get_filename()
        
        self.prmtop_names = {}
        for type in ["fullsystem", "complex", "receptor", "ligand"]:
            filename = self.wTree.get_widget("%s_prmtop_chooser" % type).get_filename()
            if filename is not None:
                self.prmtop_names[type] = filename
        
        self.guess_atom_ranges()
    
    # Hide and show the appropriate preferences for GB or PB as they are selected    
    def on_solvent_approx_changed(self, widget=None, data=None):
        selected_index = self.wTree.get_widget("solvent_approx").get_active()
        gb_params = self.wTree.get_widget("gb_params")
        pb_params = self.wTree.get_widget("pb_params")
        if selected_index == 0 or selected_index == 1: # GB
            pb_params.set_property("visible", False)
            gb_params.set_property("visible", True)
        elif selected_index == 2: # PB
            gb_params.set_property("visible", False)
            pb_params.set_property("visible", True)
        else:
            gb_params.set_property("visible", False)
            pb_params.set_property("visible", False)
            
    def get_wanted_atoms(self):
        """Extracts the atom ranges from the GUI."""
        try:
            complex_atom_start = int(self.wTree.get_widget("complex_atom_start").get_text())
            complex_atom_end = int(self.wTree.get_widget("complex_atom_end").get_text())
            complex_wanted_atoms = inclusive_range(complex_atom_start, complex_atom_end)
        except ValueError:
            # self.alert("Complex parameters messed up?")
            complex_wanted_atoms = None

        try:
            receptor_atom_start = int(self.wTree.get_widget("receptor_atom_start").get_text())
            receptor_atom_end = int(self.wTree.get_widget("receptor_atom_end").get_text())
            receptor_wanted_atoms = inclusive_range(receptor_atom_start, receptor_atom_end)
        except ValueError:
            # self.alert("Receptor parameters messed up?")
            receptor_wanted_atoms = None
            
        try:
            ligand_atom_start = int(self.wTree.get_widget("ligand_atom_start").get_text())
            ligand_atom_end = int(self.wTree.get_widget("ligand_atom_end").get_text())
            ligand_wanted_atoms = inclusive_range(ligand_atom_start, ligand_atom_end)
        except ValueError:
            # self.alert("Ligand parameters messed up?")
            ligand_wanted_atoms = None
            
        return (complex_wanted_atoms, receptor_wanted_atoms, ligand_wanted_atoms)
    
    # When the "Extract snapshots" button is clicked, actually extract the snapshots        
    def on_extract_snapshots_clicked(self, widget, data=None):
        # TODO: should be able to figure out whether it has a box from the prmtop
        traj = AmberTrajectory.from_crd_file(
            filename = self.fullsystem_trajectory,
            atoms_per_frame = guess_num_atoms_from_prmtop(self.prmtop_names["fullsystem"]),
            has_box = True)
        
        try:
            frame_stride = int(self.wTree.get_widget("frame_stride").get_text())
            wanted_frames = inclusive_range(500, 8000, frame_stride)
            info("Ignoring start and end frames for now...")
        except ValueError:
            self.alert("Which frames do you want?")
            return
            
        print "Want frames:", wanted_frames

        (complex_wanted_atoms, receptor_wanted_atoms, ligand_wanted_atoms) \
            = self.get_wanted_atoms()
        
        current_dir = os.getcwd()
        os.mkdir("snapshots", 0744)
        os.chdir("snapshots")
        if complex_wanted_atoms is not None:
            if frame_stride > 1:
                traj.get_frames(wanted_frames, complex_wanted_atoms).write_crd_files("complex")
            else:
                traj.snapshotify("complex")
        if receptor_wanted_atoms is not None:
            traj.get_frames(wanted_frames, receptor_wanted_atoms).write_crd_files("receptor")
        if ligand_wanted_atoms is not None:
            traj.get_frames(wanted_frames, ligand_wanted_atoms).write_crd_files("ligand")
        os.chdir(current_dir)
            
        print "Done writing snapshots."
        
    def show_about_box(self, widget=None, data=None):
        self.about_dialog.set_property("visible", True)
        
    def hide_about_box(self, widget=None, data=None):
        self.about_dialog.set_property("visible", False)
        
    def main(self):
        gtk.main()
        
if __name__ == "__main__":
    app = App()
    app.main()
