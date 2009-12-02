// Calculates the nonbonded energy (vdW and electrostatic) in an AMBER system.
// Assumes no cutoff. Does not calculate any other terms.
package main
import ( "encoding/binary"; "math"; "fmt"; "flag"; "os"; "bufio";
    "amber";
)

func main() {
    var prmtopFilename, rstFilename, trjFilename, outFilename string;
    var numFrames int;
    flag.StringVar(&prmtopFilename, "p", "prmtop", "Prmtop filename");
    flag.StringVar(&rstFilename, "c", "", "Inpcrd/rst filename");
    flag.StringVar(&trjFilename, "x", "", "Trajectory (in text format) filename");
    flag.IntVar(&numFrames, "n", 0, "Number of frames in trajectory to process");
    flag.StringVar(&outFilename, "o", "energies.bin", "Energy decomposition output filename");
    flag.Parse();

    mol := amber.LoadSystem(prmtopFilename);
    if mol == nil { return }
    fmt.Println("Number of atoms:", mol.NumAtoms());
    fmt.Println("Number of residues:", mol.NumResidues());
    fmt.Print("Periodic box in prmtop: ");
    if mol.GetInt("POINTERS", amber.IFBOX) > 0 {
        fmt.Println("Yes");
    } else {
        fmt.Println("No")
    }

    // Set up nonbonded parameters. We load them here so we don't have to keep
    // doing it later
    var params NonbondedParamsCache;
    params.Ntypes = mol.GetInt("POINTERS", amber.NTYPES);
    params.NBIndices = amber.VectorAsIntArray(mol.Blocks["NONBONDED_PARM_INDEX"]); // ICO
    params.AtomTypeIndices = amber.VectorAsIntArray(mol.Blocks["ATOM_TYPE_INDEX"]); // IAC
    params.LJ12 = amber.VectorAsFloat32Array(mol.Blocks["LENNARD_JONES_ACOEF"]); // CN1
    params.LJ6 = amber.VectorAsFloat32Array(mol.Blocks["LENNARD_JONES_BCOEF"]); // CN2
    params.Charges = amber.VectorAsFloat32Array(mol.Blocks["CHARGE"]);
    // If we were given a single snapshot, just do that one
    if rstFilename != "" {
        fmt.Printf("Calculating energies for single snapshot %s.\n", rstFilename);
        mol.LoadRst(rstFilename);
    
        var request EnergyCalcRequest;
        request.NBParams = params;
        request.Molecule = mol;
        request.Coords = mol.Coords[0];
        // Lookup table for bond types so we don't calculate electrostatics
        // and such between bonded atoms
        request.BondType = makeBondTypeTable(mol);
        request.ResidueMap = makeResidueMap(mol);
        request.Decomp = make([]float64, mol.NumResidues()*mol.NumResidues());
    
        fmt.Println("Electrostatic energy:", Electro(&request), "kcal/mol");
        fmt.Println("van der Waals energy:", LennardJones(&request), "kcal/mol");
        fmt.Println(amber.Status());

        amber.DumpFloat64MatrixAsText(request.Decomp, mol.NumResidues(), "decomp.txt");
    } else if trjFilename != "" {
        // Or, do the trajectory.
        trjFp, err := os.Open(trjFilename, os.O_RDONLY, 0);
        if err != nil {
            fmt.Println("Error opening", trjFilename, err);
            return;
        }
        defer trjFp.Close();
        trj := bufio.NewReader(trjFp);
        trj.ReadString('\n'); // Eat header line

        if numFrames == 0 {
            fmt.Println("Please specify the number of frames you want processed with -n.");
            return;
        }
        fmt.Printf("Calculating energies for %d frames of trajectory %s.\n", numFrames, trjFilename);

        // Lookup table for bond types so we don't calculate nonbonded energies
        // between bonded atoms
        bondType := makeBondTypeTable(mol);
        residueMap := makeResidueMap(mol);
        
        var numKids int;
        ch := make(chan int);
        decompCh := make(chan *EnergyCalcRequest, 10);
        // This goroutine will be fed the decomposition matrices made by the energy functions
        fmt.Println("Writing residue decomposition matrices to", outFilename);
        go decompProcessor(outFilename, mol.NumResidues(), numFrames, decompCh, ch);
        numAtoms := mol.NumAtoms();
        hasBox := false;
        if mol.GetInt("POINTERS", amber.IFBOX) > 0 { hasBox = true }
        
        for frame := 0; frame < numFrames; frame++ {
            //coords := amber.GetFrameFromTrajectory(trjFilename, frame, numAtoms, hasBox);
            coords := amber.GetNextFrameFromTrajectory(trj, numAtoms, hasBox);
            go calcSingleTrjFrame(mol, params, coords, frame, bondType, residueMap, decompCh, ch);
            numKids++;
        }
        
        for i := 0; i < numKids; i++ { <-ch }
        decompCh <- nil;
        <-ch; // Wait for decompProcessor to finish
    }
}

// Does something with each decomposition matrix, which is currently writing them to disk.
// This is a separate goroutine so that only one matrix is processed at a time, which is
// convenient for writing to a disk.
// XXX: Matrices are written out of order because we receive them in arbitrary order.
// That should be OK for the correlation analysis though.
func decompProcessor(filename string, numResidues, numFrames int, ch chan *EnergyCalcRequest, termCh chan int) {
    // decompTotal := make([]float64, numResidues*numResidues);
    // Output file
    outFile, _ := os.Open(filename, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0644);
    defer outFile.Close();
    tmp := make([]byte, numResidues*numResidues*4); // for converting to bytes
    for ; ; {
        request := <-ch;
        if request == nil { break }
        // Actually do stuff with the data here.
        // We could in theory do the correlation stuff now, but maybe we should
        // just write the frames to disk.
        // Dump to file. We have to explicitly convert to bytes. Yay.
        for j, n := range(request.Decomp) {
            binary.BigEndian.PutUint32(tmp[j*4:j*4+4], math.Float32bits(float32(n)));
        }
        outFile.Write(tmp);
        // Explicitly release our reference to the decomp matrix.
        // If we don't do this, the GC seems to assume that all these pointer
        // references still live inside this function and doesn't free them
        // even though they have gone out of scope and nothing else is holding them.
        // My theory is that it only checks for stuff it can free when pointer variables
        // are given new values or when all functions that touched the object have returned.
        // Who knows though...perhaps at a near-OOM state it would have tried harder
        // to free up some memory?
        request.Decomp = nil;
}
    fmt.Println("decompProcessor finished. I wrote to", filename);
    termCh <- 0; // Tell caller we're done
}

// Calculates the nonbonded energies for a single snapshot.
// Results are returned through reqOutCh.
func calcSingleTrjFrame(mol *amber.System, params NonbondedParamsCache, coords []float32, frame int,
        bondType []uint8, residueMap []int, reqOutCh chan *EnergyCalcRequest, ch chan int) {

    var request EnergyCalcRequest;
    request.Molecule = mol;
    request.Frame = frame;
    request.NBParams = params;
    request.Coords = coords;
    request.BondType = bondType;
    request.ResidueMap = residueMap;
    request.Decomp = make([]float64, mol.NumResidues()*mol.NumResidues());

/*  //DEBUG: print first few coordinates
    fmt.Printf("%d [%d]:", frame, len(coords));
    for i := 0; i < 6; i++ {
        fmt.Printf(" %f", coords[i]);
    }
    fmt.Println();
*/
    elec := Electro(&request);
    vdw := LennardJones(&request);
    if math.IsNaN(elec) || math.IsNaN(vdw) {
        fmt.Println("Weird energies. Does your trajectory have boxes but your prmtop doesn't, or vice versa?");
    }
    fmt.Printf("%d: Electrostatic: %f vdW: %f %s\n", frame, elec, vdw, amber.Status());
    
    // Send request to listening something that will probably average the decomp matrix
    // but could in theory do whatever it wants.
    reqOutCh <- &request;
    // Return frame ID through channel
    ch <- frame;
}


// This is probably a little unwieldy but I hope it's better than
// zillions of arguments to every energy calculation function
type EnergyCalcRequest struct {
    Molecule *amber.System;
    Frame int; // Frame ID
    Coords []float32;
    BondType []uint8; // Input: Bond type matrix
    ResidueMap []int; // Input: Maps atom id to residue id
    // Output: pairwise residue-residue interaction energies
    Decomp []float64; // Lots of math going on here so float64
    Energy float64;
    
    // Parameters for nonbonded energy calculations.
    // Charges, LJ coefficients, and so on
    NBParams NonbondedParamsCache;
}

// Place to stash preloaded parameters for nonbonded energy calculations
type NonbondedParamsCache struct {
    Ntypes int;
    NBIndices, AtomTypeIndices []int;
    LJ12, LJ6, Charges []float32;
}

// Computes the Lennard-Jones 6-12 energy
func LennardJones(request *EnergyCalcRequest) float64 {
    const VDW_14_SCALING_RECIP = 1/2.0;
    // NONBONDED_PARM_INDEX==ICO
    // ATOM_TYPE_INDEX=IAC
    // For NONBONDED_PARM_INDEX:
    //   index = ICO(NTYPES*(IAC(i)-1)+IAC(j))
    // If index is positive, this is an index into the
    // 6-12 parameter arrays (CN1 and CN2) otherwise it
    // is an index into the 10-12 parameter arrays (ASOL
    // and BSOL).
    // LENNARD_JONES_ACOEF=CN1
    // LENNARD_JONES_BCOEF=CN2
    /* // Check for 10-12 parameters
    for _, v := range(nbIndices) {
        if v < 0 {
            fmt.Println("10-12 L-J parameters found, which aren't supported.");
            return 0
        }
    } */
    mol := request.Molecule;
    coords := request.Coords;
    bondType := request.BondType;
    decomp := request.Decomp;
    residueMap := request.ResidueMap;
    ntypes := request.NBParams.Ntypes;
    nbIndices := request.NBParams.NBIndices; 
    atomTypeIndices := request.NBParams.AtomTypeIndices; 
    lj12 := request.NBParams.LJ12;
    lj6 := request.NBParams.LJ6;
    numAtoms := mol.NumAtoms();
    numResidues := mol.NumResidues();
    var energy float64;
    for atom_i := 0; atom_i < numAtoms; atom_i++ {
        // Get coordinates for atom i
        offs_i := atom_i*3;
        x0, y0, z0 := coords[offs_i], coords[offs_i+1], coords[offs_i+2];
        // Pulled some of the matrix indexing out of the inner loop
        nbparm_offs_i := ntypes*(atomTypeIndices[atom_i]-1);
        bondtype_offs_i := atom_i*numAtoms;
        i_res := residueMap[atom_i]; // Residue of atom i
    
        for atom_j := 0; atom_j < atom_i; atom_j++ {
            // Are these atoms connected by a bond or angle? If so, skip.
            thisBondType := bondType[bondtype_offs_i + atom_j];
            if thisBondType & (BOND | ANGLE) != 0 { continue }
            // Calculate distance reciprocals
            offs_j := atom_j*3;
            x1, y1, z1 := coords[offs_j], coords[offs_j+1], coords[offs_j+2];
            dx, dy, dz := x1-x0, y1-y0, z1-z0;
            distRecip := Invsqrt32(dx*dx + dy*dy + dz*dz);
            distRecip3 := distRecip*distRecip*distRecip;
            distRecip6 := distRecip3*distRecip3;

            // Locate L-J parameters for this atom pair
            index := nbIndices[nbparm_offs_i + atomTypeIndices[atom_j]-1] - 1;
            // A/r12 - C/r6
            thisEnergy := float64(lj12[index]*distRecip6*distRecip6 - lj6[index]*distRecip6);
            // Are these atoms 1-4 to each other? If so, divide the energy
            // by 2.0, as ff99 et al dictate.
            if thisBondType & DIHEDRAL != 0 { thisEnergy *= VDW_14_SCALING_RECIP }
            // Pairwise residue energy decomposition - symmetric
            decomp[i_res*numResidues + residueMap[atom_j]] += thisEnergy;
            decomp[i_res + residueMap[atom_j]*numResidues] += thisEnergy;
            energy += thisEnergy;
        }
    }
    request.Energy = energy;
    return energy;
}

// Calculates electrostatic interactions among all particles in an amber.System,
// according to the force field (e.g. don't include bonded atoms)
func Electro(request *EnergyCalcRequest) float64 {
    const COULOMB = 332.0636;
    const EEL_14_SCALING_RECIP = 1/1.2;
    mol := request.Molecule;
    coords := request.Coords;
    bondType := request.BondType;
    decomp := request.Decomp;
    residueMap := request.ResidueMap;
    charges := request.NBParams.Charges;
    if charges == nil {
        fmt.Println("Electro: bad CHARGE block");
        return 0
    }
    // Actually do the calculation
    var energy float64;
    numAtoms := len(coords) / 3;
    numResidues := mol.NumResidues();
    for atom_i := 0; atom_i < numAtoms; atom_i++ {
        // x y z x y z ...
        offs_i := atom_i*3;
        x0, y0, z0 := coords[offs_i], coords[offs_i+1], coords[offs_i+2];
        qi := charges[atom_i];
        i_res := residueMap[atom_i]; // Residue of atom i
        // Iterate over all atoms
        for atom_j := 0; atom_j < atom_i; atom_j++ {
            offs_j := atom_j*3;
            x1, y1, z1 := coords[offs_j], coords[offs_j+1], coords[offs_j+2];
            // Are these atoms connected by a bond or angle? If so, skip.
            thisBondType := bondType[atom_i*numAtoms+atom_j];
            if thisBondType & (BOND | ANGLE) != 0 { continue }
            // Use reciprocal sqrt to calculate distance
            dx, dy, dz := x1-x0, y1-y0, z1-z0;
            thisEnergy := float64((qi * charges[atom_j]) * Invsqrt32(dx*dx + dy*dy + dz*dz));
            // Are these atoms 1-4 to each other? If so, divide the energy
            // by 1.2, as ff99 et al dictate.
            if thisBondType & DIHEDRAL != 0 { thisEnergy *= EEL_14_SCALING_RECIP }
            // Pairwise residue energy decomposition
            decomp[i_res*numResidues + residueMap[atom_j]] += thisEnergy;
            decomp[i_res + residueMap[atom_j]*numResidues] += thisEnergy;
            energy += thisEnergy;
        }
    }
    request.Energy = energy;
    return energy; // Tell caller the resulting total energy
}

// Converts the RESIDUE_POINTER block into an array, index by atom, that
// provides the residue number (starting at 0) for each atom.
func makeResidueMap(mol *amber.System) []int {
    residueList := amber.VectorAsIntArray(mol.Blocks["RESIDUE_POINTER"]);
    residueMap := make([]int, mol.NumAtoms()); // len(residueList) == #resideus
    for res_i := 1; res_i < len(residueList); res_i++ {
        a := residueList[res_i-1]-1; // Fortran starts counting at 1
        b := residueList[res_i]-1;
        for i := a; i < b; i++ { residueMap[i] = res_i - 1 }
    }
    // RESIDUE_POINTER doesn't specify the last residue because that's implied
    numResidues := mol.NumResidues();
    for i := residueList[numResidues-1]-1; i < len(residueMap); i++ {
        residueMap[i] = numResidues-1;
    }
    return residueMap;
}

// We want to be able to quickly look up if two atoms are bonded.
// To do this, make a matrix for all atom pairs such that
// M[numAtoms*i+j] & bondtypeflag != 0
func makeBondTypeTable(mol *amber.System) []uint8 {
    numAtoms := mol.NumAtoms();
    bondType := make([]uint8, numAtoms*numAtoms);

    bondsBlocks := []string {"BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN"};
    for _, blockName := range(bondsBlocks) {
        bonds := amber.VectorAsIntArray(mol.Blocks[blockName]);
        // atom_i atom_j indexintostuff
        // These are actually coordinate array indices, not atom indices
        for i := 0; i < len(bonds); i+=3 {
            a, b := bonds[i]/3, bonds[i+1]/3;
            bondType[numAtoms*a+b] |= BOND;
            bondType[numAtoms*b+a] |= BOND;
        }
    }
    
    angleBlocks := []string {"ANGLES_WITHOUT_HYDROGEN", "ANGLES_INC_HYDROGEN"};
    for _, blockName := range(angleBlocks) {
        angles := amber.VectorAsIntArray(mol.Blocks[blockName]);
        // atom_i atom_j atom_k indexintostuff
        for i := 0; i < len(angles); i+=4 {
            a, b := angles[i]/3, angles[i+2]/3;
            bondType[numAtoms*a+b] |= ANGLE;
            bondType[numAtoms*b+a] |= ANGLE;
        }
    }
    
    dihedBlocks := []string {"DIHEDRALS_INC_HYDROGEN", "DIHEDRALS_WITHOUT_HYDROGEN"};
    for _, blockName := range(dihedBlocks) {
        diheds := amber.VectorAsIntArray(mol.Blocks[blockName]);
        // atom_i atom_j atom_k atom_l indexintostuff
        for i := 0; i < len(diheds); i+=5 {
            // Fourth coordinate index is negative if this is an improper
            a, b := amber.Abs(diheds[i]/3), amber.Abs(diheds[i+3]/3);
            bondType[numAtoms*a+b] |= DIHEDRAL;
            bondType[numAtoms*b+a] |= DIHEDRAL;
        }
    }
    return bondType;    
}

// Flags for the bondType matrix
const (
    UNBONDED = 0;
    BOND = 1;
    ANGLE = 2;
    DIHEDRAL = 4;
);