// Utility to delete water from inside lipid bilayer
// Caution! This is full of hardcoded parameters. Please read through before running.
package main

import (
    "fmt"
    "os"
    "bufio"
    "strconv"
    "log"
    "flag"
    "strings"
    "math"
    "container/list"
)

type Atom struct {
    id int
    atomType string
    x, y, z, vx, vy, vz float64
}

type Residue struct {
    name string
    num int
    atoms []Atom    
}

type Gro struct {
    description string
    numAtoms int
    residues *list.List
    dx, dy, dz float64 // box dimensions
    hasVelocities bool
}


func main() {
    log.Println("Move waters out of lipid bilayer - thomas.joseph@mountsinai.org")
    var groFile string
    var excludePore bool
    var poreX, poreY, poreRadius float64
    flag.StringVar(&groFile, "gro", "conf.gro", "the GROMACS .gro file to read")
    flag.BoolVar(&excludePore, "excludepore", false, "whether we should leave waters in the pore, as defined by a cylinder parallel to Z axis")
    flag.Float64Var(&poreX, "porex", 5.9, "X coordinate (nm) of line defining center of pore axis")
    flag.Float64Var(&poreY, "porey", 5.9, "Y coordinate (nm) of line defining center of pore axis")
    flag.Float64Var(&poreRadius, "poreradius", 1.3, "Pore radius in nm")
    flag.Parse()
    log.Printf(".gro file: %s\n", groFile)
    
    // read gro file
    f, err := os.Open(groFile)
    if err != nil {
        log.Fatal("Error opening gro file ", groFile)
    }
    gro := readGro(bufio.NewReader(f))
    
    scrubBilayer(gro, "POPC", "C2", excludePore, poreX, poreY, poreRadius)
    gro.print(os.Stdout)
}

// atomType specifies the atom in the headgroup which defines the extents of the bilayer
func scrubBilayer(gro *Gro, bilayerResname string, atomType string, excludePore bool, poreX, poreY, poreRadius float64) {
    // Find Z extents of bilayer residue to define region that solvent and ions should not be in
    minZ, maxZ := 999999999999.0, -99999999999.0
    allMinZ, allMaxZ := minZ, maxZ
    for res := gro.residues.Front(); res != nil; res = res.Next() {
        r := res.Value.(Residue)
        if r.name == bilayerResname {
            for _, a := range r.atoms {
                if a.atomType == atomType {
                    if a.z > maxZ {
                        maxZ = a.z
                    }
                    if a.z < minZ {
                        minZ = a.z
                    }
                }
            }
        }
        // Also find real extents of the system
        for _, a := range r.atoms {
            if a.z > allMaxZ {
                allMaxZ = a.z
            }
            if a.z < allMinZ {
                allMinZ = a.z
            }
        }
    }
    
    // Fudge factor to leave some waters by the headgroups
    headgroupFudge := 0.5
    maxZ -= headgroupFudge
    minZ += headgroupFudge
    
    log.Println("Lipid bilayer extents by", bilayerResname, "are", minZ, "nm to", maxZ, "nm along Z axis")
    midZ := (maxZ - minZ) / 2 + minZ
    watersToDelete, numWaters := 0, 0
    // We manually keep track of the next residue in our iteration because we might move the current residue
    // to the front of the list.
    nextRes := gro.residues.Front()
    for res := gro.residues.Front(); res != nil; res = nextRes {
        r := res.Value.(Residue)
        nextRes = res.Next()
        // Shift ions out of bilayer
        if r.name == "NA" || r.name == "CL" {
            shiftZ := 0.0
            for _, a := range r.atoms {
                if a.z >= minZ && a.z <= midZ {
                    shiftZ = (minZ - 0.7) - a.z
                } else if a.z > midZ && a.z <= maxZ {
                    shiftZ = (maxZ + 0.7) - a.z
                }
            }
            
            if shiftZ != 0.0 {
                log.Printf("Shifting %d%s by %.3f\n", r.num, r.name, shiftZ)
                for i, a := range r.atoms {
                    a.z += shiftZ
                    r.atoms[i] = a
                }
            }
        }
        
        // Delete waters from bilayer by moving them to the front of the residue list and then chopping them off.
        // We do this because apparently this doubly-linked list structure can't handle having elements removed during
        // an iteration.
        if r.name == "SOL" {
            numWaters++
            for _, a:= range r.atoms {
                if a.z>=minZ && a.z<maxZ && (excludePore==false || (excludePore && !isWithinPore(a, poreX, poreY, poreRadius))) {
                    gro.residues.MoveToFront(res)
                    watersToDelete++
                    break
                }
            }
        }
    }
    
    for i:=0; i<watersToDelete; i++ {
        gro.residues.Remove(gro.residues.Front())
    }
    gro.calcNumAtoms()
    log.Println("Deleted", watersToDelete, "solvent residues,", numWaters - watersToDelete, "remain - don't forget to edit your topol.top (or whatever it's called)")
}

func isWithinPore(a Atom, x, y, r float64) bool {
    dist := math.Hypot(x-a.x, y-a.y)
    // if dist > r {
    //     log.Println(a.x, a.y, "Distance", dist, "r", r)
    // }
    return dist <= r
}

func printResidueWithVelocities(w *os.File, r *Residue) {
    if r == nil {
        return
    }
    for _, a := range r.atoms {
        fmt.Fprintf(w, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", r.num, r.name, a.atomType, a.id, a.x, a.y, a.z, a.vx, a.vy, a.vz)
    }
}

func printResidueWithoutVelocities(w *os.File, r *Residue) {
    if r == nil {
        return
    }
    for _, a := range r.atoms {
        fmt.Fprintf(w, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", r.num, r.name, a.atomType, a.id, a.x, a.y, a.z)
    }
}

func (gro *Gro) print(w *os.File) {
    fmt.Fprintf(w, "%s\n%d\n", gro.description, gro.numAtoms)
    for res := gro.residues.Front(); res != nil; res = res.Next() {
        r := res.Value.(Residue)
        if gro.hasVelocities {
            printResidueWithVelocities(w, &r)
        } else {
            printResidueWithoutVelocities(w, &r)
        }
    }
    fmt.Fprintf(w, "%f %f %f\n", gro.dx, gro.dy, gro.dz)
}

func readGro(r *bufio.Reader) (gro *Gro) {
    gro = newGro()
    scanner := bufio.NewScanner(r)
    scanner.Scan()
    gro.description = scanner.Text()
    scanner.Scan()
    var err error
    gro.numAtoms, err = strconv.Atoi(scanner.Text())
    if err != nil {
        log.Fatal("Couldn't parse number of atoms in .gro header: ", err)
    }
    atomCount := 0
    var currentRes *Residue
    currentRes = nil
    for scanner.Scan() {
        // parse residues and atoms
        var a Atom
        var resname string
        var resnum int
        s := scanner.Text()
        resnum, err = strconv.Atoi(strings.TrimSpace(s[0:5]))
        resname = strings.TrimSpace(s[5:10])
        a.atomType = strings.TrimSpace(s[10:15])
        a.id, err = strconv.Atoi(strings.TrimSpace(s[15:20]))
        a.x, err = strconv.ParseFloat(strings.TrimSpace(s[20:28]), 64)
        a.y, err = strconv.ParseFloat(strings.TrimSpace(s[28:36]), 64)
        a.z, err = strconv.ParseFloat(strings.TrimSpace(s[36:44]), 64)
        if len(s) > 44 {
            a.vx, err = strconv.ParseFloat(strings.TrimSpace(s[44:52]), 64)
            a.vy, err = strconv.ParseFloat(strings.TrimSpace(s[52:60]), 64)
            a.vz, err = strconv.ParseFloat(strings.TrimSpace(s[60:68]), 64)
            gro.hasVelocities = true
        }
        if err != nil {
            log.Print(scanner.Text())
            log.Fatal("Failed parsing atom line: ", err)
        }
        
        if currentRes == nil || currentRes.num != resnum {
            if currentRes != nil {
                gro.residues.PushBack(*currentRes)
            }
            currentRes = new(Residue)
            currentRes.num = resnum
            currentRes.name = resname
        }
        
        currentRes.atoms = append(currentRes.atoms, a)
        
        atomCount++
        if atomCount >= gro.numAtoms {
            break
        }
    }
    // Don't forget to tack that last residue on the end there
    gro.residues.PushBack(*currentRes)

    scanner.Scan()
    box := strings.Fields(scanner.Text())
    var err2, err3 error
    gro.dx, err = strconv.ParseFloat(box[0], 64)
    gro.dy, err2 = strconv.ParseFloat(box[1], 64)
    gro.dz, err3 = strconv.ParseFloat(box[2], 64)
    if err != nil || err2 != nil || err3 != nil {
        log.Fatal("Error parsing box size")
    }
    return gro
}

func (gro *Gro) calcNumAtoms() (n int) {
    for res := gro.residues.Front(); res != nil; res = res.Next() {
        r := res.Value.(Residue)
        n += len(r.atoms)
    }
    gro.numAtoms = n
    return n
}

func newGro() (*Gro) {
    gro := new(Gro)
    gro.residues = list.New()
    gro.hasVelocities = false
    return gro
}
