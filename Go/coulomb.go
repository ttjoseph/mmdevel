package main
import (
    "fmt";
    "amber";
)

func main() {
    fmt.Println("Hello there.");

    mol := amber.LoadSystem("complex.top.x", "complex.crd.1627");
    if mol == nil { return }
    
    fmt.Println("Number of atoms:", mol.NumAtoms());
    fmt.Println("Number of residues:", mol.NumResidues());
    fmt.Println("Number of frames:", len(mol.Coords));
    fmt.Println("Total energy:", Electro(mol));
}

// Returns the Euclidean distance between two 3D points
func distance(x0, y0, z0, x1, y1, z1 float32) float32 {
    dx := x1-x0;
    dy := y1-y0;
    dz := z1-z0;
    return Sqrt32(dx*dx + dy*dy + dz*dz)
}

// Calculates electrostatic interactions among all particles in an AmberSystem
func Electro(mol *amber.System) float32 {
    // Goroutines are cheap so we can have a small blocksize...I guess
    const blockSize = 2000;
    charges := amber.VectorAsFloat32Array(mol.Blocks["CHARGE"]);
    if charges == nil {
        fmt.Println("Electro: bad CHARGE block");
        return 0
    }
    fmt.Println("Electro: Starting calculation.");
    ch := make(chan float32);
    numKids := 0;
    for i := 0; i < mol.NumAtoms(); i += blockSize {
        // Fire a goroutine to calculate part of the answer
        a, b := i, i + blockSize; 
        if b > mol.NumAtoms() { b = mol.NumAtoms() }
        go calcElecPiece(mol.Coords[0], charges, a, b, ch);
        numKids++;
    }
    
    // Collect results
    var energy float32 = 0;
    for i := 0; i < numKids; i++ { energy += <-ch }
    return energy;
}

// Returns, through the channel, the electrostatic energy of the system, only considering
// atom indices [a,b).
func calcElecPiece(coords, charges []float32, a, b int, ch chan float32) {
    const COULOMB = 332.0636;
    var energy float32 = 0.0;
    //fmt.Println("calcElecPiece", a, b);
    for atom_i := a; atom_i < b; atom_i++ {
        // x y z x y z ...
        offs_i := atom_i*3;
        x0 := coords[offs_i];
        y0 := coords[offs_i+1];
        z0 := coords[offs_i+2];
        qi := charges[atom_i];
        // Iterate over all atoms
        for atom_j := 0; atom_j < atom_i; atom_j++ {
            offs_j := atom_j*3;
            x1 := coords[offs_j];
            y1 := coords[offs_j+1];
            z1 := coords[offs_j+2];
            // Use reciprocal sqrt to calculate distance
            dx := x1-x0;
            dy := y1-y0;
            dz := z1-z0;
            energy += (qi * charges[atom_j]) * Invsqrt32(dx*dx + dy*dy + dz*dz);
        }
    }
    ch <- energy; // Tell caller the resulting total energy
}
