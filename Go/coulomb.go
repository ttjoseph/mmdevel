package main
import (
    "fmt";
    "flag";
    "amber";
)

func main() {
    var prmtopFilename, rstFilename string;
    flag.StringVar(&prmtopFilename, "p", "prmtop", "Prmtop filename");
    flag.StringVar(&rstFilename, "c", "inpcrd", "Inpcrd/rst filename");
    flag.Parse();

    mol := amber.LoadSystem(prmtopFilename, rstFilename);
    if mol == nil { return }
    // Lookup table for bond types so we don't calculate electrostatics
    // and such between bonded atoms
    bondType := makeBondTypeTable(mol);
    
    fmt.Println("Number of atoms:", mol.NumAtoms());
    fmt.Println("Number of residues:", mol.NumResidues());
    fmt.Println("Electrostatic energy:", Electro(mol, bondType), "kcal/mol");
    fmt.Println("van der Waals energy:", LennardJones(mol, bondType), "kcal/mol");
}

// Computes the Lennard-Jones 6-12 energy
func LennardJones(mol *amber.System, bondType []uint8) float32 {
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
    ntypes := mol.GetInt("POINTERS", amber.NTYPES);
    coords := mol.Coords[0];
    nbIndices := amber.VectorAsIntArray(mol.Blocks["NONBONDED_PARM_INDEX"]); // ICO
    atomTypeIndices := amber.VectorAsIntArray(mol.Blocks["ATOM_TYPE_INDEX"]); // IAC
    lj12 := amber.VectorAsFloat32Array(mol.Blocks["LENNARD_JONES_ACOEF"]); // CN1
    lj6 := amber.VectorAsFloat32Array(mol.Blocks["LENNARD_JONES_BCOEF"]); // CN2
    numAtoms := mol.NumAtoms();
    var energy float32;
    for atom_i := 0; atom_i < numAtoms; atom_i++ {
        // Get coordinates for atom i
        offs_i := atom_i*3;
        x0, y0, z0 := coords[offs_i], coords[offs_i+1], coords[offs_i+2];
    
        for atom_j := 0; atom_j < atom_i; atom_j++ {
            // Are these atoms connected by a bond or angle? If so, skip.
            thisBondType := bondType[atom_i*numAtoms+atom_j];
            if thisBondType & (BOND | ANGLE) != 0 { continue }
            // Calculate distance reciprocals
            offs_j := atom_j*3;
            x1, y1, z1 := coords[offs_j], coords[offs_j+1], coords[offs_j+2];
            dx, dy, dz := x1-x0, y1-y0, z1-z0;
            distRecip := Invsqrt32(dx*dx + dy*dy + dz*dz);
            distRecip3 := distRecip*distRecip*distRecip;
            distRecip6 := distRecip3*distRecip3;
            distRecip12 := distRecip6*distRecip6;

            // Locate L-J parameters for this atom pair
            index := nbIndices[ntypes*(atomTypeIndices[atom_i]-1) + atomTypeIndices[atom_j]-1] - 1;
            if index < 0 {
                fmt.Println("Found 10-12 Lennard-Jones parameters, which aren't supported.");
                return 0;
            }
            // A/r12 - C/r6
            k12, k6 := lj12[index], lj6[index];
            thisEnergy := k12*distRecip12 - k6*distRecip6;
            // Are these atoms 1-4 to each other? If so, divide the energy
            // by 2.0, as ff99 et al dictate.
            if thisBondType & DIHEDRAL != 0 { thisEnergy *= VDW_14_SCALING_RECIP }
            energy += thisEnergy;
        }
    }
    
    return energy;
}

// Calculates electrostatic interactions among all particles in an AmberSystem
func Electro(mol *amber.System, bondType []uint8) float32 {
    // Goroutines are cheap so we can have a small blocksize...I guess
    const blockSize = 2000;
    charges := amber.VectorAsFloat32Array(mol.Blocks["CHARGE"]);
    if charges == nil {
        fmt.Println("Electro: bad CHARGE block");
        return 0
    }
    
    // Actually do the calculation
    ch := make(chan float32);
    numKids := 0;
    for i := 0; i < mol.NumAtoms(); i += blockSize {
        // Fire a goroutine to calculate part of the answer
        a, b := i, i + blockSize; 
        if b > mol.NumAtoms() { b = mol.NumAtoms() }
        go calcElecPiece(mol.Coords[0], charges, bondType, a, b, ch);
        numKids++;
    }
    
    // Collect results
    var energy float32 = 0;
    for i := 0; i < numKids; i++ { energy += <-ch }
    return energy;
}

// Returns, through the channel, the electrostatic energy of the system, only considering
// atom indices [a,b).
func calcElecPiece(coords, charges []float32, bondType []uint8, a, b int, ch chan float32) {
    const COULOMB = 332.0636;
    const EEL_14_SCALING_RECIP = 1/1.2;
    var energy float32 = 0.0;
    numAtoms := len(coords) / 3;
    for atom_i := a; atom_i < b; atom_i++ {
        // x y z x y z ...
        offs_i := atom_i*3;
        x0, y0, z0 := coords[offs_i], coords[offs_i+1], coords[offs_i+2];
        qi := charges[atom_i];
        // Iterate over all atoms
        for atom_j := 0; atom_j < atom_i; atom_j++ {
            offs_j := atom_j*3;
            x1, y1, z1 := coords[offs_j], coords[offs_j+1], coords[offs_j+2];
            // Are these atoms connected by a bond or angle? If so, skip.
            thisBondType := bondType[atom_i*numAtoms+atom_j];
            if thisBondType & (BOND | ANGLE) != 0 { continue }
            // Use reciprocal sqrt to calculate distance
            dx, dy, dz := x1-x0, y1-y0, z1-z0;
            thisEnergy := (qi * charges[atom_j]) * Invsqrt32(dx*dx + dy*dy + dz*dz);
            // Are these atoms 1-4 to each other? If so, divide the energy
            // by 1.2, as ff99 et al dictate.
            if thisBondType & DIHEDRAL != 0 { thisEnergy *= EEL_14_SCALING_RECIP }
            
            // TODO: Add to residue-residue energy bucket for decomposition
            
            energy += thisEnergy;
        }
    }
    ch <- energy; // Tell caller the resulting total energy
}

// Returns the Euclidean distance between two 3D points
func distance(x0, y0, z0, x1, y1, z1 float32) float32 {
    dx, dy, dz := x1-x0, y1-y0, z1-z0;
    return Sqrt32(dx*dx + dy*dy + dz*dz)
}

func makeBondTypeTable(mol *amber.System) []uint8 {
    // We want to be able to quickly look up if two atoms are bonded.
    // To do this, make a matrix for all atom pairs such that
    // M[numAtoms*i+j] == bondtype (one of the constants).
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

const (
    UNBONDED = 0;
    BOND = 1;
    ANGLE = 2;
    DIHEDRAL = 4;
);