// Calculates the nonbonded energy (vdW and electrostatic) in an AMBER system.
// Assumes no cutoff. Does not calculate any other terms.
package main

import (
	"amber"
	"encoding/binary"
	"flag"
	"fmt"
	"math"
	"os"
)

func WriteInt32(file *os.File, d int) {
	tmp := make([]uint8, 4)
	binary.LittleEndian.PutUint32(tmp[0:4], uint32(d))
	file.Write(tmp)
}

func WriteInt32Array(file *os.File, d []int) {
	WriteInt32(file, len(d)) // Write size of array to file, then array itself
	tmp := make([]uint8, len(d)*4)
	for j, n := range d {
		binary.LittleEndian.PutUint32(tmp[j*4:j*4+4], uint32(n))
	}
	file.Write(tmp)
}

func WriteInt8Array(file *os.File, d []uint8) {
	WriteInt32(file, len(d)) // Write size of array to file, then array itself
	file.Write(d)
}

func WriteFloat32Array(file *os.File, d []float32) {
	WriteInt32(file, len(d)) // Write size of array to file, then array itself
	tmp := make([]uint8, len(d)*4)
	for j, n := range d {
		binary.LittleEndian.PutUint32(tmp[j*4:j*4+4], math.Float32bits(float32(n)))
	}
	file.Write(tmp)
}

func main() {
	var prmtopFilename, rstFilename, outFilename string
	var stride int
	var savePreprocessed bool

	flag.StringVar(&prmtopFilename, "p", "prmtop", "Prmtop filename (required)")
	flag.StringVar(&rstFilename, "c", "", "Inpcrd/rst filename")
	flag.IntVar(&stride, "s", 1, "Frame stride; 1 = don't skip any")
	flag.StringVar(&outFilename, "o", "energies.bin", "Energy decomposition output filename")
	flag.BoolVar(&savePreprocessed, "e", false, "Save prmtop preprocessed output (use with -c)")
	flag.Parse()

	mol := amber.LoadSystem(prmtopFilename)
	if mol == nil {
		return
	}
	fmt.Println("Number of atoms:", mol.NumAtoms())
	fmt.Println("Number of residues:", mol.NumResidues())
	hasBox := false
	fmt.Print("Periodic box in prmtop: ")
	if mol.GetInt("POINTERS", amber.IFBOX) > 0 {
		hasBox = true
		fmt.Println("Yes")
	} else {
		fmt.Println("No")
	}

	// Set up nonbonded parameters. We load them here so we don't have to keep
	// doing it later
	var params NonbondedParamsCache
	params.Ntypes = mol.GetInt("POINTERS", amber.NTYPES)
	params.NBIndices = amber.VectorAsIntArray(mol.Blocks["NONBONDED_PARM_INDEX"])  // ICO
	params.AtomTypeIndices = amber.VectorAsIntArray(mol.Blocks["ATOM_TYPE_INDEX"]) // IAC
	params.LJ12 = amber.VectorAsFloat32Array(mol.Blocks["LENNARD_JONES_ACOEF"])    // CN1
	params.LJ6 = amber.VectorAsFloat32Array(mol.Blocks["LENNARD_JONES_BCOEF"])     // CN2
	params.Charges = amber.VectorAsFloat32Array(mol.Blocks["CHARGE"])

	var request EnergyCalcRequest
	request.NBParams = params
	request.Molecule = mol
	// Lookup table for bond types so we don't calculate electrostatics
	// and such between bonded atoms
	request.BondType = makeBondTypeTable(mol)
	request.ResidueMap = makeResidueMap(mol)
	request.Decomp = make([]float64, mol.NumResidues()*mol.NumResidues())

	// Dump the preprocessed info to a file so a C version of this program can easily load and parse it
	outFile, _ := os.OpenFile("solute.top.tom", os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0644)
	defer outFile.Close()

	WriteInt32(outFile, mol.NumAtoms())
	WriteInt32(outFile, mol.NumResidues())
	WriteInt32(outFile, params.Ntypes)
	WriteInt32Array(outFile, params.NBIndices)
	WriteInt32Array(outFile, params.AtomTypeIndices)
	WriteFloat32Array(outFile, params.LJ12)
	WriteFloat32Array(outFile, params.LJ6)
	WriteFloat32Array(outFile, params.Charges)
	WriteInt8Array(outFile, request.BondType)
	WriteInt32Array(outFile, request.ResidueMap)

	if hasBox {
		WriteInt32(outFile, 1)
	} else {
		WriteInt32(outFile, 0)
	}
	fmt.Println("Wrote preprocessed prmtop data. Done!")
	os.Exit(0)
}

// Converts the RESIDUE_POINTER block into an array, index by atom, that
// provides the residue number (starting at 0) for each atom.
func makeResidueMap(mol *amber.System) []int {
	residueList := amber.VectorAsIntArray(mol.Blocks["RESIDUE_POINTER"])
	residueMap := make([]int, mol.NumAtoms()) // len(residueList) == #resideus
	for res_i := 1; res_i < len(residueList); res_i++ {
		a := residueList[res_i-1] - 1 // Fortran starts counting at 1
		b := residueList[res_i] - 1
		for i := a; i < b; i++ {
			residueMap[i] = res_i - 1
		}
	}
	// RESIDUE_POINTER doesn't specify the last residue because that's implied
	numResidues := mol.NumResidues()
	for i := residueList[numResidues-1] - 1; i < len(residueMap); i++ {
		residueMap[i] = numResidues - 1
	}
	return residueMap
}

// We want to be able to quickly look up if two atoms are bonded.
// To do this, make a matrix for all atom pairs such that
// M[numAtoms*i+j] & bondtypeflag != 0
func makeBondTypeTable(mol *amber.System) []uint8 {
	numAtoms := mol.NumAtoms()
	bondType := make([]uint8, numAtoms*numAtoms)

	bondsBlocks := []string{"BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN"}
	for _, blockName := range bondsBlocks {
		bonds := amber.VectorAsIntArray(mol.Blocks[blockName])
		// atom_i atom_j indexintostuff
		// These are actually coordinate array indices, not atom indices
		for i := 0; i < len(bonds); i += 3 {
			a, b := bonds[i]/3, bonds[i+1]/3
			bondType[numAtoms*a+b] |= BOND
			bondType[numAtoms*b+a] |= BOND
		}
	}

	angleBlocks := []string{"ANGLES_WITHOUT_HYDROGEN", "ANGLES_INC_HYDROGEN"}
	for _, blockName := range angleBlocks {
		angles := amber.VectorAsIntArray(mol.Blocks[blockName])
		// atom_i atom_j atom_k indexintostuff
		for i := 0; i < len(angles); i += 4 {
			a, b := angles[i]/3, angles[i+2]/3
			bondType[numAtoms*a+b] |= ANGLE
			bondType[numAtoms*b+a] |= ANGLE
		}
	}

	dihedBlocks := []string{"DIHEDRALS_INC_HYDROGEN", "DIHEDRALS_WITHOUT_HYDROGEN"}
	for _, blockName := range dihedBlocks {
		diheds := amber.VectorAsIntArray(mol.Blocks[blockName])
		// atom_i atom_j atom_k atom_l indexintostuff
		for i := 0; i < len(diheds); i += 5 {
			// Fourth coordinate index is negative if this is an improper
			a, b := amber.Abs(diheds[i]/3), amber.Abs(diheds[i+3]/3)
			bondType[numAtoms*a+b] |= DIHEDRAL
			bondType[numAtoms*b+a] |= DIHEDRAL
		}
	}
	return bondType
}

// Flags for the bondType matrix
const (
	UNBONDED = 0
	BOND     = 1
	ANGLE    = 2
	DIHEDRAL = 4
)

// This is probably a little unwieldy but I hope it's better than
// zillions of arguments to every energy calculation function
type EnergyCalcRequest struct {
	Molecule   *amber.System
	Frame      int // Frame ID
	Coords     []float32
	BondType   []uint8 // Input: Bond type matrix
	ResidueMap []int   // Input: Maps atom id to residue id
	// Output: pairwise residue-residue interaction energies
	Decomp []float64 // Lots of math going on here so float64
	Energy float64

	// Parameters for nonbonded energy calculations.
	// Charges, LJ coefficients, and so on
	NBParams NonbondedParamsCache
}

// Place to stash preloaded parameters for nonbonded energy calculations
type NonbondedParamsCache struct {
	Ntypes                     int
	NBIndices, AtomTypeIndices []int
	LJ12, LJ6, Charges         []float32
}
