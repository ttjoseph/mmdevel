// Functions for loading AMBER output files, as well as random other stuff
// Tom Joseph <thomas.joseph@mssm.edu>
package amber

import (
	"fmt"
	"os"
	"bufio"
	"io/ioutil"
	"regexp"
	"strconv"
	"strings"
	"container/vector"
	"compress/gzip"
	"malloc"
)
import (
	"what"
)

// Constants for POINTERS block indices in prmtops, stolen from the prmtop
// format spec on ambermd.org
const (
	NATOM  = 0  // total number of atoms
	NTYPES = 1  // total number of distinct atom types
	NBONH  = 2  // number of bonds containing hydrogen
	MBONA  = 3  // number of bonds not containing hydrogen
	NTHETH = 4  // number of angles containing hydrogen
	MTHETA = 5  // number of angles not containing hydrogen
	NPHIH  = 6  // number of dihedrals containing hydrogen
	MPHIA  = 7  // number of dihedrals not containing hydrogen
	NHPARM = 8  // currently not used
	NPARM  = 9  // currently not used
	NEXT   = 10 // number of excluded atoms
	NRES   = 11 // number of residues
	NBONA  = 12 // MBONA + number of constraint bonds
	NTHETA = 13 // MTHETA + number of constraint angles
	NPHIA  = 14 // MPHIA + number of constraint dihedrals
	NUMBND = 15 // number of unique bond types
	NUMANG = 16 // number of unique angle types
	NPTRA  = 17 // number of unique dihedral types
	NATYP  = 18 // number of atom types in parameter file, see SOLTY below
	NPHB   = 19 // number of distinct 10-12 hydrogen bond pair types
	IFPERT = 20 // set to 1 if perturbation info is to be read in
	NBPER  = 21 // number of bonds to be perturbed
	NGPER  = 22 // number of angles to be perturbed
	NDPER  = 23 // number of dihedrals to be perturbed
	MBPER  = 24 // number of bonds with atoms completely in perturbed group
	MGPER  = 25 // number of angles with atoms completely in perturbed group
	MDPER  = 26 // number of dihedrals with atoms completely in perturbed groups
	IFBOX  = 27 // set to 1 if standard periodic box, 2 when truncated octahedral
	NMXRS  = 28 // number of atoms in the largest residue
	IFCAP  = 29 // set to 1 if the CAP option from edit was specified
)

// Converts a Vector of strings to an array of float32.
func VectorAsFloat32Array(v *vector.Vector) []float32 {
	data := make([]float32, v.Len())
	for i := 0; i < v.Len(); i++ {
		x := strings.TrimSpace(v.At(i).(string))
		data[i] = what.Atof32(x)
	}
	return data
}

// Converts a Vector of strings to an array of int.
func VectorAsIntArray(v *vector.Vector) []int {
	data := make([]int, v.Len())
	for i := 0; i < v.Len(); i++ {
		x := strings.TrimSpace(v.At(i).(string))
		data[i] = what.Atoi(x)
	}
	return data
}

// Converts an array of strings to an array of float32.
func StringArrayAsFloat32Array(v []string) []float32 {
	data := make([]float32, len(v))
	for i := 0; i < len(v); i++ {
		x := strings.TrimSpace(v[i])
		data[i] = what.Atof32(x)
	}
	return data
}

// Encapsulates an AMBER prmtop/inpcrd
type System struct {
	Blocks  map[string]*vector.Vector
	Formats map[string][]string
	Coords  map[int][]float32
}

// Gets an int at specified offset in specified block.
func (mol *System) GetInt(blockName string, index int) int {
	s := mol.Blocks[blockName].At(index).(string)
	val := what.Atoi(s)
	return val
}

// Returns the number of atoms the prmtop expects there to be
func (mol *System) NumAtoms() int {
	x := mol.Blocks["POINTERS"].At(0).(string)
	val := what.Atoi(strings.TrimSpace(x))
	return val
}

// Returns the number of residues the prmtop expects there to be
func (mol *System) NumResidues() int {
	x := mol.Blocks["POINTERS"].At(11).(string)
	val := what.Atoi(strings.TrimSpace(x))
	return val
}

// Loads an AMBER system - both a prmtop and inpcrd.
func LoadSystem(prmtopFilename string) *System {
	var mol System
	mol.Blocks = make(map[string]*vector.Vector)
	mol.Formats = make(map[string][]string)

	prmtopFile, err := os.Open(prmtopFilename, os.O_RDONLY, 0)
	if err != nil {
		fmt.Println("Error opening prmtop:", err)
		return nil
	}
	// defer == awesome!
	defer prmtopFile.Close()
	prmtop := bufio.NewReader(prmtopFile)

	// Set up us the regular expression
	formatRe, _ := regexp.Compile("[(]([0-9]+)([aIE]+)([0-9.]+)[)]")

	// Eat header line
	s, err := prmtop.ReadString('\n')
	//fmt.Println(s);

	// Get first FLAG line
	s, err = prmtop.ReadString('\n')

	for err == nil {
		//fmt.Println(s);
		if len(s) < 5 {
			break
		}
		blockName := strings.TrimSpace(s[6 : len(s)-1])
		mol.Blocks[blockName] = new(vector.Vector)

		// FORMAT line
		s, err = prmtop.ReadString('\n')
		fmtSpec := formatRe.MatchStrings(s)

		var numThings, thingLen int
		if len(fmtSpec) == 4 {
			numThings, _ = strconv.Atoi(fmtSpec[1])
			//thingType = fmtSpec[2];
			lenSpec, _ := strconv.Atof32(fmtSpec[3])
			thingLen = int(lenSpec)
			mol.Formats[blockName] = fmtSpec
		} else {
			fmt.Println("Couldn't understand format specifier - bad prmtop")
			return nil
		}

		// Now that we know what type of data to expect, read it
		// until we hit another FLAG line
		for err == nil {
			s, err = prmtop.ReadString('\n')
			if len(s) >= 5 && s[0:5] == "%FLAG" {
				break
			}
			// s contains some number of things
			n := (len(s) - 1) / thingLen
			if n > numThings {
				n = numThings
			}
			for i := 0; i < (n * thingLen); i += thingLen {
				mol.Blocks[blockName].Push(s[i : i+thingLen])
			}
		}
	}
	return &mol
}

// Load an inpcrd file into this system
func (mol *System) LoadRst(inpcrdFilename string) {
	inpcrdFile, err := os.Open(inpcrdFilename, os.O_RDONLY, 0)
	if err != nil {
		fmt.Println("Error opening inpcrd:", err)
		return
	}
	// defer == awesome!
	defer inpcrdFile.Close()
	inpcrd := bufio.NewReader(inpcrdFile)
	// Header line - discard
	s, err := inpcrd.ReadString('\n')
	// Number of atoms
	s, err = inpcrd.ReadString('\n')
	n, _ := strconv.Atoi(strings.TrimSpace(s))
	numCoords := mol.NumAtoms() * 3
	if n*3 != numCoords {
		fmt.Fprintf(os.Stderr, "Inpcrd says it has %d atoms but I'm expecting %d instead.\n", n, numCoords)
		return
	}

	// Now we know how many tokens to expect. They are all of length 12,
	// with up to 6 per line.
	mol.Coords = make(map[int][]float32)
	frame := 0
	mol.Coords[frame] = make([]float32, numCoords)
	numThings := 6
	thingLen := 12
	ctr := 0

	for err == nil {
		s, err = inpcrd.ReadString('\n')
		// s contains some number of things
		n := (len(s) - 1) / thingLen
		if n > numThings {
			n = numThings
		}
		thisFrame := mol.Coords[frame]
		for i := 0; i < n*thingLen && ctr < numCoords; i += thingLen {
			thisFrame[ctr] = what.Atof32(s[i : i+thingLen])
			ctr++
		}
	}

	// TODO: if it has a box, read the box, and move on to the next frame
	if mol.GetInt("POINTERS", IFBOX) > 0 {
		fmt.Fprintf(os.Stderr, "There's a box! We should do something about this.\n")
	}
}

// Returns indices of space-delimited tokens in a string. Does not tolerate
// well strings with leading spaces, so use TrimString first.
// Useful in parsing, because strings.Split() does not like arbitrarily repeated whitespace.
// Also, this assumes that start and end are big enough (e.g. you already know
// how many tokens there are in the string).
func tokenIndices(s string, start, end []int) {
	sawWhitespace := true
	tokenCounter := 0

	for i := 0; i < len(s); i++ {
		// Increment the token count if this is not a whitespace character
		// but the previous character was
		if s[i] != ' ' && sawWhitespace {
			sawWhitespace = false
			start[tokenCounter] = i
		}
		if !sawWhitespace && (s[i] == ' ' || i == len(s)-1) {
			sawWhitespace = true
			end[tokenCounter] = i
			//fmt.Println(tokenCounter, start[tokenCounter], end[tokenCounter], s[start[tokenCounter]:end[tokenCounter]]);
			tokenCounter++
		}

	}
	return
}

// Returns, I think, a slice of a string with leading and
// trailing whitespace removed
func trimSpace(s string) string {
	var i, j int
	for i = 0; j < len(s) && s[i] == ' '; i++ {
	}
	for j = len(s) - 1; j >= 0 && s[j] == ' ' || s[j] == '\n'; j-- {
	}
	return s[i : j+1]
}

// Used to emulate ReadString. This allows us to load the entire
// file at once and read pieces of it easily. Questionable how much
// more efficient it is, and whether the additional complexity is
// worth it...
type fakeStream struct {
	data []byte
	ptr  int
}

// Acts more or less like ReadString('\n') from the Go library
func (s *fakeStream) readString() string {
	a := s.ptr
	// Advance to newline
	for ; s.ptr < len(s.data) && s.data[s.ptr] != '\n'; s.ptr++ {
	}
	b := s.ptr + 1
	// Advance past newline, or signal EOF
	if s.ptr < len(s.data)-1 {
		s.ptr++
	} else {
		s.ptr = -1
	}
	if a < 0 || b > len(s.data) {
		fmt.Println(a, b, "whee!")
	}
	return string(s.data[a:b])
}

// Loads the entire contents of a file into memory, silently
// ungzipping it if the filename ends in .gz
func inhaleFile(filename string) fakeStream {
	var fileData []byte
	var err os.Error
	// Verbose!
	if strings.HasSuffix(filename, ".gz") {
		fd, err := os.Open(filename, os.O_RDONLY, 0)
		if err != nil {
			fmt.Println("Error opening", filename, err)
			return fakeStream{nil, -1}
		}
		defer fd.Close()
		file, err := gzip.NewInflater(fd)
		fileData, err = ioutil.ReadAll(file)
	} else {
		fileData, err = ioutil.ReadFile(filename)
		if err != nil {
			fmt.Println("Error opening", filename, err)
			return fakeStream{nil, -1}
		}
	}
	return fakeStream{fileData, 0}
}

// Loads pairwise residue interaction energies from an mdout file.
// The mdout must have been generated by my hacked sander, which
// makes mdouts of the form:
// [stuff]
// TOM Start of energies
// n00 n01 n02 n03
// n10 n11 n12 n13
// n20 n21 n22 n23
// n30 n31 n32 n33
// TOM End of energies
// [stuff]
// If the filenamed ends in .gz, it is silently ungzipped during loading.
func LoadEnergiesFromMdout(filename string) ([]float32, int) {
	mdout := inhaleFile(filename)
	var s string
	var numResidues int
	var data []float32
	var start, end []int
	row := 0
	isData := false
	for mdout.ptr != -1 {
		s = mdout.readString()
		// Grab the number of residues from the mdout
		if len(s) >= 60 && isData == false && s[51:59] == "NRES   =" {
			numResidues, _ = strconv.Atoi(trimSpace(s[59:67]))
			data = make([]float32, numResidues*numResidues)
			start = make([]int, numResidues)
			end = make([]int, numResidues)
		}
		// Sentinel for start and end of residue-residue interaction energies
		if len(s) >= 3 && s[0:3] == "TOM" {
			isData = !isData
		} else if isData {
			// This is energy data
			s = trimSpace(s)
			tokenIndices(s, start, end)
			offset := row * numResidues
			for i := 0; i < numResidues; i++ {
				data[offset] = float32(Strtod(s[start[i]:end[i]]))
				offset++
			}
			row++
		}
	}
	return data, numResidues
}

// Prints to stderr what percentage of your task is done.
// Call it every time you do another piece of the total.
// It will only actually print something every few percent.
func PrintPercentDone(done, total int, description string) {
	percent := int(float(done) / float(total) * 100)
	if done%(total/20) == 0 || done == total {
		desc := " " + description
		fmt.Fprintf(os.Stderr, "%d%% done%s.\n", percent, desc)
	}
}

// Only handles space-delimited data with no extra whitespace anywhere.
// Returns data, numrows, numcols
func LoadTextAsFloat32Matrix(filename string) ([]float32, int, int) {
	contents := inhaleFile(filename)
	//nukeExtraSpaces := regexp.MustCompile("[ ]+");
	numCols := -1
	var rows vector.Vector

	for contents.ptr != -1 {
		s := trimSpace(contents.readString())
		// Go's regexp is incredibly slow
		//s = nukeExtraSpaces.ReplaceAllString(s, " ");
		thisRow := strings.Split(s, " ", 0)
		// The first non-empty row defines the number of columns for the other rows
		if numCols == -1 {
			numCols = len(thisRow)
		}
		// fmt.Fprintf(os.Stderr, "%d: %d columns in this row\n", rows.Len(), len(thisRow));
		// Except for empty lines, give up if the number of columns varies
		if len(thisRow) != numCols {
			fmt.Fprintf(os.Stderr, "Loaded %d rows and %d columns. %s\n", rows.Len(), numCols, Status())
			break
		}
		rows.Push(thisRow)
	}
	// Now we have a vector of []string
	numRows := rows.Len()
	data := make([]float32, numRows*numCols)
	for i := 0; i < rows.Len(); i++ {
		// fmt.Fprintf(os.Stderr, "row %d\n", i);
		vals := rows.At(i).([]string)
		for j := 0; j < numCols; j++ {
			data[i*numCols+j] = float32(Strtod(vals[j]))
		}
	}
	fmt.Fprintf(os.Stderr, "Finally! Done converting strings to floats. %s\n", Status())
	return data, numRows, numCols
}

// Dumps a float32 matrix (stored as a 1D array, by rows) to a text file. The number of rows
// is inferred but you do have to say how many columns it has.
func DumpFloat32MatrixAsText(data []float32, cols int, filename string) {
	fp, err := os.Open(filename, os.O_WRONLY|os.O_CREAT|os.O_TRUNC, 0644)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening %f\n", filename)
		return
	}
	defer fp.Close()
	for i := 0; i < len(data); i++ {
		// Newline after every row
		if i%cols == 0 && i > 0 {
			fmt.Fprintf(fp, "\n")
		}
		fmt.Fprintf(fp, "%f ", data[i])
	}
	fmt.Fprintf(fp, "\n")
}

// Dumps a float64 matrix (stored as a 1D array, by rows) to a text file. The number of rows
// is inferred but you do have to say how many columns it has.
func DumpFloat64MatrixAsText(data []float64, cols int, filename string) {
	fp, err := os.Open(filename, os.O_WRONLY|os.O_CREAT|os.O_TRUNC, 0644)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening %f\n", filename)
		return
	}
	defer fp.Close()
	for i := 0; i < len(data); i++ {
		// Newline after every row
		if i%cols == 0 && i > 0 {
			fmt.Fprintf(fp, "\n")
		}
		fmt.Fprintf(fp, "%f ", data[i])
	}
	fmt.Fprintf(fp, "\n")
}

// Dumps a Vector of Pairs to a text file, one pair per line.
func DumpPairVectorAsText(v *vector.Vector, filename string) {
	fp, err := os.Open(filename, os.O_WRONLY|os.O_CREAT|os.O_TRUNC, 0644)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening %f\n", filename)
		return
	}
	defer fp.Close()
	for i := 0; i < v.Len(); i++ {
		p := v.At(i).(Pair)
		fmt.Fprintf(fp, "%d %d\n", p.Row, p.Col)
	}
}

// Dumps an array of coordinates in rst format.
func DumpCoordsAsRst(coords []float32) {
	fmt.Printf("We are the champions\n%d\n", len(coords)/3)
	for i := 0; i < len(coords); i++ {
		if i > 0 && i%6 == 0 {
			fmt.Printf("\n")
		}
		fmt.Printf("%12.6f", coords[i])
	}
}

// Reads bytes from an open file up to the next newline.
func readLineFromOpenFile(fp *os.File) string {
	buf := make([]byte, 256)
	var i int
	fp.Read(buf[0:1])
	for i = 1; i < len(buf) && buf[i-1] != '\n'; i++ {
		fp.Read(buf[i : i+1])
	}
	return string(buf[0:i])
}

// Assumes file pointer is at the start of a frame
func GetNextFrameFromTrajectory(trj *bufio.Reader, numAtoms int, hasBox bool) []float32 {
	// Calculate how many lines per frame
	linesPerFrame := numAtoms * 3 / 10
	if numAtoms*3%10 != 0 {
		linesPerFrame++
	}
	// fmt.Println("Lines per frame:", linesPerFrame);

	coords := make([]float32, numAtoms*3)
	ci := 0

	for i := 0; i < linesPerFrame; i++ {
		// Read and trim a line
		line, err := trj.ReadString('\n')
		if err != nil {
			break
		}
		line = trimSpace(line)
		// Each token is ended by whitespace or eol.
		for b := 0; b < len(line); {
			var e int
			// Advance to end of token
			for e = b; e < len(line) && line[e] != ' ' && line[e] != '\n' && line[e] != '\t'; e++ {
			}
			//fmt.Println(b, e, line[b:e]);
			coords[ci] = what.Atof32(line[b:e])
			ci++
			// Now, e is either at end of line or points to whitespace.
			// Advance b to start of next token or eol
			for b = e; b < len(line) && (line[b] == ' ' || line[b] == '\n' || line[b] == '\t'); b++ {
			}
		}

	}

	// Eat box line
	if hasBox {
		trj.ReadString('\n')
	}
	return coords
}

// Prints something about the current state of the program.
// Right now just tells you how much memory is allocated, which can be different
// from the amount of memory *you* allocated. Likely because the Go memory manager
// allocates from its own pool of memory which it grows and shrinks speculatively.
func Status() string {
	return fmt.Sprintf("Allocated memory: %.1f MB", float(malloc.GetStats().Alloc)/1048576)
}

// Encapsulates the indices of a residue interaction pair
type Pair struct {
	Row, Col int
}

func (p *Pair) String() string { return fmt.Sprintf("(%d, %d)", p.Row, p.Col) }

// Absolute value of an int
func Abs(n int) int {
	if n < 0 {
		return -n
	}
	return n
}

// Absolute value of a float32
func Fabs(n float32) float32 {
	if n < 0 {
		return -n
	}
	return n
}
