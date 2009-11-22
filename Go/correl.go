package main
import ( "flag"; "fmt"; "container/vector"; "os"; "encoding/binary"; "math"; "strings"; "amber"; )

// Number of frames to process at once
const BATCH_SIZE = 200;
const NUM_RESIDUES = 709;

func main() {
    var command, filename string;
    flag.StringVar(&command, "cmd", "", "dump-mdouts,project-matrix,calc-correl");
    flag.StringVar(&filename, "f", "", "Input or output filename, depending on cmd");
    flag.Parse();
    
    switch(command) {
        case "dump-mdouts": main_dumpMdoutsToBin(filename);
        case "project-matrix": main_projectMatrix();
        case "calc-correl": main_calcCorrel(filename);
        default: flag.Usage();
    }
}

func genericErrorHandler(err os.Error) {
    if err != nil {
        fmt.Fprintf(os.Stderr, "%s\n", err.String());
        os.Exit(1);
    }
}

// Dumps pairwise residue interaction energies to a binary file.
// Doing this drastically speeds up reading it back later (Atof32 is slow as balls, looks like).
// This way we only have to deal with the mdout files once.
func main_dumpMdoutsToBin(outFilename string) {
    const dirname = "snapshots/";
    dir, err := os.Open(dirname, os.O_RDONLY, 0);
    genericErrorHandler(err);
    names, _ := dir.Readdirnames(100000);
    filenames := new(vector.Vector);
    for _, v := range(names) {
        if strings.HasSuffix(v, ".mdout.gz") { filenames.Push(dirname + v) }
    }
    // Output file
    outFile, err := os.Open(outFilename, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0644);
    genericErrorHandler(err);
    defer outFile.Close();
    // Load each mdout from gzip
    numFrames := filenames.Len();
    var numResidues int;
    var energies []float32;
    fmt.Println("Loading", numFrames, "mdout files.", amber.Status());
    for i := 0; i < numFrames; i++ {
       fname := filenames.At(i).(string);
       fmt.Println(i+1, "of", numFrames, fname, amber.Status());
       energies, numResidues = amber.LoadEnergiesFromMdout(fname);
       // Wasteful to allocate a new tmp every frame, so let's do it anyway and 
       // trust the GC to not be stupid (this way we don't have to know numResidues a priori)
       tmp := make([]byte, numResidues*numResidues*4);
       // Dump to file. We have to explicitly convert to bytes. Yay.
       for j, n := range(energies) {
           binary.BigEndian.PutUint32(tmp[j*4:j*4+4], math.Float32bits(n));
       }
       outFile.Write(tmp);
    }
    fmt.Println("Number of residues:", numResidues);
    fmt.Println(amber.Status());
}

// Project pairs interaction matrix to residue interaction matrix
func main_projectMatrix() {
    fmt.Fprintf(os.Stderr, "Projecting pairs correlation matrix to yield residue correlatin matrix.\n");
    fmt.Fprintf(os.Stderr, "Using correl.txt and pairs.txt.\n");
    correl, numPairs, _ := amber.LoadTextAsFloat32Matrix("correl.txt");
    pairsMat, numPairs2, _ := amber.LoadTextAsFloat32Matrix("pairs.txt");
    if numPairs != numPairs2 {
        fmt.Fprintf(os.Stderr, "%d and %d numPairs don't agree!", numPairs, numPairs2);
        return;
    }

    // Convert matrix to vector of Pairs
    pairs := new(vector.Vector);
    for i := 0; i < len(pairsMat); i+=2 {
        pairs.Push(amber.Pair{int(pairsMat[i]), int(pairsMat[i+1])});
    }
    
    fmt.Fprintf(os.Stderr, "Projecting matrix...\n");
    rescorrel := CorrelToResidueInteractionMatrix(correl, pairs, NUM_RESIDUES);
    amber.DumpFloat32MatrixAsText(rescorrel, NUM_RESIDUES, "rescorrel.txt");
}

// Calculates correlation in interaction energies over time
func main_calcCorrel(filename string) {
  // const filename = "energies2.bin"; // Raw binary format
  
  fmt.Println("Calculating average interaction energy matrix.");
  average := AverageEnergies(filename, NUM_RESIDUES);
  amber.DumpFloat32MatrixAsText(average, NUM_RESIDUES, "average.txt");
  //average, _, _ := amber.LoadTextAsFloat32Matrix("average.txt");
  pairs := PairsAboveCutoff(average, NUM_RESIDUES, 15);
  fmt.Println("Found", pairs.Len(), "pairs above cutoff. Dumping to file...");
  amber.DumpPairVectorAsText(pairs, "pairs.txt");
  correl := CalcCorrelations(filename, average, pairs, NUM_RESIDUES);
  correl = correl;
  fmt.Println("Done calculating correlation matrix. Dumping to file...");
  amber.DumpFloat32MatrixAsText(correl, pairs.Len(), "correl.txt");
  fmt.Println(amber.Status());
}

// "Projects" pair-pair correlation matrix to produce a residue-residue interaction energy correlation
// matrix. The resulting matrix is numResidues*numResidues.
func CorrelToResidueInteractionMatrix(correl []float32, pairs *vector.Vector, numResidues int) []float32 {
    ch := make(chan int);
    result := make([]float32, numResidues*numResidues);
    numGoroutines := 0;
    for i := 0; i < numResidues; i++ {
        for j := 0; j < i; j++ {
            // if i == j { continue }
            go projectCorrelSingleResiduePair(correl, pairs, numResidues, i, j, result, ch);
            numGoroutines++;
        }
        // if numGoroutines > 3000 { break } // DEBUG
    }
    // Wait for goroutines to finish
    for i := 0; i < numGoroutines; i++ { <-ch }
    return result;
}

// Goroutine that iterates over correlation matrix as in eqn 4 in Kong and Karplus 2007
func projectCorrelSingleResiduePair(correl []float32, pairs *vector.Vector, numResidues, i, j int, result []float32, ch chan int) {
    // fmt.Fprintf(os.Stderr, "Bonjour %d %d\n", i, j);
    numPairs := pairs.Len();
    var totalCorr float64;
    for ij := 0; ij < numPairs; ij++ {
        pij := pairs.At(ij).(amber.Pair);
        for kl := 0; kl < ij; kl++ {
            // if ij == kl { continue }
            pkl := pairs.At(kl).(amber.Pair);
            if ((pij.Row == i || pij.Col == i) && (pkl.Row == j || pkl.Col == j))
                || ((pij.Row == j || pij.Col == j) && (pkl.Row == i || pkl.Col == i)) {
                totalCorr += float64(amber.Fabs(correl[numPairs*ij+kl]))
            }
        }
    }
    result[numResidues*i+j] = float32(totalCorr);
    // Do other side of diagonal so the result matrix is symmetric
    result[numResidues*j+i] = float32(totalCorr);
    if j == 0 { fmt.Fprintf(os.Stderr, "projectCorrelSingleResiduePair done: %d. %s\n", i, amber.Status()); }
    ch <- 0;
}

// Make the raw data for correlation calculation, as follows:
// pair0 ene0 ene1 ene2 ...
// pair1 ene0 ene1 ene2 ...
// pair2 ene0 ene1 ene2 ...
//  ...
// That is, pairs.Len() rows and numFrames columns
func MakePairsEnergies(energies [][]float32, pairs *vector.Vector, numResidues, numFrames int) []float32 {
    pairsEnergies := make([]float32, pairs.Len()*numFrames);
    // Extract energies
    for i := 0; i < pairs.Len(); i++ {
        p := pairs.At(i).(amber.Pair);
        offs := numFrames*i;
        for frame := 0; frame < numFrames; frame++ {
            // Energies are stored row by row
            pairsEnergies[offs+frame] = energies[frame][numResidues*p.Row+p.Col]
        }
    }
    return pairsEnergies
}

// Determine the number of frames in the file by checking its size
// and dividing by the size of a frame
func numFramesInFile(filename string, numResidues int) int {
  stat, err := os.Stat(filename);
  if(err != nil) {
      fmt.Fprintf(os.Stderr, "Couldn't open %s: %s", filename, err);
      return -1;
  }
  // fmt.Println("Size is", stat.Size);
  return int(stat.Size / uint64(numResidues*numResidues*4));
}

// Calculates correlation matrix (numPairs*numPairs)
// Allows for processing frames in batches so we don't have to load all of them
// in memory at once.
// filename is the name of a binary file that is 
// sizeof(float32)*numPairs^2*numFrames in size
func CalcCorrelations(filename string, average []float32, pairs *vector.Vector,
        numResidues int) []float32 {
    numFrames := numFramesInFile(filename, numResidues);
    // Keep numerator, denominator matrices across calls to this function
    // Each of these matrices is numPairs*numPairs
    // We use float64 here because there could potentially be a large number of
    // frames, and we don't want too much precision error
    num := make([]float64, pairs.Len()*pairs.Len());
    denom := make([]float64, pairs.Len()*pairs.Len());
    
    // Load a batch of frames, process it, then load a new batch of frames.
    // Presumably the GC will free the old batch.
    fp, err := os.Open(filename, os.O_RDONLY, 0);
    genericErrorHandler(err);
    defer fp.Close();
    ch := make(chan int);
    for i := 0; i < numFrames; i += BATCH_SIZE {
        numFramesNow := BATCH_SIZE;
        if numFrames-i < numFramesNow { numFramesNow = numFrames - i }
        fmt.Fprintf(os.Stderr, "CalcCorrelations: frames %d-%d of %d. %s\n", i+1, i+numFramesNow, numFrames, amber.Status());
        energies := loadFrameBatch(fp, numResidues, numFramesNow);
        // fmt.Fprintf(os.Stderr, "Done loading frame batch, now processing them. %s\n", amber.Status());
        pairsEnergies := MakePairsEnergies(energies, pairs, numResidues, numFramesNow);
        const numRowsPerBatch = 200;
        numKids := 0;
        for ij := 0; ij < pairs.Len(); ij += numRowsPerBatch {
            numRowsNow := numRowsPerBatch;
            if numRowsNow > pairs.Len() - ij { numRowsNow = pairs.Len() - ij }
            go calcCorrelationsPiece(pairsEnergies, average, pairs, numResidues, numFramesNow,
                ij, ij+numRowsNow, num, denom, ch);
            numKids++;
        }
        // Wait for goroutines to finish
        for i := 0; i < numKids; i++ { <-ch }
    }
    
    // Divide numerator by denominator
    fmt.Fprintf(os.Stderr, "Done with everything. %s\n", amber.Status());
    correl := make([]float32, pairs.Len()*pairs.Len());
    for i := 0; i < len(num); i++ { correl[i] = float32(num[i] / denom[i]) }
    return correl;
}

// Goroutine to calculate some rows of the correlation matrix.
// ij = [ij_a, ij_b)
func calcCorrelationsPiece(pairsEnergies, average []float32, pairs *vector.Vector,
        numResidues, numFrames, ij_a, ij_b int, num, denom []float64, ch chan int) {
    numPairs := pairs.Len();
    for ij := ij_a; ij < ij_b; ij++ {
        for kl := 0; kl < ij; kl++ {
            // if ij == kl { continue }
            pij := pairs.At(ij).(amber.Pair);
            pkl := pairs.At(kl).(amber.Pair);
            for t := 0; t < numFrames; t++ {
                a := pairsEnergies[numFrames*ij+t] - average[numResidues*pij.Row+pij.Col];
                b := pairsEnergies[numFrames*kl+t] - average[numResidues*pkl.Row+pkl.Col];
                num[numPairs*ij+kl] += float64(a * b);
                denom[numPairs*ij+kl] += math.Sqrt(float64(a*a * b*b));
                // Do other side of diagonal so it's symmetric
                num[numPairs*kl+ij] += float64(a * b);
                denom[numPairs*kl+ij] += math.Sqrt(float64(a*a * b*b));
            }
        }
    }
    ch <- 0; // We're done
}

// Load a bunch of frames at once from a binary file of floats.
func loadFrameBatch(fp *os.File, numResidues, numFrames int) [][]float32 {
  frames := make([][]float32, BATCH_SIZE);
  for i := 0; i < BATCH_SIZE; i++ {
      frames[i] = make([]float32, numResidues*numResidues);
      loadSingleFrame(fp, numResidues, frames[i]);
  }
  return frames;
}

// Load a single frame of pairwise residue interaction energies from an open file descriptor
// that points to a binary file of floats.
func loadSingleFrame(fp *os.File, numResidues int, result []float32) {
    // A float32 will always be 4 bytes so I feel OK hardcoding it (is there a sizeof() in Go?)
    tmp := make([]byte, numResidues*numResidues*4);
    fp.Read(tmp); // Read whole frame
    for j := 0; j < numResidues*numResidues; j++ {
        result[j] = math.Float32frombits(binary.BigEndian.Uint32(tmp[j*4:j*4+4]))
    }
}

// Returns list of residue pairs whose interaction energies are above a given cutoff.
func PairsAboveCutoff(average []float32, numResidues int, cutoff float32) *vector.Vector {
    pairs := new(vector.Vector);
    for row := 0; row < numResidues; row++ {
        for col := 0; col < row-1; col++ {
            // Ignore same and adjacent residue interactions
            if amber.Abs(row-col) > 1 && amber.Fabs(average[row*numResidues+col]) > cutoff {
                // print "%d %d = %f" % (row, col, avg[row,col])
                pairs.Push(amber.Pair{row, col})
            }            
        }
    }
    return pairs;
}

// Calculates the average energy matrix by processing the file in batches.
// Turns out goroutines are not as ridiculously lightweight as I'd hoped, so
// this is single-threaded.
func AverageEnergies(filename string, numResidues int) []float32 {
  numFrames := numFramesInFile(filename, numResidues);
  fmt.Println(numFrames, "frames in file.");
  sum := make([]float64, numResidues*numResidues);
  fp, err := os.Open(filename, os.O_RDONLY, 0);
  genericErrorHandler(err);
  defer fp.Close();
  
  for frame := 0; frame < numFrames; frame += BATCH_SIZE {
    thisBatchSize := BATCH_SIZE;
    if thisBatchSize > numFrames-frame { thisBatchSize = numFrames-frame }
    fmt.Fprintf(os.Stderr, "On frames %d-%d of %d. %s\n", frame+1, frame+thisBatchSize, numFrames, amber.Status());
    energies := loadFrameBatch(fp, numResidues, thisBatchSize);
    
    // Add up all the values in this batch
    for i := 0; i < thisBatchSize; i++ {
      for j := 0; j < len(energies[i]); j++ {
        sum[j] += float64(energies[i][j])
      }
    }
  }
  
  // Divide to get average
  result := make([]float32, numResidues*numResidues);
  for i := 0; i < len(result); i++ { result[i] = float32(sum[i] / float64(numFrames)) }
  return result;
}

func min(a, b int) int {
    if a < b { return a }
    return b
}