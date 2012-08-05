// Old Go code
package oldcode

func main_oldCorrel() {
	const numFrames = 7500
	const numResidues = 709
	fp, _ := os.Open("energies-all.bin", os.O_RDONLY, 0)
	defer fp.Close()
	fmt.Println("Loading frames...")
	energies := make([][]float32, numFrames)
	for i := 0; i < numFrames; i++ {
		energies[i] = make([]float32, numResidues*numResidues)
		loadSingleFrame(fp, numResidues, energies[i])
		PrintPercentDone(i+1, numFrames, "loading frames")
	}
	fmt.Println("Done loading. Number of residues:", numResidues)
	fmt.Println(status())
	average := AverageEnergies(energies, numResidues)
	pairs := PairsAboveCutoff(average, numResidues, 20)
	fmt.Println("Found", pairs.Len(), "pairs above cutoff. Dumping to file...")
	DumpPairVectorAsText(pairs, "pairs.txt")

	pairsEnergies := MakePairsEnergies(energies, pairs, numResidues, numFrames)

	fmt.Println("Calculating correlations...")
	correl := CalcCorrelationsMonolithic(pairsEnergies, average, pairs, numResidues, numFrames)
	fmt.Println("Dumping correlation matrix to file...")
	DumpFloat32MatrixAsText(correl, pairs.Len(), "correl.txt")
	fmt.Println(status())
}

// Returns a numPairs*numPairs correlation matrix
func CalcCorrelationsMonolithic(pairsEnergies, average []float32, pairs *vector.Vector, numResidues, numFrames int) []float32 {
	correl := make([]float32, pairs.Len()*pairs.Len())
	ch := make(chan int)
	for ij := 0; ij < pairs.Len(); ij++ {
		go calcCorrelationsSingleMonolithic(pairsEnergies, average, pairs, numResidues, numFrames,
			ij, correl, ch)
	}
	// Wait for goroutines
	for i := 0; i < pairs.Len(); i++ {
		<-ch
	}
	return correl
}

// Goroutine to do part of the correlation calculation
func calcCorrelationsSingleMonolithic(pairsEnergies, average []float32, pairs *vector.Vector, numResidues, numFrames, ij int, correl []float32, ch chan int) {
	for kl := 0; kl < pairs.Len(); kl++ {
		if ij == kl {
			continue
		}
		var num, denom float32 = 0, 0
		pij := pairs.At(ij).(Pair)
		pkl := pairs.At(kl).(Pair)
		for t := 0; t < numFrames; t++ {
			a := pairsEnergies[numFrames*ij+t] - average[numResidues*pij.row+pij.col]
			b := pairsEnergies[numFrames*kl+t] - average[numResidues*pkl.row+pkl.col]
			num += a * b
			denom += float32(math.Sqrt(float64(a * a * b * b)))
		}
		correl[pairs.Len()*ij+kl] = num / denom
	}
	ch <- ij // We're done
}

func AverageEnergies(energies [][]float32, numResidues int) []float32 {
	result := make([]float32, numResidues*numResidues)
	ch := make(chan int, numResidues)
	for row := 0; row < numResidues; row++ {
		go averageEnergiesSingleRow(energies, numResidues, row, result, ch)
	}
	// Wait for goroutines to finish
	for row := 0; row < numResidues; row++ {
		<-ch
	}
	return result
}

// Goroutine to calculate average energies across a single row
func averageEnergiesSingleRow(energies [][]float32, numResidues, row int, result []float32, ch chan int) {
	numFrames := len(energies)
	offs := numResidues * row
	// Add up this row in all the frames
	for frame := 0; frame < numFrames; frame++ {
		for i := 0; i < numResidues; i++ {
			result[offs+i] += energies[frame][offs+i]
		}
	}
	// Divide by number of frames to get the average!
	for i := 0; i < numResidues; i++ {
		result[offs+i] /= float32(numFrames)
	}
	ch <- 1 // Signal we're done
}
