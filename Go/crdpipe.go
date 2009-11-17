// Extracts an arbitrary snapshot from an AMBER trajectory file.
// Ideally, you can put it on one end of a named pipe and have sander read it.
// $ mkfifo foo
// $ crdpipe -mdcrd foo.trj -frame 42 -numatoms 19401 > foo &
// $ $AMBERHOME/bin/sander -p foo.top -c foo ...
package main
import ( "fmt"; "os"; "flag"; "math"; )
import ( "what"; )

func main() {
    var trjFilename string;
    var frame, numAtoms int;
    var hasBox bool;
    flag.StringVar(&trjFilename, "mdcrd", "mdcrd", "Trajectory filename");
    flag.IntVar(&frame, "frame", 0, "Frame number you want (counting starts at 0)");
    flag.IntVar(&numAtoms, "numatoms", 0, "Number of atoms per frame in the trajectory file");
    flag.BoolVar(&hasBox, "hasbox", true, "Whether the trajectory records periodic boxes");
    flag.Parse();
    
    coords := GetFrameFromTrajectory(trjFilename, frame, numAtoms, hasBox);
    DumpCoordsAsRst(coords);
}


func DumpCoordsAsRst(coords []float32) {
    fmt.Printf("We are the champions\n%d\n", len(coords)/3);
    for i := 0; i < len(coords); i++ {
        if i > 0 && i % 6 == 0 { fmt.Printf("\n") }
        fmt.Printf("%12.6f", coords[i]);
    }
}

// Read bytes
func readLineFromOpenFile(fp *os.File) string {
    buf := make([]byte, 256);
    var i int;
    fp.Read(buf[0:1]);
    for i = 1; i < len(buf) && buf[i-1] != '\n'; i++ { fp.Read(buf[i:i+1]) }
    return string(buf[0:i]);
}

func GetFrameFromTrajectory(filename string, frame, numAtoms int, hasBox bool) []float32 {
    fp, err := os.Open(filename, os.O_RDONLY, 0);
    if err != nil {
        fmt.Fprintf(os.Stderr, "Couldn't find file %s\n", filename);
        return nil;
    }
    defer fp.Close();
    
    // Eat header
    readLineFromOpenFile(fp);
    // In an mdcrd trajectory, there are 10 coordinates per line,
    // but in a single snapshot file there are 6, for some reason.
    linesPerFrameTrj := int(math.Ceil(float64(numAtoms*3) / 10));
    // There are as many newlines as lines per frame.
    // 8 bytes per coordinate, 8*3+2+1 bytes for box
    bytesPerFrame := numAtoms*8*3 + linesPerFrameTrj;
    if hasBox { bytesPerFrame += 8*3+2+1 }
    //    fmt.Fprintf(os.Stderr, "Bytes per frame: %d, lines per frame: %d; seeking to %d\n", bytesPerFrame, 
    //        linesPerFrameTrj, frame * bytesPerFrame);
    fp.Seek(int64(frame * bytesPerFrame), 1); // Seek relative to current offset
    // fmt.Printf("%s\n%d\n", header, numAtoms);
    // Read coordinates
    coords := make([]float32, numAtoms*3);
    buf := make([]byte, 8);
    for i := 0; i < numAtoms*3; i++ {
        // Eat newline every 10 coordinates
        if i > 0 && i % 10 == 0 { fp.Read(buf[0:1]) }
        fp.Read(buf);
        coords[i] = float32(what.Atof64(string(buf)));
        // if i < 150 { fmt.Printf("%f\n", coords[i]) }
    }
    return coords;
}