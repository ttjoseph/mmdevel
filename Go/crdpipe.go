// Extracts an arbitrary snapshot from an AMBER trajectory file.
// Ideally, you can put it on one end of a named pipe and have sander read it.
// $ mkfifo foo
// $ crdpipe -mdcrd foo.trj -frame 42 -numatoms 19401 > foo &
// $ $AMBERHOME/bin/sander -p foo.top -c foo ...
package main
import ( "fmt"; "os"; "flag"; "math"; )
import ( "amber"; )

func main() {
    var trjFilename string;
    var frame, numAtoms int;
    var hasBox bool;
    flag.StringVar(&trjFilename, "mdcrd", "mdcrd", "Trajectory filename");
    flag.IntVar(&frame, "frame", 0, "Frame number you want (counting starts at 0)");
    flag.IntVar(&numAtoms, "numatoms", 0, "Number of atoms per frame in the trajectory file");
    flag.BoolVar(&hasBox, "hasbox", true, "Whether the trajectory records periodic boxes");
    flag.Parse();
    
    GetFrameFromTrajectory(trjFilename, frame, numAtoms, hasBox);
}

// Read bytes
func readLine(fp *os.File) string {
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
    
    //
    header := readLine(fp);
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
    fmt.Printf(header);
    fmt.Printf(readLine(fp));
    fmt.Printf(readLine(fp));
    fmt.Printf(readLine(fp));
    fmt.Printf(readLine(fp));
    amber.Strtod("1.2");
    return nil;
}