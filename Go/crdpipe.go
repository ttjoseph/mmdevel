// Extracts an arbitrary snapshot from an AMBER trajectory file.
// Ideally, you can put it on one end of a named pipe and have sander read it.
// $ mkfifo foo
// $ crdpipe -mdcrd foo.trj -frame 42 -numatoms 19401 > foo &
// $ $AMBERHOME/bin/sander -p foo.top -c foo ...
package main
import ( "flag"; )
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
    
    coords := amber.GetFrameFromTrajectory(trjFilename, frame, numAtoms, hasBox);
    amber.DumpCoordsAsRst(coords);
}
