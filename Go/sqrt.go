// Go's math.Sqrt is slow as balls. Here's a faster version.
// Tom Joseph <ttjoseph@gmail.com>
// Algorithm taken from http://www.azillionmonkeys.com/qed/sqroot.html
package main

import "fmt"
import "math"

/*
    double fsqrt (double y) {
    double x, z, tempf;
    unsigned long *tfptr = ((unsigned long *)&tempf) + 1;

	tempf = y;
	*tfptr = (0xbfcdd90a - *tfptr)>>1; // estimate of 1/sqrt(y)
	x =  tempf;
	z =  y*0.5;                        // hoist out the “/2”
	x = (1.5*x) - (x*x)*(x*z);         // iteration formula
	x = (1.5*x) – (x*x)*(x*z);
	x = (1.5*x) – (x*x)*(x*z);
	x = (1.5*x) – (x*x)*(x*z);
	x = (1.5*x) – (x*x)*(x*z);
	return x*y;
    }
*/

func Sqrt64(y float64) float64 {
	tempi := math.Float64bits(y)
	tempihi := tempi >> 32
	tempi &= 0x00000000ffffffff // clear high dword
	// fmt.Printf("tempi %x\n", tempi);
	tempi |= ((0xbfcdd90a - tempihi) >> 1) << 32
	// fmt.Printf("tempi %x\n", tempi);
	x := math.Float64frombits(tempi)
	// fmt.Printf("x %f\n", x);
	z := y * 0.5
	x = (1.5 * x) - (x*x)*(x*z) // iteration formula
	// fmt.Printf("x %f\n", x);
	x = (1.5 * x) - (x*x)*(x*z)
	// fmt.Printf("x %f\n", x);
	x = (1.5 * x) - (x*x)*(x*z)
	// fmt.Printf("x %f\n", x);
	x = (1.5 * x) - (x*x)*(x*z)
	// fmt.Printf("x %f\n", x);
	x = (1.5 * x) - (x*x)*(x*z)
	//fmt.Printf("x %f\n", x);
	return x * y
}

func main() {
	x := float64(42.4242424242424242)
	fmt.Printf("yee ha! %f %f\n", float64(x*x), Sqrt64(x*x))
	y := float32(13.1313131313)
	fmt.Printf("craaazy! %f %f\n", float32(y*y), Sqrt32(y*y))
}
