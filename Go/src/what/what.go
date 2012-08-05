// Miscellaneous stuff
package what

// #include <stdlib.h>
import "C"

func Atof64(s string) float64 { return float64(C.atof(C.CString(s))) }

func Atof32(s string) float32 { return float32(C.atof(C.CString(s))) }

func Atoi(s string) int { return int(C.atoi(C.CString(s))) }
