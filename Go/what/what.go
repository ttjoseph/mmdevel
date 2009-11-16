package what

// #include <stdlib.h>
import "C"

func Atof64(s string) float64 {
    return float64(C.atof(C.CString(s)))
} 