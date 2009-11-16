package simd

type Float32Q uintptr;

func Foo(a, b Float32Q) Float32Q;

// Vector addition. It's OK if a and result are the same.
func VectorAdd32(a, b, result []float32) {
	
}

func VectorSqrt32(a, b, result []float32) {

}
