// SIMD stuff for AMD64

// func Foo(a, b *float32) c *float32
TEXT main·Foo(SB),7,$0
	MOVQ a+0(FP), AX
	MOVQ b+8(FP), BX
	MOVAPS 0(AX), X0
	MOVAPS 0(BX), X1
	MULPS X1, X0
	MOVAPS X0, c+16(FP)
	RET	
