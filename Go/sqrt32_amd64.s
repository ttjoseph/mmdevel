// func Sqrt32(x float32) float32
TEXT main·Sqrt32(SB),7,$0
	MOVSS x+0(FP), X0
    SQRTSS X0, X0
	MOVSD X0, r+8(FP)
	RET

// func Invsqrt32(x float32) float32
TEXT main·Invsqrt32(SB),7,$0
	MOVSS x+0(FP), X0
	RSQRTSS X0, X0 // Inverse square root
	MOVSD X0, r+8(FP)
	RET
