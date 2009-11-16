include $(GOROOT)/src/Make.$(GOARCH)

TARG=simd

OFILES_amd64=\
	simd_amd64.$O\

OFILES=\
	$(OFILES_$(GOARCH))

GOFILES=\
	simd_decl.go

CGO_CFLAGS=
CGO_LDFLAGS=

include $(GOROOT)/src/Make.pkg
