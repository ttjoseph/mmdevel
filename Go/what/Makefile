include $(GOROOT)/src/Make.inc

TARG=what

CGOFILES=\
	what.go

include $(GOROOT)/src/Make.pkg

%: install %.go
	$(GC) $*.go
	$(LD) -o $@ $*.$O

