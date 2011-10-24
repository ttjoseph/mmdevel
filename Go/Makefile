all : coulomb crdpipe

crdpipe : crdpipe.go
	6g crdpipe.go
	6l -o crdpipe crdpipe.6

coulomb : coulomb.go sqrt32_decl.go sqrt32_amd64.s
	6g coulomb.go sqrt32_decl.go
	6a sqrt32_amd64.s
	gopack grc coulomb.a coulomb.6 sqrt32_amd64.6
	6l -o coulomb coulomb.a
	
clean :
	rm -f *.6 *.a coulomb correl
