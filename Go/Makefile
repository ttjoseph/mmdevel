all : coulomb correl crdpipe

crdpipe : crdpipe.go
	6g crdpipe.go
	6l -o crdpipe crdpipe.6

correl : correl.go sqrt32_decl.go sqrt32_amd64.s
	6g correl.go sqrt32_decl.go
	6a sqrt32_amd64.s
	6l -o correl correl.6 sqrt32_amd64.6

coulomb : coulomb.go sqrt32_decl.go sqrt32_amd64.s
	6g coulomb.go sqrt32_decl.go
	6a sqrt32_amd64.s
	6l -o coulomb coulomb.6 sqrt32_amd64.6
	
clean :
	rm -f *.6 coulomb correl
