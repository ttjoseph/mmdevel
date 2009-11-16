
#include "runtime.h"
#include "cgocall.h"

#pragma dynld initcgo initcgo "/home/josept02/build/go/pkg/linux_amd64/libcgo.so"
#pragma dynld libcgo_thread_start libcgo_thread_start "/home/josept02/build/go/pkg/linux_amd64/libcgo.so"
#pragma dynld _cgo_malloc _cgo_malloc "/home/josept02/build/go/pkg/linux_amd64/libcgo.so"
#pragma dynld _cgo_free free "/home/josept02/build/go/pkg/linux_amd64/libcgo.so"

void
what·_C_GoString(int8 *p, String s)
{
	s = gostring((byte*)p);
	FLUSH(&s);
}

void
what·_C_CString(String s, int8 *p)
{
	p = cmalloc(s.len+1);
	mcpy((byte*)p, s.str, s.len);
	p[s.len] = 0;
	FLUSH(&p);
}

#pragma dynld _cgo_atof _cgo_atof "/home/josept02/build/go/pkg/linux_amd64/what_what.so"
void (*_cgo_atof)(void*);

void
what·_C_atof(struct{uint8 x[16];}p)
{
	cgocall(_cgo_atof, &p);
}

