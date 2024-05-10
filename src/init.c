#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _coconut_em_pava(void *, void *, void *, void *, void *, void *);

//static const R_CallMethodDef CallEntries[] = {
//    {"_coconut_em_pava", (DL_FUNC) &_coconut_em_pava, 6},
//    {NULL, NULL, 0}
//};

//void R_init_coconut(DllInfo *dll)
//{
//   R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
//    R_useDynamicSymbols(dll, FALSE);
//}