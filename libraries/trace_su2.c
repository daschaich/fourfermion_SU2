/*******************  trace_su2.c  (in su2.a) ***************************
*									*
* complex trace_su2(a) su2_matrix *a;					*
* return complex trace of an su2 matrix 				*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

/* Complex trace of an su2 matrix */
void trace_su2( su2_matrix *a,complex *s) {
    complex tr;
    CADD(a->e[0][0],a->e[1][1],tr);
    s =&tr;
    //CADD(t1,a->e[2][2],t2);
    
}
