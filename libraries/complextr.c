/******************  complextr.c  (in su2.a) ****************************
*									*
* complex complextrace_su2( su2_matrix *a,*b)				*
* return Tr( A_adjoint*B )   						*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"
#include <stdio.h>

void complextrace_su2( su2_matrix *a, su2_matrix *b, complex *tr ) {
register int i,j;
register Real sumr, sumi;

    for(sumr=0.0,sumi=0.0,i=0;i<2;i++)for(j=0;j<2;j++){
     sumr+= a->e[i][j].real*b->e[i][j].real + a->e[i][j].imag*b->e[i][j].imag;
     sumi+= a->e[i][j].real*b->e[i][j].imag - a->e[i][j].imag*b->e[i][j].real;
    }
    (*tr).real= sumr; (*tr).imag=sumi; 
    
}
