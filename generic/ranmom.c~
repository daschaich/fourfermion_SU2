/*************************** ranmom.c *******************************/
/* Produce Gaussian random momenta for the gauge fields. */

#include "generic_includes.h"
//#include <defines.h>                 /* For SITERAND */

void ranmom(){
register int i,dir;
register site *s;
    FORALLSITES(i,s){
	for(dir=XUP;dir<=TUP;dir++){

	    if(dir==TUP || s->t>0){
	random_anti_hermitian( (anti_hermitmat *)&(s->mom[dir]),
		    &(s->site_prn) );
	random_anti_hermitian( (anti_hermitmat *)&(s->mom[dir]),
		    &node_prn );


/*  It would not matter if we did not set the spatial momenta at t=0
    to zero, because we never attempt to update the t=0 spatial links
*/

	    }
	    else{
		s->mom[dir].m00im = 0.0;
		s->mom[dir].m11im = 0.0;
		s->mom[dir].m01.real = 0.0;
		s->mom[dir].m01.imag = 0.0;

	    }
            
	}
    }

}


