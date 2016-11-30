#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



void exp_su2_matrix(su2_matrix *u ,su2_matrix *v){

        su2_matrix *del = malloc(sizeof(*del));
        su2_matrix *prod = malloc(sizeof(*prod));
	double fac=1.0;
        int i=1;
        for(int j=0; j < DIMF ; j++){
        for(int k=0; k < DIMF ; k++){
        if( j != k){   
        prod->e[j][k].real = 0.0 ;
        prod->e[j][k].imag = 0.0 ;
        del->e[j][k].real = 0.0 ;
        del->e[j][k].imag = 0.0 ;
        v->e[j][k].real = 0.0;
        v->e[j][k].imag = 0.0;}
        
        else{ prod->e[j][k].real = 1.0;
        prod->e[j][k].imag = 0.0 ;
        del->e[j][k].real = 1.0 ;
        del->e[j][k].imag = 0.0 ;
        v->e[j][k].real=1.0;
        v->e[j][k].imag = 0.0 ;}
        

       }

       

}
 
        register int sum=0,counter=0;
	
	do{
        fac=fac*(double)i;
        mult_su2_nn(prod,u,prod);
        //prod=prod*u;
        scalar_mult_su2_matrix(prod, (1.0/fac),del);
        //del=prod*(1.0/fac);
        scalar_mult_add_su2_matrix(v,del,1.0,v);
        //v=v+del;
        i++;}
        while(sqrt(realtrace_su2(del,del))>GAUGETOL);
        //printf("Exponentiating the matrix\n");
        sum+=i;
        counter++;
        if(counter==100000){
        printf("mean no. of terms in exp()=%.1f\n" ,(double)sum/counter);
        counter=0;sum=0;}
}
