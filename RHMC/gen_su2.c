//#include "my_gen.h"
#include "../include/complex.h"
#include "lattice.h"
#include "math.h"
#include "su2_includes.h"


//const int POS1 = (DIMF*(DIMF-1)/2);
//const int POS2 = (DIMF*(DIMF-1));



void gen_su2(){

for(int i=0 ;i < NUMGEN ; i++){

clear_su2mat(&(Lambda[i]));

}

int gen;

gen = my_gen();



if(gen == 1){printf("Computing SU(2) generators\n"); 

for(int i=0 ;i < NUMGEN ; i++){
for(int j=0 ; j< 2 ; j++){
for(int k=0 ; k < 2 ; k++){

printf("Real = %.1f Imag = %.1f\n",Lambda[i].e[j][k].real, Lambda[i].e[j][k].imag);

        }
      }
printf("------------------------------------\n");
    }

  }


}


int my_gen(void)
{
  int aa, ii, jj ;
  
  
  initgen() ;

  computegen() ;

  /* map to data structures in code */
  
  for(aa=0;aa<NUMGEN;aa++) {

    for(ii=0;ii<DIMF;ii++) {
      for(jj=0;jj<DIMF;jj++) {
	
      
      CMULREAL(genmat[aa][ii][jj], ROOT2 ,genmat[aa][ii][jj]);
      //CMUL_I(genmat[aa][ii][jj], &(Lambda[aa]->e[ii][jj]));
      Lambda[aa].e[ii][jj].real = -1.0 * genmat[aa][ii][jj].imag;
      Lambda[aa].e[ii][jj].imag = genmat[aa][ii][jj].real; 
      //Lambda[aa].set(ii,jj,(complex(0.0,1.0)*sqrt(2.0))*genmat[aa][ii][jj]);
      } 
    }

  }


  return(1) ;

}



/* Compute generator matrices */

void computegen(void) 
{
  int aa, ii, jj ;
complex adj[NUMGEN] ;
complex fund[DIMF][DIMF] ; 

  for(aa=0;aa<NUMGEN;aa++) {

    for(ii=0;ii<NUMGEN;ii++) {
    adj[ii].real =0.0;
    adj[ii].imag =0.0;}

    adj[aa].real = 1.0;
    adj[aa].imag = 0.0;
    //adj[ii] = complex(0.0,0.0) ;
    //adj[aa] = complex(1.0,0.0) ;

    vectomat(adj,fund) ;

    for(ii=0;ii<DIMF;ii++) {
      for(jj=0;jj<DIMF;jj++) {
	genmat[aa][ii][jj].real = fund[ii][jj].real ;
        genmat[aa][ii][jj].imag = fund[ii][jj].imag ;
      } 
    }

  }

}



/* SU(N) routines */

void initgen(void)
{
  int ii,jj,aa ;

  aa=0 ;

  for(ii=0;ii<DIMF;ii++) {
    for(jj=ii+1;jj<DIMF;jj++) {

      posmat[ii][jj] = aa ;

      aa+=1 ;

    }
  }

}


/* Adjoint to fundamental */

void vectomat(complex vec[NUMGEN], complex mat[DIMF][DIMF])
{
int ii,jj ;
complex mult ;
complex temp1;
complex temp[NUMGEN];

  for(ii=0;ii<DIMF;ii++) {
    for(jj=ii+1;jj<DIMF;jj++) {
     
     CMUL_I(vec[POS1+ posmat[ii][jj]],temp[POS1 + posmat[ii][jj]]);
     

     CSUB(vec[posmat[ii][jj]],temp[POS1 + posmat[ii][jj]],mat[ii][jj]);
     CMULREAL(mat[ii][jj], 0.5 , mat[ii][jj]); 
     //mat[ii][jj] = 0.5*(vec[posmat[ii][jj]]-complex(0.0,1.0)*vec[POS1+posmat[ii][jj]]) ;
      

     
     CADD(vec[posmat[ii][jj]],temp[POS1 + posmat[ii][jj]],mat[jj][ii]);
     CMULREAL(mat[jj][ii], 0.5 , mat[jj][ii]); 


      //mat[jj][ii] = 0.5*(vec[posmat[ii][jj]]+complex(0.0,1.0)*vec[POS1+posmat[ii][jj]]) ;
    }
  }
    
  for(ii=0;ii<DIMF;ii++) 
   mat[ii][ii].real = 0.0;
   mat[ii][ii].imag = 0.0; 
   //mat[ii][ii] = complex(0.0,0.0) ; 
    
  for(ii=0;ii<DIMF-1;ii++) {
    CMULREAL(vec[POS2 + ii] , (1.0/sqrt(2.0 + 2.0/(1.0 + ii)) /(1.0 + ii)) , mult);
    //mult = vec[POS2+ii]*(1.0/sqrt(2.+2./(1.+ii))/(1.+ii)) ;
    for(jj=0;jj<ii+1;jj++) {    
      
     //mat[jj][jj] = mat[jj][jj]+mult ;

     CADD(mat[jj][jj] ,mult ,mat[jj][jj]);
    }
     CMULREAL(mult,(1.0 + ii) ,temp1)
     CSUB(mat[ii+1][ii+1] , temp1 ,mat[ii+1][ii+1]);
     //mat[ii+1][ii+1] = mat[ii+1][ii+1]- 1.0*complex(1.0+ii,0.0)*mult ; 
  }
    
}


