// -----------------------------------------------------------------
// Set up SU(2) generators in an overly complicated way
#include "su2_includes.h"

#define POS1 1
#define POS2 2

#define ROOT2 1.41421356237309504880168872421
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Adjoint to fundamental
void vectomat(complex vec[NUMGEN], complex mat[DIMF][DIMF]) {
  int i, j;
  complex mult;
  complex temp1;
  complex temp[NUMGEN];

  for (i = 0; i < DIMF; i++) {
    for (j = i + 1; j < DIMF; j++) {
      CMUL_I(vec[POS1 + posmat[i][j]], temp[POS1 + posmat[i][j]]);

      CSUB(vec[posmat[i][j]], temp[POS1 + posmat[i][j]], mat[i][j]);
      CMULREAL(mat[i][j], 0.5, mat[i][j]);

      CADD(vec[posmat[i][j]], temp[POS1 + posmat[i][j]], mat[j][i]);
      CMULREAL(mat[j][i], 0.5, mat[j][i]);
    }
  }

  for (i = 0; i < DIMF; i++)
    mat[i][i] = cmplx(0.0, 0.0);

  for (i = 0; i < DIMF - 1; i++) {
    CMULREAL(vec[POS2 + i], (1.0 / sqrt(2.0 + 2.0/(1.0 + i)) / (1.0 + i)),
             mult);

    for (j = 0; j < i + 1; j++)
      CADD(mat[j][j], mult, mat[j][j]);

    CMULREAL(mult, (1.0 + i), temp1);
    CSUB(mat[i + 1][i + 1], temp1, mat[i + 1][i + 1]);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Compute generator matrices
void computegen() {
  int a, i, j;
  complex adj[NUMGEN];
  complex fund[DIMF][DIMF];

  for (a = 0; a < NUMGEN; a++) {
    for (i=0;i<NUMGEN;i++) {
      adj[i].real =0.0;
      adj[i].imag =0.0;}

    adj[a].real = 1.0;
    adj[a].imag = 0.0;
    //adj[i] = complex(0.0,0.0);
    //adj[a] = complex(1.0,0.0);

    vectomat(adj, fund);

    for (i=0;i<DIMF;i++) {
      for (j=0;j<DIMF;j++) {
        genmat[a][i][j].real = fund[i][j].real;
        genmat[a][i][j].imag = fund[i][j].imag;
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// SU(N) routines
void initgen(void) {
  int i, j, a = 0;

  for (i = 0; i < DIMF; i++) {
    for (j = i + 1; j < DIMF; j++) {
      posmat[i][j] = a;
      a++;
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void my_gen() {
  int a, i, j;

  initgen();

  computegen();

  // Map to data structures in code
  for (a = 0; a < NUMGEN; a++) {
    for (i = 0; i < DIMF; i++) {
      for (j = 0; j < DIMF; j++) {
        CMULREAL(genmat[a][i][j], ROOT2, genmat[a][i][j]);
        //CMUL_I(genmat[a][i][j], &(Lambda[a]->e[i][j]));
        Lambda[a].e[i][j].real = -1.0 * genmat[a][i][j].imag;
        Lambda[a].e[i][j].imag = genmat[a][i][j].real;
        //Lambda[a].set(i,j,(complex(0.0,1.0)*sqrt(2.0))*genmat[a][i][j]);
      }
    }
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void gen_su2() {
  int i, j, k;

  for (i = 0; i < NUMGEN; i++)
    clear_su2mat(&(Lambda[i]));

  my_gen();
  printf("Computing SU(2) generators\n");

  for (i = 0; i < NUMGEN; i++) {
    printf("Lambda[%d]:\n", i);
    for (j = 0; j < DIMF; j++) {
      for (k = 0; k < DIMF; k++) {
        node0_printf("Real = %.1f Imag = %.1f\n",
                     Lambda[i].e[j][k].real, Lambda[i].e[j][k].imag);

      }
    }
  }
}
// -----------------------------------------------------------------
