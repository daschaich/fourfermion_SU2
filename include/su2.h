// -----------------------------------------------------------------
// Defines and subroutine declarations
// for the  reduced staggered-fermion system with SU(2) gauge symmetry
#ifndef _SU2_H
#define _SU2_H

#include "../include/random.h"
#include "../include/complex.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Fermions are SU(2) vectors
//
#define DIMF 2
#define GAUGETOL 0.00000000000001

typedef struct {fcomplex e[2][2];} fsu2_matrix;
typedef struct {dcomplex e[2][2];} dsu2_matrix;
typedef struct {
  fcomplex m01;
  float m00im, m11im;
  float space;
} fanti_hermitmat;
typedef struct {
  dcomplex m01;
  double m00im, m11im;
  double space;
} danti_hermitmat;

typedef struct { fcomplex c[DIMF]; } fvector;

typedef struct { dcomplex c[DIMF]; } dvector;

#if (PRECISION == 1)
#define su2_matrix     fsu2_matrix
#define anti_hermitmat fanti_hermitmat
#define vector         fvector

#else
#define su2_matrix     dsu2_matrix
#define anti_hermitmat danti_hermitmat
#define vector         dvector
#endif

#define PLUS 1          // Flags for selecting D or D_adjoint
#define MINUS -1
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Subroutine definitions
// Vector operations
// In file clearvec.c
void clearvec(vector *v);

// In file vec_copy.c
void vec_copy(vector *a, vector *b);

// In file dumpvec.c
void dumpvec(vector *v);

// In file addvec.c
//void add_vec(vector *a, vector *b, vector *c);

// In file subvec.c
void sub_vec(vector *a, vector *b, vector *c);

// In file msq_vec.c
double magsq_vec(vector *a);

// In file dot.c
 
double dot(vector *a, vector *b);


// In file s_m_vec.c
//void scalar_mult_vec(vector *src, Real scalar, vector *dest);

// In file s_m_a_vec.c
void scalar_mult_add_vec(vector *a, vector *b, Real scalar, vector *dest);
// -----------------------------------------------------------------




// Matrix operations
// In file clear_mat.c
void clear_su2mat(su2_matrix *m);



// SU(2) generators-----------------------------------------------------------------




#define NUMGEN (DIMF*DIMF - 1) 

#define POS1 1
#define POS2 2 


#define ROOT2 1.41421356237309504880168872421


void initgen(void) ;
void vectomat(complex [NUMGEN], complex [DIMF][DIMF]) ;


void computegen(void) ;

int my_gen(void);

void gen_su2();


//------------------------------------------------------------------

// ROUTINES FOR SU(2) MATRIX OPERATIONS

void mult_su2_nn(  su2_matrix *a, su2_matrix *b, su2_matrix *c  );

void mult_su2_na( su2_matrix *a, su2_matrix *b, su2_matrix *c );


void mult_su2_an( su2_matrix *a, su2_matrix *b, su2_matrix *c );

double realtrace_su2(  su2_matrix *a, su2_matrix *b );

//void trace_su2( su2_matrix *a,complex s);
//* file "trace_su2.c"
void  complextrace_su2( su2_matrix *a, su2_matrix *b,complex *tr);


void scalar_mult_su2_matrix( su2_matrix *a, double s, su2_matrix *b );

void scalar_mult_add_su2_matrix( su2_matrix *a, su2_matrix *b,double s, su2_matrix *c);

void c_scalar_mult_sub_su2mat( su2_matrix *m1, su2_matrix *m2, complex *phase, su2_matrix *m3);
// file "cs_m_s_mat.c"

void c_scalar_mult_add_su2mat( su2_matrix *m1, su2_matrix *m2, complex *phase, su2_matrix *m3);

void c_scalar_mult_su2mat( su2_matrix *b, complex *s, su2_matrix *c );

void su2_adjoint( su2_matrix *a, su2_matrix *b );
// file "su2_adjoint.c"

void make_anti_hermitian( su2_matrix *m3,  anti_hermitmat *ah3 );
// file "make_ahmat.c"

void random_anti_hermitian( anti_hermitmat *mat_antihermit, double_prn *prn_pt );
// (prn_pt passed through to myrand())
// file "rand_ahmat.c"

void uncompress_anti_hermitian( anti_hermitmat *mat_anti, su2_matrix *mat );


void su2mat_copy(su2_matrix *a,su2_matrix *b );


void trace_su2( su2_matrix *a,complex *s);


void su2_conjug( su2_matrix *a, su2_matrix *b );

void su2_trans(su2_matrix *a, su2_matrix *b);

void det_su2(su2_matrix *a, double *sum);

void exp_su2_matrix(su2_matrix *u,su2_matrix *v);

void unit_su2mat( su2_matrix *dest );


//-------------------------------------------------------------------------------------------//

// ROUTINES FOR vector OPERATIONS ( 2 COMPONENT COMPLEX )

void scalar_mult_vec( vector *a, Real s, vector *c);

void su2_projector(vector *a,vector *b, su2_matrix *c );

void mult_su2_mat_vec( su2_matrix *a, vector *b, vector *c );

void mult_su2_mat_vec_sum( su2_matrix *a, vector *b, vector *c );


void mult_su2_mat_vec_sum_4dir( su2_matrix *a, vector *b0,
 vector *b1, vector *b2, vector *b3, vector *c );


void vec_conjug(vector *v,vector *u);


//---------------------------------------------------------------------------------------------//



// -----------------------------------------------------------------
// Miscellaneous routines
// In file gaussrand.c
Real gaussian_rand_no(double_prn *prn_pt);

#include "../include/int32type.h"
void byterevn(int32type w[], int n);
void byterevn64(int32type w[], int n);

#endif
// -----------------------------------------------------------------