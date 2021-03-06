// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
// When this routine is called the CG should already have been run,
// so that the vector **sol contains (M_adjoint*M+shift[n])^(-1) * src
#include "su2_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
Real ahmat_mag_sq(anti_hermitmat *pt);
// -----------------------------------------------------------------
double hmom_action(){
register int i ,dir;
register site *s;

double sum =0.0;
for (dir =XUP ;dir <=TUP ; dir++){
FORALLSITES(i,s)
sum +=(double) ahmat_mag_sq(&(s->mom[dir]));

}
g_doublesum(&sum);
return sum;
}

Real ahmat_mag_sq(anti_hermitmat *pt) {
  register Real x, sum;
  x = pt->m00im;      sum  = 0.5 * x * x;
  x = pt->m11im;      sum += 0.5 * x * x;
  x = pt->m01.real;   sum += x * x;
  x = pt->m01.imag;   sum += x * x;
 
  return sum;
}
// -----------------------------------------------------------------




// -----------------------------------------------------------------
// Fermion contribution to the action
// Include the ampdeg term to allow sanity check that the fermion action
// is 4*volume on average
// Since the pseudofermion src is fixed throughout the trajectory,
// ampdeg actually has no effect on Delta S (checked)
// sol, however, depends on the gauge fields through the CG
double fermion_action(vector *src, vector **sol) {
  register int i, j;
  register site *s;
  double sum = 0.0;
  double tr;
#ifdef DEBUG_CHECK
  double im = 0.0;
#endif

  FORALLSITES(i, s) {
    sum += ampdeg4 * magsq_vec(&(src[i]));
    for (j = 0; j < Norder; j++) {
      tr = dot(&(src[i]), &(sol[j][i]));   // src^dag.sol[j]
      sum += (double)(amp4[j] * tr);
    }
  }
  g_doublesum(&sum);
  
  return sum;
}
// -----------------------------------------------------------------


// Adds adjoint plaquette term
// Use tempmat for temporary storage
void plaquette_a(double *ss_plaq, double *st_plaq) {
  register int i, dir, dir2;
  register site *s;
  register su2_matrix *m1, *m4;
  
  double ss_sum = 0.0, st_sum = 0.0, td;
  complex *tc =malloc(sizeof(*tc));
  msg_tag *mtag, *mtag2;
  su2_matrix tmat;
  
FORALLSITES(i, s) {
  clear_su2mat(&(tempmat[i]));
}
 
  
  for (dir = YUP; dir <= TUP; dir++) {
    for (dir2 = XUP; dir2 < dir; dir2++) {
      mtag = start_gather_site(F_OFFSET(link[dir2]), sizeof(su2_matrix),
                               dir, EVENANDODD, gen_pt[0]);
      mtag2 = start_gather_site(F_OFFSET(link[dir]), sizeof(su2_matrix),
                                dir2, EVENANDODD, gen_pt[1]);

      FORALLSITES(i, s) {
        m1 = &(s->link[dir]);
        m4 = &(s->link[dir2]);
        
        mult_su2_an(m4, m1, &(tempmat[i]));
      }
      

      wait_gather(mtag);
      wait_gather(mtag2);
      FORALLSITES(i, s) {
        m1 = (su2_matrix *)(gen_pt[0][i]);
        m4 = (su2_matrix *)(gen_pt[1][i]);
        
        mult_su2_nn(&(tempmat[i]), m1, &tmat);
        complextrace_su2(m4, &tmat,tc);
        td = (*tc).real;  
                                  
                     
        //Needs modification here
        //td += BETA * ( (*tc).real * (*tc).real + (*tc).imag * (*tc).imag ) / 4.0;
        
        if (dir == TUP)
          st_sum += td;
          
        else
          ss_sum += td;
      }
      cleanup_gather(mtag);
      cleanup_gather(mtag2);
      
    }
  }
  
  g_doublesum(&ss_sum);
  g_doublesum(&st_sum);
  free(tc);
  
  *ss_plaq = ss_sum / ((double)volume);
  *st_plaq = st_sum / ((double)volume);

   
}


// -----------------------------------------------------------------
// Print out total action and individual contributions
double action(vector **src, vector ***sol) {
  
  int n;
  double h_act, f_act, g_action;
  double total = 0.0;
  double ssplaq ,stplaq;
  plaquette_a(&ssplaq,&stplaq);

  
  //node0_printf("action: scalar %.8g ", total);

  for (n = 0; n < Nroot; n++) {
    f_act = fermion_action(src[n], sol[n]);
    //node0_printf("fermion%d %.8g ", n, f_act);
    f_act = 0.0;
    total += f_act;
  }

  h_act = hmom_action();
  g_action = -1.0 * (BETA /2.0) * volume * (ssplaq + stplaq);
  node0_printf("mom %.8g ", h_act);
  node0_printf("plaquette %.8g ", g_action);
  total += h_act;
  total += g_action;
  node0_printf("sum %.8g\n", total);
  return total;
}
// -----------------------------------------------------------------
