// -----------------------------------------------------------------
// Update the momentum matrices
#include "su2_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
double gauge_force(Real eps) {
  register int i, dir1, dir2;
  register site *st;
  msg_tag *tag0, *tag1, *tag2;
  int start;
  su2_matrix tmat1, tmat2,tmat3,tmat4,tmat5;
  complex *ctmp = malloc(sizeof(*ctmp));
  double norm = 0;

  // Loop over directions, update mom[dir1]
  for (dir1 = XUP; dir1 <= TUP; dir1++) {
    start = 1; // Indicates staple sum not initialized
    FORALLSITES(i, st)
      clear_su2mat(&(st->staple));    // Initialize staple

    // Loop over other directions
    // Compute force from plaquettes in the dir1, dir2 plane
    for (dir2 = XUP; dir2 <= TUP; dir2++) {
      if (dir2 != dir1) {
        // Get link[dir2] from direction dir1
        tag0 = start_gather_site(F_OFFSET(link[dir2]),
                                 sizeof(su2_matrix),
                                 dir1, EVENANDODD, gen_pt[0]);

        // Start gather for the "upper staple"
       tag2 = start_gather_site(F_OFFSET(link[dir1]),
                                sizeof(su2_matrix),
                                dir2, EVENANDODD, gen_pt[2]);

        // Begin the computation "at the dir2DOWN point"
        // We will later gather the intermediate result "to the home point"
        wait_gather(tag0);
        FORALLSITES(i, st) {
          mult_su2_an(&(st->link[dir2]), &(st->link[dir1]), &tmat1);
          mult_su2_nn(&tmat1, (su2_matrix *)gen_pt[0][i],
                      (su2_matrix *)&(st->tempmat1));
        }

        // Gather lower staple "up to home site"
        tag1 = start_gather_site(F_OFFSET(tempmat1), sizeof(su2_matrix),
                                 OPP_DIR(dir2), EVENANDODD, gen_pt[1]);

        // The "upper" staple
        // One of the links has already been gathered,
        // since it was used in computing
        // the "lower" staple of the site above (in dir2)
        wait_gather(tag2);
        if (start) {  // This is the first contribution to staple
          FORALLSITES(i, st) {
            mult_su2_nn(&(st->link[dir2]), (su2_matrix *)gen_pt[2][i], &tmat1);
            mult_su2_na(&tmat1, (su2_matrix *)gen_pt[0][i], &(st->staple));
            scalar_mult_su2_matrix(&(st->staple), 1.0 ,&(st->staple));
          }
          start = 0;
        }
        else {
          FORALLSITES(i, st) {
            mult_su2_nn(&(st->link[dir2]), (su2_matrix *)gen_pt[2][i], &tmat1);
            mult_su2_na(&tmat1, (su2_matrix *)gen_pt[0][i], &tmat2);
            scalar_mult_add_su2_matrix(&(st->staple),&tmat2, 1.0 ,&(st->staple));
          }
        }

        wait_gather(tag1);
        FORALLSITES(i, st) {

          scalar_mult_add_su2_matrix(&(st->staple),(su2_matrix *)gen_pt[1][i], 1.0 ,&(st->staple));
        }
        cleanup_gather(tag0);
        cleanup_gather(tag1);
        cleanup_gather(tag2);
      }
    } // End of loop over other directions

    // Now multiply the staple sum by the link, then update momentum
    FORALLSITES(i, st) {
      mult_su2_na(&(st->staple), &(st->link[dir1]) ,&tmat3);
      su2_adjoint(&(st->staple) ,&tmat5);
      mult_su2_nn(&(st->link[dir1]), &tmat5, &tmat1);

      scalar_mult_add_su2_matrix(&tmat1 ,&tmat3 , -1.0 , &(st->staple));
      unit_su2mat(&tmat4);
      scalar_mult_su2_matrix(&tmat4 , -0.5 ,&tmat4);
      trace_su2(&(st->staple),ctmp);
      c_scalar_mult_add_su2mat(&(st->staple),&tmat4,ctmp,&(st->staple));

      scalar_mult_su2_matrix(&(st->staple) , (BETA/4.0) ,&tmat1);

      //Finally updating the momenta
      uncompress_anti_hermitian(&(st->mom[dir1]), &tmat2);
      scalar_mult_add_su2_matrix(&tmat2, &tmat1, eps, &(st->staple));
      make_anti_hermitian(&(st->staple), &(st->mom[dir1]));
      norm += (double)realtrace_su2(&tmat1, &tmat1);
    }
  } // End of loop over dir1

  g_doublesum(&norm);
  free(ctmp);
  return (eps * sqrt(norm) / volume);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Assume CG has been run and the solution is in psim[n]
double fermion_force(Real eps, vector *src, vector **sol){
  register int i, dir;
  register site *s,*st;
  int n;
  Real tr;
  double norm;
  Real ferm_epsilon = 2.0* eps;
  complex *dum = malloc(sizeof(*dum));
  vector tvec, tvec_sol;
  su2_matrix *a = malloc(sizeof(*a));
  su2_matrix *b = malloc(sizeof(*b));
  su2_matrix tmat, tmat1, tmat2, tmat3;
  msg_tag *tag0[NDIMS] ,*tag1[NDIMS];
  vector *psol= malloc(sites_on_node * sizeof(vector));

  // Zero the force collectors
  for (dir = XUP; dir <= TUP; dir++) {
    FORALLSITES(i, s)
      clear_su2mat(&sigma[dir][i]);
  }

  for(i = 0; i < NUMGEN; i++)
    clear_su2mat(&(temp_Lambda[i]));

  //Hit the solution vectors with the fermion operator
  /*for (n = 0; n < Norder; n++) {
    fermion_op(sol[n], psol, PLUS);

    for (dir = XUP; dir <= TUP; dir++) {

  //Start gathers for psol[n] = M sol[n]  , and sol[n]
  tag0[dir] = start_gather_field(psol, sizeof(vector), dir,
  EVENANDODD, gen_pt[dir]);

  tag1[dir] = start_gather_field(sol[n], sizeof(vector), dir,
  EVENANDODD, gen_pt[3+dir]);


  wait_gather(tag0[dir]);
  wait_gather(tag1[dir]);

  FORALLSITES(i,s){

  vec_copy((vector *)gen_pt[dir][i], &tvec);

  vec_copy((vector *)gen_pt[3+dir][i], &tvec_sol);

  if (dir == TUP && PBC < 0 && s->t == nt - 1)
  { scalar_mult_vec(&tvec, -1.0, &tvec);
  scalar_mult_vec(&tvec_sol, -1.0, &tvec_sol);}


  su2_projector(&tvec_sol, &(psol[i]),a);



  for(int c=0; c < NUMGEN; c++){

  if(s->parity == EVEN ){
  mult_su2_nn(&(s->link[dir]),a,&tmat);
  mult_su2_nn(&(Lambda[c]),&tmat,&tmat1);
  trace_su2(&tmat1,dum); }
  else{
  su2_conjug(&(s->link[dir]),&(s->temp_link1[dir]));
  mult_su2_nn(&(s->temp_link1[dir]),a,&tmat);
  su2_conjug(&(Lambda[c]),&(temp_Lambda[c]));
  mult_su2_nn(&(temp_Lambda[c]), &tmat,&tmat1);
  trace_su2(&tmat1,dum);}

  tr =-0.5 * (*dum).real * s->phase[dir];
  scalar_mult_add_su2_matrix(b,&(Lambda[c]),tr,b);
  }


  clear_su2mat(a);
  su2_projector(&(sol[n][i]),&tvec,a);


  for(int d=0; d < NUMGEN; d++){
  if(s->parity == EVEN ){
  mult_su2_nn(&(s->link[dir]),a,&tmat);
  mult_su2_nn(&(Lambda[d]),&tmat,&tmat1);
  trace_su2(&tmat1,dum);}
  else{
  su2_conjug(&(s->link[dir]),&(s->temp_link1[dir]));
  mult_su2_nn(&(s->temp_link1[dir]),a,&tmat);
  su2_conjug(&(Lambda[d]),&(temp_Lambda[d]));
  mult_su2_nn(&(temp_Lambda[d]),&tmat,&tmat1);
  trace_su2(&tmat1,dum);
  }

  tr = 0.5 * (*dum).real * s->phase[dir];
  scalar_mult_add_su2_matrix(b,&(Lambda[d]),tr,b);
  }

  scalar_mult_add_su2_matrix(&sigma[dir][i], b,-2.0*amp4[n],&sigma[dir][i]);

}

cleanup_gather(tag0[dir]);
cleanup_gather(tag1[dir]);


}

}

*/

  //Update the momenta with the gauge force
  for (dir = XUP; dir <= TUP; dir++) {
    FORALLSITES(i, st) {
      uncompress_anti_hermitian(&(st->mom[dir]), &tmat2);// Transform mom(an anti-hermitian matrix) to an SU(2) matrix
      scalar_mult_add_su2_matrix(&tmat2, &sigma[dir][i], ferm_epsilon, &tmat3);
      make_anti_hermitian(&tmat3, &(st->mom[dir]));
      norm += (double)realtrace_su2(&sigma[dir][i], &sigma[dir][i]);
    }
  }
  g_doublesum(&norm);

  free(psol);
  free(a);
  free(b);
  free(dum);

  return ferm_epsilon * sqrt(norm) / volume * 2;
}
// -----------------------------------------------------------------

