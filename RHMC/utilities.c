// -----------------------------------------------------------------
// Dirac operator and other helper functions for the action and force
#include "su2_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Copy gauge fields in the site struct
void gauge_field_copy(field_offset src, field_offset dest) {
register int i, dir, src2, dest2;
register site *s;

  FORALLSITES(i, s) {
    src2 = src;
    dest2 = dest;
    for (dir = XUP; dir <= TUP; dir++) {
      su2mat_copy((su2_matrix *)F_PT(s, src2),
                  (su2_matrix *)F_PT(s, dest2));
      src2 += sizeof(su2_matrix);
      dest2 += sizeof(su2_matrix);
    }
  }
}


// -----------------------------------------------------------------
// Matrix--vector operation
// Applies either the operator (sign = 1) or its adjoint (sign = -1)
// Adjoint is simply overall negative sign...
void fermion_op(vector *src, vector *dest,int sign) {
  register int i;
  register site *s;
  
  su2_matrix *m;
  su2_matrix *m1 = malloc(sizeof(*m1));
  int dir ; //L[NDIMS] = {nx,ny,nz,nt};
  //Real tr, link_mass;
  Real vev[DIMF][DIMF];
  vector tvec,tvec_dir, tvec_opp, tvec3,tvec4;
  vector *chi = malloc(sites_on_node *sizeof(*chi));
  vector *psi = malloc(sites_on_node *sizeof(*psi));
  msg_tag *tag[3 * NDIMS];
  //msg_tag *mtag[4];

  // Quick sanity check
  if (sign != 1 && sign != -1) {
    node0_printf("Error: incorrect sign in fermion_op: %d\n", sign);
    terminate(1);
  }

  vev[0][1]=1.0;  
  vev[1][0]=-1.0;
  vev[0][0]=0.0;
  vev[1][1]=0.0;
 

 if(sign == -1){ FORALLSITES(i,s){vec_copy(&(src[i]), &(chi[i])) ; vec_conjug(&(chi[i]),&(src[i]));}}
   
for (dir =XUP ;dir <=TUP ; dir++){
   FORALLSITES(i,s){
    clear_su2mat(&(s->temp_link[dir]));
      
    if(lattice[i].parity == EVEN){su2mat_copy(&(s->link[dir]),&(s->temp_link[dir])) ;}
     
    else{su2_conjug(&(s->link[dir]), &(s->temp_link[dir]));} 
    
   
   }

}
// Start gathers for kinetic term
for (dir = XUP; dir <= TUP; dir++) {

    //if(L[dir] <= 1){continue;}

    tag[dir] = start_gather_field(src, sizeof(vector), dir,
                                  EVENANDODD, gen_pt[dir]);
    tag[OPP_DIR(dir)] = start_gather_field(src, sizeof(vector), OPP_DIR(dir),
                                           EVENANDODD, gen_pt[OPP_DIR(dir)]);

    tag[11-dir] = start_gather_site(F_OFFSET(temp_link[dir]), sizeof(su2_matrix),OPP_DIR(dir), EVENANDODD, gen_pt[11-dir]);
}

//Accumulate (m * epsilon^{ab}  \ psi^{a} \psi{b}) term as gathers run
  FORALLSITES(i, s) {
  clearvec(&(dest[i]));
 
        dest[i].c[0].real = 2* site_mass * vev[0][1] * src[i].c[1].real;
        dest[i].c[0].imag = 2* site_mass * vev[0][1] * src[i].c[1].imag;

        dest[i].c[1].real = 2* site_mass * vev[1][0] * src[i].c[0].real;
        dest[i].c[1].imag = 2* site_mass * vev[1][0] * src[i].c[0].imag;
}

// Accumulate kinetic term as gathers finish
  for (dir = XUP; dir <= TUP; dir++) {
   //if (L[dir]<=1) {continue;}

    wait_gather(tag[dir]);
    wait_gather(tag[OPP_DIR(dir)]);
    wait_gather(tag[11-dir]);
FORALLSITES(i, s) {
      // Deal with BCs here 
      // Need to change the order of BC and matrix multiplication
      // Multiply by matrix and accumulate

      vec_copy((vector *)gen_pt[dir][i], &tvec_dir);
      vec_copy((vector *)gen_pt[OPP_DIR(dir)][i], &tvec_opp);
      
      if (dir == TUP && PBC < 0 && s->t == nt - 1)
     { scalar_mult_vec(&tvec_dir, -1.0, &tvec_dir); } 
      else if (dir == TUP && PBC < 0 && s->t == 0)
     { scalar_mult_vec(&tvec_opp, -1.0, &tvec_opp); }

        mult_su2_mat_vec(&(s->temp_link[dir]),&tvec_dir, &tvec3); //Multiplying the link(x) to psi( x + mu)
        m = (su2_matrix *)(gen_pt[11-dir][i]);
        su2_trans(m ,m1);//Transposing  the link(x-mu)
        mult_su2_mat_vec(m1,&tvec_opp, &tvec4); // Multiplying the link(x-mu) psi( x - mu)
        sub_vec(&tvec3, &tvec4, &tvec);
        scalar_mult_add_vec(&(dest[i]), &tvec, 0.5 * s->phase[dir], &(dest[i]));
 }
    cleanup_gather(tag[dir]);
    cleanup_gather(tag[OPP_DIR(dir)]);
    cleanup_gather(tag[11-dir]);
}
// Overall negative sign and conjugate for adjoint
if (sign == -1) {
    FORALLSITES(i, s)
      scalar_mult_vec(&(dest[i]), -1.0, &(dest[i]));
      vec_copy(&(dest[i]),&(psi[i]));
      vec_conjug(&(psi[i]),&(dest[i]));
}
free(m1);
free(chi);
free(psi);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Squared four-fermion matrix--vector operation
//   dest = D^2 src
// Use tempvec for temporary storage
void DSq(vector *src, vector *dest) {
  fermion_op(src, tempvec,PLUS);
  fermion_op(tempvec,dest,MINUS); 
}
// -----------------------------------------------------------------
