// -----------------------------------------------------------------
// Update lattice
// Leapfrog integrator
// Begin at "integral" time, with H and U evaluated at the same time

// Uncomment to print out debugging messages
#define UPDATE_DEBUG
#include "su2_includes.h"
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>         // For "finite"
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update_step(Real *fnorm, Real *snorm, vector **src, vector ***psim) {
  int step, iters = 0, n;
  Real final_rsq, eps = traj_length / (Real)nsteps[0], tr;
  node0_printf("eps %.4g\n", eps);

  // First u(t/2)

  update_u(0.5 * eps);

  for (step = 0; step < nsteps[0]; step++) {
    // Inner steps p(t) u(t)
    tr = gauge_force(eps); // Updates the gauge_momenta at the end of function
    *snorm += tr;
    if (tr > max_sf)
      max_sf = tr;

    for (n = 0; n < Nroot; n++) {
      // Do conjugate gradient to get (Mdag M)^(-1 / 4) chi
      iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);
      tr = fermion_force(eps, src[n], psim[n]);  //Updates the gauge_momenta at the end of function
      fnorm[n] += tr;
      if (tr > max_ff[n])
        max_ff[n] = tr;
    }

    if (step < nsteps[0] - 1)
      update_u(eps);
    else                // Final u(t/2)
      update_u(0.5 * eps);
  }
  return iters;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int update() {
  int j, n, iters = 0;
  Real final_rsq;
  double startaction, endaction, change;
  vector **src = malloc(Nroot * sizeof(**src));
  vector ***psim = malloc(Nroot * sizeof(***psim));

  for (n = 0; n < Nroot; n++) {
    src[n] = malloc(sites_on_node * sizeof(vector));
    psim[n] = malloc(Norder * sizeof(vector*));
    for (j = 0; j < Norder; j++)
      psim[n][j] = malloc(sites_on_node * sizeof(vector));
  }

  // Refresh the momenta
  ranmom();

  // Set up the fermion variables
  // Compute g and src = (Mdag M)^(1 / 8) g
  for (n = 0; n < Nroot; n++)
    iters += grsource(src[n]);

  // Do a CG to get psim,
  // rational approximation to (Mdag M)^(-1 / 4) src = (Mdag M)^(-1 / 8) g
  for (j = 0; j < Norder; j++)
    shift[j] = shift4[j];


#ifdef UPDATE_DEBUG
  node0_printf("Calling CG in update_leapfrog -- original action\n");
#endif
  // congrad_multi initializes psim
  for (n = 0; n < Nroot; n++)
    iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);


  node0_printf("No of iteration in CG = %.1d\n",iters);

  // Find initial action
  startaction = action(src, psim);
  snorm = 0.0;
  max_sf = 0.0;
  for (n = 0; n < Nroot; n++) {
    fnorm[n] = 0.0;
    max_ff[n] = 0.0;
  }

#ifdef HMC_ALGORITHM
  Real xrandom;   // For accept/reject test
  // Copy gauge field to old_link
  gauge_field_copy(F_OFFSET(link[0]), F_OFFSET(old_link[0]));

#endif
  // Do microcanonical updating
  iters += update_step(fnorm, &snorm, src, psim);

  // Find ending action
  // Since update_step ended on a gauge update,
  // need to do conjugate gradient to get (Mdag M)^(-1 / 4) chi
#ifdef UPDATE_DEBUG
  node0_printf("Calling CG in update_leapfrog -- new action\n");
#endif
  for (n = 0; n < Nroot; n++)
    iters += congrad_multi(src[n], psim[n], niter, rsqmin, &final_rsq);

  endaction = action(src, psim);
  change = endaction - startaction;
#ifdef HMC_ALGORITHM
  // Reject configurations giving overflow
#ifndef HAVE_IEEEFP_H
  if (fabs((double)change) > 1e20) {
#else
  if (!finite((double)change)) {
#endif
    node0_printf("WARNING: Correcting Apparent Overflow: Delta S = %.4g\n",
                 change);
    change = 1.0e20;
  }

  // Decide whether to accept, if not, copy old link field back
  // Careful -- must generate only one random number for whole lattice
  if (this_node == 0)
    xrandom = myrand(&node_prn);
  broadcast_float(&xrandom);
  if (exp(-change) < (double)xrandom) {
    if (traj_length > 0.0) {
      gauge_field_copy(F_OFFSET(old_link[0]), F_OFFSET(link[0]));
    }
    node0_printf("REJECT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startaction, endaction);
  }
  else {
    node0_printf("ACCEPT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startaction, endaction);
  }
#else
  // Only print check if not doing HMC
  node0_printf("CHECK: delta S = %.4g\n", (double)(change));
#endif // ifdef HMC

  for (n = 0; n < Nroot; n++) {
    free(src[n]);
    for (j = 0; j < Norder; j++)
      free(psim[n][j]);
    free(psim[n]);
  }
  free(src);
  free(psim);

  if (traj_length > 0) {
    node0_printf("MONITOR_FORCE_GAUGE   %.4g %.4g\n",
                 snorm / (double)(2 * nsteps[0]), max_sf);
    for (n = 0; n < Nroot; n++) {
      node0_printf("MONITOR_FORCE_FERMION%d %.4g %.4g\n",
                   n, fnorm[n] / (double)(2 * nsteps[0]), max_ff[n]);
    }
    return iters;
  }
  else
    return -99;
}
// -----------------------------------------------------------------
