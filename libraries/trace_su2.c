// -----------------------------------------------------------------
// Complex trace of a matrix
// s <-- Tr[a]
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su2.h"

void trace_su2(su2_matrix *a, complex *s) {
  CADD(a->e[0][0], a->e[1][1], *s);
}
// -----------------------------------------------------------------
