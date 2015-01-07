#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE
*/

double
hermite(double x, int n, int k)
{
	if (n == 0) {
		if (k == 0)
			return 1.0;
		else if (k >= 1)
			return 0.0;
	} else if (n == 1) {
		if (k == 0)
			return 2.0*x;
		else if (k == 1)
    			return 2.0;
		else if (k >= 2)
    			return 0.0;
	} else {
		if (k == 0)
			return 2.0*x*hermite(x, n-1, 0) - 2.0*(n-1)*hermite(x, n-2, 0);
		else
			return 2.0*k*hermite(x, n-1, k-1) + 2.0*x*hermite(x, n-1, k) - 2.0*(n-1)*hermite(x, n-2, k);
	}
	fprintf (stderr, "Bład\n");
	return -0.666;
}

void
make_spl(points_t * pts, spline_t * spl)
{

	matrix_t       *eqs = NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	eqs = make_matrix(nb, nb + 1);

	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, hermite(x[k], i, 0) * hermite(x[k], j, 0));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * hermite(x[k], j, 0));
	}

	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}

	if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = a + i*(b-a)/(spl->n-1);
			spl->x[i] = xx;
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * hermite(xx, k, 0);
				spl->f1[i] += ck * hermite(xx, k, 1);
				spl->f2[i] += ck * hermite(xx, k, 2);
				spl->f3[i] += ck * hermite(xx, k, 3);
			}
		}
	}
}
