// Host emulation of the Metal Duprat kernels, float only.
#include <cstdio>
#include <cstdlib>
#include <vector>
#define kernel
#define device
#define constant const
#include KERNEL
using std::vector;
using std::sqrt;

// 5-point Gauss-Legendre rule on [0, 1] in closed form
static void gauss5(double *x, double *w) {
  const double r = sqrt(10.0 / 7.0), s70 = sqrt(70.0);
  const double a = sqrt(5.0 - 2.0 * r) / 6.0, b = sqrt(5.0 + 2.0 * r) / 6.0;
  x[0] = 0.5 - b; x[1] = 0.5 - a; x[2] = 0.5; x[3] = 0.5 + a; x[4] = 0.5 + b;
  w[0] = w[4] = (322.0 - 13.0 * s70) / 1800.0;
  w[1] = w[3] = (322.0 + 13.0 * s70) / 1800.0;
  w[2] = 64.0 / 225.0;
}

int main(int argc, char **argv) {
  FILE *f = fopen("inputs.txt", "r");
  int n;
  if (fscanf(f, "%d", &n) != 1) return 1;
  vector<vector<float>> c(16, vector<float>(n));
  for (int i = 0; i < n; ++i)
    for (int k = 0; k < 16; ++k) {
      double d;
      if (fscanf(f, "%lf", &d) != 1) return 1;
      c[k][i] = (float)d;
    }
  fclose(f);
  vector<float> tx(n), ty(n), tz(n), al(n), d0(n, 0.123f), d1(n, 0.123f);
  double gx[5], gw[5];
  gauss5(gx, gw);
  vector<float> glx(gx, gx + 5), glw(gw, gw + 5);
  FILE *o = fopen(argv[1], "w");
  const float kappa = 0.41f, beta = 0.78f, A = 17.0f;
  for (int k = 0; k < 2; ++k) {
    if (k == 0) {
      std::fill(tx.begin(), tx.end(), 0.f); std::fill(ty.begin(), ty.end(), 0.f);
      std::fill(tz.begin(), tz.end(), 0.f);
    } else {
      tx = c[13]; ty = c[14]; tz = c[15];
    }
    const int tstep = k == 0 ? 1 : 5;
    for (int i = 0; i < n; ++i)
      duprat_compute_kernel(c[0].data(), c[1].data(), c[2].data(), c[3].data(),
                            c[4].data(), c[5].data(), c[6].data(), c[7].data(),
                            c[8].data(), c[9].data(), tx.data(), ty.data(),
                            tz.data(), al.data(), n, kappa, beta, A, glx.data(),
                            glw.data(), 5, tstep,
                            (uint)i);
    for (int i = 0; i < n; ++i)
      fprintf(o, "%.17e %.17e %.17e %.17e\n", (double)tx[i], (double)ty[i],
              (double)tz[i], (double)al[i]);
  }
  for (int i = 0; i < n; ++i)
    duprat_update_dpds_kernel(c[10].data(), c[11].data(), c[12].data(),
                              c[0].data(), c[1].data(), c[2].data(),
                              c[3].data(), c[4].data(), c[5].data(),
                              d0.data(), n, (uint)i);
  for (int i = 0; i < n; ++i)
    fprintf(o, "%.17e\n", (double)d0[i]);
  fclose(o);
  return 0;
}
