// Host emulation of the CUDA/HIP Duprat kernels: one thread, one block.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace std;
#define __device__
#define __global__
struct dim3_t { unsigned x, y, z; };
static dim3_t blockIdx{0, 0, 0}, blockDim{1, 1, 1}, threadIdx{0, 0, 0},
    gridDim{1, 1, 1};
#include KERNEL

// 5-point Gauss-Legendre rule on [0, 1] in closed form
static void gauss5(double *x, double *w) {
  const double r = sqrt(10.0 / 7.0), s70 = sqrt(70.0);
  const double a = sqrt(5.0 - 2.0 * r) / 6.0, b = sqrt(5.0 + 2.0 * r) / 6.0;
  x[0] = 0.5 - b; x[1] = 0.5 - a; x[2] = 0.5; x[3] = 0.5 + a; x[4] = 0.5 + b;
  w[0] = w[4] = (322.0 - 13.0 * s70) / 1800.0;
  w[1] = w[3] = (322.0 + 13.0 * s70) / 1800.0;
  w[2] = 64.0 / 225.0;
}

template <typename T> void run(const char *out) {
  FILE *f = fopen("inputs.txt", "r");
  int n;
  if (fscanf(f, "%d", &n) != 1) exit(1);
  vector<vector<T>> c(16, vector<T>(n));
  for (int i = 0; i < n; ++i)
    for (int k = 0; k < 16; ++k) {
      double d;
      if (fscanf(f, "%lf", &d) != 1) exit(1);
      c[k][i] = (T)d;
    }
  fclose(f);
  vector<T> tx(n), ty(n), tz(n), al(n), d0(n, (T)0.123), d1(n, (T)0.123);
  double gx[5], gw[5];
  gauss5(gx, gw);
  vector<T> glx(gx, gx + 5), glw(gw, gw + 5);
  FILE *o = fopen(out, "w");
  for (int k = 0; k < 2; ++k) {
    if (k == 0) {
      fill(tx.begin(), tx.end(), (T)0); fill(ty.begin(), ty.end(), (T)0);
      fill(tz.begin(), tz.end(), (T)0);
    } else {
      tx = c[13]; ty = c[14]; tz = c[15];
    }
    duprat_compute<T>(c[0].data(), c[1].data(), c[2].data(), c[3].data(),
                      c[4].data(), c[5].data(), c[6].data(), c[7].data(),
                      c[8].data(), c[9].data(), tx.data(), ty.data(),
                      tz.data(), al.data(), n, (T)0.41, (T)0.78, (T)17.0,
                      glx.data(), glw.data(), 5,
                      k == 0 ? 1 : 5);
    for (int i = 0; i < n; ++i)
      fprintf(o, "%.17e %.17e %.17e %.17e\n", (double)tx[i], (double)ty[i],
              (double)tz[i], (double)al[i]);
  }
  duprat_update_dpds<T>(c[10].data(), c[11].data(), c[12].data(), c[0].data(),
                        c[1].data(), c[2].data(), c[3].data(), c[4].data(),
                        c[5].data(), d0.data(), n);
  for (int i = 0; i < n; ++i)
    fprintf(o, "%.17e\n", (double)d0[i]);
  fclose(o);
}

int main(int argc, char **argv) {
  run<double>(argv[1]);
  run<float>(argv[2]);
  return 0;
}
