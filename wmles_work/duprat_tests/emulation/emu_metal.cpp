// Host emulation of the Metal Duprat kernels, float only.
#include <cstdio>
#include <cstdlib>
#include <vector>
#define kernel
#define device
#define constant const
#include KERNEL
using std::vector;

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
                            tz.data(), al.data(), n, kappa, beta, A, tstep,
                            (uint)i);
    for (int i = 0; i < n; ++i)
      fprintf(o, "%.17e %.17e %.17e %.17e\n", (double)tx[i], (double)ty[i],
              (double)tz[i], (double)al[i]);
  }
  for (int lim = 0; lim < 2; ++lim)
    for (int i = 0; i < n; ++i)
      duprat_update_dpds_kernel(c[10].data(), c[11].data(), c[12].data(),
                                c[0].data(), c[1].data(), c[2].data(),
                                c[3].data(), c[4].data(), c[5].data(),
                                c[6].data(), c[7].data(), c[8].data(), lim,
                                lim ? d1.data() : d0.data(), n, (uint)i);
  for (int i = 0; i < n; ++i)
    fprintf(o, "%.17e %.17e\n", (double)d0[i], (double)d1[i]);
  fclose(o);
  return 0;
}
