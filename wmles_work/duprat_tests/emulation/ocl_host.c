/* Run the OpenCL Duprat kernels on PoCL. argv: kernel source (already with
   real replaced), real size (4 or 8), output file. */
#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CK(x) do { cl_int e_ = (x); if (e_ != CL_SUCCESS) { fprintf(stderr, "CL error %d at %d\n", e_, __LINE__); exit(1);} } while (0)

static char *slurp(const char *p) {
  FILE *f = fopen(p, "rb"); fseek(f, 0, SEEK_END); long n = ftell(f);
  rewind(f); char *b = malloc(n + 1); if (fread(b, 1, n, f) != (size_t)n) exit(1);
  b[n] = 0; fclose(f); return b;
}

int main(int argc, char **argv) {
  const char *src = slurp(argv[1]);
  const int rs = atoi(argv[2]);
  cl_platform_id pl; cl_device_id dev; cl_int err;
  CK(clGetPlatformIDs(1, &pl, NULL));
  CK(clGetDeviceIDs(pl, CL_DEVICE_TYPE_ALL, 1, &dev, NULL));
  cl_context ctx = clCreateContext(NULL, 1, &dev, NULL, NULL, &err); CK(err);
  cl_command_queue q = clCreateCommandQueue(ctx, dev, 0, &err); CK(err);
  cl_program prg = clCreateProgramWithSource(ctx, 1, &src, NULL, &err); CK(err);
  err = clBuildProgram(prg, 1, &dev, NULL, NULL, NULL);
  if (err != CL_SUCCESS) {
    size_t ls; clGetProgramBuildInfo(prg, dev, CL_PROGRAM_BUILD_LOG, 0, NULL, &ls);
    char *log = malloc(ls); clGetProgramBuildInfo(prg, dev, CL_PROGRAM_BUILD_LOG, ls, log, NULL);
    fprintf(stderr, "%s\n", log); return 1;
  }
  FILE *f = fopen("inputs.txt", "r"); int n;
  if (fscanf(f, "%d", &n) != 1) return 1;
  double *in = malloc(sizeof(double) * 16 * n);
  for (int i = 0; i < n; ++i)
    for (int k = 0; k < 16; ++k)
      if (fscanf(f, "%lf", &in[k * n + i]) != 1) return 1;
  fclose(f);
  /* host arrays in the real type */
  void *h[16]; cl_mem d[16];
  for (int k = 0; k < 16; ++k) {
    h[k] = malloc(rs * n);
    for (int i = 0; i < n; ++i) {
      if (rs == 8) ((double *)h[k])[i] = in[k * n + i];
      else ((float *)h[k])[i] = (float)in[k * n + i];
    }
    d[k] = clCreateBuffer(ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, rs * n, h[k], &err); CK(err);
  }
  cl_mem tx = clCreateBuffer(ctx, CL_MEM_READ_WRITE, rs * n, NULL, &err);
  cl_mem ty = clCreateBuffer(ctx, CL_MEM_READ_WRITE, rs * n, NULL, &err);
  cl_mem tz = clCreateBuffer(ctx, CL_MEM_READ_WRITE, rs * n, NULL, &err);
  cl_mem al = clCreateBuffer(ctx, CL_MEM_READ_WRITE, rs * n, NULL, &err);
  cl_mem dd[2];
  void *zero = calloc(n, rs), *init = malloc(rs * n), *o4[4];
  for (int i = 0; i < n; ++i) { if (rs == 8) ((double *)init)[i] = 0.123; else ((float *)init)[i] = 0.123f; }
  for (int j = 0; j < 2; ++j) { dd[j] = clCreateBuffer(ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, rs * n, init, &err); CK(err); }
  for (int j = 0; j < 4; ++j) o4[j] = malloc(rs * n);
  double kd[3] = {0.41, 0.78, 17.0}; float kf[3] = {0.41f, 0.78f, 17.0f};
  const size_t gsz = 256 * ((n + 255) / 256), lsz = 256;
  FILE *o = fopen(argv[3], "w");
  for (int pass = 0; pass < 2; ++pass) {
    CK(clEnqueueWriteBuffer(q, tx, CL_TRUE, 0, rs * n, pass ? h[13] : zero, 0, NULL, NULL));
    CK(clEnqueueWriteBuffer(q, ty, CL_TRUE, 0, rs * n, pass ? h[14] : zero, 0, NULL, NULL));
    CK(clEnqueueWriteBuffer(q, tz, CL_TRUE, 0, rs * n, pass ? h[15] : zero, 0, NULL, NULL));
    cl_kernel kr = clCreateKernel(prg, "duprat_compute_kernel", &err); CK(err);
    for (int k = 0; k < 10; ++k) CK(clSetKernelArg(kr, k, sizeof(cl_mem), &d[k]));
    CK(clSetKernelArg(kr, 10, sizeof(cl_mem), &tx)); CK(clSetKernelArg(kr, 11, sizeof(cl_mem), &ty));
    CK(clSetKernelArg(kr, 12, sizeof(cl_mem), &tz)); CK(clSetKernelArg(kr, 13, sizeof(cl_mem), &al));
    CK(clSetKernelArg(kr, 14, sizeof(int), &n));
    for (int k = 0; k < 3; ++k) CK(clSetKernelArg(kr, 15 + k, rs, rs == 8 ? (void *)&kd[k] : (void *)&kf[k]));
    int tstep = pass ? 5 : 1; CK(clSetKernelArg(kr, 18, sizeof(int), &tstep));
    CK(clEnqueueNDRangeKernel(q, kr, 1, NULL, &gsz, &lsz, 0, NULL, NULL));
    CK(clFinish(q));
    cl_mem outs[4] = {tx, ty, tz, al};
    for (int j = 0; j < 4; ++j) CK(clEnqueueReadBuffer(q, outs[j], CL_TRUE, 0, rs * n, o4[j], 0, NULL, NULL));
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < 4; ++j) fprintf(o, "%.17e ", rs == 8 ? ((double *)o4[j])[i] : (double)((float *)o4[j])[i]);
      fprintf(o, "\n");
    }
    clReleaseKernel(kr);
  }
  for (int lim = 0; lim < 2; ++lim) {
    cl_kernel ku = clCreateKernel(prg, "duprat_update_dpds_kernel", &err); CK(err);
    int idx[12] = {10, 11, 12, 0, 1, 2, 3, 4, 5, 6, 7, 8};
    for (int k = 0; k < 12; ++k) CK(clSetKernelArg(ku, k, sizeof(cl_mem), &d[idx[k]]));
    CK(clSetKernelArg(ku, 12, sizeof(int), &lim));
    CK(clSetKernelArg(ku, 13, sizeof(cl_mem), &dd[lim]));
    CK(clSetKernelArg(ku, 14, sizeof(int), &n));
    CK(clEnqueueNDRangeKernel(q, ku, 1, NULL, &gsz, &lsz, 0, NULL, NULL));
    CK(clFinish(q));
    CK(clEnqueueReadBuffer(q, dd[lim], CL_TRUE, 0, rs * n, o4[lim], 0, NULL, NULL));
    clReleaseKernel(ku);
  }
  for (int i = 0; i < n; ++i)
    fprintf(o, "%.17e %.17e\n", rs == 8 ? ((double *)o4[0])[i] : (double)((float *)o4[0])[i],
            rs == 8 ? ((double *)o4[1])[i] : (double)((float *)o4[1])[i]);
  fclose(o);
  return 0;
}
