#ifndef CUDA_FUNCTIONS_R2_H
#define CUDA_FUNCTIONS_R2_H
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "cR2.h"
//General Functions for R2

cR2 sin(const cR2& arg);
cR2 cos(const cR2& arg);
cR2 tan(const cR2& arg);

cR2 sinh(const cR2& arg);
cR2 cosh(const cR2& arg);
cR2 tanh(const cR2& arg);

cR2 sqrt(const cR2& arg);
cR2 exp(const cR2& arg);

cR2 pow(const cR2& base, double exp);
cR2 pow(double base, const cR2& exp);
cR2 pow(const cR2& base, const cR2& exp);

cR2 asin(const cR2& arg);
cR2 acos(const cR2& arg);
cR2 atan(const cR2& arg);
cR2 atan(const cR2& y, const cR2& x);

cR2 erf(const cR2& arg);

cR2 jn(int n, const cR2& arg);

cR2 ln(const cR2& arg);
#endif // !CUDA_FUNCTIONS_R2_H