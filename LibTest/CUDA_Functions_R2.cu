#include "CUDA_Functions_R2.h"

#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "stdud.h"
#include "VecR2.h"
using namespace std;

__global__ void cR2_SIN(double* b, double* a, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		b[i] = sin(a[i]);
	}
}

__global__ void cR2_COS(double* b, double* a, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		b[i] = cos(a[i]);
	}
}

__global__ void cR2_TAN(double* b, double* a, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		if (cos(a[i]) != 0.0) {
			b[i] = sin(a[i]) / cos(a[i]);
		}
		else {
			b[i] = 0.0;
		}
	}
}


__global__ void cR2_SINH(double* b, double* a, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		b[i] = sinh(a[i]);
	}
}

__global__ void cR2_COSH(double* b, double* a, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		b[i] = cosh(a[i]);
	}
}

__global__ void cR2_TANH(double* b, double* a, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		b[i] = sinh(a[i]) / cosh(a[i]);
	}
}

__global__ void cR2_SQRT(double* b, double* a, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		b[i] = sqrt(a[i]);
	}
}

__global__ void cR2_EXP(double* b, double* a, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		b[i] = exp(a[i]);
	}
}

__global__ void cR2_POW(double* c, double* a, double b, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		c[i] = pow(a[i], b);
	}
}

__global__ void cR2_POW(double* c, double a, double* b, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		c[i] = pow(a, b[i]);
	}
}

__global__ void cR2_POW(double* c, double* a, double* b, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		c[i] = pow(a[i],b[i]);
	}
}

__global__ void cR2_ASIN(double* b, double* a, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		if (abs(a[i]) < 1.0) {
			b[i] = asin(a[i]);
		}
		else {
			b[i] = 0.0;
		}
	}
}

__global__ void cR2_ACOS(double* b, double* a, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		if (abs(a[i]) < 1.0) {
			b[i] = acos(a[i]);
		}
		else {
			b[i] = 0.0;
		}
	}
}

__global__ void cR2_ATAN(double* b, double* a, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		b[i] = atan(a[i]);
	}
}

__global__ void cR2_ATAN(double* b, double* y, double* x, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		b[i] = atan2(y[i], x[i]);
	}
}

__global__ void cR2_ERF(double* b, double* a, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		b[i] = erf(a[i]);
	}
}

__global__ void cR2_J_n(double* b, int n, double* a, int N) {
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	if (i < N) {
		b[i] = jn(n, a[i]);
	}
}

cR2 sin(const cR2 & arg) {
	cR2 temp(arg.getIndex1(), arg.getIndex2());
	temp.setDX(arg.getDX()); temp.setDY(arg.getDY());
	temp.setX1(arg.getX1()); temp.setX2(arg.getX2());
	temp.setY1(arg.getY1()); temp.setY2(arg.getY2());
	temp.setG(arg.getG()); temp.setB(arg.getB());
	//Execute GPU duals:
	cR2_SIN << <arg.getG(), arg.getB() >> > (temp.getD_ptr(), arg.getD_ptr(), arg.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 cos(const cR2 & arg) {
	cR2 temp(arg.getIndex1(), arg.getIndex2());
	temp.setDX(arg.getDX()); temp.setDY(arg.getDY());
	temp.setX1(arg.getX1()); temp.setX2(arg.getX2());
	temp.setY1(arg.getY1()); temp.setY2(arg.getY2());
	temp.setG(arg.getG()); temp.setB(arg.getB());
	//Execute GPU duals:
	cR2_COS << <arg.getG(), arg.getB() >> > (temp.getD_ptr(), arg.getD_ptr(), arg.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 tan(const cR2 & arg) {
	cR2 temp(arg.getIndex1(), arg.getIndex2());
	temp.setDX(arg.getDX()); temp.setDY(arg.getDY());
	temp.setX1(arg.getX1()); temp.setX2(arg.getX2());
	temp.setY1(arg.getY1()); temp.setY2(arg.getY2());
	temp.setG(arg.getG()); temp.setB(arg.getB());
	//Execute GPU duals:
	cR2_TAN << <arg.getG(), arg.getB() >> > (temp.getD_ptr(), arg.getD_ptr(), arg.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 sinh(const cR2 & arg) {
	cR2 temp(arg.getIndex1(), arg.getIndex2());
	temp.setDX(arg.getDX()); temp.setDY(arg.getDY());
	temp.setX1(arg.getX1()); temp.setX2(arg.getX2());
	temp.setY1(arg.getY1()); temp.setY2(arg.getY2());
	temp.setG(arg.getG()); temp.setB(arg.getB());
	//Execute GPU duals:
	cR2_SINH << <arg.getG(), arg.getB() >> > (temp.getD_ptr(), arg.getD_ptr(), arg.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 cosh(const cR2 & arg) {
	cR2 temp(arg.getIndex1(), arg.getIndex2());
	temp.setDX(arg.getDX()); temp.setDY(arg.getDY());
	temp.setX1(arg.getX1()); temp.setX2(arg.getX2());
	temp.setY1(arg.getY1()); temp.setY2(arg.getY2());
	temp.setG(arg.getG()); temp.setB(arg.getB());
	//Execute GPU duals:
	cR2_COSH << <arg.getG(), arg.getB() >> > (temp.getD_ptr(), arg.getD_ptr(), arg.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 tanh(const cR2 & arg) {
	cR2 temp(arg.getIndex1(), arg.getIndex2());
	temp.setDX(arg.getDX()); temp.setDY(arg.getDY());
	temp.setX1(arg.getX1()); temp.setX2(arg.getX2());
	temp.setY1(arg.getY1()); temp.setY2(arg.getY2());
	temp.setG(arg.getG()); temp.setB(arg.getB());
	//Execute GPU duals:
	cR2_TANH << <arg.getG(), arg.getB() >> > (temp.getD_ptr(), arg.getD_ptr(), arg.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 sqrt(const cR2 & arg) {
	cR2 temp(arg.getIndex1(), arg.getIndex2());
	temp.setDX(arg.getDX()); temp.setDY(arg.getDY());
	temp.setX1(arg.getX1()); temp.setX2(arg.getX2());
	temp.setY1(arg.getY1()); temp.setY2(arg.getY2());
	temp.setG(arg.getG()); temp.setB(arg.getB());
	//Execute GPU duals:
	cR2_SQRT << <arg.getG(), arg.getB() >> > (temp.getD_ptr(), arg.getD_ptr(), arg.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 exp(const cR2 & arg) {
	cR2 temp(arg.getIndex1(), arg.getIndex2());
	temp.setDX(arg.getDX()); temp.setDY(arg.getDY());
	temp.setX1(arg.getX1()); temp.setX2(arg.getX2());
	temp.setY1(arg.getY1()); temp.setY2(arg.getY2());
	temp.setG(arg.getG()); temp.setB(arg.getB());
	//Execute GPU duals:
	cR2_EXP << <arg.getG(), arg.getB() >> > (temp.getD_ptr(), arg.getD_ptr(), arg.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 pow(const cR2 & base, double exp) {
	cR2 temp(base.getIndex1(), base.getIndex2());
	temp.setDX(base.getDX()); temp.setDY(base.getDY());
	temp.setX1(base.getX1()); temp.setX2(base.getX2());
	temp.setY1(base.getY1()); temp.setY2(base.getY2());
	temp.setG(base.getG()); temp.setB(base.getB());
	//Execute GPU duals:
	cR2_POW << <base.getG(), base.getB() >> > (temp.getD_ptr(), base.getD_ptr(), exp, base.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 pow(double base, const cR2 & exp) {
	cR2 temp(exp.getIndex1(), exp.getIndex2());
	temp.setDX(exp.getDX()); temp.setDY(exp.getDY());
	temp.setX1(exp.getX1()); temp.setX2(exp.getX2());
	temp.setY1(exp.getY1()); temp.setY2(exp.getY2());
	temp.setG(exp.getG()); temp.setB(exp.getB());
	//Execute GPU duals:
	cR2_POW << <exp.getG(), exp.getB() >> > (temp.getD_ptr(), base, exp.getD_ptr(), exp.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 pow(const cR2 & base, const cR2 & exp) {
	cR2 temp(base.getIndex1(), base.getIndex2());
	temp.setDX(base.getDX()); temp.setDY(base.getDY());
	temp.setX1(base.getX1()); temp.setX2(base.getX2());
	temp.setY1(base.getY1()); temp.setY2(base.getY2());
	temp.setG(base.getG()); temp.setB(base.getB());
	//Execute GPU duals:
	cR2_POW << <base.getG(), base.getB() >> > (temp.getD_ptr(), base.getD_ptr(), exp.getD_ptr(), base.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 asin(const cR2 & arg) {
	cR2 temp(arg.getIndex1(), arg.getIndex2());
	temp.setDX(arg.getDX()); temp.setDY(arg.getDY());
	temp.setX1(arg.getX1()); temp.setX2(arg.getX2());
	temp.setY1(arg.getY1()); temp.setY2(arg.getY2());
	temp.setG(arg.getG()); temp.setB(arg.getB());
	//Execute GPU duals:
	cR2_ASIN << <arg.getG(), arg.getB() >> > (temp.getD_ptr(), arg.getD_ptr(), arg.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 acos(const cR2 & arg) {
	cR2 temp(arg.getIndex1(), arg.getIndex2());
	temp.setDX(arg.getDX()); temp.setDY(arg.getDY());
	temp.setX1(arg.getX1()); temp.setX2(arg.getX2());
	temp.setY1(arg.getY1()); temp.setY2(arg.getY2());
	temp.setG(arg.getG()); temp.setB(arg.getB());
	//Execute GPU duals:
	cR2_ACOS << <arg.getG(), arg.getB() >> > (temp.getD_ptr(), arg.getD_ptr(), arg.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 atan(const cR2 & arg) {
	cR2 temp(arg.getIndex1(), arg.getIndex2());
	temp.setDX(arg.getDX()); temp.setDY(arg.getDY());
	temp.setX1(arg.getX1()); temp.setX2(arg.getX2());
	temp.setY1(arg.getY1()); temp.setY2(arg.getY2());
	temp.setG(arg.getG()); temp.setB(arg.getB());
	//Execute GPU duals:
	cR2_ATAN << <arg.getG(), arg.getB() >> > (temp.getD_ptr(), arg.getD_ptr(), arg.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 atan(const cR2 & y, const cR2 & x) {
	cR2 temp(x.getIndex1(), x.getIndex2());
	temp.setDX(x.getDX()); temp.setDY(x.getDY());
	temp.setX1(x.getX1()); temp.setX2(x.getX2());
	temp.setY1(x.getY1()); temp.setY2(x.getY2());
	temp.setG(x.getG()); temp.setB(x.getB());
	//Execute GPU duals:
	cR2_ATAN << <x.getG(), x.getB() >> > (temp.getD_ptr(), y.getD_ptr(), x.getD_ptr(), x.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 erf(const cR2 & arg) {
	cR2 temp(arg.getIndex1(), arg.getIndex2());
	temp.setDX(arg.getDX()); temp.setDY(arg.getDY());
	temp.setX1(arg.getX1()); temp.setX2(arg.getX2());
	temp.setY1(arg.getY1()); temp.setY2(arg.getY2());
	temp.setG(arg.getG()); temp.setB(arg.getB());
	//Execute GPU duals:
	cR2_ERF << <arg.getG(), arg.getB() >> > (temp.getD_ptr(), arg.getD_ptr(), arg.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 jn(int n, const cR2 & arg) {
	cR2 temp(arg.getIndex1(), arg.getIndex2());
	temp.setDX(arg.getDX()); temp.setDY(arg.getDY());
	temp.setX1(arg.getX1()); temp.setX2(arg.getX2());
	temp.setY1(arg.getY1()); temp.setY2(arg.getY2());
	temp.setG(arg.getG()); temp.setB(arg.getB());
	//Execute GPU duals:
	cR2_J_n << <arg.getG(), arg.getB() >> > (temp.getD_ptr(), n, arg.getD_ptr(), arg.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}

cR2 ln(const cR2 & arg) {
	cR2 temp(arg.getIndex1(), arg.getIndex2());
	temp.setDX(arg.getDX()); temp.setDY(arg.getDY());
	temp.setX1(arg.getX1()); temp.setX2(arg.getX2());
	temp.setY1(arg.getY1()); temp.setY2(arg.getY2());
	temp.setG(arg.getG()); temp.setB(arg.getB());
	//Execute GPU duals:
	cR2_ERF << <arg.getG(), arg.getB() >> > (temp.getD_ptr(), arg.getD_ptr(), arg.getSize());

	//Optional: COPY to host, find better place!
	//cudaMemcpy(h_yy, d_yy, m*n * sizeof(double), cudaMemcpyDeviceToHost);
	return temp;
}
