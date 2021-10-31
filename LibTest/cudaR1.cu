#include "cudaR1.h"
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "stdud.h"


using namespace std;
/*
DEPENDENT ON CUDA 8.0!
*/

__global__ void cudaR1_Add(double* c, double* a, double* b) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	c[i] = a[i] + b[i];
}

__global__ void cudaR1_Sub(double* c, double* a, double* b) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	c[i] = a[i] - b[i];
}

__global__ void cudaR1_Mult(double* c, double* a, double* b) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	c[i] = a[i] * b[i];
}

__global__ void cudaR1_Divide(double* c, double* a, double* b) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (b[i] != 0.0) {
		c[i] = a[i] / b[i];
	}
	else {
		c[i] = 0.0;
	}
}

__global__ void cudaR1_Add_RHS(double* c, double* a, double rhs) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	c[i] = a[i] + rhs;
}

__global__ void cudaR1_Sub_RHS(double* c, double* a, double rhs) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	c[i] = a[i] - rhs;
}

__global__ void cudaR1_Mult_RHS(double* c, double* a, double rhs) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	c[i] = a[i] * rhs;
}

__global__ void cudaR1_Divide_RHS(double* c, double* a, double rhs) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	c[i] = a[i] / rhs;
}

__global__ void cudaR1_Add_LHS(double* c, double lhs, double* b) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	c[i] = lhs + b[i];
}

__global__ void cudaR1_Sub_LHS(double* c, double lhs, double* b) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	c[i] = lhs - b[i];
}

__global__ void cudaR1_Mult_LHS(double* c, double lhs, double* b) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	c[i] = lhs * b[i];
}

__global__ void cudaR1_Divide_LHS(double* c, double lhs, double* b) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (b[i] != 0.0) {
		c[i] = lhs / b[i];
	}
	else {
		c[i] = 0.0;
	}
}

cudaR1::cudaR1() {
	ptr = nullptr;
	x1 = 0; x2 = 0, dx = 0; m = 0;
	grid = 0; block = 0;
}

cudaR1::~cudaR1() {
	//Destructor
	//cout << "MEMORY De-Allocated: ADDRESS -> " << ptr << endl;
	delete[] ptr;
}

cudaR1::cudaR1(const cudaR1& rhs) {
	//Copy Constructor SEPERATE POINTER VALUES FOR DIFFERENT OBJECTS
	//cout << "Copy Constructor" << endl;
	m = rhs.getSize(); dx = rhs.getDX();
	x1 = rhs.getX1(); x2 = rhs.getX2();
	grid = rhs.grid; block = rhs.block;
	double* BR = rhs.ptr;
	ptr = new double[m];
	//cout << "COPIED From: ADDRESS -> " << rhs.getptr() << endl;
	//cout << "COPIED To: -> " << ptr << endl;
	for (int i = 0; i < m; ++i) {
		ptr[i] = BR[i];
	}
	//cout << "END: COPY Constructor Complete" << endl;
}

cudaR1& cudaR1::operator=(const cudaR1& rhs) {
	//cout << "Operation: ASSIGNMENT" << endl;
	//cout << "LEFT ADDRESS -> " << ptr << endl;
	//cout << "RIGHT ADDRESS -> " << rhs.getptr() << endl;
	//cout << rhs << endl;
	if ((ptr != rhs.getBase())) {
		this->m = rhs.getSize(); this->dx = rhs.getDX();
		this->x1 = rhs.getX1(); this->x2 = rhs.getX2();
		this->grid = rhs.grid; this->block = rhs.block;
		delete[] this->ptr;
		ptr = new double[m];
		for (int i = 0; i < m; ++i) {
			ptr[i] = rhs.getBase()[i];
		}
		this->ptr = ptr;
		return *this;
	}
	else {
		return *this;
	}
}

cudaR1::cudaR1(cudaR1&& inp) {
	ptr = inp.ptr;
	inp.ptr = nullptr;
	m = inp.m; dx = inp.dx;
	x1 = inp.x1; x2 = inp.x2;
	grid = inp.grid; block = inp.block;
}

cudaR1& cudaR1::operator=(cudaR1&& rhs) {
	delete[] ptr;
	ptr = rhs.ptr;
	m = rhs.m; dx = rhs.dx;
	x1 = rhs.x1; x2 = rhs.x2;
	grid = rhs.grid; block = rhs.block;
	rhs.ptr = nullptr;
	return *this;
}

cudaR1::cudaR1(int size) {
	m = size; x1 = 0; x2 = 0; dx = 0;
	ptr = new double[m];
}

cudaR1::cudaR1(double a, double b, int res) {
	//Large latency, initiate for domain.
	//Initiate as soon as possible for other operations.
	if ((res%256) == 0) {
		x1 = a; x2 = b; m = res;
		grid = (m / 256); block = 256;
		dx = (b - a) / (double(res));
		ptr = new double[res];
		for (int i = 0; i < m; ++i) {
			ptr[i] = x1 + double(i)*dx;
		}
		x2 = ptr[m - 1];
	}
	else {
		cout << "Must be a multiple of 256!" << endl;
	}
}

double cudaR1::operator()(double p0) {
	double val, p;
	int i;
	p = (p0 - x1) / dx;
	i = int(p);
	if ((i > 0) && (i < m - 1)) {
		val = ptr[i];
		return val;
	}
	else {
		return val = 0;
	}
	return val;
}

double cudaR1::operator()(int i) {
	double val;
	if ((i > 0) && (i < m - 1)) {
		val = ptr[i];
		return val;
	}
	else {
		return val = 0;
	}
}

cudaR1 cudaR1::operator+(const cudaR1 & rhs) {
	cudaR1 temp;
	if ((x1 = rhs.x1) && (x2 = rhs.x2) && (m = rhs.m) && (dx = rhs.dx)) {
		double* h_xx1 = ptr;
		double* h_xx2 = rhs.getBase();
		double* h_yy = new double[rhs.getSize()];
		temp.setBase(h_yy); temp.setDX(dx); temp.setSize(m);
		temp.setX1(x1); temp.setX2(x2);
		//GPU duals.
		double *d_xx1, *d_xx2, *d_yy;
		cudaError_t cudaStatus;
		//Check:
		cudaStatus = cudaSetDevice(0);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
			goto Error;
		}

		//Reserve:
		cudaStatus = cudaMalloc(&d_xx1, rhs.getSize() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		cudaStatus = cudaMalloc(&d_xx2, rhs.getSize() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		cudaStatus = cudaMalloc(&d_yy, rhs.getSize() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		//Copy2d:
		cudaStatus = cudaMemcpy(d_xx1, h_xx1, rhs.getSize() * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}
		cudaStatus = cudaMemcpy(d_xx2, h_xx2, rhs.getSize() * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}

		// grid has blocks, and a block has multi threads.
		//Execute: F : R1 -> R1 { xx1 + xx2 = yy }
		cudaR1_Add<<<grid, block>>>(d_yy, d_xx1, d_xx2);

		//Sync:	
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Add launch failed: %s\n", cudaGetErrorString(cudaStatus));
			goto Error;
		}

		// cudaDeviceSynchronize waits for the kernel to finish, and returns
		// any errors encountered during the launch.
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
			goto Error;
		}
		//Return2h:
		cudaStatus = cudaMemcpy(h_yy, d_yy, rhs.getSize() * sizeof(double), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}
		//Transfer back:
		//temp.setBase(h_yy);
		//ReleaseMEM:
		cudaFree(d_xx1); cudaFree(d_xx2); cudaFree(d_yy);
		cout << "Add Success!: " << false << endl;
		return temp;

	Error:
		cudaFree(d_xx1);
		cudaFree(d_xx2);
		cudaFree(d_yy);
		cout << "Error!" << true << endl;
		return temp;
	}
	else {
		cout << "Error: Dimensions must agree!" << endl;
		return temp;
	}
}

cudaR1 cudaR1::operator-(const cudaR1 & rhs) {
	cudaR1 temp;
	if ((x1 = rhs.x1) && (x2 = rhs.x2) && (m = rhs.m) && (dx = rhs.dx)) {
		double* h_xx1 = ptr;
		double* h_xx2 = rhs.getBase();
		double* h_yy = new double[rhs.getSize()];
		temp.setBase(h_yy); temp.setDX(dx); temp.setSize(m);
		temp.setX1(x1); temp.setX2(x2);
		//GPU duals.
		double *d_xx1, *d_xx2, *d_yy;
		cudaError_t cudaStatus;
		//Check:
		cudaStatus = cudaSetDevice(0);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
			goto Error;
		}

		//Reserve:
		cudaStatus = cudaMalloc(&d_xx1, rhs.getSize() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		cudaStatus = cudaMalloc(&d_xx2, rhs.getSize() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		cudaStatus = cudaMalloc(&d_yy, rhs.getSize() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		//Copy2d:
		cudaStatus = cudaMemcpy(d_xx1, h_xx1, rhs.getSize() * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}
		cudaStatus = cudaMemcpy(d_xx2, h_xx2, rhs.getSize() * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}

		// grid has blocks, and a block has multi threads.
		//Execute: F : R1 -> R1 { xx1 + xx2 = yy }
		cudaR1_Sub<<<grid, block >>>(d_yy, d_xx1, d_xx2);

		//Sync:	
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Add launch failed: %s\n", cudaGetErrorString(cudaStatus));
			goto Error;
		}

		// cudaDeviceSynchronize waits for the kernel to finish, and returns
		// any errors encountered during the launch.
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
			goto Error;
		}
		//Return2h:
		cudaStatus = cudaMemcpy(h_yy, d_yy, rhs.getSize() * sizeof(double), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}
		//Transfer back:
		//temp.setBase(h_yy);
		//ReleaseMEM:
		cudaFree(d_xx1); cudaFree(d_xx2); cudaFree(d_yy);
		cout << "Add Success!: " << false << endl;
		return temp;

	Error:
		cudaFree(d_xx1);
		cudaFree(d_xx2);
		cudaFree(d_yy);
		cout << "Error!" << true << endl;
		return temp;
	}
	else {
		cout << "Error: Dimensions must agree!" << endl;
		return temp;
	}
}

cudaR1 cudaR1::operator*(const cudaR1 & rhs) {
	cudaR1 temp;
	if ((x1 = rhs.x1) && (x2 = rhs.x2) && (m = rhs.m) && (dx = rhs.dx)) {
		double* h_xx1 = ptr;
		double* h_xx2 = rhs.getBase();
		double* h_yy = new double[rhs.getSize()];
		temp.setBase(h_yy); temp.setDX(dx); temp.setSize(m);
		temp.setX1(x1); temp.setX2(x2);
		//GPU duals.
		double *d_xx1, *d_xx2, *d_yy;
		cudaError_t cudaStatus;
		//Check:
		cudaStatus = cudaSetDevice(0);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
			goto Error;
		}

		//Reserve:
		cudaStatus = cudaMalloc(&d_xx1, rhs.getSize() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		cudaStatus = cudaMalloc(&d_xx2, rhs.getSize() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		cudaStatus = cudaMalloc(&d_yy, rhs.getSize() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		//Copy2d:
		cudaStatus = cudaMemcpy(d_xx1, h_xx1, rhs.getSize() * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}
		cudaStatus = cudaMemcpy(d_xx2, h_xx2, rhs.getSize() * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}

		// grid has blocks, and a block has multi threads.
		//Execute: F : R1 -> R1 { xx1 + xx2 = yy }
		cudaR1_Mult<<<grid, block >>>(d_yy, d_xx1, d_xx2);

		//Sync:	
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Add launch failed: %s\n", cudaGetErrorString(cudaStatus));
			goto Error;
		}

		// cudaDeviceSynchronize waits for the kernel to finish, and returns
		// any errors encountered during the launch.
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
			goto Error;
		}
		//Return2h:
		cudaStatus = cudaMemcpy(h_yy, d_yy, rhs.getSize() * sizeof(double), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}
		//Transfer back:
		//temp.setBase(h_yy);
		//ReleaseMEM:
		cudaFree(d_xx1); cudaFree(d_xx2); cudaFree(d_yy);
		cout << "Add Success!: " << false << endl;
		return temp;

	Error:
		cudaFree(d_xx1);
		cudaFree(d_xx2);
		cudaFree(d_yy);
		cout << "Error!" << true << endl;
		return temp;
	}
	else {
		cout << "Error: Dimensions must agree!" << endl;
		return temp;
	}
}

cudaR1 cudaR1::operator/(const cudaR1 & rhs) {
	cudaR1 temp;
	if ((x1 = rhs.x1) && (x2 = rhs.x2) && (m = rhs.m) && (dx = rhs.dx)) {
		double* h_xx1 = ptr;
		double* h_xx2 = rhs.getBase();
		double* h_yy = new double[rhs.getSize()];
		temp.setBase(h_yy); temp.setDX(dx); temp.setSize(m);
		temp.setX1(x1); temp.setX2(x2);
		//GPU duals.
		double *d_xx1, *d_xx2, *d_yy;
		cudaError_t cudaStatus;
		//Check:
		cudaStatus = cudaSetDevice(0);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
			goto Error;
		}

		//Reserve:
		cudaStatus = cudaMalloc(&d_xx1, rhs.getSize() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		cudaStatus = cudaMalloc(&d_xx2, rhs.getSize() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		cudaStatus = cudaMalloc(&d_yy, rhs.getSize() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		//Copy2d:
		cudaStatus = cudaMemcpy(d_xx1, h_xx1, rhs.getSize() * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}
		cudaStatus = cudaMemcpy(d_xx2, h_xx2, rhs.getSize() * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}

		// grid has blocks, and a block has multi threads.
		//Execute: F : R1 -> R1 { xx1 + xx2 = yy }
		cudaR1_Divide<<<grid, block >>>(d_yy, d_xx1, d_xx2);

		//Sync:	
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Add launch failed: %s\n", cudaGetErrorString(cudaStatus));
			goto Error;
		}

		// cudaDeviceSynchronize waits for the kernel to finish, and returns
		// any errors encountered during the launch.
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
			goto Error;
		}
		//Return2h:
		cudaStatus = cudaMemcpy(h_yy, d_yy, rhs.getSize() * sizeof(double), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}
		//Transfer back:
		//temp.setBase(h_yy);
		//ReleaseMEM:
		cudaFree(d_xx1); cudaFree(d_xx2); cudaFree(d_yy);
		cout << "Add Success!: " << false << endl;
		return temp;

	Error:
		cudaFree(d_xx1);
		cudaFree(d_xx2);
		cudaFree(d_yy);
		cout << "Error!" << true << endl;
		return temp;
	}
	else {
		cout << "Error: Dimensions must agree!" << endl;
		return temp;
	}
}

cudaR1 cudaR1::operator+(double rhs) {
	cudaR1 temp;
	double* h_xx1 = ptr;
	//double* hRHS = &rhs;
	double* h_yy = new double[m];
	temp.setBase(h_yy); temp.setDX(dx); temp.setSize(m);
	temp.setX1(x1); temp.setX2(x2);
	//GPU duals.
	double *d_xx1, *d_yy;
	cudaError_t cudaStatus;
	//Check:
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	//Reserve:
	cudaStatus = cudaMalloc(&d_xx1, m * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	/*cudaStatus = cudaMalloc(&d_rhs, sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}*/
	cudaStatus = cudaMalloc(&d_yy, m * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	//Copy2d:
	cudaStatus = cudaMemcpy(d_xx1, h_xx1, m * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	/*cudaStatus = cudaMemcpy(d_rhs, &rhs, sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}*/

	// grid has blocks, and a block has multi threads.
	//Execute: F : R1 -> R1 { xx1 + xx2 = yy }
	cudaR1_Add_RHS<<<grid, block>>>(d_yy, d_xx1, rhs);

	//Sync:	
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Add launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		goto Error;
	}
	//Return2h:
	cudaStatus = cudaMemcpy(h_yy, d_yy, m * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	//Transfer back:
	//temp.setBase(h_yy);
	//ReleaseMEM:
	cudaFree(d_xx1); cudaFree(d_yy);
	cout << "Add Success!: " << false << endl;
	return temp;

Error:
	cudaFree(d_xx1);
	//cudaFree(d_xx2);
	cudaFree(d_yy);
	cout << "Error!" << true << endl;
	return temp;
	return cudaR1();
}

cudaR1 cudaR1::operator-(double rhs) {
	cudaR1 temp;
	double* h_xx1 = ptr;
	//double* hRHS = &rhs;
	double* h_yy = new double[m];
	temp.setBase(h_yy); temp.setDX(dx); temp.setSize(m);
	temp.setX1(x1); temp.setX2(x2);
	//GPU duals.
	double *d_xx1, *d_yy;
	cudaError_t cudaStatus;
	//Check:
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	//Reserve:
	cudaStatus = cudaMalloc(&d_xx1, m * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	/*cudaStatus = cudaMalloc(&d_rhs, sizeof(double));
	if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "cudaMalloc failed!");
	goto Error;
	}*/
	cudaStatus = cudaMalloc(&d_yy, m * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	//Copy2d:
	cudaStatus = cudaMemcpy(d_xx1, h_xx1, m * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	/*cudaStatus = cudaMemcpy(d_rhs, &rhs, sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "cudaMemcpy failed!");
	goto Error;
	}*/

	// grid has blocks, and a block has multi threads.
	//Execute: F : R1 -> R1 { xx1 + xx2 = yy }
	cudaR1_Sub_RHS<<<grid, block >>>(d_yy, d_xx1, rhs);

	//Sync:	
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Add launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		goto Error;
	}
	//Return2h:
	cudaStatus = cudaMemcpy(h_yy, d_yy, m * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	//Transfer back:
	//temp.setBase(h_yy);
	//ReleaseMEM:
	cudaFree(d_xx1); cudaFree(d_yy);
	cout << "Add Success!: " << false << endl;
	return temp;

Error:
	cudaFree(d_xx1);
	//cudaFree(d_xx2);
	cudaFree(d_yy);
	cout << "Error!" << true << endl;
	return temp;
	return cudaR1();
}

cudaR1 cudaR1::operator*(double rhs) {
	cudaR1 temp;
	double* h_xx1 = ptr;
	//double* hRHS = &rhs;
	double* h_yy = new double[m];
	temp.setBase(h_yy); temp.setDX(dx); temp.setSize(m);
	temp.setX1(x1); temp.setX2(x2);
	//GPU duals.
	double *d_xx1, *d_yy;
	cudaError_t cudaStatus;
	//Check:
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	//Reserve:
	cudaStatus = cudaMalloc(&d_xx1, m * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	/*cudaStatus = cudaMalloc(&d_rhs, sizeof(double));
	if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "cudaMalloc failed!");
	goto Error;
	}*/
	cudaStatus = cudaMalloc(&d_yy, m * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	//Copy2d:
	cudaStatus = cudaMemcpy(d_xx1, h_xx1, m * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	/*cudaStatus = cudaMemcpy(d_rhs, &rhs, sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "cudaMemcpy failed!");
	goto Error;
	}*/

	// grid has blocks, and a block has multi threads.
	//Execute: F : R1 -> R1 { xx1 + xx2 = yy }
	cudaR1_Mult_RHS<<<grid, block>>>(d_yy, d_xx1, rhs);

	//Sync:	
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Add launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		goto Error;
	}
	//Return2h:
	cudaStatus = cudaMemcpy(h_yy, d_yy, m * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	//Transfer back:
	//temp.setBase(h_yy);
	//ReleaseMEM:
	cudaFree(d_xx1); cudaFree(d_yy);
	cout << "Add Success!: " << false << endl;
	return temp;

Error:
	cudaFree(d_xx1);
	//cudaFree(d_xx2);
	cudaFree(d_yy);
	cout << "Error!" << true << endl;
	return temp;
	return cudaR1();
}

cudaR1 cudaR1::operator/(double rhs) {
	cudaR1 temp;
	if (rhs != 0.0) {
		double* h_xx1 = ptr;
		//double* hRHS = &rhs;
		double* h_yy = new double[m];
		temp.setBase(h_yy); temp.setDX(dx); temp.setSize(m);
		temp.setX1(x1); temp.setX2(x2);
		//GPU duals.
		double *d_xx1, *d_yy;
		cudaError_t cudaStatus;
		//Check:
		cudaStatus = cudaSetDevice(0);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
			goto Error;
		}

		//Reserve:
		cudaStatus = cudaMalloc(&d_xx1, m * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		/*cudaStatus = cudaMalloc(&d_rhs, sizeof(double));
		if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
		}*/
		cudaStatus = cudaMalloc(&d_yy, m * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		//Copy2d:
		cudaStatus = cudaMemcpy(d_xx1, h_xx1, m * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}
		/*cudaStatus = cudaMemcpy(d_rhs, &rhs, sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
		}*/

		// grid has blocks, and a block has multi threads.
		//Execute: F : R1 -> R1 { xx1 + xx2 = yy }
		cudaR1_Divide_RHS<<<grid, block>>>(d_yy, d_xx1, rhs);

		//Sync:	
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Add launch failed: %s\n", cudaGetErrorString(cudaStatus));
			goto Error;
		}

		// cudaDeviceSynchronize waits for the kernel to finish, and returns
		// any errors encountered during the launch.
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
			goto Error;
		}
		//Return2h:
		cudaStatus = cudaMemcpy(h_yy, d_yy, m * sizeof(double), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}
		//Transfer back:
		//temp.setBase(h_yy);
		//ReleaseMEM:
		cudaFree(d_xx1); cudaFree(d_yy);
		cout << "Add Success!: " << false << endl;
		return temp;

	Error:
		cudaFree(d_xx1);
		//cudaFree(d_xx2);
		cudaFree(d_yy);
		cout << "Error!" << true << endl;
		return temp;
	}
	else {
		cout << "Error: Division by zero is not allowed!" << endl;
		return temp;
	}
}

void cudaR1::InitGB() {
	//Initialize grid/block incase forgotten.
	block = 256;
	grid = m/256;
};

cudaR1 operator+(double lhs, const cudaR1& rhs) {
	cudaR1 temp;
	//double* h_xx1 = ptr;
	double* h_xx2 = rhs.ptr;
	double* h_yy = new double[rhs.m];
	temp.setBase(h_yy); temp.setDX(rhs.dx); temp.setSize(rhs.m);
	temp.setX1(rhs.x1); temp.setX2(rhs.x2);
	//GPU duals.
	double *d_xx2, *d_yy;
	cudaError_t cudaStatus;
	//Check:
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	//Reserve:
	/*cudaStatus = cudaMalloc(&d_xx1, rhs.getSize() * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}*/
	cudaStatus = cudaMalloc(&d_xx2, rhs.m * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	cudaStatus = cudaMalloc(&d_yy, rhs.m * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	//Copy2d:
	/*cudaStatus = cudaMemcpy(d_xx1, h_xx1, rhs.getSize() * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}*/
	cudaStatus = cudaMemcpy(d_xx2, h_xx2, rhs.m * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

	// grid has blocks, and a block has multi threads.
	//Execute: F : R1 -> R1 { xx1 + xx2 = yy }
	cudaR1_Add_LHS<<<rhs.grid, rhs.block>>>(d_yy, lhs, d_xx2);

	//Sync:	
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Add launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		goto Error;
	}
	//Return2h:
	cudaStatus = cudaMemcpy(h_yy, d_yy, rhs.getSize() * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	//Transfer back:
	//temp.setBase(h_yy);
	//ReleaseMEM:
	cudaFree(d_xx2); cudaFree(d_yy);
	cout << "Add Success!: " << false << endl;
	return temp;

Error:
	//cudaFree(d_xx1);
	cudaFree(d_xx2);
	cudaFree(d_yy);
	cout << "Error!" << true << endl;
	return temp;
}

cudaR1 operator-(double lhs, const cudaR1 & rhs) {
	cudaR1 temp;
	//double* h_xx1 = ptr;
	double* h_xx2 = rhs.ptr;
	double* h_yy = new double[rhs.m];
	temp.setBase(h_yy); temp.setDX(rhs.dx); temp.setSize(rhs.m);
	temp.setX1(rhs.x1); temp.setX2(rhs.x2);
	//GPU duals.
	double *d_xx2, *d_yy;
	cudaError_t cudaStatus;
	//Check:
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	//Reserve:
	/*cudaStatus = cudaMalloc(&d_xx1, rhs.getSize() * sizeof(double));
	if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "cudaMalloc failed!");
	goto Error;
	}*/
	cudaStatus = cudaMalloc(&d_xx2, rhs.m * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	cudaStatus = cudaMalloc(&d_yy, rhs.m * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	//Copy2d:
	/*cudaStatus = cudaMemcpy(d_xx1, h_xx1, rhs.getSize() * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "cudaMemcpy failed!");
	goto Error;
	}*/
	cudaStatus = cudaMemcpy(d_xx2, h_xx2, rhs.m * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

	// grid has blocks, and a block has multi threads.
	//Execute: F : R1 -> R1 { xx1 + xx2 = yy }
	cudaR1_Sub_LHS<<<rhs.grid, rhs.block>>>(d_yy, lhs, d_xx2);

	//Sync:	
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Add launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		goto Error;
	}
	//Return2h:
	cudaStatus = cudaMemcpy(h_yy, d_yy, rhs.getSize() * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	//Transfer back:
	//temp.setBase(h_yy);
	//ReleaseMEM:
	cudaFree(d_xx2); cudaFree(d_yy);
	cout << "Add Success!: " << false << endl;
	return temp;

Error:
	//cudaFree(d_xx1);
	cudaFree(d_xx2);
	cudaFree(d_yy);
	cout << "Error!" << true << endl;
	return temp;
}

cudaR1 operator*(double lhs, const cudaR1 & rhs) {
	cudaR1 temp;
	//double* h_xx1 = ptr;
	double* h_xx2 = rhs.ptr;
	double* h_yy = new double[rhs.m];
	temp.setBase(h_yy); temp.setDX(rhs.dx); temp.setSize(rhs.m);
	temp.setX1(rhs.x1); temp.setX2(rhs.x2);
	//GPU duals.
	double *d_xx2, *d_yy;
	cudaError_t cudaStatus;
	//Check:
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	//Reserve:
	/*cudaStatus = cudaMalloc(&d_xx1, rhs.getSize() * sizeof(double));
	if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "cudaMalloc failed!");
	goto Error;
	}*/
	cudaStatus = cudaMalloc(&d_xx2, rhs.m * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	cudaStatus = cudaMalloc(&d_yy, rhs.m * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	//Copy2d:
	/*cudaStatus = cudaMemcpy(d_xx1, h_xx1, rhs.getSize() * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "cudaMemcpy failed!");
	goto Error;
	}*/
	cudaStatus = cudaMemcpy(d_xx2, h_xx2, rhs.m * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

	// grid has blocks, and a block has multi threads.
	//Execute: F : R1 -> R1 { xx1 + xx2 = yy }
	cudaR1_Mult_LHS<<<rhs.grid, rhs.block>>>(d_yy, lhs, d_xx2);

	//Sync:	
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Add launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		goto Error;
	}
	//Return2h:
	cudaStatus = cudaMemcpy(h_yy, d_yy, rhs.getSize() * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	//Transfer back:
	//temp.setBase(h_yy);
	//ReleaseMEM:
	cudaFree(d_xx2); cudaFree(d_yy);
	cout << "Add Success!: " << false << endl;
	return temp;

Error:
	//cudaFree(d_xx1);
	cudaFree(d_xx2);
	cudaFree(d_yy);
	cout << "Error!" << true << endl;
	return temp;
}

cudaR1 operator/(double lhs, const cudaR1 & rhs) {
	cudaR1 temp;
	//double* h_xx1 = ptr;
	double* h_xx2 = rhs.ptr;
	double* h_yy = new double[rhs.m];
	temp.setBase(h_yy); temp.setDX(rhs.dx); temp.setSize(rhs.m);
	temp.setX1(rhs.x1); temp.setX2(rhs.x2);
	//GPU duals.
	double *d_xx2, *d_yy;
	cudaError_t cudaStatus;
	//Check:
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		goto Error;
	}

	//Reserve:
	/*cudaStatus = cudaMalloc(&d_xx1, rhs.getSize() * sizeof(double));
	if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "cudaMalloc failed!");
	goto Error;
	}*/
	cudaStatus = cudaMalloc(&d_xx2, rhs.m * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	cudaStatus = cudaMalloc(&d_yy, rhs.m * sizeof(double));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}
	//Copy2d:
	/*cudaStatus = cudaMemcpy(d_xx1, h_xx1, rhs.getSize() * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
	fprintf(stderr, "cudaMemcpy failed!");
	goto Error;
	}*/
	cudaStatus = cudaMemcpy(d_xx2, h_xx2, rhs.m * sizeof(double), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

	// grid has blocks, and a block has multi threads.
	//Execute: F : R1 -> R1 { xx1 + xx2 = yy }
	cudaR1_Divide_LHS << <rhs.grid, rhs.block >> >(d_yy, lhs, d_xx2);

	//Sync:	
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "Add launch failed: %s\n", cudaGetErrorString(cudaStatus));
		goto Error;
	}

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
		goto Error;
	}
	//Return2h:
	cudaStatus = cudaMemcpy(h_yy, d_yy, rhs.getSize() * sizeof(double), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}
	//Transfer back:
	//temp.setBase(h_yy);
	//ReleaseMEM:
	cudaFree(d_xx2); cudaFree(d_yy);
	cout << "Add Success!: " << false << endl;
	return temp;

Error:
	//cudaFree(d_xx1);
	cudaFree(d_xx2);
	cudaFree(d_yy);
	cout << "Error!" << true << endl;
	return temp;
}

/*
cudaR1 d_dx(const cudaR1 & f, int ver) {
	cudaR1 temp;
	if (ver == 1) {
		double* h_ff = f.ptr;
		//double* h_xx2 = rhs.getBase();
		double* h_yy = new double[f.m];
		temp.setBase(h_yy); temp.setDX(f.dx); temp.setSize(f.m);
		temp.setX1(f.x1); temp.setX2(f.x2);
		//GPU duals.
		double *d_ff, *d_yy;
		cudaError_t cudaStatus;
		//Check:
		cudaStatus = cudaSetDevice(0);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
			goto Error;
		}

		//Reserve:
		cudaStatus = cudaMalloc(&d_ff, f.m* sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		cudaStatus = cudaMalloc(&d_yy, f.m * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			goto Error;
		}
		//Copy2d:
		cudaStatus = cudaMemcpy(d_ff, h_ff, f.m * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}

		// grid has blocks, and a block has multi threads.
		//Execute: F : R1 -> R1 { xx1 + xx2 = yy }
		cudaR1_d_dx_Ver1<<<grid, block>>>(d_yy, d_ff, f.dx, );

		//Sync:	
		cudaStatus = cudaGetLastError();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "Add launch failed: %s\n", cudaGetErrorString(cudaStatus));
			goto Error;
		}

		// cudaDeviceSynchronize waits for the kernel to finish, and returns
		// any errors encountered during the launch.
		cudaStatus = cudaDeviceSynchronize();
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
			goto Error;
		}
		//Return2h:
		cudaStatus = cudaMemcpy(h_yy, d_yy, f.m * sizeof(double), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			goto Error;
		}
		//Transfer back:
		//temp.setBase(h_yy);
		//ReleaseMEM:
		cudaFree(d_ff); cudaFree(d_yy);
		cout << "Add Success!: " << false << endl;
		return temp;

	Error:
		cudaFree(d_ff);
		cudaFree(d_yy);
		cout << "Error!" << true << endl;
		return temp;
	}
	else if (ver == 2) {

	}
	else if (ver == 3) {

	}
	else {
		cout << "Invalid!" << endl;
	}
}*/
