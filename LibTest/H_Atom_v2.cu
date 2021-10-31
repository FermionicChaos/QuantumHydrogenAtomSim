#include "H_Atom_v2.h"

#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include "stdud.h"
#include "extra.h"
using namespace std;
using namespace sf;

const double pi = 3.1415926535;

__global__ void Superposition(double* A, double* B, double* R, double* Theta,
	double* Psi, double* SVec,
	double tt, int I, int J, int K, int L) {

	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int k = 0;
	A[i] = 0.0; B[i] = 0.0;
	for (int n = 1; n < L + 1; ++n) {
		for (int l = 0; l < n; ++l) {
			for (int m = -l; m < (l + 1); ++m) {
				A[i] += cos(tt/double(n*n))*SVec[k] * Psi[i + I*J*k];
				B[i] += -sin(tt/double(n*n))*SVec[k] * Psi[i + I*J*k];
				k += 1;
			}
		}
	}//*/
	R[i] = sqrt((A[i] * A[i]) + (B[i] * B[i]));
	Theta[i] = atan2(B[i], A[i]);
}
//Evaluate for each eigenfunction.
// C|nlm><nlm|f> = 
__global__ void EFP(double* lambda_R, double* lambda_Theta,
	double* r, double* theta, double* R, double* Theta) {

	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int k = 0;
	lambda_R[i] = r[i] * R[i];
	lambda_Theta[i] = theta[i] - Theta[i];
}

__global__ void Superposition2(double* A, double* B, double* R, double* Theta,
	double* Psi, double* AVec, double* BVec,
	double tt, int I, int J, int K, int L) {
	//Includes complex eigen vectors.
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int k = 0;
	A[i] = 0.0; B[i] = 0.0;
	for (int n = 1; n < L + 1; ++n) {
		for (int l = 0; l < n; ++l) {
			for (int m = -l; m < (l + 1); ++m) {
				A[i] += cos(2.0*tt / double(n*n))*AVec[k] * Psi[i + I*J*k];
				B[i] += -sin(2.0*tt / double(n*n))*BVec[k] * Psi[i + I*J*k];
				k += 1;
			}
		}
	}//*/
	R[i] = sqrt((A[i] * A[i]) + (B[i] * B[i]));
	Theta[i] = atan2(B[i], A[i]);
}

__global__ void add(double* b, double* a, int I) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int j = blockDim.y * blockIdx.y + threadIdx.y;
	b[i + j*I] = a[i + j*I] * a[i + j*I];
}


H_Atom_v2::H_Atom_v2() {

}

H_Atom_v2::~H_Atom_v2() {
	//delete[] h_Ef;
	//delete[] h_vec;
	//delete[] h_pix;
	cudaFree(d_Ef);
	cudaFree(d_avec);
	cudaFree(d_bvec);
	cudaFree(d_pix);
	cudaFree(d_psia);
	cudaFree(d_psib);
	cudaFree(d_R);
	cudaFree(d_phase);
}

void H_Atom_v2::StateVec(double* psia, double* psib, cR2 & x, cR2 & y, int N) {
	tt = 0.0;
	h_avec = new double[N];
	h_bvec = new double[N];
	h_avec = psia; h_bvec = psib;
	dx = x.getDX(); dy = y.getDY();
	I = x.getIndex1(); J = y.getIndex2();
	Nh = N;
	K = N*(N + 1)*(2 * N + 1) / 6; //Develops energy eigenstate basis
	block = 1024;
	grid = I*J / 1024;
	alpha = x; phase = y;
	h_psia = new double[I*J];
	h_psib = new double[I*J];
	h_phase = new double[I*J];
	h_R = new double[I*J];
	cudaMalloc(&d_psia, I*J * sizeof(double));
	cudaMalloc(&d_psib, I*J * sizeof(double));
	cudaMalloc(&d_phase, I*J * sizeof(double));
	cudaMalloc(&d_R, I*J * sizeof(double));

	h_Ef = new double[I*J*K];
	cudaMalloc((void **)&d_Ef, I*J*K * sizeof(double));

	cR2 temp;
	int index = 0;
	cout << "<x,y,z|n,l,m>" << endl;
	for (int n = 1; n < N + 1; ++n) {
		for (int l = 0; l < n; ++l) {
			for (int m = -l; m < (l + 1); ++m) {
				cout << "|";
				cout << n << "," << l << "," << m;
				cout << ">" << endl;

				temp = Psi_nlm(n, l, m, x, y);
				temp.call2Host();
				for (int i = 0; i < I*J; ++i) {
					h_Ef[i + I*J*index] = temp.getH_ptr()[i];
				}

				cout << "Index number:  " << index << endl;
				index += 1;
			}
		}
	}
	cudaMemcpy(d_Ef, h_Ef, I*J*K * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_avec, h_avec, K * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_bvec, h_bvec, K * sizeof(double), cudaMemcpyHostToDevice);
	cout << "Data Generated!" << endl;
}
/*
What we need, State vector from hilbert space E states |n,l,m>
|Psi> = SUM(C_k * |n,l,m>);
The eigenstate vector length l;
The dimension of the screen vertical 720 x horizontal 1280.
|Psi> = C_0*|1,0,0> + C_1*|2,0,0> + C_3*|2,1,-1> + ...
Position projection
<r|Psi> = C_0*<r|1,0,0> + C_1*<r|2,0,0> + C_3*<r|2,1,-1> + ...
Time evolution:
<r|Psi> = C_0*exp(-t)*<r|1,0,0> + C_1*exp(-t/4)*<r|2,0,0> + C_3*exp(-t/4)*<r|2,1,-1> + ...
Step 1: Generate eigen functions
Step 2: Time evolve for small dt
Step 3: Convert data into pixel data.
Step 4: Draw.

A point in space will be a complex number for all x,y in the set. The phase will determine the
color of the vertex, and the magnitude will determine the alpha opacity, after which the data is normalized
to create the pixel data.

gggg
*/
void H_Atom_v2::genData(int RESX,int RESY, cR2& x, cR2& y, int N) {

	tt = 0.0;
	I = RESX; J = RESY;
	dx = x.getDX(); dy = y.getDY();
	Nh = N;
	K = N*(N + 1)*(2 * N + 1) / 6; //Develops energy eigenstate basis
	block = 512;
	grid = I*J / 512;
	//Yeah, I know, it's a lot of threads. fite me
	alpha = x; phase = y;
	h_psia = new double[RESX*RESY];
	h_psib = new double[RESX*RESY];
	h_phase = new double[RESX*RESY];
	h_R = new double[RESX*RESY];
	cudaMalloc(&d_psia, RESX*RESY * sizeof(double));
	cudaMalloc(&d_psib, RESX*RESY * sizeof(double));
	cudaMalloc(&d_phase, RESX*RESY * sizeof(double));
	cudaMalloc(&d_R, RESX*RESY * sizeof(double));

	h_Ef = new double[RESX*RESY*K];
	cudaMalloc((void **)&d_Ef, I*J*K * sizeof(double));
	h_avec = new double[K];
	cudaMalloc((void **)&d_avec, K * sizeof(double));
	h_bvec = new double[K];
	cudaMalloc((void **)&d_bvec, K * sizeof(double));

	cR2 temp;
	int index = 0;
	cout << "<x,y,z|n,l,m>" << endl;
	for (int n = 1; n < N + 1; ++n) {
		for (int l = 0; l < n; ++l) {
			for (int m = -l; m < (l + 1); ++m) {
				cout << "|";
				cout << n << "," << l << "," << m;
				cout << ">" << endl;

				temp = Psi_nlm(n, l, m, x, y);
				temp.call2Host();
				//cout << temp(0.0, 0.0) << endl;
				for (int i = 0; i < I*J; ++i) {
						h_Ef[i + I*J*index] = temp.getH_ptr()[i];
				}

				h_avec[index] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
				h_bvec[index] = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
				//h_avec[index] = 1.0;
				//h_bvec[index] = 0.0;
				cout << "Index number:  " << index << endl;
				index += 1;
			}
		}
	}
	//Update you fucking cock, nao
	cudaMemcpy(d_Ef, h_Ef, I*J*K*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_avec, h_avec, K*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_bvec, h_bvec, K * sizeof(double), cudaMemcpyHostToDevice);
	cout << "Data Generated!" << endl;
}

void H_Atom_v2::genIm(int H, int W) {
	Im = VertexArray(Points, H*W);
	h = H; w = W;
}

void H_Atom_v2::compute(float dt) {
	tt += double(dt);

	/*
	COCK << <grid, block>> > (d_psia, d_psib, d_R, d_phase, d_Ef, d_vec, tt, I, J, K, Nh);
	cudaMemcpy(h_psia, d_psia, I*J * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_psib, d_psib, I*J * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_phase, d_phase, I*J * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_R, d_R, I*J * sizeof(double), cudaMemcpyDeviceToHost);
	//*/

	Superposition2 << <grid, block >> > (d_psia, d_psib, alpha.getD_ptr(), phase.getD_ptr(), d_Ef, d_avec, d_bvec, tt, I, J, K, Nh);
	//cout << tt << endl;
	
	///*
	//cudaMemcpy(h_psia, d_psia, I*J * sizeof(double), cudaMemcpyDeviceToHost);
	//cudaMemcpy(h_psib, d_psib, I*J * sizeof(double), cudaMemcpyDeviceToHost);

	//Normalize arrays for conversion to pixel data.
	alpha.call2Host();
	phase.call2Host();
	alpha = normalize(alpha);
	alpha.call2Host();
	

	/*
	Uint8 r, g, b, a;
	double k = 255.0;
	for (int j = 0; j < J; ++j) {
		for (int i = 0; i < I; ++i) {
			r = 0; g = 0; b = 0; a = 0;
			if (cos(phase(i, j)) > 0.0) {
				r = static_cast<Uint8>(k*cos(phase(i, j)));
			}
			if (cos(phase(i, j) + 2 * pi / 3) > 0.0) {
				g = static_cast<Uint8>(k*cos(phase(i, j) + 2 * pi / 3));
			}
			if (cos(phase(i, j) - 2 * pi / 3) > 0.0) {
				b = static_cast<Uint8>(k*cos(phase(i, j) - 2 * pi / 3));
			}
			if (abs(alpha(i,j)) > 0.0) {
				a = static_cast<char>(k*alpha(i,j));
			}
			Im[i + j*w].position = Vector2f(i, j);
			Im[i + j*w].color = Color(r,g,b,a);
			Im[i + j*w].texCoords = Vector2f(i, j);
		}
	}//*/
	///*
	Uint8 r, g, b, a;
	double k = 255.0;
	for (int j = 0; j < J; ++j) {
		for (int i = 0; i < I; ++i) {
			r = 0; g = 0; b = 0; a = 0;
			r = static_cast<Uint8>(k*0.5*(cos(phase(i, j)) + 1));
			g = static_cast<Uint8>(k*0.5*(cos(phase(i, j) + 2 * pi / 3) + 1.0));
			b = static_cast<Uint8>(k*0.5*(cos(phase(i, j) - 2 * pi / 3) + 1.0));
			a = static_cast<char>(k*alpha(i, j));
			Im[i + j*w].position = Vector2f(i, j);
			Im[i + j*w].color = Color(r, g, b, a);
			Im[i + j*w].texCoords = Vector2f(i, j);
		}
	}
	//*/
}

void H_Atom_v2::draw(RenderWindow & window) {
	window.draw(Im);
}

void H_Atom_v2::ArrayTest() {
	double h_A[9];
	double *d_A;
	double h_B[9];
	double *d_B;
	cudaMalloc(&d_A, 9*sizeof(double));
	cudaMalloc(&d_B, 9*sizeof(double));
	for (int i = 0; i < 9; ++i) {
		h_A[i] = double(i);
		cout << h_A[i] << endl;
	}
	cudaMemcpy(d_A, h_A, 9 * sizeof(double), cudaMemcpyHostToDevice);
	dim3 g;
	g = dim3(3, 3);
	dim3 b;
	b = dim3(1, 1);
	add << <g, b >> > (d_B, d_A, 3);
	
	cudaMemcpy(h_B, d_B, 9 * sizeof(double), cudaMemcpyDeviceToHost);
	for (int j = 0; j < 9; ++j) {
		cout << h_A[j] << endl;
	}
	cudaFree(d_A);
	cudaFree(d_B);
}
