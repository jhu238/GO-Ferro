#pragma once
#include "matrix.h"
#include "global.h"
//#include <complex>
//#include <fftw3.h>
//#include "cufft.h"

class mathlib {
public:
	class global_parameters* pt_glb = &(global_parameters::glb);
public:
	/*mathlib(global_parameters*);*/
//#pragma acc routine(cufftPlan3d)(cufftExecD2Z)(cufftDestroy) nohost
//	void fourier_transform3D(matrix3d<double>&, matrix3d<std::complex<double>>&);
//#pragma acc routine(cufftPlan3d)(cufftExecD2Z)(cufftDestroy) nohost
//	void fourier_transform3D(matrix3d<std::complex<double>>&, matrix3d<double>&);

	void copy_to_device();

#pragma acc routine seq// nohost
	void transform_vector_glb2crt(const double&, const double&, const double&, \
		double&, double&, double&);
#pragma acc routine seq// nohost
	void transform_vector_crt2glb(const double&, const double&, const double&, \
		double&, double&, double&);
#pragma acc routine seq// nohost
	void transform_matrix_glb2crt(const double (&in)[3][3], double (&out)[3][3]);
#pragma acc routine seq// nohost
	void transform_matrix_crt2glb(const double (&in)[3][3], double (&out)[3][3]);
#pragma acc routine seq// nohost
	void strain2stress_glb(const double(&stiff)[6][6], const double(&strain)[6], double(&stress)[6]);
public:
	static mathlib mlb;
};
