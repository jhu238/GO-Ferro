#include "mathlib.h"

//mathlib::mathlib(global_parameters* pt_glb_in)
//	:pt_glb(pt_glb_in) {}
//#pragma acc routine(cufftPlan3d)(cufftExecD2Z)(cufftDestroy) nohost
//void mathlib::fourier_transform3D(matrix3d<double>& in, matrix3d<std::complex<double>>& out) {
//	//fftw_plan plan_forward;
//	//cufftHandle plan;
//
//	//fftw_complex* out_fftw;
//	//out_fftw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * out.dimz * out.dimy * out.dimx);
//	//in.transformTo1D();
//
//	//cuFFT
//	cufftHandle plan;
//	cufftPlan3d(&plan, in.dimx, in.dimy, in.dimz, CUFFT_D2Z);
//	cufftExecD2Z(plan, in.matrix, (cufftDoubleComplex*)out.matrix);
//	cufftDestroy(plan);
//
//	/*FFTW for serial*/
//	//plan_forward = fftw_plan_dft_r2c_3d(in.dimx, in.dimy, in.dimz, in.matrix, reinterpret_cast<fftw_complex*>(out.matrix), FFTW_ESTIMATE);
//	//fftw_execute(plan_forward);
//	////out.transformTo3D();
//	//fftw_destroy_plan(plan_forward);
//	//fftw_free(out_fftw);
//}
//#pragma acc routine(cufftPlan3d)(cufftExecD2Z)(cufftDestroy) nohost
//void mathlib::fourier_transform3D(matrix3d<std::complex<double>>& in, matrix3d<double>& out) {
//	//fftw_plan plan_backward;
//	//in.transformTo1D();
//	//cuFFT
//	cufftHandle plan;
//	cufftPlan3d(&plan, out.dimx, out.dimy, out.dimz, CUFFT_Z2D);
//	cufftExecZ2D(plan, (cufftDoubleComplex*)in.matrix, out.matrix);
//	cufftDestroy(plan);
//
//
//	/*FFTW for serial*/
//	//plan_backward = fftw_plan_dft_c2r_3d(out.dimx, out.dimy, out.dimz, reinterpret_cast<fftw_complex*>(in.matrix), out.matrix, FFTW_ESTIMATE);
//	//fftw_execute(plan_backward);
//	////out.transformTo3D();
//	//fftw_destroy_plan(plan_backward);
//}

#pragma acc routine seq// nohost
void mathlib::transform_vector_glb2crt(const double& vx_in, const double& vy_in, const double& vz_in, \
	double& vx_out, double& vy_out, double& vz_out) {

	if (pt_glb->if_OrientRotate == true) {
		vx_out = pt_glb->Tg2c[0][0] * vx_in + pt_glb->Tg2c[0][1] * vy_in + pt_glb->Tg2c[0][2] * vz_in;
		vy_out = pt_glb->Tg2c[1][0] * vx_in + pt_glb->Tg2c[1][1] * vy_in + pt_glb->Tg2c[1][2] * vz_in;
		vz_out = pt_glb->Tg2c[2][0] * vx_in + pt_glb->Tg2c[2][1] * vy_in + pt_glb->Tg2c[2][2] * vz_in;
	}
	else if (pt_glb->if_OrientRotate == false) {
		vx_out = vx_in;
		vy_out = vy_in;
		vz_out = vz_in;
	}
}

#pragma acc routine seq// nohost
void mathlib::transform_vector_crt2glb(const double& vx_in, const double& vy_in, const double& vz_in, \
	double& vx_out, double& vy_out, double& vz_out) {

	if (pt_glb->if_OrientRotate == true) {
		vx_out = pt_glb->Tc2g[0][0] * vx_in + pt_glb->Tc2g[0][1] * vy_in + pt_glb->Tc2g[0][2] * vz_in;
		vy_out = pt_glb->Tc2g[1][0] * vx_in + pt_glb->Tc2g[1][1] * vy_in + pt_glb->Tc2g[1][2] * vz_in;
		vz_out = pt_glb->Tc2g[2][0] * vx_in + pt_glb->Tc2g[2][1] * vy_in + pt_glb->Tc2g[2][2] * vz_in;
	}
	else if (pt_glb->if_OrientRotate == false) {
		vx_out = vx_in;
		vy_out = vy_in;
		vz_out = vz_in;
	}
}

#pragma acc routine seq// nohost
void mathlib::transform_matrix_glb2crt(const double(&in)[3][3], double(&out)[3][3]) {
	if (pt_glb->if_OrientRotate == true) {
		for (long int i = 0; i < 3; i++)
		{
			for (long int j = 0; j < 3; j++)
			{
				out[i][j] = 0.;
				for (long int k = 0; k < 3; k++)
				{
					for (long int l = 0; l < 3; l++)
					{
						out[i][j] = out[i][j] + pt_glb->Tg2c[i][k] * pt_glb->Tg2c[j][l] * in[k][l];
					}
				}
			}
		}
	}
	else if (pt_glb->if_OrientRotate == false) {
		out[0][0] = in[0][0]; out[0][1] = in[0][1]; out[0][2] = in[0][2];
		out[1][0] = in[1][0]; out[1][1] = in[1][1]; out[1][2] = in[1][2];
		out[2][0] = in[2][0]; out[2][1] = in[2][1]; out[2][2] = in[2][2];
	}
}

#pragma acc routine seq// nohost
void mathlib::transform_matrix_crt2glb(const double(&in)[3][3], double(&out)[3][3]) {
	if (pt_glb->if_OrientRotate == true) {
		for (long int i = 0; i < 3; i++)
		{
			for (long int j = 0; j < 3; j++)
			{
				out[i][j] = 0.;
				for (long int k = 0; k < 3; k++)
				{
					for (long int l = 0; l < 3; l++)
					{
						out[i][j] = out[i][j] + pt_glb->Tc2g[i][k] * pt_glb->Tc2g[j][l] * in[k][l];
					}
				}
			}
		}
	}
	else if (pt_glb->if_OrientRotate == false) {
		out[0][0] = in[0][0]; out[0][1] = in[0][1]; out[0][2] = in[0][2];
		out[1][0] = in[1][0]; out[1][1] = in[1][1]; out[1][2] = in[1][2];
		out[2][0] = in[2][0]; out[2][1] = in[2][1]; out[2][2] = in[2][2];
	}
}

#pragma acc routine seq// nohost
void mathlib::strain2stress_glb(const double(&stiff)[6][6], const double(&strain)[6], double(&stress)[6]) {
	for (long int i = 0; i < 6; i++) {
		stress[i] = 0.;
		for (long int j = 0; j < 6; j++) {
			stress[i] = stress[i] + stiff[i][j] * strain[j];
		}
	}
}

void mathlib::copy_to_device() {
#pragma acc enter data copyin(this)

#pragma acc serial default(present)
	{
		this->pt_glb = &(global_parameters::glb);
	}
}

mathlib mathlib::mlb;

