#pragma once
#include "matrix.h"
#include "geometry.h"
#include "global.h"
#include "mathlib.h"
#include <complex>
#include "elastic_system.h"
#include "EMdynamic_system.h"
#include "cufft.h"
#include "openacc.h"

class ferroelectric_system {
public:
	double e0 = 8.854187817e-12;
	std::complex<double> im_unit = { 0.0, 1.0 };
	double PI = 3.14159265358979323846;
	class geometry_parameters *pt_geo = &(geometry_parameters::geo);
	class global_parameters *pt_glb = &(global_parameters::glb);
	class mathlib *pt_math = &(mathlib::mlb);
	class elastic_system* pt_elas;
	class EMdynamic_system* pt_EM;
	long int nx, ny, nz, n, nz21;
	double scalenn;
	double dx, dy, dz;
	//double Estat0_1D;
	matrix3d<std::complex<double>> iGq_inverse;
	double er_homo;
	cufftHandle plan_electro_D2Z;
	cufftHandle plan_electro_Z2D;

	matrix3d<bool> FE_surfYZ, FE_surfXZ, FE_surfXY;

	//----------------Variables for flexo-electrics----------------//
	matrix3d<double> dpxdx, dpxdy, dpxdz;
	matrix3d<double> dpydx, dpydy, dpydz;
	matrix3d<double> dpzdx, dpzdy, dpzdz;

	//----------------Variables for electro-static solver----------------//
	matrix3d<std::complex<double>> rhofk;
	matrix3d<std::complex<double>> Dxk_inhomo, Dyk_inhomo, Dzk_inhomo;
	matrix3d<std::complex<double>> pxk, pyk, pzk;
	matrix3d<std::complex<double>> ps_terms;
	matrix3d<std::complex<double>> Edxk_n, Edyk_n, Edzk_n;

	matrix3d<double> Dx_inhomo, Dy_inhomo, Dz_inhomo;
	double e_field[4];
	double e_sum;
	double tolerance;
private:
	#pragma acc routine seq// nohost
	void get_laplacian_p_local(long int&, long int&, long int&, \
		double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&);

	void initialize_polarization();

#pragma acc routine seq
	void update_external_Efield();
#pragma acc routine seq
	void update_external_Efield_half();
#pragma acc routine seq
	void update_external_Efield_full();

	//#pragma acc routine nohost
	//double get_Edepolar_sum(/*matrix3d<double>&, matrix3d<double>&, matrix3d<double>&, \
	//	matrix3d<double>&, matrix3d<double>&, matrix3d<double>&, *//*double&*/);
public:
	matrix3d<double> px_glb, py_glb, pz_glb;
	matrix3d<double> pz_t0_glb;
//	matrix3d<double> dpx_glb, dpy_glb, dpz_glb;
	matrix3d<double> qx_glb, qy_glb, qz_glb;
//	matrix3d<double> dqx_glb, dqy_glb, dqz_glb;

	//RK4
	matrix3d<double> px_glb_store, py_glb_store, pz_glb_store;
	matrix3d<double> dpx_glb_rk1, dpy_glb_rk1, dpz_glb_rk1;
	matrix3d<double> dpx_glb_rk2, dpy_glb_rk2, dpz_glb_rk2;
	matrix3d<double> dpx_glb_rk3, dpy_glb_rk3, dpz_glb_rk3;
	matrix3d<double> dpx_glb_rk4, dpy_glb_rk4, dpz_glb_rk4;

	matrix3d<double> qx_glb_store, qy_glb_store, qz_glb_store;
	matrix3d<double> dqx_glb_rk1, dqy_glb_rk1, dqz_glb_rk1;
	matrix3d<double> dqx_glb_rk2, dqy_glb_rk2, dqz_glb_rk2;
	matrix3d<double> dqx_glb_rk3, dqy_glb_rk3, dqz_glb_rk3;
	matrix3d<double> dqx_glb_rk4, dqy_glb_rk4, dqz_glb_rk4;

	matrix3d<double> Ex_stat, Ey_stat, Ez_stat; // in global coordinate

	matrix3d<double> rhof;

	double px_ave = 0.;
	double py_ave = 0.;
	double pz_ave = 0.;
public:
	void initialize_host();
	void copy_to_device();
	void initialize_device();
	void copyp_from_device();
	void copyq_from_device();
	void copyE_from_device();

	void get_E_static();
	//void get_E_static_4RK();
#pragma acc routine gang
	void get_E_static_1D();

	void get_dq_RK1();
	void get_dq_RK2();
	void get_dq_RK3();
	void get_dq_RK4();

	void get_dp_RK1();
	void get_dp_RK2();
	void get_dp_RK3();
	void get_dp_RK4();

	void update_q_RK1();
	void update_q_RK2();
	void update_q_RK3();
	void update_q();

	void update_p_RK1();
	void update_p_RK2();
	void update_p_RK3();
	void update_p();

#pragma acc routine gang
	void get_p_gradient();

	void get_averagep();
};
