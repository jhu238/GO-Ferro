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

class magnetic_system {
private:
	double mu0 = 1.256637061e-6;
	double unit_e = 1.602176634e-19;
	double PI = 3.14159265358979323846;
public:
	class geometry_parameters* pt_geo = &(geometry_parameters::geo);
	class global_parameters* pt_glb = &(global_parameters::glb);
	class mathlib* pt_math = & (mathlib::mlb);
	class elastic_system* pt_elas;
	class EMdynamic_system* pt_EM;
	long int nx, ny, nz, n;
	double dx, dy, dz;

	matrix3d<bool> FM_surfYZ, FM_surfXZ, FM_surfXY;
	matrix3d<bool> AFM_surfYZ, AFM_surfXZ, AFM_surfXY;

	//----------Averaged parameter for inter-region exchange - Approach 1---------//
	matrix3d<double> Aij_xf, Aij_yf, Aij_zf, Aij_xb, Aij_yb, Aij_zb;
	matrix3d<double> Aij_AFM1_xf, Aij_AFM1_yf, Aij_AFM1_zf, Aij_AFM1_xb, Aij_AFM1_yb, Aij_AFM1_zb;
	matrix3d<double> Aij_AFM2_xf, Aij_AFM2_yf, Aij_AFM2_zf, Aij_AFM2_xb, Aij_AFM2_yb, Aij_AFM2_zb;
	//----------Averaged parameter for inter-region exchange - Approach 2---------//
	//matrix3d<double> AM_exchange;

	//----------Intermediate variables for Magnetostatics----------//
	matrix3d<double> Axyz, Ayzx, Azxy, Bxyz, Byzx, Bzxy, Byxz, Bzyx, Bxzy;
	matrix3d<std::complex<double>> Axyzk, Bxyzk, Bxzyk, Ayzxk, Byzxk, Byxzk, Azxyk, Bzxyk, Bzyxk;
	matrix3d<double> Mx_glb, My_glb, Mz_glb;
	matrix3d<std::complex<double>> Mx_k3D, My_k3D, Mz_k3D;
	matrix3d<std::complex<double>> Hx_stat_k, Hy_stat_k, Hz_stat_k;
	//-------------------------------------------//
	cufftHandle plan_magneto_D2Z;
	cufftHandle plan_magneto_Z2D;
	//friend class inoutput;
private:
	void initialize_intermediate_Hms();

	double A(double& x, double& y, double& z, double& dx, double& dy, double& dz);
	double F1(double& x, double& y, double& z, double& dx, double& dy, double& dz);
	double F2(double& x, double& y, double& z);
	double B(double& x, double& y, double& z, double& dx, double& dy, double& dz);
	double G1(double& x, double& y, double& z, double& dx, double& dy, double& dz);
	double G2(double& x, double& y, double& z);

#pragma acc routine seq// nohost
	void get_m_spatial_derivative_local(long int&, long int&, long int&, \
		double&, double&, \
		double&, double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&);

#pragma acc routine seq// nohost
	void get_m_AFM1_spatial_derivative_local(long int&, long int&, long int&, \
		double&, double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&);

#pragma acc routine seq// nohost
	void get_m_AFM2_spatial_derivative_local(long int&, long int&, long int&, \
		double&, double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&, \
		double&, double&, double&);
	//Laplacian is Output in the form of 
	//d2mx/dx2, d2mx/dy2, d2mx/dz2
	//d2my/dx2, d2my/dy2, d2my/dz2
	//d2mz/dx2, d2mz/dy2, d2mz/dz2

	void initialize_magnetization();
	void initialize_magnetization_AFM();

#pragma acc routine seq
	void update_external_Hfield();
#pragma acc routine seq
	void update_external_Hfield_half();
#pragma acc routine seq
	void update_external_Hfield_full();

	void get_J_ISHE_RK1();
	void get_J_ISHE_RK2();
	void get_J_ISHE_RK3();
	void get_J_ISHE_RK4();

	void get_J_ISHE();
public:
	//---------Fields_FM----------//
	matrix3d<double> mx_glb, my_glb, mz_glb;
	matrix3d<double> dmx_glb, dmy_glb, dmz_glb;

	double mx_ave = 0.; 
	double my_ave = 0.;
	double mz_ave = 0.;
	//---------Fields_AFM----------//
	matrix3d<double> mx_AFM1_glb, my_AFM1_glb, mz_AFM1_glb;
	matrix3d<double> mx_AFM2_glb, my_AFM2_glb, mz_AFM2_glb;
	matrix3d<double> dmx_AFM1_glb, dmy_AFM1_glb, dmz_AFM1_glb;
	matrix3d<double> dmx_AFM2_glb, dmy_AFM2_glb, dmz_AFM2_glb;

	double mx_AFM1_ave = 0.;
	double my_AFM1_ave = 0.;
	double mz_AFM1_ave = 0.;
	double mx_AFM2_ave = 0.;
	double my_AFM2_ave = 0.;
	double mz_AFM2_ave = 0.;

	//--------RK4 FM-------------//
	matrix3d<double> mx_glb_store, my_glb_store, mz_glb_store;
	matrix3d<double> dmx_glb_rk1, dmy_glb_rk1, dmz_glb_rk1;
	matrix3d<double> dmx_glb_rk2, dmy_glb_rk2, dmz_glb_rk2;
	matrix3d<double> dmx_glb_rk3, dmy_glb_rk3, dmz_glb_rk3;
	matrix3d<double> dmx_glb_rk4, dmy_glb_rk4, dmz_glb_rk4;

	//--------RK4 AFM-------------//
	matrix3d<double> mx_AFM1_glb_store, my_AFM1_glb_store, mz_AFM1_glb_store;
	matrix3d<double> mx_AFM2_glb_store, my_AFM2_glb_store, mz_AFM2_glb_store;

	matrix3d<double> dmx_AFM1_glb_rk1, dmy_AFM1_glb_rk1, dmz_AFM1_glb_rk1;
	matrix3d<double> dmx_AFM1_glb_rk2, dmy_AFM1_glb_rk2, dmz_AFM1_glb_rk2;
	matrix3d<double> dmx_AFM1_glb_rk3, dmy_AFM1_glb_rk3, dmz_AFM1_glb_rk3;
	matrix3d<double> dmx_AFM1_glb_rk4, dmy_AFM1_glb_rk4, dmz_AFM1_glb_rk4;

	matrix3d<double> dmx_AFM2_glb_rk1, dmy_AFM2_glb_rk1, dmz_AFM2_glb_rk1;
	matrix3d<double> dmx_AFM2_glb_rk2, dmy_AFM2_glb_rk2, dmz_AFM2_glb_rk2;
	matrix3d<double> dmx_AFM2_glb_rk3, dmy_AFM2_glb_rk3, dmz_AFM2_glb_rk3;
	matrix3d<double> dmx_AFM2_glb_rk4, dmy_AFM2_glb_rk4, dmz_AFM2_glb_rk4;

	//----------Magnetostatic field----------//
	matrix3d<double> Hx_stat, Hy_stat, Hz_stat; // in global coordinate
	//matrix3d<double> Hx_anis, Hy_anis, Hz_anis; // in global coordinate
	//matrix3d<double> Hx_exch, Hy_exch, Hz_exch; // in global coordinate
	//matrix3d<double> Hx_elas, Hy_elas, Hz_elas; // in global coordinate
	// 
	//---------Inverse Spin Hall effect induced charge current---------//
	matrix3d<double> Jx_ISHE, Jy_ISHE;

public:
	//void initialize(elastic_system* pt_elas_in, EMdynamic_system* pt_EM_in);
	void initialize_host();
	void initialize_host_supple();
	void copy_to_device();
	//void allocate_device();

	void initialize_device();

	void copym_from_device();
	void copym_AFM_from_device();
	void copyH_from_device();
	void copyJishe_from_device();

	void get_H_static(/*matrix3d<double>&, matrix3d<double>&, matrix3d<double>&*/);
	//void get_H_static_4RK();
#pragma acc routine gang
	void get_H_static_1D();

	void get_dm();
	void get_dm_AFM();

	void get_averagem();
	void get_averagem_AFM();

	//-----RK4---//
	void get_dm_RK1();
	void get_dm_RK2();
	void get_dm_RK3();
	void get_dm_RK4();

	void update_m_RK1();
	void update_m_RK2();
	void update_m_RK3();
	void update_m();

	void get_dm_AFM_RK1();
	void get_dm_AFM_RK2();
	void get_dm_AFM_RK3();
	void get_dm_AFM_RK4();

	void update_m_AFM_RK1();
	void update_m_AFM_RK2();
	void update_m_AFM_RK3();
	void update_m_AFM();

	//-------Prescribed m for test -----//
#pragma acc routine seq
	void update_prescribed_m_half();
#pragma acc routine seq
	void update_prescribed_m_full();
#pragma acc routine seq
	void update_prescribed_m();
};
