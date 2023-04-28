#pragma once
#include "matrix.h"
#include "geometry.h"
#include "global.h"
#include "mathlib.h"
#include "magnetic_system.h"
#include "ferroelectric_system.h"
#include "openacc.h"

class EMdynamic_system {
private:
	double mu0 = 1.256637061e-6;
	double e0 = 8.854187817e-12;
	double c = 2.99792458e8;
	double PI = 3.14159265358979323846;
public:
	class geometry_parameters *pt_geo = &(geometry_parameters::geo);
	class global_parameters *pt_glb = &(global_parameters::glb);
	class mathlib* pt_math = &(mathlib::mlb);
	class magnetic_system *pt_mag;
	class ferroelectric_system *pt_fe;
	long int nx, ny, nz, n;
	double dx, dy, dz;

	unsigned long planeEM_source_z;
	double sqrt_er_top, sqrt_er_bot;

	double coeff_top_half, coeff_bot_half; //For Mur ABC
	double coeff_x_half, coeff_y_half; //For Mur in-plane ABC 

	double coeff_top, coeff_bot; //For Mur ABC
	double coeff_x, coeff_y; //For Mur in-plane ABC 

	double w; //For Liao ABC

	//unsigned int Jf_z1, Jf_z2;
	//double Jf_thickness;
	unsigned long Jf_nx, Jf_ny, Jf_nz, Jf_n;

	matrix3d<double> DEx_em_t1, DEy_em_t1, DEz_em_t1;
	matrix3d<double> DEx_em_t2, DEy_em_t2, DEz_em_t2;
	matrix3d<double> DEx_em_t3, DEy_em_t3, DEz_em_t3;
	matrix3d<double> DEx_em_t4, DEy_em_t4, DEz_em_t4;

	matrix3d<double> DHx_em, DHy_em, DHz_em; // in global coordinate
	matrix3d<double> DEx_em, DEy_em, DEz_em; // in global coordinate

	matrix3d<bool> DEx_ifPEC, DEy_ifPEC, DEz_ifPEC;

	//matrix3d<double> dDHx_em, dDHy_em, dDHz_em; // in global coordinate
	//matrix3d<double> dDEx_em, dDEy_em, dDEz_em; // in global coordinate

	//RK4
	matrix3d<double> DHx_em_store, DHy_em_store, DHz_em_store; // in global coordinate
	matrix3d<double> DEx_em_store, DEy_em_store, DEz_em_store; // in global coordinate

	matrix3d<double> dDHx_em_rk1, dDHy_em_rk1, dDHz_em_rk1; // in global coordinate
	matrix3d<double> dDEx_em_rk1, dDEy_em_rk1, dDEz_em_rk1; // in global coordinate

	matrix3d<double> dDHx_em_rk2, dDHy_em_rk2, dDHz_em_rk2; // in global coordinate
	matrix3d<double> dDEx_em_rk2, dDEy_em_rk2, dDEz_em_rk2; // in global coordinate

	matrix3d<double> dDHx_em_rk3, dDHy_em_rk3, dDHz_em_rk3; // in global coordinate
	matrix3d<double> dDEx_em_rk3, dDEy_em_rk3, dDEz_em_rk3; // in global coordinate

	matrix3d<double> dDHx_em_rk4, dDHy_em_rk4, dDHz_em_rk4; // in global coordinate
	matrix3d<double> dDEx_em_rk4, dDEy_em_rk4, dDEz_em_rk4; // in global coordinate


	//friend class inoutput;
private:
	void transfer_pointer();
	void initialize_PEC_edges();

	void update_DE_Boundary_half();
	void update_DE_Boundary_full();
	void update_DE_Boundary();

	void update_DE_PEC();

#pragma acc routine gang// nohost
	void update_planeEM_half();

#pragma acc routine gang// nohost
	void update_planeEM_full();

#pragma acc routine gang// nohost
	void update_planeEM();

	void update_Jf_input();
	void update_Jf_input_half();
	void update_Jf_input_full();

	//void update_Jp();

#pragma acc routine gang// nohost
	void update_DH_cell();
#pragma acc routine gang// nohost
	void update_DE_cell();

public:
	matrix3d<double> DHx_em_cell, DHy_em_cell, DHz_em_cell; // in global coordinate
	matrix3d<double> DEx_em_cell, DEy_em_cell, DEz_em_cell; // in global coordinate

	matrix3d<double> DJfx, DJfy, DJfz;

	matrix3d<double> Jpx_n1, Jpy_n1, Jpz_n1;
	matrix3d<double> Jpx_n2, Jpy_n2, Jpz_n2;
	matrix3d<double> Jpx_n3, Jpy_n3, Jpz_n3;

	matrix3d<double> Jpx_n1_store, Jpy_n1_store, Jpz_n1_store;
	matrix3d<double> Jpx_n2_store, Jpy_n2_store, Jpz_n2_store;
	matrix3d<double> Jpx_n3_store, Jpy_n3_store, Jpz_n3_store;

	matrix3d<double> dJpx_n1_rk1, dJpy_n1_rk1, dJpz_n1_rk1;
	matrix3d<double> dJpx_n2_rk1, dJpy_n2_rk1, dJpz_n2_rk1;
	matrix3d<double> dJpx_n3_rk1, dJpy_n3_rk1, dJpz_n3_rk1;

	matrix3d<double> dJpx_n1_rk2, dJpy_n1_rk2, dJpz_n1_rk2;
	matrix3d<double> dJpx_n2_rk2, dJpy_n2_rk2, dJpz_n2_rk2;
	matrix3d<double> dJpx_n3_rk2, dJpy_n3_rk2, dJpz_n3_rk2;

	matrix3d<double> dJpx_n1_rk3, dJpy_n1_rk3, dJpz_n1_rk3;
	matrix3d<double> dJpx_n2_rk3, dJpy_n2_rk3, dJpz_n2_rk3;
	matrix3d<double> dJpx_n3_rk3, dJpy_n3_rk3, dJpz_n3_rk3;

	matrix3d<double> dJpx_n1_rk4, dJpy_n1_rk4, dJpz_n1_rk4;
	matrix3d<double> dJpx_n2_rk4, dJpy_n2_rk4, dJpz_n2_rk4;
	matrix3d<double> dJpx_n3_rk4, dJpy_n3_rk4, dJpz_n3_rk4;

public:
	void initialize_host();
	void copy_to_device();
	void initialize_device();
	void copy_from_device();
	void copyYee_from_device();
	void copyJp_from_device();

//RK4
	void get_dH_RK1();
	void get_dH_RK2();
	void get_dH_RK3();
	void get_dH_RK4();

	void get_dE_RK1();
	void get_dE_RK2();
	void get_dE_RK3();
	void get_dE_RK4();

#pragma acc routine gang
	void get_dJp_RK1();
#pragma acc routine gang
	void get_dJp_RK2();
#pragma acc routine gang
	void get_dJp_RK3();
#pragma acc routine gang
	void get_dJp_RK4();

	void update_DH_RK1();
	void update_DH_RK2();
	void update_DH_RK3();

	void update_DE_RK1();
	void update_DE_RK2();
	void update_DE_RK3();

#pragma acc routine gang
	void update_Jp_RK1();
#pragma acc routine gang
	void update_Jp_RK2();
#pragma acc routine gang
	void update_Jp_RK3();

//
	void update_DH();
	void update_DE();

#pragma acc routine gang
	void update_Jp();
};
