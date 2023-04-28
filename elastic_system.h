#pragma once
#include "matrix.h"
#include "geometry.h"
#include "global.h"
#include "mathlib.h"
#include <complex>
#include <fstream>
#include "magnetic_system.h"
#include "ferroelectric_system.h"
#include "cufft.h"
#include "openacc.h"

class elastic_system {
private:
	double PI = 3.14159265358979323846;
public:
	std::complex<double> im_unit = { 0.0, 1.0 };
	class geometry_parameters* pt_geo =&(geometry_parameters::geo);
	class global_parameters* pt_glb = &(global_parameters::glb);
	class mathlib* pt_math = &(mathlib::mlb);
	class magnetic_system *pt_mag;
	class ferroelectric_system *pt_fe;
	long int nx, ny, nz, nz21, n;
	double scalenn;
	double dx, dy, dz;
	
	double c11_1st, c44_1st, density_1st;
	double c11_last, c44_last, density_last;
	double vs_long_1st, vs_tran_1st;
	double vs_long_last, vs_tran_last;
	unsigned int source_z;

	matrix3d<bool> stress_free_surfYZ, stress_free_surfXZ, stress_free_surfXY;
	double c11_homo, c12_homo, c44_homo;
	matrix3d<std::complex<double>> Gik_11, Gik_12, Gik_13, Gik_21, Gik_22, Gik_23, Gik_31, Gik_32, Gik_33;

	matrix3d<double> duxdx_glb_surfYZ, duydx_glb_surfYZ, duzdx_glb_surfYZ;
	matrix3d<double> duxdy_glb_surfXZ, duydy_glb_surfXZ, duzdy_glb_surfXZ;
	matrix3d<double> duxdz_glb_surfXY, duydz_glb_surfXY, duzdz_glb_surfXY;

	matrix3d<double> Dux_glb_t1, Duy_glb_t1, Duz_glb_t1;
	matrix3d<double> Dux_glb_t2, Duy_glb_t2, Duz_glb_t2;
	matrix3d<double> Dux_glb_t3, Duy_glb_t3, Duz_glb_t3;
	matrix3d<double> Dux_glb_t4, Duy_glb_t4, Duz_glb_t4;

	matrix3d<double> dexxt0dx_glb, dexxt0dy_glb, dexxt0dz_glb;
	matrix3d<double> deyyt0dx_glb, deyyt0dy_glb, deyyt0dz_glb;
	matrix3d<double> dezzt0dx_glb, dezzt0dy_glb, dezzt0dz_glb;
	matrix3d<double> deyzt0dx_glb, deyzt0dy_glb, deyzt0dz_glb;
	matrix3d<double> dexzt0dx_glb, dexzt0dy_glb, dexzt0dz_glb;
	matrix3d<double> dexyt0dx_glb, dexyt0dy_glb, dexyt0dz_glb;

	matrix3d<double> dDexxdx_glb, dDexxdy_glb, dDexxdz_glb;
	matrix3d<double> dDeyydx_glb, dDeyydy_glb, dDeyydz_glb;
	matrix3d<double> dDezzdx_glb, dDezzdy_glb, dDezzdz_glb;
	matrix3d<double> dDeyzdx_glb, dDeyzdy_glb, dDeyzdz_glb;
	matrix3d<double> dDexzdx_glb, dDexzdy_glb, dDexzdz_glb;
	matrix3d<double> dDexydx_glb, dDexydy_glb, dDexydz_glb;

//----------------Variables for elasto-static solver----------------//
private:
	double e_sum, e_elas[4];
	double tolerance;
	cufftHandle plan_elasto_D2Z;
	cufftHandle plan_elasto_Z2D;

	matrix3d<double> exx_ext_crt, eyy_ext_crt, ezz_ext_crt;
	matrix3d<double> exy_ext_crt, exz_ext_crt, eyz_ext_crt;

	matrix3d<double> stress0homo_xx, stress0homo_yy, stress0homo_zz, \
		stress0homo_xy, stress0homo_xz, stress0homo_yz;

	matrix3d<double> stress0_xx, stress0_yy, stress0_zz, \
		stress0_xy, stress0_xz, stress0_yz;

	matrix3d<double> stress_xx, stress_yy, stress_zz, \
		stress_xy, stress_xz, stress_yz;

	matrix3d<double> duxdx, duxdy, duxdz, \
		duydx, duydy, duydz, \
		duzdx, duzdy, duzdz;

	matrix3d<std::complex<double>> stress0homok_xx, stress0homok_yy, stress0homok_zz, \
		stress0homok_xy, stress0homok_xz, stress0homok_yz;

	matrix3d<std::complex<double>> stressk_xx, stressk_yy, stressk_zz, \
		stressk_xy, stressk_xz, stressk_yz;

	matrix3d<std::complex<double>> duxdxk, duxdyk, duxdzk, \
		duydxk, duydyk, duydzk, \
		duzdxk, duzdyk, duzdzk;

	double* exx0_ave_film, * eyy0_ave_film, * ezz0_ave_film;
	double* eyz0_ave_film, * exz0_ave_film, * exy0_ave_film;

	//matrix3d<std::complex<double>> uxk(nx, ny, nz21), uyk(nx, ny, nz21), uzk(nx, ny, nz21);
	//matrix3d<double> ux(nx, ny, nz), uy(nx, ny, nz), uz(nx, ny, nz);
		//matrix3d<double> exx(nx, ny, nz), eyy(nx, ny, nz), ezz(nx, ny, nz), \
	//	exy(nx, ny, nz), exz(nx, ny, nz), eyz(nx, ny, nz);

private:
//#pragma acc routine nohost
//	void get_Eelas_sum(/*matrix3d<double>&, matrix3d<double>&, matrix3d<double>&, \
//		matrix3d<double>&, matrix3d<double>&, matrix3d<double>&, \
//		matrix3d<double>&, matrix3d<double>&, matrix3d<double>&, \
//		matrix3d<double>&, matrix3d<double>&, matrix3d<double>&, */double&);

#pragma acc routine gang// nohost
	void get_eigenstrain_dynamic();

	void transfer_pointer();

#pragma acc routine gang// nohost
	void update_Du_Boundary_half();
#pragma acc routine gang// nohost
	void update_Du_Boundary_full();
#pragma acc routine gang// nohost
	void update_Du_Boundary();

#pragma acc routine gang// nohost
	void update_Du_Gaussian_half();
#pragma acc routine gang// nohost
	void update_Du_Gaussian_full();
#pragma acc routine gang// nohost
	void update_Du_Gaussian();

#pragma acc routine gang// nohost
	void get_dudxyz_glb();

#pragma acc routine gang// nohost
	void get_Dstrain();

#pragma acc routine gang// nohost
	void get_Dstrain_gradient();

	void read_external_strain();
public:
	matrix3d<double> Dux_glb, Duy_glb, Duz_glb;
	//matrix3d<double> dDux_glb, dDuy_glb, dDuz_glb;

	matrix3d<double> vx_glb, vy_glb, vz_glb;
	//matrix3d<double> dvx_glb, dvy_glb, dvz_glb;

	matrix3d<double> Dexx_glb, Deyy_glb, Dezz_glb, Deyz_glb, Dexz_glb, Dexy_glb;
	matrix3d<double> exxt0_glb, eyyt0_glb, ezzt0_glb, eyzt0_glb, exzt0_glb, exyt0_glb;

	matrix3d<double> Dexx_crt, Deyy_crt, Dezz_crt, Deyz_crt, Dexz_crt, Dexy_crt;
	matrix3d<double> exxt0_crt, eyyt0_crt, ezzt0_crt, eyzt0_crt, exzt0_crt, exyt0_crt;

	matrix3d<double> Dexx0_glb, Deyy0_glb, Dezz0_glb, Deyz0_glb, Dexz0_glb, Dexy0_glb;
	//matrix3d<double> exx0t0_glb, eyy0t0_glb, ezz0t0_glb, eyz0t0_glb, exz0t0_glb, exy0t0_glb;

	matrix3d<double> Dexx0_crt, Deyy0_crt, Dezz0_crt, Deyz0_crt, Dexz0_crt, Dexy0_crt;
	matrix3d<double> exx0t0_crt, eyy0t0_crt, ezz0t0_crt, eyz0t0_crt, exz0t0_crt, exy0t0_crt;

	//RK4
	matrix3d<double> force_x, force_y, force_z;
	matrix3d<double> force_x_store, force_y_store, force_z_store;

	matrix3d<double> Dux_glb_store, Duy_glb_store, Duz_glb_store;
	matrix3d<double> dDux_glb_rk1, dDuy_glb_rk1, dDuz_glb_rk1;
	matrix3d<double> dDux_glb_rk2, dDuy_glb_rk2, dDuz_glb_rk2;
	matrix3d<double> dDux_glb_rk3, dDuy_glb_rk3, dDuz_glb_rk3;
	matrix3d<double> dDux_glb_rk4, dDuy_glb_rk4, dDuz_glb_rk4;

	matrix3d<double> vx_glb_store, vy_glb_store, vz_glb_store;
	matrix3d<double> dvx_glb_rk1, dvy_glb_rk1, dvz_glb_rk1;
	matrix3d<double> dvx_glb_rk2, dvy_glb_rk2, dvz_glb_rk2;
	matrix3d<double> dvx_glb_rk3, dvy_glb_rk3, dvz_glb_rk3;
	matrix3d<double> dvx_glb_rk4, dvy_glb_rk4, dvz_glb_rk4;

public:
	void initialize_host();
	void copy_to_device();

	void initialize_device();

	void copy_Dstrain_from_device();
	void copy_straint0_from_device();
	void copy_uandv_from_device();
	void copy_eigenstraint0_from_device();
	void copy_elastoforce_from_device();

#pragma acc routine gang// nohost
	void get_eigenstrain_static();

#pragma acc routine gang
	void get_eigenstrain_average_film();

	void get_strain_static();

#pragma acc routine gang
	void get_strain_static_1D();

#pragma acc routine gang// nohost
	void get_strain_static_glb();

#pragma acc routine gang// nohost
	void get_strain_static_glb_gradient();

	//void get_dv();
	//void get_du();

	//-----RK4---//
	void get_dv_RK1();
	void get_dv_RK2();
	void get_dv_RK3();
	void get_dv_RK4();

	void get_du_RK1();
	void get_du_RK2();
	void get_du_RK3();
	void get_du_RK4();

	void update_v_RK1();
	void update_u_RK1();
	void update_v_RK2();
	void update_u_RK2();
	void update_v_RK3();
	void update_u_RK3();
	//-----------//

	void update_v();
	void update_u();
};
