#include "elastic_system.h"

#pragma acc routine gang// nohost
void elastic_system::get_eigenstrain_dynamic() {
	unsigned int mat_type;
	material* mat;
	double mx, my, mz;
	double mx_AFM1, my_AFM1, mz_AFM1;
	double mx_AFM2, my_AFM2, mz_AFM2;
	double px, py, pz;
	double Q11, Q12, Q44;
	double F11, F12, F44;
	double lamda100, lamda111;
	double lamda100_AFM1, lamda111_AFM1;
	double lamda100_AFM2, lamda111_AFM2;
	double strain_in[3][3], strain_out[3][3];

#pragma acc loop gang vector private(mat_type,mx, my, mz,mx_AFM1, my_AFM1, mz_AFM1,mx_AFM2, my_AFM2, mz_AFM2,px, py, pz,\
Q11, Q12, Q44,F11, F12, F44,lamda100, lamda111,lamda100_AFM1, lamda111_AFM1,lamda100_AFM2, lamda111_AFM2,strain_in,strain_out,mat)
	for (long int id = 0; id < n; id++) {
		mat_type = pt_glb->material_cell(id);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[mat_type - 1]);
			Q11 = mat->Q11; Q12 = mat->Q12; Q44 = mat->Q44;
			F11 = mat->F11; F12 = mat->F12; F44 = mat->F44;
			lamda100 = mat->lamda100; lamda111 = mat->lamda111;
			lamda100_AFM1 = mat->lamda100_AFM1; lamda111_AFM1 = mat->lamda111_AFM1;
			lamda100_AFM2 = mat->lamda100_AFM2; lamda111_AFM2 = mat->lamda111_AFM2;

			pt_math->transform_vector_glb2crt(pt_fe->px_glb_store(id), pt_fe->py_glb_store(id), pt_fe->pz_glb_store(id), \
				px, py, pz);
			pt_math->transform_vector_glb2crt(pt_mag->mx_glb_store(id), pt_mag->my_glb_store(id), pt_mag->mz_glb_store(id), \
				mx, my, mz);
			pt_math->transform_vector_glb2crt(pt_mag->mx_AFM1_glb_store(id), pt_mag->my_AFM1_glb_store(id), pt_mag->mz_AFM1_glb_store(id), \
				mx_AFM1, my_AFM1, mz_AFM1);
			pt_math->transform_vector_glb2crt(pt_mag->mx_AFM2_glb_store(id), pt_mag->my_AFM2_glb_store(id), pt_mag->mz_AFM2_glb_store(id), \
				mx_AFM2, my_AFM2, mz_AFM2);

			strain_in[0][0] = 1.5 * lamda100 * (pow(mx, 2.) - 1. / 3.) + Q11 * pow(px, 2.) + Q12 * (pow(py, 2.) + pow(pz, 2.))\
				+ 0.5 * 1.5 * lamda100_AFM1 * (pow(mx_AFM1, 2.) - 1. / 3.) + 0.5 * 1.5 * lamda100_AFM2 * (pow(mx_AFM2, 2.) - 1. / 3.) \
				- exx0t0_crt(id);

			strain_in[1][1] = 1.5 * lamda100 * (pow(my, 2.) - 1. / 3.) + Q11 * pow(py, 2.) + Q12 * (pow(px, 2.) + pow(pz, 2.))\
				+ 0.5 * 1.5 * lamda100_AFM1 * (pow(my_AFM1, 2.) - 1. / 3.) + 0.5 * 1.5 * lamda100_AFM2 * (pow(my_AFM2, 2.) - 1. / 3.) \
				- eyy0t0_crt(id);

			strain_in[2][2] = 1.5 * lamda100 * (pow(mz, 2.) - 1. / 3.) + Q11 * pow(pz, 2.) + Q12 * (pow(py, 2.) + pow(px, 2.))\
				+ 0.5 * 1.5 * lamda100_AFM1 * (pow(mz_AFM1, 2.) - 1. / 3.) + 0.5 * 1.5 * lamda100_AFM2 * (pow(mz_AFM2, 2.) - 1. / 3.) \
				- ezz0t0_crt(id);

			strain_in[1][2] = 1.5 * lamda111 * my * mz + Q44 * py * pz \
				+ 0.5 * 1.5 * lamda111_AFM1 * my_AFM1 * mz_AFM1 + 0.5 * 1.5 * lamda111_AFM2 * my_AFM2 * mz_AFM2 \
				- eyz0t0_crt(id);

			strain_in[0][2] = 1.5 * lamda111 * mx * mz + Q44 * px * pz \
				+ 0.5 * 1.5 * lamda111_AFM1 * mx_AFM1 * mz_AFM1 + 0.5 * 1.5 * lamda111_AFM2 * mx_AFM2 * mz_AFM2 \
				- exz0t0_crt(id);

			strain_in[0][1] = 1.5 * lamda111 * my * mx + Q44 * py * px \
				+ 0.5 * 1.5 * lamda111_AFM1 * my_AFM1 * mx_AFM1 + 0.5 * 1.5 * lamda111_AFM2 * my_AFM2 * mx_AFM2 \
				- exy0t0_crt(id);

			if (pt_glb->if_FE_all == true && pt_glb->if_flexo == true) {
				strain_in[0][0] = strain_in[0][0] - F11 * pt_fe->dpxdx(id) - F12 * (pt_fe->dpzdz(id) + pt_fe->dpydy(id));
				strain_in[1][1] = strain_in[1][1] - F11 * pt_fe->dpydy(id) - F12 * (pt_fe->dpzdz(id) + pt_fe->dpxdx(id));
				strain_in[2][2] = strain_in[2][2] - F11 * pt_fe->dpzdz(id) - F12 * (pt_fe->dpxdx(id) + pt_fe->dpydy(id));
				strain_in[1][2] = strain_in[1][2] - F44 / 2. * (pt_fe->dpydz(id) + pt_fe->dpzdy(id));
				strain_in[0][2] = strain_in[0][2] - F44 / 2. * (pt_fe->dpxdz(id) + pt_fe->dpzdx(id));
				strain_in[0][1] = strain_in[0][1] - F44 / 2. * (pt_fe->dpydx(id) + pt_fe->dpxdy(id));
			}

			Dexx0_crt(id) = strain_in[0][0];

			Deyy0_crt(id) = strain_in[1][1];

			Dezz0_crt(id) = strain_in[2][2];

			Deyz0_crt(id) = strain_in[1][2];
			strain_in[2][1] = strain_in[1][2];

			Dexz0_crt(id) = strain_in[0][2];
			strain_in[2][0] = strain_in[0][2];

			Dexy0_crt(id) = strain_in[0][1];
			strain_in[1][0] = strain_in[0][1];

			pt_math->transform_matrix_crt2glb(strain_in, strain_out);
			Dexx0_glb(id) = strain_out[0][0];
			Deyy0_glb(id) = strain_out[1][1];
			Dezz0_glb(id) = strain_out[2][2];
			Deyz0_glb(id) = strain_out[1][2];
			Dexz0_glb(id) = strain_out[0][2];
			Dexy0_glb(id) = strain_out[0][1];
		}
		else {
			Dexx0_crt(id) = 0.;  Dexx0_glb(id) = 0.;
			Deyy0_crt(id) = 0.;  Deyy0_glb(id) = 0.;
			Dezz0_crt(id) = 0.;  Dezz0_glb(id) = 0.;
			Deyz0_crt(id) = 0.;  Deyz0_glb(id) = 0.;
			Dexz0_crt(id) = 0.;  Dexz0_glb(id) = 0.;
			Dexy0_crt(id) = 0.;  Dexy0_glb(id) = 0.;
		}
	}
}

#pragma acc routine gang// nohost
void elastic_system::get_dudxyz_glb() {
	long int i, j, k;

	//----------X direction - YZ Plane---------------//
#pragma acc loop gang vector private(i,j,k)
	for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
		i = id / (ny * nz);
		j = (id - i * (ny * nz)) / nz;
		k = id - i * (ny * nz) - j * nz;

		if (stress_free_surfYZ(i, j, k) == true) {
			duxdx_glb_surfYZ(i, j, k) = 0.;
			duydx_glb_surfYZ(i, j, k) = 0.;
			duzdx_glb_surfYZ(i, j, k) = 0.;
		}
		else {
			if (i == 0 || i == nx) {
				duxdx_glb_surfYZ(i, j, k) = (Dux_glb_store(0, j, k) - Dux_glb_store(nx - 1, j, k)) / dx;
				duydx_glb_surfYZ(i, j, k) = (Duy_glb_store(0, j, k) - Duy_glb_store(nx - 1, j, k)) / dx;
				duzdx_glb_surfYZ(i, j, k) = (Duz_glb_store(0, j, k) - Duz_glb_store(nx - 1, j, k)) / dx;
			}
			else {
				duxdx_glb_surfYZ(i, j, k) = (Dux_glb_store(i, j, k) - Dux_glb_store(i - 1, j, k)) / dx;
				duydx_glb_surfYZ(i, j, k) = (Duy_glb_store(i, j, k) - Duy_glb_store(i - 1, j, k)) / dx;
				duzdx_glb_surfYZ(i, j, k) = (Duz_glb_store(i, j, k) - Duz_glb_store(i - 1, j, k)) / dx;
			}
		}
	}

	//----------Y direction - XZ Plane---------------//
#pragma acc loop gang vector private(i,j,k)
	for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * nz);
		j = (id - i * ((ny + 1) * nz)) / nz;
		k = id - i * ((ny + 1) * nz) - j * nz;

		if (stress_free_surfXZ(i, j, k) == true) {
			duxdy_glb_surfXZ(i, j, k) = 0.;
			duydy_glb_surfXZ(i, j, k) = 0.;
			duzdy_glb_surfXZ(i, j, k) = 0.;
		}
		else {
			if (j == 0 || j == ny) {
				duxdy_glb_surfXZ(i, j, k) = (Dux_glb_store(i, 0, k) - Dux_glb_store(i, ny - 1, k)) / dy;
				duydy_glb_surfXZ(i, j, k) = (Duy_glb_store(i, 0, k) - Duy_glb_store(i, ny - 1, k)) / dy;
				duzdy_glb_surfXZ(i, j, k) = (Duz_glb_store(i, 0, k) - Duz_glb_store(i, ny - 1, k)) / dy;
			}
			else {
				duxdy_glb_surfXZ(i, j, k) = (Dux_glb_store(i, j, k) - Dux_glb_store(i, j - 1, k)) / dy;
				duydy_glb_surfXZ(i, j, k) = (Duy_glb_store(i, j, k) - Duy_glb_store(i, j - 1, k)) / dy;
				duzdy_glb_surfXZ(i, j, k) = (Duz_glb_store(i, j, k) - Duz_glb_store(i, j - 1, k)) / dy;
			}
		}
	}

	//----------Z direction - XY Plane---------------//
#pragma acc loop gang vector private(i,j,k)
	for (long int id = 0; id < nx * ny * (nz + 1); id++) {
		i = id / (ny * (nz + 1));
		j = (id - i * (ny * (nz + 1))) / (nz + 1);
		k = id - i * (ny * (nz + 1)) - j * (nz + 1);

		if (stress_free_surfXY(i, j, k) == true) {
			duxdz_glb_surfXY(i, j, k) = 0.;
			duydz_glb_surfXY(i, j, k) = 0.;
			duzdz_glb_surfXY(i, j, k) = 0.;
		}
		else {
			if (k == 0 || k == nz) {
				duxdz_glb_surfXY(i, j, k) = (Dux_glb_store(i, j, 0) - Dux_glb_store(i, j, nz - 1)) / dz;
				duydz_glb_surfXY(i, j, k) = (Duy_glb_store(i, j, 0) - Duy_glb_store(i, j, nz - 1)) / dz;
				duzdz_glb_surfXY(i, j, k) = (Duz_glb_store(i, j, 0) - Duz_glb_store(i, j, nz - 1)) / dz;
			}
			else {
				duxdz_glb_surfXY(i, j, k) = (Dux_glb_store(i, j, k) - Dux_glb_store(i, j, k - 1)) / dz;
				duydz_glb_surfXY(i, j, k) = (Duy_glb_store(i, j, k) - Duy_glb_store(i, j, k - 1)) / dz;
				duzdz_glb_surfXY(i, j, k) = (Duz_glb_store(i, j, k) - Duz_glb_store(i, j, k - 1)) / dz;
			}
		}
	}
}

#pragma acc routine gang// nohost
void elastic_system::get_Dstrain() {
	long int i, j, k;
	double strain_in[3][3];
	double strain_out[3][3];
	//#pragma acc declare device_resident(strain_in,strain_out)

#pragma acc loop gang vector private(strain_in,strain_out,i,j,k)
	for (long int id = 0; id < n; id++) {
		i = id / (ny * nz);
		j = (id - i * (ny * nz)) / nz;
		k = id - i * (ny * nz) - j * nz;

		strain_in[0][0] = (duxdx_glb_surfYZ(i, j, k) + duxdx_glb_surfYZ(i + 1, j, k)) / 2.;
		strain_in[1][1] = (duydy_glb_surfXZ(i, j, k) + duydy_glb_surfXZ(i, j + 1, k)) / 2.;
		strain_in[2][2] = (duzdz_glb_surfXY(i, j, k) + duzdz_glb_surfXY(i, j, k + 1)) / 2.;

		strain_in[1][2] = (duydz_glb_surfXY(i, j, k) + duydz_glb_surfXY(i, j, k + 1) + \
			duzdy_glb_surfXZ(i, j, k) + duzdy_glb_surfXZ(i, j + 1, k)) / 4.;
		strain_in[2][1] = strain_in[1][2];

		strain_in[0][2] = (duxdz_glb_surfXY(i, j, k) + duxdz_glb_surfXY(i, j, k + 1) + \
			duzdx_glb_surfYZ(i, j, k) + duzdx_glb_surfYZ(i + 1, j, k)) / 4.;
		strain_in[2][0] = strain_in[0][2];

		strain_in[0][1] = (duxdy_glb_surfXZ(i, j, k) + duxdy_glb_surfXZ(i, j + 1, k) + \
			duydx_glb_surfYZ(i, j, k) + duydx_glb_surfYZ(i + 1, j, k)) / 4.;
		strain_in[1][0] = strain_in[0][1];

		pt_math->transform_matrix_glb2crt(strain_in, strain_out);

		Dexx_glb(id) = strain_in[0][0]; Dexx_crt(id) = strain_out[0][0];
		Deyy_glb(id) = strain_in[1][1]; Deyy_crt(id) = strain_out[1][1];
		Dezz_glb(id) = strain_in[2][2]; Dezz_crt(id) = strain_out[2][2];
		Deyz_glb(id) = strain_in[1][2]; Deyz_crt(id) = strain_out[1][2];
		Dexz_glb(id) = strain_in[0][2]; Dexz_crt(id) = strain_out[0][2];
		Dexy_glb(id) = strain_in[0][1]; Dexy_crt(id) = strain_out[0][1];
	}
}

#pragma acc routine gang// nohost
void elastic_system::get_Dstrain_gradient() {
	long i, j, k;
	unsigned int mat_type;
	long fwd, bwd;
	double scale;

#pragma acc loop gang vector private(i,j,k,fwd,bwd,scale,mat_type)
	for (long id = 0; id < n; id++) {
		mat_type = pt_glb->material_cell(id);
		if (mat_type == 0) {
			dDexxdx_glb(id) = 0.; dDexxdy_glb(id) = 0.; dDexxdz_glb(id) = 0.;
			dDeyydx_glb(id) = 0.; dDeyydy_glb(id) = 0.; dDeyydz_glb(id) = 0.;
			dDezzdx_glb(id) = 0.; dDezzdy_glb(id) = 0.; dDezzdz_glb(id) = 0.;
			dDeyzdx_glb(id) = 0.; dDeyzdy_glb(id) = 0.; dDeyzdz_glb(id) = 0.;
			dDexzdx_glb(id) = 0.; dDexzdy_glb(id) = 0.; dDexzdz_glb(id) = 0.;
			dDexydx_glb(id) = 0.; dDexydy_glb(id) = 0.; dDexydz_glb(id) = 0.;
		}
		else {
			i = id / (ny * nz);
			j = (id - i * (ny * nz)) / nz;
			k = id - i * (ny * nz) - j * nz;
			//------------------spatial derivative along X--------------//
			if (stress_free_surfYZ(i, j, k) == stress_free_surfYZ(i + 1, j, k)) {
				if (stress_free_surfYZ(i, j, k) == true) {
					bwd = i; fwd = i; scale = 1.;
				}
				else {
					scale = 2.;
					if (i == 0) {
						bwd = nx - 1;
					}
					else {
						bwd = i - 1;
					}
					if (i == nx - 1) {
						fwd = 0;
					}
					else {
						fwd = i + 1;
					}
				}
			}
			else {
				if (stress_free_surfYZ(i, j, k) == true) {
					bwd = i; scale = 1.;
					if (i == nx - 1) {
						fwd = 0;
					}
					else {
						fwd = i + 1;
					}
				}
				else {
					fwd = i; scale = 1.;
					if (i == 0) {
						bwd = nx - 1;
					}
					else {
						bwd = i - 1;
					}
				}
			}
			dDexxdx_glb(id) = (duxdx_glb_surfYZ(i + 1, j, k) - duxdx_glb_surfYZ(i, j, k)) / dx;
			dDeyydx_glb(id) = (Deyy_glb(fwd, j, k) - Deyy_glb(bwd, j, k)) / dx / scale;
			dDezzdx_glb(id) = (Dezz_glb(fwd, j, k) - Dezz_glb(bwd, j, k)) / dx / scale;

			dDeyzdx_glb(id) = (Deyz_glb(fwd, j, k) - Deyz_glb(bwd, j, k)) / dx / scale;
			dDexzdx_glb(id) = (Dexz_glb(fwd, j, k) - Dexz_glb(bwd, j, k)) / dx / scale;
			dDexydx_glb(id) = (Dexy_glb(fwd, j, k) - Dexy_glb(bwd, j, k)) / dx / scale;

			//------------------spatial derivative along Y--------------//
			if (stress_free_surfXZ(i, j, k) == stress_free_surfXZ(i, j + 1, k)) {
				if (stress_free_surfXZ(i, j, k) == true) {
					bwd = j; fwd = j; scale = 1.;
				}
				else {
					scale = 2.;
					if (j == 0) {
						bwd = ny - 1;
					}
					else {
						bwd = j - 1;
					}
					if (j == ny - 1) {
						fwd = 0;
					}
					else {
						fwd = j + 1;
					}
				}
			}
			else {
				if (stress_free_surfXZ(i, j, k) == true) {
					bwd = j; scale = 1.;
					if (j == ny - 1) {
						fwd = 0;
					}
					else {
						fwd = j + 1;
					}
				}
				else {
					fwd = j; scale = 1.;
					if (j == 0) {
						bwd = ny - 1;
					}
					else {
						bwd = j - 1;
					}
				}
			}
			dDexxdy_glb(id) = (Dexx_glb(i, fwd, k) - Dexx_glb(i, bwd, k)) / dy / scale;
			dDeyydy_glb(id) = (duydy_glb_surfXZ(i, j + 1, k) - duydy_glb_surfXZ(i, j, k)) / dy;
			dDezzdy_glb(id) = (Dezz_glb(i, fwd, k) - Dezz_glb(i, bwd, k)) / dy / scale;

			dDeyzdy_glb(id) = (Deyz_glb(i, fwd, k) - Deyz_glb(i, bwd, k)) / dy / scale;
			dDexzdy_glb(id) = (Dexz_glb(i, fwd, k) - Dexz_glb(i, bwd, k)) / dy / scale;
			dDexydy_glb(id) = (Dexy_glb(i, fwd, k) - Dexy_glb(i, bwd, k)) / dy / scale;

			//------------------spatial derivative along Z--------------//
			if (stress_free_surfXY(i, j, k) == stress_free_surfXY(i, j, k + 1)) {
				if (stress_free_surfXY(i, j, k) == true) {
					bwd = k; fwd = k; scale = 1.;
				}
				else {
					scale = 2.;
					if (k == 0) {
						bwd = nz - 1;
					}
					else {
						bwd = k - 1;
					}
					if (k == nz - 1) {
						fwd = 0;
					}
					else {
						fwd = k + 1;
					}
				}
			}
			else {
				if (stress_free_surfXY(i, j, k) == true) {
					bwd = k; scale = 1.;
					if (k == nz - 1) {
						fwd = 0;
					}
					else {
						fwd = k + 1;
					}
				}
				else {
					fwd = k; scale = 1.;
					if (k == 0) {
						bwd = nz - 1;
					}
					else {
						bwd = k - 1;
					}
				}
			}
			dDexxdz_glb(id) = (Dexx_glb(i, j, fwd) - Dexx_glb(i, j, bwd)) / dz / scale;
			dDeyydz_glb(id) = (Deyy_glb(i, j, fwd) - Deyy_glb(i, j, bwd)) / dz / scale;
			dDezzdz_glb(id) = (duzdz_glb_surfXY(i, j, k + 1) - duzdz_glb_surfXY(i, j, k)) / dz;

			dDeyzdz_glb(id) = (Deyz_glb(i, j, fwd) - Deyz_glb(i, j, bwd)) / dz / scale;
			dDexzdz_glb(id) = (Dexz_glb(i, j, fwd) - Dexz_glb(i, j, bwd)) / dz / scale;
			dDexydz_glb(id) = (Dexy_glb(i, j, fwd) - Dexy_glb(i, j, bwd)) / dz / scale;
		}
	}
}

//void elastic_system::get_dv() {
//	double stress1_surfYZ[6], stress2_surfYZ[6];
//	double stress1_surfXZ[6], stress2_surfXZ[6];
//	double stress1_surfXY[6], stress2_surfXY[6];
//	double strain[6], stress[6];
//	unsigned int mat_type;
//	material* mat;
//	double density, damping;
//	long int i, j, k;
//
//	if (pt_glb->if_elasto_backaction_from_pORm == true) {
//#pragma acc parallel default(present)
//		{
//			get_eigenstrain_dynamic();
//		}
//	}
//
//#pragma acc parallel default(present)
//	{
//#pragma acc loop gang vector private(stress1_surfYZ,stress2_surfYZ,stress1_surfXZ,stress2_surfXZ,stress1_surfXY,stress2_surfXY,\
//strain,stress,mat_type,density,damping, mat,i,j,k)
//		for (long int id = 0; id < n; id++) {
//			i = id / (ny * nz);
//			j = (id - i * (ny * nz)) / nz;
//			k = id - i * (ny * nz) - j * nz;
//
//			mat_type = (pt_glb->material_cell(id));
//			if (mat_type == 0) {
//				dvx_glb(id) = 0.;
//				dvy_glb(id) = 0.;
//				dvz_glb(id) = 0.;
//			}
//			else {
//				//----------YZ surface in -X direction-----//
//				if (stress_free_surfYZ(i, j, k) == true) {
//					stress1_surfYZ[0] = 0.; stress1_surfYZ[1] = 0.; stress1_surfYZ[2] = 0.;
//					stress1_surfYZ[3] = 0.; stress1_surfYZ[4] = 0.; stress1_surfYZ[5] = 0.;
//				}
//				else {
//					mat_type = (pt_glb->material_cell(i, j, k));
//					mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//					strain[0] = duxdx_glb_surfYZ(i, j, k) - Dexx0_glb(i, j, k);
//					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
//					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);
//					strain[3] = (Deyz_glb(i, j, k) - Deyz0_glb(i, j, k)) * 2.;
//					strain[4] = ((duxdz_glb_surfXY(i, j, k) + duxdz_glb_surfXY(i, j, k + 1) + duzdx_glb_surfYZ(i, j, k) * 2.) / 4. \
//						- Dexz0_glb(i, j, k)) * 2.;
//					strain[5] = ((duxdy_glb_surfXZ(i, j, k) + duxdy_glb_surfXZ(i, j + 1, k) + duydx_glb_surfYZ(i, j, k) * 2.) / 4. \
//						- Dexy0_glb(i, j, k)) * 2.;
//
//					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress1_surfYZ);
//
//					if (i == 0) {
//						mat_type = (pt_glb->material_cell(nx - 1, j, k));
//						mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//						strain[0] = duxdx_glb_surfYZ(i, j, k) - Dexx0_glb(nx - 1, j, k);
//						strain[1] = Deyy_glb(nx - 1, j, k) - Deyy0_glb(nx - 1, j, k);
//						strain[2] = Dezz_glb(nx - 1, j, k) - Dezz0_glb(nx - 1, j, k);
//						strain[3] = (Deyz_glb(nx - 1, j, k) - Deyz0_glb(nx - 1, j, k)) * 2.;
//						strain[4] = ((duxdz_glb_surfXY(nx - 1, j, k) + duxdz_glb_surfXY(nx - 1, j, k + 1) + duzdx_glb_surfYZ(i, j, k) * 2.) / 4. \
//							- Dexz0_glb(nx - 1, j, k)) * 2.;
//						strain[5] = ((duxdy_glb_surfXZ(nx - 1, j, k) + duxdy_glb_surfXZ(nx - 1, j + 1, k) + duydx_glb_surfYZ(i, j, k) * 2.) / 4. \
//							- Dexy0_glb(nx - 1, j, k)) * 2.;
//
//						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
//					}
//					else {
//						mat_type = (pt_glb->material_cell(i - 1, j, k));
//						mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//						strain[0] = duxdx_glb_surfYZ(i, j, k) - Dexx0_glb(i - 1, j, k);
//						strain[1] = Deyy_glb(i - 1, j, k) - Deyy0_glb(i - 1, j, k);
//						strain[2] = Dezz_glb(i - 1, j, k) - Dezz0_glb(i - 1, j, k);
//						strain[3] = (Deyz_glb(i - 1, j, k) - Deyz0_glb(i - 1, j, k)) * 2.;
//						strain[4] = ((duxdz_glb_surfXY(i - 1, j, k) + duxdz_glb_surfXY(i - 1, j, k + 1) + duzdx_glb_surfYZ(i, j, k) * 2.) / 4. \
//							- Dexz0_glb(i - 1, j, k)) * 2.;
//						strain[5] = ((duxdy_glb_surfXZ(i - 1, j, k) + duxdy_glb_surfXZ(i - 1, j + 1, k) + duydx_glb_surfYZ(i, j, k) * 2.) / 4. \
//							- Dexy0_glb(i - 1, j, k)) * 2.;
//
//						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
//					}
//					stress1_surfYZ[0] = (stress1_surfYZ[0] + stress[0]) / 2.;
//					stress1_surfYZ[1] = (stress1_surfYZ[1] + stress[1]) / 2.;
//					stress1_surfYZ[2] = (stress1_surfYZ[2] + stress[2]) / 2.;
//					stress1_surfYZ[3] = (stress1_surfYZ[3] + stress[3]) / 2.;
//					stress1_surfYZ[4] = (stress1_surfYZ[4] + stress[4]) / 2.;
//					stress1_surfYZ[5] = (stress1_surfYZ[5] + stress[5]) / 2.;
//				}
//
//				//----------YZ surface in +X direction-----//
//				if (stress_free_surfYZ(i + 1, j, k) == true) {
//					stress2_surfYZ[0] = 0.; stress2_surfYZ[1] = 0.; stress2_surfYZ[2] = 0.;
//					stress2_surfYZ[3] = 0.; stress2_surfYZ[4] = 0.; stress2_surfYZ[5] = 0.;
//				}
//				else {
//					mat_type = (pt_glb->material_cell(i, j, k));
//					mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//					strain[0] = duxdx_glb_surfYZ(i + 1, j, k) - Dexx0_glb(i, j, k);
//					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
//					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);
//					strain[3] = (Deyz_glb(i, j, k) - Deyz0_glb(i, j, k)) * 2.;
//					strain[4] = ((duxdz_glb_surfXY(i, j, k) + duxdz_glb_surfXY(i, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
//						- Dexz0_glb(i, j, k)) * 2.;
//					strain[5] = ((duxdy_glb_surfXZ(i, j, k) + duxdy_glb_surfXZ(i, j + 1, k) + duydx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
//						- Dexy0_glb(i, j, k)) * 2.;
//
//					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress2_surfYZ);
//
//					if (i == nx - 1) {
//						mat_type = (pt_glb->material_cell(0, j, k));
//						mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//						strain[0] = duxdx_glb_surfYZ(i + 1, j, k) - Dexx0_glb(0, j, k);
//						strain[1] = Deyy_glb(0, j, k) - Deyy0_glb(0, j, k);
//						strain[2] = Dezz_glb(0, j, k) - Dezz0_glb(0, j, k);
//						strain[3] = (Deyz_glb(0, j, k) - Deyz0_glb(0, j, k)) * 2.;
//						strain[4] = ((duxdz_glb_surfXY(0, j, k) + duxdz_glb_surfXY(0, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
//							- Dexz0_glb(0, j, k)) * 2.;
//						strain[5] = ((duxdy_glb_surfXZ(0, j, k) + duxdy_glb_surfXZ(0, j + 1, k) + duydx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
//							- Dexy0_glb(0, j, k)) * 2.;
//
//						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
//					}
//					else {
//						mat_type = (pt_glb->material_cell(i + 1, j, k));
//						mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//						strain[0] = duxdx_glb_surfYZ(i + 1, j, k) - Dexx0_glb(i + 1, j, k);
//						strain[1] = Deyy_glb(i + 1, j, k) - Deyy0_glb(i + 1, j, k);
//						strain[2] = Dezz_glb(i + 1, j, k) - Dezz0_glb(i + 1, j, k);
//						strain[3] = (Deyz_glb(i + 1, j, k) - Deyz0_glb(i + 1, j, k)) * 2.;
//						strain[4] = ((duxdz_glb_surfXY(i + 1, j, k) + duxdz_glb_surfXY(i + 1, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
//							- Dexz0_glb(i + 1, j, k)) * 2.;
//						strain[5] = ((duxdy_glb_surfXZ(i + 1, j, k) + duxdy_glb_surfXZ(i + 1, j + 1, k) + duydx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
//							- Dexy0_glb(i + 1, j, k)) * 2.;
//
//						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
//					}
//					stress2_surfYZ[0] = (stress2_surfYZ[0] + stress[0]) / 2.;
//					stress2_surfYZ[1] = (stress2_surfYZ[1] + stress[1]) / 2.;
//					stress2_surfYZ[2] = (stress2_surfYZ[2] + stress[2]) / 2.;
//					stress2_surfYZ[3] = (stress2_surfYZ[3] + stress[3]) / 2.;
//					stress2_surfYZ[4] = (stress2_surfYZ[4] + stress[4]) / 2.;
//					stress2_surfYZ[5] = (stress2_surfYZ[5] + stress[5]) / 2.;
//				}
//				//----------END YZ surface in +-X direction-----//
//
//				//----------XZ surface in -Y direction-----//
//				if (stress_free_surfXZ(i, j, k) == true) {
//					stress1_surfXZ[0] = 0.; stress1_surfXZ[1] = 0.; stress1_surfXZ[2] = 0.;
//					stress1_surfXZ[3] = 0.; stress1_surfXZ[4] = 0.; stress1_surfXZ[5] = 0.;
//				}
//				else {
//					mat_type = (pt_glb->material_cell(i, j, k));
//					mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
//					strain[1] = duydy_glb_surfXZ(i, j, k) - Deyy0_glb(i, j, k);
//					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);
//
//					strain[3] = ((duydz_glb_surfXY(i, j, k) + duydz_glb_surfXY(i, j, k + 1) + \
//						duzdy_glb_surfXZ(i, j, k) * 2.) / 4. - Deyz0_glb(i, j, k)) * 2.;
//					strain[4] = (Dexz_glb(i, j, k) - Dexz0_glb(i, j, k)) * 2.;
//					strain[5] = ((duxdy_glb_surfXZ(i, j, k) * 2. + \
//						duydx_glb_surfYZ(i, j, k) + duydx_glb_surfYZ(i + 1, j, k)) / 4. - Dexy0_glb(i, j, k)) * 2.;
//
//					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress1_surfXZ);
//
//					if (j == 0) {
//						mat_type = (pt_glb->material_cell(i, ny - 1, k));
//						mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//						strain[0] = Dexx_glb(i, ny - 1, k) - Dexx0_glb(i, ny - 1, k);
//						strain[1] = duydy_glb_surfXZ(i, j, k) - Deyy0_glb(i, ny - 1, k);
//						strain[2] = Dezz_glb(i, ny - 1, k) - Dezz0_glb(i, ny - 1, k);
//
//						strain[3] = ((duydz_glb_surfXY(i, ny - 1, k) + duydz_glb_surfXY(i, ny - 1, k + 1) + \
//							duzdy_glb_surfXZ(i, j, k) * 2.) / 4. - Deyz0_glb(i, ny - 1, k)) * 2.;
//
//						strain[4] = (Dexz_glb(i, ny - 1, k) - Dexz0_glb(i, ny - 1, k)) * 2.;
//
//						strain[5] = ((duxdy_glb_surfXZ(i, j, k) * 2. + \
//							duydx_glb_surfYZ(i, ny - 1, k) + duydx_glb_surfYZ(i + 1, ny - 1, k)) / 4. - Dexy0_glb(i, ny - 1, k)) * 2.;
//
//						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
//					}
//					else {
//						mat_type = (pt_glb->material_cell(i, j - 1, k));
//						mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//						strain[0] = Dexx_glb(i, j - 1, k) - Dexx0_glb(i, j - 1, k);
//						strain[1] = duydy_glb_surfXZ(i, j, k) - Deyy0_glb(i, j - 1, k);
//						strain[2] = Dezz_glb(i, j - 1, k) - Dezz0_glb(i, j - 1, k);
//
//						strain[3] = ((duydz_glb_surfXY(i, j - 1, k) + duydz_glb_surfXY(i, j - 1, k + 1) + \
//							duzdy_glb_surfXZ(i, j, k) * 2.) / 4. - Deyz0_glb(i, j - 1, k)) * 2.;
//
//						strain[4] = (Dexz_glb(i, j - 1, k) - Dexz0_glb(i, j - 1, k)) * 2.;
//
//						strain[5] = ((duxdy_glb_surfXZ(i, j, k) * 2. + \
//							duydx_glb_surfYZ(i, j - 1, k) + duydx_glb_surfYZ(i + 1, j - 1, k)) / 4. - Dexy0_glb(i, j - 1, k)) * 2.;
//
//						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
//					}
//					stress1_surfXZ[0] = (stress1_surfXZ[0] + stress[0]) / 2.;
//					stress1_surfXZ[1] = (stress1_surfXZ[1] + stress[1]) / 2.;
//					stress1_surfXZ[2] = (stress1_surfXZ[2] + stress[2]) / 2.;
//					stress1_surfXZ[3] = (stress1_surfXZ[3] + stress[3]) / 2.;
//					stress1_surfXZ[4] = (stress1_surfXZ[4] + stress[4]) / 2.;
//					stress1_surfXZ[5] = (stress1_surfXZ[5] + stress[5]) / 2.;
//				}
//
//				//----------XZ surface in +Y direction-----//
//				if (stress_free_surfXZ(i, j + 1, k) == true) {
//					stress2_surfXZ[0] = 0.; stress2_surfXZ[1] = 0.; stress2_surfXZ[2] = 0.;
//					stress2_surfXZ[3] = 0.; stress2_surfXZ[4] = 0.; stress2_surfXZ[5] = 0.;
//				}
//				else {
//					mat_type = (pt_glb->material_cell(i, j, k));
//					mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
//					strain[1] = duydy_glb_surfXZ(i, j + 1, k) - Deyy0_glb(i, j, k);
//					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);
//
//					strain[3] = ((duydz_glb_surfXY(i, j, k) + duydz_glb_surfXY(i, j, k + 1) + \
//						duzdy_glb_surfXZ(i, j + 1, k) * 2.) / 4. - Deyz0_glb(i, j, k)) * 2.;
//
//					strain[4] = (Dexz_glb(i, j, k) - Dexz0_glb(i, j, k)) * 2.;
//
//					strain[5] = ((duxdy_glb_surfXZ(i, j + 1, k) * 2. + \
//						duydx_glb_surfYZ(i, j, k) + duydx_glb_surfYZ(i + 1, j, k)) / 4. - Dexy0_glb(i, j, k)) * 2.;
//
//					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress2_surfXZ);
//
//					if (j == ny - 1) {
//						mat_type = (pt_glb->material_cell(i, 0, k));
//						mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//						strain[0] = Dexx_glb(i, 0, k) - Dexx0_glb(i, 0, k);
//						strain[1] = duydy_glb_surfXZ(i, j + 1, k) - Deyy0_glb(i, 0, k);
//						strain[2] = Dezz_glb(i, 0, k) - Dezz0_glb(i, 0, k);
//
//						strain[3] = ((duydz_glb_surfXY(i, 0, k) + duydz_glb_surfXY(i, 0, k + 1) + \
//							duzdy_glb_surfXZ(i, j + 1, k) * 2.) / 4. - Deyz0_glb(i, 0, k)) * 2.;
//
//						strain[4] = (Dexz_glb(i, 0, k) - Dexz0_glb(i, 0, k)) * 2.;
//
//						strain[5] = ((duxdy_glb_surfXZ(i, j + 1, k) * 2. + \
//							duydx_glb_surfYZ(i, 0, k) + duydx_glb_surfYZ(i + 1, 0, k)) / 4. - Dexy0_glb(i, 0, k)) * 2.;
//
//						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
//					}
//					else {
//						mat_type = (pt_glb->material_cell(i, j + 1, k));
//						mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//						strain[0] = Dexx_glb(i, j + 1, k) - Dexx0_glb(i, j + 1, k);
//						strain[1] = duydy_glb_surfXZ(i, j + 1, k) - Deyy0_glb(i, j + 1, k);
//						strain[2] = Dezz_glb(i, j + 1, k) - Dezz0_glb(i, j + 1, k);
//
//						strain[3] = ((duydz_glb_surfXY(i, j + 1, k) + duydz_glb_surfXY(i, j + 1, k + 1) + \
//							duzdy_glb_surfXZ(i, j + 1, k) * 2.) / 4. - Deyz0_glb(i, j + 1, k)) * 2.;
//
//						strain[4] = (Dexz_glb(i, j + 1, k) - Dexz0_glb(i, j + 1, k)) * 2.;
//
//						strain[5] = ((duxdy_glb_surfXZ(i, j + 1, k) * 2. + \
//							duydx_glb_surfYZ(i, j + 1, k) + duydx_glb_surfYZ(i + 1, j + 1, k)) / 4. - Dexy0_glb(i, j + 1, k)) * 2.;
//
//						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
//					}
//					stress2_surfXZ[0] = (stress2_surfXZ[0] + stress[0]) / 2.;
//					stress2_surfXZ[1] = (stress2_surfXZ[1] + stress[1]) / 2.;
//					stress2_surfXZ[2] = (stress2_surfXZ[2] + stress[2]) / 2.;
//					stress2_surfXZ[3] = (stress2_surfXZ[3] + stress[3]) / 2.;
//					stress2_surfXZ[4] = (stress2_surfXZ[4] + stress[4]) / 2.;
//					stress2_surfXZ[5] = (stress2_surfXZ[5] + stress[5]) / 2.;
//				}
//				//----------END XZ surface in +-Y direction-----//
//
//				//----------XY surface in -Z direction-----//
//				if (stress_free_surfXY(i, j, k) == true) {
//					stress1_surfXY[0] = 0.; stress1_surfXY[1] = 0.; stress1_surfXY[2] = 0.;
//					stress1_surfXY[3] = 0.; stress1_surfXY[4] = 0.; stress1_surfXY[5] = 0.;
//				}
//				else {
//					mat_type = (pt_glb->material_cell(i, j, k));
//					mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
//					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
//					strain[2] = duzdz_glb_surfXY(i, j, k) - Dezz0_glb(i, j, k);
//
//					strain[3] = ((duydz_glb_surfXY(i, j, k) * 2. + \
//						duzdy_glb_surfXZ(i, j, k) + duzdy_glb_surfXZ(i, j + 1, k)) / 4. - Deyz0_glb(i, j, k)) * 2.;
//
//					strain[4] = ((duxdz_glb_surfXY(i, j, k) * 2. + \
//						duzdx_glb_surfYZ(i, j, k) + duzdx_glb_surfYZ(i + 1, j, k)) / 4. - Dexz0_glb(i, j, k)) * 2.;
//
//					strain[5] = (Dexy_glb(i, j, k) - Dexy0_glb(i, j, k)) * 2.;
//
//					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress1_surfXY);
//
//					if (k == 0) {
//						mat_type = (pt_glb->material_cell(i, j, nz - 1));
//						mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//						strain[0] = Dexx_glb(i, j, nz - 1) - Dexx0_glb(i, j, nz - 1);
//						strain[1] = Deyy_glb(i, j, nz - 1) - Deyy0_glb(i, j, nz - 1);
//						strain[2] = duzdz_glb_surfXY(i, j, k) - Dezz0_glb(i, j, nz - 1);
//
//						strain[3] = ((duydz_glb_surfXY(i, j, k) * 2. + \
//							duzdy_glb_surfXZ(i, j, nz - 1) + duzdy_glb_surfXZ(i, j + 1, nz - 1)) / 4. - Deyz0_glb(i, j, nz - 1)) * 2.;
//
//						strain[4] = ((duxdz_glb_surfXY(i, j, k) * 2. + \
//							duzdx_glb_surfYZ(i, j, nz - 1) + duzdx_glb_surfYZ(i + 1, j, nz - 1)) / 4. - Dexz0_glb(i, j, nz - 1)) * 2.;
//
//						strain[5] = (Dexy_glb(i, j, nz - 1) - Dexy0_glb(i, j, nz - 1)) * 2.;
//
//						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
//					}
//					else {
//						mat_type = (pt_glb->material_cell(i, j, k - 1));
//						mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//						strain[0] = Dexx_glb(i, j, k - 1) - Dexx0_glb(i, j, k - 1);
//						strain[1] = Deyy_glb(i, j, k - 1) - Deyy0_glb(i, j, k - 1);
//						strain[2] = duzdz_glb_surfXY(i, j, k) - Dezz0_glb(i, j, k - 1);
//
//						strain[3] = ((duydz_glb_surfXY(i, j, k) * 2. + \
//							duzdy_glb_surfXZ(i, j, k - 1) + duzdy_glb_surfXZ(i, j + 1, k - 1)) / 4. - Deyz0_glb(i, j, k - 1)) * 2.;
//
//						strain[4] = ((duxdz_glb_surfXY(i, j, k) * 2. + \
//							duzdx_glb_surfYZ(i, j, k - 1) + duzdx_glb_surfYZ(i + 1, j, k - 1)) / 4. - Dexz0_glb(i, j, k - 1)) * 2.;
//
//						strain[5] = (Dexy_glb(i, j, k - 1) - Dexy0_glb(i, j, k - 1)) * 2.;
//
//						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
//					}
//					stress1_surfXY[0] = (stress1_surfXY[0] + stress[0]) / 2.;
//					stress1_surfXY[1] = (stress1_surfXY[1] + stress[1]) / 2.;
//					stress1_surfXY[2] = (stress1_surfXY[2] + stress[2]) / 2.;
//					stress1_surfXY[3] = (stress1_surfXY[3] + stress[3]) / 2.;
//					stress1_surfXY[4] = (stress1_surfXY[4] + stress[4]) / 2.;
//					stress1_surfXY[5] = (stress1_surfXY[5] + stress[5]) / 2.;
//				}
//
//				//----------XY surface in +Z direction-----//
//				if (stress_free_surfXY(i, j, k + 1) == true) {
//					stress2_surfXY[0] = 0.; stress2_surfXY[1] = 0.; stress2_surfXY[2] = 0.;
//					stress2_surfXY[3] = 0.; stress2_surfXY[4] = 0.; stress2_surfXY[5] = 0.;
//				}
//				else {
//					mat_type = (pt_glb->material_cell(i, j, k));
//					mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
//					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
//					strain[2] = duzdz_glb_surfXY(i, j, k + 1) - Dezz0_glb(i, j, k);
//
//					strain[3] = ((duydz_glb_surfXY(i, j, k + 1) * 2. + \
//						duzdy_glb_surfXZ(i, j, k) + duzdy_glb_surfXZ(i, j + 1, k)) / 4. - Deyz0_glb(i, j, k)) * 2.;
//
//					strain[4] = ((duxdz_glb_surfXY(i, j, k + 1) * 2. + \
//						duzdx_glb_surfYZ(i, j, k) + duzdx_glb_surfYZ(i + 1, j, k)) / 4. - Dexz0_glb(i, j, k)) * 2.;
//
//					strain[5] = (Dexy_glb(i, j, k) - Dexy0_glb(i, j, k)) * 2.;
//
//					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress2_surfXY);
//
//					if (k == ny - 1) {
//						mat_type = (pt_glb->material_cell(i, j, 0));
//						mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//						strain[0] = Dexx_glb(i, j, 0) - Dexx0_glb(i, j, 0);
//						strain[1] = Deyy_glb(i, j, 0) - Deyy0_glb(i, j, 0);
//						strain[2] = duzdz_glb_surfXY(i, j, k + 1) - Dezz0_glb(i, j, 0);
//
//						strain[3] = ((duydz_glb_surfXY(i, j, k + 1) * 2. + \
//							duzdy_glb_surfXZ(i, j, 0) + duzdy_glb_surfXZ(i, j + 1, 0)) / 4. - Deyz0_glb(i, j, 0)) * 2.;
//
//						strain[4] = ((duxdz_glb_surfXY(i, j, k + 1) * 2. + \
//							duzdx_glb_surfYZ(i, j, 0) + duzdx_glb_surfYZ(i + 1, j, 0)) / 4. - Dexz0_glb(i, j, 0)) * 2.;
//
//						strain[5] = (Dexy_glb(i, j, 0) - Dexy0_glb(i, j, 0)) * 2.;
//
//						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
//					}
//					else {
//						mat_type = (pt_glb->material_cell(i, j, k + 1));
//						mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//						strain[0] = Dexx_glb(i, j, k + 1) - Dexx0_glb(i, j, k + 1);
//						strain[1] = Deyy_glb(i, j, k + 1) - Deyy0_glb(i, j, k + 1);
//						strain[2] = duzdz_glb_surfXY(i, j, k + 1) - Dezz0_glb(i, j, k + 1);
//
//						strain[3] = ((duydz_glb_surfXY(i, j, k + 1) * 2. + \
//							duzdy_glb_surfXZ(i, j, k + 1) + duzdy_glb_surfXZ(i, j + 1, k + 1)) / 4. - Deyz0_glb(i, j, k + 1)) * 2.;
//
//						strain[4] = ((duxdz_glb_surfXY(i, j, k + 1) * 2. + \
//							duzdx_glb_surfYZ(i, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k + 1)) / 4. - Dexz0_glb(i, j, k + 1)) * 2.;
//
//						strain[5] = (Dexy_glb(i, j, k + 1) - Dexy0_glb(i, j, k + 1)) * 2.;
//
//						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
//					}
//					stress2_surfXY[0] = (stress2_surfXY[0] + stress[0]) / 2.;
//					stress2_surfXY[1] = (stress2_surfXY[1] + stress[1]) / 2.;
//					stress2_surfXY[2] = (stress2_surfXY[2] + stress[2]) / 2.;
//					stress2_surfXY[3] = (stress2_surfXY[3] + stress[3]) / 2.;
//					stress2_surfXY[4] = (stress2_surfXY[4] + stress[4]) / 2.;
//					stress2_surfXY[5] = (stress2_surfXY[5] + stress[5]) / 2.;
//				}
//				//----------END XY surface in +-Z direction-----//
//
//				mat_type = (pt_glb->material_cell(i, j, k));
//				mat = &(pt_glb->material_parameters[(mat_type)-1]);
//
//				density = pt_glb->dt / mat->density;
//				damping = -density * mat->elast_mass_damping;
//
//				dvx_glb(i, j, k) = ((stress2_surfYZ[0] - stress1_surfYZ[0]) / dx + \
//					(stress2_surfXZ[5] - stress1_surfXZ[5]) / dy + \
//					(stress2_surfXY[4] - stress1_surfXY[4]) / dz) * density + damping * vx_glb(i, j, k);
//
//				dvy_glb(i, j, k) = ((stress2_surfYZ[5] - stress1_surfYZ[5]) / dx + \
//					(stress2_surfXZ[1] - stress1_surfXZ[1]) / dy + \
//					(stress2_surfXY[3] - stress1_surfXY[3]) / dz) * density + damping * vy_glb(i, j, k);
//
//				dvz_glb(i, j, k) = ((stress2_surfYZ[4] - stress1_surfYZ[4]) / dx + \
//					(stress2_surfXZ[3] - stress1_surfXZ[3]) / dy + \
//					(stress2_surfXY[2] - stress1_surfXY[2]) / dz) * density + damping * vz_glb(i, j, k);
//			}
//		}
//	}
//}

void elastic_system::get_dv_RK1() {
	double stress1_surfYZ[6], stress2_surfYZ[6];
	double stress1_surfXZ[6], stress2_surfXZ[6];
	double stress1_surfXY[6], stress2_surfXY[6];
	double strain[6], stress[6];
	unsigned int mat_type;
	material* mat;
	double density, mass_damping, stiff_damping;
	long int i, j, k;
	double temporal_var;

	if (pt_glb->if_elasto_backaction_from_pORm == true) {
#pragma acc parallel default(present) async(1)
		{
			get_eigenstrain_dynamic();
		}
	}

#pragma acc wait(1) async(2)
#pragma acc parallel default(present) async(2)
	{
#pragma acc loop gang vector private(stress1_surfYZ,stress2_surfYZ,stress1_surfXZ,stress2_surfXZ,stress1_surfXY,stress2_surfXY,\
strain,stress,mat_type,density,mass_damping, stiff_damping, mat,i,j,k,temporal_var)
		for (long int id = 0; id < n; id++) {
			i = id / (ny * nz);
			j = (id - i * (ny * nz)) / nz;
			k = id - i * (ny * nz) - j * nz;

			mat_type = (pt_glb->material_cell(id));
			if (mat_type == 0) {
				dvx_glb_rk1(id) = 0.;
				dvy_glb_rk1(id) = 0.;
				dvz_glb_rk1(id) = 0.;
			}
			else {
				//----------YZ surface in -X direction-----//
				if (stress_free_surfYZ(i, j, k) == true) {
					stress1_surfYZ[0] = 0.; stress1_surfYZ[1] = 0.; stress1_surfYZ[2] = 0.;
					stress1_surfYZ[3] = 0.; stress1_surfYZ[4] = 0.; stress1_surfYZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = duxdx_glb_surfYZ(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);
					strain[3] = (Deyz_glb(i, j, k) - Deyz0_glb(i, j, k)) * 2.;
					strain[4] = ((duxdz_glb_surfXY(i, j, k) + duxdz_glb_surfXY(i, j, k + 1) + duzdx_glb_surfYZ(i, j, k) * 2.) / 4. \
						- Dexz0_glb(i, j, k)) * 2.;
					strain[5] = ((duxdy_glb_surfXZ(i, j, k) + duxdy_glb_surfXZ(i, j + 1, k) + duydx_glb_surfYZ(i, j, k) * 2.) / 4. \
						- Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress1_surfYZ);

					if (i == 0) {
						mat_type = (pt_glb->material_cell(nx - 1, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i, j, k) - Dexx0_glb(nx - 1, j, k);
						strain[1] = Deyy_glb(nx - 1, j, k) - Deyy0_glb(nx - 1, j, k);
						strain[2] = Dezz_glb(nx - 1, j, k) - Dezz0_glb(nx - 1, j, k);
						strain[3] = (Deyz_glb(nx - 1, j, k) - Deyz0_glb(nx - 1, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(nx - 1, j, k) + duxdz_glb_surfXY(nx - 1, j, k + 1) + duzdx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexz0_glb(nx - 1, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(nx - 1, j, k) + duxdy_glb_surfXZ(nx - 1, j + 1, k) + duydx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexy0_glb(nx - 1, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i - 1, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i, j, k) - Dexx0_glb(i - 1, j, k);
						strain[1] = Deyy_glb(i - 1, j, k) - Deyy0_glb(i - 1, j, k);
						strain[2] = Dezz_glb(i - 1, j, k) - Dezz0_glb(i - 1, j, k);
						strain[3] = (Deyz_glb(i - 1, j, k) - Deyz0_glb(i - 1, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(i - 1, j, k) + duxdz_glb_surfXY(i - 1, j, k + 1) + duzdx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexz0_glb(i - 1, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(i - 1, j, k) + duxdy_glb_surfXZ(i - 1, j + 1, k) + duydx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexy0_glb(i - 1, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress1_surfYZ[0] = (stress1_surfYZ[0] + stress[0]) / 2.;
					stress1_surfYZ[1] = (stress1_surfYZ[1] + stress[1]) / 2.;
					stress1_surfYZ[2] = (stress1_surfYZ[2] + stress[2]) / 2.;
					stress1_surfYZ[3] = (stress1_surfYZ[3] + stress[3]) / 2.;
					stress1_surfYZ[4] = (stress1_surfYZ[4] + stress[4]) / 2.;
					stress1_surfYZ[5] = (stress1_surfYZ[5] + stress[5]) / 2.;
				}

				//----------YZ surface in +X direction-----//
				if (stress_free_surfYZ(i + 1, j, k) == true) {
					stress2_surfYZ[0] = 0.; stress2_surfYZ[1] = 0.; stress2_surfYZ[2] = 0.;
					stress2_surfYZ[3] = 0.; stress2_surfYZ[4] = 0.; stress2_surfYZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = duxdx_glb_surfYZ(i + 1, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);
					strain[3] = (Deyz_glb(i, j, k) - Deyz0_glb(i, j, k)) * 2.;
					strain[4] = ((duxdz_glb_surfXY(i, j, k) + duxdz_glb_surfXY(i, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
						- Dexz0_glb(i, j, k)) * 2.;
					strain[5] = ((duxdy_glb_surfXZ(i, j, k) + duxdy_glb_surfXZ(i, j + 1, k) + duydx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
						- Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress2_surfYZ);

					if (i == nx - 1) {
						mat_type = (pt_glb->material_cell(0, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i + 1, j, k) - Dexx0_glb(0, j, k);
						strain[1] = Deyy_glb(0, j, k) - Deyy0_glb(0, j, k);
						strain[2] = Dezz_glb(0, j, k) - Dezz0_glb(0, j, k);
						strain[3] = (Deyz_glb(0, j, k) - Deyz0_glb(0, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(0, j, k) + duxdz_glb_surfXY(0, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexz0_glb(0, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(0, j, k) + duxdy_glb_surfXZ(0, j + 1, k) + duydx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexy0_glb(0, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i + 1, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i + 1, j, k) - Dexx0_glb(i + 1, j, k);
						strain[1] = Deyy_glb(i + 1, j, k) - Deyy0_glb(i + 1, j, k);
						strain[2] = Dezz_glb(i + 1, j, k) - Dezz0_glb(i + 1, j, k);
						strain[3] = (Deyz_glb(i + 1, j, k) - Deyz0_glb(i + 1, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(i + 1, j, k) + duxdz_glb_surfXY(i + 1, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexz0_glb(i + 1, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(i + 1, j, k) + duxdy_glb_surfXZ(i + 1, j + 1, k) + duydx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexy0_glb(i + 1, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress2_surfYZ[0] = (stress2_surfYZ[0] + stress[0]) / 2.;
					stress2_surfYZ[1] = (stress2_surfYZ[1] + stress[1]) / 2.;
					stress2_surfYZ[2] = (stress2_surfYZ[2] + stress[2]) / 2.;
					stress2_surfYZ[3] = (stress2_surfYZ[3] + stress[3]) / 2.;
					stress2_surfYZ[4] = (stress2_surfYZ[4] + stress[4]) / 2.;
					stress2_surfYZ[5] = (stress2_surfYZ[5] + stress[5]) / 2.;
				}
				//----------END YZ surface in +-X direction-----//

				//----------XZ surface in -Y direction-----//
				if (stress_free_surfXZ(i, j, k) == true) {
					stress1_surfXZ[0] = 0.; stress1_surfXZ[1] = 0.; stress1_surfXZ[2] = 0.;
					stress1_surfXZ[3] = 0.; stress1_surfXZ[4] = 0.; stress1_surfXZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = duydy_glb_surfXZ(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k) + duydz_glb_surfXY(i, j, k + 1) + \
						duzdy_glb_surfXZ(i, j, k) * 2.) / 4. - Deyz0_glb(i, j, k)) * 2.;
					strain[4] = (Dexz_glb(i, j, k) - Dexz0_glb(i, j, k)) * 2.;
					strain[5] = ((duxdy_glb_surfXZ(i, j, k) * 2. + \
						duydx_glb_surfYZ(i, j, k) + duydx_glb_surfYZ(i + 1, j, k)) / 4. - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress1_surfXZ);

					if (j == 0) {
						mat_type = (pt_glb->material_cell(i, ny - 1, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, ny - 1, k) - Dexx0_glb(i, ny - 1, k);
						strain[1] = duydy_glb_surfXZ(i, j, k) - Deyy0_glb(i, ny - 1, k);
						strain[2] = Dezz_glb(i, ny - 1, k) - Dezz0_glb(i, ny - 1, k);

						strain[3] = ((duydz_glb_surfXY(i, ny - 1, k) + duydz_glb_surfXY(i, ny - 1, k + 1) + \
							duzdy_glb_surfXZ(i, j, k) * 2.) / 4. - Deyz0_glb(i, ny - 1, k)) * 2.;

						strain[4] = (Dexz_glb(i, ny - 1, k) - Dexz0_glb(i, ny - 1, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j, k) * 2. + \
							duydx_glb_surfYZ(i, ny - 1, k) + duydx_glb_surfYZ(i + 1, ny - 1, k)) / 4. - Dexy0_glb(i, ny - 1, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j - 1, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j - 1, k) - Dexx0_glb(i, j - 1, k);
						strain[1] = duydy_glb_surfXZ(i, j, k) - Deyy0_glb(i, j - 1, k);
						strain[2] = Dezz_glb(i, j - 1, k) - Dezz0_glb(i, j - 1, k);

						strain[3] = ((duydz_glb_surfXY(i, j - 1, k) + duydz_glb_surfXY(i, j - 1, k + 1) + \
							duzdy_glb_surfXZ(i, j, k) * 2.) / 4. - Deyz0_glb(i, j - 1, k)) * 2.;

						strain[4] = (Dexz_glb(i, j - 1, k) - Dexz0_glb(i, j - 1, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j, k) * 2. + \
							duydx_glb_surfYZ(i, j - 1, k) + duydx_glb_surfYZ(i + 1, j - 1, k)) / 4. - Dexy0_glb(i, j - 1, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress1_surfXZ[0] = (stress1_surfXZ[0] + stress[0]) / 2.;
					stress1_surfXZ[1] = (stress1_surfXZ[1] + stress[1]) / 2.;
					stress1_surfXZ[2] = (stress1_surfXZ[2] + stress[2]) / 2.;
					stress1_surfXZ[3] = (stress1_surfXZ[3] + stress[3]) / 2.;
					stress1_surfXZ[4] = (stress1_surfXZ[4] + stress[4]) / 2.;
					stress1_surfXZ[5] = (stress1_surfXZ[5] + stress[5]) / 2.;
				}

				//----------XZ surface in +Y direction-----//
				if (stress_free_surfXZ(i, j + 1, k) == true) {
					stress2_surfXZ[0] = 0.; stress2_surfXZ[1] = 0.; stress2_surfXZ[2] = 0.;
					stress2_surfXZ[3] = 0.; stress2_surfXZ[4] = 0.; stress2_surfXZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = duydy_glb_surfXZ(i, j + 1, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k) + duydz_glb_surfXY(i, j, k + 1) + \
						duzdy_glb_surfXZ(i, j + 1, k) * 2.) / 4. - Deyz0_glb(i, j, k)) * 2.;

					strain[4] = (Dexz_glb(i, j, k) - Dexz0_glb(i, j, k)) * 2.;

					strain[5] = ((duxdy_glb_surfXZ(i, j + 1, k) * 2. + \
						duydx_glb_surfYZ(i, j, k) + duydx_glb_surfYZ(i + 1, j, k)) / 4. - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress2_surfXZ);

					if (j == ny - 1) {
						mat_type = (pt_glb->material_cell(i, 0, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, 0, k) - Dexx0_glb(i, 0, k);
						strain[1] = duydy_glb_surfXZ(i, j + 1, k) - Deyy0_glb(i, 0, k);
						strain[2] = Dezz_glb(i, 0, k) - Dezz0_glb(i, 0, k);

						strain[3] = ((duydz_glb_surfXY(i, 0, k) + duydz_glb_surfXY(i, 0, k + 1) + \
							duzdy_glb_surfXZ(i, j + 1, k) * 2.) / 4. - Deyz0_glb(i, 0, k)) * 2.;

						strain[4] = (Dexz_glb(i, 0, k) - Dexz0_glb(i, 0, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j + 1, k) * 2. + \
							duydx_glb_surfYZ(i, 0, k) + duydx_glb_surfYZ(i + 1, 0, k)) / 4. - Dexy0_glb(i, 0, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j + 1, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j + 1, k) - Dexx0_glb(i, j + 1, k);
						strain[1] = duydy_glb_surfXZ(i, j + 1, k) - Deyy0_glb(i, j + 1, k);
						strain[2] = Dezz_glb(i, j + 1, k) - Dezz0_glb(i, j + 1, k);

						strain[3] = ((duydz_glb_surfXY(i, j + 1, k) + duydz_glb_surfXY(i, j + 1, k + 1) + \
							duzdy_glb_surfXZ(i, j + 1, k) * 2.) / 4. - Deyz0_glb(i, j + 1, k)) * 2.;

						strain[4] = (Dexz_glb(i, j + 1, k) - Dexz0_glb(i, j + 1, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j + 1, k) * 2. + \
							duydx_glb_surfYZ(i, j + 1, k) + duydx_glb_surfYZ(i + 1, j + 1, k)) / 4. - Dexy0_glb(i, j + 1, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress2_surfXZ[0] = (stress2_surfXZ[0] + stress[0]) / 2.;
					stress2_surfXZ[1] = (stress2_surfXZ[1] + stress[1]) / 2.;
					stress2_surfXZ[2] = (stress2_surfXZ[2] + stress[2]) / 2.;
					stress2_surfXZ[3] = (stress2_surfXZ[3] + stress[3]) / 2.;
					stress2_surfXZ[4] = (stress2_surfXZ[4] + stress[4]) / 2.;
					stress2_surfXZ[5] = (stress2_surfXZ[5] + stress[5]) / 2.;
				}
				//----------END XZ surface in +-Y direction-----//

				//----------XY surface in -Z direction-----//
				if (stress_free_surfXY(i, j, k) == true) {
					stress1_surfXY[0] = 0.; stress1_surfXY[1] = 0.; stress1_surfXY[2] = 0.;
					stress1_surfXY[3] = 0.; stress1_surfXY[4] = 0.; stress1_surfXY[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = duzdz_glb_surfXY(i, j, k) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k) * 2. + \
						duzdy_glb_surfXZ(i, j, k) + duzdy_glb_surfXZ(i, j + 1, k)) / 4. - Deyz0_glb(i, j, k)) * 2.;

					strain[4] = ((duxdz_glb_surfXY(i, j, k) * 2. + \
						duzdx_glb_surfYZ(i, j, k) + duzdx_glb_surfYZ(i + 1, j, k)) / 4. - Dexz0_glb(i, j, k)) * 2.;

					strain[5] = (Dexy_glb(i, j, k) - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress1_surfXY);

					if (k == 0) {
						mat_type = (pt_glb->material_cell(i, j, nz - 1));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, nz - 1) - Dexx0_glb(i, j, nz - 1);
						strain[1] = Deyy_glb(i, j, nz - 1) - Deyy0_glb(i, j, nz - 1);
						strain[2] = duzdz_glb_surfXY(i, j, k) - Dezz0_glb(i, j, nz - 1);

						strain[3] = ((duydz_glb_surfXY(i, j, k) * 2. + \
							duzdy_glb_surfXZ(i, j, nz - 1) + duzdy_glb_surfXZ(i, j + 1, nz - 1)) / 4. - Deyz0_glb(i, j, nz - 1)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k) * 2. + \
							duzdx_glb_surfYZ(i, j, nz - 1) + duzdx_glb_surfYZ(i + 1, j, nz - 1)) / 4. - Dexz0_glb(i, j, nz - 1)) * 2.;

						strain[5] = (Dexy_glb(i, j, nz - 1) - Dexy0_glb(i, j, nz - 1)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j, k - 1));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, k - 1) - Dexx0_glb(i, j, k - 1);
						strain[1] = Deyy_glb(i, j, k - 1) - Deyy0_glb(i, j, k - 1);
						strain[2] = duzdz_glb_surfXY(i, j, k) - Dezz0_glb(i, j, k - 1);

						strain[3] = ((duydz_glb_surfXY(i, j, k) * 2. + \
							duzdy_glb_surfXZ(i, j, k - 1) + duzdy_glb_surfXZ(i, j + 1, k - 1)) / 4. - Deyz0_glb(i, j, k - 1)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k) * 2. + \
							duzdx_glb_surfYZ(i, j, k - 1) + duzdx_glb_surfYZ(i + 1, j, k - 1)) / 4. - Dexz0_glb(i, j, k - 1)) * 2.;

						strain[5] = (Dexy_glb(i, j, k - 1) - Dexy0_glb(i, j, k - 1)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress1_surfXY[0] = (stress1_surfXY[0] + stress[0]) / 2.;
					stress1_surfXY[1] = (stress1_surfXY[1] + stress[1]) / 2.;
					stress1_surfXY[2] = (stress1_surfXY[2] + stress[2]) / 2.;
					stress1_surfXY[3] = (stress1_surfXY[3] + stress[3]) / 2.;
					stress1_surfXY[4] = (stress1_surfXY[4] + stress[4]) / 2.;
					stress1_surfXY[5] = (stress1_surfXY[5] + stress[5]) / 2.;
				}

				//----------XY surface in +Z direction-----//
				if (stress_free_surfXY(i, j, k + 1) == true) {
					if (pt_glb->if_gaussian_strain_pulse == true && pt_glb->input_strain_type == 3 && k == source_z) {

						temporal_var = pt_glb->time_device - pt_glb->dt - 5.0 * pt_glb->sigma_gauss;
						if (pt_glb->input_strain_component == 'x') {
							stress2_surfXY[0] = 0.;
							stress2_surfXY[1] = 0.;
							stress2_surfXY[2] = 0.;
							stress2_surfXY[3] = 0.;
							stress2_surfXY[4] = pt_glb->amplitude_gauss * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->sigma_gauss, 2.)));
							stress2_surfXY[5] = 0.;
						}
						else if (pt_glb->input_strain_component == 'y') {
							stress2_surfXY[0] = 0.;
							stress2_surfXY[1] = 0.;
							stress2_surfXY[2] = 0.;
							stress2_surfXY[3] = pt_glb->amplitude_gauss * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->sigma_gauss, 2.)));
							stress2_surfXY[4] = 0.;
							stress2_surfXY[5] = 0.;
						}
						else if (pt_glb->input_strain_component == 'z') {
							stress2_surfXY[0] = 0.; 
							stress2_surfXY[1] = 0.;
							stress2_surfXY[2] = pt_glb->amplitude_gauss * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->sigma_gauss, 2.)));
							stress2_surfXY[3] = 0.;
							stress2_surfXY[4] = 0.; 
							stress2_surfXY[5] = 0.;
						}

					}
					else {
						stress2_surfXY[0] = 0.; stress2_surfXY[1] = 0.; stress2_surfXY[2] = 0.;
						stress2_surfXY[3] = 0.; stress2_surfXY[4] = 0.; stress2_surfXY[5] = 0.;
					}
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = duzdz_glb_surfXY(i, j, k + 1) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k + 1) * 2. + \
						duzdy_glb_surfXZ(i, j, k) + duzdy_glb_surfXZ(i, j + 1, k)) / 4. - Deyz0_glb(i, j, k)) * 2.;

					strain[4] = ((duxdz_glb_surfXY(i, j, k + 1) * 2. + \
						duzdx_glb_surfYZ(i, j, k) + duzdx_glb_surfYZ(i + 1, j, k)) / 4. - Dexz0_glb(i, j, k)) * 2.;

					strain[5] = (Dexy_glb(i, j, k) - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress2_surfXY);

					if (k == ny - 1) {
						mat_type = (pt_glb->material_cell(i, j, 0));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, 0) - Dexx0_glb(i, j, 0);
						strain[1] = Deyy_glb(i, j, 0) - Deyy0_glb(i, j, 0);
						strain[2] = duzdz_glb_surfXY(i, j, k + 1) - Dezz0_glb(i, j, 0);

						strain[3] = ((duydz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdy_glb_surfXZ(i, j, 0) + duzdy_glb_surfXZ(i, j + 1, 0)) / 4. - Deyz0_glb(i, j, 0)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdx_glb_surfYZ(i, j, 0) + duzdx_glb_surfYZ(i + 1, j, 0)) / 4. - Dexz0_glb(i, j, 0)) * 2.;

						strain[5] = (Dexy_glb(i, j, 0) - Dexy0_glb(i, j, 0)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j, k + 1));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, k + 1) - Dexx0_glb(i, j, k + 1);
						strain[1] = Deyy_glb(i, j, k + 1) - Deyy0_glb(i, j, k + 1);
						strain[2] = duzdz_glb_surfXY(i, j, k + 1) - Dezz0_glb(i, j, k + 1);

						strain[3] = ((duydz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdy_glb_surfXZ(i, j, k + 1) + duzdy_glb_surfXZ(i, j + 1, k + 1)) / 4. - Deyz0_glb(i, j, k + 1)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdx_glb_surfYZ(i, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k + 1)) / 4. - Dexz0_glb(i, j, k + 1)) * 2.;

						strain[5] = (Dexy_glb(i, j, k + 1) - Dexy0_glb(i, j, k + 1)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress2_surfXY[0] = (stress2_surfXY[0] + stress[0]) / 2.;
					stress2_surfXY[1] = (stress2_surfXY[1] + stress[1]) / 2.;
					stress2_surfXY[2] = (stress2_surfXY[2] + stress[2]) / 2.;
					stress2_surfXY[3] = (stress2_surfXY[3] + stress[3]) / 2.;
					stress2_surfXY[4] = (stress2_surfXY[4] + stress[4]) / 2.;
					stress2_surfXY[5] = (stress2_surfXY[5] + stress[5]) / 2.;
				}
				//----------END XY surface in +-Z direction-----//

				mat_type = (pt_glb->material_cell(i, j, k));
				mat = &(pt_glb->material_parameters[(mat_type)-1]);

				density = pt_glb->dt / mat->density;
				mass_damping = -pt_glb->dt * mat->elast_mass_damping;
				stiff_damping = mat->elast_stiff_damping / mat->density;

				force_x(i, j, k) = (stress2_surfYZ[0] - stress1_surfYZ[0]) / dx + \
					(stress2_surfXZ[5] - stress1_surfXZ[5]) / dy + \
					(stress2_surfXY[4] - stress1_surfXY[4]) / dz;

				force_y(i, j, k) = (stress2_surfYZ[5] - stress1_surfYZ[5]) / dx + \
					(stress2_surfXZ[1] - stress1_surfXZ[1]) / dy + \
					(stress2_surfXY[3] - stress1_surfXY[3]) / dz;

				force_z(i, j, k) = (stress2_surfYZ[4] - stress1_surfYZ[4]) / dx + \
					(stress2_surfXZ[3] - stress1_surfXZ[3]) / dy + \
					(stress2_surfXY[2] - stress1_surfXY[2]) / dz;

				dvx_glb_rk1(i, j, k) = force_x(i, j, k) * density \
					+ mass_damping * vx_glb(i, j, k) \
					+ stiff_damping * (force_x(i, j, k) - force_x_store(i, j, k));

				dvy_glb_rk1(i, j, k) = force_y(i, j, k) * density \
					+ mass_damping * vy_glb(i, j, k) \
					+ stiff_damping * (force_y(i, j, k) - force_y_store(i, j, k));

				dvz_glb_rk1(i, j, k) = force_z(i, j, k) * density \
					+ mass_damping * vz_glb(i, j, k) \
					+ stiff_damping * (force_z(i, j, k) - force_z_store(i, j, k));

				force_x_store(i, j, k) = force_x(i, j, k);
				force_y_store(i, j, k) = force_y(i, j, k);
				force_z_store(i, j, k) = force_z(i, j, k);
			}
		}
	}
}

void elastic_system::get_dv_RK2() {
	double stress1_surfYZ[6], stress2_surfYZ[6];
	double stress1_surfXZ[6], stress2_surfXZ[6];
	double stress1_surfXY[6], stress2_surfXY[6];
	double strain[6], stress[6];
	unsigned int mat_type;
	material* mat;
	double density, mass_damping, stiff_damping;
	long int i, j, k;
	double temporal_var;

	if (pt_glb->if_elasto_backaction_from_pORm == true) {
#pragma acc parallel default(present) async(1)
		{
			get_eigenstrain_dynamic();
		}
	}

#pragma acc wait(1) async(2)
#pragma acc parallel default(present) async(2)
	{
#pragma acc loop gang vector private(stress1_surfYZ,stress2_surfYZ,stress1_surfXZ,stress2_surfXZ,stress1_surfXY,stress2_surfXY,\
strain,stress,mat_type,density, mass_damping, stiff_damping, mat,i,j,k,temporal_var)
		for (long int id = 0; id < n; id++) {
			i = id / (ny * nz);
			j = (id - i * (ny * nz)) / nz;
			k = id - i * (ny * nz) - j * nz;

			mat_type = (pt_glb->material_cell(id));
			if (mat_type == 0) {
				dvx_glb_rk2(id) = 0.;
				dvy_glb_rk2(id) = 0.;
				dvz_glb_rk2(id) = 0.;
			}
			else {
				//----------YZ surface in -X direction-----//
				if (stress_free_surfYZ(i, j, k) == true) {
					stress1_surfYZ[0] = 0.; stress1_surfYZ[1] = 0.; stress1_surfYZ[2] = 0.;
					stress1_surfYZ[3] = 0.; stress1_surfYZ[4] = 0.; stress1_surfYZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = duxdx_glb_surfYZ(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);
					strain[3] = (Deyz_glb(i, j, k) - Deyz0_glb(i, j, k)) * 2.;
					strain[4] = ((duxdz_glb_surfXY(i, j, k) + duxdz_glb_surfXY(i, j, k + 1) + duzdx_glb_surfYZ(i, j, k) * 2.) / 4. \
						- Dexz0_glb(i, j, k)) * 2.;
					strain[5] = ((duxdy_glb_surfXZ(i, j, k) + duxdy_glb_surfXZ(i, j + 1, k) + duydx_glb_surfYZ(i, j, k) * 2.) / 4. \
						- Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress1_surfYZ);

					if (i == 0) {
						mat_type = (pt_glb->material_cell(nx - 1, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i, j, k) - Dexx0_glb(nx - 1, j, k);
						strain[1] = Deyy_glb(nx - 1, j, k) - Deyy0_glb(nx - 1, j, k);
						strain[2] = Dezz_glb(nx - 1, j, k) - Dezz0_glb(nx - 1, j, k);
						strain[3] = (Deyz_glb(nx - 1, j, k) - Deyz0_glb(nx - 1, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(nx - 1, j, k) + duxdz_glb_surfXY(nx - 1, j, k + 1) + duzdx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexz0_glb(nx - 1, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(nx - 1, j, k) + duxdy_glb_surfXZ(nx - 1, j + 1, k) + duydx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexy0_glb(nx - 1, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i - 1, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i, j, k) - Dexx0_glb(i - 1, j, k);
						strain[1] = Deyy_glb(i - 1, j, k) - Deyy0_glb(i - 1, j, k);
						strain[2] = Dezz_glb(i - 1, j, k) - Dezz0_glb(i - 1, j, k);
						strain[3] = (Deyz_glb(i - 1, j, k) - Deyz0_glb(i - 1, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(i - 1, j, k) + duxdz_glb_surfXY(i - 1, j, k + 1) + duzdx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexz0_glb(i - 1, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(i - 1, j, k) + duxdy_glb_surfXZ(i - 1, j + 1, k) + duydx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexy0_glb(i - 1, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress1_surfYZ[0] = (stress1_surfYZ[0] + stress[0]) / 2.;
					stress1_surfYZ[1] = (stress1_surfYZ[1] + stress[1]) / 2.;
					stress1_surfYZ[2] = (stress1_surfYZ[2] + stress[2]) / 2.;
					stress1_surfYZ[3] = (stress1_surfYZ[3] + stress[3]) / 2.;
					stress1_surfYZ[4] = (stress1_surfYZ[4] + stress[4]) / 2.;
					stress1_surfYZ[5] = (stress1_surfYZ[5] + stress[5]) / 2.;
				}

				//----------YZ surface in +X direction-----//
				if (stress_free_surfYZ(i + 1, j, k) == true) {
					stress2_surfYZ[0] = 0.; stress2_surfYZ[1] = 0.; stress2_surfYZ[2] = 0.;
					stress2_surfYZ[3] = 0.; stress2_surfYZ[4] = 0.; stress2_surfYZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = duxdx_glb_surfYZ(i + 1, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);
					strain[3] = (Deyz_glb(i, j, k) - Deyz0_glb(i, j, k)) * 2.;
					strain[4] = ((duxdz_glb_surfXY(i, j, k) + duxdz_glb_surfXY(i, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
						- Dexz0_glb(i, j, k)) * 2.;
					strain[5] = ((duxdy_glb_surfXZ(i, j, k) + duxdy_glb_surfXZ(i, j + 1, k) + duydx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
						- Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress2_surfYZ);

					if (i == nx - 1) {
						mat_type = (pt_glb->material_cell(0, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i + 1, j, k) - Dexx0_glb(0, j, k);
						strain[1] = Deyy_glb(0, j, k) - Deyy0_glb(0, j, k);
						strain[2] = Dezz_glb(0, j, k) - Dezz0_glb(0, j, k);
						strain[3] = (Deyz_glb(0, j, k) - Deyz0_glb(0, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(0, j, k) + duxdz_glb_surfXY(0, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexz0_glb(0, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(0, j, k) + duxdy_glb_surfXZ(0, j + 1, k) + duydx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexy0_glb(0, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i + 1, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i + 1, j, k) - Dexx0_glb(i + 1, j, k);
						strain[1] = Deyy_glb(i + 1, j, k) - Deyy0_glb(i + 1, j, k);
						strain[2] = Dezz_glb(i + 1, j, k) - Dezz0_glb(i + 1, j, k);
						strain[3] = (Deyz_glb(i + 1, j, k) - Deyz0_glb(i + 1, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(i + 1, j, k) + duxdz_glb_surfXY(i + 1, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexz0_glb(i + 1, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(i + 1, j, k) + duxdy_glb_surfXZ(i + 1, j + 1, k) + duydx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexy0_glb(i + 1, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress2_surfYZ[0] = (stress2_surfYZ[0] + stress[0]) / 2.;
					stress2_surfYZ[1] = (stress2_surfYZ[1] + stress[1]) / 2.;
					stress2_surfYZ[2] = (stress2_surfYZ[2] + stress[2]) / 2.;
					stress2_surfYZ[3] = (stress2_surfYZ[3] + stress[3]) / 2.;
					stress2_surfYZ[4] = (stress2_surfYZ[4] + stress[4]) / 2.;
					stress2_surfYZ[5] = (stress2_surfYZ[5] + stress[5]) / 2.;
				}
				//----------END YZ surface in +-X direction-----//

				//----------XZ surface in -Y direction-----//
				if (stress_free_surfXZ(i, j, k) == true) {
					stress1_surfXZ[0] = 0.; stress1_surfXZ[1] = 0.; stress1_surfXZ[2] = 0.;
					stress1_surfXZ[3] = 0.; stress1_surfXZ[4] = 0.; stress1_surfXZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = duydy_glb_surfXZ(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k) + duydz_glb_surfXY(i, j, k + 1) + \
						duzdy_glb_surfXZ(i, j, k) * 2.) / 4. - Deyz0_glb(i, j, k)) * 2.;
					strain[4] = (Dexz_glb(i, j, k) - Dexz0_glb(i, j, k)) * 2.;
					strain[5] = ((duxdy_glb_surfXZ(i, j, k) * 2. + \
						duydx_glb_surfYZ(i, j, k) + duydx_glb_surfYZ(i + 1, j, k)) / 4. - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress1_surfXZ);

					if (j == 0) {
						mat_type = (pt_glb->material_cell(i, ny - 1, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, ny - 1, k) - Dexx0_glb(i, ny - 1, k);
						strain[1] = duydy_glb_surfXZ(i, j, k) - Deyy0_glb(i, ny - 1, k);
						strain[2] = Dezz_glb(i, ny - 1, k) - Dezz0_glb(i, ny - 1, k);

						strain[3] = ((duydz_glb_surfXY(i, ny - 1, k) + duydz_glb_surfXY(i, ny - 1, k + 1) + \
							duzdy_glb_surfXZ(i, j, k) * 2.) / 4. - Deyz0_glb(i, ny - 1, k)) * 2.;

						strain[4] = (Dexz_glb(i, ny - 1, k) - Dexz0_glb(i, ny - 1, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j, k) * 2. + \
							duydx_glb_surfYZ(i, ny - 1, k) + duydx_glb_surfYZ(i + 1, ny - 1, k)) / 4. - Dexy0_glb(i, ny - 1, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j - 1, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j - 1, k) - Dexx0_glb(i, j - 1, k);
						strain[1] = duydy_glb_surfXZ(i, j, k) - Deyy0_glb(i, j - 1, k);
						strain[2] = Dezz_glb(i, j - 1, k) - Dezz0_glb(i, j - 1, k);

						strain[3] = ((duydz_glb_surfXY(i, j - 1, k) + duydz_glb_surfXY(i, j - 1, k + 1) + \
							duzdy_glb_surfXZ(i, j, k) * 2.) / 4. - Deyz0_glb(i, j - 1, k)) * 2.;

						strain[4] = (Dexz_glb(i, j - 1, k) - Dexz0_glb(i, j - 1, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j, k) * 2. + \
							duydx_glb_surfYZ(i, j - 1, k) + duydx_glb_surfYZ(i + 1, j - 1, k)) / 4. - Dexy0_glb(i, j - 1, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress1_surfXZ[0] = (stress1_surfXZ[0] + stress[0]) / 2.;
					stress1_surfXZ[1] = (stress1_surfXZ[1] + stress[1]) / 2.;
					stress1_surfXZ[2] = (stress1_surfXZ[2] + stress[2]) / 2.;
					stress1_surfXZ[3] = (stress1_surfXZ[3] + stress[3]) / 2.;
					stress1_surfXZ[4] = (stress1_surfXZ[4] + stress[4]) / 2.;
					stress1_surfXZ[5] = (stress1_surfXZ[5] + stress[5]) / 2.;
				}

				//----------XZ surface in +Y direction-----//
				if (stress_free_surfXZ(i, j + 1, k) == true) {
					stress2_surfXZ[0] = 0.; stress2_surfXZ[1] = 0.; stress2_surfXZ[2] = 0.;
					stress2_surfXZ[3] = 0.; stress2_surfXZ[4] = 0.; stress2_surfXZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = duydy_glb_surfXZ(i, j + 1, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k) + duydz_glb_surfXY(i, j, k + 1) + \
						duzdy_glb_surfXZ(i, j + 1, k) * 2.) / 4. - Deyz0_glb(i, j, k)) * 2.;

					strain[4] = (Dexz_glb(i, j, k) - Dexz0_glb(i, j, k)) * 2.;

					strain[5] = ((duxdy_glb_surfXZ(i, j + 1, k) * 2. + \
						duydx_glb_surfYZ(i, j, k) + duydx_glb_surfYZ(i + 1, j, k)) / 4. - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress2_surfXZ);

					if (j == ny - 1) {
						mat_type = (pt_glb->material_cell(i, 0, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, 0, k) - Dexx0_glb(i, 0, k);
						strain[1] = duydy_glb_surfXZ(i, j + 1, k) - Deyy0_glb(i, 0, k);
						strain[2] = Dezz_glb(i, 0, k) - Dezz0_glb(i, 0, k);

						strain[3] = ((duydz_glb_surfXY(i, 0, k) + duydz_glb_surfXY(i, 0, k + 1) + \
							duzdy_glb_surfXZ(i, j + 1, k) * 2.) / 4. - Deyz0_glb(i, 0, k)) * 2.;

						strain[4] = (Dexz_glb(i, 0, k) - Dexz0_glb(i, 0, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j + 1, k) * 2. + \
							duydx_glb_surfYZ(i, 0, k) + duydx_glb_surfYZ(i + 1, 0, k)) / 4. - Dexy0_glb(i, 0, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j + 1, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j + 1, k) - Dexx0_glb(i, j + 1, k);
						strain[1] = duydy_glb_surfXZ(i, j + 1, k) - Deyy0_glb(i, j + 1, k);
						strain[2] = Dezz_glb(i, j + 1, k) - Dezz0_glb(i, j + 1, k);

						strain[3] = ((duydz_glb_surfXY(i, j + 1, k) + duydz_glb_surfXY(i, j + 1, k + 1) + \
							duzdy_glb_surfXZ(i, j + 1, k) * 2.) / 4. - Deyz0_glb(i, j + 1, k)) * 2.;

						strain[4] = (Dexz_glb(i, j + 1, k) - Dexz0_glb(i, j + 1, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j + 1, k) * 2. + \
							duydx_glb_surfYZ(i, j + 1, k) + duydx_glb_surfYZ(i + 1, j + 1, k)) / 4. - Dexy0_glb(i, j + 1, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress2_surfXZ[0] = (stress2_surfXZ[0] + stress[0]) / 2.;
					stress2_surfXZ[1] = (stress2_surfXZ[1] + stress[1]) / 2.;
					stress2_surfXZ[2] = (stress2_surfXZ[2] + stress[2]) / 2.;
					stress2_surfXZ[3] = (stress2_surfXZ[3] + stress[3]) / 2.;
					stress2_surfXZ[4] = (stress2_surfXZ[4] + stress[4]) / 2.;
					stress2_surfXZ[5] = (stress2_surfXZ[5] + stress[5]) / 2.;
				}
				//----------END XZ surface in +-Y direction-----//

				//----------XY surface in -Z direction-----//
				if (stress_free_surfXY(i, j, k) == true) {
					stress1_surfXY[0] = 0.; stress1_surfXY[1] = 0.; stress1_surfXY[2] = 0.;
					stress1_surfXY[3] = 0.; stress1_surfXY[4] = 0.; stress1_surfXY[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = duzdz_glb_surfXY(i, j, k) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k) * 2. + \
						duzdy_glb_surfXZ(i, j, k) + duzdy_glb_surfXZ(i, j + 1, k)) / 4. - Deyz0_glb(i, j, k)) * 2.;

					strain[4] = ((duxdz_glb_surfXY(i, j, k) * 2. + \
						duzdx_glb_surfYZ(i, j, k) + duzdx_glb_surfYZ(i + 1, j, k)) / 4. - Dexz0_glb(i, j, k)) * 2.;

					strain[5] = (Dexy_glb(i, j, k) - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress1_surfXY);

					if (k == 0) {
						mat_type = (pt_glb->material_cell(i, j, nz - 1));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, nz - 1) - Dexx0_glb(i, j, nz - 1);
						strain[1] = Deyy_glb(i, j, nz - 1) - Deyy0_glb(i, j, nz - 1);
						strain[2] = duzdz_glb_surfXY(i, j, k) - Dezz0_glb(i, j, nz - 1);

						strain[3] = ((duydz_glb_surfXY(i, j, k) * 2. + \
							duzdy_glb_surfXZ(i, j, nz - 1) + duzdy_glb_surfXZ(i, j + 1, nz - 1)) / 4. - Deyz0_glb(i, j, nz - 1)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k) * 2. + \
							duzdx_glb_surfYZ(i, j, nz - 1) + duzdx_glb_surfYZ(i + 1, j, nz - 1)) / 4. - Dexz0_glb(i, j, nz - 1)) * 2.;

						strain[5] = (Dexy_glb(i, j, nz - 1) - Dexy0_glb(i, j, nz - 1)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j, k - 1));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, k - 1) - Dexx0_glb(i, j, k - 1);
						strain[1] = Deyy_glb(i, j, k - 1) - Deyy0_glb(i, j, k - 1);
						strain[2] = duzdz_glb_surfXY(i, j, k) - Dezz0_glb(i, j, k - 1);

						strain[3] = ((duydz_glb_surfXY(i, j, k) * 2. + \
							duzdy_glb_surfXZ(i, j, k - 1) + duzdy_glb_surfXZ(i, j + 1, k - 1)) / 4. - Deyz0_glb(i, j, k - 1)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k) * 2. + \
							duzdx_glb_surfYZ(i, j, k - 1) + duzdx_glb_surfYZ(i + 1, j, k - 1)) / 4. - Dexz0_glb(i, j, k - 1)) * 2.;

						strain[5] = (Dexy_glb(i, j, k - 1) - Dexy0_glb(i, j, k - 1)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress1_surfXY[0] = (stress1_surfXY[0] + stress[0]) / 2.;
					stress1_surfXY[1] = (stress1_surfXY[1] + stress[1]) / 2.;
					stress1_surfXY[2] = (stress1_surfXY[2] + stress[2]) / 2.;
					stress1_surfXY[3] = (stress1_surfXY[3] + stress[3]) / 2.;
					stress1_surfXY[4] = (stress1_surfXY[4] + stress[4]) / 2.;
					stress1_surfXY[5] = (stress1_surfXY[5] + stress[5]) / 2.;
				}

				//----------XY surface in +Z direction-----//
				if (stress_free_surfXY(i, j, k + 1) == true) {
					if (pt_glb->if_gaussian_strain_pulse == true && pt_glb->input_strain_type == 3 && k == source_z) {

						temporal_var = pt_glb->time_device - 0.5 * pt_glb->dt - 5.0 * pt_glb->sigma_gauss;
						if (pt_glb->input_strain_component == 'x') {
							stress2_surfXY[0] = 0.;
							stress2_surfXY[1] = 0.;
							stress2_surfXY[2] = 0.;
							stress2_surfXY[3] = 0.;
							stress2_surfXY[4] = pt_glb->amplitude_gauss * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->sigma_gauss, 2.)));
							stress2_surfXY[5] = 0.;
						}
						else if (pt_glb->input_strain_component == 'y') {
							stress2_surfXY[0] = 0.;
							stress2_surfXY[1] = 0.;
							stress2_surfXY[2] = 0.;
							stress2_surfXY[3] = pt_glb->amplitude_gauss * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->sigma_gauss, 2.)));
							stress2_surfXY[4] = 0.;
							stress2_surfXY[5] = 0.;
						}
						else if (pt_glb->input_strain_component == 'z') {
							stress2_surfXY[0] = 0.;
							stress2_surfXY[1] = 0.;
							stress2_surfXY[2] = pt_glb->amplitude_gauss * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->sigma_gauss, 2.)));
							stress2_surfXY[3] = 0.;
							stress2_surfXY[4] = 0.;
							stress2_surfXY[5] = 0.;
						}

					}
					else {
						stress2_surfXY[0] = 0.; stress2_surfXY[1] = 0.; stress2_surfXY[2] = 0.;
						stress2_surfXY[3] = 0.; stress2_surfXY[4] = 0.; stress2_surfXY[5] = 0.;
					}
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = duzdz_glb_surfXY(i, j, k + 1) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k + 1) * 2. + \
						duzdy_glb_surfXZ(i, j, k) + duzdy_glb_surfXZ(i, j + 1, k)) / 4. - Deyz0_glb(i, j, k)) * 2.;

					strain[4] = ((duxdz_glb_surfXY(i, j, k + 1) * 2. + \
						duzdx_glb_surfYZ(i, j, k) + duzdx_glb_surfYZ(i + 1, j, k)) / 4. - Dexz0_glb(i, j, k)) * 2.;

					strain[5] = (Dexy_glb(i, j, k) - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress2_surfXY);

					if (k == ny - 1) {
						mat_type = (pt_glb->material_cell(i, j, 0));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, 0) - Dexx0_glb(i, j, 0);
						strain[1] = Deyy_glb(i, j, 0) - Deyy0_glb(i, j, 0);
						strain[2] = duzdz_glb_surfXY(i, j, k + 1) - Dezz0_glb(i, j, 0);

						strain[3] = ((duydz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdy_glb_surfXZ(i, j, 0) + duzdy_glb_surfXZ(i, j + 1, 0)) / 4. - Deyz0_glb(i, j, 0)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdx_glb_surfYZ(i, j, 0) + duzdx_glb_surfYZ(i + 1, j, 0)) / 4. - Dexz0_glb(i, j, 0)) * 2.;

						strain[5] = (Dexy_glb(i, j, 0) - Dexy0_glb(i, j, 0)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j, k + 1));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, k + 1) - Dexx0_glb(i, j, k + 1);
						strain[1] = Deyy_glb(i, j, k + 1) - Deyy0_glb(i, j, k + 1);
						strain[2] = duzdz_glb_surfXY(i, j, k + 1) - Dezz0_glb(i, j, k + 1);

						strain[3] = ((duydz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdy_glb_surfXZ(i, j, k + 1) + duzdy_glb_surfXZ(i, j + 1, k + 1)) / 4. - Deyz0_glb(i, j, k + 1)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdx_glb_surfYZ(i, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k + 1)) / 4. - Dexz0_glb(i, j, k + 1)) * 2.;

						strain[5] = (Dexy_glb(i, j, k + 1) - Dexy0_glb(i, j, k + 1)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress2_surfXY[0] = (stress2_surfXY[0] + stress[0]) / 2.;
					stress2_surfXY[1] = (stress2_surfXY[1] + stress[1]) / 2.;
					stress2_surfXY[2] = (stress2_surfXY[2] + stress[2]) / 2.;
					stress2_surfXY[3] = (stress2_surfXY[3] + stress[3]) / 2.;
					stress2_surfXY[4] = (stress2_surfXY[4] + stress[4]) / 2.;
					stress2_surfXY[5] = (stress2_surfXY[5] + stress[5]) / 2.;
				}
				//----------END XY surface in +-Z direction-----//

				mat_type = (pt_glb->material_cell(i, j, k));
				mat = &(pt_glb->material_parameters[(mat_type)-1]);

				density = pt_glb->dt / mat->density;
				mass_damping = -pt_glb->dt * mat->elast_mass_damping;
				stiff_damping = mat->elast_stiff_damping / mat->density;

				force_x(i, j, k) = (stress2_surfYZ[0] - stress1_surfYZ[0]) / dx + \
					(stress2_surfXZ[5] - stress1_surfXZ[5]) / dy + \
					(stress2_surfXY[4] - stress1_surfXY[4]) / dz;

				force_y(i, j, k) = (stress2_surfYZ[5] - stress1_surfYZ[5]) / dx + \
					(stress2_surfXZ[1] - stress1_surfXZ[1]) / dy + \
					(stress2_surfXY[3] - stress1_surfXY[3]) / dz;

				force_z(i, j, k) = (stress2_surfYZ[4] - stress1_surfYZ[4]) / dx + \
					(stress2_surfXZ[3] - stress1_surfXZ[3]) / dy + \
					(stress2_surfXY[2] - stress1_surfXY[2]) / dz;

				dvx_glb_rk2(i, j, k) = force_x(i, j, k) * density \
					+ mass_damping * vx_glb(i, j, k) \
					+ stiff_damping * (force_x(i, j, k) - force_x_store(i, j, k)) * 2.;

				dvy_glb_rk2(i, j, k) = force_y(i, j, k) * density \
					+ mass_damping * vy_glb(i, j, k) \
					+ stiff_damping * (force_y(i, j, k) - force_y_store(i, j, k)) * 2.;

				dvz_glb_rk2(i, j, k) = force_z(i, j, k) * density \
					+ mass_damping * vz_glb(i, j, k) \
					+ stiff_damping * (force_z(i, j, k) - force_z_store(i, j, k)) * 2.;
			}
		}
	}
}

void elastic_system::get_dv_RK3() {
	double stress1_surfYZ[6], stress2_surfYZ[6];
	double stress1_surfXZ[6], stress2_surfXZ[6];
	double stress1_surfXY[6], stress2_surfXY[6];
	double strain[6], stress[6];
	unsigned int mat_type;
	material* mat;
	double density, mass_damping, stiff_damping;
	long int i, j, k;
	double temporal_var;

	if (pt_glb->if_elasto_backaction_from_pORm == true) {
#pragma acc parallel default(present) async(1)
		{
			get_eigenstrain_dynamic();
		}
	}

#pragma acc wait(1) async(2)
#pragma acc parallel default(present) async(2)
	{
#pragma acc loop gang vector private(stress1_surfYZ,stress2_surfYZ,stress1_surfXZ,stress2_surfXZ,stress1_surfXY,stress2_surfXY,\
strain,stress,mat_type,density, mass_damping, stiff_damping, mat,i,j,k,temporal_var)
		for (long int id = 0; id < n; id++) {
			i = id / (ny * nz);
			j = (id - i * (ny * nz)) / nz;
			k = id - i * (ny * nz) - j * nz;

			mat_type = (pt_glb->material_cell(id));
			if (mat_type == 0) {
				dvx_glb_rk3(id) = 0.;
				dvy_glb_rk3(id) = 0.;
				dvz_glb_rk3(id) = 0.;
			}
			else {
				//----------YZ surface in -X direction-----//
				if (stress_free_surfYZ(i, j, k) == true) {
					stress1_surfYZ[0] = 0.; stress1_surfYZ[1] = 0.; stress1_surfYZ[2] = 0.;
					stress1_surfYZ[3] = 0.; stress1_surfYZ[4] = 0.; stress1_surfYZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = duxdx_glb_surfYZ(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);
					strain[3] = (Deyz_glb(i, j, k) - Deyz0_glb(i, j, k)) * 2.;
					strain[4] = ((duxdz_glb_surfXY(i, j, k) + duxdz_glb_surfXY(i, j, k + 1) + duzdx_glb_surfYZ(i, j, k) * 2.) / 4. \
						- Dexz0_glb(i, j, k)) * 2.;
					strain[5] = ((duxdy_glb_surfXZ(i, j, k) + duxdy_glb_surfXZ(i, j + 1, k) + duydx_glb_surfYZ(i, j, k) * 2.) / 4. \
						- Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress1_surfYZ);

					if (i == 0) {
						mat_type = (pt_glb->material_cell(nx - 1, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i, j, k) - Dexx0_glb(nx - 1, j, k);
						strain[1] = Deyy_glb(nx - 1, j, k) - Deyy0_glb(nx - 1, j, k);
						strain[2] = Dezz_glb(nx - 1, j, k) - Dezz0_glb(nx - 1, j, k);
						strain[3] = (Deyz_glb(nx - 1, j, k) - Deyz0_glb(nx - 1, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(nx - 1, j, k) + duxdz_glb_surfXY(nx - 1, j, k + 1) + duzdx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexz0_glb(nx - 1, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(nx - 1, j, k) + duxdy_glb_surfXZ(nx - 1, j + 1, k) + duydx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexy0_glb(nx - 1, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i - 1, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i, j, k) - Dexx0_glb(i - 1, j, k);
						strain[1] = Deyy_glb(i - 1, j, k) - Deyy0_glb(i - 1, j, k);
						strain[2] = Dezz_glb(i - 1, j, k) - Dezz0_glb(i - 1, j, k);
						strain[3] = (Deyz_glb(i - 1, j, k) - Deyz0_glb(i - 1, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(i - 1, j, k) + duxdz_glb_surfXY(i - 1, j, k + 1) + duzdx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexz0_glb(i - 1, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(i - 1, j, k) + duxdy_glb_surfXZ(i - 1, j + 1, k) + duydx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexy0_glb(i - 1, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress1_surfYZ[0] = (stress1_surfYZ[0] + stress[0]) / 2.;
					stress1_surfYZ[1] = (stress1_surfYZ[1] + stress[1]) / 2.;
					stress1_surfYZ[2] = (stress1_surfYZ[2] + stress[2]) / 2.;
					stress1_surfYZ[3] = (stress1_surfYZ[3] + stress[3]) / 2.;
					stress1_surfYZ[4] = (stress1_surfYZ[4] + stress[4]) / 2.;
					stress1_surfYZ[5] = (stress1_surfYZ[5] + stress[5]) / 2.;
				}

				//----------YZ surface in +X direction-----//
				if (stress_free_surfYZ(i + 1, j, k) == true) {
					stress2_surfYZ[0] = 0.; stress2_surfYZ[1] = 0.; stress2_surfYZ[2] = 0.;
					stress2_surfYZ[3] = 0.; stress2_surfYZ[4] = 0.; stress2_surfYZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = duxdx_glb_surfYZ(i + 1, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);
					strain[3] = (Deyz_glb(i, j, k) - Deyz0_glb(i, j, k)) * 2.;
					strain[4] = ((duxdz_glb_surfXY(i, j, k) + duxdz_glb_surfXY(i, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
						- Dexz0_glb(i, j, k)) * 2.;
					strain[5] = ((duxdy_glb_surfXZ(i, j, k) + duxdy_glb_surfXZ(i, j + 1, k) + duydx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
						- Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress2_surfYZ);

					if (i == nx - 1) {
						mat_type = (pt_glb->material_cell(0, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i + 1, j, k) - Dexx0_glb(0, j, k);
						strain[1] = Deyy_glb(0, j, k) - Deyy0_glb(0, j, k);
						strain[2] = Dezz_glb(0, j, k) - Dezz0_glb(0, j, k);
						strain[3] = (Deyz_glb(0, j, k) - Deyz0_glb(0, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(0, j, k) + duxdz_glb_surfXY(0, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexz0_glb(0, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(0, j, k) + duxdy_glb_surfXZ(0, j + 1, k) + duydx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexy0_glb(0, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i + 1, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i + 1, j, k) - Dexx0_glb(i + 1, j, k);
						strain[1] = Deyy_glb(i + 1, j, k) - Deyy0_glb(i + 1, j, k);
						strain[2] = Dezz_glb(i + 1, j, k) - Dezz0_glb(i + 1, j, k);
						strain[3] = (Deyz_glb(i + 1, j, k) - Deyz0_glb(i + 1, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(i + 1, j, k) + duxdz_glb_surfXY(i + 1, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexz0_glb(i + 1, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(i + 1, j, k) + duxdy_glb_surfXZ(i + 1, j + 1, k) + duydx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexy0_glb(i + 1, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress2_surfYZ[0] = (stress2_surfYZ[0] + stress[0]) / 2.;
					stress2_surfYZ[1] = (stress2_surfYZ[1] + stress[1]) / 2.;
					stress2_surfYZ[2] = (stress2_surfYZ[2] + stress[2]) / 2.;
					stress2_surfYZ[3] = (stress2_surfYZ[3] + stress[3]) / 2.;
					stress2_surfYZ[4] = (stress2_surfYZ[4] + stress[4]) / 2.;
					stress2_surfYZ[5] = (stress2_surfYZ[5] + stress[5]) / 2.;
				}
				//----------END YZ surface in +-X direction-----//

				//----------XZ surface in -Y direction-----//
				if (stress_free_surfXZ(i, j, k) == true) {
					stress1_surfXZ[0] = 0.; stress1_surfXZ[1] = 0.; stress1_surfXZ[2] = 0.;
					stress1_surfXZ[3] = 0.; stress1_surfXZ[4] = 0.; stress1_surfXZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = duydy_glb_surfXZ(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k) + duydz_glb_surfXY(i, j, k + 1) + \
						duzdy_glb_surfXZ(i, j, k) * 2.) / 4. - Deyz0_glb(i, j, k)) * 2.;
					strain[4] = (Dexz_glb(i, j, k) - Dexz0_glb(i, j, k)) * 2.;
					strain[5] = ((duxdy_glb_surfXZ(i, j, k) * 2. + \
						duydx_glb_surfYZ(i, j, k) + duydx_glb_surfYZ(i + 1, j, k)) / 4. - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress1_surfXZ);

					if (j == 0) {
						mat_type = (pt_glb->material_cell(i, ny - 1, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, ny - 1, k) - Dexx0_glb(i, ny - 1, k);
						strain[1] = duydy_glb_surfXZ(i, j, k) - Deyy0_glb(i, ny - 1, k);
						strain[2] = Dezz_glb(i, ny - 1, k) - Dezz0_glb(i, ny - 1, k);

						strain[3] = ((duydz_glb_surfXY(i, ny - 1, k) + duydz_glb_surfXY(i, ny - 1, k + 1) + \
							duzdy_glb_surfXZ(i, j, k) * 2.) / 4. - Deyz0_glb(i, ny - 1, k)) * 2.;

						strain[4] = (Dexz_glb(i, ny - 1, k) - Dexz0_glb(i, ny - 1, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j, k) * 2. + \
							duydx_glb_surfYZ(i, ny - 1, k) + duydx_glb_surfYZ(i + 1, ny - 1, k)) / 4. - Dexy0_glb(i, ny - 1, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j - 1, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j - 1, k) - Dexx0_glb(i, j - 1, k);
						strain[1] = duydy_glb_surfXZ(i, j, k) - Deyy0_glb(i, j - 1, k);
						strain[2] = Dezz_glb(i, j - 1, k) - Dezz0_glb(i, j - 1, k);

						strain[3] = ((duydz_glb_surfXY(i, j - 1, k) + duydz_glb_surfXY(i, j - 1, k + 1) + \
							duzdy_glb_surfXZ(i, j, k) * 2.) / 4. - Deyz0_glb(i, j - 1, k)) * 2.;

						strain[4] = (Dexz_glb(i, j - 1, k) - Dexz0_glb(i, j - 1, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j, k) * 2. + \
							duydx_glb_surfYZ(i, j - 1, k) + duydx_glb_surfYZ(i + 1, j - 1, k)) / 4. - Dexy0_glb(i, j - 1, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress1_surfXZ[0] = (stress1_surfXZ[0] + stress[0]) / 2.;
					stress1_surfXZ[1] = (stress1_surfXZ[1] + stress[1]) / 2.;
					stress1_surfXZ[2] = (stress1_surfXZ[2] + stress[2]) / 2.;
					stress1_surfXZ[3] = (stress1_surfXZ[3] + stress[3]) / 2.;
					stress1_surfXZ[4] = (stress1_surfXZ[4] + stress[4]) / 2.;
					stress1_surfXZ[5] = (stress1_surfXZ[5] + stress[5]) / 2.;
				}

				//----------XZ surface in +Y direction-----//
				if (stress_free_surfXZ(i, j + 1, k) == true) {
					stress2_surfXZ[0] = 0.; stress2_surfXZ[1] = 0.; stress2_surfXZ[2] = 0.;
					stress2_surfXZ[3] = 0.; stress2_surfXZ[4] = 0.; stress2_surfXZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = duydy_glb_surfXZ(i, j + 1, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k) + duydz_glb_surfXY(i, j, k + 1) + \
						duzdy_glb_surfXZ(i, j + 1, k) * 2.) / 4. - Deyz0_glb(i, j, k)) * 2.;

					strain[4] = (Dexz_glb(i, j, k) - Dexz0_glb(i, j, k)) * 2.;

					strain[5] = ((duxdy_glb_surfXZ(i, j + 1, k) * 2. + \
						duydx_glb_surfYZ(i, j, k) + duydx_glb_surfYZ(i + 1, j, k)) / 4. - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress2_surfXZ);

					if (j == ny - 1) {
						mat_type = (pt_glb->material_cell(i, 0, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, 0, k) - Dexx0_glb(i, 0, k);
						strain[1] = duydy_glb_surfXZ(i, j + 1, k) - Deyy0_glb(i, 0, k);
						strain[2] = Dezz_glb(i, 0, k) - Dezz0_glb(i, 0, k);

						strain[3] = ((duydz_glb_surfXY(i, 0, k) + duydz_glb_surfXY(i, 0, k + 1) + \
							duzdy_glb_surfXZ(i, j + 1, k) * 2.) / 4. - Deyz0_glb(i, 0, k)) * 2.;

						strain[4] = (Dexz_glb(i, 0, k) - Dexz0_glb(i, 0, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j + 1, k) * 2. + \
							duydx_glb_surfYZ(i, 0, k) + duydx_glb_surfYZ(i + 1, 0, k)) / 4. - Dexy0_glb(i, 0, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j + 1, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j + 1, k) - Dexx0_glb(i, j + 1, k);
						strain[1] = duydy_glb_surfXZ(i, j + 1, k) - Deyy0_glb(i, j + 1, k);
						strain[2] = Dezz_glb(i, j + 1, k) - Dezz0_glb(i, j + 1, k);

						strain[3] = ((duydz_glb_surfXY(i, j + 1, k) + duydz_glb_surfXY(i, j + 1, k + 1) + \
							duzdy_glb_surfXZ(i, j + 1, k) * 2.) / 4. - Deyz0_glb(i, j + 1, k)) * 2.;

						strain[4] = (Dexz_glb(i, j + 1, k) - Dexz0_glb(i, j + 1, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j + 1, k) * 2. + \
							duydx_glb_surfYZ(i, j + 1, k) + duydx_glb_surfYZ(i + 1, j + 1, k)) / 4. - Dexy0_glb(i, j + 1, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress2_surfXZ[0] = (stress2_surfXZ[0] + stress[0]) / 2.;
					stress2_surfXZ[1] = (stress2_surfXZ[1] + stress[1]) / 2.;
					stress2_surfXZ[2] = (stress2_surfXZ[2] + stress[2]) / 2.;
					stress2_surfXZ[3] = (stress2_surfXZ[3] + stress[3]) / 2.;
					stress2_surfXZ[4] = (stress2_surfXZ[4] + stress[4]) / 2.;
					stress2_surfXZ[5] = (stress2_surfXZ[5] + stress[5]) / 2.;
				}
				//----------END XZ surface in +-Y direction-----//

				//----------XY surface in -Z direction-----//
				if (stress_free_surfXY(i, j, k) == true) {
					stress1_surfXY[0] = 0.; stress1_surfXY[1] = 0.; stress1_surfXY[2] = 0.;
					stress1_surfXY[3] = 0.; stress1_surfXY[4] = 0.; stress1_surfXY[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = duzdz_glb_surfXY(i, j, k) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k) * 2. + \
						duzdy_glb_surfXZ(i, j, k) + duzdy_glb_surfXZ(i, j + 1, k)) / 4. - Deyz0_glb(i, j, k)) * 2.;

					strain[4] = ((duxdz_glb_surfXY(i, j, k) * 2. + \
						duzdx_glb_surfYZ(i, j, k) + duzdx_glb_surfYZ(i + 1, j, k)) / 4. - Dexz0_glb(i, j, k)) * 2.;

					strain[5] = (Dexy_glb(i, j, k) - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress1_surfXY);

					if (k == 0) {
						mat_type = (pt_glb->material_cell(i, j, nz - 1));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, nz - 1) - Dexx0_glb(i, j, nz - 1);
						strain[1] = Deyy_glb(i, j, nz - 1) - Deyy0_glb(i, j, nz - 1);
						strain[2] = duzdz_glb_surfXY(i, j, k) - Dezz0_glb(i, j, nz - 1);

						strain[3] = ((duydz_glb_surfXY(i, j, k) * 2. + \
							duzdy_glb_surfXZ(i, j, nz - 1) + duzdy_glb_surfXZ(i, j + 1, nz - 1)) / 4. - Deyz0_glb(i, j, nz - 1)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k) * 2. + \
							duzdx_glb_surfYZ(i, j, nz - 1) + duzdx_glb_surfYZ(i + 1, j, nz - 1)) / 4. - Dexz0_glb(i, j, nz - 1)) * 2.;

						strain[5] = (Dexy_glb(i, j, nz - 1) - Dexy0_glb(i, j, nz - 1)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j, k - 1));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, k - 1) - Dexx0_glb(i, j, k - 1);
						strain[1] = Deyy_glb(i, j, k - 1) - Deyy0_glb(i, j, k - 1);
						strain[2] = duzdz_glb_surfXY(i, j, k) - Dezz0_glb(i, j, k - 1);

						strain[3] = ((duydz_glb_surfXY(i, j, k) * 2. + \
							duzdy_glb_surfXZ(i, j, k - 1) + duzdy_glb_surfXZ(i, j + 1, k - 1)) / 4. - Deyz0_glb(i, j, k - 1)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k) * 2. + \
							duzdx_glb_surfYZ(i, j, k - 1) + duzdx_glb_surfYZ(i + 1, j, k - 1)) / 4. - Dexz0_glb(i, j, k - 1)) * 2.;

						strain[5] = (Dexy_glb(i, j, k - 1) - Dexy0_glb(i, j, k - 1)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress1_surfXY[0] = (stress1_surfXY[0] + stress[0]) / 2.;
					stress1_surfXY[1] = (stress1_surfXY[1] + stress[1]) / 2.;
					stress1_surfXY[2] = (stress1_surfXY[2] + stress[2]) / 2.;
					stress1_surfXY[3] = (stress1_surfXY[3] + stress[3]) / 2.;
					stress1_surfXY[4] = (stress1_surfXY[4] + stress[4]) / 2.;
					stress1_surfXY[5] = (stress1_surfXY[5] + stress[5]) / 2.;
				}

				//----------XY surface in +Z direction-----//
				if (stress_free_surfXY(i, j, k + 1) == true) {
					if (pt_glb->if_gaussian_strain_pulse == true && pt_glb->input_strain_type == 3 && k == source_z) {

						temporal_var = pt_glb->time_device - 0.5 * pt_glb->dt - 5.0 * pt_glb->sigma_gauss;
						if (pt_glb->input_strain_component == 'x') {
							stress2_surfXY[0] = 0.;
							stress2_surfXY[1] = 0.;
							stress2_surfXY[2] = 0.;
							stress2_surfXY[3] = 0.;
							stress2_surfXY[4] = pt_glb->amplitude_gauss * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->sigma_gauss, 2.)));
							stress2_surfXY[5] = 0.;
						}
						else if (pt_glb->input_strain_component == 'y') {
							stress2_surfXY[0] = 0.;
							stress2_surfXY[1] = 0.;
							stress2_surfXY[2] = 0.;
							stress2_surfXY[3] = pt_glb->amplitude_gauss * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->sigma_gauss, 2.)));
							stress2_surfXY[4] = 0.;
							stress2_surfXY[5] = 0.;
						}
						else if (pt_glb->input_strain_component == 'z') {
							stress2_surfXY[0] = 0.;
							stress2_surfXY[1] = 0.;
							stress2_surfXY[2] = pt_glb->amplitude_gauss * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->sigma_gauss, 2.)));
							stress2_surfXY[3] = 0.;
							stress2_surfXY[4] = 0.;
							stress2_surfXY[5] = 0.;
						}

					}
					else {
						stress2_surfXY[0] = 0.; stress2_surfXY[1] = 0.; stress2_surfXY[2] = 0.;
						stress2_surfXY[3] = 0.; stress2_surfXY[4] = 0.; stress2_surfXY[5] = 0.;
					}
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = duzdz_glb_surfXY(i, j, k + 1) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k + 1) * 2. + \
						duzdy_glb_surfXZ(i, j, k) + duzdy_glb_surfXZ(i, j + 1, k)) / 4. - Deyz0_glb(i, j, k)) * 2.;

					strain[4] = ((duxdz_glb_surfXY(i, j, k + 1) * 2. + \
						duzdx_glb_surfYZ(i, j, k) + duzdx_glb_surfYZ(i + 1, j, k)) / 4. - Dexz0_glb(i, j, k)) * 2.;

					strain[5] = (Dexy_glb(i, j, k) - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress2_surfXY);

					if (k == ny - 1) {
						mat_type = (pt_glb->material_cell(i, j, 0));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, 0) - Dexx0_glb(i, j, 0);
						strain[1] = Deyy_glb(i, j, 0) - Deyy0_glb(i, j, 0);
						strain[2] = duzdz_glb_surfXY(i, j, k + 1) - Dezz0_glb(i, j, 0);

						strain[3] = ((duydz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdy_glb_surfXZ(i, j, 0) + duzdy_glb_surfXZ(i, j + 1, 0)) / 4. - Deyz0_glb(i, j, 0)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdx_glb_surfYZ(i, j, 0) + duzdx_glb_surfYZ(i + 1, j, 0)) / 4. - Dexz0_glb(i, j, 0)) * 2.;

						strain[5] = (Dexy_glb(i, j, 0) - Dexy0_glb(i, j, 0)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j, k + 1));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, k + 1) - Dexx0_glb(i, j, k + 1);
						strain[1] = Deyy_glb(i, j, k + 1) - Deyy0_glb(i, j, k + 1);
						strain[2] = duzdz_glb_surfXY(i, j, k + 1) - Dezz0_glb(i, j, k + 1);

						strain[3] = ((duydz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdy_glb_surfXZ(i, j, k + 1) + duzdy_glb_surfXZ(i, j + 1, k + 1)) / 4. - Deyz0_glb(i, j, k + 1)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdx_glb_surfYZ(i, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k + 1)) / 4. - Dexz0_glb(i, j, k + 1)) * 2.;

						strain[5] = (Dexy_glb(i, j, k + 1) - Dexy0_glb(i, j, k + 1)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress2_surfXY[0] = (stress2_surfXY[0] + stress[0]) / 2.;
					stress2_surfXY[1] = (stress2_surfXY[1] + stress[1]) / 2.;
					stress2_surfXY[2] = (stress2_surfXY[2] + stress[2]) / 2.;
					stress2_surfXY[3] = (stress2_surfXY[3] + stress[3]) / 2.;
					stress2_surfXY[4] = (stress2_surfXY[4] + stress[4]) / 2.;
					stress2_surfXY[5] = (stress2_surfXY[5] + stress[5]) / 2.;
				}
				//----------END XY surface in +-Z direction-----//

				mat_type = (pt_glb->material_cell(i, j, k));
				mat = &(pt_glb->material_parameters[(mat_type)-1]);

				density = pt_glb->dt / mat->density;
				mass_damping = -pt_glb->dt * mat->elast_mass_damping;
				stiff_damping = mat->elast_stiff_damping / mat->density;

				force_x(i, j, k) = (stress2_surfYZ[0] - stress1_surfYZ[0]) / dx + \
					(stress2_surfXZ[5] - stress1_surfXZ[5]) / dy + \
					(stress2_surfXY[4] - stress1_surfXY[4]) / dz;

				force_y(i, j, k) = (stress2_surfYZ[5] - stress1_surfYZ[5]) / dx + \
					(stress2_surfXZ[1] - stress1_surfXZ[1]) / dy + \
					(stress2_surfXY[3] - stress1_surfXY[3]) / dz;

				force_z(i, j, k) = (stress2_surfYZ[4] - stress1_surfYZ[4]) / dx + \
					(stress2_surfXZ[3] - stress1_surfXZ[3]) / dy + \
					(stress2_surfXY[2] - stress1_surfXY[2]) / dz;

				dvx_glb_rk3(i, j, k) = force_x(i, j, k) * density \
					+ mass_damping * vx_glb(i, j, k) \
					+ stiff_damping * (force_x(i, j, k) - force_x_store(i, j, k)) * 2.;

				dvy_glb_rk3(i, j, k) = force_y(i, j, k) * density \
					+ mass_damping * vy_glb(i, j, k) \
					+ stiff_damping * (force_y(i, j, k) - force_y_store(i, j, k)) * 2.;

				dvz_glb_rk3(i, j, k) = force_z(i, j, k) * density \
					+ mass_damping * vz_glb(i, j, k) \
					+ stiff_damping * (force_z(i, j, k) - force_z_store(i, j, k)) * 2.;
			}
		}
	}
}

void elastic_system::get_dv_RK4() {
	double stress1_surfYZ[6], stress2_surfYZ[6];
	double stress1_surfXZ[6], stress2_surfXZ[6];
	double stress1_surfXY[6], stress2_surfXY[6];
	double strain[6], stress[6];
	unsigned int mat_type;
	material* mat;
	double density, mass_damping, stiff_damping;
	long int i, j, k;
	double temporal_var;

	if (pt_glb->if_elasto_backaction_from_pORm == true) {
#pragma acc parallel default(present) async(1)
		{
			get_eigenstrain_dynamic();
		}
	}

#pragma acc wait(1) async(2)
#pragma acc parallel default(present) async(2)
	{
#pragma acc loop gang vector private(stress1_surfYZ,stress2_surfYZ,stress1_surfXZ,stress2_surfXZ,stress1_surfXY,stress2_surfXY,\
strain,stress,mat_type,density, mass_damping, stiff_damping, mat,i,j,k,temporal_var)
		for (long int id = 0; id < n; id++) {
			i = id / (ny * nz);
			j = (id - i * (ny * nz)) / nz;
			k = id - i * (ny * nz) - j * nz;

			mat_type = (pt_glb->material_cell(id));
			if (mat_type == 0) {
				dvx_glb_rk4(id) = 0.;
				dvy_glb_rk4(id) = 0.;
				dvz_glb_rk4(id) = 0.;
			}
			else {
				//----------YZ surface in -X direction-----//
				if (stress_free_surfYZ(i, j, k) == true) {
					stress1_surfYZ[0] = 0.; stress1_surfYZ[1] = 0.; stress1_surfYZ[2] = 0.;
					stress1_surfYZ[3] = 0.; stress1_surfYZ[4] = 0.; stress1_surfYZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = duxdx_glb_surfYZ(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);
					strain[3] = (Deyz_glb(i, j, k) - Deyz0_glb(i, j, k)) * 2.;
					strain[4] = ((duxdz_glb_surfXY(i, j, k) + duxdz_glb_surfXY(i, j, k + 1) + duzdx_glb_surfYZ(i, j, k) * 2.) / 4. \
						- Dexz0_glb(i, j, k)) * 2.;
					strain[5] = ((duxdy_glb_surfXZ(i, j, k) + duxdy_glb_surfXZ(i, j + 1, k) + duydx_glb_surfYZ(i, j, k) * 2.) / 4. \
						- Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress1_surfYZ);

					if (i == 0) {
						mat_type = (pt_glb->material_cell(nx - 1, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i, j, k) - Dexx0_glb(nx - 1, j, k);
						strain[1] = Deyy_glb(nx - 1, j, k) - Deyy0_glb(nx - 1, j, k);
						strain[2] = Dezz_glb(nx - 1, j, k) - Dezz0_glb(nx - 1, j, k);
						strain[3] = (Deyz_glb(nx - 1, j, k) - Deyz0_glb(nx - 1, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(nx - 1, j, k) + duxdz_glb_surfXY(nx - 1, j, k + 1) + duzdx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexz0_glb(nx - 1, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(nx - 1, j, k) + duxdy_glb_surfXZ(nx - 1, j + 1, k) + duydx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexy0_glb(nx - 1, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i - 1, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i, j, k) - Dexx0_glb(i - 1, j, k);
						strain[1] = Deyy_glb(i - 1, j, k) - Deyy0_glb(i - 1, j, k);
						strain[2] = Dezz_glb(i - 1, j, k) - Dezz0_glb(i - 1, j, k);
						strain[3] = (Deyz_glb(i - 1, j, k) - Deyz0_glb(i - 1, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(i - 1, j, k) + duxdz_glb_surfXY(i - 1, j, k + 1) + duzdx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexz0_glb(i - 1, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(i - 1, j, k) + duxdy_glb_surfXZ(i - 1, j + 1, k) + duydx_glb_surfYZ(i, j, k) * 2.) / 4. \
							- Dexy0_glb(i - 1, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress1_surfYZ[0] = (stress1_surfYZ[0] + stress[0]) / 2.;
					stress1_surfYZ[1] = (stress1_surfYZ[1] + stress[1]) / 2.;
					stress1_surfYZ[2] = (stress1_surfYZ[2] + stress[2]) / 2.;
					stress1_surfYZ[3] = (stress1_surfYZ[3] + stress[3]) / 2.;
					stress1_surfYZ[4] = (stress1_surfYZ[4] + stress[4]) / 2.;
					stress1_surfYZ[5] = (stress1_surfYZ[5] + stress[5]) / 2.;
				}

				//----------YZ surface in +X direction-----//
				if (stress_free_surfYZ(i + 1, j, k) == true) {
					stress2_surfYZ[0] = 0.; stress2_surfYZ[1] = 0.; stress2_surfYZ[2] = 0.;
					stress2_surfYZ[3] = 0.; stress2_surfYZ[4] = 0.; stress2_surfYZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = duxdx_glb_surfYZ(i + 1, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);
					strain[3] = (Deyz_glb(i, j, k) - Deyz0_glb(i, j, k)) * 2.;
					strain[4] = ((duxdz_glb_surfXY(i, j, k) + duxdz_glb_surfXY(i, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
						- Dexz0_glb(i, j, k)) * 2.;
					strain[5] = ((duxdy_glb_surfXZ(i, j, k) + duxdy_glb_surfXZ(i, j + 1, k) + duydx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
						- Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress2_surfYZ);

					if (i == nx - 1) {
						mat_type = (pt_glb->material_cell(0, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i + 1, j, k) - Dexx0_glb(0, j, k);
						strain[1] = Deyy_glb(0, j, k) - Deyy0_glb(0, j, k);
						strain[2] = Dezz_glb(0, j, k) - Dezz0_glb(0, j, k);
						strain[3] = (Deyz_glb(0, j, k) - Deyz0_glb(0, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(0, j, k) + duxdz_glb_surfXY(0, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexz0_glb(0, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(0, j, k) + duxdy_glb_surfXZ(0, j + 1, k) + duydx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexy0_glb(0, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i + 1, j, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = duxdx_glb_surfYZ(i + 1, j, k) - Dexx0_glb(i + 1, j, k);
						strain[1] = Deyy_glb(i + 1, j, k) - Deyy0_glb(i + 1, j, k);
						strain[2] = Dezz_glb(i + 1, j, k) - Dezz0_glb(i + 1, j, k);
						strain[3] = (Deyz_glb(i + 1, j, k) - Deyz0_glb(i + 1, j, k)) * 2.;
						strain[4] = ((duxdz_glb_surfXY(i + 1, j, k) + duxdz_glb_surfXY(i + 1, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexz0_glb(i + 1, j, k)) * 2.;
						strain[5] = ((duxdy_glb_surfXZ(i + 1, j, k) + duxdy_glb_surfXZ(i + 1, j + 1, k) + duydx_glb_surfYZ(i + 1, j, k) * 2.) / 4. \
							- Dexy0_glb(i + 1, j, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress2_surfYZ[0] = (stress2_surfYZ[0] + stress[0]) / 2.;
					stress2_surfYZ[1] = (stress2_surfYZ[1] + stress[1]) / 2.;
					stress2_surfYZ[2] = (stress2_surfYZ[2] + stress[2]) / 2.;
					stress2_surfYZ[3] = (stress2_surfYZ[3] + stress[3]) / 2.;
					stress2_surfYZ[4] = (stress2_surfYZ[4] + stress[4]) / 2.;
					stress2_surfYZ[5] = (stress2_surfYZ[5] + stress[5]) / 2.;
				}
				//----------END YZ surface in +-X direction-----//

				//----------XZ surface in -Y direction-----//
				if (stress_free_surfXZ(i, j, k) == true) {
					stress1_surfXZ[0] = 0.; stress1_surfXZ[1] = 0.; stress1_surfXZ[2] = 0.;
					stress1_surfXZ[3] = 0.; stress1_surfXZ[4] = 0.; stress1_surfXZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = duydy_glb_surfXZ(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k) + duydz_glb_surfXY(i, j, k + 1) + \
						duzdy_glb_surfXZ(i, j, k) * 2.) / 4. - Deyz0_glb(i, j, k)) * 2.;
					strain[4] = (Dexz_glb(i, j, k) - Dexz0_glb(i, j, k)) * 2.;
					strain[5] = ((duxdy_glb_surfXZ(i, j, k) * 2. + \
						duydx_glb_surfYZ(i, j, k) + duydx_glb_surfYZ(i + 1, j, k)) / 4. - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress1_surfXZ);

					if (j == 0) {
						mat_type = (pt_glb->material_cell(i, ny - 1, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, ny - 1, k) - Dexx0_glb(i, ny - 1, k);
						strain[1] = duydy_glb_surfXZ(i, j, k) - Deyy0_glb(i, ny - 1, k);
						strain[2] = Dezz_glb(i, ny - 1, k) - Dezz0_glb(i, ny - 1, k);

						strain[3] = ((duydz_glb_surfXY(i, ny - 1, k) + duydz_glb_surfXY(i, ny - 1, k + 1) + \
							duzdy_glb_surfXZ(i, j, k) * 2.) / 4. - Deyz0_glb(i, ny - 1, k)) * 2.;

						strain[4] = (Dexz_glb(i, ny - 1, k) - Dexz0_glb(i, ny - 1, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j, k) * 2. + \
							duydx_glb_surfYZ(i, ny - 1, k) + duydx_glb_surfYZ(i + 1, ny - 1, k)) / 4. - Dexy0_glb(i, ny - 1, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j - 1, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j - 1, k) - Dexx0_glb(i, j - 1, k);
						strain[1] = duydy_glb_surfXZ(i, j, k) - Deyy0_glb(i, j - 1, k);
						strain[2] = Dezz_glb(i, j - 1, k) - Dezz0_glb(i, j - 1, k);

						strain[3] = ((duydz_glb_surfXY(i, j - 1, k) + duydz_glb_surfXY(i, j - 1, k + 1) + \
							duzdy_glb_surfXZ(i, j, k) * 2.) / 4. - Deyz0_glb(i, j - 1, k)) * 2.;

						strain[4] = (Dexz_glb(i, j - 1, k) - Dexz0_glb(i, j - 1, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j, k) * 2. + \
							duydx_glb_surfYZ(i, j - 1, k) + duydx_glb_surfYZ(i + 1, j - 1, k)) / 4. - Dexy0_glb(i, j - 1, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress1_surfXZ[0] = (stress1_surfXZ[0] + stress[0]) / 2.;
					stress1_surfXZ[1] = (stress1_surfXZ[1] + stress[1]) / 2.;
					stress1_surfXZ[2] = (stress1_surfXZ[2] + stress[2]) / 2.;
					stress1_surfXZ[3] = (stress1_surfXZ[3] + stress[3]) / 2.;
					stress1_surfXZ[4] = (stress1_surfXZ[4] + stress[4]) / 2.;
					stress1_surfXZ[5] = (stress1_surfXZ[5] + stress[5]) / 2.;
				}

				//----------XZ surface in +Y direction-----//
				if (stress_free_surfXZ(i, j + 1, k) == true) {
					stress2_surfXZ[0] = 0.; stress2_surfXZ[1] = 0.; stress2_surfXZ[2] = 0.;
					stress2_surfXZ[3] = 0.; stress2_surfXZ[4] = 0.; stress2_surfXZ[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = duydy_glb_surfXZ(i, j + 1, k) - Deyy0_glb(i, j, k);
					strain[2] = Dezz_glb(i, j, k) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k) + duydz_glb_surfXY(i, j, k + 1) + \
						duzdy_glb_surfXZ(i, j + 1, k) * 2.) / 4. - Deyz0_glb(i, j, k)) * 2.;

					strain[4] = (Dexz_glb(i, j, k) - Dexz0_glb(i, j, k)) * 2.;

					strain[5] = ((duxdy_glb_surfXZ(i, j + 1, k) * 2. + \
						duydx_glb_surfYZ(i, j, k) + duydx_glb_surfYZ(i + 1, j, k)) / 4. - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress2_surfXZ);

					if (j == ny - 1) {
						mat_type = (pt_glb->material_cell(i, 0, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, 0, k) - Dexx0_glb(i, 0, k);
						strain[1] = duydy_glb_surfXZ(i, j + 1, k) - Deyy0_glb(i, 0, k);
						strain[2] = Dezz_glb(i, 0, k) - Dezz0_glb(i, 0, k);

						strain[3] = ((duydz_glb_surfXY(i, 0, k) + duydz_glb_surfXY(i, 0, k + 1) + \
							duzdy_glb_surfXZ(i, j + 1, k) * 2.) / 4. - Deyz0_glb(i, 0, k)) * 2.;

						strain[4] = (Dexz_glb(i, 0, k) - Dexz0_glb(i, 0, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j + 1, k) * 2. + \
							duydx_glb_surfYZ(i, 0, k) + duydx_glb_surfYZ(i + 1, 0, k)) / 4. - Dexy0_glb(i, 0, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j + 1, k));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j + 1, k) - Dexx0_glb(i, j + 1, k);
						strain[1] = duydy_glb_surfXZ(i, j + 1, k) - Deyy0_glb(i, j + 1, k);
						strain[2] = Dezz_glb(i, j + 1, k) - Dezz0_glb(i, j + 1, k);

						strain[3] = ((duydz_glb_surfXY(i, j + 1, k) + duydz_glb_surfXY(i, j + 1, k + 1) + \
							duzdy_glb_surfXZ(i, j + 1, k) * 2.) / 4. - Deyz0_glb(i, j + 1, k)) * 2.;

						strain[4] = (Dexz_glb(i, j + 1, k) - Dexz0_glb(i, j + 1, k)) * 2.;

						strain[5] = ((duxdy_glb_surfXZ(i, j + 1, k) * 2. + \
							duydx_glb_surfYZ(i, j + 1, k) + duydx_glb_surfYZ(i + 1, j + 1, k)) / 4. - Dexy0_glb(i, j + 1, k)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress2_surfXZ[0] = (stress2_surfXZ[0] + stress[0]) / 2.;
					stress2_surfXZ[1] = (stress2_surfXZ[1] + stress[1]) / 2.;
					stress2_surfXZ[2] = (stress2_surfXZ[2] + stress[2]) / 2.;
					stress2_surfXZ[3] = (stress2_surfXZ[3] + stress[3]) / 2.;
					stress2_surfXZ[4] = (stress2_surfXZ[4] + stress[4]) / 2.;
					stress2_surfXZ[5] = (stress2_surfXZ[5] + stress[5]) / 2.;
				}
				//----------END XZ surface in +-Y direction-----//

				//----------XY surface in -Z direction-----//
				if (stress_free_surfXY(i, j, k) == true) {
					stress1_surfXY[0] = 0.; stress1_surfXY[1] = 0.; stress1_surfXY[2] = 0.;
					stress1_surfXY[3] = 0.; stress1_surfXY[4] = 0.; stress1_surfXY[5] = 0.;
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = duzdz_glb_surfXY(i, j, k) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k) * 2. + \
						duzdy_glb_surfXZ(i, j, k) + duzdy_glb_surfXZ(i, j + 1, k)) / 4. - Deyz0_glb(i, j, k)) * 2.;

					strain[4] = ((duxdz_glb_surfXY(i, j, k) * 2. + \
						duzdx_glb_surfYZ(i, j, k) + duzdx_glb_surfYZ(i + 1, j, k)) / 4. - Dexz0_glb(i, j, k)) * 2.;

					strain[5] = (Dexy_glb(i, j, k) - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress1_surfXY);

					if (k == 0) {
						mat_type = (pt_glb->material_cell(i, j, nz - 1));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, nz - 1) - Dexx0_glb(i, j, nz - 1);
						strain[1] = Deyy_glb(i, j, nz - 1) - Deyy0_glb(i, j, nz - 1);
						strain[2] = duzdz_glb_surfXY(i, j, k) - Dezz0_glb(i, j, nz - 1);

						strain[3] = ((duydz_glb_surfXY(i, j, k) * 2. + \
							duzdy_glb_surfXZ(i, j, nz - 1) + duzdy_glb_surfXZ(i, j + 1, nz - 1)) / 4. - Deyz0_glb(i, j, nz - 1)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k) * 2. + \
							duzdx_glb_surfYZ(i, j, nz - 1) + duzdx_glb_surfYZ(i + 1, j, nz - 1)) / 4. - Dexz0_glb(i, j, nz - 1)) * 2.;

						strain[5] = (Dexy_glb(i, j, nz - 1) - Dexy0_glb(i, j, nz - 1)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j, k - 1));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, k - 1) - Dexx0_glb(i, j, k - 1);
						strain[1] = Deyy_glb(i, j, k - 1) - Deyy0_glb(i, j, k - 1);
						strain[2] = duzdz_glb_surfXY(i, j, k) - Dezz0_glb(i, j, k - 1);

						strain[3] = ((duydz_glb_surfXY(i, j, k) * 2. + \
							duzdy_glb_surfXZ(i, j, k - 1) + duzdy_glb_surfXZ(i, j + 1, k - 1)) / 4. - Deyz0_glb(i, j, k - 1)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k) * 2. + \
							duzdx_glb_surfYZ(i, j, k - 1) + duzdx_glb_surfYZ(i + 1, j, k - 1)) / 4. - Dexz0_glb(i, j, k - 1)) * 2.;

						strain[5] = (Dexy_glb(i, j, k - 1) - Dexy0_glb(i, j, k - 1)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress1_surfXY[0] = (stress1_surfXY[0] + stress[0]) / 2.;
					stress1_surfXY[1] = (stress1_surfXY[1] + stress[1]) / 2.;
					stress1_surfXY[2] = (stress1_surfXY[2] + stress[2]) / 2.;
					stress1_surfXY[3] = (stress1_surfXY[3] + stress[3]) / 2.;
					stress1_surfXY[4] = (stress1_surfXY[4] + stress[4]) / 2.;
					stress1_surfXY[5] = (stress1_surfXY[5] + stress[5]) / 2.;
				}

				//----------XY surface in +Z direction-----//
				if (stress_free_surfXY(i, j, k + 1) == true) {
					if (pt_glb->if_gaussian_strain_pulse == true && pt_glb->input_strain_type == 3 && k == source_z) {

						temporal_var = pt_glb->time_device - 5.0 * pt_glb->sigma_gauss;
						if (pt_glb->input_strain_component == 'x') {
							stress2_surfXY[0] = 0.;
							stress2_surfXY[1] = 0.;
							stress2_surfXY[2] = 0.;
							stress2_surfXY[3] = 0.;
							stress2_surfXY[4] = pt_glb->amplitude_gauss * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->sigma_gauss, 2.)));
							stress2_surfXY[5] = 0.;
						}
						else if (pt_glb->input_strain_component == 'y') {
							stress2_surfXY[0] = 0.;
							stress2_surfXY[1] = 0.;
							stress2_surfXY[2] = 0.;
							stress2_surfXY[3] = pt_glb->amplitude_gauss * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->sigma_gauss, 2.)));
							stress2_surfXY[4] = 0.;
							stress2_surfXY[5] = 0.;
						}
						else if (pt_glb->input_strain_component == 'z') {
							stress2_surfXY[0] = 0.;
							stress2_surfXY[1] = 0.;
							stress2_surfXY[2] = pt_glb->amplitude_gauss * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->sigma_gauss, 2.)));
							stress2_surfXY[3] = 0.;
							stress2_surfXY[4] = 0.;
							stress2_surfXY[5] = 0.;
						}

					}
					else {
						stress2_surfXY[0] = 0.; stress2_surfXY[1] = 0.; stress2_surfXY[2] = 0.;
						stress2_surfXY[3] = 0.; stress2_surfXY[4] = 0.; stress2_surfXY[5] = 0.;
					}
				}
				else {
					mat_type = (pt_glb->material_cell(i, j, k));
					mat = &(pt_glb->material_parameters[(mat_type)-1]);

					strain[0] = Dexx_glb(i, j, k) - Dexx0_glb(i, j, k);
					strain[1] = Deyy_glb(i, j, k) - Deyy0_glb(i, j, k);
					strain[2] = duzdz_glb_surfXY(i, j, k + 1) - Dezz0_glb(i, j, k);

					strain[3] = ((duydz_glb_surfXY(i, j, k + 1) * 2. + \
						duzdy_glb_surfXZ(i, j, k) + duzdy_glb_surfXZ(i, j + 1, k)) / 4. - Deyz0_glb(i, j, k)) * 2.;

					strain[4] = ((duxdz_glb_surfXY(i, j, k + 1) * 2. + \
						duzdx_glb_surfYZ(i, j, k) + duzdx_glb_surfYZ(i + 1, j, k)) / 4. - Dexz0_glb(i, j, k)) * 2.;

					strain[5] = (Dexy_glb(i, j, k) - Dexy0_glb(i, j, k)) * 2.;

					pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress2_surfXY);

					if (k == ny - 1) {
						mat_type = (pt_glb->material_cell(i, j, 0));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, 0) - Dexx0_glb(i, j, 0);
						strain[1] = Deyy_glb(i, j, 0) - Deyy0_glb(i, j, 0);
						strain[2] = duzdz_glb_surfXY(i, j, k + 1) - Dezz0_glb(i, j, 0);

						strain[3] = ((duydz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdy_glb_surfXZ(i, j, 0) + duzdy_glb_surfXZ(i, j + 1, 0)) / 4. - Deyz0_glb(i, j, 0)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdx_glb_surfYZ(i, j, 0) + duzdx_glb_surfYZ(i + 1, j, 0)) / 4. - Dexz0_glb(i, j, 0)) * 2.;

						strain[5] = (Dexy_glb(i, j, 0) - Dexy0_glb(i, j, 0)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					else {
						mat_type = (pt_glb->material_cell(i, j, k + 1));
						mat = &(pt_glb->material_parameters[(mat_type)-1]);

						strain[0] = Dexx_glb(i, j, k + 1) - Dexx0_glb(i, j, k + 1);
						strain[1] = Deyy_glb(i, j, k + 1) - Deyy0_glb(i, j, k + 1);
						strain[2] = duzdz_glb_surfXY(i, j, k + 1) - Dezz0_glb(i, j, k + 1);

						strain[3] = ((duydz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdy_glb_surfXZ(i, j, k + 1) + duzdy_glb_surfXZ(i, j + 1, k + 1)) / 4. - Deyz0_glb(i, j, k + 1)) * 2.;

						strain[4] = ((duxdz_glb_surfXY(i, j, k + 1) * 2. + \
							duzdx_glb_surfYZ(i, j, k + 1) + duzdx_glb_surfYZ(i + 1, j, k + 1)) / 4. - Dexz0_glb(i, j, k + 1)) * 2.;

						strain[5] = (Dexy_glb(i, j, k + 1) - Dexy0_glb(i, j, k + 1)) * 2.;

						pt_math->strain2stress_glb(mat->cijkl_glb, strain, stress);
					}
					stress2_surfXY[0] = (stress2_surfXY[0] + stress[0]) / 2.;
					stress2_surfXY[1] = (stress2_surfXY[1] + stress[1]) / 2.;
					stress2_surfXY[2] = (stress2_surfXY[2] + stress[2]) / 2.;
					stress2_surfXY[3] = (stress2_surfXY[3] + stress[3]) / 2.;
					stress2_surfXY[4] = (stress2_surfXY[4] + stress[4]) / 2.;
					stress2_surfXY[5] = (stress2_surfXY[5] + stress[5]) / 2.;
				}
				//----------END XY surface in +-Z direction-----//

				mat_type = (pt_glb->material_cell(i, j, k));
				mat = &(pt_glb->material_parameters[(mat_type)-1]);

				density = pt_glb->dt / mat->density;
				mass_damping = -pt_glb->dt * mat->elast_mass_damping;
				stiff_damping = mat->elast_stiff_damping / mat->density;

				force_x(i, j, k) = (stress2_surfYZ[0] - stress1_surfYZ[0]) / dx + \
					(stress2_surfXZ[5] - stress1_surfXZ[5]) / dy + \
					(stress2_surfXY[4] - stress1_surfXY[4]) / dz;

				force_y(i, j, k) = (stress2_surfYZ[5] - stress1_surfYZ[5]) / dx + \
					(stress2_surfXZ[1] - stress1_surfXZ[1]) / dy + \
					(stress2_surfXY[3] - stress1_surfXY[3]) / dz;

				force_z(i, j, k) = (stress2_surfYZ[4] - stress1_surfYZ[4]) / dx + \
					(stress2_surfXZ[3] - stress1_surfXZ[3]) / dy + \
					(stress2_surfXY[2] - stress1_surfXY[2]) / dz;

				dvx_glb_rk4(i, j, k) = force_x(i, j, k) * density \
					+ mass_damping * vx_glb(i, j, k) \
					+ stiff_damping * (force_x(i, j, k) - force_x_store(i, j, k));

				dvy_glb_rk4(i, j, k) = force_y(i, j, k) * density \
					+ mass_damping * vy_glb(i, j, k) \
					+ stiff_damping * (force_y(i, j, k) - force_y_store(i, j, k));

				dvz_glb_rk4(i, j, k) = force_z(i, j, k) * density \
					+ mass_damping * vz_glb(i, j, k) \
					+ stiff_damping * (force_z(i, j, k) - force_z_store(i, j, k));
			}
		}
	}
}

//void elastic_system::get_du() {
//#pragma acc parallel default(present)
//	{
//#pragma acc loop gang vector
//		for (long int id = 0; id < n; id++) {
//			dDux_glb(id) = vx_glb(id) * pt_glb->dt;
//			dDuy_glb(id) = vy_glb(id) * pt_glb->dt;
//			dDuz_glb(id) = vz_glb(id) * pt_glb->dt;
//		}
//	}
//}

void elastic_system::get_du_RK1() {
#pragma acc parallel default(present) async(3)
	{
#pragma acc loop gang vector
		for (long int id = 0; id < n; id++) {
			dDux_glb_rk1(id) = vx_glb_store(id) * pt_glb->dt;
			dDuy_glb_rk1(id) = vy_glb_store(id) * pt_glb->dt;
			dDuz_glb_rk1(id) = vz_glb_store(id) * pt_glb->dt;
		}
	}
}

void elastic_system::get_du_RK2() {
#pragma acc parallel default(present) async(3)
	{
#pragma acc loop gang vector
		for (long int id = 0; id < n; id++) {
			dDux_glb_rk2(id) = vx_glb_store(id) * pt_glb->dt;
			dDuy_glb_rk2(id) = vy_glb_store(id) * pt_glb->dt;
			dDuz_glb_rk2(id) = vz_glb_store(id) * pt_glb->dt;
		}
	}
}

void elastic_system::get_du_RK3() {
#pragma acc parallel default(present) async(3)
	{
#pragma acc loop gang vector
		for (long int id = 0; id < n; id++) {
			dDux_glb_rk3(id) = vx_glb_store(id) * pt_glb->dt;
			dDuy_glb_rk3(id) = vy_glb_store(id) * pt_glb->dt;
			dDuz_glb_rk3(id) = vz_glb_store(id) * pt_glb->dt;
		}
	}
}

void elastic_system::get_du_RK4() {
#pragma acc parallel default(present) async(3)
	{
#pragma acc loop gang vector
		for (long int id = 0; id < n; id++) {
			dDux_glb_rk4(id) = vx_glb_store(id) * pt_glb->dt;
			dDuy_glb_rk4(id) = vy_glb_store(id) * pt_glb->dt;
			dDuz_glb_rk4(id) = vz_glb_store(id) * pt_glb->dt;
		}
	}
}

void elastic_system::update_v_RK1() {
#pragma acc parallel default(present) async(1)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			vx_glb_store(id) = vx_glb(id) + dvx_glb_rk1(id) * 0.5;
			vy_glb_store(id) = vy_glb(id) + dvy_glb_rk1(id) * 0.5;
			vz_glb_store(id) = vz_glb(id) + dvz_glb_rk1(id) * 0.5;
		}
	}
}

void elastic_system::update_v_RK2() {
#pragma acc parallel default(present) async(1)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			vx_glb_store(id) = vx_glb(id) + dvx_glb_rk2(id) * 0.5;
			vy_glb_store(id) = vy_glb(id) + dvy_glb_rk2(id) * 0.5;
			vz_glb_store(id) = vz_glb(id) + dvz_glb_rk2(id) * 0.5;
		}
	}
}

void elastic_system::update_v_RK3() {
#pragma acc parallel default(present) async(1)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			vx_glb_store(id) = vx_glb(id) + dvx_glb_rk3(id);
			vy_glb_store(id) = vy_glb(id) + dvy_glb_rk3(id);
			vz_glb_store(id) = vz_glb(id) + dvz_glb_rk3(id);
		}
	}
}

void elastic_system::update_v() {
#pragma acc parallel default(present) async(1)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			vx_glb(id) = vx_glb(id) + dvx_glb_rk1(id) / 6. + dvx_glb_rk2(id) / 3. + dvx_glb_rk3(id) / 3. + dvx_glb_rk4(id) / 6.;
			vx_glb_store(id) = vx_glb(id);

			vy_glb(id) = vy_glb(id) + dvy_glb_rk1(id) / 6. + dvy_glb_rk2(id) / 3. + dvy_glb_rk3(id) / 3. + dvy_glb_rk4(id) / 6.;
			vy_glb_store(id) = vy_glb(id);

			vz_glb(id) = vz_glb(id) + dvz_glb_rk1(id) / 6. + dvz_glb_rk2(id) / 3. + dvz_glb_rk3(id) / 3. + dvz_glb_rk4(id) / 6.;
			vz_glb_store(id) = vz_glb(id);
		}
	}
}

void elastic_system::update_u_RK1() {

#pragma acc parallel default(present) async(2)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			Dux_glb_store(id) = Dux_glb(id) + dDux_glb_rk1(id) * 0.5;
			Duy_glb_store(id) = Duy_glb(id) + dDuy_glb_rk1(id) * 0.5;
			Duz_glb_store(id) = Duz_glb(id) + dDuz_glb_rk1(id) * 0.5;
		}
	}

#pragma acc parallel default(present) async(2)
	{
		update_Du_Boundary_half();
		if (pt_glb->if_gaussian_strain_pulse == true) {
			update_Du_Gaussian_half();
		}
	}

#pragma acc parallel default(present) async(2)
	{
		get_dudxyz_glb();
	}
#pragma acc parallel default(present) async(2)
	{
		get_Dstrain();
	}

	if (pt_glb->if_FE_all == true && pt_glb->if_flexo == true) {
#pragma acc parallel default(present) async(2)
		{
			get_Dstrain_gradient();
		}
	}
}

void elastic_system::update_u_RK2() {

#pragma acc parallel default(present) async(2)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			Dux_glb_store(id) = Dux_glb(id) + dDux_glb_rk2(id) * 0.5;
			Duy_glb_store(id) = Duy_glb(id) + dDuy_glb_rk2(id) * 0.5;
			Duz_glb_store(id) = Duz_glb(id) + dDuz_glb_rk2(id) * 0.5;
		}
	}

#pragma acc parallel default(present) async(2)
	{
		update_Du_Boundary_half();
		if (pt_glb->if_gaussian_strain_pulse == true) {
			update_Du_Gaussian_half();
		}
	}

#pragma acc parallel default(present) async(2)
	{
		get_dudxyz_glb();
	}
#pragma acc parallel default(present) async(2)
	{
		get_Dstrain();
	}

	if (pt_glb->if_FE_all == true && pt_glb->if_flexo == true) {
#pragma acc parallel default(present) async(2)
		{
			get_Dstrain_gradient();
		}
	}
}

void elastic_system::update_u_RK3() {

#pragma acc parallel default(present) async(2)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			Dux_glb_store(id) = Dux_glb(id) + dDux_glb_rk3(id);
			Duy_glb_store(id) = Duy_glb(id) + dDuy_glb_rk3(id);
			Duz_glb_store(id) = Duz_glb(id) + dDuz_glb_rk3(id);
		}
	}

#pragma acc parallel default(present) async(2)
	{
		update_Du_Boundary_full();
		if (pt_glb->if_gaussian_strain_pulse == true) {
			update_Du_Gaussian_full();
		}
	}

#pragma acc parallel default(present) async(2)
	{
		get_dudxyz_glb();
	}
#pragma acc parallel default(present) async(2)
	{
		get_Dstrain();
	}

	if (pt_glb->if_FE_all == true && pt_glb->if_flexo == true) {
#pragma acc parallel default(present) async(2)
		{
			get_Dstrain_gradient();
		}
	}
}

void elastic_system::update_u() {

	transfer_pointer();

#pragma acc parallel default(present) async(2)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			Dux_glb(id) = Dux_glb(id) + dDux_glb_rk1(id) / 6. + dDux_glb_rk2(id) / 3. + dDux_glb_rk3(id) / 3. + dDux_glb_rk4(id) / 6.;
			Duy_glb(id) = Duy_glb(id) + dDuy_glb_rk1(id) / 6. + dDuy_glb_rk2(id) / 3. + dDuy_glb_rk3(id) / 3. + dDuy_glb_rk4(id) / 6.;
			Duz_glb(id) = Duz_glb(id) + dDuz_glb_rk1(id) / 6. + dDuz_glb_rk2(id) / 3. + dDuz_glb_rk3(id) / 3. + dDuz_glb_rk4(id) / 6.;
		}
	}

#pragma acc parallel default(present) async(2)
	{
		update_Du_Boundary();
		if (pt_glb->if_gaussian_strain_pulse == true) {
			update_Du_Gaussian();
		}
	}

#pragma acc parallel default(present) async(2)
	{
		Dux_glb_store = Dux_glb;
		Duy_glb_store = Duy_glb;
		Duz_glb_store = Duz_glb;
	}

#pragma acc parallel default(present) async(2)
	{
		get_dudxyz_glb();
	}
#pragma acc parallel default(present) async(2)
	{
		get_Dstrain();
	}

	if (pt_glb->if_FE_all == true && pt_glb->if_flexo == true) {
#pragma acc parallel default(present) async(2)
		{
			get_Dstrain_gradient();
		}
	}
}

#pragma acc routine gang// nohost
void elastic_system::update_Du_Gaussian_half() {
	unsigned int i, j;

#pragma acc loop gang vector private(i,j)
	for (long int id = 0; id < nx * ny; id++) {
		i = id / ny;
		j = id - i * ny;
		if (pt_glb->material_cell(i, j, source_z) != 0) {
			if (pt_glb->input_strain_type == 1 && pt_glb->time_device < 12.5 * pt_glb->sigma_gauss) {

				if (pt_glb->input_strain_component == 'x') {
					Dux_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * 1.6487210 * pt_glb->sigma_gauss * \
						vs_tran_1st * exp(-1. * pow((pt_glb->time_device - 0.5 * pt_glb->dt - 5.0 * pt_glb->sigma_gauss), 2.) / (2. * pt_glb->sigma_gauss * pt_glb->sigma_gauss));
				}
				else if (pt_glb->input_strain_component == 'y') {
					Duy_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * 1.6487210 * pt_glb->sigma_gauss * \
						vs_tran_1st * exp(-1. * pow((pt_glb->time_device - 0.5 * pt_glb->dt - 5.0 * pt_glb->sigma_gauss), 2.) / (2. * pt_glb->sigma_gauss * pt_glb->sigma_gauss));
				}
				else if (pt_glb->input_strain_component == 'z') {
					Duz_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * 1.6487210 * pt_glb->sigma_gauss * \
						vs_long_1st * exp(-1. * pow((pt_glb->time_device - 0.5 * pt_glb->dt - 5.0 * pt_glb->sigma_gauss), 2.) / (2. * pt_glb->sigma_gauss * pt_glb->sigma_gauss));
				}

			}
			else if (pt_glb->input_strain_type == 2) {
				if (pt_glb->num_strain_cycles == 0) {
					if (pt_glb->input_strain_component == 'x') {
						Dux_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * vs_tran_1st / 2. / PI / pt_glb->sigma_gauss * \
							(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device - 0.5 * pt_glb->dt)) - 1);
					}
					else if (pt_glb->input_strain_component == 'y') {
						Duy_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * vs_tran_1st / 2. / PI / pt_glb->sigma_gauss * \
							(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device - 0.5 * pt_glb->dt)) - 1);
					}
					else if (pt_glb->input_strain_component == 'z') {
						Duz_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * vs_long_1st / 2. / PI / pt_glb->sigma_gauss * \
							(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device - 0.5 * pt_glb->dt)) - 1);
					}
				}
				else {
					if (pt_glb->input_strain_component == 'x') {
						if ((pt_glb->time_device - 0.5 * pt_glb->dt) <= static_cast<double>(pt_glb->num_strain_cycles) / pt_glb->sigma_gauss) {
							Dux_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * vs_tran_1st / 2. / PI / pt_glb->sigma_gauss * \
								(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device - 0.5 * pt_glb->dt)) - 1);
						}
						else {
							Dux_glb_store(i, j, source_z) = 0.;
						}
					}
					else if (pt_glb->input_strain_component == 'y') {
						if ((pt_glb->time_device - 0.5 * pt_glb->dt) <= static_cast<double>(pt_glb->num_strain_cycles) / pt_glb->sigma_gauss) {
							Duy_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * vs_tran_1st / 2. / PI / pt_glb->sigma_gauss * \
								(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device - 0.5 * pt_glb->dt)) - 1);
						}
						else {
							Duy_glb_store(i, j, source_z) = 0.;
						}
					}
					else if (pt_glb->input_strain_component == 'z') {
						if ((pt_glb->time_device - 0.5 * pt_glb->dt) <= static_cast<double>(pt_glb->num_strain_cycles) / pt_glb->sigma_gauss) {
							Duz_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * vs_long_1st / 2. / PI / pt_glb->sigma_gauss * \
								(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device - 0.5 * pt_glb->dt)) - 1);
						}
						else {
							Duz_glb_store(i, j, source_z) = 0.;
						}
					}
				}

			}
		}
	}

}

#pragma acc routine gang// nohost
void elastic_system::update_Du_Gaussian_full() {
	unsigned int i, j;

#pragma acc loop gang vector private(i,j)
	for (long int id = 0; id < nx * ny; id++) {
		i = id / ny;
		j = id - i * ny;
		if (pt_glb->material_cell(i, j, source_z) != 0) {
			if (pt_glb->input_strain_type == 1 && pt_glb->time_device < 12.5 * pt_glb->sigma_gauss) {

				if (pt_glb->input_strain_component == 'x') {
					Dux_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * 1.6487210 * pt_glb->sigma_gauss * \
						vs_tran_1st * exp(-1. * pow((pt_glb->time_device - 5.0 * pt_glb->sigma_gauss), 2.) / (2. * pt_glb->sigma_gauss * pt_glb->sigma_gauss));
				}
				else if (pt_glb->input_strain_component == 'y') {
					Duy_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * 1.6487210 * pt_glb->sigma_gauss * \
						vs_tran_1st * exp(-1. * pow((pt_glb->time_device - 5.0 * pt_glb->sigma_gauss), 2.) / (2. * pt_glb->sigma_gauss * pt_glb->sigma_gauss));
				}
				else if (pt_glb->input_strain_component == 'z') {
					Duz_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * 1.6487210 * pt_glb->sigma_gauss * \
						vs_long_1st * exp(-1. * pow((pt_glb->time_device - 5.0 * pt_glb->sigma_gauss), 2.) / (2. * pt_glb->sigma_gauss * pt_glb->sigma_gauss));
				}

			}
			else if (pt_glb->input_strain_type == 2) {
				if (pt_glb->num_strain_cycles > 0) {
					if (pt_glb->input_strain_component == 'x') {
						if ((pt_glb->time_device) <= static_cast<double>(pt_glb->num_strain_cycles) / pt_glb->sigma_gauss) {
							Dux_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * vs_tran_1st / 2. / PI / pt_glb->sigma_gauss * \
								(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device)) - 1);
						}
						else {
							Dux_glb_store(i, j, source_z) = 0.;
						}
					}
					else if (pt_glb->input_strain_component == 'y') {
						if ((pt_glb->time_device) <= static_cast<double>(pt_glb->num_strain_cycles) / pt_glb->sigma_gauss) {
							Duy_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * vs_tran_1st / 2. / PI / pt_glb->sigma_gauss * \
								(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device)) - 1);
						}
						else {
							Duy_glb_store(i, j, source_z) = 0.;
						}
					}
					else if (pt_glb->input_strain_component == 'z') {
						if ((pt_glb->time_device) <= static_cast<double>(pt_glb->num_strain_cycles) / pt_glb->sigma_gauss) {
							Duz_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * vs_long_1st / 2. / PI / pt_glb->sigma_gauss * \
								(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device)) - 1);
						}
						else {
							Duz_glb_store(i, j, source_z) = 0.;
						}
					}
				}
				else {
					if (pt_glb->input_strain_component == 'x') {
						Dux_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * vs_tran_1st / 2. / PI / pt_glb->sigma_gauss * \
							(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device)) - 1);
					}
					else if (pt_glb->input_strain_component == 'y') {
						Duy_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * vs_tran_1st / 2. / PI / pt_glb->sigma_gauss * \
							(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device)) - 1);
					}
					else if (pt_glb->input_strain_component == 'z') {
						Duz_glb_store(i, j, source_z) = pt_glb->amplitude_gauss * vs_long_1st / 2. / PI / pt_glb->sigma_gauss * \
							(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device)) - 1);
					}
				}
			}
		}
	}

}

#pragma acc routine gang// nohost
void elastic_system::update_Du_Gaussian() {
	unsigned int i, j;

#pragma acc loop gang vector private(i,j)
	for (long int id = 0; id < nx * ny; id++) {
		i = id / ny;
		j = id - i * ny;
		if (pt_glb->material_cell(i, j, source_z) != 0) {
			if (pt_glb->input_strain_type == 1 && pt_glb->time_device < 12.5 * pt_glb->sigma_gauss) {

				if (pt_glb->input_strain_component == 'x') {
					Dux_glb(i, j, source_z) = pt_glb->amplitude_gauss * 1.6487210 * pt_glb->sigma_gauss * \
						vs_tran_1st * exp(-1. * pow((pt_glb->time_device - 5.0 * pt_glb->sigma_gauss), 2.) / (2. * pt_glb->sigma_gauss * pt_glb->sigma_gauss));
				}
				else if (pt_glb->input_strain_component == 'y') {
					Duy_glb(i, j, source_z) = pt_glb->amplitude_gauss * 1.6487210 * pt_glb->sigma_gauss * \
						vs_tran_1st * exp(-1. * pow((pt_glb->time_device - 5.0 * pt_glb->sigma_gauss), 2.) / (2. * pt_glb->sigma_gauss * pt_glb->sigma_gauss));
				}
				else if (pt_glb->input_strain_component == 'z') {
					Duz_glb(i, j, source_z) = pt_glb->amplitude_gauss * 1.6487210 * pt_glb->sigma_gauss * \
						vs_long_1st * exp(-1. * pow((pt_glb->time_device - 5.0 * pt_glb->sigma_gauss), 2.) / (2. * pt_glb->sigma_gauss * pt_glb->sigma_gauss));
				}

			}
			else if (pt_glb->input_strain_type == 2) {
				if (pt_glb->num_strain_cycles == 0) {
					if (pt_glb->input_strain_component == 'x') {
						Dux_glb(i, j, source_z) = pt_glb->amplitude_gauss * vs_tran_1st / 2. / PI / pt_glb->sigma_gauss * \
							(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device)) - 1);
					}
					else if (pt_glb->input_strain_component == 'y') {
						Duy_glb(i, j, source_z) = pt_glb->amplitude_gauss * vs_tran_1st / 2. / PI / pt_glb->sigma_gauss * \
							(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device)) - 1);
					}
					else if (pt_glb->input_strain_component == 'z') {
						Duz_glb(i, j, source_z) = pt_glb->amplitude_gauss * vs_long_1st / 2. / PI / pt_glb->sigma_gauss * \
							(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device)) - 1);
					}
				}
				else {
					if (pt_glb->input_strain_component == 'x') {
						if ((pt_glb->time_device) <= static_cast<double>(pt_glb->num_strain_cycles) / pt_glb->sigma_gauss) {
							Dux_glb(i, j, source_z) = pt_glb->amplitude_gauss * vs_tran_1st / 2. / PI / pt_glb->sigma_gauss * \
								(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device)) - 1);
						}
						else {
							Dux_glb(i, j, source_z) = 0.;
						}
					}
					else if (pt_glb->input_strain_component == 'y') {
						if ((pt_glb->time_device) <= static_cast<double>(pt_glb->num_strain_cycles) / pt_glb->sigma_gauss) {
							Duy_glb(i, j, source_z) = pt_glb->amplitude_gauss * vs_tran_1st / 2. / PI / pt_glb->sigma_gauss * \
								(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device)) - 1);
						}
						else {
							Duy_glb(i, j, source_z) = 0.;
						}
					}
					else if (pt_glb->input_strain_component == 'z') {
						if ((pt_glb->time_device) <= static_cast<double>(pt_glb->num_strain_cycles) / pt_glb->sigma_gauss) {
							Duz_glb(i, j, source_z) = pt_glb->amplitude_gauss * vs_long_1st / 2. / PI / pt_glb->sigma_gauss * \
								(cos(2. * PI * pt_glb->sigma_gauss * (pt_glb->time_device)) - 1);
						}
						else {
							Duz_glb(i, j, source_z) = 0.;
						}
					}
				}

			}
		}
	}

}
