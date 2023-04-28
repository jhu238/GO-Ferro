#include "magnetic_system.h"

//void magnetic_system::initialize(elastic_system* pt_elas_in, EMdynamic_system* pt_EM_in) {
void magnetic_system::initialize_host() {
	//pt_elas = pt_elas_in;
	//pt_EM = pt_EM_in;

	nx = pt_geo->nx_system;
	ny = pt_geo->ny_system;
	nz = pt_geo->nz_system;

	n = nx * ny * nz;

	dx = pt_geo->dx;
	dy = pt_geo->dy;
	dz = pt_geo->dz;

	mx_glb.initialize(nx, ny, nz); my_glb.initialize(nx, ny, nz); mz_glb.initialize(nx, ny, nz);

	mx_AFM1_glb.initialize(nx, ny, nz); my_AFM1_glb.initialize(nx, ny, nz); mz_AFM1_glb.initialize(nx, ny, nz);
	mx_AFM2_glb.initialize(nx, ny, nz); my_AFM2_glb.initialize(nx, ny, nz); mz_AFM2_glb.initialize(nx, ny, nz);

	if (pt_glb->if_FM_all == true && pt_glb->if_input_m == false) {
		initialize_magnetization();
	}
	if (pt_glb->if_AFM_all == true && pt_glb->if_input_AFMm == false) {
		initialize_magnetization_AFM();
	}

	Hx_stat.initialize(nx, ny, nz); Hy_stat.initialize(nx, ny, nz); Hz_stat.initialize(nx, ny, nz);
	if (pt_glb->if_mag_1Dmodel == false && pt_glb->if_magnetostatics == true && pt_glb->if_demag_factor == false) {
		initialize_intermediate_Hms();
	}

	dmx_glb.initialize(nx, ny, nz); dmy_glb.initialize(nx, ny, nz); dmz_glb.initialize(nx, ny, nz);
	//Hx_anis.initialize(nx, ny, nz); Hy_anis.initialize(nx, ny, nz); Hz_anis.initialize(nx, ny, nz);
	//Hx_exch.initialize(nx, ny, nz); Hy_exch.initialize(nx, ny, nz); Hz_exch.initialize(nx, ny, nz);
	//Hx_elas.initialize(nx, ny, nz); Hy_elas.initialize(nx, ny, nz); Hz_elas.initialize(nx, ny, nz);

	dmx_AFM1_glb.initialize(nx, ny, nz); dmy_AFM1_glb.initialize(nx, ny, nz); dmz_AFM1_glb.initialize(nx, ny, nz);
	dmx_AFM2_glb.initialize(nx, ny, nz); dmy_AFM2_glb.initialize(nx, ny, nz); dmz_AFM2_glb.initialize(nx, ny, nz);

	//-----------------RK4-------------------//
	mx_glb_store.initialize(nx, ny, nz); my_glb_store.initialize(nx, ny, nz); mz_glb_store.initialize(nx, ny, nz);
	dmx_glb_rk1.initialize(nx, ny, nz); dmy_glb_rk1.initialize(nx, ny, nz); dmz_glb_rk1.initialize(nx, ny, nz);
	dmx_glb_rk2.initialize(nx, ny, nz); dmy_glb_rk2.initialize(nx, ny, nz); dmz_glb_rk2.initialize(nx, ny, nz);
	dmx_glb_rk3.initialize(nx, ny, nz); dmy_glb_rk3.initialize(nx, ny, nz); dmz_glb_rk3.initialize(nx, ny, nz);
	dmx_glb_rk4.initialize(nx, ny, nz); dmy_glb_rk4.initialize(nx, ny, nz); dmz_glb_rk4.initialize(nx, ny, nz);

	mx_AFM1_glb_store.initialize(nx, ny, nz); my_AFM1_glb_store.initialize(nx, ny, nz); mz_AFM1_glb_store.initialize(nx, ny, nz);
	mx_AFM2_glb_store.initialize(nx, ny, nz); my_AFM2_glb_store.initialize(nx, ny, nz); mz_AFM2_glb_store.initialize(nx, ny, nz);

	dmx_AFM1_glb_rk1.initialize(nx, ny, nz); dmy_AFM1_glb_rk1.initialize(nx, ny, nz); dmz_AFM1_glb_rk1.initialize(nx, ny, nz);
	dmx_AFM1_glb_rk2.initialize(nx, ny, nz); dmy_AFM1_glb_rk2.initialize(nx, ny, nz); dmz_AFM1_glb_rk2.initialize(nx, ny, nz);
	dmx_AFM1_glb_rk3.initialize(nx, ny, nz); dmy_AFM1_glb_rk3.initialize(nx, ny, nz); dmz_AFM1_glb_rk3.initialize(nx, ny, nz);
	dmx_AFM1_glb_rk4.initialize(nx, ny, nz); dmy_AFM1_glb_rk4.initialize(nx, ny, nz); dmz_AFM1_glb_rk4.initialize(nx, ny, nz);

	dmx_AFM2_glb_rk1.initialize(nx, ny, nz); dmy_AFM2_glb_rk1.initialize(nx, ny, nz); dmz_AFM2_glb_rk1.initialize(nx, ny, nz);
	dmx_AFM2_glb_rk2.initialize(nx, ny, nz); dmy_AFM2_glb_rk2.initialize(nx, ny, nz); dmz_AFM2_glb_rk2.initialize(nx, ny, nz);
	dmx_AFM2_glb_rk3.initialize(nx, ny, nz); dmy_AFM2_glb_rk3.initialize(nx, ny, nz); dmz_AFM2_glb_rk3.initialize(nx, ny, nz);
	dmx_AFM2_glb_rk4.initialize(nx, ny, nz); dmy_AFM2_glb_rk4.initialize(nx, ny, nz); dmz_AFM2_glb_rk4.initialize(nx, ny, nz);

	if (pt_glb->if_FM_all == true) {
		FM_surfYZ.initialize(nx + 1, ny, nz);
		FM_surfXZ.initialize(nx, ny + 1, nz);
		FM_surfXY.initialize(nx, ny, nz + 1);

		//----------Averaged parameter for inter-region exchange - Approach 1---------//
		Aij_xf.initialize(nx, ny, nz);
		Aij_yf.initialize(nx, ny, nz);
		Aij_zf.initialize(nx, ny, nz);
		Aij_xb.initialize(nx, ny, nz);
		Aij_yb.initialize(nx, ny, nz);
		Aij_zb.initialize(nx, ny, nz);
		//----------Averaged parameter for inter-region exchange - Approach 2---------//
		//AM_exchange.initialize(nx, ny, nz);
	}

	if (pt_glb->if_AFM_all == true) {
		AFM_surfYZ.initialize(nx + 1, ny, nz);
		AFM_surfXZ.initialize(nx, ny + 1, nz);
		AFM_surfXY.initialize(nx, ny, nz + 1);

		Aij_AFM1_xf.initialize(nx, ny, nz); Aij_AFM2_xf.initialize(nx, ny, nz);
		Aij_AFM1_yf.initialize(nx, ny, nz);	Aij_AFM2_yf.initialize(nx, ny, nz);
		Aij_AFM1_zf.initialize(nx, ny, nz);	Aij_AFM2_zf.initialize(nx, ny, nz);
		Aij_AFM1_xb.initialize(nx, ny, nz);	Aij_AFM2_xb.initialize(nx, ny, nz);
		Aij_AFM1_yb.initialize(nx, ny, nz);	Aij_AFM2_yb.initialize(nx, ny, nz);
		Aij_AFM1_zb.initialize(nx, ny, nz);	Aij_AFM2_zb.initialize(nx, ny, nz);
	}

	if (pt_glb->if_mag_1Dmodel == false && pt_glb->if_magnetostatics == true && pt_glb->if_demag_factor == false) {
		Axyzk.initialize(nx, ny, nz / 2 + 1);
		Bxyzk.initialize(nx, ny, nz / 2 + 1);
		Bxzyk.initialize(nx, ny, nz / 2 + 1);
		Ayzxk.initialize(nx, ny, nz / 2 + 1);
		Byzxk.initialize(nx, ny, nz / 2 + 1);
		Byxzk.initialize(nx, ny, nz / 2 + 1);
		Azxyk.initialize(nx, ny, nz / 2 + 1);
		Bzxyk.initialize(nx, ny, nz / 2 + 1);
		Bzyxk.initialize(nx, ny, nz / 2 + 1);

		Mx_glb.initialize(nx, ny, nz); My_glb.initialize(nx, ny, nz); Mz_glb.initialize(nx, ny, nz);
		Mx_k3D.initialize(nx, ny, nz / 2 + 1); My_k3D.initialize(nx, ny, nz / 2 + 1); Mz_k3D.initialize(nx, ny, nz / 2 + 1);
		Hx_stat_k.initialize(nx, ny, nz / 2 + 1); Hy_stat_k.initialize(nx, ny, nz / 2 + 1); Hz_stat_k.initialize(nx, ny, nz / 2 + 1);
	}

	Jx_ISHE.initialize(nx, ny, nz);
	Jy_ISHE.initialize(nx, ny, nz);

	initialize_host_supple();

	//----------Initialize cuHandle for cuFFT--------------//
	cufftPlan3d(&plan_magneto_D2Z, nx, ny, nz, CUFFT_D2Z);
	cufftPlan3d(&plan_magneto_Z2D, nx, ny, nz, CUFFT_Z2D);
}

void magnetic_system::initialize_host_supple() {
	unsigned int material_type1, material_type2;
	//unsigned int type_neighbor[6];
	//double MA_temp, count;
	bool if_F1, if_F2;
	unsigned int mat_type, mat_typen;
	material* mat, * matn;
	long x, y, z;
	double Aex1, Aex2;
	double Aex1_AFM1, Aex2_AFM1;
	double Aex1_AFM2, Aex2_AFM2;
	//double Ms1, Ms2;
	//#pragma acc declare device_resident(material_type1, material_type2,if_F1, if_F2)

	//------------------------------------------//
	//	Determine if FM/Non-FM boundaries		//
	//------------------------------------------//

	if (pt_glb->if_FM_all == true) {
		//#pragma acc parallel default(present)
		{
			//------------X Direction (YZ plane)----------------//
//#pragma acc loop gang
			for (long int j = 0; j < ny; j++) {
				//#pragma acc loop vector private(material_type1, material_type2,if_F1, if_F2)
				for (long int k = 0; k < nz; k++) {
					material_type1 = pt_glb->material_cell(0, j, k);
					if (material_type1 == 0) {
						if_F1 = false;
					}
					else {
						if_F1 = pt_glb->material_parameters[(material_type1)-1].if_FM;
					}

					material_type2 = pt_glb->material_cell(nx - 1, j, k);
					if (material_type2 == 0) {
						if_F2 = false;
					}
					else {
						if_F2 = pt_glb->material_parameters[(material_type2)-1].if_FM;
					}

					if (pt_geo->periodicX == true) {
						if (if_F1 ^ if_F2) {
							FM_surfYZ(0, j, k) = true;
							FM_surfYZ(nx, j, k) = true;
						}
						else {
							FM_surfYZ(0, j, k) = false;
							FM_surfYZ(nx, j, k) = false;
						}
					}
					else if (pt_geo->periodicX == false) {
						if (if_F1 == true) {
							FM_surfYZ(0, j, k) = true;
						}
						else {
							FM_surfYZ(0, j, k) = false;
						}

						if (if_F2 == true) {
							FM_surfYZ(nx, j, k) = true;
						}
						else {
							FM_surfYZ(nx, j, k) = false;
						}
					}
				}
			}

			//#pragma acc loop gang
			for (long int i = 1; i < nx; i++) {
				//#pragma acc loop worker
				for (long int j = 0; j < ny; j++) {
					//#pragma acc loop vector private(material_type1, material_type2,if_F1, if_F2)
					for (long int k = 0; k < nz; k++) {
						material_type1 = pt_glb->material_cell(i - 1, j, k);
						if (material_type1 == 0) {
							if_F1 = false;
						}
						else {
							if_F1 = pt_glb->material_parameters[(material_type1)-1].if_FM;
						}

						material_type2 = pt_glb->material_cell(i, j, k);
						if (material_type2 == 0) {
							if_F2 = false;
						}
						else {
							if_F2 = pt_glb->material_parameters[(material_type2)-1].if_FM;
						}

						if (if_F1 ^ if_F2) {
							FM_surfYZ(i, j, k) = true;
						}
						else {
							FM_surfYZ(i, j, k) = false;
						}
					}
				}
			}

			//------------Y Direction (XZ plane)----------------//
//#pragma acc loop gang
			for (long int i = 0; i < nx; i++) {
				//#pragma acc loop vector private(material_type1, material_type2,if_F1, if_F2)
				for (long int k = 0; k < nz; k++) {
					material_type1 = pt_glb->material_cell(i, 0, k);
					if (material_type1 == 0) {
						if_F1 = false;
					}
					else {
						if_F1 = pt_glb->material_parameters[(material_type1)-1].if_FM;
					}

					material_type2 = pt_glb->material_cell(i, ny - 1, k);
					if (material_type2 == 0) {
						if_F2 = false;
					}
					else {
						if_F2 = pt_glb->material_parameters[(material_type2)-1].if_FM;
					}

					if (pt_geo->periodicY == true) {
						if (if_F1 ^ if_F2) {
							FM_surfXZ(i, 0, k) = true;
							FM_surfXZ(i, ny, k) = true;
						}
						else {
							FM_surfXZ(i, 0, k) = false;
							FM_surfXZ(i, ny, k) = false;
						}
					}
					else if (pt_geo->periodicY == false) {
						if (if_F1 == true) {
							FM_surfXZ(i, 0, k) = true;
						}
						else {
							FM_surfXZ(i, 0, k) = false;
						}

						if (if_F2 == true) {
							FM_surfXZ(i, ny, k) = true;
						}
						else {
							FM_surfXZ(i, ny, k) = false;
						}
					}
				}
			}
			//#pragma acc loop gang
			for (long int i = 0; i < nx; i++) {
				//#pragma acc loop worker
				for (long int j = 1; j < ny; j++) {
					//#pragma acc loop vector private(material_type1, material_type2,if_F1, if_F2)
					for (long int k = 0; k < nz; k++) {
						material_type1 = pt_glb->material_cell(i, j - 1, k);
						if (material_type1 == 0) {
							if_F1 = false;
						}
						else {
							if_F1 = pt_glb->material_parameters[(material_type1)-1].if_FM;
						}

						material_type2 = pt_glb->material_cell(i, j, k);
						if (material_type2 == 0) {
							if_F2 = false;
						}
						else {
							if_F2 = pt_glb->material_parameters[(material_type2)-1].if_FM;
						}

						if (if_F1 ^ if_F2) {
							FM_surfXZ(i, j, k) = true;
						}
						else {
							FM_surfXZ(i, j, k) = false;
						}
					}
				}
			}

			//------------Z Direction (XY plane)----------------//
//#pragma acc loop gang
			for (long int i = 0; i < nx; i++) {
				//#pragma acc loop vector private(material_type1, material_type2,if_F1, if_F2)
				for (long int j = 0; j < ny; j++) {
					material_type1 = pt_glb->material_cell(i, j, 0);
					if (material_type1 == 0) {
						if_F1 = false;
					}
					else {
						if_F1 = pt_glb->material_parameters[(material_type1)-1].if_FM;
					}

					material_type2 = pt_glb->material_cell(i, j, nz - 1);
					if (material_type2 == 0) {
						if_F2 = false;
					}
					else {
						if_F2 = pt_glb->material_parameters[(material_type2)-1].if_FM;
					}

					if (pt_geo->periodicZ == true) {
						if (if_F1 ^ if_F2) {
							FM_surfXY(i, j, 0) = true;
							FM_surfXY(i, j, nz) = true;
						}
						else {
							FM_surfXY(i, j, 0) = false;
							FM_surfXY(i, j, nz) = false;
						}
					}
					else if (pt_geo->periodicZ == false) {
						if (if_F1 == true) {
							FM_surfXY(i, j, 0) = true;
						}
						else {
							FM_surfXY(i, j, 0) = false;
						}

						if (if_F2 == true) {
							FM_surfXY(i, j, nz) = true;
						}
						else {
							FM_surfXY(i, j, nz) = false;
						}
					}
				}
			}
			//#pragma acc loop gang
			for (long int i = 0; i < nx; i++) {
				//#pragma acc loop worker
				for (long int j = 0; j < ny; j++) {
					//#pragma acc loop vector private(material_type1, material_type2,if_F1, if_F2)
					for (long int k = 1; k < nz; k++) {
						material_type1 = pt_glb->material_cell(i, j, k - 1);
						if (material_type1 == 0) {
							if_F1 = false;
						}
						else {
							if_F1 = pt_glb->material_parameters[(material_type1)-1].if_FM;
						}

						material_type2 = pt_glb->material_cell(i, j, k);
						if (material_type2 == 0) {
							if_F2 = false;
						}
						else {
							if_F2 = pt_glb->material_parameters[(material_type2)-1].if_FM;
						}

						if (if_F1 ^ if_F2) {
							FM_surfXY(i, j, k) = true;
						}
						else {
							FM_surfXY(i, j, k) = false;
						}
					}
				}
			}
		}

		//----------Averaged parameter for inter-region exchange - Approach 1---------//
		for (long int id = 0; id < n; id++) {
			mat_type = pt_glb->material_cell(id);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_FM == true) {
					x = id / (ny * nz);
					y = (id - x * (ny * nz)) / nz;
					z = id - x * (ny * nz) - y * nz;

					if (FM_surfYZ(x + 1, y, z) == true) {
						Aij_xf(id) = mat->Aex;
					}
					else {
						if (x == nx - 1) {
							mat_typen = pt_glb->material_cell(0, y, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1 = mat->Aex; Aex2 = matn->Aex;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//Aij_xf(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_xf(id) = 2. * Aex1 * Aex2 / (Aex1 + Aex2);
						}
						else {
							mat_typen = pt_glb->material_cell(x + 1, y, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1 = mat->Aex; Aex2 = matn->Aex;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//Aij_xf(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_xf(id) = 2. * Aex1 * Aex2 / (Aex1 + Aex2);
						}
					}

					if (FM_surfYZ(x, y, z) == true) {
						Aij_xb(id) = mat->Aex;
					}
					else {
						if (x == 0) {
							mat_typen = pt_glb->material_cell(nx - 1, y, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1 = mat->Aex; Aex2 = matn->Aex;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_xb(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_xb(id) = 2. * Aex1 * Aex2 / (Aex1 + Aex2);
						}
						else {
							mat_typen = pt_glb->material_cell(x - 1, y, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1 = mat->Aex; Aex2 = matn->Aex;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_xb(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_xb(id) = 2. * Aex1 * Aex2 / (Aex1 + Aex2);
						}
					}

					if (FM_surfXZ(x, y + 1, z) == true) {
						Aij_yf(id) = mat->Aex;
					}
					else {
						if (y == ny - 1) {
							mat_typen = pt_glb->material_cell(x, 0, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1 = mat->Aex; Aex2 = matn->Aex;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_yf(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_yf(id) = 2. * Aex1 * Aex2 / (Aex1 + Aex2);
						}
						else {
							mat_typen = pt_glb->material_cell(x, y + 1, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1 = mat->Aex; Aex2 = matn->Aex;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_yf(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_yf(id) = 2. * Aex1 * Aex2 / (Aex1 + Aex2);
						}
					}

					if (FM_surfXZ(x, y, z) == true) {
						Aij_yb(id) = mat->Aex;
					}
					else {
						if (y == 0) {
							mat_typen = pt_glb->material_cell(x, ny - 1, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1 = mat->Aex; Aex2 = matn->Aex;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_yb(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_yb(id) = 2. * Aex1 * Aex2 / (Aex1 + Aex2);
						}
						else {
							mat_typen = pt_glb->material_cell(x, y - 1, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1 = mat->Aex; Aex2 = matn->Aex;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_yb(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_yb(id) = 2. * Aex1 * Aex2 / (Aex1 + Aex2);
						}
					}

					if (FM_surfXY(x, y, z + 1) == true) {
						Aij_zf(id) = mat->Aex;
					}
					else {
						if (z == nz - 1) {
							mat_typen = pt_glb->material_cell(x, y, 0);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1 = mat->Aex; Aex2 = matn->Aex;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_zf(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_zf(id) = 2. * Aex1 * Aex2 / (Aex1 + Aex2);
						}
						else {
							mat_typen = pt_glb->material_cell(x, y, z + 1);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1 = mat->Aex; Aex2 = matn->Aex;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_zf(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_zf(id) = 2. * Aex1 * Aex2 / (Aex1 + Aex2);
						}
					}

					if (FM_surfXY(x, y, z) == true) {
						Aij_zb(id) = mat->Aex;
					}
					else {
						if (z == 0) {
							mat_typen = pt_glb->material_cell(x, y, nz - 1);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1 = mat->Aex; Aex2 = matn->Aex;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_zb(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_zb(id) = 2. * Aex1 * Aex2 / (Aex1 + Aex2);
						}
						else {
							mat_typen = pt_glb->material_cell(x, y, z - 1);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1 = mat->Aex; Aex2 = matn->Aex;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_zb(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_zb(id) = 2. * Aex1 * Aex2 / (Aex1 + Aex2);
						}
					}
				}
			}
		}

		//----------Averaged parameter for inter-region exchange - Approach 2---------//
				//for (long id = 0; id < n; id++) {
				//	mat_type = pt_glb->material_cell(id);
				//	if (mat_type == 0) {
				//		AM_exchange(id) = 0.;
				//	}
				//	else {
				//		x = id / (ny * nz);
				//		y = (id - x * (ny * nz)) / nz;
				//		z = id - x * (ny * nz) - y * nz;
				//		mat= &(pt_glb->material_parameters[(mat_type)-1]);

				//		//-X
				//		if (FM_surfYZ(x, y, z) == true) {
				//			type_neighbor[0] = 0;
				//		}
				//		else {
				//			if (x == 0) {
				//				type_neighbor[0] = pt_glb->material_cell(nx - 1, y, z);
				//			}
				//			else {
				//				type_neighbor[0] = pt_glb->material_cell(x - 1, y, z);
				//			}
				//		}

				//		//+X
				//		if (FM_surfYZ(x + 1, y, z) == true) {
				//			type_neighbor[1] = 0;
				//		}
				//		else {
				//			if (x == nx - 1) {
				//				type_neighbor[1] = pt_glb->material_cell(0, y, z);
				//			}
				//			else {
				//				type_neighbor[1] = pt_glb->material_cell(x + 1, y, z);
				//			}
				//		}

				//		//-Y
				//		if (FM_surfXZ(x, y, z) == true) {
				//			type_neighbor[2] = 0;
				//		}
				//		else {
				//			if (y == 0) {
				//				type_neighbor[2] = pt_glb->material_cell(x, ny-1, z);
				//			}
				//			else {
				//				type_neighbor[2] = pt_glb->material_cell(x, y-1, z);
				//			}
				//		}

				//		//+Y
				//		if (FM_surfXZ(x, y+1, z) == true) {
				//			type_neighbor[3] = 0;
				//		}
				//		else {
				//			if (y == ny - 1) {
				//				type_neighbor[3] = pt_glb->material_cell(x, 0, z);
				//			}
				//			else {
				//				type_neighbor[3] = pt_glb->material_cell(x, y+1, z);
				//			}
				//		}

				//		//-Z
				//		if (FM_surfXY(x, y, z) == true) {
				//			type_neighbor[4] = 0;
				//		}
				//		else {
				//			if (z == 0) {
				//				type_neighbor[4] = pt_glb->material_cell(x, y, nz-1);
				//			}
				//			else {
				//				type_neighbor[4] = pt_glb->material_cell(x, y, z-1);
				//			}
				//		}

				//		//+Z
				//		if (FM_surfXY(x, y, z+1) == true) {
				//			type_neighbor[5] = 0;
				//		}
				//		else {
				//			if (z == nz - 1) {
				//				type_neighbor[5] = pt_glb->material_cell(x, y, 0);
				//			}
				//			else {
				//				type_neighbor[5] = pt_glb->material_cell(x, y, z+1);
				//			}
				//		}

				//		MA_temp = mat->Ms / mat->Aex; count = 1.;
				//		for (long i = 0; i < 6; i++) {
				//			if (type_neighbor[i] != 0 && type_neighbor[i] != mat_type) {
				//				matn = &(pt_glb->material_parameters[type_neighbor[i] - 1]);
				//				MA_temp = MA_temp + matn->Ms / matn->Aex;
				//				count = count + 1.;
				//			}
				//		}

				//		AM_exchange(id) = count / MA_temp;

				//	}
				//}
	}

	//------------------------------------------//
	//	Determine if AFM/Non-AFM boundaries		//
	//------------------------------------------//

	if (pt_glb->if_AFM_all == true) {
		{
			//------------X Direction (YZ plane)----------------//
			for (long int j = 0; j < ny; j++) {
				for (long int k = 0; k < nz; k++) {
					material_type1 = pt_glb->material_cell(0, j, k);
					if (material_type1 == 0) {
						if_F1 = false;
					}
					else {
						if_F1 = pt_glb->material_parameters[(material_type1)-1].if_AFM;
					}

					material_type2 = pt_glb->material_cell(nx - 1, j, k);
					if (material_type2 == 0) {
						if_F2 = false;
					}
					else {
						if_F2 = pt_glb->material_parameters[(material_type2)-1].if_AFM;
					}

					if (pt_geo->periodicX == true) {
						if (if_F1 ^ if_F2) {
							AFM_surfYZ(0, j, k) = true;
							AFM_surfYZ(nx, j, k) = true;
						}
						else {
							AFM_surfYZ(0, j, k) = false;
							AFM_surfYZ(nx, j, k) = false;
						}
					}
					else if (pt_geo->periodicX == false) {
						if (if_F1 == true) {
							AFM_surfYZ(0, j, k) = true;
						}
						else {
							AFM_surfYZ(0, j, k) = false;
						}

						if (if_F2 == true) {
							AFM_surfYZ(nx, j, k) = true;
						}
						else {
							AFM_surfYZ(nx, j, k) = false;
						}
					}
				}
			}

			for (long int i = 1; i < nx; i++) {
				for (long int j = 0; j < ny; j++) {
					for (long int k = 0; k < nz; k++) {
						material_type1 = pt_glb->material_cell(i - 1, j, k);
						if (material_type1 == 0) {
							if_F1 = false;
						}
						else {
							if_F1 = pt_glb->material_parameters[(material_type1)-1].if_AFM;
						}

						material_type2 = pt_glb->material_cell(i, j, k);
						if (material_type2 == 0) {
							if_F2 = false;
						}
						else {
							if_F2 = pt_glb->material_parameters[(material_type2)-1].if_AFM;
						}

						if (if_F1 ^ if_F2) {
							AFM_surfYZ(i, j, k) = true;
						}
						else {
							AFM_surfYZ(i, j, k) = false;
						}
					}
				}
			}

			//------------Y Direction (XZ plane)----------------//
			for (long int i = 0; i < nx; i++) {
				for (long int k = 0; k < nz; k++) {
					material_type1 = pt_glb->material_cell(i, 0, k);
					if (material_type1 == 0) {
						if_F1 = false;
					}
					else {
						if_F1 = pt_glb->material_parameters[(material_type1)-1].if_AFM;
					}

					material_type2 = pt_glb->material_cell(i, ny - 1, k);
					if (material_type2 == 0) {
						if_F2 = false;
					}
					else {
						if_F2 = pt_glb->material_parameters[(material_type2)-1].if_AFM;
					}

					if (pt_geo->periodicY == true) {
						if (if_F1 ^ if_F2) {
							AFM_surfXZ(i, 0, k) = true;
							AFM_surfXZ(i, ny, k) = true;
						}
						else {
							AFM_surfXZ(i, 0, k) = false;
							AFM_surfXZ(i, ny, k) = false;
						}
					}
					else if (pt_geo->periodicY == false) {
						if (if_F1 == true) {
							AFM_surfXZ(i, 0, k) = true;
						}
						else {
							AFM_surfXZ(i, 0, k) = false;
						}

						if (if_F2 == true) {
							AFM_surfXZ(i, ny, k) = true;
						}
						else {
							AFM_surfXZ(i, ny, k) = false;
						}
					}
				}
			}

			for (long int i = 0; i < nx; i++) {
				for (long int j = 1; j < ny; j++) {
					for (long int k = 0; k < nz; k++) {
						material_type1 = pt_glb->material_cell(i, j - 1, k);
						if (material_type1 == 0) {
							if_F1 = false;
						}
						else {
							if_F1 = pt_glb->material_parameters[(material_type1)-1].if_AFM;
						}

						material_type2 = pt_glb->material_cell(i, j, k);
						if (material_type2 == 0) {
							if_F2 = false;
						}
						else {
							if_F2 = pt_glb->material_parameters[(material_type2)-1].if_AFM;
						}

						if (if_F1 ^ if_F2) {
							AFM_surfXZ(i, j, k) = true;
						}
						else {
							AFM_surfXZ(i, j, k) = false;
						}
					}
				}
			}

			//------------Z Direction (XY plane)----------------//
			for (long int i = 0; i < nx; i++) {
				for (long int j = 0; j < ny; j++) {
					material_type1 = pt_glb->material_cell(i, j, 0);
					if (material_type1 == 0) {
						if_F1 = false;
					}
					else {
						if_F1 = pt_glb->material_parameters[(material_type1)-1].if_AFM;
					}

					material_type2 = pt_glb->material_cell(i, j, nz - 1);
					if (material_type2 == 0) {
						if_F2 = false;
					}
					else {
						if_F2 = pt_glb->material_parameters[(material_type2)-1].if_AFM;
					}

					if (pt_geo->periodicZ == true) {
						if (if_F1 ^ if_F2) {
							AFM_surfXY(i, j, 0) = true;
							AFM_surfXY(i, j, nz) = true;
						}
						else {
							AFM_surfXY(i, j, 0) = false;
							AFM_surfXY(i, j, nz) = false;
						}
					}
					else if (pt_geo->periodicZ == false) {
						if (if_F1 == true) {
							AFM_surfXY(i, j, 0) = true;
						}
						else {
							AFM_surfXY(i, j, 0) = false;
						}

						if (if_F2 == true) {
							AFM_surfXY(i, j, nz) = true;
						}
						else {
							AFM_surfXY(i, j, nz) = false;
						}
					}
				}
			}

			for (long int i = 0; i < nx; i++) {
				for (long int j = 0; j < ny; j++) {
					for (long int k = 1; k < nz; k++) {
						material_type1 = pt_glb->material_cell(i, j, k - 1);
						if (material_type1 == 0) {
							if_F1 = false;
						}
						else {
							if_F1 = pt_glb->material_parameters[(material_type1)-1].if_AFM;
						}

						material_type2 = pt_glb->material_cell(i, j, k);
						if (material_type2 == 0) {
							if_F2 = false;
						}
						else {
							if_F2 = pt_glb->material_parameters[(material_type2)-1].if_AFM;
						}

						if (if_F1 ^ if_F2) {
							AFM_surfXY(i, j, k) = true;
						}
						else {
							AFM_surfXY(i, j, k) = false;
						}
					}
				}
			}
		}

		//----------Averaged parameter for inter-region exchange - Approach 1---------//
		for (long int id = 0; id < n; id++) {
			mat_type = pt_glb->material_cell(id);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_AFM == true) {
					x = id / (ny * nz);
					y = (id - x * (ny * nz)) / nz;
					z = id - x * (ny * nz) - y * nz;

					if (AFM_surfYZ(x + 1, y, z) == true) {
						Aij_AFM1_xf(id) = mat->Aex_AFM1;
						Aij_AFM2_xf(id) = mat->Aex_AFM2;
					}
					else {
						if (x == nx - 1) {
							mat_typen = pt_glb->material_cell(0, y, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1_AFM1 = mat->Aex_AFM1; Aex2_AFM1 = matn->Aex_AFM1;
							Aex1_AFM2 = mat->Aex_AFM2; Aex2_AFM2 = matn->Aex_AFM2;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//Aij_xf(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_AFM1_xf(id) = 2. * Aex1_AFM1 * Aex2_AFM1 / (Aex1_AFM1 + Aex2_AFM1);
							Aij_AFM2_xf(id) = 2. * Aex1_AFM2 * Aex2_AFM2 / (Aex1_AFM2 + Aex2_AFM2);
						}
						else {
							mat_typen = pt_glb->material_cell(x + 1, y, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1_AFM1 = mat->Aex_AFM1; Aex2_AFM1 = matn->Aex_AFM1;
							Aex1_AFM2 = mat->Aex_AFM2; Aex2_AFM2 = matn->Aex_AFM2;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//Aij_xf(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_AFM1_xf(id) = 2. * Aex1_AFM1 * Aex2_AFM1 / (Aex1_AFM1 + Aex2_AFM1);
							Aij_AFM2_xf(id) = 2. * Aex1_AFM2 * Aex2_AFM2 / (Aex1_AFM2 + Aex2_AFM2);
						}
					}

					if (AFM_surfYZ(x, y, z) == true) {
						Aij_AFM1_xb(id) = mat->Aex_AFM1;
						Aij_AFM2_xb(id) = mat->Aex_AFM2;
					}
					else {
						if (x == 0) {
							mat_typen = pt_glb->material_cell(nx - 1, y, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1_AFM1 = mat->Aex_AFM1; Aex2_AFM1 = matn->Aex_AFM1;
							Aex1_AFM2 = mat->Aex_AFM2; Aex2_AFM2 = matn->Aex_AFM2;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_xb(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_AFM1_xb(id) = 2. * Aex1_AFM1 * Aex2_AFM1 / (Aex1_AFM1 + Aex2_AFM1);
							Aij_AFM2_xb(id) = 2. * Aex1_AFM2 * Aex2_AFM2 / (Aex1_AFM2 + Aex2_AFM2);
						}
						else {
							mat_typen = pt_glb->material_cell(x - 1, y, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1_AFM1 = mat->Aex_AFM1; Aex2_AFM1 = matn->Aex_AFM1;
							Aex1_AFM2 = mat->Aex_AFM2; Aex2_AFM2 = matn->Aex_AFM2;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_xb(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_AFM1_xb(id) = 2. * Aex1_AFM1 * Aex2_AFM1 / (Aex1_AFM1 + Aex2_AFM1);
							Aij_AFM2_xb(id) = 2. * Aex1_AFM2 * Aex2_AFM2 / (Aex1_AFM2 + Aex2_AFM2);
						}
					}

					if (AFM_surfXZ(x, y + 1, z) == true) {
						Aij_AFM1_yf(id) = mat->Aex_AFM1;
						Aij_AFM2_yf(id) = mat->Aex_AFM2;
					}
					else {
						if (y == ny - 1) {
							mat_typen = pt_glb->material_cell(x, 0, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1_AFM1 = mat->Aex_AFM1; Aex2_AFM1 = matn->Aex_AFM1;
							Aex1_AFM2 = mat->Aex_AFM2; Aex2_AFM2 = matn->Aex_AFM2;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_yf(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_AFM1_yf(id) = 2. * Aex1_AFM1 * Aex2_AFM1 / (Aex1_AFM1 + Aex2_AFM1);
							Aij_AFM2_yf(id) = 2. * Aex1_AFM2 * Aex2_AFM2 / (Aex1_AFM2 + Aex2_AFM2);
						}
						else {
							mat_typen = pt_glb->material_cell(x, y + 1, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1_AFM1 = mat->Aex_AFM1; Aex2_AFM1 = matn->Aex_AFM1;
							Aex1_AFM2 = mat->Aex_AFM2; Aex2_AFM2 = matn->Aex_AFM2;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_yf(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_AFM1_yf(id) = 2. * Aex1_AFM1 * Aex2_AFM1 / (Aex1_AFM1 + Aex2_AFM1);
							Aij_AFM2_yf(id) = 2. * Aex1_AFM2 * Aex2_AFM2 / (Aex1_AFM2 + Aex2_AFM2);
						}
					}

					if (AFM_surfXZ(x, y, z) == true) {
						Aij_AFM1_yb(id) = mat->Aex_AFM1;
						Aij_AFM2_yb(id) = mat->Aex_AFM2;
					}
					else {
						if (y == 0) {
							mat_typen = pt_glb->material_cell(x, ny - 1, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1_AFM1 = mat->Aex_AFM1; Aex2_AFM1 = matn->Aex_AFM1;
							Aex1_AFM2 = mat->Aex_AFM2; Aex2_AFM2 = matn->Aex_AFM2;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_yb(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_AFM1_yb(id) = 2. * Aex1_AFM1 * Aex2_AFM1 / (Aex1_AFM1 + Aex2_AFM1);
							Aij_AFM2_yb(id) = 2. * Aex1_AFM2 * Aex2_AFM2 / (Aex1_AFM2 + Aex2_AFM2);
						}
						else {
							mat_typen = pt_glb->material_cell(x, y - 1, z);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1_AFM1 = mat->Aex_AFM1; Aex2_AFM1 = matn->Aex_AFM1;
							Aex1_AFM2 = mat->Aex_AFM2; Aex2_AFM2 = matn->Aex_AFM2;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_yb(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_AFM1_yb(id) = 2. * Aex1_AFM1 * Aex2_AFM1 / (Aex1_AFM1 + Aex2_AFM1);
							Aij_AFM2_yb(id) = 2. * Aex1_AFM2 * Aex2_AFM2 / (Aex1_AFM2 + Aex2_AFM2);
						}
					}

					if (AFM_surfXY(x, y, z + 1) == true) {
						Aij_AFM1_zf(id) = mat->Aex_AFM1;
						Aij_AFM2_zf(id) = mat->Aex_AFM2;
					}
					else {
						if (z == nz - 1) {
							mat_typen = pt_glb->material_cell(x, y, 0);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1_AFM1 = mat->Aex_AFM1; Aex2_AFM1 = matn->Aex_AFM1;
							Aex1_AFM2 = mat->Aex_AFM2; Aex2_AFM2 = matn->Aex_AFM2;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_zf(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_AFM1_zf(id) = 2. * Aex1_AFM1 * Aex2_AFM1 / (Aex1_AFM1 + Aex2_AFM1);
							Aij_AFM2_zf(id) = 2. * Aex1_AFM2 * Aex2_AFM2 / (Aex1_AFM2 + Aex2_AFM2);
						}
						else {
							mat_typen = pt_glb->material_cell(x, y, z + 1);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1_AFM1 = mat->Aex_AFM1; Aex2_AFM1 = matn->Aex_AFM1;
							Aex1_AFM2 = mat->Aex_AFM2; Aex2_AFM2 = matn->Aex_AFM2;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_zf(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_AFM1_zf(id) = 2. * Aex1_AFM1 * Aex2_AFM1 / (Aex1_AFM1 + Aex2_AFM1);
							Aij_AFM2_zf(id) = 2. * Aex1_AFM2 * Aex2_AFM2 / (Aex1_AFM2 + Aex2_AFM2);
						}
					}

					if (AFM_surfXY(x, y, z) == true) {
						Aij_AFM1_zb(id) = mat->Aex_AFM1;
						Aij_AFM2_zb(id) = mat->Aex_AFM2;
					}
					else {
						if (z == 0) {
							mat_typen = pt_glb->material_cell(x, y, nz - 1);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1_AFM1 = mat->Aex_AFM1; Aex2_AFM1 = matn->Aex_AFM1;
							Aex1_AFM2 = mat->Aex_AFM2; Aex2_AFM2 = matn->Aex_AFM2;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_zb(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_AFM1_zb(id) = 2. * Aex1_AFM1 * Aex2_AFM1 / (Aex1_AFM1 + Aex2_AFM1);
							Aij_AFM2_zb(id) = 2. * Aex1_AFM2 * Aex2_AFM2 / (Aex1_AFM2 + Aex2_AFM2);
						}
						else {
							mat_typen = pt_glb->material_cell(x, y, z - 1);
							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
							Aex1_AFM1 = mat->Aex_AFM1; Aex2_AFM1 = matn->Aex_AFM1;
							Aex1_AFM2 = mat->Aex_AFM2; Aex2_AFM2 = matn->Aex_AFM2;
							//Ms1 = mat->Ms; Ms2 = matn->Ms;
							//AM_zb(id) = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
							Aij_AFM1_zb(id) = 2. * Aex1_AFM1 * Aex2_AFM1 / (Aex1_AFM1 + Aex2_AFM1);
							Aij_AFM2_zb(id) = 2. * Aex1_AFM2 * Aex2_AFM2 / (Aex1_AFM2 + Aex2_AFM2);
						}
					}
				}
			}
		}
	}
	//	if (pt_glb->if_mag_1Dmodel == false) {
	//
	//		cufftPlan3d(&plan, nx, ny, nz, CUFFT_D2Z);
	//
	//		cufftExecD2Z(plan, Axyz.matrix, (cufftDoubleComplex*)(Axyzk.matrix));
	//		cufftExecD2Z(plan, Ayzx.matrix, (cufftDoubleComplex*)(Ayzxk.matrix));
	//		cufftExecD2Z(plan, Azxy.matrix, (cufftDoubleComplex*)(Azxyk.matrix));
	//		cufftExecD2Z(plan, Bxyz.matrix, (cufftDoubleComplex*)(Bxyzk.matrix));
	//		cufftExecD2Z(plan, Byzx.matrix, (cufftDoubleComplex*)(Byzxk.matrix));
	//		cufftExecD2Z(plan, Bzxy.matrix, (cufftDoubleComplex*)(Bzxyk.matrix));
	//		cufftExecD2Z(plan, Byxz.matrix, (cufftDoubleComplex*)(Byxzk.matrix));
	//		cufftExecD2Z(plan, Bzyx.matrix, (cufftDoubleComplex*)(Bzyxk.matrix));
	//		cufftExecD2Z(plan, Bxzy.matrix, (cufftDoubleComplex*)(Bxzyk.matrix));
	//
	//		cufftDestroy(plan);
	//
	//#pragma acc wait
	//
	//		Axyzk *= (1. / static_cast<double>(nx * ny * nz));
	//		Ayzxk *= (1. / static_cast<double>(nx * ny * nz));
	//		Azxyk *= (1. / static_cast<double>(nx * ny * nz));
	//
	//		Bxyzk *= (1. / static_cast<double>(nx * ny * nz));
	//		Byzxk *= (1. / static_cast<double>(nx * ny * nz));
	//		Bzxyk *= (1. / static_cast<double>(nx * ny * nz));
	//		Byxzk *= (1. / static_cast<double>(nx * ny * nz));
	//		Bzyxk *= (1. / static_cast<double>(nx * ny * nz));
	//		Bxzyk *= (1. / static_cast<double>(nx * ny * nz));
	//	}
}

void magnetic_system::initialize_device()
{
	if (pt_glb->if_prescribed_m == true) {
#pragma acc serial default(present) async(1)
		{
			update_prescribed_m();
		}
	}

#pragma acc parallel default(present) async(1)
	{
		mx_glb_store = mx_glb;
		my_glb_store = my_glb;
		mz_glb_store = mz_glb;
	}

#pragma acc parallel default(present) async(1)
	{
		mx_AFM1_glb_store = mx_AFM1_glb;
		my_AFM1_glb_store = my_AFM1_glb;
		mz_AFM1_glb_store = mz_AFM1_glb;

		mx_AFM2_glb_store = mx_AFM2_glb;
		my_AFM2_glb_store = my_AFM2_glb;
		mz_AFM2_glb_store = mz_AFM2_glb;
	}

	if (pt_glb->if_mag_1Dmodel == false && pt_glb->if_magnetostatics == true && pt_glb->if_demag_factor == false) {
#pragma acc host_data use_device(Axyz.matrix, Ayzx.matrix, Azxy.matrix, Bxyz.matrix, Byzx.matrix, Bzxy.matrix, Byxz.matrix, Bzyx.matrix, Bxzy.matrix, \
								 Axyzk.matrix,Ayzxk.matrix,Azxyk.matrix,Bxyzk.matrix,Byzxk.matrix,Bzxyk.matrix,Byxzk.matrix,Bzyxk.matrix,Bxzyk.matrix)
	{
		cufftExecD2Z(plan_magneto_D2Z, Axyz.matrix, (cufftDoubleComplex*)(Axyzk.matrix));
		cufftExecD2Z(plan_magneto_D2Z, Ayzx.matrix, (cufftDoubleComplex*)(Ayzxk.matrix));
		cufftExecD2Z(plan_magneto_D2Z, Azxy.matrix, (cufftDoubleComplex*)(Azxyk.matrix));
		cufftExecD2Z(plan_magneto_D2Z, Bxyz.matrix, (cufftDoubleComplex*)(Bxyzk.matrix));
		cufftExecD2Z(plan_magneto_D2Z, Byzx.matrix, (cufftDoubleComplex*)(Byzxk.matrix));
		cufftExecD2Z(plan_magneto_D2Z, Bzxy.matrix, (cufftDoubleComplex*)(Bzxyk.matrix));
		cufftExecD2Z(plan_magneto_D2Z, Byxz.matrix, (cufftDoubleComplex*)(Byxzk.matrix));
		cufftExecD2Z(plan_magneto_D2Z, Bzyx.matrix, (cufftDoubleComplex*)(Bzyxk.matrix));
		cufftExecD2Z(plan_magneto_D2Z, Bxzy.matrix, (cufftDoubleComplex*)(Bxzyk.matrix));
	}
#pragma acc wait

//#pragma acc parallel default(present) async(1)
//	{
//		Axyzk *= (1. / static_cast<double>(n));
//		Ayzxk *= (1. / static_cast<double>(n));
//		Azxyk *= (1. / static_cast<double>(n));
//
//		Bxyzk *= (1. / static_cast<double>(n));
//		Byzxk *= (1. / static_cast<double>(n));
//		Bzxyk *= (1. / static_cast<double>(n));
//		Byxzk *= (1. / static_cast<double>(n));
//		Bzyxk *= (1. / static_cast<double>(n));
//		Bxzyk *= (1. / static_cast<double>(n));
//	}
	}
}

void magnetic_system::initialize_magnetization() {
	unsigned int material_type;
	double theta = (pt_glb->theta_mag) / 180. * PI;
	double phi = (pt_glb->phi_mag) / 180 * PI;
	double mx = sin(theta) * cos(phi);
	double my = sin(theta) * sin(phi);
	double mz = cos(theta);

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				material_type = pt_glb->material_cell(i, j, k);
				if (material_type != 0) {
					if (pt_glb->material_parameters[(material_type - 1)].if_FM == true) {
						mx_glb(i, j, k) = mx; my_glb(i, j, k) = my; mz_glb(i, j, k) = mz;
					}
				}
			}
		}
	}
}

void magnetic_system::initialize_magnetization_AFM() {
	unsigned int material_type;
	double theta = (pt_glb->theta_mag) / 180. * PI;
	double phi = (pt_glb->phi_mag) / 180 * PI;
	double mx = sin(theta) * cos(phi);
	double my = sin(theta) * sin(phi);
	double mz = cos(theta);

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				material_type = pt_glb->material_cell(i, j, k);
				if (material_type != 0) {
					if (pt_glb->material_parameters[(material_type - 1)].if_AFM == true) {
						mx_AFM1_glb(i, j, k) = mx;  my_AFM1_glb(i, j, k) = my;  mz_AFM1_glb(i, j, k) = mz;
						mx_AFM2_glb(i, j, k) = -mx; my_AFM2_glb(i, j, k) = -my; mz_AFM2_glb(i, j, k) = -mz;
					}
				}
			}
		}
	}
}

void magnetic_system::initialize_intermediate_Hms() {
	double x, y, z;
	double twopi = 2. * PI;
	long int shiftx = nx / 2 - 1;
	long int shifty = ny / 2 - 1;
	long int shiftz = nz / 2 - 1;
	matrix3d<double> WAxyz, WAyzx, WAzxy, WBxyz, WByzx, WBzxy, WByxz, WBzyx, WBxzy;

	//Axyzk.initialize(nx, ny, nz / 2 + 1);
	//Bxyzk.initialize(nx, ny, nz / 2 + 1);
	//Bxzyk.initialize(nx, ny, nz / 2 + 1);
	//Ayzxk.initialize(nx, ny, nz / 2 + 1);
	//Byzxk.initialize(nx, ny, nz / 2 + 1);
	//Byxzk.initialize(nx, ny, nz / 2 + 1);
	//Azxyk.initialize(nx, ny, nz / 2 + 1);
	//Bzxyk.initialize(nx, ny, nz / 2 + 1);
	//Bzyxk.initialize(nx, ny, nz / 2 + 1);

	WAxyz.initialize(nx - 1, ny - 1, nz - 1); Axyz.initialize(nx, ny, nz);
	WAyzx.initialize(nx - 1, ny - 1, nz - 1); Ayzx.initialize(nx, ny, nz);
	WAzxy.initialize(nx - 1, ny - 1, nz - 1); Azxy.initialize(nx, ny, nz);
	WBxyz.initialize(nx - 1, ny - 1, nz - 1); Bxyz.initialize(nx, ny, nz);
	WByzx.initialize(nx - 1, ny - 1, nz - 1); Byzx.initialize(nx, ny, nz);
	WBzxy.initialize(nx - 1, ny - 1, nz - 1); Bzxy.initialize(nx, ny, nz);
	WByxz.initialize(nx - 1, ny - 1, nz - 1); Byxz.initialize(nx, ny, nz);
	WBzyx.initialize(nx - 1, ny - 1, nz - 1); Bzyx.initialize(nx, ny, nz);
	WBxzy.initialize(nx - 1, ny - 1, nz - 1); Bxzy.initialize(nx, ny, nz);

	for (long int i = 1 - (nx) / 2; i < (nx) / 2; i++) {
		for (long int j = 1 - (ny) / 2; j < (ny) / 2; j++) {
			for (long int k = 1 - (nz) / 2; k < (nz) / 2; k++) {
				x = static_cast<double>(i) * dx;
				y = static_cast<double>(j) * dy;
				z = static_cast<double>(k) * dz;

				WAxyz(i + shiftx, j + shifty, k + shiftz) = A(x, y, z, dx, dy, dz) / (2. * twopi);
				WAyzx(i + shiftx, j + shifty, k + shiftz) = A(y, z, x, dy, dz, dx) / (2. * twopi);
				WAzxy(i + shiftx, j + shifty, k + shiftz) = A(z, x, y, dz, dx, dy) / (2. * twopi);

				WBxyz(i + shiftx, j + shifty, k + shiftz) = B(x, y, z, dx, dy, dz) / (2. * twopi);
				WByzx(i + shiftx, j + shifty, k + shiftz) = B(y, z, x, dy, dz, dx) / (2. * twopi);
				WBzxy(i + shiftx, j + shifty, k + shiftz) = B(z, x, y, dz, dx, dy) / (2. * twopi);
				WByxz(i + shiftx, j + shifty, k + shiftz) = B(y, x, z, dy, dx, dz) / (2. * twopi);
				WBzyx(i + shiftx, j + shifty, k + shiftz) = B(z, y, x, dz, dy, dx) / (2. * twopi);
				WBxzy(i + shiftx, j + shifty, k + shiftz) = B(x, z, y, dx, dz, dy) / (2. * twopi);
			}
		}
	}

	for (long i = 1; i < nx / 2 + 1; i++) {
		for (long j = 1; j < ny / 2 + 1; j++) {
			for (long k = 1; k < nz / 2 + 1; k++) {
				Axyz(i - 1, j - 1, k - 1) = WAxyz(i + shiftx - 1, j + shifty - 1, k + shiftz - 1);
				Ayzx(i - 1, j - 1, k - 1) = WAyzx(i + shiftx - 1, j + shifty - 1, k + shiftz - 1);
				Azxy(i - 1, j - 1, k - 1) = WAzxy(i + shiftx - 1, j + shifty - 1, k + shiftz - 1);

				Bxyz(i - 1, j - 1, k - 1) = WBxyz(i + shiftx - 1, j + shifty - 1, k + shiftz - 1);
				Byzx(i - 1, j - 1, k - 1) = WByzx(i + shiftx - 1, j + shifty - 1, k + shiftz - 1);
				Bzxy(i - 1, j - 1, k - 1) = WBzxy(i + shiftx - 1, j + shifty - 1, k + shiftz - 1);
				Byxz(i - 1, j - 1, k - 1) = WByxz(i + shiftx - 1, j + shifty - 1, k + shiftz - 1);
				Bzyx(i - 1, j - 1, k - 1) = WBzyx(i + shiftx - 1, j + shifty - 1, k + shiftz - 1);
				Bxzy(i - 1, j - 1, k - 1) = WBxzy(i + shiftx - 1, j + shifty - 1, k + shiftz - 1);
			}
		}
	}

	for (long i = 1; i < nx / 2; i++) {
		for (long j = 1; j < ny / 2 + 1; j++) {
			for (long k = 1; k < nz / 2 + 1; k++) {
				Axyz(nx - i, j - 1, k - 1) = WAxyz(shiftx - i, j - 1 + shifty, k - 1 + shiftz);
				Ayzx(nx - i, j - 1, k - 1) = WAyzx(shiftx - i, j - 1 + shifty, k - 1 + shiftz);
				Azxy(nx - i, j - 1, k - 1) = WAzxy(shiftx - i, j - 1 + shifty, k - 1 + shiftz);
				Bxyz(nx - i, j - 1, k - 1) = WBxyz(shiftx - i, j - 1 + shifty, k - 1 + shiftz);
				Byzx(nx - i, j - 1, k - 1) = WByzx(shiftx - i, j - 1 + shifty, k - 1 + shiftz);
				Bzxy(nx - i, j - 1, k - 1) = WBzxy(shiftx - i, j - 1 + shifty, k - 1 + shiftz);
				Byxz(nx - i, j - 1, k - 1) = WByxz(shiftx - i, j - 1 + shifty, k - 1 + shiftz);
				Bzyx(nx - i, j - 1, k - 1) = WBzyx(shiftx - i, j - 1 + shifty, k - 1 + shiftz);
				Bxzy(nx - i, j - 1, k - 1) = WBxzy(shiftx - i, j - 1 + shifty, k - 1 + shiftz);
			}
		}
	}

	for (long i = 1; i < nx / 2; i++) {
		for (long j = 1; j < ny / 2; j++) {
			for (long k = 1; k < nz / 2 + 1; k++) {
				Axyz(nx - i, ny - j, k - 1) = WAxyz(shiftx - i, shifty - j, k - 1 + shiftz);
				Ayzx(nx - i, ny - j, k - 1) = WAyzx(shiftx - i, shifty - j, k - 1 + shiftz);
				Azxy(nx - i, ny - j, k - 1) = WAzxy(shiftx - i, shifty - j, k - 1 + shiftz);
				Bxyz(nx - i, ny - j, k - 1) = WBxyz(shiftx - i, shifty - j, k - 1 + shiftz);
				Byzx(nx - i, ny - j, k - 1) = WByzx(shiftx - i, shifty - j, k - 1 + shiftz);
				Bzxy(nx - i, ny - j, k - 1) = WBzxy(shiftx - i, shifty - j, k - 1 + shiftz);
				Byxz(nx - i, ny - j, k - 1) = WByxz(shiftx - i, shifty - j, k - 1 + shiftz);
				Bzyx(nx - i, ny - j, k - 1) = WBzyx(shiftx - i, shifty - j, k - 1 + shiftz);
				Bxzy(nx - i, ny - j, k - 1) = WBxzy(shiftx - i, shifty - j, k - 1 + shiftz);
			}
		}
	}

	for (long i = 1; i < nx / 2 + 1; i++) {
		for (long j = 1; j < ny / 2; j++) {
			for (long k = 1; k < nz / 2 + 1; k++) {
				Axyz(i - 1, ny - j, k - 1) = WAxyz(i - 1 + shiftx, shifty - j, k - 1 + shiftz);
				Ayzx(i - 1, ny - j, k - 1) = WAyzx(i - 1 + shiftx, shifty - j, k - 1 + shiftz);
				Azxy(i - 1, ny - j, k - 1) = WAzxy(i - 1 + shiftx, shifty - j, k - 1 + shiftz);

				Bxyz(i - 1, ny - j, k - 1) = WBxyz(i - 1 + shiftx, shifty - j, k - 1 + shiftz);
				Byzx(i - 1, ny - j, k - 1) = WByzx(i - 1 + shiftx, shifty - j, k - 1 + shiftz);
				Bzxy(i - 1, ny - j, k - 1) = WBzxy(i - 1 + shiftx, shifty - j, k - 1 + shiftz);
				Byxz(i - 1, ny - j, k - 1) = WByxz(i - 1 + shiftx, shifty - j, k - 1 + shiftz);
				Bzyx(i - 1, ny - j, k - 1) = WBzyx(i - 1 + shiftx, shifty - j, k - 1 + shiftz);
				Bxzy(i - 1, ny - j, k - 1) = WBxzy(i - 1 + shiftx, shifty - j, k - 1 + shiftz);
			}
		}
	}

	for (long i = 1; i < nx / 2; i++) {
		for (long j = 1; j < ny / 2 + 1; j++) {
			for (long k = 1; k < nz / 2; k++) {
				Axyz(nx - i, j - 1, nz - k) = WAxyz(shiftx - i, j - 1 + shifty, shiftz - k);
				Ayzx(nx - i, j - 1, nz - k) = WAyzx(shiftx - i, j - 1 + shifty, shiftz - k);
				Azxy(nx - i, j - 1, nz - k) = WAzxy(shiftx - i, j - 1 + shifty, shiftz - k);

				Bxyz(nx - i, j - 1, nz - k) = WBxyz(shiftx - i, j - 1 + shifty, shiftz - k);
				Byzx(nx - i, j - 1, nz - k) = WByzx(shiftx - i, j - 1 + shifty, shiftz - k);
				Bzxy(nx - i, j - 1, nz - k) = WBzxy(shiftx - i, j - 1 + shifty, shiftz - k);
				Byxz(nx - i, j - 1, nz - k) = WByxz(shiftx - i, j - 1 + shifty, shiftz - k);
				Bzyx(nx - i, j - 1, nz - k) = WBzyx(shiftx - i, j - 1 + shifty, shiftz - k);
				Bxzy(nx - i, j - 1, nz - k) = WBxzy(shiftx - i, j - 1 + shifty, shiftz - k);
			}
		}
	}

	for (long i = 1; i < nx / 2 + 1; i++) {
		for (long j = 1; j < ny / 2 + 1; j++) {
			for (long k = 1; k < nz / 2; k++) {
				Axyz(i - 1, j - 1, nz - k) = WAxyz(i - 1 + shiftx, j - 1 + shifty, shiftz - k);
				Ayzx(i - 1, j - 1, nz - k) = WAyzx(i - 1 + shiftx, j - 1 + shifty, shiftz - k);
				Azxy(i - 1, j - 1, nz - k) = WAzxy(i - 1 + shiftx, j - 1 + shifty, shiftz - k);

				Bxyz(i - 1, j - 1, nz - k) = WBxyz(i - 1 + shiftx, j - 1 + shifty, shiftz - k);
				Byzx(i - 1, j - 1, nz - k) = WByzx(i - 1 + shiftx, j - 1 + shifty, shiftz - k);
				Bzxy(i - 1, j - 1, nz - k) = WBzxy(i - 1 + shiftx, j - 1 + shifty, shiftz - k);
				Byxz(i - 1, j - 1, nz - k) = WByxz(i - 1 + shiftx, j - 1 + shifty, shiftz - k);
				Bzyx(i - 1, j - 1, nz - k) = WBzyx(i - 1 + shiftx, j - 1 + shifty, shiftz - k);
				Bxzy(i - 1, j - 1, nz - k) = WBxzy(i - 1 + shiftx, j - 1 + shifty, shiftz - k);
			}
		}
	}

	for (long i = 1; i < nx / 2; i++) {
		for (long j = 1; j < ny / 2; j++) {
			for (long k = 1; k < nz / 2; k++) {
				Axyz(nx - i, ny - j, nz - k) = WAxyz(shiftx - i, shifty - j, shiftz - k);
				Ayzx(nx - i, ny - j, nz - k) = WAyzx(shiftx - i, shifty - j, shiftz - k);
				Azxy(nx - i, ny - j, nz - k) = WAzxy(shiftx - i, shifty - j, shiftz - k);

				Bxyz(nx - i, ny - j, nz - k) = WBxyz(shiftx - i, shifty - j, shiftz - k);
				Byzx(nx - i, ny - j, nz - k) = WByzx(shiftx - i, shifty - j, shiftz - k);
				Bzxy(nx - i, ny - j, nz - k) = WBzxy(shiftx - i, shifty - j, shiftz - k);
				Byxz(nx - i, ny - j, nz - k) = WByxz(shiftx - i, shifty - j, shiftz - k);
				Bzyx(nx - i, ny - j, nz - k) = WBzyx(shiftx - i, shifty - j, shiftz - k);
				Bxzy(nx - i, ny - j, nz - k) = WBxzy(shiftx - i, shifty - j, shiftz - k);
			}
		}
	}

	for (long i = 1; i < nx / 2 + 1; i++) {
		for (long j = 1; j < ny / 2; j++) {
			for (long k = 1; k < nz / 2; k++) {
				Axyz(i - 1, ny - j, nz - k) = WAxyz(i - 1 + shiftx, shifty - j, shiftz - k);
				Ayzx(i - 1, ny - j, nz - k) = WAyzx(i - 1 + shiftx, shifty - j, shiftz - k);
				Azxy(i - 1, ny - j, nz - k) = WAzxy(i - 1 + shiftx, shifty - j, shiftz - k);

				Bxyz(i - 1, ny - j, nz - k) = WBxyz(i - 1 + shiftx, shifty - j, shiftz - k);
				Byzx(i - 1, ny - j, nz - k) = WByzx(i - 1 + shiftx, shifty - j, shiftz - k);
				Bzxy(i - 1, ny - j, nz - k) = WBzxy(i - 1 + shiftx, shifty - j, shiftz - k);
				Byxz(i - 1, ny - j, nz - k) = WByxz(i - 1 + shiftx, shifty - j, shiftz - k);
				Bzyx(i - 1, ny - j, nz - k) = WBzyx(i - 1 + shiftx, shifty - j, shiftz - k);
				Bxzy(i - 1, ny - j, nz - k) = WBxzy(i - 1 + shiftx, shifty - j, shiftz - k);
			}
		}
	}

	for (long i = 0; i < nx; i++) {
		for (long j = 0; j < ny; j++) {
			for (long k = 0; k < nz; k++) {
				Axyz(i, j, k) = Axyz(i, j, k) / static_cast<double>(n);
				Ayzx(i, j, k) = Ayzx(i, j, k) / static_cast<double>(n);
				Azxy(i, j, k) = Azxy(i, j, k) / static_cast<double>(n);
				Bxyz(i, j, k) = Bxyz(i, j, k) / static_cast<double>(n);
				Byzx(i, j, k) = Byzx(i, j, k) / static_cast<double>(n);
				Bzxy(i, j, k) = Bzxy(i, j, k) / static_cast<double>(n);
				Byxz(i, j, k) = Byxz(i, j, k) / static_cast<double>(n);
				Bzyx(i, j, k) = Bzyx(i, j, k) / static_cast<double>(n);
				Bxzy(i, j, k) = Bxzy(i, j, k) / static_cast<double>(n);
			}
		}
	}

	//pt_math->fourier_transform3D(Axyz, Axyzk);
	//pt_math->fourier_transform3D(Ayzx, Ayzxk);
	//pt_math->fourier_transform3D(Azxy, Azxyk);
	//pt_math->fourier_transform3D(Bxyz, Bxyzk);
	//pt_math->fourier_transform3D(Byzx, Byzxk);
	//pt_math->fourier_transform3D(Bzxy, Bzxyk);
	//pt_math->fourier_transform3D(Byxz, Byxzk);
	//pt_math->fourier_transform3D(Bzyx, Bzyxk);
	//pt_math->fourier_transform3D(Bxzy, Bxzyk);

	//Axyzk *= (1. / static_cast<double>(nx * ny * nz));
	//Ayzxk *= (1. / static_cast<double>(nx * ny * nz));
	//Azxyk *= (1. / static_cast<double>(nx * ny * nz));

	//Bxyzk *= (1. / static_cast<double>(nx * ny * nz));
	//Byzxk *= (1. / static_cast<double>(nx * ny * nz));
	//Bzxyk *= (1. / static_cast<double>(nx * ny * nz));
	//Byxzk *= (1. / static_cast<double>(nx * ny * nz));
	//Bzyxk *= (1. / static_cast<double>(nx * ny * nz));
	//Bxzyk *= (1. / static_cast<double>(nx * ny * nz));

	WAxyz.nullify(); WAyzx.nullify(); WAzxy.nullify();
	WBxyz.nullify(); WByzx.nullify(); WBzxy.nullify();
	WByxz.nullify(); WBzyx.nullify(); WBxzy.nullify();;
	//Axyz.nullify(); Ayzx.nullify(); Azxy.nullify();
	//Bxyz.nullify(); Byzx.nullify(); Bzxy.nullify();
	//Byxz.nullify(); Bzyx.nullify(); Bxzy.nullify();;
}

double magnetic_system::A(double& x, double& y, double& z, double& dx, double& dy, double& dz) {
	double a;
	double xp = x + dx;
	double xm = x - dx;
	double zp = z + dz;
	double zm = z - dz;

	a = 4. * F1(x, y, z, dx, dy, dz) - 2. * F1(x, y, zp, dx, dy, dz) \
		- 2. * F1(x, y, zm, dx, dy, dz) - 2. * F1(xp, y, z, dx, dy, dz) \
		+ F1(xp, y, zp, dx, dy, dz) \
		+ F1(xp, y, zm, dx, dy, dz) \
		- 2. * F1(xm, y, z, dx, dy, dz) + F1(xm, y, zp, dx, dy, dz) + F1(xm, y, zm, dx, dy, dz);
	return a;
}

double magnetic_system::F1(double& x, double& y, double& z, double& dx, double& dy, double& dz) {
	double f1;
	double ym = y - dy;
	double yp = y + dy;
	f1 = 2. * F2(x, y, z) - F2(x, ym, z) - F2(x, yp, z);
	return f1;
}

double magnetic_system::F2(double& x, double& y, double& z) {
	double f2 = 0.;
	double r;
	double u, v, w;
	double sigma;

	sigma = 1.0e-33;
	u = x * x;
	v = y * y;
	w = z * z;
	r = sqrt(u + v + w);


	if (abs(x) > sigma && abs(y) > sigma && abs(z) > sigma) {
		f2 = z * (u - v) * log(r - z) / 2. \
			+ y * (u - w) * log(r - y) / 2. \
			- x * y * z * atan(y * z / (x * r)) \
			+ (2. * u - v - w) * r / 6.;
	}

	else if (abs(x) < sigma && abs(y) > sigma && abs(z) > sigma) {
		f2 = -(z * v / 2.) * log(r - z) - (y * w / 2.) * log(r - y) - r * r * r / 6.;
	}

	else if (abs(x) > sigma && abs(y) < sigma && abs(z) > sigma) {
		f2 = (z * u / 2.) * log(r - z) + (2. * u - w) * r / 6.;
	}

	else if (abs(x) > sigma && abs(y) > sigma && abs(z) < sigma) {
		f2 = (y * u / 2.) * log(r - y) + (2. * u - v) * r / 6.;
	}

	else if (abs(x) < sigma && abs(y) < sigma && abs(z) > sigma) {
		f2 = -abs(z) * abs(z) * abs(z) / 6.;
	}

	else if (abs(x) < sigma && abs(y) > sigma && abs(z) < sigma) {
		f2 = -abs(y) * abs(y) * abs(y) / 6.;
	}

	else if (abs(x) > sigma && abs(y) < sigma && abs(z) < sigma) {
		f2 = abs(x) * abs(x) * abs(x) / 3.;
	}

	else if (abs(x) < sigma && abs(y) < sigma && abs(z) < sigma) {
		f2 = 0.;
	}

	return f2;
}

double magnetic_system::B(double& x, double& y, double& z, double& dx, double& dy, double& dz) {
	double b;
	double yp = y + dy / 2.;
	double ym = y - dy / 2.;
	b = G1(x, yp, z, dx, dy, dz) - G1(x, ym, z, dx, dy, dz);
	return b;
}

double magnetic_system::G1(double& x, double& y, double& z, double& dx, double& dy, double& dz) {
	double g1;
	double xp = x + dx;
	double xm = x - dx;
	double yp = y + dy / 2.;
	double ym = y - dy / 2.;
	double zp = z + dz;
	double zm = z - dz;

	g1 = 4. * G2(x, yp, z) - 4. * G2(x, ym, z) \
		- 2. * G2(x, yp, zp) + 2. * G2(x, ym, zp) \
		- 2. * G2(x, yp, zm) + 2. * G2(x, ym, zm) \
		- 2. * G2(xm, yp, z) + 2. * G2(xm, ym, z) \
		+ G2(xm, yp, zp) - G2(xm, ym, zp) \
		+ G2(xm, yp, zm) - G2(xm, ym, zm) \
		- 2. * G2(xp, yp, z) + 2. * G2(xp, ym, z) \
		+ G2(xp, yp, zp) - G2(xp, ym, zp) \
		+ G2(xp, yp, zm) - G2(xp, ym, zm);
	return g1;
}

double magnetic_system::G2(double& x, double& y, double& z) {
	double g2;
	double r;
	double u, v, w;
	double sigma;

	sigma = 1.0e-33;

	u = x * x;
	v = y * y;
	w = z * z;
	r = sqrt(u + v + w);

	if (abs(x * y * z) > sigma) {
		g2 = (y * y * y / 6. - y * w / 2.) * log(r + x)\
			+ (x / 6.) * (3. * w - u) * log(r - y)\
			+ x * y * r / 3. + x * y * z * log(r - z)\
			+ (v * z / 2.) * atan(x * z / (r * y))\
			+ (z * z * z / 6.) * atan(x * y / (r * z))\
			+ (u * z / 2.) * atan(y * z / (r * x));
	}

	else if (abs(x) < sigma && abs(y * z) > sigma) {
		g2 = (y * y * y / 6. - y * w / 2.) * log(r);
	}

	else if (abs(x * z) > sigma && abs(y) < sigma) {
		g2 = (x * w / 2. - x * x * x / 6.) * log(r);
	}

	else if (abs(x * y) > sigma && abs(z) < sigma) {
		g2 = ((y * y * y) / 6.) * log(r + x) - (x * x * x / 6.) * log(r - y) + x * y * r / 3.;
	}

	else if ((abs(x) + abs(y)) < sigma && abs(z) > sigma) {
		g2 = 0.;
	}

	else if ((abs(x) + abs(z)) < sigma && abs(y) > sigma) {
		g2 = (y * y * y) / 12. * log(v);
	}

	else if (abs(x) > sigma && (abs(y) + abs(z)) < sigma) {
		g2 = -((x * x * x) / 12.) * log(u);
	}

	else if (abs(x) < sigma && abs(y) < sigma && abs(z) < sigma) {
		g2 = 0.;
	}
	return g2;
}

void magnetic_system::copy_to_device() {
	long int n21 = nx * ny * (nz / 2 + 1);

#pragma acc enter data copyin(this)

#pragma acc serial default(present)
	{
		this->pt_geo = &(geometry_parameters::geo);
		this->pt_glb = &(global_parameters::glb);
		this->pt_math = &(mathlib::mlb);
	}
	//present(geometry_parameters::geo,global_parameters::glb,mathlib::mlb) deviceptr(this->pt_geo,this->pt_glb,this->pt_math)

	if (pt_glb->if_mag_1Dmodel == false && pt_glb->if_magnetostatics == true && pt_glb->if_demag_factor == false) {
#pragma acc enter data copyin(this->Axyz)
#pragma acc enter data copyin(this->Bxyz)
#pragma acc enter data copyin(this->Bxzy)
#pragma acc enter data copyin(this->Ayzx)
#pragma acc enter data copyin(this->Byzx)
#pragma acc enter data copyin(this->Byxz)
#pragma acc enter data copyin(this->Azxy)
#pragma acc enter data copyin(this->Bzxy)
#pragma acc enter data copyin(this->Bzyx)

#pragma acc enter data copyin(this->Axyz.matrix[0:n])
#pragma acc enter data copyin(this->Bxyz.matrix[0:n])
#pragma acc enter data copyin(this->Bxzy.matrix[0:n])
#pragma acc enter data copyin(this->Ayzx.matrix[0:n])
#pragma acc enter data copyin(this->Byzx.matrix[0:n])
#pragma acc enter data copyin(this->Byxz.matrix[0:n])
#pragma acc enter data copyin(this->Azxy.matrix[0:n])
#pragma acc enter data copyin(this->Bzxy.matrix[0:n])
#pragma acc enter data copyin(this->Bzyx.matrix[0:n])

#pragma acc enter data copyin(this->Axyzk)
#pragma acc enter data copyin(this->Bxyzk)
#pragma acc enter data copyin(this->Bxzyk)
#pragma acc enter data copyin(this->Ayzxk)
#pragma acc enter data copyin(this->Byzxk)
#pragma acc enter data copyin(this->Byxzk)
#pragma acc enter data copyin(this->Azxyk)
#pragma acc enter data copyin(this->Bzxyk)
#pragma acc enter data copyin(this->Bzyxk)

#pragma acc enter data copyin(this->Axyzk.matrix[0:n21])
#pragma acc enter data copyin(this->Bxyzk.matrix[0:n21])
#pragma acc enter data copyin(this->Bxzyk.matrix[0:n21])
#pragma acc enter data copyin(this->Ayzxk.matrix[0:n21])
#pragma acc enter data copyin(this->Byzxk.matrix[0:n21])
#pragma acc enter data copyin(this->Byxzk.matrix[0:n21])
#pragma acc enter data copyin(this->Azxyk.matrix[0:n21])
#pragma acc enter data copyin(this->Bzxyk.matrix[0:n21])
#pragma acc enter data copyin(this->Bzyxk.matrix[0:n21])

#pragma acc enter data copyin(this->Mx_glb,this->My_glb,this->Mz_glb)
#pragma acc enter data copyin(this->Mx_k3D,this->My_k3D,this->Mz_k3D)	
#pragma acc enter data copyin(this->Hx_stat_k,this->Hy_stat_k,this->Hz_stat_k)	

#pragma acc enter data copyin(this->Mx_glb.matrix[0:n],this->My_glb.matrix[0:n],this->Mz_glb.matrix[0:n])
#pragma acc enter data copyin(this->Mx_k3D.matrix[0:n21],this->My_k3D.matrix[0:n21],this->Mz_k3D.matrix[0:n21])	
#pragma acc enter data copyin(this->Hx_stat_k.matrix[0:n21],this->Hy_stat_k.matrix[0:n21],this->Hz_stat_k.matrix[0:n21])	
	}

	if (pt_glb->if_FM_all == true) {
#pragma acc enter data copyin(this->FM_surfYZ)
#pragma acc enter data copyin(this->FM_surfXZ)
#pragma acc enter data copyin(this->FM_surfXY)

#pragma acc enter data copyin(this->FM_surfYZ.matrix[0:(nx + 1)*ny*nz])
#pragma acc enter data copyin(this->FM_surfXZ.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->FM_surfXY.matrix[0:nx*ny*(nz+1)])

		//----------Averaged parameter for inter-region exchange - Approach 1---------//
#pragma acc enter data copyin(this->Aij_xf)
#pragma acc enter data copyin(this->Aij_yf)
#pragma acc enter data copyin(this->Aij_zf)
#pragma acc enter data copyin(this->Aij_xb)
#pragma acc enter data copyin(this->Aij_yb)
#pragma acc enter data copyin(this->Aij_zb)

#pragma acc enter data copyin(this->Aij_xf.matrix[0:n])
#pragma acc enter data copyin(this->Aij_yf.matrix[0:n])
#pragma acc enter data copyin(this->Aij_zf.matrix[0:n])
#pragma acc enter data copyin(this->Aij_xb.matrix[0:n])
#pragma acc enter data copyin(this->Aij_yb.matrix[0:n])
#pragma acc enter data copyin(this->Aij_zb.matrix[0:n])

//----------Averaged parameter for inter-region exchange - Approach 2---------//
//#pragma acc enter data copyin(this->AM_exchange)
//#pragma acc enter data copyin(this->AM_exchange.matrix[0:n])
	}

	if (pt_glb->if_AFM_all == true) {
#pragma acc enter data copyin(this->AFM_surfYZ)
#pragma acc enter data copyin(this->AFM_surfXZ)
#pragma acc enter data copyin(this->AFM_surfXY)

#pragma acc enter data copyin(this->AFM_surfYZ.matrix[0:(nx + 1)*ny*nz])
#pragma acc enter data copyin(this->AFM_surfXZ.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->AFM_surfXY.matrix[0:nx*ny*(nz+1)])

		//----------Averaged parameter for inter-region exchange - Approach 1---------//
#pragma acc enter data copyin(this->Aij_AFM1_xf,this->Aij_AFM2_xf)
#pragma acc enter data copyin(this->Aij_AFM1_yf,this->Aij_AFM2_yf)
#pragma acc enter data copyin(this->Aij_AFM1_zf,this->Aij_AFM2_zf)
#pragma acc enter data copyin(this->Aij_AFM1_xb,this->Aij_AFM2_xb)
#pragma acc enter data copyin(this->Aij_AFM1_yb,this->Aij_AFM2_yb)
#pragma acc enter data copyin(this->Aij_AFM1_zb,this->Aij_AFM2_zb)

#pragma acc enter data copyin(this->Aij_AFM1_xf.matrix[0:n],this->Aij_AFM2_xf.matrix[0:n])
#pragma acc enter data copyin(this->Aij_AFM1_yf.matrix[0:n],this->Aij_AFM2_yf.matrix[0:n])
#pragma acc enter data copyin(this->Aij_AFM1_zf.matrix[0:n],this->Aij_AFM2_zf.matrix[0:n])
#pragma acc enter data copyin(this->Aij_AFM1_xb.matrix[0:n],this->Aij_AFM2_xb.matrix[0:n])
#pragma acc enter data copyin(this->Aij_AFM1_yb.matrix[0:n],this->Aij_AFM2_yb.matrix[0:n])
#pragma acc enter data copyin(this->Aij_AFM1_zb.matrix[0:n],this->Aij_AFM2_zb.matrix[0:n])
	}
	//
//------------------FM------------------//
#pragma acc enter data copyin(this->mx_glb)
#pragma acc enter data copyin(this->my_glb)
#pragma acc enter data copyin(this->mz_glb)

#pragma acc enter data copyin(this->mx_glb.matrix[0:n])
#pragma acc enter data copyin(this->my_glb.matrix[0:n])
#pragma acc enter data copyin(this->mz_glb.matrix[0:n])

#pragma acc enter data copyin(this->dmx_glb)
#pragma acc enter data copyin(this->dmy_glb)
#pragma acc enter data copyin(this->dmz_glb)

#pragma acc enter data copyin(this->dmx_glb.matrix[0:n])
#pragma acc enter data copyin(this->dmy_glb.matrix[0:n])
#pragma acc enter data copyin(this->dmz_glb.matrix[0:n])

#pragma acc enter data copyin(this->mx_glb_store, this->my_glb_store, this->mz_glb_store)
#pragma acc enter data copyin(this->mx_glb_store.matrix[0:n], this->my_glb_store.matrix[0:n], this->mz_glb_store.matrix[0:n])

#pragma acc enter data copyin(this->dmx_glb_rk1, this->dmy_glb_rk1, this->dmz_glb_rk1)
#pragma acc enter data copyin(this->dmx_glb_rk2, this->dmy_glb_rk2, this->dmz_glb_rk2)
#pragma acc enter data copyin(this->dmx_glb_rk3, this->dmy_glb_rk3, this->dmz_glb_rk3)
#pragma acc enter data copyin(this->dmx_glb_rk4, this->dmy_glb_rk4, this->dmz_glb_rk4)

#pragma acc enter data copyin(this->dmx_glb_rk1.matrix[0:n], this->dmy_glb_rk1.matrix[0:n], this->dmz_glb_rk1.matrix[0:n])
#pragma acc enter data copyin(this->dmx_glb_rk2.matrix[0:n], this->dmy_glb_rk2.matrix[0:n], this->dmz_glb_rk2.matrix[0:n])
#pragma acc enter data copyin(this->dmx_glb_rk3.matrix[0:n], this->dmy_glb_rk3.matrix[0:n], this->dmz_glb_rk3.matrix[0:n])
#pragma acc enter data copyin(this->dmx_glb_rk4.matrix[0:n], this->dmy_glb_rk4.matrix[0:n], this->dmz_glb_rk4.matrix[0:n])

	//------------------AFM------------------//
#pragma acc enter data copyin(this->mx_AFM1_glb,this->mx_AFM2_glb)
#pragma acc enter data copyin(this->my_AFM1_glb,this->my_AFM2_glb)
#pragma acc enter data copyin(this->mz_AFM1_glb,this->mz_AFM2_glb)

#pragma acc enter data copyin(this->mx_AFM1_glb.matrix[0:n],this->mx_AFM2_glb.matrix[0:n])
#pragma acc enter data copyin(this->my_AFM1_glb.matrix[0:n],this->my_AFM2_glb.matrix[0:n])
#pragma acc enter data copyin(this->mz_AFM1_glb.matrix[0:n],this->mz_AFM2_glb.matrix[0:n])

#pragma acc enter data copyin(this->dmx_AFM1_glb,this->dmx_AFM2_glb)
#pragma acc enter data copyin(this->dmy_AFM1_glb,this->dmy_AFM2_glb)
#pragma acc enter data copyin(this->dmz_AFM1_glb,this->dmz_AFM2_glb)

#pragma acc enter data copyin(this->dmx_AFM1_glb.matrix[0:n],this->dmx_AFM2_glb.matrix[0:n])
#pragma acc enter data copyin(this->dmy_AFM1_glb.matrix[0:n],this->dmy_AFM2_glb.matrix[0:n])
#pragma acc enter data copyin(this->dmz_AFM1_glb.matrix[0:n],this->dmz_AFM2_glb.matrix[0:n])

#pragma acc enter data copyin(this->mx_AFM1_glb_store, this->my_AFM1_glb_store, this->mz_AFM1_glb_store)
#pragma acc enter data copyin(this->mx_AFM2_glb_store, this->my_AFM2_glb_store, this->mz_AFM2_glb_store)

#pragma acc enter data copyin(this->mx_AFM1_glb_store.matrix[0:n], this->my_AFM1_glb_store.matrix[0:n], this->mz_AFM1_glb_store.matrix[0:n])
#pragma acc enter data copyin(this->mx_AFM2_glb_store.matrix[0:n], this->my_AFM2_glb_store.matrix[0:n], this->mz_AFM2_glb_store.matrix[0:n])

#pragma acc enter data copyin(this->dmx_AFM1_glb_rk1, this->dmy_AFM1_glb_rk1, this->dmz_AFM1_glb_rk1)
#pragma acc enter data copyin(this->dmx_AFM1_glb_rk2, this->dmy_AFM1_glb_rk2, this->dmz_AFM1_glb_rk2)
#pragma acc enter data copyin(this->dmx_AFM1_glb_rk3, this->dmy_AFM1_glb_rk3, this->dmz_AFM1_glb_rk3)
#pragma acc enter data copyin(this->dmx_AFM1_glb_rk4, this->dmy_AFM1_glb_rk4, this->dmz_AFM1_glb_rk4)

#pragma acc enter data copyin(this->dmx_AFM1_glb_rk1.matrix[0:n], this->dmy_AFM1_glb_rk1.matrix[0:n], this->dmz_AFM1_glb_rk1.matrix[0:n])
#pragma acc enter data copyin(this->dmx_AFM1_glb_rk2.matrix[0:n], this->dmy_AFM1_glb_rk2.matrix[0:n], this->dmz_AFM1_glb_rk2.matrix[0:n])
#pragma acc enter data copyin(this->dmx_AFM1_glb_rk3.matrix[0:n], this->dmy_AFM1_glb_rk3.matrix[0:n], this->dmz_AFM1_glb_rk3.matrix[0:n])
#pragma acc enter data copyin(this->dmx_AFM1_glb_rk4.matrix[0:n], this->dmy_AFM1_glb_rk4.matrix[0:n], this->dmz_AFM1_glb_rk4.matrix[0:n])

#pragma acc enter data copyin(this->dmx_AFM2_glb_rk1, this->dmy_AFM2_glb_rk1, this->dmz_AFM2_glb_rk1)
#pragma acc enter data copyin(this->dmx_AFM2_glb_rk2, this->dmy_AFM2_glb_rk2, this->dmz_AFM2_glb_rk2)
#pragma acc enter data copyin(this->dmx_AFM2_glb_rk3, this->dmy_AFM2_glb_rk3, this->dmz_AFM2_glb_rk3)
#pragma acc enter data copyin(this->dmx_AFM2_glb_rk4, this->dmy_AFM2_glb_rk4, this->dmz_AFM2_glb_rk4)

#pragma acc enter data copyin(this->dmx_AFM2_glb_rk1.matrix[0:n], this->dmy_AFM2_glb_rk1.matrix[0:n], this->dmz_AFM2_glb_rk1.matrix[0:n])
#pragma acc enter data copyin(this->dmx_AFM2_glb_rk2.matrix[0:n], this->dmy_AFM2_glb_rk2.matrix[0:n], this->dmz_AFM2_glb_rk2.matrix[0:n])
#pragma acc enter data copyin(this->dmx_AFM2_glb_rk3.matrix[0:n], this->dmy_AFM2_glb_rk3.matrix[0:n], this->dmz_AFM2_glb_rk3.matrix[0:n])
#pragma acc enter data copyin(this->dmx_AFM2_glb_rk4.matrix[0:n], this->dmy_AFM2_glb_rk4.matrix[0:n], this->dmz_AFM2_glb_rk4.matrix[0:n])
	//----------------------------------------------------//

#pragma acc enter data copyin(this->Hx_stat)
#pragma acc enter data copyin(this->Hy_stat)
#pragma acc enter data copyin(this->Hz_stat)

#pragma acc enter data copyin(this->Hx_stat.matrix[0:n])
#pragma acc enter data copyin(this->Hy_stat.matrix[0:n])
#pragma acc enter data copyin(this->Hz_stat.matrix[0:n])

#pragma acc enter data copyin(this->Jx_ISHE)
#pragma acc enter data copyin(this->Jy_ISHE)

#pragma acc enter data copyin(this->Jx_ISHE.matrix[0:n])
#pragma acc enter data copyin(this->Jy_ISHE.matrix[0:n])

	//#pragma acc enter data copyin(this->Hx_anis)
	//#pragma acc enter data copyin(this->Hy_anis)
	//#pragma acc enter data copyin(this->Hz_anis)
	//
	//#pragma acc enter data copyin(this->Hx_anis.matrix[0:n])
	//#pragma acc enter data copyin(this->Hy_anis.matrix[0:n])
	//#pragma acc enter data copyin(this->Hz_anis.matrix[0:n])
	//
	//#pragma acc enter data copyin(this->Hx_exch)
	//#pragma acc enter data copyin(this->Hy_exch)
	//#pragma acc enter data copyin(this->Hz_exch)
	//
	//#pragma acc enter data copyin(this->Hx_exch.matrix[0:n])
	//#pragma acc enter data copyin(this->Hy_exch.matrix[0:n])
	//#pragma acc enter data copyin(this->Hz_exch.matrix[0:n])
	//
	//#pragma acc enter data copyin(this->Hx_elas)
	//#pragma acc enter data copyin(this->Hy_elas)
	//#pragma acc enter data copyin(this->Hz_elas)
	//
	//#pragma acc enter data copyin(this->Hx_elas.matrix[0:n])
	//#pragma acc enter data copyin(this->Hy_elas.matrix[0:n])
	//#pragma acc enter data copyin(this->Hz_elas.matrix[0:n])
}

void magnetic_system::copym_from_device() {
	//#pragma acc update self(mx_glb, my_glb, mz_glb)
	//	;
#pragma acc update self(mx_glb.matrix[0:n], my_glb.matrix[0:n], mz_glb.matrix[0:n])
	;
}

void magnetic_system::copym_AFM_from_device() {
	//#pragma acc update self(mx_glb, my_glb, mz_glb)
	//	;
#pragma acc update self(mx_AFM1_glb.matrix[0:n], my_AFM1_glb.matrix[0:n], mz_AFM1_glb.matrix[0:n])
#pragma acc update self(mx_AFM2_glb.matrix[0:n], my_AFM2_glb.matrix[0:n], mz_AFM2_glb.matrix[0:n])
	;
}

void magnetic_system::copyH_from_device() {
	//#pragma acc update self(Hx_stat, Hy_stat, Hz_stat)
	//	;
#pragma acc update self(Hx_stat.matrix[0:n], Hy_stat.matrix[0:n], Hz_stat.matrix[0:n])
	;
}

void magnetic_system::copyJishe_from_device() {
#pragma acc update self(Jx_ISHE.matrix[0:n],Jy_ISHE.matrix[0:n])
	;
}
