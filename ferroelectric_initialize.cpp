#include "ferroelectric_system.h"
#include<iostream>
void ferroelectric_system::initialize_host() {
	unsigned int material_type1, material_type2;
	bool if_F1, if_F2;
	//double kq[3];
	double temporary_var;

	//pt_elas=pt_elas_in; 
	//pt_EM=pt_EM_in;

	nx = pt_geo->nx_system;
	ny = pt_geo->ny_system;
	nz = pt_geo->nz_system;

	n = nx * ny * nz;
	nz21 = nz / 2 + 1;
	scalenn = 1. / static_cast<double>(n);

	dx = pt_geo->dx;
	dy = pt_geo->dy;
	dz = pt_geo->dz;

	px_glb.initialize(nx, ny, nz);  py_glb.initialize(nx, ny, nz);  pz_glb.initialize(nx, ny, nz);
	if (pt_glb->if_elec_1D_compensate == true) {
		pz_t0_glb.initialize(nx, ny, nz);
	}
	//dpx_glb.initialize(nx, ny, nz); dpy_glb.initialize(nx, ny, nz); dpz_glb.initialize(nx, ny, nz);
	qx_glb.initialize(nx, ny, nz);  qy_glb.initialize(nx, ny, nz);  qz_glb.initialize(nx, ny, nz);
	//dqx_glb.initialize(nx, ny, nz); dqy_glb.initialize(nx, ny, nz); dqz_glb.initialize(nx, ny, nz);
	// 
	//RK4
	px_glb_store.initialize(nx, ny, nz); py_glb_store.initialize(nx, ny, nz); pz_glb_store.initialize(nx, ny, nz);
	dpx_glb_rk1.initialize(nx, ny, nz); dpy_glb_rk1.initialize(nx, ny, nz); dpz_glb_rk1.initialize(nx, ny, nz);
	dpx_glb_rk2.initialize(nx, ny, nz); dpy_glb_rk2.initialize(nx, ny, nz); dpz_glb_rk2.initialize(nx, ny, nz);
	dpx_glb_rk3.initialize(nx, ny, nz); dpy_glb_rk3.initialize(nx, ny, nz); dpz_glb_rk3.initialize(nx, ny, nz);
	dpx_glb_rk4.initialize(nx, ny, nz); dpy_glb_rk4.initialize(nx, ny, nz); dpz_glb_rk4.initialize(nx, ny, nz);

	qx_glb_store.initialize(nx, ny, nz); qy_glb_store.initialize(nx, ny, nz); qz_glb_store.initialize(nx, ny, nz);
	dqx_glb_rk1.initialize(nx, ny, nz); dqy_glb_rk1.initialize(nx, ny, nz); dqz_glb_rk1.initialize(nx, ny, nz);
	dqx_glb_rk2.initialize(nx, ny, nz); dqy_glb_rk2.initialize(nx, ny, nz); dqz_glb_rk2.initialize(nx, ny, nz);
	dqx_glb_rk3.initialize(nx, ny, nz); dqy_glb_rk3.initialize(nx, ny, nz); dqz_glb_rk3.initialize(nx, ny, nz);
	dqx_glb_rk4.initialize(nx, ny, nz); dqy_glb_rk4.initialize(nx, ny, nz); dqz_glb_rk4.initialize(nx, ny, nz);
	//----------------------//

	if (pt_glb->if_FE_all == true && pt_glb->if_flexo == true) {
		dpxdx.initialize(nx, ny, nz); dpxdy.initialize(nx, ny, nz); dpxdz.initialize(nx, ny, nz);
		dpydx.initialize(nx, ny, nz); dpydy.initialize(nx, ny, nz); dpydz.initialize(nx, ny, nz);
		dpzdx.initialize(nx, ny, nz); dpzdy.initialize(nx, ny, nz); dpzdz.initialize(nx, ny, nz);
	}

	if (pt_glb->if_FE_all == true && pt_glb->if_input_pandq == false) {
		initialize_polarization();
	}

	//Electrostatics
	//Estat0_1D = pt_glb->total_free_charge / -2. / e0;

	Ex_stat.initialize(nx, ny, nz); Ey_stat.initialize(nx, ny, nz); Ez_stat.initialize(nx, ny, nz);

	rhof.initialize(nx, ny, nz);

	rhofk.initialize(nx, ny, nz / 2 + 1);
	Dxk_inhomo.initialize(nx, ny, nz / 2 + 1); Dyk_inhomo.initialize(nx, ny, nz / 2 + 1); Dzk_inhomo.initialize(nx, ny, nz / 2 + 1);
	pxk.initialize(nx, ny, nz / 2 + 1); pyk.initialize(nx, ny, nz / 2 + 1); pzk.initialize(nx, ny, nz / 2 + 1);
	ps_terms.initialize(nx, ny, nz / 2 + 1);
	Edxk_n.initialize(nx, ny, nz / 2 + 1); Edyk_n.initialize(nx, ny, nz / 2 + 1); Edzk_n.initialize(nx, ny, nz / 2 + 1);

	Dx_inhomo.initialize(nx, ny, nz); Dy_inhomo.initialize(nx, ny, nz); Dz_inhomo.initialize(nx, ny, nz);

	//------------------------------------------//
	//	Determine if FE/Non-FE boundaries		//
	//------------------------------------------//
	FE_surfYZ.initialize(nx + 1, ny, nz);
	FE_surfXZ.initialize(nx, ny + 1, nz);
	FE_surfXY.initialize(nx, ny, nz + 1);

	if (pt_glb->if_FE_all == true) {
		//------------X Direction (YZ plane)----------------//
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				material_type1 = pt_glb->material_cell(0, j, k);
				if (material_type1 == 0) {
					if_F1 = false;
				}
				else {
					if_F1 = pt_glb->material_parameters[static_cast<long int>(material_type1) - 1].if_FE;
				}

				material_type2 = pt_glb->material_cell(nx - 1, j, k);
				if (material_type2 == 0) {
					if_F2 = false;
				}
				else {
					if_F2 = pt_glb->material_parameters[static_cast<long int>(material_type2) - 1].if_FE;
				}

				if (pt_geo->periodicX == true) {
					if (if_F1 ^ if_F2) {
						FE_surfYZ(0, j, k) = true;
						FE_surfYZ(nx, j, k) = true;
					}
					else {
						FE_surfYZ(0, j, k) = false;
						FE_surfYZ(nx, j, k) = false;
					}
				}
				else if (pt_geo->periodicX == false) {
					if (if_F1 == true) {
						FE_surfYZ(0, j, k) = true;
					}
					else {
						FE_surfYZ(0, j, k) = false;
					}

					if (if_F2 == true) {
						FE_surfYZ(nx, j, k) = true;
					}
					else {
						FE_surfYZ(nx, j, k) = false;
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
						if_F1 = pt_glb->material_parameters[static_cast<long int>(material_type1) - 1].if_FE;
					}

					material_type2 = pt_glb->material_cell(i, j, k);
					if (material_type2 == 0) {
						if_F2 = false;
					}
					else {
						if_F2 = pt_glb->material_parameters[static_cast<long int>(material_type2) - 1].if_FE;
					}

					if (if_F1 ^ if_F2) {
						FE_surfYZ(i, j, k) = true;
					}
					else {
						FE_surfYZ(i, j, k) = false;
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
					if_F1 = pt_glb->material_parameters[static_cast<long int>(material_type1) - 1].if_FE;
				}

				material_type2 = pt_glb->material_cell(i, ny - 1, k);
				if (material_type2 == 0) {
					if_F2 = false;
				}
				else {
					if_F2 = pt_glb->material_parameters[static_cast<long int>(material_type2) - 1].if_FE;
				}

				if (pt_geo->periodicY == true) {
					if (if_F1 ^ if_F2) {
						FE_surfXZ(i, 0, k) = true;
						FE_surfXZ(i, ny, k) = true;
					}
					else {
						FE_surfXZ(i, 0, k) = false;
						FE_surfXZ(i, ny, k) = false;
					}
				}
				else if (pt_geo->periodicY == false) {
					if (if_F1 == true) {
						FE_surfXZ(i, 0, k) = true;
					}
					else {
						FE_surfXZ(i, 0, k) = false;
					}

					if (if_F2 == true) {
						FE_surfXZ(i, ny, k) = true;
					}
					else {
						FE_surfXZ(i, ny, k) = false;
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
						if_F1 = pt_glb->material_parameters[static_cast<long int>(material_type1) - 1].if_FE;
					}

					material_type2 = pt_glb->material_cell(i, j, k);
					if (material_type2 == 0) {
						if_F2 = false;
					}
					else {
						if_F2 = pt_glb->material_parameters[static_cast<long int>(material_type2) - 1].if_FE;
					}

					if (if_F1 ^ if_F2) {
						FE_surfXZ(i, j, k) = true;
					}
					else {
						FE_surfXZ(i, j, k) = false;
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
					if_F1 = pt_glb->material_parameters[static_cast<long int>(material_type1) - 1].if_FE;
				}

				material_type2 = pt_glb->material_cell(i, j, nz - 1);
				if (material_type2 == 0) {
					if_F2 = false;
				}
				else {
					if_F2 = pt_glb->material_parameters[static_cast<long int>(material_type2) - 1].if_FE;
				}

				if (pt_geo->periodicZ == true) {
					if (if_F1 ^ if_F2) {
						FE_surfXY(i, j, 0) = true;
						FE_surfXY(i, j, nz) = true;
					}
					else {
						FE_surfXY(i, j, 0) = false;
						FE_surfXY(i, j, nz) = false;
					}
				}
				else if (pt_geo->periodicZ == false) {
					if (if_F1 == true) {
						FE_surfXY(i, j, 0) = true;
					}
					else {
						FE_surfXY(i, j, 0) = false;
					}

					if (if_F2 == true) {
						FE_surfXY(i, j, nz) = true;
					}
					else {
						FE_surfXY(i, j, nz) = false;
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
						if_F1 = pt_glb->material_parameters[static_cast<long int>(material_type1) - 1].if_FE;
					}

					material_type2 = pt_glb->material_cell(i, j, k);
					if (material_type2 == 0) {
						if_F2 = false;
					}
					else {
						if_F2 = pt_glb->material_parameters[static_cast<long int>(material_type2) - 1].if_FE;
					}

					if (if_F1 ^ if_F2) {
						FE_surfXY(i, j, k) = true;
					}
					else {
						FE_surfXY(i, j, k) = false;
					}
				}
			}
		}
	}

	//	Calculate homogeneous electrical permittivity for static solver//
	er_homo = pt_glb->material_parameters[0].r_permittivity;

	for (long int i = 1; i < pt_glb->num_materials; i++) {
		if (pt_glb->material_parameters[i].r_permittivity > er_homo) {
			er_homo = pt_glb->material_parameters[i].r_permittivity;
		}
	}
	er_homo = er_homo * pt_glb->scale_elec;

	//--------Initialize intermediate variables for electrostatic solver--------//	
	iGq_inverse.initialize(nx, ny, nz / 2 + 1);

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz / 2 + 1; k++) {
				temporary_var = er_homo * pt_glb->k_norm(i, j, k);

				if (abs(temporary_var) > 1.e-8) {
					temporary_var = -1. / temporary_var;
				}
				else
				{
					temporary_var = -1.e30;
				}
				iGq_inverse(i, j, k) = im_unit * temporary_var;
			}
		}
	}

	//for (long i = 0; i < nx; i++) {
	//	for (long j = 0; j < ny; j++) {
	//		for (long k = pt_glb->free_charge_nzi - 1; k < pt_glb->free_charge_nzf; k++) {
	//			rhof(i, j, k) = pt_glb->free_charge;
	//		}
	//	}
	//}
	e_field[0] = 0.; e_field[1] = 0.; e_field[2] = 0.; e_field[3] = 0.;
	e_sum = 0.;
	tolerance = 0.;
	//----------Initialize cuHandle for cuFFT--------------//
	cufftPlan3d(&plan_electro_D2Z, nx, ny, nz, CUFFT_D2Z);
	cufftPlan3d(&plan_electro_Z2D, nx, ny, nz, CUFFT_Z2D);
}

void ferroelectric_system::initialize_device() {
//	if (pt_glb->if_input_pandq == true) {
#pragma acc parallel default(present) async(2)
		{
			px_glb_store = px_glb;
			py_glb_store = py_glb;
			pz_glb_store = pz_glb;

			qx_glb_store = qx_glb;
			qy_glb_store = qy_glb;
			qz_glb_store = qz_glb;

			if (pt_glb->if_elec_1D_compensate == true) {
				pz_t0_glb = pz_glb;
			}
		}

//		if (pt_glb->if_flexo == true) {
//#pragma acc parallel default(present) async(2)
//			{
//				get_p_gradient();
//			}
//		}
//	}
}

void ferroelectric_system::initialize_polarization() {
	unsigned int material_type;
	double px = pt_glb->px_init;
	double py = pt_glb->py_init;
	double pz = pt_glb->pz_init;

	for (long int id = 0; id < n; id++) {
		material_type = pt_glb->material_cell(id);
		if (material_type != 0) {
			if (pt_glb->material_parameters[(material_type - 1)].if_FE == true) {
				px_glb(id) = px; py_glb(id) = py; pz_glb(id) = pz;
			}
		}
	}
}

void ferroelectric_system::copy_to_device() {
	long n21 = nx * ny * (nz / 2 + 1);

#pragma acc enter data copyin(this)

#pragma acc serial default(present)
	{
		this->pt_geo = &(geometry_parameters::geo);
		this->pt_glb = &(global_parameters::glb);
		this->pt_math = &(mathlib::mlb);
	}

#pragma acc enter data copyin(this->iGq_inverse)
#pragma acc enter data copyin(this->iGq_inverse.matrix[0:n21])

	if (pt_glb->if_FE_all == true) {
#pragma acc enter data copyin(this->FE_surfYZ)
#pragma acc enter data copyin(this->FE_surfXZ)
#pragma acc enter data copyin(this->FE_surfXY)
#pragma acc enter data copyin(this->FE_surfYZ.matrix[0:(nx + 1)*ny*nz])
#pragma acc enter data copyin(this->FE_surfXZ.matrix[0:nx*(ny + 1)*nz])
#pragma acc enter data copyin(this->FE_surfXY.matrix[0:nx*ny*(nz + 1)])
	}

#pragma acc enter data copyin(this->px_glb, this->px_glb_store)
#pragma acc enter data copyin(this->py_glb, this->py_glb_store)
#pragma acc enter data copyin(this->pz_glb, this->pz_glb_store)
#pragma acc enter data copyin(this->px_glb.matrix[0:n], this->px_glb_store.matrix[0:n])
#pragma acc enter data copyin(this->py_glb.matrix[0:n], this->py_glb_store.matrix[0:n])
#pragma acc enter data copyin(this->pz_glb.matrix[0:n], this->pz_glb_store.matrix[0:n])

	if (pt_glb->if_elec_1D_compensate == true) {
#pragma acc enter data copyin(this->pz_t0_glb)
#pragma acc enter data copyin(this->pz_t0_glb.matrix[0:n])
	}

	//#pragma acc enter data copyin(this->dpx_glb)
	//#pragma acc enter data copyin(this->dpy_glb)
	//#pragma acc enter data copyin(this->dpz_glb)
	//#pragma acc enter data copyin(this->dpx_glb.matrix[0:n])
	//#pragma acc enter data copyin(this->dpy_glb.matrix[0:n])
	//#pragma acc enter data copyin(this->dpz_glb.matrix[0:n])

#pragma acc enter data copyin(this->dpx_glb_rk1,this->dpx_glb_rk2,this->dpx_glb_rk3,this->dpx_glb_rk4)
#pragma acc enter data copyin(this->dpy_glb_rk1,this->dpy_glb_rk2,this->dpy_glb_rk3,this->dpy_glb_rk4)
#pragma acc enter data copyin(this->dpz_glb_rk1,this->dpz_glb_rk2,this->dpz_glb_rk3,this->dpz_glb_rk4)
#pragma acc enter data copyin(this->dpx_glb_rk1.matrix[0:n])
#pragma acc enter data copyin(this->dpy_glb_rk1.matrix[0:n])
#pragma acc enter data copyin(this->dpz_glb_rk1.matrix[0:n])
#pragma acc enter data copyin(this->dpx_glb_rk2.matrix[0:n])
#pragma acc enter data copyin(this->dpy_glb_rk2.matrix[0:n])
#pragma acc enter data copyin(this->dpz_glb_rk2.matrix[0:n])
#pragma acc enter data copyin(this->dpx_glb_rk3.matrix[0:n])
#pragma acc enter data copyin(this->dpy_glb_rk3.matrix[0:n])
#pragma acc enter data copyin(this->dpz_glb_rk3.matrix[0:n])
#pragma acc enter data copyin(this->dpx_glb_rk4.matrix[0:n])
#pragma acc enter data copyin(this->dpy_glb_rk4.matrix[0:n])
#pragma acc enter data copyin(this->dpz_glb_rk4.matrix[0:n])

#pragma acc enter data copyin(this->qx_glb, this->qx_glb_store)
#pragma acc enter data copyin(this->qy_glb, this->qy_glb_store)
#pragma acc enter data copyin(this->qz_glb, this->qz_glb_store)
#pragma acc enter data copyin(this->qx_glb.matrix[0:n], this->qx_glb_store.matrix[0:n])
#pragma acc enter data copyin(this->qy_glb.matrix[0:n], this->qy_glb_store.matrix[0:n])
#pragma acc enter data copyin(this->qz_glb.matrix[0:n], this->qz_glb_store.matrix[0:n])

//#pragma acc enter data copyin(this->dqx_glb)
//#pragma acc enter data copyin(this->dqy_glb)
//#pragma acc enter data copyin(this->dqz_glb)
//#pragma acc enter data copyin(this->dqx_glb.matrix[0:n])
//#pragma acc enter data copyin(this->dqy_glb.matrix[0:n])
//#pragma acc enter data copyin(this->dqz_glb.matrix[0:n])
#pragma acc enter data copyin(this->dqx_glb_rk1,this->dqx_glb_rk2,this->dqx_glb_rk3,this->dqx_glb_rk4)
#pragma acc enter data copyin(this->dqy_glb_rk1,this->dqy_glb_rk2,this->dqy_glb_rk3,this->dqy_glb_rk4)
#pragma acc enter data copyin(this->dqz_glb_rk1,this->dqz_glb_rk2,this->dqz_glb_rk3,this->dqz_glb_rk4)
#pragma acc enter data copyin(this->dqx_glb_rk1.matrix[0:n])
#pragma acc enter data copyin(this->dqy_glb_rk1.matrix[0:n])
#pragma acc enter data copyin(this->dqz_glb_rk1.matrix[0:n])
#pragma acc enter data copyin(this->dqx_glb_rk2.matrix[0:n])
#pragma acc enter data copyin(this->dqy_glb_rk2.matrix[0:n])
#pragma acc enter data copyin(this->dqz_glb_rk2.matrix[0:n])
#pragma acc enter data copyin(this->dqx_glb_rk3.matrix[0:n])
#pragma acc enter data copyin(this->dqy_glb_rk3.matrix[0:n])
#pragma acc enter data copyin(this->dqz_glb_rk3.matrix[0:n])
#pragma acc enter data copyin(this->dqx_glb_rk4.matrix[0:n])
#pragma acc enter data copyin(this->dqy_glb_rk4.matrix[0:n])
#pragma acc enter data copyin(this->dqz_glb_rk4.matrix[0:n])

	if (pt_glb->if_FE_all == true && pt_glb->if_flexo == true) {
#pragma acc enter data copyin(dpxdx, dpxdy, dpxdz)
#pragma acc enter data copyin(dpydx, dpydy, dpydz)
#pragma acc enter data copyin(dpzdx, dpzdy, dpzdz)

#pragma acc enter data copyin(dpxdx.matrix[0:n], dpxdy.matrix[0:n], dpxdz.matrix[0:n])
#pragma acc enter data copyin(dpydx.matrix[0:n], dpydy.matrix[0:n], dpydz.matrix[0:n])
#pragma acc enter data copyin(dpzdx.matrix[0:n], dpzdy.matrix[0:n], dpzdz.matrix[0:n])
	}

#pragma acc enter data copyin(this->Ex_stat)
#pragma acc enter data copyin(this->Ey_stat)
#pragma acc enter data copyin(this->Ez_stat)
#pragma acc enter data copyin(this->Ex_stat.matrix[0:n])
#pragma acc enter data copyin(this->Ey_stat.matrix[0:n])
#pragma acc enter data copyin(this->Ez_stat.matrix[0:n])

#pragma acc enter data copyin(this->rhof)
#pragma acc enter data copyin(this->rhof.matrix[0:n])

#pragma acc enter data copyin(this->rhofk)
#pragma acc enter data copyin(this->rhofk.matrix[0:n21])

#pragma acc host_data use_device(rhof.matrix, rhofk.matrix)
	{
		cufftExecD2Z(plan_electro_D2Z, rhof.matrix, (cufftDoubleComplex*)(rhofk.matrix));
	}

#pragma acc enter data copyin(this->Dxk_inhomo)
#pragma acc enter data copyin(this->Dyk_inhomo)
#pragma acc enter data copyin(this->Dzk_inhomo)
#pragma acc enter data copyin(this->Dxk_inhomo.matrix[0:n21])
#pragma acc enter data copyin(this->Dyk_inhomo.matrix[0:n21])
#pragma acc enter data copyin(this->Dzk_inhomo.matrix[0:n21])

#pragma acc enter data copyin(this->pxk)
#pragma acc enter data copyin(this->pyk)
#pragma acc enter data copyin(this->pzk)
#pragma acc enter data copyin(this->pxk.matrix[0:n21])
#pragma acc enter data copyin(this->pyk.matrix[0:n21])
#pragma acc enter data copyin(this->pzk.matrix[0:n21])

#pragma acc enter data copyin(this->ps_terms)
#pragma acc enter data copyin(this->ps_terms.matrix[0:n21])

#pragma acc enter data copyin(this->Edxk_n)
#pragma acc enter data copyin(this->Edyk_n)
#pragma acc enter data copyin(this->Edzk_n)
#pragma acc enter data copyin(this->Edxk_n.matrix[0:n21])
#pragma acc enter data copyin(this->Edyk_n.matrix[0:n21])
#pragma acc enter data copyin(this->Edzk_n.matrix[0:n21])

#pragma acc enter data copyin(this->Dx_inhomo)
#pragma acc enter data copyin(this->Dy_inhomo)
#pragma acc enter data copyin(this->Dz_inhomo)
#pragma acc enter data copyin(this->Dx_inhomo.matrix[0:n])
#pragma acc enter data copyin(this->Dy_inhomo.matrix[0:n])
#pragma acc enter data copyin(this->Dz_inhomo.matrix[0:n])

#pragma acc enter data copyin(this->e_field[0:4])
}

void ferroelectric_system::copyp_from_device() {
#pragma acc update host(px_glb.matrix[0:n],py_glb.matrix[0:n],pz_glb.matrix[0:n])
	;
}

void ferroelectric_system::copyq_from_device() {
#pragma acc update host(qx_glb.matrix[0:n],qy_glb.matrix[0:n],qz_glb.matrix[0:n])
	;
}

void ferroelectric_system::copyE_from_device() {
#pragma acc update host(Ex_stat.matrix[0:n],Ey_stat.matrix[0:n],Ez_stat.matrix[0:n])
	;
}
