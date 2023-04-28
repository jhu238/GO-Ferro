#include "ferroelectric_system.h"

#pragma acc routine seq// nohost
void ferroelectric_system::get_laplacian_p_local\
(long int& i, long int& j, long int& k, \
	double& d2pxdx2, double& d2pxdy2, double& d2pxdz2, \
	double& d2pydx2, double& d2pydy2, double& d2pydz2, \
	double& d2pzdx2, double& d2pzdy2, double& d2pzdz2) {

	double px_fwd, px_center, px_bwd;
	double py_fwd, py_center, py_bwd;
	double pz_fwd, pz_center, pz_bwd;

	px_center = px_glb_store(i, j, k);
	py_center = py_glb_store(i, j, k);
	pz_center = pz_glb_store(i, j, k);

	//--------------X direction (YZ surface)------------//
	//------------- negative X direction---------------//
	if (FE_surfYZ(i, j, k) == true) {
		px_bwd = px_center; py_bwd = py_center; pz_bwd = pz_center;
	}
	else if (FE_surfYZ(i, j, k) == false) {
		if (i == 0) {
			px_bwd = px_glb_store(nx - 1, j, k);
			py_bwd = py_glb_store(nx - 1, j, k);
			pz_bwd = pz_glb_store(nx - 1, j, k);
		}
		else {
			px_bwd = px_glb_store(i - 1, j, k);
			py_bwd = py_glb_store(i - 1, j, k);
			pz_bwd = pz_glb_store(i - 1, j, k);
		}
	}
	//------------- positive X direction---------------//
	if (FE_surfYZ(i + 1, j, k) == true) {
		px_fwd = px_center; py_fwd = py_center; pz_fwd = pz_center;
	}
	else if (FE_surfYZ(i + 1, j, k) == false) {
		if (i == nx - 1) {
			px_fwd = px_glb_store(0, j, k);
			py_fwd = py_glb_store(0, j, k);
			pz_fwd = pz_glb_store(0, j, k);
		}
		else {
			px_fwd = px_glb_store(i + 1, j, k);
			py_fwd = py_glb_store(i + 1, j, k);
			pz_fwd = pz_glb_store(i + 1, j, k);
		}
	}
	d2pxdx2 = (px_fwd - 2. * px_center + px_bwd) / dx / dx;
	d2pydx2 = (py_fwd - 2. * py_center + py_bwd) / dx / dx;
	d2pzdx2 = (pz_fwd - 2. * pz_center + pz_bwd) / dx / dx;
	//--------------END X direction (YZ surface)------------//

	//--------------Y direction (XZ surface)------------//
	//------------- negative Y direction---------------//
	if (FE_surfXZ(i, j, k) == true) {
		px_bwd = px_center; py_bwd = py_center; pz_bwd = pz_center;
	}
	else if (FE_surfXZ(i, j, k) == false) {
		if (j == 0) {
			px_bwd = px_glb_store(i, ny - 1, k);
			py_bwd = py_glb_store(i, ny - 1, k);
			pz_bwd = pz_glb_store(i, ny - 1, k);
		}
		else {
			px_bwd = px_glb_store(i, j - 1, k);
			py_bwd = py_glb_store(i, j - 1, k);
			pz_bwd = pz_glb_store(i, j - 1, k);
		}
	}
	//------------- positive Y direction---------------//
	if (FE_surfXZ(i, j + 1, k) == true) {
		px_fwd = px_center; py_fwd = py_center; pz_fwd = pz_center;
	}
	else if (FE_surfXZ(i, j + 1, k) == false) {
		if (j == ny - 1) {
			px_fwd = px_glb_store(i, 0, k);
			py_fwd = py_glb_store(i, 0, k);
			pz_fwd = pz_glb_store(i, 0, k);
		}
		else {
			px_fwd = px_glb_store(i, j + 1, k);
			py_fwd = py_glb_store(i, j + 1, k);
			pz_fwd = pz_glb_store(i, j + 1, k);
		}
	}
	d2pxdy2 = (px_fwd - 2. * px_center + px_bwd) / dy / dy;
	d2pydy2 = (py_fwd - 2. * py_center + py_bwd) / dy / dy;
	d2pzdy2 = (pz_fwd - 2. * pz_center + pz_bwd) / dy / dy;
	//--------------END Y direction (XZ surface)------------//	

	//--------------Z direction (XY surface)------------//
	//------------- negative Z direction---------------//
	if (FE_surfXY(i, j, k) == true) {
		px_bwd = px_center; py_bwd = py_center; pz_bwd = pz_center;
	}
	else if (FE_surfXY(i, j, k) == false) {
		if (k == 0) {
			px_bwd = px_glb_store(i, j, nz - 1);
			py_bwd = py_glb_store(i, j, nz - 1);
			pz_bwd = pz_glb_store(i, j, nz - 1);
		}
		else {
			px_bwd = px_glb_store(i, j, k - 1);
			py_bwd = py_glb_store(i, j, k - 1);
			pz_bwd = pz_glb_store(i, j, k - 1);
		}
	}
	//------------- positive Z direction---------------//
	if (FE_surfXY(i, j, k + 1) == true) {
		px_fwd = px_center; py_fwd = py_center; pz_fwd = pz_center;
	}
	else if (FE_surfXY(i, j, k + 1) == false) {
		if (k == nz - 1) {
			px_fwd = px_glb_store(i, j, 0);
			py_fwd = py_glb_store(i, j, 0);
			pz_fwd = pz_glb_store(i, j, 0);
		}
		else {
			px_fwd = px_glb_store(i, j, k + 1);
			py_fwd = py_glb_store(i, j, k + 1);
			pz_fwd = pz_glb_store(i, j, k + 1);
		}
	}
	d2pxdz2 = (px_fwd - 2. * px_center + px_bwd) / dz / dz;
	d2pydz2 = (py_fwd - 2. * py_center + py_bwd) / dz / dz;
	d2pzdz2 = (pz_fwd - 2. * pz_center + pz_bwd) / dz / dz;
	//--------------END Z direction (XY surface)------------//
}

#pragma acc routine gang
void ferroelectric_system::get_p_gradient() {
	long i, j, k;
	unsigned int mat_type;
	material* mat;
	double px_fwd, px_bwd;
	double py_fwd, py_bwd;
	double pz_fwd, pz_bwd;

#pragma acc loop gang vector private(i,j,k,px_fwd,px_bwd,py_fwd,py_bwd,pz_fwd,pz_bwd, mat, mat_type)
	for (long id = 0; id < n; id++) {
		mat_type = pt_glb->material_cell(id);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[mat_type - 1]);
			if (mat->if_FE == true) {
				i = id / (ny * nz);
				j = (id - i * (ny * nz)) / nz;
				k = id - i * (ny * nz) - j * nz;

				//--------------X direction (YZ surface)------------//
				//------------- negative X direction---------------//
				if (FE_surfYZ(i, j, k) == true) {
					px_bwd = px_glb_store(id);
					py_bwd = py_glb_store(id);
					pz_bwd = pz_glb_store(id);
				}
				else if (FE_surfYZ(i, j, k) == false) {
					if (i == 0) {
						px_bwd = px_glb_store(nx - 1, j, k);
						py_bwd = py_glb_store(nx - 1, j, k);
						pz_bwd = pz_glb_store(nx - 1, j, k);
					}
					else {
						px_bwd = px_glb_store(i - 1, j, k);
						py_bwd = py_glb_store(i - 1, j, k);
						pz_bwd = pz_glb_store(i - 1, j, k);
					}
				}
				//------------- positive X direction---------------//
				if (FE_surfYZ(i + 1, j, k) == true) {
					px_fwd = px_glb_store(id);
					py_fwd = py_glb_store(id);
					pz_fwd = pz_glb_store(id);
				}
				else if (FE_surfYZ(i + 1, j, k) == false) {
					if (i == nx - 1) {
						px_fwd = px_glb_store(0, j, k);
						py_fwd = py_glb_store(0, j, k);
						pz_fwd = pz_glb_store(0, j, k);
					}
					else {
						px_fwd = px_glb_store(i + 1, j, k);
						py_fwd = py_glb_store(i + 1, j, k);
						pz_fwd = pz_glb_store(i + 1, j, k);
					}
				}
				dpxdx(id) = (px_fwd - px_bwd) / dx / 2.;
				dpydx(id) = (py_fwd - py_bwd) / dx / 2.;
				dpzdx(id) = (pz_fwd - pz_bwd) / dx / 2.;
				//--------------END X direction (YZ surface)------------//

				//--------------Y direction (XZ surface)------------//
				//------------- negative Y direction---------------//
				if (FE_surfXZ(i, j, k) == true) {
					px_bwd = px_glb_store(id);
					py_bwd = py_glb_store(id);
					pz_bwd = pz_glb_store(id);
				}
				else if (FE_surfXZ(i, j, k) == false) {
					if (j == 0) {
						px_bwd = px_glb_store(i, ny - 1, k);
						py_bwd = py_glb_store(i, ny - 1, k);
						pz_bwd = pz_glb_store(i, ny - 1, k);
					}
					else {
						px_bwd = px_glb_store(i, j - 1, k);
						py_bwd = py_glb_store(i, j - 1, k);
						pz_bwd = pz_glb_store(i, j - 1, k);
					}
				}
				//------------- positive Y direction---------------//
				if (FE_surfXZ(i, j + 1, k) == true) {
					px_fwd = px_glb_store(id);
					py_fwd = py_glb_store(id);
					pz_fwd = pz_glb_store(id);
				}
				else if (FE_surfXZ(i, j + 1, k) == false) {
					if (j == ny - 1) {
						px_fwd = px_glb_store(i, 0, k);
						py_fwd = py_glb_store(i, 0, k);
						pz_fwd = pz_glb_store(i, 0, k);
					}
					else {
						px_fwd = px_glb_store(i, j + 1, k);
						py_fwd = py_glb_store(i, j + 1, k);
						pz_fwd = pz_glb_store(i, j + 1, k);
					}
				}
				dpxdy(id) = (px_fwd - px_bwd) / dy / 2.;
				dpydy(id) = (py_fwd - py_bwd) / dy / 2.;
				dpzdy(id) = (pz_fwd - pz_bwd) / dy / 2.;
				//--------------END Y direction (XZ surface)------------//	

				//--------------Z direction (XY surface)------------//
				//------------- negative Z direction---------------//
				if (FE_surfXY(i, j, k) == true) {
					px_bwd = px_glb_store(id);
					py_bwd = py_glb_store(id);
					pz_bwd = pz_glb_store(id);
				}
				else if (FE_surfXY(i, j, k) == false) {
					if (k == 0) {
						px_bwd = px_glb_store(i, j, nz - 1);
						py_bwd = py_glb_store(i, j, nz - 1);
						pz_bwd = pz_glb_store(i, j, nz - 1);
					}
					else {
						px_bwd = px_glb_store(i, j, k - 1);
						py_bwd = py_glb_store(i, j, k - 1);
						pz_bwd = pz_glb_store(i, j, k - 1);
					}
				}
				//------------- positive Z direction---------------//
				if (FE_surfXY(i, j, k + 1) == true) {
					px_fwd = px_glb_store(id);
					py_fwd = py_glb_store(id);
					pz_fwd = pz_glb_store(id);
				}
				else if (FE_surfXY(i, j, k + 1) == false) {
					if (k == nz - 1) {
						px_fwd = px_glb_store(i, j, 0);
						py_fwd = py_glb_store(i, j, 0);
						pz_fwd = pz_glb_store(i, j, 0);
					}
					else {
						px_fwd = px_glb_store(i, j, k + 1);
						py_fwd = py_glb_store(i, j, k + 1);
						pz_fwd = pz_glb_store(i, j, k + 1);
					}
				}
				dpxdz(id) = (px_fwd - px_bwd) / dz / 2.;
				dpydz(id) = (py_fwd - py_bwd) / dz / 2.;
				dpzdz(id) = (pz_fwd - pz_bwd) / dz / 2.;
				//--------------END Z direction (XY surface)------------//
			}
		}
	}

}

//void ferroelectric_system::get_dq() {
//	double dqx_glb_local, dqy_glb_local, dqz_glb_local;
//	double dqx_crt_local, dqy_crt_local, dqz_crt_local;
//	double px_glb_local, py_glb_local, pz_glb_local;
//	double px, py, pz;
//
//	double d2pxdx2_glb, d2pxdy2_glb, d2pxdz2_glb;
//	double d2pydx2_glb, d2pydy2_glb, d2pydz2_glb;
//	double d2pzdx2_glb, d2pzdy2_glb, d2pzdz2_glb;
//
//	double dexxdx, dexxdy, dexxdz;
//	double deyydx, deyydy, deyydz;
//	double dezzdx, dezzdy, dezzdz;
//	double deyzdx, deyzdy, deyzdz;
//	double dexzdx, dexzdy, dexzdz;
//	double dexydx, dexydy, dexydz;
//
//	unsigned int mat_type;
//	material* mat;
//	double exx, eyy, ezz, eyz, exz, exy;
//	double exx0, eyy0, ezz0, eyz0, exz0, exy0;
//
//	double a1, \
//		a11, a12, \
//		a111, a112, a123, \
//		a1111, a1112, a1122, a1123;
//	double tv1, tv2, tv3, tv4, tv5, tv6, tv7, tv8;
//	double f11, f12, f44;
//	double G11;
//
//	double FE_mass, FE_damping;
//	long i, j, k;
//
//#pragma acc parallel default(present)
//	{
//#pragma acc loop gang vector private(\
//	dqx_glb_local, dqy_glb_local, dqz_glb_local,\
//	dqx_crt_local, dqy_crt_local, dqz_crt_local,\
//	px_glb_local, py_glb_local, pz_glb_local,\
//	px, py, pz,\
//	d2pxdx2_glb, d2pxdy2_glb, d2pxdz2_glb,\
//	d2pydx2_glb, d2pydy2_glb, d2pydz2_glb,\
//	d2pzdx2_glb, d2pzdy2_glb, d2pzdz2_glb,\
//	mat_type,\
//	exx, eyy, ezz, eyz, exz, exy,\
//	exx0, eyy0, ezz0, eyz0, exz0, exy0,\
//	a1, \
//	a11, a12, \
//	a111, a112, a123, \
//	a1111, a1112, a1122, a1123,\
//	tv1, tv2, tv3, tv4, tv5, tv6, tv7, tv8,\
//	f11, f12, f44,\
//	G11,\
//	FE_mass, FE_damping,mat,i,j,k,\
//	dexxdx, dexxdy, dexxdz,\
//	deyydx, deyydy, deyydz,\
//	dezzdx, dezzdy, dezzdz,\
//	deyzdx, deyzdy, deyzdz,\
//	dexzdx, dexzdy, dexzdz,\
//	dexydx, dexydy, dexydz)
//		for (long int id = 0; id < n; id++) {
//			i = id / (ny * nz);
//			j = (id - i * (ny * nz)) / nz;
//			k = id - i * (ny * nz) - j * nz;
//
//			mat_type = pt_glb->material_cell(i, j, k);
//			if (mat_type != 0) {
//				mat = &(pt_glb->material_parameters[mat_type - 1]);
//				if (mat->if_FE == true) {
//					px_glb_local = px_glb(i, j, k);
//					py_glb_local = py_glb(i, j, k);
//					pz_glb_local = pz_glb(i, j, k);
//
//					get_laplacian_p_local \
//						(i, j, k, \
//							d2pxdx2_glb, d2pxdy2_glb, d2pxdz2_glb, \
//							d2pydx2_glb, d2pydy2_glb, d2pydz2_glb, \
//							d2pzdx2_glb, d2pzdy2_glb, d2pzdz2_glb);
//
//					pt_math->transform_vector_glb2crt(px_glb_local, py_glb_local, pz_glb_local, \
//						px, py, pz);
//
//					exx = pt_elas->exxt0_crt(i, j, k) + pt_elas->Dexx_crt(i, j, k);
//					eyy = pt_elas->eyyt0_crt(i, j, k) + pt_elas->Deyy_crt(i, j, k);
//					ezz = pt_elas->ezzt0_crt(i, j, k) + pt_elas->Dezz_crt(i, j, k);
//					eyz = pt_elas->eyzt0_crt(i, j, k) + pt_elas->Deyz_crt(i, j, k);
//					exz = pt_elas->exzt0_crt(i, j, k) + pt_elas->Dexz_crt(i, j, k);
//					exy = pt_elas->exyt0_crt(i, j, k) + pt_elas->Dexy_crt(i, j, k);
//
//					exx0 = pt_elas->exx0t0_crt(i, j, k) + pt_elas->Dexx0_crt(i, j, k);
//					eyy0 = pt_elas->eyy0t0_crt(i, j, k) + pt_elas->Deyy0_crt(i, j, k);
//					ezz0 = pt_elas->ezz0t0_crt(i, j, k) + pt_elas->Dezz0_crt(i, j, k);
//					eyz0 = pt_elas->eyz0t0_crt(i, j, k) + pt_elas->Deyz0_crt(i, j, k);
//					exz0 = pt_elas->exz0t0_crt(i, j, k) + pt_elas->Dexz0_crt(i, j, k);
//					exy0 = pt_elas->exy0t0_crt(i, j, k) + pt_elas->Dexy0_crt(i, j, k);
//
//					if (pt_glb->if_flexo == true) {
//						dexxdx = pt_elas->dexxt0dx_glb(i, j, k) + pt_elas->dDexxdx_glb(i, j, k);
//						dexxdy = pt_elas->dexxt0dy_glb(i, j, k) + pt_elas->dDexxdy_glb(i, j, k);
//						dexxdz = pt_elas->dexxt0dz_glb(i, j, k) + pt_elas->dDexxdz_glb(i, j, k);
//						deyydx = pt_elas->deyyt0dx_glb(i, j, k) + pt_elas->dDeyydx_glb(i, j, k);
//						deyydy = pt_elas->deyyt0dy_glb(i, j, k) + pt_elas->dDeyydy_glb(i, j, k);
//						deyydz = pt_elas->deyyt0dz_glb(i, j, k) + pt_elas->dDeyydz_glb(i, j, k);
//						dezzdx = pt_elas->dezzt0dx_glb(i, j, k) + pt_elas->dDezzdx_glb(i, j, k);
//						dezzdy = pt_elas->dezzt0dy_glb(i, j, k) + pt_elas->dDezzdy_glb(i, j, k);
//						dezzdz = pt_elas->dezzt0dz_glb(i, j, k) + pt_elas->dDezzdz_glb(i, j, k);
//						deyzdx = pt_elas->deyzt0dx_glb(i, j, k) + pt_elas->dDeyzdx_glb(i, j, k);
//						deyzdy = pt_elas->deyzt0dy_glb(i, j, k) + pt_elas->dDeyzdy_glb(i, j, k);
//						deyzdz = pt_elas->deyzt0dz_glb(i, j, k) + pt_elas->dDeyzdz_glb(i, j, k);
//						dexzdx = pt_elas->dexzt0dx_glb(i, j, k) + pt_elas->dDexzdx_glb(i, j, k);
//						dexzdy = pt_elas->dexzt0dy_glb(i, j, k) + pt_elas->dDexzdy_glb(i, j, k);
//						dexzdz = pt_elas->dexzt0dz_glb(i, j, k) + pt_elas->dDexzdz_glb(i, j, k);
//						dexydx = pt_elas->dexyt0dx_glb(i, j, k) + pt_elas->dDexydx_glb(i, j, k);
//						dexydy = pt_elas->dexyt0dy_glb(i, j, k) + pt_elas->dDexydy_glb(i, j, k);
//						dexydz = pt_elas->dexyt0dz_glb(i, j, k) + pt_elas->dDexydz_glb(i, j, k);
//					}
//
//					a1 = mat->a1;
//					a11 = mat->a11; a12 = mat->a12;
//					a111 = mat->a111; a112 = mat->a112; a123 = mat->a123;
//					a1111 = mat->a1111; a1112 = mat->a1112; a1122 = mat->a1122; a1123 = mat->a1123;
//
//					tv1 = mat->tv1;	tv2 = mat->tv2; tv3 = mat->tv3; tv4 = mat->tv4;
//					tv5 = mat->tv5; tv6 = mat->tv6; tv7 = mat->tv7; tv8 = mat->tv8;
//
//					f11 = mat->f11; f12 = mat->f12; f44 = mat->f44;
//
//					G11 = mat->G11;
//
//					FE_mass = -1. * pt_glb->dt / mat->FE_mass;
//					FE_damping = FE_mass * mat->FE_damping;
//
//					//--------------------X compoennt-------------------//
//					//---------Landau energy effective field------------//
//					dqx_crt_local = FE_mass * (\
//						2. * a1 * px + 4. * a11 * pow(px, 3.) \
//						+ 2. * a12 * (pow(py, 2.) + pow(pz, 2.)) * px \
//						+ 6. * a111 * pow(px, 5.) \
//						+ 4. * a112 * pow(px, 3.) * (pow(py, 2.) + pow(pz, 2.)) \
//						+ 2. * a112 * (pow(py, 4.) + pow(pz, 4.)) * px \
//						+ 2. * a123 * px * pow(py, 2.) * pow(pz, 2.) \
//						+ 8. * a1111 * pow(px, 7.) \
//						+ 2. * a1112 * px * (pow(py, 6.) + pow(pz, 6.) \
//							+ 3. * pow(px, 4.) * (pow(py, 2.) + pow(pz, 2.))) \
//						+ 4. * a1122 * pow(px, 3.) * (pow(py, 4.) + pow(pz, 4.)) \
//						+ 2. * a1123 * px * pow(py, 2.) * pow(pz, 2.) * \
//						(2. * pow(px, 2.) + pow(py, 2.) + pow(pz, 2.)) \
//						//--------Elastic energy effective field----------//
//						+ (tv1 * (exx0 - exx) + tv2 * (eyy0 - eyy + ezz0 - ezz)) * px \
//						+ tv3 * ((exy0 - exy) * py + (exz0 - exz) * pz));
//
//					//--------------------Y compoennt-------------------//	
//					//---------Landau energy effective field------------//
//					dqy_crt_local = FE_mass * (\
//						2. * a1 * py + 4. * a11 * pow(py, 3.) \
//						+ 2. * a12 * (pow(pz, 2.) + pow(px, 2.)) * py \
//						+ 6. * a111 * pow(py, 5.) \
//						+ 4. * a112 * pow(py, 3.) * (pow(px, 2.) + pow(pz, 2.)) \
//						+ 2. * a112 * (pow(px, 4.) + pow(pz, 4.)) * py \
//						+ 2. * a123 * py * pow(px, 2.) * pow(pz, 2.) \
//						+ 8. * a1111 * pow(py, 7) \
//						+ 2. * a1112 * py * (pow(px, 6.) + pow(pz, 6.) \
//							+ 3. * pow(py, 4.) * (pow(px, 2.) + pow(pz, 2.))) \
//						+ 4. * a1122 * pow(py, 3.) * (pow(px, 4.) + pow(pz, 4.)) \
//						+ 2. * a1123 * py * pow(px, 2.) * pow(pz, 2.) * \
//						(pow(px, 2.) + 2 * pow(py, 2.) + pow(pz, 2.))\
//						//--------Elastic energy effective field----------//
//						+ (tv1 * (eyy0 - eyy) + tv2 * (exx0 - exx + ezz0 - ezz)) * py \
//						+ tv3 * ((exy0 - exy) * px + (eyz0 - eyz) * pz));
//
//					//--------------------Z compoennt-------------------//	
//					//---------Landau energy effective field------------//
//					dqz_crt_local = FE_mass * (\
//						2. * a1 * pz + 4. * a11 * pow(pz, 3.) \
//						+ 2. * a12 * (pow(py, 2.) + pow(px, 2.)) * pz \
//						+ 6. * a111 * pow(pz, 5.) \
//						+ 4. * a112 * pow(pz, 3.) * (pow(py, 2.) + pow(px, 2.)) \
//						+ 2. * a112 * (pow(py, 4.) + pow(px, 4.)) * pz \
//						+ 2. * a123 * pz * pow(py, 2.) * pow(px, 2.) \
//						+ 8. * a1111 * pow(pz, 7.) \
//						+ 2. * a1112 * pz * (pow(py, 6.) + pow(px, 6.) \
//							+ 3. * pow(pz, 4.) * (pow(py, 2.) + pow(px, 2.))) \
//						+ 4. * a1122 * pow(pz, 3.) * (pow(py, 4.) + pow(px, 4.)) \
//						+ 2. * a1123 * pz * pow(px, 2.) * pow(py, 2.) * \
//						(pow(px, 2.) + pow(py, 2.) + 2. * pow(pz, 2.)) \
//						//--------Elastic energy effective field----------//
//						+ (tv1 * (ezz0 - ezz) + tv2 * (exx0 - exx + eyy0 - eyy)) * pz \
//						+ tv3 * ((exz0 - exz) * px + (eyz0 - eyz) * py));
//
//					pt_math->transform_vector_crt2glb(dqx_crt_local, dqy_crt_local, dqz_crt_local, \
//						dqx_glb_local, dqy_glb_local, dqz_glb_local);
//
//					dqx_glb_local = dqx_glb_local + FE_mass * (\
//						- 1. * G11 * (d2pxdx2_glb + d2pxdy2_glb + d2pxdz2_glb) \
//						- Ex_stat(i, j, k) - pt_EM->DEx_em_cell(i, j, k) - pt_glb->Eext[0]) \
//						+ FE_damping * qx_glb(i, j, k);
//
//					dqy_glb_local = dqy_glb_local + FE_mass * (\
//						- 1. * G11 * (d2pydx2_glb + d2pydy2_glb + d2pydz2_glb) \
//						- Ey_stat(i, j, k) - pt_EM->DEy_em_cell(i, j, k) - pt_glb->Eext[1]) \
//						+ FE_damping * qy_glb(i, j, k);
//
//					dqz_glb_local = dqz_glb_local + FE_mass * (\
//						- 1. * G11 * (d2pzdx2_glb + d2pzdy2_glb + d2pzdz2_glb) \
//						- Ez_stat(i, j, k) - pt_EM->DEz_em_cell(i, j, k) - pt_glb->Eext[2]) \
//						+ FE_damping * qz_glb(i, j, k);
//
//					if (pt_glb->if_flexo == true) {
//						//--------------------X compoennt-------------------//
//						dqx_glb_local = dqx_glb_local + FE_mass * (\
//							+ tv4 * px_glb_local * dpxdx(i,j,k) \
//							+ tv5 * (py_glb_local * dpydx(i, j, k) + pz_glb_local * dpzdx(i, j, k)) \
//							+ tv6 * (px_glb_local * dpydy(i, j, k) + px_glb_local * dpzdz(i, j, k) + py_glb_local * dpxdy(i, j, k) + pz_glb_local * dpxdz(i, j, k)) \
//							-tv7 * d2pxdx2_glb \
//							-tv8 * (d2pxdz2_glb + d2pxdy2_glb) \
//							//-------------Flexoelectric energy effectiv field--------//
//							-f11 * dexxdx - f12 * (deyydx + dezzdx) \
//							- 2. * f44 * (dexydy + dexzdz));
//
//						//--------------------Y compoennt-------------------//
//						dqy_glb_local = dqy_glb_local + FE_mass * (\
//							+ tv4 * py_glb_local * dpydy(i, j, k) \
//							+tv5 * (px_glb_local * dpxdy(i,j,k) + pz_glb_local * dpzdy(i,j,k)) \
//							+tv6 * (py_glb_local * dpxdx(i,j,k) + py_glb_local * dpzdz(i,j,k) + px_glb_local * dpydx(i, j, k) + pz_glb_local * dpydz(i, j, k)) \
//							-tv7 * d2pydy2_glb \
//							-tv8 * (d2pydx2_glb + d2pydz2_glb) \
//							//-------------Flexoelectric energy effectiv field--------//
//							-f11 * deyydy - f12 * (dexxdy + dezzdy) \
//							- 2. * f44 * (dexydx + deyzdz));
//
//						//--------------------Z compoennt-------------------//
//						dqz_glb_local = dqz_glb_local + FE_mass * (\
//							+ tv4 * pz_glb_local * dpzdz(i, j, k) \
//							+ tv5 * (px_glb_local * dpxdz(i,j,k) + py_glb_local * dpydz(i,j,k)) \
//							+ tv6 * (pz_glb_local * dpxdx(i,j,k) + pz_glb_local * dpydy(i,j,k) + px_glb_local * dpzdx(i, j, k) + py_glb_local * dpzdy(i, j, k)) \
//							- tv7 * d2pzdz2_glb\
//							- tv8 * (d2pzdx2_glb + d2pzdy2_glb) \
//							//-------------Flexoelectric energy effectiv field--------//
//							-f11 * dezzdz - f12 * (dexxdz + deyydz) \
//							- 2. * f44 * (deyzdy + dexzdx));
//					}
//
//					dqx_glb(i, j, k) = dqx_glb_local;
//					dqy_glb(i, j, k) = dqy_glb_local;
//					dqz_glb(i, j, k) = dqz_glb_local;
//				}
//			}
//		}
//	}
//}

void ferroelectric_system::get_dq_RK1() {
	double dqx_glb_local, dqy_glb_local, dqz_glb_local;
	double dqx_crt_local, dqy_crt_local, dqz_crt_local;
	double px_glb_local, py_glb_local, pz_glb_local;
	double px, py, pz;

	double d2pxdx2_glb, d2pxdy2_glb, d2pxdz2_glb;
	double d2pydx2_glb, d2pydy2_glb, d2pydz2_glb;
	double d2pzdx2_glb, d2pzdy2_glb, d2pzdz2_glb;

	double dexxdx, dexxdy, dexxdz;
	double deyydx, deyydy, deyydz;
	double dezzdx, dezzdy, dezzdz;
	double deyzdx, deyzdy, deyzdz;
	double dexzdx, dexzdy, dexzdz;
	double dexydx, dexydy, dexydz;

	unsigned int mat_type;
	material* mat;
	double exx, eyy, ezz, eyz, exz, exy;
	double exx0, eyy0, ezz0, eyz0, exz0, exy0;

	double a1, \
		a11, a12, \
		a111, a112, a123, \
		a1111, a1112, a1122, a1123;
	double tv1, tv2, tv3, tv4, tv5, tv6, tv7, tv8;
	double f11, f12, f44;
	double G11;

	double FE_mass, FE_damping;
	long i, j, k;

	if (pt_glb->if_prescribe_Eext == false) {
#pragma acc serial default(present) async(5)
		{
			update_external_Efield();
		}
	}

#pragma acc wait(1) async(5)
#pragma acc parallel default(present) async(5)
	{
#pragma acc loop gang vector private(\
	dqx_glb_local, dqy_glb_local, dqz_glb_local,\
	dqx_crt_local, dqy_crt_local, dqz_crt_local,\
	px_glb_local, py_glb_local, pz_glb_local,\
	px, py, pz,\
	d2pxdx2_glb, d2pxdy2_glb, d2pxdz2_glb,\
	d2pydx2_glb, d2pydy2_glb, d2pydz2_glb,\
	d2pzdx2_glb, d2pzdy2_glb, d2pzdz2_glb,\
	mat_type,\
	exx, eyy, ezz, eyz, exz, exy,\
	exx0, eyy0, ezz0, eyz0, exz0, exy0,\
	a1, \
	a11, a12, \
	a111, a112, a123, \
	a1111, a1112, a1122, a1123,\
	tv1, tv2, tv3, tv4, tv5, tv6, tv7, tv8,\
	f11, f12, f44,\
	G11,\
	FE_mass, FE_damping,mat,i,j,k,\
	dexxdx, dexxdy, dexxdz,\
	deyydx, deyydy, deyydz,\
	dezzdx, dezzdy, dezzdz,\
	deyzdx, deyzdy, deyzdz,\
	dexzdx, dexzdy, dexzdz,\
	dexydx, dexydy, dexydz)
		for (long int id = 0; id < n; id++) {
			i = id / (ny * nz);
			j = (id - i * (ny * nz)) / nz;
			k = id - i * (ny * nz) - j * nz;

			mat_type = pt_glb->material_cell(i, j, k);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[mat_type - 1]);
				if (mat->if_FE == true) {
					px_glb_local = px_glb_store(i, j, k);
					py_glb_local = py_glb_store(i, j, k);
					pz_glb_local = pz_glb_store(i, j, k);

					get_laplacian_p_local \
						(i, j, k, \
							d2pxdx2_glb, d2pxdy2_glb, d2pxdz2_glb, \
							d2pydx2_glb, d2pydy2_glb, d2pydz2_glb, \
							d2pzdx2_glb, d2pzdy2_glb, d2pzdz2_glb);

					pt_math->transform_vector_glb2crt(px_glb_local, py_glb_local, pz_glb_local, \
						px, py, pz);

					if (pt_glb->if_elasto_on_pORm == true) {
						exx = pt_elas->Dexx_crt(i, j, k);
						eyy = pt_elas->Deyy_crt(i, j, k);
						ezz = pt_elas->Dezz_crt(i, j, k);
						eyz = pt_elas->Deyz_crt(i, j, k);
						exz = pt_elas->Dexz_crt(i, j, k);
						exy = pt_elas->Dexy_crt(i, j, k);

						exx0 = pt_elas->Dexx0_crt(i, j, k);
						eyy0 = pt_elas->Deyy0_crt(i, j, k);
						ezz0 = pt_elas->Dezz0_crt(i, j, k);
						eyz0 = pt_elas->Deyz0_crt(i, j, k);
						exz0 = pt_elas->Dexz0_crt(i, j, k);
						exy0 = pt_elas->Dexy0_crt(i, j, k);
						if (pt_glb->if_elastostatic == true) {
							exx = exx + pt_elas->exxt0_crt(i, j, k);
							eyy = eyy + pt_elas->eyyt0_crt(i, j, k);
							ezz = ezz + pt_elas->ezzt0_crt(i, j, k);
							eyz = eyz + pt_elas->eyzt0_crt(i, j, k);
							exz = exz + pt_elas->exzt0_crt(i, j, k);
							exy = exy + pt_elas->exyt0_crt(i, j, k);

							exx0 = exx0 + pt_elas->exx0t0_crt(i, j, k);
							eyy0 = eyy0 + pt_elas->eyy0t0_crt(i, j, k);
							ezz0 = ezz0 + pt_elas->ezz0t0_crt(i, j, k);
							eyz0 = eyz0 + pt_elas->eyz0t0_crt(i, j, k);
							exz0 = exz0 + pt_elas->exz0t0_crt(i, j, k);
							exy0 = exy0 + pt_elas->exy0t0_crt(i, j, k);
						}
					}
					else {
						exx = 0.; exx0 = 0;
						eyy = 0.; eyy0 = 0;
						ezz = 0.; ezz0 = 0;
						eyz = 0.; eyz0 = 0;
						exz = 0.; exz0 = 0;
						exy = 0.; exy0 = 0;
					}

					if (pt_glb->if_flexo == true) {
						dexxdx = pt_elas->dDexxdx_glb(i, j, k);
						dexxdy = pt_elas->dDexxdy_glb(i, j, k);
						dexxdz = pt_elas->dDexxdz_glb(i, j, k);
						deyydx = pt_elas->dDeyydx_glb(i, j, k);
						deyydy = pt_elas->dDeyydy_glb(i, j, k);
						deyydz = pt_elas->dDeyydz_glb(i, j, k);
						dezzdx = pt_elas->dDezzdx_glb(i, j, k);
						dezzdy = pt_elas->dDezzdy_glb(i, j, k);
						dezzdz = pt_elas->dDezzdz_glb(i, j, k);
						deyzdx = pt_elas->dDeyzdx_glb(i, j, k);
						deyzdy = pt_elas->dDeyzdy_glb(i, j, k);
						deyzdz = pt_elas->dDeyzdz_glb(i, j, k);
						dexzdx = pt_elas->dDexzdx_glb(i, j, k);
						dexzdy = pt_elas->dDexzdy_glb(i, j, k);
						dexzdz = pt_elas->dDexzdz_glb(i, j, k);
						dexydx = pt_elas->dDexydx_glb(i, j, k);
						dexydy = pt_elas->dDexydy_glb(i, j, k);
						dexydz = pt_elas->dDexydz_glb(i, j, k);
						if (pt_glb->if_elastostatic == true) {
							dexxdx = dexxdx + pt_elas->dexxt0dx_glb(i, j, k);
							dexxdy = dexxdy + pt_elas->dexxt0dy_glb(i, j, k);
							dexxdz = dexxdz + pt_elas->dexxt0dz_glb(i, j, k);
							deyydx = deyydx + pt_elas->deyyt0dx_glb(i, j, k);
							deyydy = deyydy + pt_elas->deyyt0dy_glb(i, j, k);
							deyydz = deyydz + pt_elas->deyyt0dz_glb(i, j, k);
							dezzdx = dezzdx + pt_elas->dezzt0dx_glb(i, j, k);
							dezzdy = dezzdy + pt_elas->dezzt0dy_glb(i, j, k);
							dezzdz = dezzdz + pt_elas->dezzt0dz_glb(i, j, k);
							deyzdx = deyzdx + pt_elas->deyzt0dx_glb(i, j, k);
							deyzdy = deyzdy + pt_elas->deyzt0dy_glb(i, j, k);
							deyzdz = deyzdz + pt_elas->deyzt0dz_glb(i, j, k);
							dexzdx = dexzdx + pt_elas->dexzt0dx_glb(i, j, k);
							dexzdy = dexzdy + pt_elas->dexzt0dy_glb(i, j, k);
							dexzdz = dexzdz + pt_elas->dexzt0dz_glb(i, j, k);
							dexydx = dexydx + pt_elas->dexyt0dx_glb(i, j, k);
							dexydy = dexydy + pt_elas->dexyt0dy_glb(i, j, k);
							dexydz = dexydz + pt_elas->dexyt0dz_glb(i, j, k);
						}
					}

					a1 = mat->a1;
					a11 = mat->a11; a12 = mat->a12;
					a111 = mat->a111; a112 = mat->a112; a123 = mat->a123;
					a1111 = mat->a1111; a1112 = mat->a1112; a1122 = mat->a1122; a1123 = mat->a1123;

					tv1 = mat->tv1;	tv2 = mat->tv2; tv3 = mat->tv3; tv4 = mat->tv4;
					tv5 = mat->tv5; tv6 = mat->tv6; tv7 = mat->tv7; tv8 = mat->tv8;

					f11 = mat->f11; f12 = mat->f12; f44 = mat->f44;

					G11 = mat->G11;

					FE_mass = -1. * pt_glb->dt / mat->FE_mass;
					FE_damping = FE_mass * mat->FE_damping;

					//--------------------X compoennt-------------------//
					//---------Landau energy effective field------------//
					dqx_crt_local = FE_mass * (\
						2. * a1 * px + 4. * a11 * pow(px, 3.) \
						+ 2. * a12 * (pow(py, 2.) + pow(pz, 2.)) * px \
						+ 6. * a111 * pow(px, 5.) \
						+ 4. * a112 * pow(px, 3.) * (pow(py, 2.) + pow(pz, 2.)) \
						+ 2. * a112 * (pow(py, 4.) + pow(pz, 4.)) * px \
						+ 2. * a123 * px * pow(py, 2.) * pow(pz, 2.) \
						+ 8. * a1111 * pow(px, 7.) \
						+ 2. * a1112 * px * (pow(py, 6.) + pow(pz, 6.) \
							+ 3. * pow(px, 4.) * (pow(py, 2.) + pow(pz, 2.))) \
						+ 4. * a1122 * pow(px, 3.) * (pow(py, 4.) + pow(pz, 4.)) \
						+ 2. * a1123 * px * pow(py, 2.) * pow(pz, 2.) * \
						(2. * pow(px, 2.) + pow(py, 2.) + pow(pz, 2.)) \
						//--------Elastic energy effective field----------//
						+ (tv1 * (exx0 - exx) + tv2 * (eyy0 - eyy + ezz0 - ezz)) * px \
						+ tv3 * ((exy0 - exy) * py + (exz0 - exz) * pz));

					//--------------------Y compoennt-------------------//	
					//---------Landau energy effective field------------//
					dqy_crt_local = FE_mass * (\
						2. * a1 * py + 4. * a11 * pow(py, 3.) \
						+ 2. * a12 * (pow(pz, 2.) + pow(px, 2.)) * py \
						+ 6. * a111 * pow(py, 5.) \
						+ 4. * a112 * pow(py, 3.) * (pow(px, 2.) + pow(pz, 2.)) \
						+ 2. * a112 * (pow(px, 4.) + pow(pz, 4.)) * py \
						+ 2. * a123 * py * pow(px, 2.) * pow(pz, 2.) \
						+ 8. * a1111 * pow(py, 7) \
						+ 2. * a1112 * py * (pow(px, 6.) + pow(pz, 6.) \
							+ 3. * pow(py, 4.) * (pow(px, 2.) + pow(pz, 2.))) \
						+ 4. * a1122 * pow(py, 3.) * (pow(px, 4.) + pow(pz, 4.)) \
						+ 2. * a1123 * py * pow(px, 2.) * pow(pz, 2.) * \
						(pow(px, 2.) + 2 * pow(py, 2.) + pow(pz, 2.))\
						//--------Elastic energy effective field----------//
						+ (tv1 * (eyy0 - eyy) + tv2 * (exx0 - exx + ezz0 - ezz)) * py \
						+ tv3 * ((exy0 - exy) * px + (eyz0 - eyz) * pz));

					//--------------------Z compoennt-------------------//	
					//---------Landau energy effective field------------//
					dqz_crt_local = FE_mass * (\
						2. * a1 * pz + 4. * a11 * pow(pz, 3.) \
						+ 2. * a12 * (pow(py, 2.) + pow(px, 2.)) * pz \
						+ 6. * a111 * pow(pz, 5.) \
						+ 4. * a112 * pow(pz, 3.) * (pow(py, 2.) + pow(px, 2.)) \
						+ 2. * a112 * (pow(py, 4.) + pow(px, 4.)) * pz \
						+ 2. * a123 * pz * pow(py, 2.) * pow(px, 2.) \
						+ 8. * a1111 * pow(pz, 7.) \
						+ 2. * a1112 * pz * (pow(py, 6.) + pow(px, 6.) \
							+ 3. * pow(pz, 4.) * (pow(py, 2.) + pow(px, 2.))) \
						+ 4. * a1122 * pow(pz, 3.) * (pow(py, 4.) + pow(px, 4.)) \
						+ 2. * a1123 * pz * pow(px, 2.) * pow(py, 2.) * \
						(pow(px, 2.) + pow(py, 2.) + 2. * pow(pz, 2.)) \
						//--------Elastic energy effective field----------//
						+ (tv1 * (ezz0 - ezz) + tv2 * (exx0 - exx + eyy0 - eyy)) * pz \
						+ tv3 * ((exz0 - exz) * px + (eyz0 - eyz) * py));

					pt_math->transform_vector_crt2glb(dqx_crt_local, dqy_crt_local, dqz_crt_local, \
						dqx_glb_local, dqy_glb_local, dqz_glb_local);

					dqx_glb_local = dqx_glb_local + FE_mass * (\
						- 1. * G11 * (d2pxdx2_glb + d2pxdy2_glb + d2pxdz2_glb) \
						- Ex_stat(i, j, k) - pt_glb->Eext[0]) \
						+ FE_damping * qx_glb_store(i, j, k);

					dqy_glb_local = dqy_glb_local + FE_mass * (\
						- 1. * G11 * (d2pydx2_glb + d2pydy2_glb + d2pydz2_glb) \
						- Ey_stat(i, j, k) - pt_glb->Eext[1]) \
						+ FE_damping * qy_glb_store(i, j, k);

					dqz_glb_local = dqz_glb_local + FE_mass * (\
						- 1. * G11 * (d2pzdx2_glb + d2pzdy2_glb + d2pzdz2_glb) \
						- Ez_stat(i, j, k) - pt_glb->Eext[2]) \
						+ FE_damping * qz_glb_store(i, j, k);

					if (pt_glb->if_EM_backaction_on_P == true) {
						dqx_glb_local = dqx_glb_local - FE_mass * pt_EM->DEx_em_cell(i, j, k);
						dqy_glb_local = dqy_glb_local - FE_mass * pt_EM->DEy_em_cell(i, j, k);
						dqz_glb_local = dqz_glb_local - FE_mass * pt_EM->DEz_em_cell(i, j, k);
					}

					if (pt_glb->if_flexo == true) {
						//--------------------X compoennt-------------------//
						dqx_glb_local = dqx_glb_local + FE_mass * (\
							+ tv4 * px_glb_local * dpxdx(i, j, k) \
							+ tv5 * (py_glb_local * dpydx(i, j, k) + pz_glb_local * dpzdx(i, j, k)) \
							+ tv6 * (px_glb_local * dpydy(i, j, k) + px_glb_local * dpzdz(i, j, k) + py_glb_local * dpxdy(i, j, k) + pz_glb_local * dpxdz(i, j, k)) \
							- tv7 * d2pxdx2_glb \
							- tv8 * (d2pxdz2_glb + d2pxdy2_glb) \
							//-------------Flexoelectric energy effectiv field--------//
							- f11 * dexxdx - f12 * (deyydx + dezzdx) \
							- 2. * f44 * (dexydy + dexzdz));

						//--------------------Y compoennt-------------------//
						dqy_glb_local = dqy_glb_local + FE_mass * (\
							+ tv4 * py_glb_local * dpydy(i, j, k) \
							+ tv5 * (px_glb_local * dpxdy(i, j, k) + pz_glb_local * dpzdy(i, j, k)) \
							+ tv6 * (py_glb_local * dpxdx(i, j, k) + py_glb_local * dpzdz(i, j, k) + px_glb_local * dpydx(i, j, k) + pz_glb_local * dpydz(i, j, k)) \
							- tv7 * d2pydy2_glb \
							- tv8 * (d2pydx2_glb + d2pydz2_glb) \
							//-------------Flexoelectric energy effectiv field--------//
							- f11 * deyydy - f12 * (dexxdy + dezzdy) \
							- 2. * f44 * (dexydx + deyzdz));

						//--------------------Z compoennt-------------------//
						dqz_glb_local = dqz_glb_local + FE_mass * (\
							+ tv4 * pz_glb_local * dpzdz(i, j, k) \
							+ tv5 * (px_glb_local * dpxdz(i, j, k) + py_glb_local * dpydz(i, j, k)) \
							+ tv6 * (pz_glb_local * dpxdx(i, j, k) + pz_glb_local * dpydy(i, j, k) + px_glb_local * dpzdx(i, j, k) + py_glb_local * dpzdy(i, j, k)) \
							- tv7 * d2pzdz2_glb\
							- tv8 * (d2pzdx2_glb + d2pzdy2_glb) \
							//-------------Flexoelectric energy effectiv field--------//
							- f11 * dezzdz - f12 * (dexxdz + deyydz) \
							- 2. * f44 * (deyzdy + dexzdx));
					}

					dqx_glb_rk1(i, j, k) = dqx_glb_local; //RK
					dqy_glb_rk1(i, j, k) = dqy_glb_local; //RK
					dqz_glb_rk1(i, j, k) = dqz_glb_local; //RK
				}
			}
		}
	}
}

void ferroelectric_system::get_dq_RK2() {
	double dqx_glb_local, dqy_glb_local, dqz_glb_local;
	double dqx_crt_local, dqy_crt_local, dqz_crt_local;
	double px_glb_local, py_glb_local, pz_glb_local;
	double px, py, pz;

	double d2pxdx2_glb, d2pxdy2_glb, d2pxdz2_glb;
	double d2pydx2_glb, d2pydy2_glb, d2pydz2_glb;
	double d2pzdx2_glb, d2pzdy2_glb, d2pzdz2_glb;

	double dexxdx, dexxdy, dexxdz;
	double deyydx, deyydy, deyydz;
	double dezzdx, dezzdy, dezzdz;
	double deyzdx, deyzdy, deyzdz;
	double dexzdx, dexzdy, dexzdz;
	double dexydx, dexydy, dexydz;

	unsigned int mat_type;
	material* mat;
	double exx, eyy, ezz, eyz, exz, exy;
	double exx0, eyy0, ezz0, eyz0, exz0, exy0;

	double a1, \
		a11, a12, \
		a111, a112, a123, \
		a1111, a1112, a1122, a1123;
	double tv1, tv2, tv3, tv4, tv5, tv6, tv7, tv8;
	double f11, f12, f44;
	double G11;

	double FE_mass, FE_damping;
	long i, j, k;

	if (pt_glb->if_prescribe_Eext == false) {
#pragma acc serial default(present) async(5)
		{
			update_external_Efield_half();
		}
	}

#pragma acc wait(1) async(5)
#pragma acc parallel default(present) async(5)
	{
#pragma acc loop gang vector private(\
	dqx_glb_local, dqy_glb_local, dqz_glb_local,\
	dqx_crt_local, dqy_crt_local, dqz_crt_local,\
	px_glb_local, py_glb_local, pz_glb_local,\
	px, py, pz,\
	d2pxdx2_glb, d2pxdy2_glb, d2pxdz2_glb,\
	d2pydx2_glb, d2pydy2_glb, d2pydz2_glb,\
	d2pzdx2_glb, d2pzdy2_glb, d2pzdz2_glb,\
	mat_type,\
	exx, eyy, ezz, eyz, exz, exy,\
	exx0, eyy0, ezz0, eyz0, exz0, exy0,\
	a1, \
	a11, a12, \
	a111, a112, a123, \
	a1111, a1112, a1122, a1123,\
	tv1, tv2, tv3, tv4, tv5, tv6, tv7, tv8,\
	f11, f12, f44,\
	G11,\
	FE_mass, FE_damping,mat,i,j,k,\
	dexxdx, dexxdy, dexxdz,\
	deyydx, deyydy, deyydz,\
	dezzdx, dezzdy, dezzdz,\
	deyzdx, deyzdy, deyzdz,\
	dexzdx, dexzdy, dexzdz,\
	dexydx, dexydy, dexydz)
		for (long int id = 0; id < n; id++) {
			i = id / (ny * nz);
			j = (id - i * (ny * nz)) / nz;
			k = id - i * (ny * nz) - j * nz;

			mat_type = pt_glb->material_cell(i, j, k);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[mat_type - 1]);
				if (mat->if_FE == true) {
					px_glb_local = px_glb_store(i, j, k);
					py_glb_local = py_glb_store(i, j, k);
					pz_glb_local = pz_glb_store(i, j, k);

					get_laplacian_p_local \
						(i, j, k, \
							d2pxdx2_glb, d2pxdy2_glb, d2pxdz2_glb, \
							d2pydx2_glb, d2pydy2_glb, d2pydz2_glb, \
							d2pzdx2_glb, d2pzdy2_glb, d2pzdz2_glb);

					pt_math->transform_vector_glb2crt(px_glb_local, py_glb_local, pz_glb_local, \
						px, py, pz);

					if (pt_glb->if_elasto_on_pORm == true) {
						exx = pt_elas->Dexx_crt(i, j, k);
						eyy = pt_elas->Deyy_crt(i, j, k);
						ezz = pt_elas->Dezz_crt(i, j, k);
						eyz = pt_elas->Deyz_crt(i, j, k);
						exz = pt_elas->Dexz_crt(i, j, k);
						exy = pt_elas->Dexy_crt(i, j, k);

						exx0 = pt_elas->Dexx0_crt(i, j, k);
						eyy0 = pt_elas->Deyy0_crt(i, j, k);
						ezz0 = pt_elas->Dezz0_crt(i, j, k);
						eyz0 = pt_elas->Deyz0_crt(i, j, k);
						exz0 = pt_elas->Dexz0_crt(i, j, k);
						exy0 = pt_elas->Dexy0_crt(i, j, k);
						if (pt_glb->if_elastostatic == true) {
							exx = exx + pt_elas->exxt0_crt(i, j, k);
							eyy = eyy + pt_elas->eyyt0_crt(i, j, k);
							ezz = ezz + pt_elas->ezzt0_crt(i, j, k);
							eyz = eyz + pt_elas->eyzt0_crt(i, j, k);
							exz = exz + pt_elas->exzt0_crt(i, j, k);
							exy = exy + pt_elas->exyt0_crt(i, j, k);

							exx0 = exx0 + pt_elas->exx0t0_crt(i, j, k);
							eyy0 = eyy0 + pt_elas->eyy0t0_crt(i, j, k);
							ezz0 = ezz0 + pt_elas->ezz0t0_crt(i, j, k);
							eyz0 = eyz0 + pt_elas->eyz0t0_crt(i, j, k);
							exz0 = exz0 + pt_elas->exz0t0_crt(i, j, k);
							exy0 = exy0 + pt_elas->exy0t0_crt(i, j, k);
						}
					}
					else {
						exx = 0.; exx0 = 0;
						eyy = 0.; eyy0 = 0;
						ezz = 0.; ezz0 = 0;
						eyz = 0.; eyz0 = 0;
						exz = 0.; exz0 = 0;
						exy = 0.; exy0 = 0;
					}

					if (pt_glb->if_flexo == true) {
						dexxdx = pt_elas->dDexxdx_glb(i, j, k);
						dexxdy = pt_elas->dDexxdy_glb(i, j, k);
						dexxdz = pt_elas->dDexxdz_glb(i, j, k);
						deyydx = pt_elas->dDeyydx_glb(i, j, k);
						deyydy = pt_elas->dDeyydy_glb(i, j, k);
						deyydz = pt_elas->dDeyydz_glb(i, j, k);
						dezzdx = pt_elas->dDezzdx_glb(i, j, k);
						dezzdy = pt_elas->dDezzdy_glb(i, j, k);
						dezzdz = pt_elas->dDezzdz_glb(i, j, k);
						deyzdx = pt_elas->dDeyzdx_glb(i, j, k);
						deyzdy = pt_elas->dDeyzdy_glb(i, j, k);
						deyzdz = pt_elas->dDeyzdz_glb(i, j, k);
						dexzdx = pt_elas->dDexzdx_glb(i, j, k);
						dexzdy = pt_elas->dDexzdy_glb(i, j, k);
						dexzdz = pt_elas->dDexzdz_glb(i, j, k);
						dexydx = pt_elas->dDexydx_glb(i, j, k);
						dexydy = pt_elas->dDexydy_glb(i, j, k);
						dexydz = pt_elas->dDexydz_glb(i, j, k);
						if (pt_glb->if_elastostatic == true) {
							dexxdx = dexxdx + pt_elas->dexxt0dx_glb(i, j, k);
							dexxdy = dexxdy + pt_elas->dexxt0dy_glb(i, j, k);
							dexxdz = dexxdz + pt_elas->dexxt0dz_glb(i, j, k);
							deyydx = deyydx + pt_elas->deyyt0dx_glb(i, j, k);
							deyydy = deyydy + pt_elas->deyyt0dy_glb(i, j, k);
							deyydz = deyydz + pt_elas->deyyt0dz_glb(i, j, k);
							dezzdx = dezzdx + pt_elas->dezzt0dx_glb(i, j, k);
							dezzdy = dezzdy + pt_elas->dezzt0dy_glb(i, j, k);
							dezzdz = dezzdz + pt_elas->dezzt0dz_glb(i, j, k);
							deyzdx = deyzdx + pt_elas->deyzt0dx_glb(i, j, k);
							deyzdy = deyzdy + pt_elas->deyzt0dy_glb(i, j, k);
							deyzdz = deyzdz + pt_elas->deyzt0dz_glb(i, j, k);
							dexzdx = dexzdx + pt_elas->dexzt0dx_glb(i, j, k);
							dexzdy = dexzdy + pt_elas->dexzt0dy_glb(i, j, k);
							dexzdz = dexzdz + pt_elas->dexzt0dz_glb(i, j, k);
							dexydx = dexydx + pt_elas->dexyt0dx_glb(i, j, k);
							dexydy = dexydy + pt_elas->dexyt0dy_glb(i, j, k);
							dexydz = dexydz + pt_elas->dexyt0dz_glb(i, j, k);
						}
					}

					a1 = mat->a1;
					a11 = mat->a11; a12 = mat->a12;
					a111 = mat->a111; a112 = mat->a112; a123 = mat->a123;
					a1111 = mat->a1111; a1112 = mat->a1112; a1122 = mat->a1122; a1123 = mat->a1123;

					tv1 = mat->tv1;	tv2 = mat->tv2; tv3 = mat->tv3; tv4 = mat->tv4;
					tv5 = mat->tv5; tv6 = mat->tv6; tv7 = mat->tv7; tv8 = mat->tv8;

					f11 = mat->f11; f12 = mat->f12; f44 = mat->f44;

					G11 = mat->G11;

					FE_mass = -1. * pt_glb->dt / mat->FE_mass;
					FE_damping = FE_mass * mat->FE_damping;

					//--------------------X compoennt-------------------//
					//---------Landau energy effective field------------//
					dqx_crt_local = FE_mass * (\
						2. * a1 * px + 4. * a11 * pow(px, 3.) \
						+ 2. * a12 * (pow(py, 2.) + pow(pz, 2.)) * px \
						+ 6. * a111 * pow(px, 5.) \
						+ 4. * a112 * pow(px, 3.) * (pow(py, 2.) + pow(pz, 2.)) \
						+ 2. * a112 * (pow(py, 4.) + pow(pz, 4.)) * px \
						+ 2. * a123 * px * pow(py, 2.) * pow(pz, 2.) \
						+ 8. * a1111 * pow(px, 7.) \
						+ 2. * a1112 * px * (pow(py, 6.) + pow(pz, 6.) \
							+ 3. * pow(px, 4.) * (pow(py, 2.) + pow(pz, 2.))) \
						+ 4. * a1122 * pow(px, 3.) * (pow(py, 4.) + pow(pz, 4.)) \
						+ 2. * a1123 * px * pow(py, 2.) * pow(pz, 2.) * \
						(2. * pow(px, 2.) + pow(py, 2.) + pow(pz, 2.)) \
						//--------Elastic energy effective field----------//
						+ (tv1 * (exx0 - exx) + tv2 * (eyy0 - eyy + ezz0 - ezz)) * px \
						+ tv3 * ((exy0 - exy) * py + (exz0 - exz) * pz));

					//--------------------Y compoennt-------------------//	
					//---------Landau energy effective field------------//
					dqy_crt_local = FE_mass * (\
						2. * a1 * py + 4. * a11 * pow(py, 3.) \
						+ 2. * a12 * (pow(pz, 2.) + pow(px, 2.)) * py \
						+ 6. * a111 * pow(py, 5.) \
						+ 4. * a112 * pow(py, 3.) * (pow(px, 2.) + pow(pz, 2.)) \
						+ 2. * a112 * (pow(px, 4.) + pow(pz, 4.)) * py \
						+ 2. * a123 * py * pow(px, 2.) * pow(pz, 2.) \
						+ 8. * a1111 * pow(py, 7) \
						+ 2. * a1112 * py * (pow(px, 6.) + pow(pz, 6.) \
							+ 3. * pow(py, 4.) * (pow(px, 2.) + pow(pz, 2.))) \
						+ 4. * a1122 * pow(py, 3.) * (pow(px, 4.) + pow(pz, 4.)) \
						+ 2. * a1123 * py * pow(px, 2.) * pow(pz, 2.) * \
						(pow(px, 2.) + 2 * pow(py, 2.) + pow(pz, 2.))\
						//--------Elastic energy effective field----------//
						+ (tv1 * (eyy0 - eyy) + tv2 * (exx0 - exx + ezz0 - ezz)) * py \
						+ tv3 * ((exy0 - exy) * px + (eyz0 - eyz) * pz));

					//--------------------Z compoennt-------------------//	
					//---------Landau energy effective field------------//
					dqz_crt_local = FE_mass * (\
						2. * a1 * pz + 4. * a11 * pow(pz, 3.) \
						+ 2. * a12 * (pow(py, 2.) + pow(px, 2.)) * pz \
						+ 6. * a111 * pow(pz, 5.) \
						+ 4. * a112 * pow(pz, 3.) * (pow(py, 2.) + pow(px, 2.)) \
						+ 2. * a112 * (pow(py, 4.) + pow(px, 4.)) * pz \
						+ 2. * a123 * pz * pow(py, 2.) * pow(px, 2.) \
						+ 8. * a1111 * pow(pz, 7.) \
						+ 2. * a1112 * pz * (pow(py, 6.) + pow(px, 6.) \
							+ 3. * pow(pz, 4.) * (pow(py, 2.) + pow(px, 2.))) \
						+ 4. * a1122 * pow(pz, 3.) * (pow(py, 4.) + pow(px, 4.)) \
						+ 2. * a1123 * pz * pow(px, 2.) * pow(py, 2.) * \
						(pow(px, 2.) + pow(py, 2.) + 2. * pow(pz, 2.)) \
						//--------Elastic energy effective field----------//
						+ (tv1 * (ezz0 - ezz) + tv2 * (exx0 - exx + eyy0 - eyy)) * pz \
						+ tv3 * ((exz0 - exz) * px + (eyz0 - eyz) * py));

					pt_math->transform_vector_crt2glb(dqx_crt_local, dqy_crt_local, dqz_crt_local, \
						dqx_glb_local, dqy_glb_local, dqz_glb_local);

					dqx_glb_local = dqx_glb_local + FE_mass * (\
						- 1. * G11 * (d2pxdx2_glb + d2pxdy2_glb + d2pxdz2_glb) \
						- Ex_stat(i, j, k) - pt_glb->Eext[0]) \
						+ FE_damping * qx_glb_store(i, j, k);

					dqy_glb_local = dqy_glb_local + FE_mass * (\
						- 1. * G11 * (d2pydx2_glb + d2pydy2_glb + d2pydz2_glb) \
						- Ey_stat(i, j, k) - pt_glb->Eext[1]) \
						+ FE_damping * qy_glb_store(i, j, k);

					dqz_glb_local = dqz_glb_local + FE_mass * (\
						- 1. * G11 * (d2pzdx2_glb + d2pzdy2_glb + d2pzdz2_glb) \
						- Ez_stat(i, j, k) - pt_glb->Eext[2]) \
						+ FE_damping * qz_glb_store(i, j, k);

					if (pt_glb->if_EM_backaction_on_P == true) {
						dqx_glb_local = dqx_glb_local - FE_mass * pt_EM->DEx_em_cell(i, j, k);
						dqy_glb_local = dqy_glb_local - FE_mass * pt_EM->DEy_em_cell(i, j, k);
						dqz_glb_local = dqz_glb_local - FE_mass * pt_EM->DEz_em_cell(i, j, k);
					}

					if (pt_glb->if_flexo == true) {
						//--------------------X compoennt-------------------//
						dqx_glb_local = dqx_glb_local + FE_mass * (\
							+ tv4 * px_glb_local * dpxdx(i, j, k) \
							+ tv5 * (py_glb_local * dpydx(i, j, k) + pz_glb_local * dpzdx(i, j, k)) \
							+ tv6 * (px_glb_local * dpydy(i, j, k) + px_glb_local * dpzdz(i, j, k) + py_glb_local * dpxdy(i, j, k) + pz_glb_local * dpxdz(i, j, k)) \
							- tv7 * d2pxdx2_glb \
							- tv8 * (d2pxdz2_glb + d2pxdy2_glb) \
							//-------------Flexoelectric energy effectiv field--------//
							- f11 * dexxdx - f12 * (deyydx + dezzdx) \
							- 2. * f44 * (dexydy + dexzdz));

						//--------------------Y compoennt-------------------//
						dqy_glb_local = dqy_glb_local + FE_mass * (\
							+ tv4 * py_glb_local * dpydy(i, j, k) \
							+ tv5 * (px_glb_local * dpxdy(i, j, k) + pz_glb_local * dpzdy(i, j, k)) \
							+ tv6 * (py_glb_local * dpxdx(i, j, k) + py_glb_local * dpzdz(i, j, k) + px_glb_local * dpydx(i, j, k) + pz_glb_local * dpydz(i, j, k)) \
							- tv7 * d2pydy2_glb \
							- tv8 * (d2pydx2_glb + d2pydz2_glb) \
							//-------------Flexoelectric energy effectiv field--------//
							- f11 * deyydy - f12 * (dexxdy + dezzdy) \
							- 2. * f44 * (dexydx + deyzdz));

						//--------------------Z compoennt-------------------//
						dqz_glb_local = dqz_glb_local + FE_mass * (\
							+ tv4 * pz_glb_local * dpzdz(i, j, k) \
							+ tv5 * (px_glb_local * dpxdz(i, j, k) + py_glb_local * dpydz(i, j, k)) \
							+ tv6 * (pz_glb_local * dpxdx(i, j, k) + pz_glb_local * dpydy(i, j, k) + px_glb_local * dpzdx(i, j, k) + py_glb_local * dpzdy(i, j, k)) \
							- tv7 * d2pzdz2_glb\
							- tv8 * (d2pzdx2_glb + d2pzdy2_glb) \
							//-------------Flexoelectric energy effectiv field--------//
							- f11 * dezzdz - f12 * (dexxdz + deyydz) \
							- 2. * f44 * (deyzdy + dexzdx));
					}

					dqx_glb_rk2(i, j, k) = dqx_glb_local; //RK
					dqy_glb_rk2(i, j, k) = dqy_glb_local; //RK
					dqz_glb_rk2(i, j, k) = dqz_glb_local; //RK
				}
			}
		}
	}
}

void ferroelectric_system::get_dq_RK3() {
	double dqx_glb_local, dqy_glb_local, dqz_glb_local;
	double dqx_crt_local, dqy_crt_local, dqz_crt_local;
	double px_glb_local, py_glb_local, pz_glb_local;
	double px, py, pz;

	double d2pxdx2_glb, d2pxdy2_glb, d2pxdz2_glb;
	double d2pydx2_glb, d2pydy2_glb, d2pydz2_glb;
	double d2pzdx2_glb, d2pzdy2_glb, d2pzdz2_glb;

	double dexxdx, dexxdy, dexxdz;
	double deyydx, deyydy, deyydz;
	double dezzdx, dezzdy, dezzdz;
	double deyzdx, deyzdy, deyzdz;
	double dexzdx, dexzdy, dexzdz;
	double dexydx, dexydy, dexydz;

	unsigned int mat_type;
	material* mat;
	double exx, eyy, ezz, eyz, exz, exy;
	double exx0, eyy0, ezz0, eyz0, exz0, exy0;

	double a1, \
		a11, a12, \
		a111, a112, a123, \
		a1111, a1112, a1122, a1123;
	double tv1, tv2, tv3, tv4, tv5, tv6, tv7, tv8;
	double f11, f12, f44;
	double G11;

	double FE_mass, FE_damping;
	long i, j, k;

	if (pt_glb->if_prescribe_Eext == false) {
#pragma acc serial default(present) async(5)
		{
			update_external_Efield_half();
		}
	}

#pragma acc wait(1) async(5)
#pragma acc parallel default(present) async(5)
	{
#pragma acc loop gang vector private(\
	dqx_glb_local, dqy_glb_local, dqz_glb_local,\
	dqx_crt_local, dqy_crt_local, dqz_crt_local,\
	px_glb_local, py_glb_local, pz_glb_local,\
	px, py, pz,\
	d2pxdx2_glb, d2pxdy2_glb, d2pxdz2_glb,\
	d2pydx2_glb, d2pydy2_glb, d2pydz2_glb,\
	d2pzdx2_glb, d2pzdy2_glb, d2pzdz2_glb,\
	mat_type,\
	exx, eyy, ezz, eyz, exz, exy,\
	exx0, eyy0, ezz0, eyz0, exz0, exy0,\
	a1, \
	a11, a12, \
	a111, a112, a123, \
	a1111, a1112, a1122, a1123,\
	tv1, tv2, tv3, tv4, tv5, tv6, tv7, tv8,\
	f11, f12, f44,\
	G11,\
	FE_mass, FE_damping,mat,i,j,k,\
	dexxdx, dexxdy, dexxdz,\
	deyydx, deyydy, deyydz,\
	dezzdx, dezzdy, dezzdz,\
	deyzdx, deyzdy, deyzdz,\
	dexzdx, dexzdy, dexzdz,\
	dexydx, dexydy, dexydz)
		for (long int id = 0; id < n; id++) {
			i = id / (ny * nz);
			j = (id - i * (ny * nz)) / nz;
			k = id - i * (ny * nz) - j * nz;

			mat_type = pt_glb->material_cell(i, j, k);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[mat_type - 1]);
				if (mat->if_FE == true) {
					px_glb_local = px_glb_store(i, j, k);
					py_glb_local = py_glb_store(i, j, k);
					pz_glb_local = pz_glb_store(i, j, k);

					get_laplacian_p_local \
						(i, j, k, \
							d2pxdx2_glb, d2pxdy2_glb, d2pxdz2_glb, \
							d2pydx2_glb, d2pydy2_glb, d2pydz2_glb, \
							d2pzdx2_glb, d2pzdy2_glb, d2pzdz2_glb);

					pt_math->transform_vector_glb2crt(px_glb_local, py_glb_local, pz_glb_local, \
						px, py, pz);

					if (pt_glb->if_elasto_on_pORm == true) {
						exx = pt_elas->Dexx_crt(i, j, k);
						eyy = pt_elas->Deyy_crt(i, j, k);
						ezz = pt_elas->Dezz_crt(i, j, k);
						eyz = pt_elas->Deyz_crt(i, j, k);
						exz = pt_elas->Dexz_crt(i, j, k);
						exy = pt_elas->Dexy_crt(i, j, k);

						exx0 = pt_elas->Dexx0_crt(i, j, k);
						eyy0 = pt_elas->Deyy0_crt(i, j, k);
						ezz0 = pt_elas->Dezz0_crt(i, j, k);
						eyz0 = pt_elas->Deyz0_crt(i, j, k);
						exz0 = pt_elas->Dexz0_crt(i, j, k);
						exy0 = pt_elas->Dexy0_crt(i, j, k);
						if (pt_glb->if_elastostatic == true) {
							exx = exx + pt_elas->exxt0_crt(i, j, k);
							eyy = eyy + pt_elas->eyyt0_crt(i, j, k);
							ezz = ezz + pt_elas->ezzt0_crt(i, j, k);
							eyz = eyz + pt_elas->eyzt0_crt(i, j, k);
							exz = exz + pt_elas->exzt0_crt(i, j, k);
							exy = exy + pt_elas->exyt0_crt(i, j, k);

							exx0 = exx0 + pt_elas->exx0t0_crt(i, j, k);
							eyy0 = eyy0 + pt_elas->eyy0t0_crt(i, j, k);
							ezz0 = ezz0 + pt_elas->ezz0t0_crt(i, j, k);
							eyz0 = eyz0 + pt_elas->eyz0t0_crt(i, j, k);
							exz0 = exz0 + pt_elas->exz0t0_crt(i, j, k);
							exy0 = exy0 + pt_elas->exy0t0_crt(i, j, k);
						}
					}
					else {
						exx = 0.; exx0 = 0;
						eyy = 0.; eyy0 = 0;
						ezz = 0.; ezz0 = 0;
						eyz = 0.; eyz0 = 0;
						exz = 0.; exz0 = 0;
						exy = 0.; exy0 = 0;
					}

					if (pt_glb->if_flexo == true) {
						dexxdx = pt_elas->dDexxdx_glb(i, j, k);
						dexxdy = pt_elas->dDexxdy_glb(i, j, k);
						dexxdz = pt_elas->dDexxdz_glb(i, j, k);
						deyydx = pt_elas->dDeyydx_glb(i, j, k);
						deyydy = pt_elas->dDeyydy_glb(i, j, k);
						deyydz = pt_elas->dDeyydz_glb(i, j, k);
						dezzdx = pt_elas->dDezzdx_glb(i, j, k);
						dezzdy = pt_elas->dDezzdy_glb(i, j, k);
						dezzdz = pt_elas->dDezzdz_glb(i, j, k);
						deyzdx = pt_elas->dDeyzdx_glb(i, j, k);
						deyzdy = pt_elas->dDeyzdy_glb(i, j, k);
						deyzdz = pt_elas->dDeyzdz_glb(i, j, k);
						dexzdx = pt_elas->dDexzdx_glb(i, j, k);
						dexzdy = pt_elas->dDexzdy_glb(i, j, k);
						dexzdz = pt_elas->dDexzdz_glb(i, j, k);
						dexydx = pt_elas->dDexydx_glb(i, j, k);
						dexydy = pt_elas->dDexydy_glb(i, j, k);
						dexydz = pt_elas->dDexydz_glb(i, j, k);
						if (pt_glb->if_elastostatic == true) {
							dexxdx = dexxdx + pt_elas->dexxt0dx_glb(i, j, k);
							dexxdy = dexxdy + pt_elas->dexxt0dy_glb(i, j, k);
							dexxdz = dexxdz + pt_elas->dexxt0dz_glb(i, j, k);
							deyydx = deyydx + pt_elas->deyyt0dx_glb(i, j, k);
							deyydy = deyydy + pt_elas->deyyt0dy_glb(i, j, k);
							deyydz = deyydz + pt_elas->deyyt0dz_glb(i, j, k);
							dezzdx = dezzdx + pt_elas->dezzt0dx_glb(i, j, k);
							dezzdy = dezzdy + pt_elas->dezzt0dy_glb(i, j, k);
							dezzdz = dezzdz + pt_elas->dezzt0dz_glb(i, j, k);
							deyzdx = deyzdx + pt_elas->deyzt0dx_glb(i, j, k);
							deyzdy = deyzdy + pt_elas->deyzt0dy_glb(i, j, k);
							deyzdz = deyzdz + pt_elas->deyzt0dz_glb(i, j, k);
							dexzdx = dexzdx + pt_elas->dexzt0dx_glb(i, j, k);
							dexzdy = dexzdy + pt_elas->dexzt0dy_glb(i, j, k);
							dexzdz = dexzdz + pt_elas->dexzt0dz_glb(i, j, k);
							dexydx = dexydx + pt_elas->dexyt0dx_glb(i, j, k);
							dexydy = dexydy + pt_elas->dexyt0dy_glb(i, j, k);
							dexydz = dexydz + pt_elas->dexyt0dz_glb(i, j, k);
						}
					}

					a1 = mat->a1;
					a11 = mat->a11; a12 = mat->a12;
					a111 = mat->a111; a112 = mat->a112; a123 = mat->a123;
					a1111 = mat->a1111; a1112 = mat->a1112; a1122 = mat->a1122; a1123 = mat->a1123;

					tv1 = mat->tv1;	tv2 = mat->tv2; tv3 = mat->tv3; tv4 = mat->tv4;
					tv5 = mat->tv5; tv6 = mat->tv6; tv7 = mat->tv7; tv8 = mat->tv8;

					f11 = mat->f11; f12 = mat->f12; f44 = mat->f44;

					G11 = mat->G11;

					FE_mass = -1. * pt_glb->dt / mat->FE_mass;
					FE_damping = FE_mass * mat->FE_damping;

					//--------------------X compoennt-------------------//
					//---------Landau energy effective field------------//
					dqx_crt_local = FE_mass * (\
						2. * a1 * px + 4. * a11 * pow(px, 3.) \
						+ 2. * a12 * (pow(py, 2.) + pow(pz, 2.)) * px \
						+ 6. * a111 * pow(px, 5.) \
						+ 4. * a112 * pow(px, 3.) * (pow(py, 2.) + pow(pz, 2.)) \
						+ 2. * a112 * (pow(py, 4.) + pow(pz, 4.)) * px \
						+ 2. * a123 * px * pow(py, 2.) * pow(pz, 2.) \
						+ 8. * a1111 * pow(px, 7.) \
						+ 2. * a1112 * px * (pow(py, 6.) + pow(pz, 6.) \
							+ 3. * pow(px, 4.) * (pow(py, 2.) + pow(pz, 2.))) \
						+ 4. * a1122 * pow(px, 3.) * (pow(py, 4.) + pow(pz, 4.)) \
						+ 2. * a1123 * px * pow(py, 2.) * pow(pz, 2.) * \
						(2. * pow(px, 2.) + pow(py, 2.) + pow(pz, 2.)) \
						//--------Elastic energy effective field----------//
						+ (tv1 * (exx0 - exx) + tv2 * (eyy0 - eyy + ezz0 - ezz)) * px \
						+ tv3 * ((exy0 - exy) * py + (exz0 - exz) * pz));

					//--------------------Y compoennt-------------------//	
					//---------Landau energy effective field------------//
					dqy_crt_local = FE_mass * (\
						2. * a1 * py + 4. * a11 * pow(py, 3.) \
						+ 2. * a12 * (pow(pz, 2.) + pow(px, 2.)) * py \
						+ 6. * a111 * pow(py, 5.) \
						+ 4. * a112 * pow(py, 3.) * (pow(px, 2.) + pow(pz, 2.)) \
						+ 2. * a112 * (pow(px, 4.) + pow(pz, 4.)) * py \
						+ 2. * a123 * py * pow(px, 2.) * pow(pz, 2.) \
						+ 8. * a1111 * pow(py, 7) \
						+ 2. * a1112 * py * (pow(px, 6.) + pow(pz, 6.) \
							+ 3. * pow(py, 4.) * (pow(px, 2.) + pow(pz, 2.))) \
						+ 4. * a1122 * pow(py, 3.) * (pow(px, 4.) + pow(pz, 4.)) \
						+ 2. * a1123 * py * pow(px, 2.) * pow(pz, 2.) * \
						(pow(px, 2.) + 2 * pow(py, 2.) + pow(pz, 2.))\
						//--------Elastic energy effective field----------//
						+ (tv1 * (eyy0 - eyy) + tv2 * (exx0 - exx + ezz0 - ezz)) * py \
						+ tv3 * ((exy0 - exy) * px + (eyz0 - eyz) * pz));

					//--------------------Z compoennt-------------------//	
					//---------Landau energy effective field------------//
					dqz_crt_local = FE_mass * (\
						2. * a1 * pz + 4. * a11 * pow(pz, 3.) \
						+ 2. * a12 * (pow(py, 2.) + pow(px, 2.)) * pz \
						+ 6. * a111 * pow(pz, 5.) \
						+ 4. * a112 * pow(pz, 3.) * (pow(py, 2.) + pow(px, 2.)) \
						+ 2. * a112 * (pow(py, 4.) + pow(px, 4.)) * pz \
						+ 2. * a123 * pz * pow(py, 2.) * pow(px, 2.) \
						+ 8. * a1111 * pow(pz, 7.) \
						+ 2. * a1112 * pz * (pow(py, 6.) + pow(px, 6.) \
							+ 3. * pow(pz, 4.) * (pow(py, 2.) + pow(px, 2.))) \
						+ 4. * a1122 * pow(pz, 3.) * (pow(py, 4.) + pow(px, 4.)) \
						+ 2. * a1123 * pz * pow(px, 2.) * pow(py, 2.) * \
						(pow(px, 2.) + pow(py, 2.) + 2. * pow(pz, 2.)) \
						//--------Elastic energy effective field----------//
						+ (tv1 * (ezz0 - ezz) + tv2 * (exx0 - exx + eyy0 - eyy)) * pz \
						+ tv3 * ((exz0 - exz) * px + (eyz0 - eyz) * py));

					pt_math->transform_vector_crt2glb(dqx_crt_local, dqy_crt_local, dqz_crt_local, \
						dqx_glb_local, dqy_glb_local, dqz_glb_local);

					dqx_glb_local = dqx_glb_local + FE_mass * (\
						- 1. * G11 * (d2pxdx2_glb + d2pxdy2_glb + d2pxdz2_glb) \
						- Ex_stat(i, j, k) - pt_glb->Eext[0]) \
						+ FE_damping * qx_glb_store(i, j, k);

					dqy_glb_local = dqy_glb_local + FE_mass * (\
						- 1. * G11 * (d2pydx2_glb + d2pydy2_glb + d2pydz2_glb) \
						- Ey_stat(i, j, k) - pt_glb->Eext[1]) \
						+ FE_damping * qy_glb_store(i, j, k);

					dqz_glb_local = dqz_glb_local + FE_mass * (\
						- 1. * G11 * (d2pzdx2_glb + d2pzdy2_glb + d2pzdz2_glb) \
						- Ez_stat(i, j, k) - pt_glb->Eext[2]) \
						+ FE_damping * qz_glb_store(i, j, k);

					if (pt_glb->if_EM_backaction_on_P == true) {
						dqx_glb_local = dqx_glb_local - FE_mass * pt_EM->DEx_em_cell(i, j, k);
						dqy_glb_local = dqy_glb_local - FE_mass * pt_EM->DEy_em_cell(i, j, k);
						dqz_glb_local = dqz_glb_local - FE_mass * pt_EM->DEz_em_cell(i, j, k);
					}

					if (pt_glb->if_flexo == true) {
						//--------------------X compoennt-------------------//
						dqx_glb_local = dqx_glb_local + FE_mass * (\
							+ tv4 * px_glb_local * dpxdx(i, j, k) \
							+ tv5 * (py_glb_local * dpydx(i, j, k) + pz_glb_local * dpzdx(i, j, k)) \
							+ tv6 * (px_glb_local * dpydy(i, j, k) + px_glb_local * dpzdz(i, j, k) + py_glb_local * dpxdy(i, j, k) + pz_glb_local * dpxdz(i, j, k)) \
							- tv7 * d2pxdx2_glb \
							- tv8 * (d2pxdz2_glb + d2pxdy2_glb) \
							//-------------Flexoelectric energy effectiv field--------//
							- f11 * dexxdx - f12 * (deyydx + dezzdx) \
							- 2. * f44 * (dexydy + dexzdz));

						//--------------------Y compoennt-------------------//
						dqy_glb_local = dqy_glb_local + FE_mass * (\
							+ tv4 * py_glb_local * dpydy(i, j, k) \
							+ tv5 * (px_glb_local * dpxdy(i, j, k) + pz_glb_local * dpzdy(i, j, k)) \
							+ tv6 * (py_glb_local * dpxdx(i, j, k) + py_glb_local * dpzdz(i, j, k) + px_glb_local * dpydx(i, j, k) + pz_glb_local * dpydz(i, j, k)) \
							- tv7 * d2pydy2_glb \
							- tv8 * (d2pydx2_glb + d2pydz2_glb) \
							//-------------Flexoelectric energy effectiv field--------//
							- f11 * deyydy - f12 * (dexxdy + dezzdy) \
							- 2. * f44 * (dexydx + deyzdz));

						//--------------------Z compoennt-------------------//
						dqz_glb_local = dqz_glb_local + FE_mass * (\
							+ tv4 * pz_glb_local * dpzdz(i, j, k) \
							+ tv5 * (px_glb_local * dpxdz(i, j, k) + py_glb_local * dpydz(i, j, k)) \
							+ tv6 * (pz_glb_local * dpxdx(i, j, k) + pz_glb_local * dpydy(i, j, k) + px_glb_local * dpzdx(i, j, k) + py_glb_local * dpzdy(i, j, k)) \
							- tv7 * d2pzdz2_glb\
							- tv8 * (d2pzdx2_glb + d2pzdy2_glb) \
							//-------------Flexoelectric energy effectiv field--------//
							- f11 * dezzdz - f12 * (dexxdz + deyydz) \
							- 2. * f44 * (deyzdy + dexzdx));
					}

					dqx_glb_rk3(i, j, k) = dqx_glb_local; //RK
					dqy_glb_rk3(i, j, k) = dqy_glb_local; //RK
					dqz_glb_rk3(i, j, k) = dqz_glb_local; //RK
				}
			}
		}
	}
}

void ferroelectric_system::get_dq_RK4() {
	double dqx_glb_local, dqy_glb_local, dqz_glb_local;
	double dqx_crt_local, dqy_crt_local, dqz_crt_local;
	double px_glb_local, py_glb_local, pz_glb_local;
	double px, py, pz;

	double d2pxdx2_glb, d2pxdy2_glb, d2pxdz2_glb;
	double d2pydx2_glb, d2pydy2_glb, d2pydz2_glb;
	double d2pzdx2_glb, d2pzdy2_glb, d2pzdz2_glb;

	double dexxdx, dexxdy, dexxdz;
	double deyydx, deyydy, deyydz;
	double dezzdx, dezzdy, dezzdz;
	double deyzdx, deyzdy, deyzdz;
	double dexzdx, dexzdy, dexzdz;
	double dexydx, dexydy, dexydz;

	unsigned int mat_type;
	material* mat;
	double exx, eyy, ezz, eyz, exz, exy;
	double exx0, eyy0, ezz0, eyz0, exz0, exy0;

	double a1, \
		a11, a12, \
		a111, a112, a123, \
		a1111, a1112, a1122, a1123;
	double tv1, tv2, tv3, tv4, tv5, tv6, tv7, tv8;
	double f11, f12, f44;
	double G11;

	double FE_mass, FE_damping;
	long i, j, k;

	if (pt_glb->if_prescribe_Eext == false) {
#pragma acc serial default(present) async(5)
		{
			update_external_Efield_full();
		}
	}

#pragma acc wait(1) async(5)
#pragma acc parallel default(present) async(5)
	{
#pragma acc loop gang vector private(\
	dqx_glb_local, dqy_glb_local, dqz_glb_local,\
	dqx_crt_local, dqy_crt_local, dqz_crt_local,\
	px_glb_local, py_glb_local, pz_glb_local,\
	px, py, pz,\
	d2pxdx2_glb, d2pxdy2_glb, d2pxdz2_glb,\
	d2pydx2_glb, d2pydy2_glb, d2pydz2_glb,\
	d2pzdx2_glb, d2pzdy2_glb, d2pzdz2_glb,\
	mat_type,\
	exx, eyy, ezz, eyz, exz, exy,\
	exx0, eyy0, ezz0, eyz0, exz0, exy0,\
	a1, \
	a11, a12, \
	a111, a112, a123, \
	a1111, a1112, a1122, a1123,\
	tv1, tv2, tv3, tv4, tv5, tv6, tv7, tv8,\
	f11, f12, f44,\
	G11,\
	FE_mass, FE_damping,mat,i,j,k,\
	dexxdx, dexxdy, dexxdz,\
	deyydx, deyydy, deyydz,\
	dezzdx, dezzdy, dezzdz,\
	deyzdx, deyzdy, deyzdz,\
	dexzdx, dexzdy, dexzdz,\
	dexydx, dexydy, dexydz)
		for (long int id = 0; id < n; id++) {
			i = id / (ny * nz);
			j = (id - i * (ny * nz)) / nz;
			k = id - i * (ny * nz) - j * nz;

			mat_type = pt_glb->material_cell(i, j, k);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[mat_type - 1]);
				if (mat->if_FE == true) {
					px_glb_local = px_glb_store(i, j, k);
					py_glb_local = py_glb_store(i, j, k);
					pz_glb_local = pz_glb_store(i, j, k);

					get_laplacian_p_local \
						(i, j, k, \
							d2pxdx2_glb, d2pxdy2_glb, d2pxdz2_glb, \
							d2pydx2_glb, d2pydy2_glb, d2pydz2_glb, \
							d2pzdx2_glb, d2pzdy2_glb, d2pzdz2_glb);

					pt_math->transform_vector_glb2crt(px_glb_local, py_glb_local, pz_glb_local, \
						px, py, pz);

					if (pt_glb->if_elasto_on_pORm == true) {
						exx = pt_elas->Dexx_crt(i, j, k);
						eyy = pt_elas->Deyy_crt(i, j, k);
						ezz = pt_elas->Dezz_crt(i, j, k);
						eyz = pt_elas->Deyz_crt(i, j, k);
						exz = pt_elas->Dexz_crt(i, j, k);
						exy = pt_elas->Dexy_crt(i, j, k);

						exx0 = pt_elas->Dexx0_crt(i, j, k);
						eyy0 = pt_elas->Deyy0_crt(i, j, k);
						ezz0 = pt_elas->Dezz0_crt(i, j, k);
						eyz0 = pt_elas->Deyz0_crt(i, j, k);
						exz0 = pt_elas->Dexz0_crt(i, j, k);
						exy0 = pt_elas->Dexy0_crt(i, j, k);
						if (pt_glb->if_elastostatic == true) {
							exx = exx + pt_elas->exxt0_crt(i, j, k);
							eyy = eyy + pt_elas->eyyt0_crt(i, j, k);
							ezz = ezz + pt_elas->ezzt0_crt(i, j, k);
							eyz = eyz + pt_elas->eyzt0_crt(i, j, k);
							exz = exz + pt_elas->exzt0_crt(i, j, k);
							exy = exy + pt_elas->exyt0_crt(i, j, k);

							exx0 = exx0 + pt_elas->exx0t0_crt(i, j, k);
							eyy0 = eyy0 + pt_elas->eyy0t0_crt(i, j, k);
							ezz0 = ezz0 + pt_elas->ezz0t0_crt(i, j, k);
							eyz0 = eyz0 + pt_elas->eyz0t0_crt(i, j, k);
							exz0 = exz0 + pt_elas->exz0t0_crt(i, j, k);
							exy0 = exy0 + pt_elas->exy0t0_crt(i, j, k);
						}
					}
					else {
						exx = 0.; exx0 = 0;
						eyy = 0.; eyy0 = 0;
						ezz = 0.; ezz0 = 0;
						eyz = 0.; eyz0 = 0;
						exz = 0.; exz0 = 0;
						exy = 0.; exy0 = 0;
					}

					if (pt_glb->if_flexo == true) {
						dexxdx = pt_elas->dDexxdx_glb(i, j, k);
						dexxdy = pt_elas->dDexxdy_glb(i, j, k);
						dexxdz = pt_elas->dDexxdz_glb(i, j, k);
						deyydx = pt_elas->dDeyydx_glb(i, j, k);
						deyydy = pt_elas->dDeyydy_glb(i, j, k);
						deyydz = pt_elas->dDeyydz_glb(i, j, k);
						dezzdx = pt_elas->dDezzdx_glb(i, j, k);
						dezzdy = pt_elas->dDezzdy_glb(i, j, k);
						dezzdz = pt_elas->dDezzdz_glb(i, j, k);
						deyzdx = pt_elas->dDeyzdx_glb(i, j, k);
						deyzdy = pt_elas->dDeyzdy_glb(i, j, k);
						deyzdz = pt_elas->dDeyzdz_glb(i, j, k);
						dexzdx = pt_elas->dDexzdx_glb(i, j, k);
						dexzdy = pt_elas->dDexzdy_glb(i, j, k);
						dexzdz = pt_elas->dDexzdz_glb(i, j, k);
						dexydx = pt_elas->dDexydx_glb(i, j, k);
						dexydy = pt_elas->dDexydy_glb(i, j, k);
						dexydz = pt_elas->dDexydz_glb(i, j, k);
						if (pt_glb->if_elastostatic == true) {
							dexxdx = dexxdx + pt_elas->dexxt0dx_glb(i, j, k);
							dexxdy = dexxdy + pt_elas->dexxt0dy_glb(i, j, k);
							dexxdz = dexxdz + pt_elas->dexxt0dz_glb(i, j, k);
							deyydx = deyydx + pt_elas->deyyt0dx_glb(i, j, k);
							deyydy = deyydy + pt_elas->deyyt0dy_glb(i, j, k);
							deyydz = deyydz + pt_elas->deyyt0dz_glb(i, j, k);
							dezzdx = dezzdx + pt_elas->dezzt0dx_glb(i, j, k);
							dezzdy = dezzdy + pt_elas->dezzt0dy_glb(i, j, k);
							dezzdz = dezzdz + pt_elas->dezzt0dz_glb(i, j, k);
							deyzdx = deyzdx + pt_elas->deyzt0dx_glb(i, j, k);
							deyzdy = deyzdy + pt_elas->deyzt0dy_glb(i, j, k);
							deyzdz = deyzdz + pt_elas->deyzt0dz_glb(i, j, k);
							dexzdx = dexzdx + pt_elas->dexzt0dx_glb(i, j, k);
							dexzdy = dexzdy + pt_elas->dexzt0dy_glb(i, j, k);
							dexzdz = dexzdz + pt_elas->dexzt0dz_glb(i, j, k);
							dexydx = dexydx + pt_elas->dexyt0dx_glb(i, j, k);
							dexydy = dexydy + pt_elas->dexyt0dy_glb(i, j, k);
							dexydz = dexydz + pt_elas->dexyt0dz_glb(i, j, k);
						}
					}

					a1 = mat->a1;
					a11 = mat->a11; a12 = mat->a12;
					a111 = mat->a111; a112 = mat->a112; a123 = mat->a123;
					a1111 = mat->a1111; a1112 = mat->a1112; a1122 = mat->a1122; a1123 = mat->a1123;

					tv1 = mat->tv1;	tv2 = mat->tv2; tv3 = mat->tv3; tv4 = mat->tv4;
					tv5 = mat->tv5; tv6 = mat->tv6; tv7 = mat->tv7; tv8 = mat->tv8;

					f11 = mat->f11; f12 = mat->f12; f44 = mat->f44;

					G11 = mat->G11;

					FE_mass = -1. * pt_glb->dt / mat->FE_mass;
					FE_damping = FE_mass * mat->FE_damping;

					//--------------------X compoennt-------------------//
					//---------Landau energy effective field------------//
					dqx_crt_local = FE_mass * (\
						2. * a1 * px + 4. * a11 * pow(px, 3.) \
						+ 2. * a12 * (pow(py, 2.) + pow(pz, 2.)) * px \
						+ 6. * a111 * pow(px, 5.) \
						+ 4. * a112 * pow(px, 3.) * (pow(py, 2.) + pow(pz, 2.)) \
						+ 2. * a112 * (pow(py, 4.) + pow(pz, 4.)) * px \
						+ 2. * a123 * px * pow(py, 2.) * pow(pz, 2.) \
						+ 8. * a1111 * pow(px, 7.) \
						+ 2. * a1112 * px * (pow(py, 6.) + pow(pz, 6.) \
							+ 3. * pow(px, 4.) * (pow(py, 2.) + pow(pz, 2.))) \
						+ 4. * a1122 * pow(px, 3.) * (pow(py, 4.) + pow(pz, 4.)) \
						+ 2. * a1123 * px * pow(py, 2.) * pow(pz, 2.) * \
						(2. * pow(px, 2.) + pow(py, 2.) + pow(pz, 2.)) \
						//--------Elastic energy effective field----------//
						+ (tv1 * (exx0 - exx) + tv2 * (eyy0 - eyy + ezz0 - ezz)) * px \
						+ tv3 * ((exy0 - exy) * py + (exz0 - exz) * pz));

					//--------------------Y compoennt-------------------//	
					//---------Landau energy effective field------------//
					dqy_crt_local = FE_mass * (\
						2. * a1 * py + 4. * a11 * pow(py, 3.) \
						+ 2. * a12 * (pow(pz, 2.) + pow(px, 2.)) * py \
						+ 6. * a111 * pow(py, 5.) \
						+ 4. * a112 * pow(py, 3.) * (pow(px, 2.) + pow(pz, 2.)) \
						+ 2. * a112 * (pow(px, 4.) + pow(pz, 4.)) * py \
						+ 2. * a123 * py * pow(px, 2.) * pow(pz, 2.) \
						+ 8. * a1111 * pow(py, 7) \
						+ 2. * a1112 * py * (pow(px, 6.) + pow(pz, 6.) \
							+ 3. * pow(py, 4.) * (pow(px, 2.) + pow(pz, 2.))) \
						+ 4. * a1122 * pow(py, 3.) * (pow(px, 4.) + pow(pz, 4.)) \
						+ 2. * a1123 * py * pow(px, 2.) * pow(pz, 2.) * \
						(pow(px, 2.) + 2 * pow(py, 2.) + pow(pz, 2.))\
						//--------Elastic energy effective field----------//
						+ (tv1 * (eyy0 - eyy) + tv2 * (exx0 - exx + ezz0 - ezz)) * py \
						+ tv3 * ((exy0 - exy) * px + (eyz0 - eyz) * pz));

					//--------------------Z compoennt-------------------//	
					//---------Landau energy effective field------------//
					dqz_crt_local = FE_mass * (\
						2. * a1 * pz + 4. * a11 * pow(pz, 3.) \
						+ 2. * a12 * (pow(py, 2.) + pow(px, 2.)) * pz \
						+ 6. * a111 * pow(pz, 5.) \
						+ 4. * a112 * pow(pz, 3.) * (pow(py, 2.) + pow(px, 2.)) \
						+ 2. * a112 * (pow(py, 4.) + pow(px, 4.)) * pz \
						+ 2. * a123 * pz * pow(py, 2.) * pow(px, 2.) \
						+ 8. * a1111 * pow(pz, 7.) \
						+ 2. * a1112 * pz * (pow(py, 6.) + pow(px, 6.) \
							+ 3. * pow(pz, 4.) * (pow(py, 2.) + pow(px, 2.))) \
						+ 4. * a1122 * pow(pz, 3.) * (pow(py, 4.) + pow(px, 4.)) \
						+ 2. * a1123 * pz * pow(px, 2.) * pow(py, 2.) * \
						(pow(px, 2.) + pow(py, 2.) + 2. * pow(pz, 2.)) \
						//--------Elastic energy effective field----------//
						+ (tv1 * (ezz0 - ezz) + tv2 * (exx0 - exx + eyy0 - eyy)) * pz \
						+ tv3 * ((exz0 - exz) * px + (eyz0 - eyz) * py));

					pt_math->transform_vector_crt2glb(dqx_crt_local, dqy_crt_local, dqz_crt_local, \
						dqx_glb_local, dqy_glb_local, dqz_glb_local);

					dqx_glb_local = dqx_glb_local + FE_mass * (\
						- 1. * G11 * (d2pxdx2_glb + d2pxdy2_glb + d2pxdz2_glb) \
						- Ex_stat(i, j, k) - pt_glb->Eext[0]) \
						+ FE_damping * qx_glb_store(i, j, k);

					dqy_glb_local = dqy_glb_local + FE_mass * (\
						- 1. * G11 * (d2pydx2_glb + d2pydy2_glb + d2pydz2_glb) \
						- Ey_stat(i, j, k) - pt_glb->Eext[1]) \
						+ FE_damping * qy_glb_store(i, j, k);

					dqz_glb_local = dqz_glb_local + FE_mass * (\
						- 1. * G11 * (d2pzdx2_glb + d2pzdy2_glb + d2pzdz2_glb) \
						- Ez_stat(i, j, k) - pt_glb->Eext[2]) \
						+ FE_damping * qz_glb_store(i, j, k);

					if (pt_glb->if_EM_backaction_on_P == true) {
						dqx_glb_local = dqx_glb_local - FE_mass * pt_EM->DEx_em_cell(i, j, k);
						dqy_glb_local = dqy_glb_local - FE_mass * pt_EM->DEy_em_cell(i, j, k);
						dqz_glb_local = dqz_glb_local - FE_mass * pt_EM->DEz_em_cell(i, j, k);
					}

					if (pt_glb->if_flexo == true) {
						//--------------------X compoennt-------------------//
						dqx_glb_local = dqx_glb_local + FE_mass * (\
							+ tv4 * px_glb_local * dpxdx(i, j, k) \
							+ tv5 * (py_glb_local * dpydx(i, j, k) + pz_glb_local * dpzdx(i, j, k)) \
							+ tv6 * (px_glb_local * dpydy(i, j, k) + px_glb_local * dpzdz(i, j, k) + py_glb_local * dpxdy(i, j, k) + pz_glb_local * dpxdz(i, j, k)) \
							- tv7 * d2pxdx2_glb \
							- tv8 * (d2pxdz2_glb + d2pxdy2_glb) \
							//-------------Flexoelectric energy effectiv field--------//
							- f11 * dexxdx - f12 * (deyydx + dezzdx) \
							- 2. * f44 * (dexydy + dexzdz));

						//--------------------Y compoennt-------------------//
						dqy_glb_local = dqy_glb_local + FE_mass * (\
							+ tv4 * py_glb_local * dpydy(i, j, k) \
							+ tv5 * (px_glb_local * dpxdy(i, j, k) + pz_glb_local * dpzdy(i, j, k)) \
							+ tv6 * (py_glb_local * dpxdx(i, j, k) + py_glb_local * dpzdz(i, j, k) + px_glb_local * dpydx(i, j, k) + pz_glb_local * dpydz(i, j, k)) \
							- tv7 * d2pydy2_glb \
							- tv8 * (d2pydx2_glb + d2pydz2_glb) \
							//-------------Flexoelectric energy effectiv field--------//
							- f11 * deyydy - f12 * (dexxdy + dezzdy) \
							- 2. * f44 * (dexydx + deyzdz));

						//--------------------Z compoennt-------------------//
						dqz_glb_local = dqz_glb_local + FE_mass * (\
							+ tv4 * pz_glb_local * dpzdz(i, j, k) \
							+ tv5 * (px_glb_local * dpxdz(i, j, k) + py_glb_local * dpydz(i, j, k)) \
							+ tv6 * (pz_glb_local * dpxdx(i, j, k) + pz_glb_local * dpydy(i, j, k) + px_glb_local * dpzdx(i, j, k) + py_glb_local * dpzdy(i, j, k)) \
							- tv7 * d2pzdz2_glb\
							- tv8 * (d2pzdx2_glb + d2pzdy2_glb) \
							//-------------Flexoelectric energy effectiv field--------//
							- f11 * dezzdz - f12 * (dexxdz + deyydz) \
							- 2. * f44 * (deyzdy + dexzdx));
					}

					dqx_glb_rk4(i, j, k) = dqx_glb_local; //RK
					dqy_glb_rk4(i, j, k) = dqy_glb_local; //RK
					dqz_glb_rk4(i, j, k) = dqz_glb_local; //RK
				}
			}
		}
	}
}

//void ferroelectric_system::get_dp() {
//#pragma acc parallel default(present)
//	{
//#pragma acc loop gang vector
//		for (long int id = 0; id < n; id++) {
//			dpx_glb(id) = qx_glb(id) * pt_glb->dt;
//			dpy_glb(id) = qy_glb(id) * pt_glb->dt;
//			dpz_glb(id) = qz_glb(id) * pt_glb->dt;
//		}
//	}
//}

void ferroelectric_system::get_dp_RK1() {
#pragma acc parallel default(present) async(6)
	{
#pragma acc loop gang vector
		for (long int id = 0; id < n; id++) {
			dpx_glb_rk1(id) = qx_glb_store(id) * pt_glb->dt;
			dpy_glb_rk1(id) = qy_glb_store(id) * pt_glb->dt;
			dpz_glb_rk1(id) = qz_glb_store(id) * pt_glb->dt;
		}
	}
}

void ferroelectric_system::get_dp_RK2() {
#pragma acc parallel default(present) async(6)
	{
#pragma acc loop gang vector
		for (long int id = 0; id < n; id++) {
			dpx_glb_rk2(id) = qx_glb_store(id) * pt_glb->dt;
			dpy_glb_rk2(id) = qy_glb_store(id) * pt_glb->dt;
			dpz_glb_rk2(id) = qz_glb_store(id) * pt_glb->dt;
		}
	}
}

void ferroelectric_system::get_dp_RK3() {
#pragma acc parallel default(present) async(6)
	{
#pragma acc loop gang vector
		for (long int id = 0; id < n; id++) {
			dpx_glb_rk3(id) = qx_glb_store(id) * pt_glb->dt;
			dpy_glb_rk3(id) = qy_glb_store(id) * pt_glb->dt;
			dpz_glb_rk3(id) = qz_glb_store(id) * pt_glb->dt;
		}
	}
}

void ferroelectric_system::get_dp_RK4() {
#pragma acc parallel default(present) async(6)
	{
#pragma acc loop gang vector
		for (long int id = 0; id < n; id++) {
			dpx_glb_rk4(id) = qx_glb_store(id) * pt_glb->dt;
			dpy_glb_rk4(id) = qy_glb_store(id) * pt_glb->dt;
			dpz_glb_rk4(id) = qz_glb_store(id) * pt_glb->dt;
		}
	}
}


void ferroelectric_system::update_q_RK1() {
#pragma acc parallel default(present) async(4)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			qx_glb_store(id) = qx_glb(id) + dqx_glb_rk1(id) * 0.5;
			qy_glb_store(id) = qy_glb(id) + dqy_glb_rk1(id) * 0.5;
			qz_glb_store(id) = qz_glb(id) + dqz_glb_rk1(id) * 0.5;
		}
	}
}

void ferroelectric_system::update_q_RK2() {
#pragma acc parallel default(present) async(4)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			qx_glb_store(id) = qx_glb(id) + dqx_glb_rk2(id) * 0.5;
			qy_glb_store(id) = qy_glb(id) + dqy_glb_rk2(id) * 0.5;
			qz_glb_store(id) = qz_glb(id) + dqz_glb_rk2(id) * 0.5;
		}
	}
}

void ferroelectric_system::update_q_RK3() {
#pragma acc parallel default(present) async(4)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			qx_glb_store(id) = qx_glb(id) + dqx_glb_rk3(id);
			qy_glb_store(id) = qy_glb(id) + dqy_glb_rk3(id);
			qz_glb_store(id) = qz_glb(id) + dqz_glb_rk3(id);
		}
	}
}

void ferroelectric_system::update_q() {
#pragma acc parallel default(present) async(4)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			qx_glb(id) = qx_glb(id) + dqx_glb_rk1(id) / 6. + dqx_glb_rk2(id) / 3. + dqx_glb_rk3(id) / 3. + dqx_glb_rk4(id) / 6.;
			qx_glb_store(id) = qx_glb(id);
			qy_glb(id) = qy_glb(id) + dqy_glb_rk1(id) / 6. + dqy_glb_rk2(id) / 3. + dqy_glb_rk3(id) / 3. + dqy_glb_rk4(id) / 6.;
			qy_glb_store(id) = qy_glb(id);
			qz_glb(id) = qz_glb(id) + dqz_glb_rk1(id) / 6. + dqz_glb_rk2(id) / 3. + dqz_glb_rk3(id) / 3. + dqz_glb_rk4(id) / 6.;
			qz_glb_store(id) = qz_glb(id);
		}
	}
}

void ferroelectric_system::update_p_RK1() {
#pragma acc parallel default(present) async(5)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			px_glb_store(id) = px_glb(id) + dpx_glb_rk1(id) * 0.5;
			py_glb_store(id) = py_glb(id) + dpy_glb_rk1(id) * 0.5;
			pz_glb_store(id) = pz_glb(id) + dpz_glb_rk1(id) * 0.5;
		}
	}

	if (pt_glb->if_EMdynamic == false && pt_glb->if_elec_1Dmodel == true) {
#pragma acc wait(5) async(6)
#pragma acc parallel default(present) async(6)
		{
			get_E_static_1D();
		}
	}

	if (pt_glb->if_flexo == true) {
#pragma acc wait(5) async(7)
#pragma acc parallel default(present) async(7)
		{
			get_p_gradient();
		}
	}
}

void ferroelectric_system::update_p_RK2() {
#pragma acc parallel default(present) async(5)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			px_glb_store(id) = px_glb(id) + dpx_glb_rk2(id) * 0.5;
			py_glb_store(id) = py_glb(id) + dpy_glb_rk2(id) * 0.5;
			pz_glb_store(id) = pz_glb(id) + dpz_glb_rk2(id) * 0.5;
		}
	}

	if (pt_glb->if_EMdynamic == false && pt_glb->if_elec_1Dmodel == true) {
#pragma acc wait(5) async(6)
#pragma acc parallel default(present) async(6)
		{
			get_E_static_1D();
		}
	}

	if (pt_glb->if_flexo == true) {
#pragma acc wait(5) async(7)
#pragma acc parallel default(present) async(7)
		{
			get_p_gradient();
		}
	}
}

void ferroelectric_system::update_p_RK3() {
#pragma acc parallel default(present) async(5)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			px_glb_store(id) = px_glb(id) + dpx_glb_rk3(id);
			py_glb_store(id) = py_glb(id) + dpy_glb_rk3(id);
			pz_glb_store(id) = pz_glb(id) + dpz_glb_rk3(id);
		}
	}

	if (pt_glb->if_EMdynamic == false && pt_glb->if_elec_1Dmodel == true) {
#pragma acc wait(5) async(6)
#pragma acc parallel default(present) async(6)
		{
			get_E_static_1D();
		}
	}

	if (pt_glb->if_flexo == true) {
#pragma acc wait(5) async(7)
#pragma acc parallel default(present) async(7)
		{
			get_p_gradient();
		}
	}
}

void ferroelectric_system::update_p() {
#pragma acc parallel default(present) async(5)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			px_glb(id) = px_glb(id) + dpx_glb_rk1(id) / 6. + dpx_glb_rk2(id) / 3. + dpx_glb_rk3(id) / 3. + dpx_glb_rk4(id) / 6.;
			px_glb_store(id) = px_glb(id);
			py_glb(id) = py_glb(id) + dpy_glb_rk1(id) / 6. + dpy_glb_rk2(id) / 3. + dpy_glb_rk3(id) / 3. + dpy_glb_rk4(id) / 6.;
			py_glb_store(id) = py_glb(id);
			pz_glb(id) = pz_glb(id) + dpz_glb_rk1(id) / 6. + dpz_glb_rk2(id) / 3. + dpz_glb_rk3(id) / 3. + dpz_glb_rk4(id) / 6.;
			pz_glb_store(id) = pz_glb(id);
		}
	}

	if (pt_glb->if_EMdynamic == false && pt_glb->if_elec_1Dmodel == true) {
#pragma acc wait(5) async(6)
#pragma acc parallel default(present) async(6)
		{
			get_E_static_1D();
		}
	}

	if (pt_glb->if_flexo == true) {
#pragma acc wait(5) async(7)
#pragma acc parallel default(present) async(7)
		{
			get_p_gradient();
		}
	}
}

#pragma acc routine seq
void ferroelectric_system::update_external_Efield() {
	double local_time;

	local_time = pt_glb->time_device - pt_glb->dt;
	pt_glb->Eext[0] = pt_glb->Eext_stat[0] + pt_glb->Eext_altr[0] * sin(2. * PI * pt_glb->E_altr_freq[0] * local_time);
	pt_glb->Eext[1] = pt_glb->Eext_stat[1] + pt_glb->Eext_altr[1] * sin(2. * PI * pt_glb->E_altr_freq[1] * local_time);
	pt_glb->Eext[2] = pt_glb->Eext_stat[2] + pt_glb->Eext_altr[2] * sin(2. * PI * pt_glb->E_altr_freq[2] * local_time);
}

#pragma acc routine seq
void ferroelectric_system::update_external_Efield_half() {
	double local_time;

	local_time = pt_glb->time_device - pt_glb->dt * 0.5;
	pt_glb->Eext[0] = pt_glb->Eext_stat[0] + pt_glb->Eext_altr[0] * sin(2. * PI * pt_glb->E_altr_freq[0] * local_time);
	pt_glb->Eext[1] = pt_glb->Eext_stat[1] + pt_glb->Eext_altr[1] * sin(2. * PI * pt_glb->E_altr_freq[1] * local_time);
	pt_glb->Eext[2] = pt_glb->Eext_stat[2] + pt_glb->Eext_altr[2] * sin(2. * PI * pt_glb->E_altr_freq[2] * local_time);
}

#pragma acc routine seq
void ferroelectric_system::update_external_Efield_full() {
	double local_time;

	local_time = pt_glb->time_device;
	pt_glb->Eext[0] = pt_glb->Eext_stat[0] + pt_glb->Eext_altr[0] * sin(2. * PI * pt_glb->E_altr_freq[0] * local_time);
	pt_glb->Eext[1] = pt_glb->Eext_stat[1] + pt_glb->Eext_altr[1] * sin(2. * PI * pt_glb->E_altr_freq[1] * local_time);
	pt_glb->Eext[2] = pt_glb->Eext_stat[2] + pt_glb->Eext_altr[2] * sin(2. * PI * pt_glb->E_altr_freq[2] * local_time);
}

void ferroelectric_system::get_averagep() {
	double ax, ay, az;
	ax = 0.; ay = 0.; az = 0.;
#pragma acc parallel default(present)
	{
#pragma acc loop gang vector reduction(+:ax)
		for (long int id = 0; id < n; id++) {
			ax = ax + px_glb(id);
		}

#pragma acc loop gang vector reduction(+:ay)
		for (long int id = 0; id < n; id++) {
			ay = ay + py_glb(id);
		}

#pragma acc loop gang vector reduction(+:az)
		for (long int id = 0; id < n; id++) {
			az = az + pz_glb(id);
		}
	}

	px_ave = ax / static_cast<double>(pt_glb->NFE);
	py_ave = ay / static_cast<double>(pt_glb->NFE);
	pz_ave = az / static_cast<double>(pt_glb->NFE);
}


