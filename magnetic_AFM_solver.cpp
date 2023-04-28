#include "magnetic_system.h"

void magnetic_system::get_dm_AFM_RK1() { //RK
	double Hx_anis_crt_AFM1, Hy_anis_crt_AFM1, Hz_anis_crt_AFM1;
	double Hx_anis_glb_AFM1, Hy_anis_glb_AFM1, Hz_anis_glb_AFM1;
	double Hx_elas_crt_AFM1, Hy_elas_crt_AFM1, Hz_elas_crt_AFM1;
	double Hx_elas_glb_AFM1, Hy_elas_glb_AFM1, Hz_elas_glb_AFM1;
	double Hx_exch_glb_AFM1, Hy_exch_glb_AFM1, Hz_exch_glb_AFM1;
	double Hx_iAFM_glb_AFM1, Hy_iAFM_glb_AFM1, Hz_iAFM_glb_AFM1;
	double Hx_eff_glb_AFM1, Hy_eff_glb_AFM1, Hz_eff_glb_AFM1;
	double precess_termx_AFM1, precess_termy_AFM1, precess_termz_AFM1;
	double damping_termx_AFM1, damping_termy_AFM1, damping_termz_AFM1;
	//double Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb;

	double Hx_anis_crt_AFM2, Hy_anis_crt_AFM2, Hz_anis_crt_AFM2;
	double Hx_anis_glb_AFM2, Hy_anis_glb_AFM2, Hz_anis_glb_AFM2;
	double Hx_elas_crt_AFM2, Hy_elas_crt_AFM2, Hz_elas_crt_AFM2;
	double Hx_elas_glb_AFM2, Hy_elas_glb_AFM2, Hz_elas_glb_AFM2;
	double Hx_exch_glb_AFM2, Hy_exch_glb_AFM2, Hz_exch_glb_AFM2;
	double Hx_iAFM_glb_AFM2, Hy_iAFM_glb_AFM2, Hz_iAFM_glb_AFM2;
	double Hx_eff_glb_AFM2, Hy_eff_glb_AFM2, Hz_eff_glb_AFM2;
	double precess_termx_AFM2, precess_termy_AFM2, precess_termz_AFM2;
	double damping_termx_AFM2, damping_termy_AFM2, damping_termz_AFM2;
	//double Hx_iDMI_glb_AFM2, Hy_iDMI_glb_AFM2, Hz_iDMI_glb_AFM2;

	long int x, y, z;
	unsigned int mat_type;
	material* mat;
	double Ms_AFM1, K1_AFM1, K2_AFM1, B1_AFM1, B2_AFM1;
	double Ms_AFM2, K1_AFM2, K2_AFM2, B1_AFM2, B2_AFM2;
	double ux, uy, uz;
	double udotm_AFM1, udotm_AFM2;
	//double Aexi;//, iDMIi;
	double JAFM_AFM1, JAFM_AFM2;
	double alpha_AFM1, gyro_AFM1, precess_AFM1;
	double alpha_AFM2, gyro_AFM2, precess_AFM2;
	double mx_glb_local_AFM1, my_glb_local_AFM1, mz_glb_local_AFM1;
	double mx_glb_local_AFM2, my_glb_local_AFM2, mz_glb_local_AFM2;
	double mx_AFM1, my_AFM1, mz_AFM1;
	double mx_AFM2, my_AFM2, mz_AFM2;
	double exx, eyy, ezz, eyz, exz, exy;
	//double d2mxdx2, d2mxdy2, d2mxdz2;
	//double d2mydx2, d2mydy2, d2mydz2;
	//double d2mzdx2, d2mzdy2, d2mzdz2;
	double dmxdx_AFM1, dmydy_AFM1, dmzdx_AFM1, dmzdy_AFM1;
	double dmxdx2_f_AFM1, dmydx2_f_AFM1, dmzdx2_f_AFM1;
	double dmxdx2_b_AFM1, dmydx2_b_AFM1, dmzdx2_b_AFM1;
	double dmxdy2_f_AFM1, dmydy2_f_AFM1, dmzdy2_f_AFM1;
	double dmxdy2_b_AFM1, dmydy2_b_AFM1, dmzdy2_b_AFM1;
	double dmxdz2_f_AFM1, dmydz2_f_AFM1, dmzdz2_f_AFM1;
	double dmxdz2_b_AFM1, dmydz2_b_AFM1, dmzdz2_b_AFM1;

	double dmxdx_AFM2, dmydy_AFM2, dmzdx_AFM2, dmzdy_AFM2;
	double dmxdx2_f_AFM2, dmydx2_f_AFM2, dmzdx2_f_AFM2;
	double dmxdx2_b_AFM2, dmydx2_b_AFM2, dmzdx2_b_AFM2;
	double dmxdy2_f_AFM2, dmydy2_f_AFM2, dmzdy2_f_AFM2;
	double dmxdy2_b_AFM2, dmydy2_b_AFM2, dmzdy2_b_AFM2;
	double dmxdz2_f_AFM2, dmydz2_f_AFM2, dmzdz2_f_AFM2;
	double dmxdz2_b_AFM2, dmydz2_b_AFM2, dmzdz2_b_AFM2;

#pragma acc serial default(present) async(4)
	{
		update_external_Hfield();  //RK
	}

#pragma acc parallel default(present) async(4)
	{
		//------------------Sublattice 1 ---------------------//
#pragma acc loop gang vector private(\
	Hx_anis_crt_AFM1, Hy_anis_crt_AFM1, Hz_anis_crt_AFM1,\
	Hx_anis_glb_AFM1, Hy_anis_glb_AFM1, Hz_anis_glb_AFM1,\
	Hx_elas_crt_AFM1, Hy_elas_crt_AFM1, Hz_elas_crt_AFM1,\
	Hx_elas_glb_AFM1, Hy_elas_glb_AFM1, Hz_elas_glb_AFM1,\
	Hx_exch_glb_AFM1, Hy_exch_glb_AFM1, Hz_exch_glb_AFM1,\
	Hx_iAFM_glb_AFM1, Hy_iAFM_glb_AFM1, Hz_iAFM_glb_AFM1,\
	Hx_eff_glb_AFM1, Hy_eff_glb_AFM1, Hz_eff_glb_AFM1,\
	precess_termx_AFM1, precess_termy_AFM1, precess_termz_AFM1,\
	damping_termx_AFM1, damping_termy_AFM1, damping_termz_AFM1,\
	x,y,z,mat_type,mat,\
	Ms_AFM1, K1_AFM1, K2_AFM1, B1_AFM1, B2_AFM1, \
	ux, uy, uz, udotm_AFM1, \
	JAFM_AFM1,\
	alpha_AFM1, gyro_AFM1, precess_AFM1,\
	mx_glb_local_AFM1, my_glb_local_AFM1, mz_glb_local_AFM1,\
	mx_AFM1, my_AFM1, mz_AFM1,\
	exx, eyy, ezz, eyz, exz, exy,\
	dmxdx_AFM1, dmydy_AFM1, dmzdx_AFM1, dmzdy_AFM1,\
	dmxdx2_f_AFM1, dmydx2_f_AFM1, dmzdx2_f_AFM1,\
	dmxdx2_b_AFM1, dmydx2_b_AFM1, dmzdx2_b_AFM1,\
	dmxdy2_f_AFM1, dmydy2_f_AFM1, dmzdy2_f_AFM1,\
	dmxdy2_b_AFM1, dmydy2_b_AFM1, dmzdy2_b_AFM1,\
	dmxdz2_f_AFM1, dmydz2_f_AFM1, dmzdz2_f_AFM1,\
	dmxdz2_b_AFM1, dmydz2_b_AFM1, dmzdz2_b_AFM1)
//Aexi, Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb, iDMIi)
		for (long int i = 0; i < n; i++) {
			mat_type = pt_glb->material_cell(i);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_AFM == true) {
					x = i / (ny * nz);
					y = (i - x * (ny * nz)) / nz;
					z = i - x * (ny * nz) - y * nz;

					ux = mat->uniaxis_AFM[0]; uy = mat->uniaxis_AFM[1]; uz = mat->uniaxis_AFM[2];

					if (pt_glb->if_elasto_on_pORm == true) {
						exx = pt_elas->exxt0_crt(i) + pt_elas->Dexx_crt(i);
						eyy = pt_elas->eyyt0_crt(i) + pt_elas->Deyy_crt(i);
						ezz = pt_elas->ezzt0_crt(i) + pt_elas->Dezz_crt(i);
						eyz = pt_elas->eyzt0_crt(i) + pt_elas->Deyz_crt(i);
						exz = pt_elas->exzt0_crt(i) + pt_elas->Dexz_crt(i);
						exy = pt_elas->exyt0_crt(i) + pt_elas->Dexy_crt(i);
					}
					else {
						exx = 0.;
						eyy = 0.;
						ezz = 0.;
						eyz = 0.;
						exz = 0.;
						exy = 0.;
					}

					mx_glb_local_AFM1 = mx_AFM1_glb_store(i);
					my_glb_local_AFM1 = my_AFM1_glb_store(i);
					mz_glb_local_AFM1 = mz_AFM1_glb_store(i);

					pt_math->transform_vector_glb2crt(mx_glb_local_AFM1, my_glb_local_AFM1, mz_glb_local_AFM1, \
						mx_AFM1, my_AFM1, mz_AFM1);

					Ms_AFM1 = mat->Ms_AFM1;
					K1_AFM1 = mat->K1_AFM1 / mu0 / Ms_AFM1; //K2_AFM1 = mat->K2 / mu0 / Ms;					
					udotm_AFM1 = mx_AFM1 * ux + my_AFM1 * uy + mz_AFM1 * uz;
					B1_AFM1 = mat->B1_AFM1 / mu0 / Ms_AFM1; B2_AFM1 = mat->B2_AFM1 / mu0 / Ms_AFM1;
					//Aexi_AFM1 = mat->Aex_AFM1; //iDMIi = mat->iDMI;
					JAFM_AFM1 = -1. * mat->J_AFM / mu0 / Ms_AFM1;

					get_m_AFM1_spatial_derivative_local \
						(x, y, z, \
							dmxdx_AFM1, dmydy_AFM1, dmzdx_AFM1, dmzdy_AFM1, \
							dmxdx2_f_AFM1, dmydx2_f_AFM1, dmzdx2_f_AFM1, \
							dmxdx2_b_AFM1, dmydx2_b_AFM1, dmzdx2_b_AFM1, \
							dmxdy2_f_AFM1, dmydy2_f_AFM1, dmzdy2_f_AFM1, \
							dmxdy2_b_AFM1, dmydy2_b_AFM1, dmzdy2_b_AFM1, \
							dmxdz2_f_AFM1, dmydz2_f_AFM1, dmzdz2_f_AFM1, \
							dmxdz2_b_AFM1, dmydz2_b_AFM1, dmzdz2_b_AFM1);

					//get_laplacian_m_local \
						//	(x, y, z, \
						//		d2mxdx2, d2mxdy2, d2mxdz2, \
						//		d2mydx2, d2mydy2, d2mydz2, \
						//		d2mzdx2, d2mzdy2, d2mzdz2);

						//------------Magnetocrystalline anisotropy field--------------//
					Hx_anis_crt_AFM1 = 2. * K1_AFM1 * udotm_AFM1 * ux;
					Hy_anis_crt_AFM1 = 2. * K1_AFM1 * udotm_AFM1 * uy;
					Hz_anis_crt_AFM1 = 2. * K1_AFM1 * udotm_AFM1 * uz;

					pt_math->transform_vector_crt2glb(Hx_anis_crt_AFM1, Hy_anis_crt_AFM1, Hz_anis_crt_AFM1, \
						Hx_anis_glb_AFM1, Hy_anis_glb_AFM1, Hz_anis_glb_AFM1);

					//Hx_anis(i) = Hx_anis_glb;
					//Hy_anis(i) = Hy_anis_glb;
					//Hz_anis(i) = Hz_anis_glb;

					//------------------Magnetoelastic effective field---------------//
					Hx_elas_crt_AFM1 = -1. * (B1_AFM1 * (exx)*mx_AFM1 + \
						B2_AFM1 * ((exy)*my_AFM1 + (exz)*mz_AFM1));
					Hy_elas_crt_AFM1 = -1. * (B1_AFM1 * (eyy)*my_AFM1 + \
						B2_AFM1 * ((exy)*mx_AFM1 + (eyz)*mz_AFM1));
					Hz_elas_crt_AFM1 = -1. * (B1_AFM1 * (ezz)*mz_AFM1 + \
						B2_AFM1 * ((exz)*mx_AFM1 + (eyz)*my_AFM1));

					pt_math->transform_vector_crt2glb(Hx_elas_crt_AFM1, Hy_elas_crt_AFM1, Hz_elas_crt_AFM1, \
						Hx_elas_glb_AFM1, Hy_elas_glb_AFM1, Hz_elas_glb_AFM1);

					//Hx_elas(i) = Hx_elas_glb;
					//Hy_elas(i) = Hy_elas_glb;
					//Hz_elas(i) = Hz_elas_glb;

					//-----------------Magnetic exchange field----------------------//
					//----------Averaged parameter for inter-region exchange - Approach 1---------//
					Hx_exch_glb_AFM1 = 2. / mu0 / Ms_AFM1 * \
							 (Aij_AFM1_xf(i) * dmxdx2_f_AFM1 + Aij_AFM1_xb(i) * dmxdx2_b_AFM1 \
							+ Aij_AFM1_yf(i) * dmxdy2_f_AFM1 + Aij_AFM1_yb(i) * dmxdy2_b_AFM1 \
							+ Aij_AFM1_zf(i) * dmxdz2_f_AFM1 + Aij_AFM1_zb(i) * dmxdz2_b_AFM1);
					Hy_exch_glb_AFM1 = 2. / mu0 / Ms_AFM1 * \
							 (Aij_AFM1_xf(i) * dmydx2_f_AFM1 + Aij_AFM1_xb(i) * dmydx2_b_AFM1 \
							+ Aij_AFM1_yf(i) * dmydy2_f_AFM1 + Aij_AFM1_yb(i) * dmydy2_b_AFM1 \
							+ Aij_AFM1_zf(i) * dmydz2_f_AFM1 + Aij_AFM1_zb(i) * dmydz2_b_AFM1);
					Hz_exch_glb_AFM1 = 2. / mu0 / Ms_AFM1 * \
							 (Aij_AFM1_xf(i) * dmzdx2_f_AFM1 + Aij_AFM1_xb(i) * dmzdx2_b_AFM1 \
							+ Aij_AFM1_yf(i) * dmzdy2_f_AFM1 + Aij_AFM1_yb(i) * dmzdy2_b_AFM1 \
							+ Aij_AFM1_zf(i) * dmzdz2_f_AFM1 + Aij_AFM1_zb(i) * dmzdz2_b_AFM1);

					//Hx_exch(i) = 2. * Aex * (d2mxdx2 + d2mxdy2 + d2mxdz2);
					//Hy_exch(i) = 2. * Aex * (d2mydx2 + d2mydy2 + d2mydz2);
					//Hz_exch(i) = 2. * Aex * (d2mzdx2 + d2mzdy2 + d2mzdz2);
					// 
					// -----------------Inter-sublattice AFM exchange field---------------//
					Hx_iAFM_glb_AFM1 = JAFM_AFM1 * mx_AFM2_glb_store(i);
					Hy_iAFM_glb_AFM1 = JAFM_AFM1 * my_AFM2_glb_store(i);
					Hz_iAFM_glb_AFM1 = JAFM_AFM1 * mz_AFM2_glb_store(i);

					//------------------Interfacial DMI effective field-------------------//
					//Hx_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdx;
					//Hy_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdy;
					//Hz_iDMI_glb = -2. * iDMIi / mu0 / Ms * (dmxdx + dmydy);

					//-------------------LLG Equation------------------------------//
					Hx_eff_glb_AFM1 = Hx_stat(i) + \
						Hx_anis_glb_AFM1 + Hx_elas_glb_AFM1 + \
						Hx_exch_glb_AFM1 + Hx_iAFM_glb_AFM1 /*+ Hx_iDMI_glb*/ + pt_glb->Hext[0];

					Hy_eff_glb_AFM1 = Hy_stat(i) + \
						Hy_anis_glb_AFM1 + Hy_elas_glb_AFM1 + \
						Hy_exch_glb_AFM1 + Hy_iAFM_glb_AFM1/*+ Hy_iDMI_glb*/ + pt_glb->Hext[1];

					Hz_eff_glb_AFM1 = Hz_stat(i) + \
						Hz_anis_glb_AFM1 + Hz_elas_glb_AFM1 + \
						Hz_exch_glb_AFM1 + Hz_iAFM_glb_AFM1/*+ Hz_iDMI_glb*/ + pt_glb->Hext[2];

					if (pt_glb->if_EM_backaction_on_M == true) {
						Hx_eff_glb_AFM1 = Hx_eff_glb_AFM1 + pt_EM->DHx_em_cell(i);
						Hy_eff_glb_AFM1 = Hy_eff_glb_AFM1 + pt_EM->DHy_em_cell(i);
						Hz_eff_glb_AFM1 = Hz_eff_glb_AFM1 + pt_EM->DHz_em_cell(i);
					}

					alpha_AFM1 = mat->damping_AFM1;
					gyro_AFM1 = mat->gyro_AFM1;
					precess_AFM1 = -gyro_AFM1 / (1. + alpha_AFM1 * alpha_AFM1) * pt_glb->dt;

					precess_termx_AFM1 = my_glb_local_AFM1 * Hz_eff_glb_AFM1 - mz_glb_local_AFM1 * Hy_eff_glb_AFM1;
					precess_termy_AFM1 = mz_glb_local_AFM1 * Hx_eff_glb_AFM1 - mx_glb_local_AFM1 * Hz_eff_glb_AFM1;
					precess_termz_AFM1 = mx_glb_local_AFM1 * Hy_eff_glb_AFM1 - my_glb_local_AFM1 * Hx_eff_glb_AFM1;

					damping_termx_AFM1 = my_glb_local_AFM1 * precess_termz_AFM1 - mz_glb_local_AFM1 * precess_termy_AFM1;
					damping_termy_AFM1 = mz_glb_local_AFM1 * precess_termx_AFM1 - mx_glb_local_AFM1 * precess_termz_AFM1;
					damping_termz_AFM1 = mx_glb_local_AFM1 * precess_termy_AFM1 - my_glb_local_AFM1 * precess_termx_AFM1;

					dmx_AFM1_glb_rk1(i) = precess_AFM1 * precess_termx_AFM1 + alpha_AFM1 * precess_AFM1 * damping_termx_AFM1; //RK
					dmy_AFM1_glb_rk1(i) = precess_AFM1 * precess_termy_AFM1 + alpha_AFM1 * precess_AFM1 * damping_termy_AFM1; //RK
					dmz_AFM1_glb_rk1(i) = precess_AFM1 * precess_termz_AFM1 + alpha_AFM1 * precess_AFM1 * damping_termz_AFM1; //RK
				}
			}
		}

		//------------------Sublattice 2 ---------------------//
#pragma acc loop gang vector private(\
	Hx_anis_crt_AFM2, Hy_anis_crt_AFM2, Hz_anis_crt_AFM2,\
	Hx_anis_glb_AFM2, Hy_anis_glb_AFM2, Hz_anis_glb_AFM2,\
	Hx_elas_crt_AFM2, Hy_elas_crt_AFM2, Hz_elas_crt_AFM2,\
	Hx_elas_glb_AFM2, Hy_elas_glb_AFM2, Hz_elas_glb_AFM2,\
	Hx_exch_glb_AFM2, Hy_exch_glb_AFM2, Hz_exch_glb_AFM2,\
	Hx_iAFM_glb_AFM2, Hy_iAFM_glb_AFM2, Hz_iAFM_glb_AFM2,\
	Hx_eff_glb_AFM2, Hy_eff_glb_AFM2, Hz_eff_glb_AFM2,\
	precess_termx_AFM2, precess_termy_AFM2, precess_termz_AFM2,\
	damping_termx_AFM2, damping_termy_AFM2, damping_termz_AFM2,\
	x,y,z,mat_type,mat,\
	Ms_AFM2, K1_AFM2, K2_AFM2, B1_AFM2, B2_AFM2, \
	ux, uy, uz, udotm_AFM2, \
	JAFM_AFM2,\
	alpha_AFM2, gyro_AFM2, precess_AFM2,\
	mx_glb_local_AFM2, my_glb_local_AFM2, mz_glb_local_AFM2,\
	mx_AFM2, my_AFM2, mz_AFM2,\
	exx, eyy, ezz, eyz, exz, exy,\
	dmxdx_AFM2, dmydy_AFM2, dmzdx_AFM2, dmzdy_AFM2,\
	dmxdx2_f_AFM2, dmydx2_f_AFM2, dmzdx2_f_AFM2,\
	dmxdx2_b_AFM2, dmydx2_b_AFM2, dmzdx2_b_AFM2,\
	dmxdy2_f_AFM2, dmydy2_f_AFM2, dmzdy2_f_AFM2,\
	dmxdy2_b_AFM2, dmydy2_b_AFM2, dmzdy2_b_AFM2,\
	dmxdz2_f_AFM2, dmydz2_f_AFM2, dmzdz2_f_AFM2,\
	dmxdz2_b_AFM2, dmydz2_b_AFM2, dmzdz2_b_AFM2)
	//Aexi, Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb, iDMIi)
		for (long int i = 0; i < n; i++) {
			mat_type = pt_glb->material_cell(i);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_AFM == true) {
					x = i / (ny * nz);
					y = (i - x * (ny * nz)) / nz;
					z = i - x * (ny * nz) - y * nz;

					ux = mat->uniaxis_AFM[0]; uy = mat->uniaxis_AFM[1]; uz = mat->uniaxis_AFM[2];

					if (pt_glb->if_elasto_on_pORm == true) {
						exx = pt_elas->exxt0_crt(i) + pt_elas->Dexx_crt(i);
						eyy = pt_elas->eyyt0_crt(i) + pt_elas->Deyy_crt(i);
						ezz = pt_elas->ezzt0_crt(i) + pt_elas->Dezz_crt(i);
						eyz = pt_elas->eyzt0_crt(i) + pt_elas->Deyz_crt(i);
						exz = pt_elas->exzt0_crt(i) + pt_elas->Dexz_crt(i);
						exy = pt_elas->exyt0_crt(i) + pt_elas->Dexy_crt(i);
					}
					else {
						exx = 0.;
						eyy = 0.;
						ezz = 0.;
						eyz = 0.;
						exz = 0.;
						exy = 0.;
					}
				
					mx_glb_local_AFM2 = mx_AFM2_glb_store(i);
					my_glb_local_AFM2 = my_AFM2_glb_store(i);
					mz_glb_local_AFM2 = mz_AFM2_glb_store(i);

					pt_math->transform_vector_glb2crt(mx_glb_local_AFM2, my_glb_local_AFM2, mz_glb_local_AFM2, \
						mx_AFM2, my_AFM2, mz_AFM2);

					Ms_AFM2 = mat->Ms_AFM2;
					K1_AFM2 = mat->K1_AFM2 / mu0 / Ms_AFM2; //K2 = mat->K2_AFM2 / mu0 / Ms;					
					udotm_AFM2 = mx_AFM2 * ux + my_AFM2 * uy + mz_AFM2 * uz;
					B1_AFM2 = mat->B1_AFM2 / mu0 / Ms_AFM2; B2_AFM2 = mat->B2_AFM2 / mu0 / Ms_AFM2;
					//Aexi_AFM2 = mat->Aex_AFM2; //iDMIi = mat->iDMI;
					JAFM_AFM2 = -1. * mat->J_AFM / mu0 / Ms_AFM2;

					get_m_AFM2_spatial_derivative_local \
						(x, y, z, \
							dmxdx_AFM2, dmydy_AFM2, dmzdx_AFM2, dmzdy_AFM2, \
							dmxdx2_f_AFM2, dmydx2_f_AFM2, dmzdx2_f_AFM2, \
							dmxdx2_b_AFM2, dmydx2_b_AFM2, dmzdx2_b_AFM2, \
							dmxdy2_f_AFM2, dmydy2_f_AFM2, dmzdy2_f_AFM2, \
							dmxdy2_b_AFM2, dmydy2_b_AFM2, dmzdy2_b_AFM2, \
							dmxdz2_f_AFM2, dmydz2_f_AFM2, dmzdz2_f_AFM2, \
							dmxdz2_b_AFM2, dmydz2_b_AFM2, dmzdz2_b_AFM2);

					//get_laplacian_m_local \
						//	(x, y, z, \
						//		d2mxdx2, d2mxdy2, d2mxdz2, \
						//		d2mydx2, d2mydy2, d2mydz2, \
						//		d2mzdx2, d2mzdy2, d2mzdz2);

						//------------Magnetocrystalline anisotropy field--------------//
					Hx_anis_crt_AFM2 = 2. * K1_AFM2 * udotm_AFM2 * ux;
					Hy_anis_crt_AFM2 = 2. * K1_AFM2 * udotm_AFM2 * uy;
					Hz_anis_crt_AFM2 = 2. * K1_AFM2 * udotm_AFM2 * uz;

					pt_math->transform_vector_crt2glb(Hx_anis_crt_AFM2, Hy_anis_crt_AFM2, Hz_anis_crt_AFM2, \
						Hx_anis_glb_AFM2, Hy_anis_glb_AFM2, Hz_anis_glb_AFM2);

					//Hx_anis(i) = Hx_anis_glb;
					//Hy_anis(i) = Hy_anis_glb;
					//Hz_anis(i) = Hz_anis_glb;

					//------------------Magnetoelastic effective field---------------//
					Hx_elas_crt_AFM2 = -1. * (B1_AFM2 * (exx)*mx_AFM2 + \
						B2_AFM2 * ((exy)*my_AFM2 + (exz)*mz_AFM2));
					Hy_elas_crt_AFM2 = -1. * (B1_AFM2 * (eyy)*my_AFM2 + \
						B2_AFM2 * ((exy)*mx_AFM2 + (eyz)*mz_AFM2));
					Hz_elas_crt_AFM2 = -1. * (B1_AFM2 * (ezz)*mz_AFM2 + \
						B2_AFM2 * ((exz)*mx_AFM2 + (eyz)*my_AFM2));

					pt_math->transform_vector_crt2glb(Hx_elas_crt_AFM2, Hy_elas_crt_AFM2, Hz_elas_crt_AFM2, \
						Hx_elas_glb_AFM2, Hy_elas_glb_AFM2, Hz_elas_glb_AFM2);

					//Hx_elas(i) = Hx_elas_glb;
					//Hy_elas(i) = Hy_elas_glb;
					//Hz_elas(i) = Hz_elas_glb;

					//-----------------Magnetic exchange field----------------------//
					//----------Averaged parameter for inter-region exchange - Approach 1---------//
					Hx_exch_glb_AFM2 = 2. / mu0 / Ms_AFM2 * \
							 (Aij_AFM2_xf(i) * dmxdx2_f_AFM2 + Aij_AFM2_xb(i) * dmxdx2_b_AFM2 \
							+ Aij_AFM2_yf(i) * dmxdy2_f_AFM2 + Aij_AFM2_yb(i) * dmxdy2_b_AFM2 \
							+ Aij_AFM2_zf(i) * dmxdz2_f_AFM2 + Aij_AFM2_zb(i) * dmxdz2_b_AFM2);
					Hy_exch_glb_AFM2 = 2. / mu0 / Ms_AFM2 * \
							 (Aij_AFM2_xf(i) * dmydx2_f_AFM2 + Aij_AFM2_xb(i) * dmydx2_b_AFM2 \
							+ Aij_AFM2_yf(i) * dmydy2_f_AFM2 + Aij_AFM2_yb(i) * dmydy2_b_AFM2 \
							+ Aij_AFM2_zf(i) * dmydz2_f_AFM2 + Aij_AFM2_zb(i) * dmydz2_b_AFM2);
					Hz_exch_glb_AFM2 = 2. / mu0 / Ms_AFM2 * \
							 (Aij_AFM2_xf(i) * dmzdx2_f_AFM2 + Aij_AFM2_xb(i) * dmzdx2_b_AFM2 \
							+ Aij_AFM2_yf(i) * dmzdy2_f_AFM2 + Aij_AFM2_yb(i) * dmzdy2_b_AFM2 \
							+ Aij_AFM2_zf(i) * dmzdz2_f_AFM2 + Aij_AFM2_zb(i) * dmzdz2_b_AFM2);

					//Hx_exch(i) = 2. * Aex * (d2mxdx2 + d2mxdy2 + d2mxdz2);
					//Hy_exch(i) = 2. * Aex * (d2mydx2 + d2mydy2 + d2mydz2);
					//Hz_exch(i) = 2. * Aex * (d2mzdx2 + d2mzdy2 + d2mzdz2);
					// 
					// -----------------Inter-sublattice AFM exchange field---------------//
					Hx_iAFM_glb_AFM2 = JAFM_AFM2 * mx_AFM1_glb_store(i);
					Hy_iAFM_glb_AFM2 = JAFM_AFM2 * my_AFM1_glb_store(i);
					Hz_iAFM_glb_AFM2 = JAFM_AFM2 * mz_AFM1_glb_store(i);
					//------------------Interfacial DMI effective field-------------------//
					//Hx_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdx;
					//Hy_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdy;
					//Hz_iDMI_glb = -2. * iDMIi / mu0 / Ms * (dmxdx + dmydy);

					//-------------------LLG Equation------------------------------//
					Hx_eff_glb_AFM2 = Hx_stat(i) + \
						Hx_anis_glb_AFM2 + Hx_elas_glb_AFM2 + \
						Hx_exch_glb_AFM2 + Hx_iAFM_glb_AFM2 /*+ Hx_iDMI_glb*/ + pt_glb->Hext[0];

					Hy_eff_glb_AFM2 = Hy_stat(i) + \
						Hy_anis_glb_AFM2 + Hy_elas_glb_AFM2 + \
						Hy_exch_glb_AFM2 + Hy_iAFM_glb_AFM2 /*+ Hy_iDMI_glb*/ + pt_glb->Hext[1];

					Hz_eff_glb_AFM2 = Hz_stat(i) + \
						Hz_anis_glb_AFM2 + Hz_elas_glb_AFM2 + \
						Hz_exch_glb_AFM2 + Hz_iAFM_glb_AFM2 /*+ Hz_iDMI_glb*/ + pt_glb->Hext[2];

					if (pt_glb->if_EM_backaction_on_M == true) {
						Hx_eff_glb_AFM2 = Hx_eff_glb_AFM2 + pt_EM->DHx_em_cell(i);
						Hy_eff_glb_AFM2 = Hy_eff_glb_AFM2 + pt_EM->DHy_em_cell(i);
						Hz_eff_glb_AFM2 = Hz_eff_glb_AFM2 + pt_EM->DHz_em_cell(i);
					}

					alpha_AFM2 = mat->damping_AFM2;
					gyro_AFM2 = mat->gyro_AFM2;
					precess_AFM2 = -gyro_AFM2 / (1 + alpha_AFM2 * alpha_AFM2) * pt_glb->dt;

					precess_termx_AFM2 = my_glb_local_AFM2 * Hz_eff_glb_AFM2 - mz_glb_local_AFM2 * Hy_eff_glb_AFM2;
					precess_termy_AFM2 = mz_glb_local_AFM2 * Hx_eff_glb_AFM2 - mx_glb_local_AFM2 * Hz_eff_glb_AFM2;
					precess_termz_AFM2 = mx_glb_local_AFM2 * Hy_eff_glb_AFM2 - my_glb_local_AFM2 * Hx_eff_glb_AFM2;

					damping_termx_AFM2 = my_glb_local_AFM2 * precess_termz_AFM2 - mz_glb_local_AFM2 * precess_termy_AFM2;
					damping_termy_AFM2 = mz_glb_local_AFM2 * precess_termx_AFM2 - mx_glb_local_AFM2 * precess_termz_AFM2;
					damping_termz_AFM2 = mx_glb_local_AFM2 * precess_termy_AFM2 - my_glb_local_AFM2 * precess_termx_AFM2;

					dmx_AFM2_glb_rk1(i) = precess_AFM2 * precess_termx_AFM2 + alpha_AFM2 * precess_AFM2 * damping_termx_AFM2; //RK
					dmy_AFM2_glb_rk1(i) = precess_AFM2 * precess_termy_AFM2 + alpha_AFM2 * precess_AFM2 * damping_termy_AFM2; //RK
					dmz_AFM2_glb_rk1(i) = precess_AFM2 * precess_termz_AFM2 + alpha_AFM2 * precess_AFM2 * damping_termz_AFM2; //RK
				}
			}
		}
	}

	if (pt_glb->if_spin_pumping == true) {
		get_J_ISHE_RK1(); //RK
	}
}

void magnetic_system::get_dm_AFM_RK2() { //RK
	double Hx_anis_crt_AFM1, Hy_anis_crt_AFM1, Hz_anis_crt_AFM1;
	double Hx_anis_glb_AFM1, Hy_anis_glb_AFM1, Hz_anis_glb_AFM1;
	double Hx_elas_crt_AFM1, Hy_elas_crt_AFM1, Hz_elas_crt_AFM1;
	double Hx_elas_glb_AFM1, Hy_elas_glb_AFM1, Hz_elas_glb_AFM1;
	double Hx_exch_glb_AFM1, Hy_exch_glb_AFM1, Hz_exch_glb_AFM1;
	double Hx_iAFM_glb_AFM1, Hy_iAFM_glb_AFM1, Hz_iAFM_glb_AFM1;
	double Hx_eff_glb_AFM1, Hy_eff_glb_AFM1, Hz_eff_glb_AFM1;
	double precess_termx_AFM1, precess_termy_AFM1, precess_termz_AFM1;
	double damping_termx_AFM1, damping_termy_AFM1, damping_termz_AFM1;
	//double Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb;

	double Hx_anis_crt_AFM2, Hy_anis_crt_AFM2, Hz_anis_crt_AFM2;
	double Hx_anis_glb_AFM2, Hy_anis_glb_AFM2, Hz_anis_glb_AFM2;
	double Hx_elas_crt_AFM2, Hy_elas_crt_AFM2, Hz_elas_crt_AFM2;
	double Hx_elas_glb_AFM2, Hy_elas_glb_AFM2, Hz_elas_glb_AFM2;
	double Hx_exch_glb_AFM2, Hy_exch_glb_AFM2, Hz_exch_glb_AFM2;
	double Hx_iAFM_glb_AFM2, Hy_iAFM_glb_AFM2, Hz_iAFM_glb_AFM2;
	double Hx_eff_glb_AFM2, Hy_eff_glb_AFM2, Hz_eff_glb_AFM2;
	double precess_termx_AFM2, precess_termy_AFM2, precess_termz_AFM2;
	double damping_termx_AFM2, damping_termy_AFM2, damping_termz_AFM2;
	//double Hx_iDMI_glb_AFM2, Hy_iDMI_glb_AFM2, Hz_iDMI_glb_AFM2;

	long int x, y, z;
	unsigned int mat_type;
	material* mat;
	double Ms_AFM1, K1_AFM1, K2_AFM1, B1_AFM1, B2_AFM1;
	double Ms_AFM2, K1_AFM2, K2_AFM2, B1_AFM2, B2_AFM2;
	double ux, uy, uz;
	double udotm_AFM1, udotm_AFM2;
	//double Aexi;//, iDMIi;
	double JAFM_AFM1, JAFM_AFM2;
	double alpha_AFM1, gyro_AFM1, precess_AFM1;
	double alpha_AFM2, gyro_AFM2, precess_AFM2;
	double mx_glb_local_AFM1, my_glb_local_AFM1, mz_glb_local_AFM1;
	double mx_glb_local_AFM2, my_glb_local_AFM2, mz_glb_local_AFM2;
	double mx_AFM1, my_AFM1, mz_AFM1;
	double mx_AFM2, my_AFM2, mz_AFM2;
	double exx, eyy, ezz, eyz, exz, exy;
	//double d2mxdx2, d2mxdy2, d2mxdz2;
	//double d2mydx2, d2mydy2, d2mydz2;
	//double d2mzdx2, d2mzdy2, d2mzdz2;
	double dmxdx_AFM1, dmydy_AFM1, dmzdx_AFM1, dmzdy_AFM1;
	double dmxdx2_f_AFM1, dmydx2_f_AFM1, dmzdx2_f_AFM1;
	double dmxdx2_b_AFM1, dmydx2_b_AFM1, dmzdx2_b_AFM1;
	double dmxdy2_f_AFM1, dmydy2_f_AFM1, dmzdy2_f_AFM1;
	double dmxdy2_b_AFM1, dmydy2_b_AFM1, dmzdy2_b_AFM1;
	double dmxdz2_f_AFM1, dmydz2_f_AFM1, dmzdz2_f_AFM1;
	double dmxdz2_b_AFM1, dmydz2_b_AFM1, dmzdz2_b_AFM1;

	double dmxdx_AFM2, dmydy_AFM2, dmzdx_AFM2, dmzdy_AFM2;
	double dmxdx2_f_AFM2, dmydx2_f_AFM2, dmzdx2_f_AFM2;
	double dmxdx2_b_AFM2, dmydx2_b_AFM2, dmzdx2_b_AFM2;
	double dmxdy2_f_AFM2, dmydy2_f_AFM2, dmzdy2_f_AFM2;
	double dmxdy2_b_AFM2, dmydy2_b_AFM2, dmzdy2_b_AFM2;
	double dmxdz2_f_AFM2, dmydz2_f_AFM2, dmzdz2_f_AFM2;
	double dmxdz2_b_AFM2, dmydz2_b_AFM2, dmzdz2_b_AFM2;

#pragma acc serial default(present) async(4)
	{
		update_external_Hfield_half();  //RK
	}

#pragma acc parallel default(present) async(4)
	{
		//------------------Sublattice 1 ---------------------//
#pragma acc loop gang vector private(\
	Hx_anis_crt_AFM1, Hy_anis_crt_AFM1, Hz_anis_crt_AFM1,\
	Hx_anis_glb_AFM1, Hy_anis_glb_AFM1, Hz_anis_glb_AFM1,\
	Hx_elas_crt_AFM1, Hy_elas_crt_AFM1, Hz_elas_crt_AFM1,\
	Hx_elas_glb_AFM1, Hy_elas_glb_AFM1, Hz_elas_glb_AFM1,\
	Hx_exch_glb_AFM1, Hy_exch_glb_AFM1, Hz_exch_glb_AFM1,\
	Hx_iAFM_glb_AFM1, Hy_iAFM_glb_AFM1, Hz_iAFM_glb_AFM1,\
	Hx_eff_glb_AFM1, Hy_eff_glb_AFM1, Hz_eff_glb_AFM1,\
	precess_termx_AFM1, precess_termy_AFM1, precess_termz_AFM1,\
	damping_termx_AFM1, damping_termy_AFM1, damping_termz_AFM1,\
	x,y,z,mat_type,mat,\
	Ms_AFM1, K1_AFM1, K2_AFM1, B1_AFM1, B2_AFM1, \
	ux, uy, uz, udotm_AFM1, \
	JAFM_AFM1,\
	alpha_AFM1, gyro_AFM1, precess_AFM1,\
	mx_glb_local_AFM1, my_glb_local_AFM1, mz_glb_local_AFM1,\
	mx_AFM1, my_AFM1, mz_AFM1,\
	exx, eyy, ezz, eyz, exz, exy,\
	dmxdx_AFM1, dmydy_AFM1, dmzdx_AFM1, dmzdy_AFM1,\
	dmxdx2_f_AFM1, dmydx2_f_AFM1, dmzdx2_f_AFM1,\
	dmxdx2_b_AFM1, dmydx2_b_AFM1, dmzdx2_b_AFM1,\
	dmxdy2_f_AFM1, dmydy2_f_AFM1, dmzdy2_f_AFM1,\
	dmxdy2_b_AFM1, dmydy2_b_AFM1, dmzdy2_b_AFM1,\
	dmxdz2_f_AFM1, dmydz2_f_AFM1, dmzdz2_f_AFM1,\
	dmxdz2_b_AFM1, dmydz2_b_AFM1, dmzdz2_b_AFM1)
//Aexi, Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb, iDMIi)
		for (long int i = 0; i < n; i++) {
			mat_type = pt_glb->material_cell(i);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_AFM == true) {
					x = i / (ny * nz);
					y = (i - x * (ny * nz)) / nz;
					z = i - x * (ny * nz) - y * nz;

					ux = mat->uniaxis_AFM[0]; uy = mat->uniaxis_AFM[1]; uz = mat->uniaxis_AFM[2];

					if (pt_glb->if_elasto_on_pORm == true) {
						exx = pt_elas->exxt0_crt(i) + pt_elas->Dexx_crt(i);
						eyy = pt_elas->eyyt0_crt(i) + pt_elas->Deyy_crt(i);
						ezz = pt_elas->ezzt0_crt(i) + pt_elas->Dezz_crt(i);
						eyz = pt_elas->eyzt0_crt(i) + pt_elas->Deyz_crt(i);
						exz = pt_elas->exzt0_crt(i) + pt_elas->Dexz_crt(i);
						exy = pt_elas->exyt0_crt(i) + pt_elas->Dexy_crt(i);
					}
					else {
						exx = 0.;
						eyy = 0.;
						ezz = 0.;
						eyz = 0.;
						exz = 0.;
						exy = 0.;
					}

					mx_glb_local_AFM1 = mx_AFM1_glb_store(i);
					my_glb_local_AFM1 = my_AFM1_glb_store(i);
					mz_glb_local_AFM1 = mz_AFM1_glb_store(i);

					pt_math->transform_vector_glb2crt(mx_glb_local_AFM1, my_glb_local_AFM1, mz_glb_local_AFM1, \
						mx_AFM1, my_AFM1, mz_AFM1);

					Ms_AFM1 = mat->Ms_AFM1;
					K1_AFM1 = mat->K1_AFM1 / mu0 / Ms_AFM1; //K2_AFM1 = mat->K2 / mu0 / Ms;					
					udotm_AFM1 = mx_AFM1 * ux + my_AFM1 * uy + mz_AFM1 * uz;
					B1_AFM1 = mat->B1_AFM1 / mu0 / Ms_AFM1; B2_AFM1 = mat->B2_AFM1 / mu0 / Ms_AFM1;
					//Aexi_AFM1 = mat->Aex_AFM1; //iDMIi = mat->iDMI;
					JAFM_AFM1 = -1. * mat->J_AFM / mu0 / Ms_AFM1;

					get_m_AFM1_spatial_derivative_local \
						(x, y, z, \
							dmxdx_AFM1, dmydy_AFM1, dmzdx_AFM1, dmzdy_AFM1, \
							dmxdx2_f_AFM1, dmydx2_f_AFM1, dmzdx2_f_AFM1, \
							dmxdx2_b_AFM1, dmydx2_b_AFM1, dmzdx2_b_AFM1, \
							dmxdy2_f_AFM1, dmydy2_f_AFM1, dmzdy2_f_AFM1, \
							dmxdy2_b_AFM1, dmydy2_b_AFM1, dmzdy2_b_AFM1, \
							dmxdz2_f_AFM1, dmydz2_f_AFM1, dmzdz2_f_AFM1, \
							dmxdz2_b_AFM1, dmydz2_b_AFM1, dmzdz2_b_AFM1);

					//get_laplacian_m_local \
						//	(x, y, z, \
						//		d2mxdx2, d2mxdy2, d2mxdz2, \
						//		d2mydx2, d2mydy2, d2mydz2, \
						//		d2mzdx2, d2mzdy2, d2mzdz2);

						//------------Magnetocrystalline anisotropy field--------------//
					Hx_anis_crt_AFM1 = 2. * K1_AFM1 * udotm_AFM1 * ux;
					Hy_anis_crt_AFM1 = 2. * K1_AFM1 * udotm_AFM1 * uy;
					Hz_anis_crt_AFM1 = 2. * K1_AFM1 * udotm_AFM1 * uz;

					pt_math->transform_vector_crt2glb(Hx_anis_crt_AFM1, Hy_anis_crt_AFM1, Hz_anis_crt_AFM1, \
						Hx_anis_glb_AFM1, Hy_anis_glb_AFM1, Hz_anis_glb_AFM1);

					//Hx_anis(i) = Hx_anis_glb;
					//Hy_anis(i) = Hy_anis_glb;
					//Hz_anis(i) = Hz_anis_glb;

					//------------------Magnetoelastic effective field---------------//
					Hx_elas_crt_AFM1 = -1. * (B1_AFM1 * (exx)*mx_AFM1 + \
						B2_AFM1 * ((exy)*my_AFM1 + (exz)*mz_AFM1));
					Hy_elas_crt_AFM1 = -1. * (B1_AFM1 * (eyy)*my_AFM1 + \
						B2_AFM1 * ((exy)*mx_AFM1 + (eyz)*mz_AFM1));
					Hz_elas_crt_AFM1 = -1. * (B1_AFM1 * (ezz)*mz_AFM1 + \
						B2_AFM1 * ((exz)*mx_AFM1 + (eyz)*my_AFM1));

					pt_math->transform_vector_crt2glb(Hx_elas_crt_AFM1, Hy_elas_crt_AFM1, Hz_elas_crt_AFM1, \
						Hx_elas_glb_AFM1, Hy_elas_glb_AFM1, Hz_elas_glb_AFM1);

					//Hx_elas(i) = Hx_elas_glb;
					//Hy_elas(i) = Hy_elas_glb;
					//Hz_elas(i) = Hz_elas_glb;

					//-----------------Magnetic exchange field----------------------//
					//----------Averaged parameter for inter-region exchange - Approach 1---------//
					Hx_exch_glb_AFM1 = 2. / mu0 / Ms_AFM1 * \
						(Aij_AFM1_xf(i) * dmxdx2_f_AFM1 + Aij_AFM1_xb(i) * dmxdx2_b_AFM1 \
							+ Aij_AFM1_yf(i) * dmxdy2_f_AFM1 + Aij_AFM1_yb(i) * dmxdy2_b_AFM1 \
							+ Aij_AFM1_zf(i) * dmxdz2_f_AFM1 + Aij_AFM1_zb(i) * dmxdz2_b_AFM1);
					Hy_exch_glb_AFM1 = 2. / mu0 / Ms_AFM1 * \
						(Aij_AFM1_xf(i) * dmydx2_f_AFM1 + Aij_AFM1_xb(i) * dmydx2_b_AFM1 \
							+ Aij_AFM1_yf(i) * dmydy2_f_AFM1 + Aij_AFM1_yb(i) * dmydy2_b_AFM1 \
							+ Aij_AFM1_zf(i) * dmydz2_f_AFM1 + Aij_AFM1_zb(i) * dmydz2_b_AFM1);
					Hz_exch_glb_AFM1 = 2. / mu0 / Ms_AFM1 * \
						(Aij_AFM1_xf(i) * dmzdx2_f_AFM1 + Aij_AFM1_xb(i) * dmzdx2_b_AFM1 \
							+ Aij_AFM1_yf(i) * dmzdy2_f_AFM1 + Aij_AFM1_yb(i) * dmzdy2_b_AFM1 \
							+ Aij_AFM1_zf(i) * dmzdz2_f_AFM1 + Aij_AFM1_zb(i) * dmzdz2_b_AFM1);

					//Hx_exch(i) = 2. * Aex * (d2mxdx2 + d2mxdy2 + d2mxdz2);
					//Hy_exch(i) = 2. * Aex * (d2mydx2 + d2mydy2 + d2mydz2);
					//Hz_exch(i) = 2. * Aex * (d2mzdx2 + d2mzdy2 + d2mzdz2);
					// 
					// -----------------Inter-sublattice AFM exchange field---------------//
					Hx_iAFM_glb_AFM1 = JAFM_AFM1 * mx_AFM2_glb_store(i);
					Hy_iAFM_glb_AFM1 = JAFM_AFM1 * my_AFM2_glb_store(i);
					Hz_iAFM_glb_AFM1 = JAFM_AFM1 * mz_AFM2_glb_store(i);

					//------------------Interfacial DMI effective field-------------------//
					//Hx_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdx;
					//Hy_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdy;
					//Hz_iDMI_glb = -2. * iDMIi / mu0 / Ms * (dmxdx + dmydy);

					//-------------------LLG Equation------------------------------//
					Hx_eff_glb_AFM1 = Hx_stat(i) + \
						Hx_anis_glb_AFM1 + Hx_elas_glb_AFM1 + \
						Hx_exch_glb_AFM1 + Hx_iAFM_glb_AFM1 /*+ Hx_iDMI_glb*/ + pt_glb->Hext[0];

					Hy_eff_glb_AFM1 = Hy_stat(i) + \
						Hy_anis_glb_AFM1 + Hy_elas_glb_AFM1 + \
						Hy_exch_glb_AFM1 + Hy_iAFM_glb_AFM1/*+ Hy_iDMI_glb*/ + pt_glb->Hext[1];

					Hz_eff_glb_AFM1 = Hz_stat(i) + \
						Hz_anis_glb_AFM1 + Hz_elas_glb_AFM1 + \
						Hz_exch_glb_AFM1 + Hz_iAFM_glb_AFM1/*+ Hz_iDMI_glb*/ + pt_glb->Hext[2];

					if (pt_glb->if_EM_backaction_on_M == true) {
						Hx_eff_glb_AFM1 = Hx_eff_glb_AFM1 + pt_EM->DHx_em_cell(i);
						Hy_eff_glb_AFM1 = Hy_eff_glb_AFM1 + pt_EM->DHy_em_cell(i);
						Hz_eff_glb_AFM1 = Hz_eff_glb_AFM1 + pt_EM->DHz_em_cell(i);
					}

					alpha_AFM1 = mat->damping_AFM1;
					gyro_AFM1 = mat->gyro_AFM1;
					precess_AFM1 = -gyro_AFM1 / (1. + alpha_AFM1 * alpha_AFM1) * pt_glb->dt;

					precess_termx_AFM1 = my_glb_local_AFM1 * Hz_eff_glb_AFM1 - mz_glb_local_AFM1 * Hy_eff_glb_AFM1;
					precess_termy_AFM1 = mz_glb_local_AFM1 * Hx_eff_glb_AFM1 - mx_glb_local_AFM1 * Hz_eff_glb_AFM1;
					precess_termz_AFM1 = mx_glb_local_AFM1 * Hy_eff_glb_AFM1 - my_glb_local_AFM1 * Hx_eff_glb_AFM1;

					damping_termx_AFM1 = my_glb_local_AFM1 * precess_termz_AFM1 - mz_glb_local_AFM1 * precess_termy_AFM1;
					damping_termy_AFM1 = mz_glb_local_AFM1 * precess_termx_AFM1 - mx_glb_local_AFM1 * precess_termz_AFM1;
					damping_termz_AFM1 = mx_glb_local_AFM1 * precess_termy_AFM1 - my_glb_local_AFM1 * precess_termx_AFM1;

					dmx_AFM1_glb_rk2(i) = precess_AFM1 * precess_termx_AFM1 + alpha_AFM1 * precess_AFM1 * damping_termx_AFM1; //RK
					dmy_AFM1_glb_rk2(i) = precess_AFM1 * precess_termy_AFM1 + alpha_AFM1 * precess_AFM1 * damping_termy_AFM1; //RK
					dmz_AFM1_glb_rk2(i) = precess_AFM1 * precess_termz_AFM1 + alpha_AFM1 * precess_AFM1 * damping_termz_AFM1; //RK
				}
			}
		}

		//------------------Sublattice 2 ---------------------//
#pragma acc loop gang vector private(\
	Hx_anis_crt_AFM2, Hy_anis_crt_AFM2, Hz_anis_crt_AFM2,\
	Hx_anis_glb_AFM2, Hy_anis_glb_AFM2, Hz_anis_glb_AFM2,\
	Hx_elas_crt_AFM2, Hy_elas_crt_AFM2, Hz_elas_crt_AFM2,\
	Hx_elas_glb_AFM2, Hy_elas_glb_AFM2, Hz_elas_glb_AFM2,\
	Hx_exch_glb_AFM2, Hy_exch_glb_AFM2, Hz_exch_glb_AFM2,\
	Hx_iAFM_glb_AFM2, Hy_iAFM_glb_AFM2, Hz_iAFM_glb_AFM2,\
	Hx_eff_glb_AFM2, Hy_eff_glb_AFM2, Hz_eff_glb_AFM2,\
	precess_termx_AFM2, precess_termy_AFM2, precess_termz_AFM2,\
	damping_termx_AFM2, damping_termy_AFM2, damping_termz_AFM2,\
	x,y,z,mat_type,mat,\
	Ms_AFM2, K1_AFM2, K2_AFM2, B1_AFM2, B2_AFM2, \
	ux, uy, uz, udotm_AFM2, \
	JAFM_AFM2,\
	alpha_AFM2, gyro_AFM2, precess_AFM2,\
	mx_glb_local_AFM2, my_glb_local_AFM2, mz_glb_local_AFM2,\
	mx_AFM2, my_AFM2, mz_AFM2,\
	exx, eyy, ezz, eyz, exz, exy,\
	dmxdx_AFM2, dmydy_AFM2, dmzdx_AFM2, dmzdy_AFM2,\
	dmxdx2_f_AFM2, dmydx2_f_AFM2, dmzdx2_f_AFM2,\
	dmxdx2_b_AFM2, dmydx2_b_AFM2, dmzdx2_b_AFM2,\
	dmxdy2_f_AFM2, dmydy2_f_AFM2, dmzdy2_f_AFM2,\
	dmxdy2_b_AFM2, dmydy2_b_AFM2, dmzdy2_b_AFM2,\
	dmxdz2_f_AFM2, dmydz2_f_AFM2, dmzdz2_f_AFM2,\
	dmxdz2_b_AFM2, dmydz2_b_AFM2, dmzdz2_b_AFM2)
	//Aexi, Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb, iDMIi)
		for (long int i = 0; i < n; i++) {
			mat_type = pt_glb->material_cell(i);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_AFM == true) {
					x = i / (ny * nz);
					y = (i - x * (ny * nz)) / nz;
					z = i - x * (ny * nz) - y * nz;

					ux = mat->uniaxis_AFM[0]; uy = mat->uniaxis_AFM[1]; uz = mat->uniaxis_AFM[2];

					if (pt_glb->if_elasto_on_pORm == true) {
						exx = pt_elas->exxt0_crt(i) + pt_elas->Dexx_crt(i);
						eyy = pt_elas->eyyt0_crt(i) + pt_elas->Deyy_crt(i);
						ezz = pt_elas->ezzt0_crt(i) + pt_elas->Dezz_crt(i);
						eyz = pt_elas->eyzt0_crt(i) + pt_elas->Deyz_crt(i);
						exz = pt_elas->exzt0_crt(i) + pt_elas->Dexz_crt(i);
						exy = pt_elas->exyt0_crt(i) + pt_elas->Dexy_crt(i);
					}
					else {
						exx = 0.;
						eyy = 0.;
						ezz = 0.;
						eyz = 0.;
						exz = 0.;
						exy = 0.;
					}

					mx_glb_local_AFM2 = mx_AFM2_glb_store(i);
					my_glb_local_AFM2 = my_AFM2_glb_store(i);
					mz_glb_local_AFM2 = mz_AFM2_glb_store(i);

					pt_math->transform_vector_glb2crt(mx_glb_local_AFM2, my_glb_local_AFM2, mz_glb_local_AFM2, \
						mx_AFM2, my_AFM2, mz_AFM2);

					Ms_AFM2 = mat->Ms_AFM2;
					K1_AFM2 = mat->K1_AFM2 / mu0 / Ms_AFM2; //K2 = mat->K2_AFM2 / mu0 / Ms;					
					udotm_AFM2 = mx_AFM2 * ux + my_AFM2 * uy + mz_AFM2 * uz;
					B1_AFM2 = mat->B1_AFM2 / mu0 / Ms_AFM2; B2_AFM2 = mat->B2_AFM2 / mu0 / Ms_AFM2;
					//Aexi_AFM2 = mat->Aex_AFM2; //iDMIi = mat->iDMI;
					JAFM_AFM2 = -1. * mat->J_AFM / mu0 / Ms_AFM2;

					get_m_AFM2_spatial_derivative_local \
						(x, y, z, \
							dmxdx_AFM2, dmydy_AFM2, dmzdx_AFM2, dmzdy_AFM2, \
							dmxdx2_f_AFM2, dmydx2_f_AFM2, dmzdx2_f_AFM2, \
							dmxdx2_b_AFM2, dmydx2_b_AFM2, dmzdx2_b_AFM2, \
							dmxdy2_f_AFM2, dmydy2_f_AFM2, dmzdy2_f_AFM2, \
							dmxdy2_b_AFM2, dmydy2_b_AFM2, dmzdy2_b_AFM2, \
							dmxdz2_f_AFM2, dmydz2_f_AFM2, dmzdz2_f_AFM2, \
							dmxdz2_b_AFM2, dmydz2_b_AFM2, dmzdz2_b_AFM2);

					//get_laplacian_m_local \
						//	(x, y, z, \
						//		d2mxdx2, d2mxdy2, d2mxdz2, \
						//		d2mydx2, d2mydy2, d2mydz2, \
						//		d2mzdx2, d2mzdy2, d2mzdz2);

						//------------Magnetocrystalline anisotropy field--------------//
					Hx_anis_crt_AFM2 = 2. * K1_AFM2 * udotm_AFM2 * ux;
					Hy_anis_crt_AFM2 = 2. * K1_AFM2 * udotm_AFM2 * uy;
					Hz_anis_crt_AFM2 = 2. * K1_AFM2 * udotm_AFM2 * uz;

					pt_math->transform_vector_crt2glb(Hx_anis_crt_AFM2, Hy_anis_crt_AFM2, Hz_anis_crt_AFM2, \
						Hx_anis_glb_AFM2, Hy_anis_glb_AFM2, Hz_anis_glb_AFM2);

					//Hx_anis(i) = Hx_anis_glb;
					//Hy_anis(i) = Hy_anis_glb;
					//Hz_anis(i) = Hz_anis_glb;

					//------------------Magnetoelastic effective field---------------//
					Hx_elas_crt_AFM2 = -1. * (B1_AFM2 * (exx)*mx_AFM2 + \
						B2_AFM2 * ((exy)*my_AFM2 + (exz)*mz_AFM2));
					Hy_elas_crt_AFM2 = -1. * (B1_AFM2 * (eyy)*my_AFM2 + \
						B2_AFM2 * ((exy)*mx_AFM2 + (eyz)*mz_AFM2));
					Hz_elas_crt_AFM2 = -1. * (B1_AFM2 * (ezz)*mz_AFM2 + \
						B2_AFM2 * ((exz)*mx_AFM2 + (eyz)*my_AFM2));

					pt_math->transform_vector_crt2glb(Hx_elas_crt_AFM2, Hy_elas_crt_AFM2, Hz_elas_crt_AFM2, \
						Hx_elas_glb_AFM2, Hy_elas_glb_AFM2, Hz_elas_glb_AFM2);

					//Hx_elas(i) = Hx_elas_glb;
					//Hy_elas(i) = Hy_elas_glb;
					//Hz_elas(i) = Hz_elas_glb;

					//-----------------Magnetic exchange field----------------------//
					//----------Averaged parameter for inter-region exchange - Approach 1---------//
					Hx_exch_glb_AFM2 = 2. / mu0 / Ms_AFM2 * \
						(Aij_AFM2_xf(i) * dmxdx2_f_AFM2 + Aij_AFM2_xb(i) * dmxdx2_b_AFM2 \
							+ Aij_AFM2_yf(i) * dmxdy2_f_AFM2 + Aij_AFM2_yb(i) * dmxdy2_b_AFM2 \
							+ Aij_AFM2_zf(i) * dmxdz2_f_AFM2 + Aij_AFM2_zb(i) * dmxdz2_b_AFM2);
					Hy_exch_glb_AFM2 = 2. / mu0 / Ms_AFM2 * \
						(Aij_AFM2_xf(i) * dmydx2_f_AFM2 + Aij_AFM2_xb(i) * dmydx2_b_AFM2 \
							+ Aij_AFM2_yf(i) * dmydy2_f_AFM2 + Aij_AFM2_yb(i) * dmydy2_b_AFM2 \
							+ Aij_AFM2_zf(i) * dmydz2_f_AFM2 + Aij_AFM2_zb(i) * dmydz2_b_AFM2);
					Hz_exch_glb_AFM2 = 2. / mu0 / Ms_AFM2 * \
						(Aij_AFM2_xf(i) * dmzdx2_f_AFM2 + Aij_AFM2_xb(i) * dmzdx2_b_AFM2 \
							+ Aij_AFM2_yf(i) * dmzdy2_f_AFM2 + Aij_AFM2_yb(i) * dmzdy2_b_AFM2 \
							+ Aij_AFM2_zf(i) * dmzdz2_f_AFM2 + Aij_AFM2_zb(i) * dmzdz2_b_AFM2);

					//Hx_exch(i) = 2. * Aex * (d2mxdx2 + d2mxdy2 + d2mxdz2);
					//Hy_exch(i) = 2. * Aex * (d2mydx2 + d2mydy2 + d2mydz2);
					//Hz_exch(i) = 2. * Aex * (d2mzdx2 + d2mzdy2 + d2mzdz2);
					// 
					// -----------------Inter-sublattice AFM exchange field---------------//
					Hx_iAFM_glb_AFM2 = JAFM_AFM2 * mx_AFM1_glb_store(i);
					Hy_iAFM_glb_AFM2 = JAFM_AFM2 * my_AFM1_glb_store(i);
					Hz_iAFM_glb_AFM2 = JAFM_AFM2 * mz_AFM1_glb_store(i);
					//------------------Interfacial DMI effective field-------------------//
					//Hx_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdx;
					//Hy_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdy;
					//Hz_iDMI_glb = -2. * iDMIi / mu0 / Ms * (dmxdx + dmydy);

					//-------------------LLG Equation------------------------------//
					Hx_eff_glb_AFM2 = Hx_stat(i) + \
						Hx_anis_glb_AFM2 + Hx_elas_glb_AFM2 + \
						Hx_exch_glb_AFM2 + Hx_iAFM_glb_AFM2 /*+ Hx_iDMI_glb*/ + pt_glb->Hext[0];

					Hy_eff_glb_AFM2 = Hy_stat(i) + \
						Hy_anis_glb_AFM2 + Hy_elas_glb_AFM2 + \
						Hy_exch_glb_AFM2 + Hy_iAFM_glb_AFM2 /*+ Hy_iDMI_glb*/ + pt_glb->Hext[1];

					Hz_eff_glb_AFM2 = Hz_stat(i) + \
						Hz_anis_glb_AFM2 + Hz_elas_glb_AFM2 + \
						Hz_exch_glb_AFM2 + Hz_iAFM_glb_AFM2 /*+ Hz_iDMI_glb*/ + pt_glb->Hext[2];

					if (pt_glb->if_EM_backaction_on_M == true) {
						Hx_eff_glb_AFM2 = Hx_eff_glb_AFM2 + pt_EM->DHx_em_cell(i);
						Hy_eff_glb_AFM2 = Hy_eff_glb_AFM2 + pt_EM->DHy_em_cell(i);
						Hz_eff_glb_AFM2 = Hz_eff_glb_AFM2 + pt_EM->DHz_em_cell(i);
					}

					alpha_AFM2 = mat->damping_AFM2;
					gyro_AFM2 = mat->gyro_AFM2;
					precess_AFM2 = -gyro_AFM2 / (1 + alpha_AFM2 * alpha_AFM2) * pt_glb->dt;

					precess_termx_AFM2 = my_glb_local_AFM2 * Hz_eff_glb_AFM2 - mz_glb_local_AFM2 * Hy_eff_glb_AFM2;
					precess_termy_AFM2 = mz_glb_local_AFM2 * Hx_eff_glb_AFM2 - mx_glb_local_AFM2 * Hz_eff_glb_AFM2;
					precess_termz_AFM2 = mx_glb_local_AFM2 * Hy_eff_glb_AFM2 - my_glb_local_AFM2 * Hx_eff_glb_AFM2;

					damping_termx_AFM2 = my_glb_local_AFM2 * precess_termz_AFM2 - mz_glb_local_AFM2 * precess_termy_AFM2;
					damping_termy_AFM2 = mz_glb_local_AFM2 * precess_termx_AFM2 - mx_glb_local_AFM2 * precess_termz_AFM2;
					damping_termz_AFM2 = mx_glb_local_AFM2 * precess_termy_AFM2 - my_glb_local_AFM2 * precess_termx_AFM2;

					dmx_AFM2_glb_rk2(i) = precess_AFM2 * precess_termx_AFM2 + alpha_AFM2 * precess_AFM2 * damping_termx_AFM2; //RK
					dmy_AFM2_glb_rk2(i) = precess_AFM2 * precess_termy_AFM2 + alpha_AFM2 * precess_AFM2 * damping_termy_AFM2; //RK
					dmz_AFM2_glb_rk2(i) = precess_AFM2 * precess_termz_AFM2 + alpha_AFM2 * precess_AFM2 * damping_termz_AFM2; //RK
				}
			}
		}
	}

	if (pt_glb->if_spin_pumping == true) {
		get_J_ISHE_RK2(); //RK
	}
}

void magnetic_system::get_dm_AFM_RK3() { //RK
	double Hx_anis_crt_AFM1, Hy_anis_crt_AFM1, Hz_anis_crt_AFM1;
	double Hx_anis_glb_AFM1, Hy_anis_glb_AFM1, Hz_anis_glb_AFM1;
	double Hx_elas_crt_AFM1, Hy_elas_crt_AFM1, Hz_elas_crt_AFM1;
	double Hx_elas_glb_AFM1, Hy_elas_glb_AFM1, Hz_elas_glb_AFM1;
	double Hx_exch_glb_AFM1, Hy_exch_glb_AFM1, Hz_exch_glb_AFM1;
	double Hx_iAFM_glb_AFM1, Hy_iAFM_glb_AFM1, Hz_iAFM_glb_AFM1;
	double Hx_eff_glb_AFM1, Hy_eff_glb_AFM1, Hz_eff_glb_AFM1;
	double precess_termx_AFM1, precess_termy_AFM1, precess_termz_AFM1;
	double damping_termx_AFM1, damping_termy_AFM1, damping_termz_AFM1;
	//double Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb;

	double Hx_anis_crt_AFM2, Hy_anis_crt_AFM2, Hz_anis_crt_AFM2;
	double Hx_anis_glb_AFM2, Hy_anis_glb_AFM2, Hz_anis_glb_AFM2;
	double Hx_elas_crt_AFM2, Hy_elas_crt_AFM2, Hz_elas_crt_AFM2;
	double Hx_elas_glb_AFM2, Hy_elas_glb_AFM2, Hz_elas_glb_AFM2;
	double Hx_exch_glb_AFM2, Hy_exch_glb_AFM2, Hz_exch_glb_AFM2;
	double Hx_iAFM_glb_AFM2, Hy_iAFM_glb_AFM2, Hz_iAFM_glb_AFM2;
	double Hx_eff_glb_AFM2, Hy_eff_glb_AFM2, Hz_eff_glb_AFM2;
	double precess_termx_AFM2, precess_termy_AFM2, precess_termz_AFM2;
	double damping_termx_AFM2, damping_termy_AFM2, damping_termz_AFM2;
	//double Hx_iDMI_glb_AFM2, Hy_iDMI_glb_AFM2, Hz_iDMI_glb_AFM2;

	long int x, y, z;
	unsigned int mat_type;
	material* mat;
	double Ms_AFM1, K1_AFM1, K2_AFM1, B1_AFM1, B2_AFM1;
	double Ms_AFM2, K1_AFM2, K2_AFM2, B1_AFM2, B2_AFM2;
	double ux, uy, uz;
	double udotm_AFM1, udotm_AFM2;
	//double Aexi;//, iDMIi;
	double JAFM_AFM1, JAFM_AFM2;
	double alpha_AFM1, gyro_AFM1, precess_AFM1;
	double alpha_AFM2, gyro_AFM2, precess_AFM2;
	double mx_glb_local_AFM1, my_glb_local_AFM1, mz_glb_local_AFM1;
	double mx_glb_local_AFM2, my_glb_local_AFM2, mz_glb_local_AFM2;
	double mx_AFM1, my_AFM1, mz_AFM1;
	double mx_AFM2, my_AFM2, mz_AFM2;
	double exx, eyy, ezz, eyz, exz, exy;
	//double d2mxdx2, d2mxdy2, d2mxdz2;
	//double d2mydx2, d2mydy2, d2mydz2;
	//double d2mzdx2, d2mzdy2, d2mzdz2;
	double dmxdx_AFM1, dmydy_AFM1, dmzdx_AFM1, dmzdy_AFM1;
	double dmxdx2_f_AFM1, dmydx2_f_AFM1, dmzdx2_f_AFM1;
	double dmxdx2_b_AFM1, dmydx2_b_AFM1, dmzdx2_b_AFM1;
	double dmxdy2_f_AFM1, dmydy2_f_AFM1, dmzdy2_f_AFM1;
	double dmxdy2_b_AFM1, dmydy2_b_AFM1, dmzdy2_b_AFM1;
	double dmxdz2_f_AFM1, dmydz2_f_AFM1, dmzdz2_f_AFM1;
	double dmxdz2_b_AFM1, dmydz2_b_AFM1, dmzdz2_b_AFM1;

	double dmxdx_AFM2, dmydy_AFM2, dmzdx_AFM2, dmzdy_AFM2;
	double dmxdx2_f_AFM2, dmydx2_f_AFM2, dmzdx2_f_AFM2;
	double dmxdx2_b_AFM2, dmydx2_b_AFM2, dmzdx2_b_AFM2;
	double dmxdy2_f_AFM2, dmydy2_f_AFM2, dmzdy2_f_AFM2;
	double dmxdy2_b_AFM2, dmydy2_b_AFM2, dmzdy2_b_AFM2;
	double dmxdz2_f_AFM2, dmydz2_f_AFM2, dmzdz2_f_AFM2;
	double dmxdz2_b_AFM2, dmydz2_b_AFM2, dmzdz2_b_AFM2;

#pragma acc serial default(present) async(4)
	{
		update_external_Hfield_half();  //RK
	}

#pragma acc parallel default(present) async(4)
	{
		//------------------Sublattice 1 ---------------------//
#pragma acc loop gang vector private(\
	Hx_anis_crt_AFM1, Hy_anis_crt_AFM1, Hz_anis_crt_AFM1,\
	Hx_anis_glb_AFM1, Hy_anis_glb_AFM1, Hz_anis_glb_AFM1,\
	Hx_elas_crt_AFM1, Hy_elas_crt_AFM1, Hz_elas_crt_AFM1,\
	Hx_elas_glb_AFM1, Hy_elas_glb_AFM1, Hz_elas_glb_AFM1,\
	Hx_exch_glb_AFM1, Hy_exch_glb_AFM1, Hz_exch_glb_AFM1,\
	Hx_iAFM_glb_AFM1, Hy_iAFM_glb_AFM1, Hz_iAFM_glb_AFM1,\
	Hx_eff_glb_AFM1, Hy_eff_glb_AFM1, Hz_eff_glb_AFM1,\
	precess_termx_AFM1, precess_termy_AFM1, precess_termz_AFM1,\
	damping_termx_AFM1, damping_termy_AFM1, damping_termz_AFM1,\
	x,y,z,mat_type,mat,\
	Ms_AFM1, K1_AFM1, K2_AFM1, B1_AFM1, B2_AFM1, \
	ux, uy, uz, udotm_AFM1, \
	JAFM_AFM1,\
	alpha_AFM1, gyro_AFM1, precess_AFM1,\
	mx_glb_local_AFM1, my_glb_local_AFM1, mz_glb_local_AFM1,\
	mx_AFM1, my_AFM1, mz_AFM1,\
	exx, eyy, ezz, eyz, exz, exy,\
	dmxdx_AFM1, dmydy_AFM1, dmzdx_AFM1, dmzdy_AFM1,\
	dmxdx2_f_AFM1, dmydx2_f_AFM1, dmzdx2_f_AFM1,\
	dmxdx2_b_AFM1, dmydx2_b_AFM1, dmzdx2_b_AFM1,\
	dmxdy2_f_AFM1, dmydy2_f_AFM1, dmzdy2_f_AFM1,\
	dmxdy2_b_AFM1, dmydy2_b_AFM1, dmzdy2_b_AFM1,\
	dmxdz2_f_AFM1, dmydz2_f_AFM1, dmzdz2_f_AFM1,\
	dmxdz2_b_AFM1, dmydz2_b_AFM1, dmzdz2_b_AFM1)
//Aexi, Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb, iDMIi)
		for (long int i = 0; i < n; i++) {
			mat_type = pt_glb->material_cell(i);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_AFM == true) {
					x = i / (ny * nz);
					y = (i - x * (ny * nz)) / nz;
					z = i - x * (ny * nz) - y * nz;

					ux = mat->uniaxis_AFM[0]; uy = mat->uniaxis_AFM[1]; uz = mat->uniaxis_AFM[2];

					if (pt_glb->if_elasto_on_pORm == true) {
						exx = pt_elas->exxt0_crt(i) + pt_elas->Dexx_crt(i);
						eyy = pt_elas->eyyt0_crt(i) + pt_elas->Deyy_crt(i);
						ezz = pt_elas->ezzt0_crt(i) + pt_elas->Dezz_crt(i);
						eyz = pt_elas->eyzt0_crt(i) + pt_elas->Deyz_crt(i);
						exz = pt_elas->exzt0_crt(i) + pt_elas->Dexz_crt(i);
						exy = pt_elas->exyt0_crt(i) + pt_elas->Dexy_crt(i);
					}
					else {
						exx = 0.;
						eyy = 0.;
						ezz = 0.;
						eyz = 0.;
						exz = 0.;
						exy = 0.;
					}

					mx_glb_local_AFM1 = mx_AFM1_glb_store(i);
					my_glb_local_AFM1 = my_AFM1_glb_store(i);
					mz_glb_local_AFM1 = mz_AFM1_glb_store(i);

					pt_math->transform_vector_glb2crt(mx_glb_local_AFM1, my_glb_local_AFM1, mz_glb_local_AFM1, \
						mx_AFM1, my_AFM1, mz_AFM1);

					Ms_AFM1 = mat->Ms_AFM1;
					K1_AFM1 = mat->K1_AFM1 / mu0 / Ms_AFM1; //K2_AFM1 = mat->K2 / mu0 / Ms;					
					udotm_AFM1 = mx_AFM1 * ux + my_AFM1 * uy + mz_AFM1 * uz;
					B1_AFM1 = mat->B1_AFM1 / mu0 / Ms_AFM1; B2_AFM1 = mat->B2_AFM1 / mu0 / Ms_AFM1;
					//Aexi_AFM1 = mat->Aex_AFM1; //iDMIi = mat->iDMI;
					JAFM_AFM1 = -1. * mat->J_AFM / mu0 / Ms_AFM1;

					get_m_AFM1_spatial_derivative_local \
						(x, y, z, \
							dmxdx_AFM1, dmydy_AFM1, dmzdx_AFM1, dmzdy_AFM1, \
							dmxdx2_f_AFM1, dmydx2_f_AFM1, dmzdx2_f_AFM1, \
							dmxdx2_b_AFM1, dmydx2_b_AFM1, dmzdx2_b_AFM1, \
							dmxdy2_f_AFM1, dmydy2_f_AFM1, dmzdy2_f_AFM1, \
							dmxdy2_b_AFM1, dmydy2_b_AFM1, dmzdy2_b_AFM1, \
							dmxdz2_f_AFM1, dmydz2_f_AFM1, dmzdz2_f_AFM1, \
							dmxdz2_b_AFM1, dmydz2_b_AFM1, dmzdz2_b_AFM1);

					//get_laplacian_m_local \
						//	(x, y, z, \
						//		d2mxdx2, d2mxdy2, d2mxdz2, \
						//		d2mydx2, d2mydy2, d2mydz2, \
						//		d2mzdx2, d2mzdy2, d2mzdz2);

						//------------Magnetocrystalline anisotropy field--------------//
					Hx_anis_crt_AFM1 = 2. * K1_AFM1 * udotm_AFM1 * ux;
					Hy_anis_crt_AFM1 = 2. * K1_AFM1 * udotm_AFM1 * uy;
					Hz_anis_crt_AFM1 = 2. * K1_AFM1 * udotm_AFM1 * uz;

					pt_math->transform_vector_crt2glb(Hx_anis_crt_AFM1, Hy_anis_crt_AFM1, Hz_anis_crt_AFM1, \
						Hx_anis_glb_AFM1, Hy_anis_glb_AFM1, Hz_anis_glb_AFM1);

					//Hx_anis(i) = Hx_anis_glb;
					//Hy_anis(i) = Hy_anis_glb;
					//Hz_anis(i) = Hz_anis_glb;

					//------------------Magnetoelastic effective field---------------//
					Hx_elas_crt_AFM1 = -1. * (B1_AFM1 * (exx)*mx_AFM1 + \
						B2_AFM1 * ((exy)*my_AFM1 + (exz)*mz_AFM1));
					Hy_elas_crt_AFM1 = -1. * (B1_AFM1 * (eyy)*my_AFM1 + \
						B2_AFM1 * ((exy)*mx_AFM1 + (eyz)*mz_AFM1));
					Hz_elas_crt_AFM1 = -1. * (B1_AFM1 * (ezz)*mz_AFM1 + \
						B2_AFM1 * ((exz)*mx_AFM1 + (eyz)*my_AFM1));

					pt_math->transform_vector_crt2glb(Hx_elas_crt_AFM1, Hy_elas_crt_AFM1, Hz_elas_crt_AFM1, \
						Hx_elas_glb_AFM1, Hy_elas_glb_AFM1, Hz_elas_glb_AFM1);

					//Hx_elas(i) = Hx_elas_glb;
					//Hy_elas(i) = Hy_elas_glb;
					//Hz_elas(i) = Hz_elas_glb;

					//-----------------Magnetic exchange field----------------------//
					//----------Averaged parameter for inter-region exchange - Approach 1---------//
					Hx_exch_glb_AFM1 = 2. / mu0 / Ms_AFM1 * \
						(Aij_AFM1_xf(i) * dmxdx2_f_AFM1 + Aij_AFM1_xb(i) * dmxdx2_b_AFM1 \
							+ Aij_AFM1_yf(i) * dmxdy2_f_AFM1 + Aij_AFM1_yb(i) * dmxdy2_b_AFM1 \
							+ Aij_AFM1_zf(i) * dmxdz2_f_AFM1 + Aij_AFM1_zb(i) * dmxdz2_b_AFM1);
					Hy_exch_glb_AFM1 = 2. / mu0 / Ms_AFM1 * \
						(Aij_AFM1_xf(i) * dmydx2_f_AFM1 + Aij_AFM1_xb(i) * dmydx2_b_AFM1 \
							+ Aij_AFM1_yf(i) * dmydy2_f_AFM1 + Aij_AFM1_yb(i) * dmydy2_b_AFM1 \
							+ Aij_AFM1_zf(i) * dmydz2_f_AFM1 + Aij_AFM1_zb(i) * dmydz2_b_AFM1);
					Hz_exch_glb_AFM1 = 2. / mu0 / Ms_AFM1 * \
						(Aij_AFM1_xf(i) * dmzdx2_f_AFM1 + Aij_AFM1_xb(i) * dmzdx2_b_AFM1 \
							+ Aij_AFM1_yf(i) * dmzdy2_f_AFM1 + Aij_AFM1_yb(i) * dmzdy2_b_AFM1 \
							+ Aij_AFM1_zf(i) * dmzdz2_f_AFM1 + Aij_AFM1_zb(i) * dmzdz2_b_AFM1);

					//Hx_exch(i) = 2. * Aex * (d2mxdx2 + d2mxdy2 + d2mxdz2);
					//Hy_exch(i) = 2. * Aex * (d2mydx2 + d2mydy2 + d2mydz2);
					//Hz_exch(i) = 2. * Aex * (d2mzdx2 + d2mzdy2 + d2mzdz2);
					// 
					// -----------------Inter-sublattice AFM exchange field---------------//
					Hx_iAFM_glb_AFM1 = JAFM_AFM1 * mx_AFM2_glb_store(i);
					Hy_iAFM_glb_AFM1 = JAFM_AFM1 * my_AFM2_glb_store(i);
					Hz_iAFM_glb_AFM1 = JAFM_AFM1 * mz_AFM2_glb_store(i);

					//------------------Interfacial DMI effective field-------------------//
					//Hx_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdx;
					//Hy_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdy;
					//Hz_iDMI_glb = -2. * iDMIi / mu0 / Ms * (dmxdx + dmydy);

					//-------------------LLG Equation------------------------------//
					Hx_eff_glb_AFM1 = Hx_stat(i) + \
						Hx_anis_glb_AFM1 + Hx_elas_glb_AFM1 + \
						Hx_exch_glb_AFM1 + Hx_iAFM_glb_AFM1 /*+ Hx_iDMI_glb*/ + pt_glb->Hext[0];

					Hy_eff_glb_AFM1 = Hy_stat(i) + \
						Hy_anis_glb_AFM1 + Hy_elas_glb_AFM1 + \
						Hy_exch_glb_AFM1 + Hy_iAFM_glb_AFM1/*+ Hy_iDMI_glb*/ + pt_glb->Hext[1];

					Hz_eff_glb_AFM1 = Hz_stat(i) + \
						Hz_anis_glb_AFM1 + Hz_elas_glb_AFM1 + \
						Hz_exch_glb_AFM1 + Hz_iAFM_glb_AFM1/*+ Hz_iDMI_glb*/ + pt_glb->Hext[2];

					if (pt_glb->if_EM_backaction_on_M == true) {
						Hx_eff_glb_AFM1 = Hx_eff_glb_AFM1 + pt_EM->DHx_em_cell(i);
						Hy_eff_glb_AFM1 = Hy_eff_glb_AFM1 + pt_EM->DHy_em_cell(i);
						Hz_eff_glb_AFM1 = Hz_eff_glb_AFM1 + pt_EM->DHz_em_cell(i);
					}

					alpha_AFM1 = mat->damping_AFM1;
					gyro_AFM1 = mat->gyro_AFM1;
					precess_AFM1 = -gyro_AFM1 / (1. + alpha_AFM1 * alpha_AFM1) * pt_glb->dt;

					precess_termx_AFM1 = my_glb_local_AFM1 * Hz_eff_glb_AFM1 - mz_glb_local_AFM1 * Hy_eff_glb_AFM1;
					precess_termy_AFM1 = mz_glb_local_AFM1 * Hx_eff_glb_AFM1 - mx_glb_local_AFM1 * Hz_eff_glb_AFM1;
					precess_termz_AFM1 = mx_glb_local_AFM1 * Hy_eff_glb_AFM1 - my_glb_local_AFM1 * Hx_eff_glb_AFM1;

					damping_termx_AFM1 = my_glb_local_AFM1 * precess_termz_AFM1 - mz_glb_local_AFM1 * precess_termy_AFM1;
					damping_termy_AFM1 = mz_glb_local_AFM1 * precess_termx_AFM1 - mx_glb_local_AFM1 * precess_termz_AFM1;
					damping_termz_AFM1 = mx_glb_local_AFM1 * precess_termy_AFM1 - my_glb_local_AFM1 * precess_termx_AFM1;

					dmx_AFM1_glb_rk3(i) = precess_AFM1 * precess_termx_AFM1 + alpha_AFM1 * precess_AFM1 * damping_termx_AFM1; //RK
					dmy_AFM1_glb_rk3(i) = precess_AFM1 * precess_termy_AFM1 + alpha_AFM1 * precess_AFM1 * damping_termy_AFM1; //RK
					dmz_AFM1_glb_rk3(i) = precess_AFM1 * precess_termz_AFM1 + alpha_AFM1 * precess_AFM1 * damping_termz_AFM1; //RK
				}
			}
		}

		//------------------Sublattice 2 ---------------------//
#pragma acc loop gang vector private(\
	Hx_anis_crt_AFM2, Hy_anis_crt_AFM2, Hz_anis_crt_AFM2,\
	Hx_anis_glb_AFM2, Hy_anis_glb_AFM2, Hz_anis_glb_AFM2,\
	Hx_elas_crt_AFM2, Hy_elas_crt_AFM2, Hz_elas_crt_AFM2,\
	Hx_elas_glb_AFM2, Hy_elas_glb_AFM2, Hz_elas_glb_AFM2,\
	Hx_exch_glb_AFM2, Hy_exch_glb_AFM2, Hz_exch_glb_AFM2,\
	Hx_iAFM_glb_AFM2, Hy_iAFM_glb_AFM2, Hz_iAFM_glb_AFM2,\
	Hx_eff_glb_AFM2, Hy_eff_glb_AFM2, Hz_eff_glb_AFM2,\
	precess_termx_AFM2, precess_termy_AFM2, precess_termz_AFM2,\
	damping_termx_AFM2, damping_termy_AFM2, damping_termz_AFM2,\
	x,y,z,mat_type,mat,\
	Ms_AFM2, K1_AFM2, K2_AFM2, B1_AFM2, B2_AFM2, \
	ux, uy, uz, udotm_AFM2, \
	JAFM_AFM2,\
	alpha_AFM2, gyro_AFM2, precess_AFM2,\
	mx_glb_local_AFM2, my_glb_local_AFM2, mz_glb_local_AFM2,\
	mx_AFM2, my_AFM2, mz_AFM2,\
	exx, eyy, ezz, eyz, exz, exy,\
	dmxdx_AFM2, dmydy_AFM2, dmzdx_AFM2, dmzdy_AFM2,\
	dmxdx2_f_AFM2, dmydx2_f_AFM2, dmzdx2_f_AFM2,\
	dmxdx2_b_AFM2, dmydx2_b_AFM2, dmzdx2_b_AFM2,\
	dmxdy2_f_AFM2, dmydy2_f_AFM2, dmzdy2_f_AFM2,\
	dmxdy2_b_AFM2, dmydy2_b_AFM2, dmzdy2_b_AFM2,\
	dmxdz2_f_AFM2, dmydz2_f_AFM2, dmzdz2_f_AFM2,\
	dmxdz2_b_AFM2, dmydz2_b_AFM2, dmzdz2_b_AFM2)
	//Aexi, Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb, iDMIi)
		for (long int i = 0; i < n; i++) {
			mat_type = pt_glb->material_cell(i);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_AFM == true) {
					x = i / (ny * nz);
					y = (i - x * (ny * nz)) / nz;
					z = i - x * (ny * nz) - y * nz;

					ux = mat->uniaxis_AFM[0]; uy = mat->uniaxis_AFM[1]; uz = mat->uniaxis_AFM[2];

					if (pt_glb->if_elasto_on_pORm == true) {
						exx = pt_elas->exxt0_crt(i) + pt_elas->Dexx_crt(i);
						eyy = pt_elas->eyyt0_crt(i) + pt_elas->Deyy_crt(i);
						ezz = pt_elas->ezzt0_crt(i) + pt_elas->Dezz_crt(i);
						eyz = pt_elas->eyzt0_crt(i) + pt_elas->Deyz_crt(i);
						exz = pt_elas->exzt0_crt(i) + pt_elas->Dexz_crt(i);
						exy = pt_elas->exyt0_crt(i) + pt_elas->Dexy_crt(i);
					}
					else {
						exx = 0.;
						eyy = 0.;
						ezz = 0.;
						eyz = 0.;
						exz = 0.;
						exy = 0.;
					}

					mx_glb_local_AFM2 = mx_AFM2_glb_store(i);
					my_glb_local_AFM2 = my_AFM2_glb_store(i);
					mz_glb_local_AFM2 = mz_AFM2_glb_store(i);

					pt_math->transform_vector_glb2crt(mx_glb_local_AFM2, my_glb_local_AFM2, mz_glb_local_AFM2, \
						mx_AFM2, my_AFM2, mz_AFM2);

					Ms_AFM2 = mat->Ms_AFM2;
					K1_AFM2 = mat->K1_AFM2 / mu0 / Ms_AFM2; //K2 = mat->K2_AFM2 / mu0 / Ms;					
					udotm_AFM2 = mx_AFM2 * ux + my_AFM2 * uy + mz_AFM2 * uz;
					B1_AFM2 = mat->B1_AFM2 / mu0 / Ms_AFM2; B2_AFM2 = mat->B2_AFM2 / mu0 / Ms_AFM2;
					//Aexi_AFM2 = mat->Aex_AFM2; //iDMIi = mat->iDMI;
					JAFM_AFM2 = -1. * mat->J_AFM / mu0 / Ms_AFM2;

					get_m_AFM2_spatial_derivative_local \
						(x, y, z, \
							dmxdx_AFM2, dmydy_AFM2, dmzdx_AFM2, dmzdy_AFM2, \
							dmxdx2_f_AFM2, dmydx2_f_AFM2, dmzdx2_f_AFM2, \
							dmxdx2_b_AFM2, dmydx2_b_AFM2, dmzdx2_b_AFM2, \
							dmxdy2_f_AFM2, dmydy2_f_AFM2, dmzdy2_f_AFM2, \
							dmxdy2_b_AFM2, dmydy2_b_AFM2, dmzdy2_b_AFM2, \
							dmxdz2_f_AFM2, dmydz2_f_AFM2, dmzdz2_f_AFM2, \
							dmxdz2_b_AFM2, dmydz2_b_AFM2, dmzdz2_b_AFM2);

					//get_laplacian_m_local \
						//	(x, y, z, \
						//		d2mxdx2, d2mxdy2, d2mxdz2, \
						//		d2mydx2, d2mydy2, d2mydz2, \
						//		d2mzdx2, d2mzdy2, d2mzdz2);

						//------------Magnetocrystalline anisotropy field--------------//
					Hx_anis_crt_AFM2 = 2. * K1_AFM2 * udotm_AFM2 * ux;
					Hy_anis_crt_AFM2 = 2. * K1_AFM2 * udotm_AFM2 * uy;
					Hz_anis_crt_AFM2 = 2. * K1_AFM2 * udotm_AFM2 * uz;

					pt_math->transform_vector_crt2glb(Hx_anis_crt_AFM2, Hy_anis_crt_AFM2, Hz_anis_crt_AFM2, \
						Hx_anis_glb_AFM2, Hy_anis_glb_AFM2, Hz_anis_glb_AFM2);

					//Hx_anis(i) = Hx_anis_glb;
					//Hy_anis(i) = Hy_anis_glb;
					//Hz_anis(i) = Hz_anis_glb;

					//------------------Magnetoelastic effective field---------------//
					Hx_elas_crt_AFM2 = -1. * (B1_AFM2 * (exx)*mx_AFM2 + \
						B2_AFM2 * ((exy)*my_AFM2 + (exz)*mz_AFM2));
					Hy_elas_crt_AFM2 = -1. * (B1_AFM2 * (eyy)*my_AFM2 + \
						B2_AFM2 * ((exy)*mx_AFM2 + (eyz)*mz_AFM2));
					Hz_elas_crt_AFM2 = -1. * (B1_AFM2 * (ezz)*mz_AFM2 + \
						B2_AFM2 * ((exz)*mx_AFM2 + (eyz)*my_AFM2));

					pt_math->transform_vector_crt2glb(Hx_elas_crt_AFM2, Hy_elas_crt_AFM2, Hz_elas_crt_AFM2, \
						Hx_elas_glb_AFM2, Hy_elas_glb_AFM2, Hz_elas_glb_AFM2);

					//Hx_elas(i) = Hx_elas_glb;
					//Hy_elas(i) = Hy_elas_glb;
					//Hz_elas(i) = Hz_elas_glb;

					//-----------------Magnetic exchange field----------------------//
					//----------Averaged parameter for inter-region exchange - Approach 1---------//
					Hx_exch_glb_AFM2 = 2. / mu0 / Ms_AFM2 * \
						(Aij_AFM2_xf(i) * dmxdx2_f_AFM2 + Aij_AFM2_xb(i) * dmxdx2_b_AFM2 \
							+ Aij_AFM2_yf(i) * dmxdy2_f_AFM2 + Aij_AFM2_yb(i) * dmxdy2_b_AFM2 \
							+ Aij_AFM2_zf(i) * dmxdz2_f_AFM2 + Aij_AFM2_zb(i) * dmxdz2_b_AFM2);
					Hy_exch_glb_AFM2 = 2. / mu0 / Ms_AFM2 * \
						(Aij_AFM2_xf(i) * dmydx2_f_AFM2 + Aij_AFM2_xb(i) * dmydx2_b_AFM2 \
							+ Aij_AFM2_yf(i) * dmydy2_f_AFM2 + Aij_AFM2_yb(i) * dmydy2_b_AFM2 \
							+ Aij_AFM2_zf(i) * dmydz2_f_AFM2 + Aij_AFM2_zb(i) * dmydz2_b_AFM2);
					Hz_exch_glb_AFM2 = 2. / mu0 / Ms_AFM2 * \
						(Aij_AFM2_xf(i) * dmzdx2_f_AFM2 + Aij_AFM2_xb(i) * dmzdx2_b_AFM2 \
							+ Aij_AFM2_yf(i) * dmzdy2_f_AFM2 + Aij_AFM2_yb(i) * dmzdy2_b_AFM2 \
							+ Aij_AFM2_zf(i) * dmzdz2_f_AFM2 + Aij_AFM2_zb(i) * dmzdz2_b_AFM2);

					//Hx_exch(i) = 2. * Aex * (d2mxdx2 + d2mxdy2 + d2mxdz2);
					//Hy_exch(i) = 2. * Aex * (d2mydx2 + d2mydy2 + d2mydz2);
					//Hz_exch(i) = 2. * Aex * (d2mzdx2 + d2mzdy2 + d2mzdz2);
					// 
					// -----------------Inter-sublattice AFM exchange field---------------//
					Hx_iAFM_glb_AFM2 = JAFM_AFM2 * mx_AFM1_glb_store(i);
					Hy_iAFM_glb_AFM2 = JAFM_AFM2 * my_AFM1_glb_store(i);
					Hz_iAFM_glb_AFM2 = JAFM_AFM2 * mz_AFM1_glb_store(i);
					//------------------Interfacial DMI effective field-------------------//
					//Hx_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdx;
					//Hy_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdy;
					//Hz_iDMI_glb = -2. * iDMIi / mu0 / Ms * (dmxdx + dmydy);

					//-------------------LLG Equation------------------------------//
					Hx_eff_glb_AFM2 = Hx_stat(i) + \
						Hx_anis_glb_AFM2 + Hx_elas_glb_AFM2 + \
						Hx_exch_glb_AFM2 + Hx_iAFM_glb_AFM2 /*+ Hx_iDMI_glb*/ + pt_glb->Hext[0];

					Hy_eff_glb_AFM2 = Hy_stat(i) + \
						Hy_anis_glb_AFM2 + Hy_elas_glb_AFM2 + \
						Hy_exch_glb_AFM2 + Hy_iAFM_glb_AFM2 /*+ Hy_iDMI_glb*/ + pt_glb->Hext[1];

					Hz_eff_glb_AFM2 = Hz_stat(i) + \
						Hz_anis_glb_AFM2 + Hz_elas_glb_AFM2 + \
						Hz_exch_glb_AFM2 + Hz_iAFM_glb_AFM2 /*+ Hz_iDMI_glb*/ + pt_glb->Hext[2];

					if (pt_glb->if_EM_backaction_on_M == true) {
						Hx_eff_glb_AFM2 = Hx_eff_glb_AFM2 + pt_EM->DHx_em_cell(i);
						Hy_eff_glb_AFM2 = Hy_eff_glb_AFM2 + pt_EM->DHy_em_cell(i);
						Hz_eff_glb_AFM2 = Hz_eff_glb_AFM2 + pt_EM->DHz_em_cell(i);
					}

					alpha_AFM2 = mat->damping_AFM2;
					gyro_AFM2 = mat->gyro_AFM2;
					precess_AFM2 = -gyro_AFM2 / (1 + alpha_AFM2 * alpha_AFM2) * pt_glb->dt;

					precess_termx_AFM2 = my_glb_local_AFM2 * Hz_eff_glb_AFM2 - mz_glb_local_AFM2 * Hy_eff_glb_AFM2;
					precess_termy_AFM2 = mz_glb_local_AFM2 * Hx_eff_glb_AFM2 - mx_glb_local_AFM2 * Hz_eff_glb_AFM2;
					precess_termz_AFM2 = mx_glb_local_AFM2 * Hy_eff_glb_AFM2 - my_glb_local_AFM2 * Hx_eff_glb_AFM2;

					damping_termx_AFM2 = my_glb_local_AFM2 * precess_termz_AFM2 - mz_glb_local_AFM2 * precess_termy_AFM2;
					damping_termy_AFM2 = mz_glb_local_AFM2 * precess_termx_AFM2 - mx_glb_local_AFM2 * precess_termz_AFM2;
					damping_termz_AFM2 = mx_glb_local_AFM2 * precess_termy_AFM2 - my_glb_local_AFM2 * precess_termx_AFM2;

					dmx_AFM2_glb_rk3(i) = precess_AFM2 * precess_termx_AFM2 + alpha_AFM2 * precess_AFM2 * damping_termx_AFM2; //RK
					dmy_AFM2_glb_rk3(i) = precess_AFM2 * precess_termy_AFM2 + alpha_AFM2 * precess_AFM2 * damping_termy_AFM2; //RK
					dmz_AFM2_glb_rk3(i) = precess_AFM2 * precess_termz_AFM2 + alpha_AFM2 * precess_AFM2 * damping_termz_AFM2; //RK
				}
			}
		}
	}

	if (pt_glb->if_spin_pumping == true) {
		get_J_ISHE_RK3(); //RK
	}
}

void magnetic_system::get_dm_AFM_RK4() { //RK
	double Hx_anis_crt_AFM1, Hy_anis_crt_AFM1, Hz_anis_crt_AFM1;
	double Hx_anis_glb_AFM1, Hy_anis_glb_AFM1, Hz_anis_glb_AFM1;
	double Hx_elas_crt_AFM1, Hy_elas_crt_AFM1, Hz_elas_crt_AFM1;
	double Hx_elas_glb_AFM1, Hy_elas_glb_AFM1, Hz_elas_glb_AFM1;
	double Hx_exch_glb_AFM1, Hy_exch_glb_AFM1, Hz_exch_glb_AFM1;
	double Hx_iAFM_glb_AFM1, Hy_iAFM_glb_AFM1, Hz_iAFM_glb_AFM1;
	double Hx_eff_glb_AFM1, Hy_eff_glb_AFM1, Hz_eff_glb_AFM1;
	double precess_termx_AFM1, precess_termy_AFM1, precess_termz_AFM1;
	double damping_termx_AFM1, damping_termy_AFM1, damping_termz_AFM1;
	//double Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb;

	double Hx_anis_crt_AFM2, Hy_anis_crt_AFM2, Hz_anis_crt_AFM2;
	double Hx_anis_glb_AFM2, Hy_anis_glb_AFM2, Hz_anis_glb_AFM2;
	double Hx_elas_crt_AFM2, Hy_elas_crt_AFM2, Hz_elas_crt_AFM2;
	double Hx_elas_glb_AFM2, Hy_elas_glb_AFM2, Hz_elas_glb_AFM2;
	double Hx_exch_glb_AFM2, Hy_exch_glb_AFM2, Hz_exch_glb_AFM2;
	double Hx_iAFM_glb_AFM2, Hy_iAFM_glb_AFM2, Hz_iAFM_glb_AFM2;
	double Hx_eff_glb_AFM2, Hy_eff_glb_AFM2, Hz_eff_glb_AFM2;
	double precess_termx_AFM2, precess_termy_AFM2, precess_termz_AFM2;
	double damping_termx_AFM2, damping_termy_AFM2, damping_termz_AFM2;
	//double Hx_iDMI_glb_AFM2, Hy_iDMI_glb_AFM2, Hz_iDMI_glb_AFM2;

	long int x, y, z;
	unsigned int mat_type;
	material* mat;
	double Ms_AFM1, K1_AFM1, K2_AFM1, B1_AFM1, B2_AFM1;
	double Ms_AFM2, K1_AFM2, K2_AFM2, B1_AFM2, B2_AFM2;
	double ux, uy, uz;
	double udotm_AFM1, udotm_AFM2;
	//double Aexi;//, iDMIi;
	double JAFM_AFM1, JAFM_AFM2;
	double alpha_AFM1, gyro_AFM1, precess_AFM1;
	double alpha_AFM2, gyro_AFM2, precess_AFM2;
	double mx_glb_local_AFM1, my_glb_local_AFM1, mz_glb_local_AFM1;
	double mx_glb_local_AFM2, my_glb_local_AFM2, mz_glb_local_AFM2;
	double mx_AFM1, my_AFM1, mz_AFM1;
	double mx_AFM2, my_AFM2, mz_AFM2;
	double exx, eyy, ezz, eyz, exz, exy;
	//double d2mxdx2, d2mxdy2, d2mxdz2;
	//double d2mydx2, d2mydy2, d2mydz2;
	//double d2mzdx2, d2mzdy2, d2mzdz2;
	double dmxdx_AFM1, dmydy_AFM1, dmzdx_AFM1, dmzdy_AFM1;
	double dmxdx2_f_AFM1, dmydx2_f_AFM1, dmzdx2_f_AFM1;
	double dmxdx2_b_AFM1, dmydx2_b_AFM1, dmzdx2_b_AFM1;
	double dmxdy2_f_AFM1, dmydy2_f_AFM1, dmzdy2_f_AFM1;
	double dmxdy2_b_AFM1, dmydy2_b_AFM1, dmzdy2_b_AFM1;
	double dmxdz2_f_AFM1, dmydz2_f_AFM1, dmzdz2_f_AFM1;
	double dmxdz2_b_AFM1, dmydz2_b_AFM1, dmzdz2_b_AFM1;

	double dmxdx_AFM2, dmydy_AFM2, dmzdx_AFM2, dmzdy_AFM2;
	double dmxdx2_f_AFM2, dmydx2_f_AFM2, dmzdx2_f_AFM2;
	double dmxdx2_b_AFM2, dmydx2_b_AFM2, dmzdx2_b_AFM2;
	double dmxdy2_f_AFM2, dmydy2_f_AFM2, dmzdy2_f_AFM2;
	double dmxdy2_b_AFM2, dmydy2_b_AFM2, dmzdy2_b_AFM2;
	double dmxdz2_f_AFM2, dmydz2_f_AFM2, dmzdz2_f_AFM2;
	double dmxdz2_b_AFM2, dmydz2_b_AFM2, dmzdz2_b_AFM2;

#pragma acc serial default(present) async(4)
	{
		update_external_Hfield_full();  //RK
	}

#pragma acc parallel default(present) async(4)
	{
		//------------------Sublattice 1 ---------------------//
#pragma acc loop gang vector private(\
	Hx_anis_crt_AFM1, Hy_anis_crt_AFM1, Hz_anis_crt_AFM1,\
	Hx_anis_glb_AFM1, Hy_anis_glb_AFM1, Hz_anis_glb_AFM1,\
	Hx_elas_crt_AFM1, Hy_elas_crt_AFM1, Hz_elas_crt_AFM1,\
	Hx_elas_glb_AFM1, Hy_elas_glb_AFM1, Hz_elas_glb_AFM1,\
	Hx_exch_glb_AFM1, Hy_exch_glb_AFM1, Hz_exch_glb_AFM1,\
	Hx_iAFM_glb_AFM1, Hy_iAFM_glb_AFM1, Hz_iAFM_glb_AFM1,\
	Hx_eff_glb_AFM1, Hy_eff_glb_AFM1, Hz_eff_glb_AFM1,\
	precess_termx_AFM1, precess_termy_AFM1, precess_termz_AFM1,\
	damping_termx_AFM1, damping_termy_AFM1, damping_termz_AFM1,\
	x,y,z,mat_type,mat,\
	Ms_AFM1, K1_AFM1, K2_AFM1, B1_AFM1, B2_AFM1, \
	ux, uy, uz, udotm_AFM1, \
	JAFM_AFM1,\
	alpha_AFM1, gyro_AFM1, precess_AFM1,\
	mx_glb_local_AFM1, my_glb_local_AFM1, mz_glb_local_AFM1,\
	mx_AFM1, my_AFM1, mz_AFM1,\
	exx, eyy, ezz, eyz, exz, exy,\
	dmxdx_AFM1, dmydy_AFM1, dmzdx_AFM1, dmzdy_AFM1,\
	dmxdx2_f_AFM1, dmydx2_f_AFM1, dmzdx2_f_AFM1,\
	dmxdx2_b_AFM1, dmydx2_b_AFM1, dmzdx2_b_AFM1,\
	dmxdy2_f_AFM1, dmydy2_f_AFM1, dmzdy2_f_AFM1,\
	dmxdy2_b_AFM1, dmydy2_b_AFM1, dmzdy2_b_AFM1,\
	dmxdz2_f_AFM1, dmydz2_f_AFM1, dmzdz2_f_AFM1,\
	dmxdz2_b_AFM1, dmydz2_b_AFM1, dmzdz2_b_AFM1)
//Aexi, Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb, iDMIi)
		for (long int i = 0; i < n; i++) {
			mat_type = pt_glb->material_cell(i);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_AFM == true) {
					x = i / (ny * nz);
					y = (i - x * (ny * nz)) / nz;
					z = i - x * (ny * nz) - y * nz;

					ux = mat->uniaxis_AFM[0]; uy = mat->uniaxis_AFM[1]; uz = mat->uniaxis_AFM[2];

					if (pt_glb->if_elasto_on_pORm == true) {
						exx = pt_elas->exxt0_crt(i) + pt_elas->Dexx_crt(i);
						eyy = pt_elas->eyyt0_crt(i) + pt_elas->Deyy_crt(i);
						ezz = pt_elas->ezzt0_crt(i) + pt_elas->Dezz_crt(i);
						eyz = pt_elas->eyzt0_crt(i) + pt_elas->Deyz_crt(i);
						exz = pt_elas->exzt0_crt(i) + pt_elas->Dexz_crt(i);
						exy = pt_elas->exyt0_crt(i) + pt_elas->Dexy_crt(i);
					}
					else {
						exx = 0.;
						eyy = 0.;
						ezz = 0.;
						eyz = 0.;
						exz = 0.;
						exy = 0.;
					}

					mx_glb_local_AFM1 = mx_AFM1_glb_store(i);
					my_glb_local_AFM1 = my_AFM1_glb_store(i);
					mz_glb_local_AFM1 = mz_AFM1_glb_store(i);

					pt_math->transform_vector_glb2crt(mx_glb_local_AFM1, my_glb_local_AFM1, mz_glb_local_AFM1, \
						mx_AFM1, my_AFM1, mz_AFM1);

					Ms_AFM1 = mat->Ms_AFM1;
					K1_AFM1 = mat->K1_AFM1 / mu0 / Ms_AFM1; //K2_AFM1 = mat->K2 / mu0 / Ms;					
					udotm_AFM1 = mx_AFM1 * ux + my_AFM1 * uy + mz_AFM1 * uz;
					B1_AFM1 = mat->B1_AFM1 / mu0 / Ms_AFM1; B2_AFM1 = mat->B2_AFM1 / mu0 / Ms_AFM1;
					//Aexi_AFM1 = mat->Aex_AFM1; //iDMIi = mat->iDMI;
					JAFM_AFM1 = -1. * mat->J_AFM / mu0 / Ms_AFM1;

					get_m_AFM1_spatial_derivative_local \
						(x, y, z, \
							dmxdx_AFM1, dmydy_AFM1, dmzdx_AFM1, dmzdy_AFM1, \
							dmxdx2_f_AFM1, dmydx2_f_AFM1, dmzdx2_f_AFM1, \
							dmxdx2_b_AFM1, dmydx2_b_AFM1, dmzdx2_b_AFM1, \
							dmxdy2_f_AFM1, dmydy2_f_AFM1, dmzdy2_f_AFM1, \
							dmxdy2_b_AFM1, dmydy2_b_AFM1, dmzdy2_b_AFM1, \
							dmxdz2_f_AFM1, dmydz2_f_AFM1, dmzdz2_f_AFM1, \
							dmxdz2_b_AFM1, dmydz2_b_AFM1, dmzdz2_b_AFM1);

					//get_laplacian_m_local \
						//	(x, y, z, \
						//		d2mxdx2, d2mxdy2, d2mxdz2, \
						//		d2mydx2, d2mydy2, d2mydz2, \
						//		d2mzdx2, d2mzdy2, d2mzdz2);

						//------------Magnetocrystalline anisotropy field--------------//
					Hx_anis_crt_AFM1 = 2. * K1_AFM1 * udotm_AFM1 * ux;
					Hy_anis_crt_AFM1 = 2. * K1_AFM1 * udotm_AFM1 * uy;
					Hz_anis_crt_AFM1 = 2. * K1_AFM1 * udotm_AFM1 * uz;

					pt_math->transform_vector_crt2glb(Hx_anis_crt_AFM1, Hy_anis_crt_AFM1, Hz_anis_crt_AFM1, \
						Hx_anis_glb_AFM1, Hy_anis_glb_AFM1, Hz_anis_glb_AFM1);

					//Hx_anis(i) = Hx_anis_glb;
					//Hy_anis(i) = Hy_anis_glb;
					//Hz_anis(i) = Hz_anis_glb;

					//------------------Magnetoelastic effective field---------------//
					Hx_elas_crt_AFM1 = -1. * (B1_AFM1 * (exx)*mx_AFM1 + \
						B2_AFM1 * ((exy)*my_AFM1 + (exz)*mz_AFM1));
					Hy_elas_crt_AFM1 = -1. * (B1_AFM1 * (eyy)*my_AFM1 + \
						B2_AFM1 * ((exy)*mx_AFM1 + (eyz)*mz_AFM1));
					Hz_elas_crt_AFM1 = -1. * (B1_AFM1 * (ezz)*mz_AFM1 + \
						B2_AFM1 * ((exz)*mx_AFM1 + (eyz)*my_AFM1));

					pt_math->transform_vector_crt2glb(Hx_elas_crt_AFM1, Hy_elas_crt_AFM1, Hz_elas_crt_AFM1, \
						Hx_elas_glb_AFM1, Hy_elas_glb_AFM1, Hz_elas_glb_AFM1);

					//Hx_elas(i) = Hx_elas_glb;
					//Hy_elas(i) = Hy_elas_glb;
					//Hz_elas(i) = Hz_elas_glb;

					//-----------------Magnetic exchange field----------------------//
					//----------Averaged parameter for inter-region exchange - Approach 1---------//
					Hx_exch_glb_AFM1 = 2. / mu0 / Ms_AFM1 * \
						(Aij_AFM1_xf(i) * dmxdx2_f_AFM1 + Aij_AFM1_xb(i) * dmxdx2_b_AFM1 \
							+ Aij_AFM1_yf(i) * dmxdy2_f_AFM1 + Aij_AFM1_yb(i) * dmxdy2_b_AFM1 \
							+ Aij_AFM1_zf(i) * dmxdz2_f_AFM1 + Aij_AFM1_zb(i) * dmxdz2_b_AFM1);
					Hy_exch_glb_AFM1 = 2. / mu0 / Ms_AFM1 * \
						(Aij_AFM1_xf(i) * dmydx2_f_AFM1 + Aij_AFM1_xb(i) * dmydx2_b_AFM1 \
							+ Aij_AFM1_yf(i) * dmydy2_f_AFM1 + Aij_AFM1_yb(i) * dmydy2_b_AFM1 \
							+ Aij_AFM1_zf(i) * dmydz2_f_AFM1 + Aij_AFM1_zb(i) * dmydz2_b_AFM1);
					Hz_exch_glb_AFM1 = 2. / mu0 / Ms_AFM1 * \
						(Aij_AFM1_xf(i) * dmzdx2_f_AFM1 + Aij_AFM1_xb(i) * dmzdx2_b_AFM1 \
							+ Aij_AFM1_yf(i) * dmzdy2_f_AFM1 + Aij_AFM1_yb(i) * dmzdy2_b_AFM1 \
							+ Aij_AFM1_zf(i) * dmzdz2_f_AFM1 + Aij_AFM1_zb(i) * dmzdz2_b_AFM1);

					//Hx_exch(i) = 2. * Aex * (d2mxdx2 + d2mxdy2 + d2mxdz2);
					//Hy_exch(i) = 2. * Aex * (d2mydx2 + d2mydy2 + d2mydz2);
					//Hz_exch(i) = 2. * Aex * (d2mzdx2 + d2mzdy2 + d2mzdz2);
					// 
					// -----------------Inter-sublattice AFM exchange field---------------//
					Hx_iAFM_glb_AFM1 = JAFM_AFM1 * mx_AFM2_glb_store(i);
					Hy_iAFM_glb_AFM1 = JAFM_AFM1 * my_AFM2_glb_store(i);
					Hz_iAFM_glb_AFM1 = JAFM_AFM1 * mz_AFM2_glb_store(i);

					//------------------Interfacial DMI effective field-------------------//
					//Hx_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdx;
					//Hy_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdy;
					//Hz_iDMI_glb = -2. * iDMIi / mu0 / Ms * (dmxdx + dmydy);

					//-------------------LLG Equation------------------------------//
					Hx_eff_glb_AFM1 = Hx_stat(i) + \
						Hx_anis_glb_AFM1 + Hx_elas_glb_AFM1 + \
						Hx_exch_glb_AFM1 + Hx_iAFM_glb_AFM1 /*+ Hx_iDMI_glb*/ + pt_glb->Hext[0];

					Hy_eff_glb_AFM1 = Hy_stat(i) + \
						Hy_anis_glb_AFM1 + Hy_elas_glb_AFM1 + \
						Hy_exch_glb_AFM1 + Hy_iAFM_glb_AFM1/*+ Hy_iDMI_glb*/ + pt_glb->Hext[1];

					Hz_eff_glb_AFM1 = Hz_stat(i) + \
						Hz_anis_glb_AFM1 + Hz_elas_glb_AFM1 + \
						Hz_exch_glb_AFM1 + Hz_iAFM_glb_AFM1/*+ Hz_iDMI_glb*/ + pt_glb->Hext[2];

					if (pt_glb->if_EM_backaction_on_M == true) {
						Hx_eff_glb_AFM1 = Hx_eff_glb_AFM1 + pt_EM->DHx_em_cell(i);
						Hy_eff_glb_AFM1 = Hy_eff_glb_AFM1 + pt_EM->DHy_em_cell(i);
						Hz_eff_glb_AFM1 = Hz_eff_glb_AFM1 + pt_EM->DHz_em_cell(i);
					}

					alpha_AFM1 = mat->damping_AFM1;
					gyro_AFM1 = mat->gyro_AFM1;
					precess_AFM1 = -gyro_AFM1 / (1. + alpha_AFM1 * alpha_AFM1) * pt_glb->dt;

					precess_termx_AFM1 = my_glb_local_AFM1 * Hz_eff_glb_AFM1 - mz_glb_local_AFM1 * Hy_eff_glb_AFM1;
					precess_termy_AFM1 = mz_glb_local_AFM1 * Hx_eff_glb_AFM1 - mx_glb_local_AFM1 * Hz_eff_glb_AFM1;
					precess_termz_AFM1 = mx_glb_local_AFM1 * Hy_eff_glb_AFM1 - my_glb_local_AFM1 * Hx_eff_glb_AFM1;

					damping_termx_AFM1 = my_glb_local_AFM1 * precess_termz_AFM1 - mz_glb_local_AFM1 * precess_termy_AFM1;
					damping_termy_AFM1 = mz_glb_local_AFM1 * precess_termx_AFM1 - mx_glb_local_AFM1 * precess_termz_AFM1;
					damping_termz_AFM1 = mx_glb_local_AFM1 * precess_termy_AFM1 - my_glb_local_AFM1 * precess_termx_AFM1;

					dmx_AFM1_glb_rk4(i) = precess_AFM1 * precess_termx_AFM1 + alpha_AFM1 * precess_AFM1 * damping_termx_AFM1; //RK
					dmy_AFM1_glb_rk4(i) = precess_AFM1 * precess_termy_AFM1 + alpha_AFM1 * precess_AFM1 * damping_termy_AFM1; //RK
					dmz_AFM1_glb_rk4(i) = precess_AFM1 * precess_termz_AFM1 + alpha_AFM1 * precess_AFM1 * damping_termz_AFM1; //RK
				}
			}
		}

		//------------------Sublattice 2 ---------------------//
#pragma acc loop gang vector private(\
	Hx_anis_crt_AFM2, Hy_anis_crt_AFM2, Hz_anis_crt_AFM2,\
	Hx_anis_glb_AFM2, Hy_anis_glb_AFM2, Hz_anis_glb_AFM2,\
	Hx_elas_crt_AFM2, Hy_elas_crt_AFM2, Hz_elas_crt_AFM2,\
	Hx_elas_glb_AFM2, Hy_elas_glb_AFM2, Hz_elas_glb_AFM2,\
	Hx_exch_glb_AFM2, Hy_exch_glb_AFM2, Hz_exch_glb_AFM2,\
	Hx_iAFM_glb_AFM2, Hy_iAFM_glb_AFM2, Hz_iAFM_glb_AFM2,\
	Hx_eff_glb_AFM2, Hy_eff_glb_AFM2, Hz_eff_glb_AFM2,\
	precess_termx_AFM2, precess_termy_AFM2, precess_termz_AFM2,\
	damping_termx_AFM2, damping_termy_AFM2, damping_termz_AFM2,\
	x,y,z,mat_type,mat,\
	Ms_AFM2, K1_AFM2, K2_AFM2, B1_AFM2, B2_AFM2, \
	ux, uy, uz, udotm_AFM2, \
	JAFM_AFM2,\
	alpha_AFM2, gyro_AFM2, precess_AFM2,\
	mx_glb_local_AFM2, my_glb_local_AFM2, mz_glb_local_AFM2,\
	mx_AFM2, my_AFM2, mz_AFM2,\
	exx, eyy, ezz, eyz, exz, exy,\
	dmxdx_AFM2, dmydy_AFM2, dmzdx_AFM2, dmzdy_AFM2,\
	dmxdx2_f_AFM2, dmydx2_f_AFM2, dmzdx2_f_AFM2,\
	dmxdx2_b_AFM2, dmydx2_b_AFM2, dmzdx2_b_AFM2,\
	dmxdy2_f_AFM2, dmydy2_f_AFM2, dmzdy2_f_AFM2,\
	dmxdy2_b_AFM2, dmydy2_b_AFM2, dmzdy2_b_AFM2,\
	dmxdz2_f_AFM2, dmydz2_f_AFM2, dmzdz2_f_AFM2,\
	dmxdz2_b_AFM2, dmydz2_b_AFM2, dmzdz2_b_AFM2)
	//Aexi, Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb, iDMIi)
		for (long int i = 0; i < n; i++) {
			mat_type = pt_glb->material_cell(i);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_AFM == true) {
					x = i / (ny * nz);
					y = (i - x * (ny * nz)) / nz;
					z = i - x * (ny * nz) - y * nz;

					ux = mat->uniaxis_AFM[0]; uy = mat->uniaxis_AFM[1]; uz = mat->uniaxis_AFM[2];

					if (pt_glb->if_elasto_on_pORm == true) {
						exx = pt_elas->exxt0_crt(i) + pt_elas->Dexx_crt(i);
						eyy = pt_elas->eyyt0_crt(i) + pt_elas->Deyy_crt(i);
						ezz = pt_elas->ezzt0_crt(i) + pt_elas->Dezz_crt(i);
						eyz = pt_elas->eyzt0_crt(i) + pt_elas->Deyz_crt(i);
						exz = pt_elas->exzt0_crt(i) + pt_elas->Dexz_crt(i);
						exy = pt_elas->exyt0_crt(i) + pt_elas->Dexy_crt(i);
					}
					else {
						exx = 0.;
						eyy = 0.;
						ezz = 0.;
						eyz = 0.;
						exz = 0.;
						exy = 0.;
					}

					mx_glb_local_AFM2 = mx_AFM2_glb_store(i);
					my_glb_local_AFM2 = my_AFM2_glb_store(i);
					mz_glb_local_AFM2 = mz_AFM2_glb_store(i);

					pt_math->transform_vector_glb2crt(mx_glb_local_AFM2, my_glb_local_AFM2, mz_glb_local_AFM2, \
						mx_AFM2, my_AFM2, mz_AFM2);

					Ms_AFM2 = mat->Ms_AFM2;
					K1_AFM2 = mat->K1_AFM2 / mu0 / Ms_AFM2; //K2 = mat->K2_AFM2 / mu0 / Ms;					
					udotm_AFM2 = mx_AFM2 * ux + my_AFM2 * uy + mz_AFM2 * uz;
					B1_AFM2 = mat->B1_AFM2 / mu0 / Ms_AFM2; B2_AFM2 = mat->B2_AFM2 / mu0 / Ms_AFM2;
					//Aexi_AFM2 = mat->Aex_AFM2; //iDMIi = mat->iDMI;
					JAFM_AFM2 = -1. * mat->J_AFM / mu0 / Ms_AFM2;

					get_m_AFM2_spatial_derivative_local \
						(x, y, z, \
							dmxdx_AFM2, dmydy_AFM2, dmzdx_AFM2, dmzdy_AFM2, \
							dmxdx2_f_AFM2, dmydx2_f_AFM2, dmzdx2_f_AFM2, \
							dmxdx2_b_AFM2, dmydx2_b_AFM2, dmzdx2_b_AFM2, \
							dmxdy2_f_AFM2, dmydy2_f_AFM2, dmzdy2_f_AFM2, \
							dmxdy2_b_AFM2, dmydy2_b_AFM2, dmzdy2_b_AFM2, \
							dmxdz2_f_AFM2, dmydz2_f_AFM2, dmzdz2_f_AFM2, \
							dmxdz2_b_AFM2, dmydz2_b_AFM2, dmzdz2_b_AFM2);

					//get_laplacian_m_local \
						//	(x, y, z, \
						//		d2mxdx2, d2mxdy2, d2mxdz2, \
						//		d2mydx2, d2mydy2, d2mydz2, \
						//		d2mzdx2, d2mzdy2, d2mzdz2);

						//------------Magnetocrystalline anisotropy field--------------//
					Hx_anis_crt_AFM2 = 2. * K1_AFM2 * udotm_AFM2 * ux;
					Hy_anis_crt_AFM2 = 2. * K1_AFM2 * udotm_AFM2 * uy;
					Hz_anis_crt_AFM2 = 2. * K1_AFM2 * udotm_AFM2 * uz;

					pt_math->transform_vector_crt2glb(Hx_anis_crt_AFM2, Hy_anis_crt_AFM2, Hz_anis_crt_AFM2, \
						Hx_anis_glb_AFM2, Hy_anis_glb_AFM2, Hz_anis_glb_AFM2);

					//Hx_anis(i) = Hx_anis_glb;
					//Hy_anis(i) = Hy_anis_glb;
					//Hz_anis(i) = Hz_anis_glb;

					//------------------Magnetoelastic effective field---------------//
					Hx_elas_crt_AFM2 = -1. * (B1_AFM2 * (exx)*mx_AFM2 + \
						B2_AFM2 * ((exy)*my_AFM2 + (exz)*mz_AFM2));
					Hy_elas_crt_AFM2 = -1. * (B1_AFM2 * (eyy)*my_AFM2 + \
						B2_AFM2 * ((exy)*mx_AFM2 + (eyz)*mz_AFM2));
					Hz_elas_crt_AFM2 = -1. * (B1_AFM2 * (ezz)*mz_AFM2 + \
						B2_AFM2 * ((exz)*mx_AFM2 + (eyz)*my_AFM2));

					pt_math->transform_vector_crt2glb(Hx_elas_crt_AFM2, Hy_elas_crt_AFM2, Hz_elas_crt_AFM2, \
						Hx_elas_glb_AFM2, Hy_elas_glb_AFM2, Hz_elas_glb_AFM2);

					//Hx_elas(i) = Hx_elas_glb;
					//Hy_elas(i) = Hy_elas_glb;
					//Hz_elas(i) = Hz_elas_glb;

					//-----------------Magnetic exchange field----------------------//
					//----------Averaged parameter for inter-region exchange - Approach 1---------//
					Hx_exch_glb_AFM2 = 2. / mu0 / Ms_AFM2 * \
						(Aij_AFM2_xf(i) * dmxdx2_f_AFM2 + Aij_AFM2_xb(i) * dmxdx2_b_AFM2 \
							+ Aij_AFM2_yf(i) * dmxdy2_f_AFM2 + Aij_AFM2_yb(i) * dmxdy2_b_AFM2 \
							+ Aij_AFM2_zf(i) * dmxdz2_f_AFM2 + Aij_AFM2_zb(i) * dmxdz2_b_AFM2);
					Hy_exch_glb_AFM2 = 2. / mu0 / Ms_AFM2 * \
						(Aij_AFM2_xf(i) * dmydx2_f_AFM2 + Aij_AFM2_xb(i) * dmydx2_b_AFM2 \
							+ Aij_AFM2_yf(i) * dmydy2_f_AFM2 + Aij_AFM2_yb(i) * dmydy2_b_AFM2 \
							+ Aij_AFM2_zf(i) * dmydz2_f_AFM2 + Aij_AFM2_zb(i) * dmydz2_b_AFM2);
					Hz_exch_glb_AFM2 = 2. / mu0 / Ms_AFM2 * \
						(Aij_AFM2_xf(i) * dmzdx2_f_AFM2 + Aij_AFM2_xb(i) * dmzdx2_b_AFM2 \
							+ Aij_AFM2_yf(i) * dmzdy2_f_AFM2 + Aij_AFM2_yb(i) * dmzdy2_b_AFM2 \
							+ Aij_AFM2_zf(i) * dmzdz2_f_AFM2 + Aij_AFM2_zb(i) * dmzdz2_b_AFM2);

					//Hx_exch(i) = 2. * Aex * (d2mxdx2 + d2mxdy2 + d2mxdz2);
					//Hy_exch(i) = 2. * Aex * (d2mydx2 + d2mydy2 + d2mydz2);
					//Hz_exch(i) = 2. * Aex * (d2mzdx2 + d2mzdy2 + d2mzdz2);
					// 
					// -----------------Inter-sublattice AFM exchange field---------------//
					Hx_iAFM_glb_AFM2 = JAFM_AFM2 * mx_AFM1_glb_store(i);
					Hy_iAFM_glb_AFM2 = JAFM_AFM2 * my_AFM1_glb_store(i);
					Hz_iAFM_glb_AFM2 = JAFM_AFM2 * mz_AFM1_glb_store(i);
					//------------------Interfacial DMI effective field-------------------//
					//Hx_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdx;
					//Hy_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdy;
					//Hz_iDMI_glb = -2. * iDMIi / mu0 / Ms * (dmxdx + dmydy);

					//-------------------LLG Equation------------------------------//
					Hx_eff_glb_AFM2 = Hx_stat(i) + \
						Hx_anis_glb_AFM2 + Hx_elas_glb_AFM2 + \
						Hx_exch_glb_AFM2 + Hx_iAFM_glb_AFM2 /*+ Hx_iDMI_glb*/ + pt_glb->Hext[0];

					Hy_eff_glb_AFM2 = Hy_stat(i) + \
						Hy_anis_glb_AFM2 + Hy_elas_glb_AFM2 + \
						Hy_exch_glb_AFM2 + Hy_iAFM_glb_AFM2 /*+ Hy_iDMI_glb*/ + pt_glb->Hext[1];

					Hz_eff_glb_AFM2 = Hz_stat(i) + \
						Hz_anis_glb_AFM2 + Hz_elas_glb_AFM2 + \
						Hz_exch_glb_AFM2 + Hz_iAFM_glb_AFM2 /*+ Hz_iDMI_glb*/ + pt_glb->Hext[2];

					if (pt_glb->if_EM_backaction_on_M == true) {
						Hx_eff_glb_AFM2 = Hx_eff_glb_AFM2 + pt_EM->DHx_em_cell(i);
						Hy_eff_glb_AFM2 = Hy_eff_glb_AFM2 + pt_EM->DHy_em_cell(i);
						Hz_eff_glb_AFM2 = Hz_eff_glb_AFM2 + pt_EM->DHz_em_cell(i);
					}

					alpha_AFM2 = mat->damping_AFM2;
					gyro_AFM2 = mat->gyro_AFM2;
					precess_AFM2 = -gyro_AFM2 / (1 + alpha_AFM2 * alpha_AFM2) * pt_glb->dt;

					precess_termx_AFM2 = my_glb_local_AFM2 * Hz_eff_glb_AFM2 - mz_glb_local_AFM2 * Hy_eff_glb_AFM2;
					precess_termy_AFM2 = mz_glb_local_AFM2 * Hx_eff_glb_AFM2 - mx_glb_local_AFM2 * Hz_eff_glb_AFM2;
					precess_termz_AFM2 = mx_glb_local_AFM2 * Hy_eff_glb_AFM2 - my_glb_local_AFM2 * Hx_eff_glb_AFM2;

					damping_termx_AFM2 = my_glb_local_AFM2 * precess_termz_AFM2 - mz_glb_local_AFM2 * precess_termy_AFM2;
					damping_termy_AFM2 = mz_glb_local_AFM2 * precess_termx_AFM2 - mx_glb_local_AFM2 * precess_termz_AFM2;
					damping_termz_AFM2 = mx_glb_local_AFM2 * precess_termy_AFM2 - my_glb_local_AFM2 * precess_termx_AFM2;

					dmx_AFM2_glb_rk4(i) = precess_AFM2 * precess_termx_AFM2 + alpha_AFM2 * precess_AFM2 * damping_termx_AFM2; //RK
					dmy_AFM2_glb_rk4(i) = precess_AFM2 * precess_termy_AFM2 + alpha_AFM2 * precess_AFM2 * damping_termy_AFM2; //RK
					dmz_AFM2_glb_rk4(i) = precess_AFM2 * precess_termz_AFM2 + alpha_AFM2 * precess_AFM2 * damping_termz_AFM2; //RK
				}
			}
		}
	}

	if (pt_glb->if_spin_pumping == true) {
		get_J_ISHE_RK4(); //RK
	}
}

void magnetic_system::get_dm_AFM() {
#pragma acc parallel default(present) async(4)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			dmx_AFM1_glb(id) = dmx_AFM1_glb_rk1(id) / 6. + dmx_AFM1_glb_rk2(id) / 3. + dmx_AFM1_glb_rk3(id) / 3. + dmx_AFM1_glb_rk4(id) / 6.;
			dmy_AFM1_glb(id) = dmy_AFM1_glb_rk1(id) / 6. + dmy_AFM1_glb_rk2(id) / 3. + dmy_AFM1_glb_rk3(id) / 3. + dmy_AFM1_glb_rk4(id) / 6.;
			dmz_AFM1_glb(id) = dmz_AFM1_glb_rk1(id) / 6. + dmz_AFM1_glb_rk2(id) / 3. + dmz_AFM1_glb_rk3(id) / 3. + dmz_AFM1_glb_rk4(id) / 6.;

			dmx_AFM2_glb(id) = dmx_AFM2_glb_rk1(id) / 6. + dmx_AFM2_glb_rk2(id) / 3. + dmx_AFM2_glb_rk3(id) / 3. + dmx_AFM2_glb_rk4(id) / 6.;
			dmy_AFM2_glb(id) = dmy_AFM2_glb_rk1(id) / 6. + dmy_AFM2_glb_rk2(id) / 3. + dmy_AFM2_glb_rk3(id) / 3. + dmy_AFM2_glb_rk4(id) / 6.;
			dmz_AFM2_glb(id) = dmz_AFM2_glb_rk1(id) / 6. + dmz_AFM2_glb_rk2(id) / 3. + dmz_AFM2_glb_rk3(id) / 3. + dmz_AFM2_glb_rk4(id) / 6.;
		}
	}
}

void magnetic_system::update_m_AFM_RK1() {
#pragma acc parallel default(present) async(3)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			mx_AFM1_glb_store(id) = mx_AFM1_glb(id) + dmx_AFM1_glb_rk1(id) * 0.5;
			my_AFM1_glb_store(id) = my_AFM1_glb(id) + dmy_AFM1_glb_rk1(id) * 0.5;
			mz_AFM1_glb_store(id) = mz_AFM1_glb(id) + dmz_AFM1_glb_rk1(id) * 0.5;

			mx_AFM2_glb_store(id) = mx_AFM2_glb(id) + dmx_AFM2_glb_rk1(id) * 0.5;
			my_AFM2_glb_store(id) = my_AFM2_glb(id) + dmy_AFM2_glb_rk1(id) * 0.5;
			mz_AFM2_glb_store(id) = mz_AFM2_glb(id) + dmz_AFM2_glb_rk1(id) * 0.5;
		}
	}

	if (pt_glb->if_EMdynamic == false && pt_glb->if_mag_1Dmodel == true) {
#pragma acc parallel default(present) async(3)
		{
			get_H_static_1D();
		}
	}
}

void magnetic_system::update_m_AFM_RK2() {
#pragma acc parallel default(present) async(3)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			mx_AFM1_glb_store(id) = mx_AFM1_glb(id) + dmx_AFM1_glb_rk2(id) * 0.5;
			my_AFM1_glb_store(id) = my_AFM1_glb(id) + dmy_AFM1_glb_rk2(id) * 0.5;
			mz_AFM1_glb_store(id) = mz_AFM1_glb(id) + dmz_AFM1_glb_rk2(id) * 0.5;

			mx_AFM2_glb_store(id) = mx_AFM2_glb(id) + dmx_AFM2_glb_rk2(id) * 0.5;
			my_AFM2_glb_store(id) = my_AFM2_glb(id) + dmy_AFM2_glb_rk2(id) * 0.5;
			mz_AFM2_glb_store(id) = mz_AFM2_glb(id) + dmz_AFM2_glb_rk2(id) * 0.5;
		}
	}

	if (pt_glb->if_EMdynamic == false && pt_glb->if_mag_1Dmodel == true) {
#pragma acc parallel default(present) async(3)
		{
			get_H_static_1D();
		}
	}
}

void magnetic_system::update_m_AFM_RK3() {
#pragma acc parallel default(present) async(3)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			mx_AFM1_glb_store(id) = mx_AFM1_glb(id) + dmx_AFM1_glb_rk3(id);
			my_AFM1_glb_store(id) = my_AFM1_glb(id) + dmy_AFM1_glb_rk3(id);
			mz_AFM1_glb_store(id) = mz_AFM1_glb(id) + dmz_AFM1_glb_rk3(id);

			mx_AFM2_glb_store(id) = mx_AFM2_glb(id) + dmx_AFM2_glb_rk3(id);
			my_AFM2_glb_store(id) = my_AFM2_glb(id) + dmy_AFM2_glb_rk3(id);
			mz_AFM2_glb_store(id) = mz_AFM2_glb(id) + dmz_AFM2_glb_rk3(id);
		}
	}

	if (pt_glb->if_EMdynamic == false && pt_glb->if_mag_1Dmodel == true) {
#pragma acc parallel default(present) async(3)
		{
			get_H_static_1D();
		}
	}
}

void magnetic_system::update_m_AFM() {

	if (pt_glb->if_spin_pumping == true) {
		get_J_ISHE();
	}

#pragma acc parallel default(present) async(3)
	{
		mx_AFM1_glb += dmx_AFM1_glb;
		my_AFM1_glb += dmy_AFM1_glb;
		mz_AFM1_glb += dmz_AFM1_glb;

		mx_AFM2_glb += dmx_AFM2_glb;
		my_AFM2_glb += dmy_AFM2_glb;
		mz_AFM2_glb += dmz_AFM2_glb;
	}

#pragma acc parallel default(present) async(3)
	{
		mx_AFM1_glb_store = mx_AFM1_glb;
		my_AFM1_glb_store = my_AFM1_glb;
		mz_AFM1_glb_store = mz_AFM1_glb;

		mx_AFM2_glb_store = mx_AFM2_glb;
		my_AFM2_glb_store = my_AFM2_glb;
		mz_AFM2_glb_store = mz_AFM2_glb;
	}

	if (pt_glb->if_EMdynamic == false && pt_glb->if_mag_1Dmodel == true) {
#pragma acc parallel default(present) async(3)
		{
			get_H_static_1D();
		}
	}
}

#pragma acc routine seq// nohost
void magnetic_system::get_m_AFM1_spatial_derivative_local(long int& i, long int& j, long int& k, \
	double& dmxdx, double& dmydy, double& dmzdx, double& dmzdy, \
	double& dmxdx2_f, double& dmydx2_f, double& dmzdx2_f, \
	double& dmxdx2_b, double& dmydx2_b, double& dmzdx2_b, \
	double& dmxdy2_f, double& dmydy2_f, double& dmzdy2_f, \
	double& dmxdy2_b, double& dmydy2_b, double& dmzdy2_b, \
	double& dmxdz2_f, double& dmydz2_f, double& dmzdz2_f, \
	double& dmxdz2_b, double& dmydz2_b, double& dmzdz2_b) {
	double mx_fwd, mx_center, mx_bwd;
	double my_fwd, my_center, my_bwd;
	double mz_fwd, mz_center, mz_bwd;

	double dmxdx_fwd, dmydx_fwd, dmzdx_fwd;
	double dmxdx_bwd, dmydx_bwd, dmzdx_bwd;
	double dmxdy_fwd, dmydy_fwd, dmzdy_fwd;
	double dmxdy_bwd, dmydy_bwd, dmzdy_bwd;
	double dmxdz_fwd, dmydz_fwd, dmzdz_fwd;
	double dmxdz_bwd, dmydz_bwd, dmzdz_bwd;

	mx_center = mx_AFM1_glb_store(i, j, k);
	my_center = my_AFM1_glb_store(i, j, k);
	mz_center = mz_AFM1_glb_store(i, j, k);

	//--------------X direction (YZ surface)------------//
	//------------- negative X direction---------------//
	if (AFM_surfYZ(i, j, k) == true) {
		//mx_bwd = mx_center; my_bwd = my_center; mz_bwd = mz_center;
		dmxdx_bwd = 0.; // -D / 2. / A * mz_center;
		dmydx_bwd = 0.;
		dmzdx_bwd = 0.; // D / 2. / A * mx_center;
	}
	else if (AFM_surfYZ(i, j, k) == false) {
		if (i == 0) {
			mx_bwd = mx_AFM1_glb_store(nx - 1, j, k);
			my_bwd = my_AFM1_glb_store(nx - 1, j, k);
			mz_bwd = mz_AFM1_glb_store(nx - 1, j, k);
		}
		else {
			mx_bwd = mx_AFM1_glb_store(i - 1, j, k);
			my_bwd = my_AFM1_glb_store(i - 1, j, k);
			mz_bwd = mz_AFM1_glb_store(i - 1, j, k);
		}

		dmxdx_bwd = (mx_center - mx_bwd) / dx;
		dmydx_bwd = (my_center - my_bwd) / dx;
		dmzdx_bwd = (mz_center - mz_bwd) / dx;
	}
	//------------- positive X direction---------------//
	if (AFM_surfYZ(i + 1, j, k) == true) {
		//mx_fwd = mx_center; my_fwd = my_center; mz_fwd = mz_center;
		dmxdx_fwd = 0.;	// -D / 2. / A * mz_center;
		dmydx_fwd = 0.;
		dmzdx_fwd = 0.; // D / 2. / A * mx_center;
	}
	else if (AFM_surfYZ(i + 1, j, k) == false) {
		if (i == nx - 1) {
			mx_fwd = mx_AFM1_glb_store(0, j, k);
			my_fwd = my_AFM1_glb_store(0, j, k);
			mz_fwd = mz_AFM1_glb_store(0, j, k);
		}
		else {
			mx_fwd = mx_AFM1_glb_store(i + 1, j, k);
			my_fwd = my_AFM1_glb_store(i + 1, j, k);
			mz_fwd = mz_AFM1_glb_store(i + 1, j, k);
		}
		dmxdx_fwd = (mx_fwd - mx_center) / dx;
		dmydx_fwd = (my_fwd - my_center) / dx;
		dmzdx_fwd = (mz_fwd - mz_center) / dx;
	}

	dmxdx = (dmxdx_bwd + dmxdx_fwd) * 0.5;
	dmzdx = (dmzdx_bwd + dmzdx_fwd) * 0.5;
	dmxdx2_f = dmxdx_fwd / dx;
	dmydx2_f = dmydx_fwd / dx;
	dmzdx2_f = dmzdx_fwd / dx;
	dmxdx2_b = -dmxdx_bwd / dx;
	dmydx2_b = -dmydx_bwd / dx;
	dmzdx2_b = -dmzdx_bwd / dx;
	//--------------END X direction (YZ surface)------------//

	//--------------Y direction (XZ surface)------------//
	//------------- negative Y direction---------------//
	if (AFM_surfXZ(i, j, k) == true) {
		//mx_bwd = mx_center; my_bwd = my_center; mz_bwd = mz_center;
		dmxdy_bwd = 0.;
		dmydy_bwd = 0.; // -D / 2. / A * mz_center;
		dmzdy_bwd = 0.; // D / 2. / A * my_center;
	}
	else if (AFM_surfXZ(i, j, k) == false) {
		if (j == 0) {
			mx_bwd = mx_AFM1_glb_store(i, ny - 1, k);
			my_bwd = my_AFM1_glb_store(i, ny - 1, k);
			mz_bwd = mz_AFM1_glb_store(i, ny - 1, k);
		}
		else {
			mx_bwd = mx_AFM1_glb_store(i, j - 1, k);
			my_bwd = my_AFM1_glb_store(i, j - 1, k);
			mz_bwd = mz_AFM1_glb_store(i, j - 1, k);
		}

		dmxdy_bwd = (mx_center - mx_bwd) / dy;
		dmydy_bwd = (my_center - my_bwd) / dy;
		dmzdy_bwd = (mz_center - mz_bwd) / dy;
	}
	//------------- positive Y direction---------------//
	if (AFM_surfXZ(i, j + 1, k) == true) {
		//mx_fwd = mx_center; my_fwd = my_center; mz_fwd = mz_center;
		dmxdy_fwd = 0.;
		dmydy_fwd = 0.; // -D / 2. / A * mz_center;
		dmzdy_fwd = 0.; // D / 2. / A * my_center;
	}
	else if (AFM_surfXZ(i, j + 1, k) == false) {
		if (j == ny - 1) {
			mx_fwd = mx_AFM1_glb_store(i, 0, k);
			my_fwd = my_AFM1_glb_store(i, 0, k);
			mz_fwd = mz_AFM1_glb_store(i, 0, k);
		}
		else {
			mx_fwd = mx_AFM1_glb_store(i, j + 1, k);
			my_fwd = my_AFM1_glb_store(i, j + 1, k);
			mz_fwd = mz_AFM1_glb_store(i, j + 1, k);
		}

		dmxdy_fwd = (mx_fwd - mx_center) / dy;
		dmydy_fwd = (my_fwd - my_center) / dy;
		dmzdy_fwd = (mz_fwd - mz_center) / dy;
	}

	dmydy = (dmydy_bwd + dmydy_fwd) * 0.5;
	dmzdy = (dmzdy_bwd + dmzdy_fwd) * 0.5;
	dmxdy2_f = dmxdy_fwd / dy;
	dmydy2_f = dmydy_fwd / dy;
	dmzdy2_f = dmzdy_fwd / dy;
	dmxdy2_b = -dmxdy_bwd / dy;
	dmydy2_b = -dmydy_bwd / dy;
	dmzdy2_b = -dmzdy_bwd / dy;
	//--------------END Y direction (XZ surface)------------//	

	//--------------Z direction (XY surface)------------//
	//------------- negative Z direction---------------//
	if (AFM_surfXY(i, j, k) == true) {
		//mx_bwd = mx_center; my_bwd = my_center; mz_bwd = mz_center;
		dmxdz_bwd = 0.; dmydz_bwd = 0.; dmzdz_bwd = 0.;
	}
	else if (AFM_surfXY(i, j, k) == false) {
		if (k == 0) {
			mx_bwd = mx_AFM1_glb_store(i, j, nz - 1);
			my_bwd = my_AFM1_glb_store(i, j, nz - 1);
			mz_bwd = mz_AFM1_glb_store(i, j, nz - 1);
		}
		else {
			mx_bwd = mx_AFM1_glb_store(i, j, k - 1);
			my_bwd = my_AFM1_glb_store(i, j, k - 1);
			mz_bwd = mz_AFM1_glb_store(i, j, k - 1);
		}

		dmxdz_bwd = (mx_center - mx_bwd) / dz;
		dmydz_bwd = (my_center - my_bwd) / dz;
		dmzdz_bwd = (mz_center - mz_bwd) / dz;
	}
	//------------- positive Z direction---------------//
	if (AFM_surfXY(i, j, k + 1) == true) {
		//mx_fwd = mx_center; my_fwd = my_center; mz_fwd = mz_center;
		dmxdz_fwd = 0.; dmydz_fwd = 0.; dmzdz_fwd = 0.;
	}
	else if (AFM_surfXY(i, j, k + 1) == false) {
		if (k == nz - 1) {
			mx_fwd = mx_AFM1_glb_store(i, j, 0);
			my_fwd = my_AFM1_glb_store(i, j, 0);
			mz_fwd = mz_AFM1_glb_store(i, j, 0);
		}
		else {
			mx_fwd = mx_AFM1_glb_store(i, j, k + 1);
			my_fwd = my_AFM1_glb_store(i, j, k + 1);
			mz_fwd = mz_AFM1_glb_store(i, j, k + 1);
		}

		dmxdz_fwd = (mx_fwd - mx_center) / dz;
		dmydz_fwd = (my_fwd - my_center) / dz;
		dmzdz_fwd = (mz_fwd - mz_center) / dz;
	}

	dmxdz2_f = dmxdz_fwd / dz;
	dmydz2_f = dmydz_fwd / dz;
	dmzdz2_f = dmzdz_fwd / dz;
	dmxdz2_b = -dmxdz_bwd / dz;
	dmydz2_b = -dmydz_bwd / dz;
	dmzdz2_b = -dmzdz_bwd / dz;
	//--------------END Z direction (XY surface)------------//
}

#pragma acc routine seq// nohost
void magnetic_system::get_m_AFM2_spatial_derivative_local(long int& i, long int& j, long int& k, \
	double& dmxdx, double& dmydy, double& dmzdx, double& dmzdy, \
	double& dmxdx2_f, double& dmydx2_f, double& dmzdx2_f, \
	double& dmxdx2_b, double& dmydx2_b, double& dmzdx2_b, \
	double& dmxdy2_f, double& dmydy2_f, double& dmzdy2_f, \
	double& dmxdy2_b, double& dmydy2_b, double& dmzdy2_b, \
	double& dmxdz2_f, double& dmydz2_f, double& dmzdz2_f, \
	double& dmxdz2_b, double& dmydz2_b, double& dmzdz2_b) {
	double mx_fwd, mx_center, mx_bwd;
	double my_fwd, my_center, my_bwd;
	double mz_fwd, mz_center, mz_bwd;

	double dmxdx_fwd, dmydx_fwd, dmzdx_fwd;
	double dmxdx_bwd, dmydx_bwd, dmzdx_bwd;
	double dmxdy_fwd, dmydy_fwd, dmzdy_fwd;
	double dmxdy_bwd, dmydy_bwd, dmzdy_bwd;
	double dmxdz_fwd, dmydz_fwd, dmzdz_fwd;
	double dmxdz_bwd, dmydz_bwd, dmzdz_bwd;

	mx_center = mx_AFM2_glb_store(i, j, k);
	my_center = my_AFM2_glb_store(i, j, k);
	mz_center = mz_AFM2_glb_store(i, j, k);

	//--------------X direction (YZ surface)------------//
	//------------- negative X direction---------------//
	if (AFM_surfYZ(i, j, k) == true) {
		//mx_bwd = mx_center; my_bwd = my_center; mz_bwd = mz_center;
		dmxdx_bwd = 0.; // -D / 2. / A * mz_center;
		dmydx_bwd = 0.;
		dmzdx_bwd = 0.; // D / 2. / A * mx_center;
	}
	else if (AFM_surfYZ(i, j, k) == false) {
		if (i == 0) {
			mx_bwd = mx_AFM2_glb_store(nx - 1, j, k);
			my_bwd = my_AFM2_glb_store(nx - 1, j, k);
			mz_bwd = mz_AFM2_glb_store(nx - 1, j, k);
		}
		else {
			mx_bwd = mx_AFM2_glb_store(i - 1, j, k);
			my_bwd = my_AFM2_glb_store(i - 1, j, k);
			mz_bwd = mz_AFM2_glb_store(i - 1, j, k);
		}

		dmxdx_bwd = (mx_center - mx_bwd) / dx;
		dmydx_bwd = (my_center - my_bwd) / dx;
		dmzdx_bwd = (mz_center - mz_bwd) / dx;
	}
	//------------- positive X direction---------------//
	if (AFM_surfYZ(i + 1, j, k) == true) {
		//mx_fwd = mx_center; my_fwd = my_center; mz_fwd = mz_center;
		dmxdx_fwd = 0.;	// -D / 2. / A * mz_center;
		dmydx_fwd = 0.;
		dmzdx_fwd = 0.; // D / 2. / A * mx_center;
	}
	else if (AFM_surfYZ(i + 1, j, k) == false) {
		if (i == nx - 1) {
			mx_fwd = mx_AFM2_glb_store(0, j, k);
			my_fwd = my_AFM2_glb_store(0, j, k);
			mz_fwd = mz_AFM2_glb_store(0, j, k);
		}
		else {
			mx_fwd = mx_AFM2_glb_store(i + 1, j, k);
			my_fwd = my_AFM2_glb_store(i + 1, j, k);
			mz_fwd = mz_AFM2_glb_store(i + 1, j, k);
		}
		dmxdx_fwd = (mx_fwd - mx_center) / dx;
		dmydx_fwd = (my_fwd - my_center) / dx;
		dmzdx_fwd = (mz_fwd - mz_center) / dx;
	}

	dmxdx = (dmxdx_bwd + dmxdx_fwd) * 0.5;
	dmzdx = (dmzdx_bwd + dmzdx_fwd) * 0.5;
	dmxdx2_f = dmxdx_fwd / dx;
	dmydx2_f = dmydx_fwd / dx;
	dmzdx2_f = dmzdx_fwd / dx;
	dmxdx2_b = -dmxdx_bwd / dx;
	dmydx2_b = -dmydx_bwd / dx;
	dmzdx2_b = -dmzdx_bwd / dx;
	//--------------END X direction (YZ surface)------------//

	//--------------Y direction (XZ surface)------------//
	//------------- negative Y direction---------------//
	if (AFM_surfXZ(i, j, k) == true) {
		//mx_bwd = mx_center; my_bwd = my_center; mz_bwd = mz_center;
		dmxdy_bwd = 0.;
		dmydy_bwd = 0.; // -D / 2. / A * mz_center;
		dmzdy_bwd = 0.; // D / 2. / A * my_center;
	}
	else if (AFM_surfXZ(i, j, k) == false) {
		if (j == 0) {
			mx_bwd = mx_AFM2_glb_store(i, ny - 1, k);
			my_bwd = my_AFM2_glb_store(i, ny - 1, k);
			mz_bwd = mz_AFM2_glb_store(i, ny - 1, k);
		}
		else {
			mx_bwd = mx_AFM2_glb_store(i, j - 1, k);
			my_bwd = my_AFM2_glb_store(i, j - 1, k);
			mz_bwd = mz_AFM2_glb_store(i, j - 1, k);
		}

		dmxdy_bwd = (mx_center - mx_bwd) / dy;
		dmydy_bwd = (my_center - my_bwd) / dy;
		dmzdy_bwd = (mz_center - mz_bwd) / dy;
	}
	//------------- positive Y direction---------------//
	if (AFM_surfXZ(i, j + 1, k) == true) {
		//mx_fwd = mx_center; my_fwd = my_center; mz_fwd = mz_center;
		dmxdy_fwd = 0.;
		dmydy_fwd = 0.; // -D / 2. / A * mz_center;
		dmzdy_fwd = 0.; // D / 2. / A * my_center;
	}
	else if (AFM_surfXZ(i, j + 1, k) == false) {
		if (j == ny - 1) {
			mx_fwd = mx_AFM2_glb_store(i, 0, k);
			my_fwd = my_AFM2_glb_store(i, 0, k);
			mz_fwd = mz_AFM2_glb_store(i, 0, k);
		}
		else {
			mx_fwd = mx_AFM2_glb_store(i, j + 1, k);
			my_fwd = my_AFM2_glb_store(i, j + 1, k);
			mz_fwd = mz_AFM2_glb_store(i, j + 1, k);
		}

		dmxdy_fwd = (mx_fwd - mx_center) / dy;
		dmydy_fwd = (my_fwd - my_center) / dy;
		dmzdy_fwd = (mz_fwd - mz_center) / dy;
	}

	dmydy = (dmydy_bwd + dmydy_fwd) * 0.5;
	dmzdy = (dmzdy_bwd + dmzdy_fwd) * 0.5;
	dmxdy2_f = dmxdy_fwd / dy;
	dmydy2_f = dmydy_fwd / dy;
	dmzdy2_f = dmzdy_fwd / dy;
	dmxdy2_b = -dmxdy_bwd / dy;
	dmydy2_b = -dmydy_bwd / dy;
	dmzdy2_b = -dmzdy_bwd / dy;
	//--------------END Y direction (XZ surface)------------//	

	//--------------Z direction (XY surface)------------//
	//------------- negative Z direction---------------//
	if (AFM_surfXY(i, j, k) == true) {
		//mx_bwd = mx_center; my_bwd = my_center; mz_bwd = mz_center;
		dmxdz_bwd = 0.; dmydz_bwd = 0.; dmzdz_bwd = 0.;
	}
	else if (AFM_surfXY(i, j, k) == false) {
		if (k == 0) {
			mx_bwd = mx_AFM2_glb_store(i, j, nz - 1);
			my_bwd = my_AFM2_glb_store(i, j, nz - 1);
			mz_bwd = mz_AFM2_glb_store(i, j, nz - 1);
		}
		else {
			mx_bwd = mx_AFM2_glb_store(i, j, k - 1);
			my_bwd = my_AFM2_glb_store(i, j, k - 1);
			mz_bwd = mz_AFM2_glb_store(i, j, k - 1);
		}

		dmxdz_bwd = (mx_center - mx_bwd) / dz;
		dmydz_bwd = (my_center - my_bwd) / dz;
		dmzdz_bwd = (mz_center - mz_bwd) / dz;
	}
	//------------- positive Z direction---------------//
	if (AFM_surfXY(i, j, k + 1) == true) {
		//mx_fwd = mx_center; my_fwd = my_center; mz_fwd = mz_center;
		dmxdz_fwd = 0.; dmydz_fwd = 0.; dmzdz_fwd = 0.;
	}
	else if (AFM_surfXY(i, j, k + 1) == false) {
		if (k == nz - 1) {
			mx_fwd = mx_AFM2_glb_store(i, j, 0);
			my_fwd = my_AFM2_glb_store(i, j, 0);
			mz_fwd = mz_AFM2_glb_store(i, j, 0);
		}
		else {
			mx_fwd = mx_AFM2_glb_store(i, j, k + 1);
			my_fwd = my_AFM2_glb_store(i, j, k + 1);
			mz_fwd = mz_AFM2_glb_store(i, j, k + 1);
		}

		dmxdz_fwd = (mx_fwd - mx_center) / dz;
		dmydz_fwd = (my_fwd - my_center) / dz;
		dmzdz_fwd = (mz_fwd - mz_center) / dz;
	}

	dmxdz2_f = dmxdz_fwd / dz;
	dmydz2_f = dmydz_fwd / dz;
	dmzdz2_f = dmzdz_fwd / dz;
	dmxdz2_b = -dmxdz_bwd / dz;
	dmydz2_b = -dmydz_bwd / dz;
	dmzdz2_b = -dmzdz_bwd / dz;
	//--------------END Z direction (XY surface)------------//

}

void magnetic_system::get_averagem_AFM() {
	double ax1, ay1, az1;
	double ax2, ay2, az2;
	ax1 = 0.; ay1 = 0.; az1 = 0.;
	ax2 = 0.; ay2 = 0.; az2 = 0.;

#pragma acc parallel default(present)
	{
#pragma acc loop gang vector reduction(+:ax1)
		for (long int id = 0; id < n; id++) {
			ax1 = ax1 + mx_AFM1_glb(id);
		}

#pragma acc loop gang vector reduction(+:ay1)
		for (long int id = 0; id < n; id++) {
			ay1 = ay1 + my_AFM1_glb(id);
		}

#pragma acc loop gang vector reduction(+:az1)
		for (long int id = 0; id < n; id++) {
			az1 = az1 + mz_AFM1_glb(id);
		}

#pragma acc loop gang vector reduction(+:ax2)
		for (long int id = 0; id < n; id++) {
			ax2 = ax2 + mx_AFM2_glb(id);
		}

#pragma acc loop gang vector reduction(+:ay2)
		for (long int id = 0; id < n; id++) {
			ay2 = ay2 + my_AFM2_glb(id);
		}

#pragma acc loop gang vector reduction(+:az2)
		for (long int id = 0; id < n; id++) {
			az2 = az2 + mz_AFM2_glb(id);
		}
	}

	mx_AFM1_ave = ax1 / static_cast<double>(pt_glb->NAFM);
	my_AFM1_ave = ay1 / static_cast<double>(pt_glb->NAFM);
	mz_AFM1_ave = az1 / static_cast<double>(pt_glb->NAFM);
	mx_AFM2_ave = ax2 / static_cast<double>(pt_glb->NAFM);
	my_AFM2_ave = ay2 / static_cast<double>(pt_glb->NAFM);
	mz_AFM2_ave = az2 / static_cast<double>(pt_glb->NAFM);
}
