#include "magnetic_system.h"

void magnetic_system::get_dm_RK1() {
	double Hx_anis_crt, Hy_anis_crt, Hz_anis_crt;
	double Hx_anis_glb, Hy_anis_glb, Hz_anis_glb;
	double Hx_elas_crt, Hy_elas_crt, Hz_elas_crt;
	double Hx_elas_glb, Hy_elas_glb, Hz_elas_glb;
	double Hx_exch_glb, Hy_exch_glb, Hz_exch_glb;
	double Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb;

	double Hx_eff_glb, Hy_eff_glb, Hz_eff_glb;
	double precess_termx, precess_termy, precess_termz;
	double damping_termx, damping_termy, damping_termz;

	long int x, y, z;
	unsigned int mat_type;
	material* mat;
	double Ms, K1, K2, B1, B2;
	double Aexi, iDMIi;
	double alpha, gyro, precess;
	double mx_glb_local, my_glb_local, mz_glb_local;
	double mx, my, mz;
	double exx, eyy, ezz, eyz, exz, exy;
	//double d2mxdx2, d2mxdy2, d2mxdz2;
	//double d2mydx2, d2mydy2, d2mydz2;
	//double d2mzdx2, d2mzdy2, d2mzdz2;
	double dmxdx, dmydy, dmzdx, dmzdy;
	double dmxdx2_f, dmydx2_f, dmzdx2_f;
	double dmxdx2_b, dmydx2_b, dmzdx2_b;
	double dmxdy2_f, dmydy2_f, dmzdy2_f;
	double dmxdy2_b, dmydy2_b, dmzdy2_b;
	double dmxdz2_f, dmydz2_f, dmzdz2_f;
	double dmxdz2_b, dmydz2_b, dmzdz2_b;

	if (pt_glb->if_prescribe_Hext == false) {
#pragma acc serial default(present) async(4)
		{
			update_external_Hfield();
		}
	}

#pragma acc parallel default(present) async(4)
	{
#pragma acc loop gang vector private(Hx_anis_crt, Hy_anis_crt, Hz_anis_crt,\
	Hx_anis_glb, Hy_anis_glb, Hz_anis_glb,\
	Hx_elas_crt, Hy_elas_crt, Hz_elas_crt,\
	Hx_elas_glb, Hy_elas_glb, Hz_elas_glb,\
	Hx_exch_glb, Hy_exch_glb, Hz_exch_glb,\
	Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb,\
	Hx_eff_glb, Hy_eff_glb, Hz_eff_glb,\
	precess_termx, precess_termy, precess_termz,\
	damping_termx, damping_termy, damping_termz,\
	x,y,z,\
	mat_type,\
	Ms, K1, K2, B1, B2, \
	Aexi, iDMIi, \
	alpha, gyro, precess,\
	mx_glb_local, my_glb_local, mz_glb_local,\
	mx, my, mz,\
	exx, eyy, ezz, eyz, exz, exy, mat,\
	dmxdx, dmydy, dmzdx, dmzdy, \
	dmxdx2_f, dmydx2_f, dmzdx2_f,\
	dmxdx2_b, dmydx2_b, dmzdx2_b,\
	dmxdy2_f, dmydy2_f, dmzdy2_f,\
	dmxdy2_b, dmydy2_b, dmzdy2_b,\
	dmxdz2_f, dmydz2_f, dmzdz2_f,\
	dmxdz2_b, dmydz2_b, dmzdz2_b)
		for (long int i = 0; i < n; i++) {
			mat_type = pt_glb->material_cell(i);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_FM == true) {
					x = i / (ny * nz);
					y = (i - x * (ny * nz)) / nz;
					z = i - x * (ny * nz) - y * nz;

					mx_glb_local = mx_glb_store(i);
					my_glb_local = my_glb_store(i);
					mz_glb_local = mz_glb_store(i);

					pt_math->transform_vector_glb2crt(mx_glb_local, my_glb_local, mz_glb_local, \
						mx, my, mz);

					Ms = mat->Ms;
					K1 = mat->K1 / mu0 / Ms; K2 = mat->K2 / mu0 / Ms;
					B1 = mat->B1 / mu0 / Ms; B2 = mat->B2 / mu0 / Ms;
					Aexi = mat->Aex; iDMIi = mat->iDMI;

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

					get_m_spatial_derivative_local \
						(x, y, z, \
							Aexi, iDMIi, \
							dmxdx, dmydy, dmzdx, dmzdy, \
							dmxdx2_f, dmydx2_f, dmzdx2_f, \
							dmxdx2_b, dmydx2_b, dmzdx2_b, \
							dmxdy2_f, dmydy2_f, dmzdy2_f, \
							dmxdy2_b, dmydy2_b, dmzdy2_b, \
							dmxdz2_f, dmydz2_f, dmzdz2_f, \
							dmxdz2_b, dmydz2_b, dmzdz2_b);

					//get_laplacian_m_local \
						//	(x, y, z, \
						//		d2mxdx2, d2mxdy2, d2mxdz2, \
						//		d2mydx2, d2mydy2, d2mydz2, \
						//		d2mzdx2, d2mzdy2, d2mzdz2);

						//------------Magnetocrystalline anisotropy field--------------//
					if (mat->anisotropy_type == 1) {
						Hx_anis_crt = -2. * \
							(K1 * (my * my + mz * mz) \
								+ K2 * my * my * mz * mz) * mx;
						Hy_anis_crt = -2. * \
							(K1 * (mx * mx + mz * mz) \
								+ K2 * mx * mx * mz * mz) * my;
						Hz_anis_crt = -2. * \
							(K1 * (mx * mx + my * my) \
								+ K2 * mx * mx * my * my) * mz;
					}
					else if (mat->anisotropy_type == 2) {
						Hx_anis_crt = 0.;
						Hy_anis_crt = 0.;
						Hz_anis_crt = 2. * K1 * mz;
					}
					pt_math->transform_vector_crt2glb(Hx_anis_crt, Hy_anis_crt, Hz_anis_crt, \
						Hx_anis_glb, Hy_anis_glb, Hz_anis_glb);

					//Hx_anis(i) = Hx_anis_glb;
					//Hy_anis(i) = Hy_anis_glb;
					//Hz_anis(i) = Hz_anis_glb;

					//------------------Magnetoelastic effective field---------------//
					Hx_elas_crt = -2. * (B1 * (exx)*mx + \
						B2 * ((exy)*my + (exz)*mz));
					Hy_elas_crt = -2. * (B1 * (eyy)*my + \
						B2 * ((exy)*mx + (eyz)*mz));
					Hz_elas_crt = -2. * (B1 * (ezz)*mz + \
						B2 * ((exz)*mx + (eyz)*my));

					pt_math->transform_vector_crt2glb(Hx_elas_crt, Hy_elas_crt, Hz_elas_crt, \
						Hx_elas_glb, Hy_elas_glb, Hz_elas_glb);

					//Hx_elas(i) = Hx_elas_glb;
					//Hy_elas(i) = Hy_elas_glb;
					//Hz_elas(i) = Hz_elas_glb;

					//-----------------Magnetic exchange field----------------------//
					//----------Averaged parameter for inter-region exchange - Approach 1---------//
					Hx_exch_glb = 2. / mu0 / Ms * \
						(Aij_xf(i) * dmxdx2_f + Aij_xb(i) * dmxdx2_b + Aij_yf(i) * dmxdy2_f + Aij_yb(i) * dmxdy2_b + Aij_zf(i) * dmxdz2_f + Aij_zb(i) * dmxdz2_b);
					Hy_exch_glb = 2. / mu0 / Ms * \
						(Aij_xf(i) * dmydx2_f + Aij_xb(i) * dmydx2_b + Aij_yf(i) * dmydy2_f + Aij_yb(i) * dmydy2_b + Aij_zf(i) * dmydz2_f + Aij_zb(i) * dmydz2_b);
					Hz_exch_glb = 2. / mu0 / Ms * \
						(Aij_xf(i) * dmzdx2_f + Aij_xb(i) * dmzdx2_b + Aij_yf(i) * dmzdy2_f + Aij_yb(i) * dmzdy2_b + Aij_zf(i) * dmzdz2_f + Aij_zb(i) * dmzdz2_b);
					//----------Averaged parameter for inter-region exchange - Approach 2---------//
					//Hx_exch_glb = 2. / mu0 * \
					//	AM_exchange(i) * (dmxdx2_f + dmxdx2_b + dmxdy2_f + dmxdy2_b + dmxdz2_f + dmxdz2_b);
					//Hy_exch_glb = 2. / mu0 * \
					//	AM_exchange(i) * (dmydx2_f + dmydx2_b + dmydy2_f + dmydy2_b + dmydz2_f + dmydz2_b);
					//Hz_exch_glb = 2. / mu0 * \
					//	AM_exchange(i) * (dmzdx2_f + dmzdx2_b + dmzdy2_f + dmzdy2_b + dmzdz2_f + dmzdz2_b);

					//Hx_exch(i) = 2. * Aex * (d2mxdx2 + d2mxdy2 + d2mxdz2);
					//Hy_exch(i) = 2. * Aex * (d2mydx2 + d2mydy2 + d2mydz2);
					//Hz_exch(i) = 2. * Aex * (d2mzdx2 + d2mzdy2 + d2mzdz2);
					// 
					//------------------Interfacial DMI effective field-------------------//
					Hx_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdx;
					Hy_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdy;
					Hz_iDMI_glb = -2. * iDMIi / mu0 / Ms * (dmxdx + dmydy);

					//-------------------LLG Equation------------------------------//
					Hx_eff_glb = Hx_stat(i) + \
						Hx_anis_glb + Hx_elas_glb + \
						Hx_exch_glb + Hx_iDMI_glb + pt_glb->Hext[0] + pt_glb->Hextx_nonunif(i);

					Hy_eff_glb = Hy_stat(i) + \
						Hy_anis_glb + Hy_elas_glb + \
						Hy_exch_glb + Hy_iDMI_glb + pt_glb->Hext[1] + pt_glb->Hexty_nonunif(i);

					Hz_eff_glb = Hz_stat(i) + \
						Hz_anis_glb + Hz_elas_glb + \
						Hz_exch_glb + Hz_iDMI_glb + pt_glb->Hext[2] + pt_glb->Hextz_nonunif(i);

					if (pt_glb->if_EM_backaction_on_M == true) {
						Hx_eff_glb = Hx_eff_glb + pt_EM->DHx_em_cell(i);
						Hy_eff_glb = Hy_eff_glb + pt_EM->DHy_em_cell(i);
						Hz_eff_glb = Hz_eff_glb + pt_EM->DHz_em_cell(i);
					}

					alpha = mat->FM_damping;
					gyro = mat->gyro;
					precess = -gyro / (1 + alpha * alpha) * pt_glb->dt;

					precess_termx = my_glb_local * Hz_eff_glb - mz_glb_local * Hy_eff_glb;
					precess_termy = mz_glb_local * Hx_eff_glb - mx_glb_local * Hz_eff_glb;
					precess_termz = mx_glb_local * Hy_eff_glb - my_glb_local * Hx_eff_glb;

					damping_termx = my_glb_local * precess_termz - mz_glb_local * precess_termy;
					damping_termy = mz_glb_local * precess_termx - mx_glb_local * precess_termz;
					damping_termz = mx_glb_local * precess_termy - my_glb_local * precess_termx;

					dmx_glb_rk1(i) = precess * precess_termx + alpha * precess * damping_termx;
					dmy_glb_rk1(i) = precess * precess_termy + alpha * precess * damping_termy;
					dmz_glb_rk1(i) = precess * precess_termz + alpha * precess * damping_termz;
				}
			}
		}
	}

	if (pt_glb->if_spin_pumping == true) {
		get_J_ISHE_RK1();
	}
}

void magnetic_system::get_J_ISHE_RK1() {
	long i, j, sz;
	unsigned int mat_type;
	material* mat;
	double coeff_Js, coeff_L;
	double distance;

#pragma acc parallel default(present) async(4)
	{
#pragma acc loop gang vector private(i,j,mat_type, mat, coeff_Js,coeff_L, sz, distance)
		for (long k = 0; k < nz; k++) {

			if (pt_glb->if_J_ISHE[k] == true) {
				mat_type = pt_glb->material_cell(0, 0, k);
				mat = &(pt_glb->material_parameters[(mat_type)-1]);

				sz = pt_glb->source_z[k];

				coeff_Js = mat->spin_hall_angle * unit_e / 2. / PI * mat->spin_mix_cond;

				distance = (abs(static_cast<double>((k - sz))) - 0.5) * dz;

				if (k > sz) {
					coeff_L = sinh((pt_glb->thickness_pumping_layer[k] - distance) / mat->spin_diffus_length) \
						/ sinh(pt_glb->thickness_pumping_layer[k] / mat->spin_diffus_length);
				}
				else {
					coeff_L = -1. * sinh((pt_glb->thickness_pumping_layer[k] - distance) / mat->spin_diffus_length) \
						/ sinh(pt_glb->thickness_pumping_layer[k] / mat->spin_diffus_length);
				}

#pragma acc loop seq
				for (i = 0; i < nx; i++) {
#pragma acc loop seq
					for (j = 0; j < ny; j++) {
						Jx_ISHE(i, j, k) = coeff_Js * coeff_L * \
							(mx_glb_store(i, j, sz) * dmz_glb_rk1(i, j, sz) - mz_glb_store(i, j, sz) * dmx_glb_rk1(i, j, sz)\
								+ mx_AFM1_glb_store(i, j, sz) * dmz_AFM1_glb_rk1(i, j, sz) - mz_AFM1_glb_store(i, j, sz) * dmx_AFM1_glb_rk1(i, j, sz)\
								+ mx_AFM2_glb_store(i, j, sz) * dmz_AFM2_glb_rk1(i, j, sz) - mz_AFM2_glb_store(i, j, sz) * dmx_AFM2_glb_rk1(i, j, sz)) \
							/ pt_glb->dt;

						Jy_ISHE(i, j, k) = coeff_Js * coeff_L * \
							(my_glb_store(i, j, sz) * dmz_glb_rk1(i, j, sz) - mz_glb_store(i, j, sz) * dmy_glb_rk1(i, j, sz)\
								+ my_AFM1_glb_store(i, j, sz) * dmz_AFM1_glb_rk1(i, j, sz) - mz_AFM1_glb_store(i, j, sz) * dmy_AFM1_glb_rk1(i, j, sz)\
								+ my_AFM2_glb_store(i, j, sz) * dmz_AFM2_glb_rk1(i, j, sz) - mz_AFM2_glb_store(i, j, sz) * dmy_AFM2_glb_rk1(i, j, sz)) \
							/ pt_glb->dt;
					}
				}
			}

		}
	}
}

void magnetic_system::get_dm_RK2() {
	double Hx_anis_crt, Hy_anis_crt, Hz_anis_crt;
	double Hx_anis_glb, Hy_anis_glb, Hz_anis_glb;
	double Hx_elas_crt, Hy_elas_crt, Hz_elas_crt;
	double Hx_elas_glb, Hy_elas_glb, Hz_elas_glb;
	double Hx_exch_glb, Hy_exch_glb, Hz_exch_glb;
	double Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb;

	double Hx_eff_glb, Hy_eff_glb, Hz_eff_glb;
	double precess_termx, precess_termy, precess_termz;
	double damping_termx, damping_termy, damping_termz;

	long int x, y, z;
	unsigned int mat_type;
	material* mat;
	double Ms, K1, K2, B1, B2;
	double Aexi, iDMIi;
	double alpha, gyro, precess;
	double mx_glb_local, my_glb_local, mz_glb_local;
	double mx, my, mz;
	double exx, eyy, ezz, eyz, exz, exy;
	//double d2mxdx2, d2mxdy2, d2mxdz2;
	//double d2mydx2, d2mydy2, d2mydz2;
	//double d2mzdx2, d2mzdy2, d2mzdz2;
	double dmxdx, dmydy, dmzdx, dmzdy;
	double dmxdx2_f, dmydx2_f, dmzdx2_f;
	double dmxdx2_b, dmydx2_b, dmzdx2_b;
	double dmxdy2_f, dmydy2_f, dmzdy2_f;
	double dmxdy2_b, dmydy2_b, dmzdy2_b;
	double dmxdz2_f, dmydz2_f, dmzdz2_f;
	double dmxdz2_b, dmydz2_b, dmzdz2_b;

	if (pt_glb->if_prescribe_Hext == false) {
#pragma acc serial default(present) async(4)
		{
			update_external_Hfield_half();
		}
	}

#pragma acc parallel default(present) async(4)
	{
#pragma acc loop gang vector private(Hx_anis_crt, Hy_anis_crt, Hz_anis_crt,\
	Hx_anis_glb, Hy_anis_glb, Hz_anis_glb,\
	Hx_elas_crt, Hy_elas_crt, Hz_elas_crt,\
	Hx_elas_glb, Hy_elas_glb, Hz_elas_glb,\
	Hx_exch_glb, Hy_exch_glb, Hz_exch_glb,\
	Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb,\
	Hx_eff_glb, Hy_eff_glb, Hz_eff_glb,\
	precess_termx, precess_termy, precess_termz,\
	damping_termx, damping_termy, damping_termz,\
	x,y,z,\
	mat_type,\
	Ms, K1, K2, B1, B2, \
	Aexi, iDMIi, \
	alpha, gyro, precess,\
	mx_glb_local, my_glb_local, mz_glb_local,\
	mx, my, mz,\
	exx, eyy, ezz, eyz, exz, exy, mat,\
	dmxdx, dmydy, dmzdx, dmzdy, \
	dmxdx2_f, dmydx2_f, dmzdx2_f,\
	dmxdx2_b, dmydx2_b, dmzdx2_b,\
	dmxdy2_f, dmydy2_f, dmzdy2_f,\
	dmxdy2_b, dmydy2_b, dmzdy2_b,\
	dmxdz2_f, dmydz2_f, dmzdz2_f,\
	dmxdz2_b, dmydz2_b, dmzdz2_b)
		for (long int i = 0; i < n; i++) {
			mat_type = pt_glb->material_cell(i);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_FM == true) {
					x = i / (ny * nz);
					y = (i - x * (ny * nz)) / nz;
					z = i - x * (ny * nz) - y * nz;

					mx_glb_local = mx_glb_store(i);
					my_glb_local = my_glb_store(i);
					mz_glb_local = mz_glb_store(i);

					pt_math->transform_vector_glb2crt(mx_glb_local, my_glb_local, mz_glb_local, \
						mx, my, mz);

					Ms = mat->Ms;
					K1 = mat->K1 / mu0 / Ms; K2 = mat->K2 / mu0 / Ms;
					B1 = mat->B1 / mu0 / Ms; B2 = mat->B2 / mu0 / Ms;
					Aexi = mat->Aex; iDMIi = mat->iDMI;

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

					get_m_spatial_derivative_local \
						(x, y, z, \
							Aexi, iDMIi, \
							dmxdx, dmydy, dmzdx, dmzdy, \
							dmxdx2_f, dmydx2_f, dmzdx2_f, \
							dmxdx2_b, dmydx2_b, dmzdx2_b, \
							dmxdy2_f, dmydy2_f, dmzdy2_f, \
							dmxdy2_b, dmydy2_b, dmzdy2_b, \
							dmxdz2_f, dmydz2_f, dmzdz2_f, \
							dmxdz2_b, dmydz2_b, dmzdz2_b);

					//get_laplacian_m_local \
						//	(x, y, z, \
						//		d2mxdx2, d2mxdy2, d2mxdz2, \
						//		d2mydx2, d2mydy2, d2mydz2, \
						//		d2mzdx2, d2mzdy2, d2mzdz2);

						//------------Magnetocrystalline anisotropy field--------------//
					if (mat->anisotropy_type == 1) {
						Hx_anis_crt = -2. * \
							(K1 * (my * my + mz * mz) \
								+ K2 * my * my * mz * mz) * mx;
						Hy_anis_crt = -2. * \
							(K1 * (mx * mx + mz * mz) \
								+ K2 * mx * mx * mz * mz) * my;
						Hz_anis_crt = -2. * \
							(K1 * (mx * mx + my * my) \
								+ K2 * mx * mx * my * my) * mz;
					}
					else if (mat->anisotropy_type == 2) {
						Hx_anis_crt = 0.;
						Hy_anis_crt = 0.;
						Hz_anis_crt = 2. * K1 * mz;
					}
					pt_math->transform_vector_crt2glb(Hx_anis_crt, Hy_anis_crt, Hz_anis_crt, \
						Hx_anis_glb, Hy_anis_glb, Hz_anis_glb);

					//Hx_anis(i) = Hx_anis_glb;
					//Hy_anis(i) = Hy_anis_glb;
					//Hz_anis(i) = Hz_anis_glb;

					//------------------Magnetoelastic effective field---------------//
					Hx_elas_crt = -2. * (B1 * (exx)*mx + \
						B2 * ((exy)*my + (exz)*mz));
					Hy_elas_crt = -2. * (B1 * (eyy)*my + \
						B2 * ((exy)*mx + (eyz)*mz));
					Hz_elas_crt = -2. * (B1 * (ezz)*mz + \
						B2 * ((exz)*mx + (eyz)*my));

					pt_math->transform_vector_crt2glb(Hx_elas_crt, Hy_elas_crt, Hz_elas_crt, \
						Hx_elas_glb, Hy_elas_glb, Hz_elas_glb);

					//Hx_elas(i) = Hx_elas_glb;
					//Hy_elas(i) = Hy_elas_glb;
					//Hz_elas(i) = Hz_elas_glb;

					//-----------------Magnetic exchange field----------------------//
					//----------Averaged parameter for inter-region exchange - Approach 1---------//
					Hx_exch_glb = 2. / mu0 / Ms * \
						(Aij_xf(i) * dmxdx2_f + Aij_xb(i) * dmxdx2_b + Aij_yf(i) * dmxdy2_f + Aij_yb(i) * dmxdy2_b + Aij_zf(i) * dmxdz2_f + Aij_zb(i) * dmxdz2_b);
					Hy_exch_glb = 2. / mu0 / Ms * \
						(Aij_xf(i) * dmydx2_f + Aij_xb(i) * dmydx2_b + Aij_yf(i) * dmydy2_f + Aij_yb(i) * dmydy2_b + Aij_zf(i) * dmydz2_f + Aij_zb(i) * dmydz2_b);
					Hz_exch_glb = 2. / mu0 / Ms * \
						(Aij_xf(i) * dmzdx2_f + Aij_xb(i) * dmzdx2_b + Aij_yf(i) * dmzdy2_f + Aij_yb(i) * dmzdy2_b + Aij_zf(i) * dmzdz2_f + Aij_zb(i) * dmzdz2_b);

					//----------Averaged parameter for inter-region exchange - Approach 2---------//
					//Hx_exch_glb = 2. / mu0 * \
					//	AM_exchange(i) * (dmxdx2_f + dmxdx2_b + dmxdy2_f + dmxdy2_b + dmxdz2_f + dmxdz2_b);
					//Hy_exch_glb = 2. / mu0 * \
					//	AM_exchange(i) * (dmydx2_f + dmydx2_b + dmydy2_f + dmydy2_b + dmydz2_f + dmydz2_b);
					//Hz_exch_glb = 2. / mu0 * \
					//	AM_exchange(i) * (dmzdx2_f + dmzdx2_b + dmzdy2_f + dmzdy2_b + dmzdz2_f + dmzdz2_b);

					//Hx_exch(i) = 2. * Aex * (d2mxdx2 + d2mxdy2 + d2mxdz2);
					//Hy_exch(i) = 2. * Aex * (d2mydx2 + d2mydy2 + d2mydz2);
					//Hz_exch(i) = 2. * Aex * (d2mzdx2 + d2mzdy2 + d2mzdz2);

					//------------------Interfacial DMI effective field-------------------//
					Hx_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdx;
					Hy_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdy;
					Hz_iDMI_glb = -2. * iDMIi / mu0 / Ms * (dmxdx + dmydy);

					//-------------------LLG Equation------------------------------//
					Hx_eff_glb = Hx_stat(i) + \
						Hx_anis_glb + Hx_elas_glb + \
						Hx_exch_glb + Hx_iDMI_glb + pt_glb->Hext[0] + pt_glb->Hextx_nonunif(i);

					Hy_eff_glb = Hy_stat(i) + \
						Hy_anis_glb + Hy_elas_glb + \
						Hy_exch_glb + Hy_iDMI_glb + pt_glb->Hext[1] + pt_glb->Hexty_nonunif(i);

					Hz_eff_glb = Hz_stat(i) + \
						Hz_anis_glb + Hz_elas_glb + \
						Hz_exch_glb + Hz_iDMI_glb + pt_glb->Hext[2] + pt_glb->Hextz_nonunif(i);

					if (pt_glb->if_EM_backaction_on_M == true) {
						Hx_eff_glb = Hx_eff_glb + pt_EM->DHx_em_cell(i);
						Hy_eff_glb = Hy_eff_glb + pt_EM->DHy_em_cell(i);
						Hz_eff_glb = Hz_eff_glb + pt_EM->DHz_em_cell(i);
					}

					alpha = mat->FM_damping;
					gyro = mat->gyro;
					precess = -gyro / (1 + alpha * alpha) * pt_glb->dt;

					precess_termx = my_glb_local * Hz_eff_glb - mz_glb_local * Hy_eff_glb;
					precess_termy = mz_glb_local * Hx_eff_glb - mx_glb_local * Hz_eff_glb;
					precess_termz = mx_glb_local * Hy_eff_glb - my_glb_local * Hx_eff_glb;

					damping_termx = my_glb_local * precess_termz - mz_glb_local * precess_termy;
					damping_termy = mz_glb_local * precess_termx - mx_glb_local * precess_termz;
					damping_termz = mx_glb_local * precess_termy - my_glb_local * precess_termx;

					dmx_glb_rk2(i) = precess * precess_termx + alpha * precess * damping_termx;
					dmy_glb_rk2(i) = precess * precess_termy + alpha * precess * damping_termy;
					dmz_glb_rk2(i) = precess * precess_termz + alpha * precess * damping_termz;
				}
			}
		}
	}

	if (pt_glb->if_spin_pumping == true) {
		get_J_ISHE_RK2();
	}
}

void magnetic_system::get_J_ISHE_RK2() {
	long i, j, sz;
	unsigned int mat_type;
	material* mat;
	double coeff_Js, coeff_L;
	double distance;

#pragma acc parallel default(present) async(4)
	{
#pragma acc loop gang vector private(i,j,mat_type, mat, coeff_Js,coeff_L, sz, distance)
		for (long k = 0; k < nz; k++) {

			if (pt_glb->if_J_ISHE[k] == true) {
				mat_type = pt_glb->material_cell(0, 0, k);
				mat = &(pt_glb->material_parameters[(mat_type)-1]);

				sz = pt_glb->source_z[k];

				coeff_Js = mat->spin_hall_angle * unit_e / 2. / PI * mat->spin_mix_cond;

				distance = (abs(static_cast<double>((k - sz))) - 0.5) * dz;

				if (k > sz) {
					coeff_L = sinh((pt_glb->thickness_pumping_layer[k] - distance) / mat->spin_diffus_length) \
						/ sinh(pt_glb->thickness_pumping_layer[k] / mat->spin_diffus_length);
				}
				else {
					coeff_L = -1. * sinh((pt_glb->thickness_pumping_layer[k] - distance) / mat->spin_diffus_length) \
						/ sinh(pt_glb->thickness_pumping_layer[k] / mat->spin_diffus_length);
				}

#pragma acc loop seq
				for (i = 0; i < nx; i++) {
#pragma acc loop seq
					for (j = 0; j < ny; j++) {
						Jx_ISHE(i, j, k) = coeff_Js * coeff_L * \
							(mx_glb_store(i, j, sz) * dmz_glb_rk2(i, j, sz) - mz_glb_store(i, j, sz) * dmx_glb_rk2(i, j, sz)\
								+ mx_AFM1_glb_store(i, j, sz) * dmz_AFM1_glb_rk2(i, j, sz) - mz_AFM1_glb_store(i, j, sz) * dmx_AFM1_glb_rk2(i, j, sz)\
								+ mx_AFM2_glb_store(i, j, sz) * dmz_AFM2_glb_rk2(i, j, sz) - mz_AFM2_glb_store(i, j, sz) * dmx_AFM2_glb_rk2(i, j, sz)) \
							/ pt_glb->dt;

						Jy_ISHE(i, j, k) = coeff_Js * coeff_L * \
							(my_glb_store(i, j, sz) * dmz_glb_rk2(i, j, sz) - mz_glb_store(i, j, sz) * dmy_glb_rk2(i, j, sz)\
								+ my_AFM1_glb_store(i, j, sz) * dmz_AFM1_glb_rk2(i, j, sz) - mz_AFM1_glb_store(i, j, sz) * dmy_AFM1_glb_rk2(i, j, sz)\
								+ my_AFM2_glb_store(i, j, sz) * dmz_AFM2_glb_rk2(i, j, sz) - mz_AFM2_glb_store(i, j, sz) * dmy_AFM2_glb_rk2(i, j, sz)) \
							/ pt_glb->dt;
					}
				}
			}

		}
	}
}

void magnetic_system::get_dm_RK3() {
	double Hx_anis_crt, Hy_anis_crt, Hz_anis_crt;
	double Hx_anis_glb, Hy_anis_glb, Hz_anis_glb;
	double Hx_elas_crt, Hy_elas_crt, Hz_elas_crt;
	double Hx_elas_glb, Hy_elas_glb, Hz_elas_glb;
	double Hx_exch_glb, Hy_exch_glb, Hz_exch_glb;
	double Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb;

	double Hx_eff_glb, Hy_eff_glb, Hz_eff_glb;
	double precess_termx, precess_termy, precess_termz;
	double damping_termx, damping_termy, damping_termz;

	long int x, y, z;
	unsigned int mat_type;
	material* mat;
	double Ms, K1, K2, B1, B2;
	double Aexi, iDMIi;
	double alpha, gyro, precess;
	double mx_glb_local, my_glb_local, mz_glb_local;
	double mx, my, mz;
	double exx, eyy, ezz, eyz, exz, exy;
	//double d2mxdx2, d2mxdy2, d2mxdz2;
	//double d2mydx2, d2mydy2, d2mydz2;
	//double d2mzdx2, d2mzdy2, d2mzdz2;
	double dmxdx, dmydy, dmzdx, dmzdy;
	double dmxdx2_f, dmydx2_f, dmzdx2_f;
	double dmxdx2_b, dmydx2_b, dmzdx2_b;
	double dmxdy2_f, dmydy2_f, dmzdy2_f;
	double dmxdy2_b, dmydy2_b, dmzdy2_b;
	double dmxdz2_f, dmydz2_f, dmzdz2_f;
	double dmxdz2_b, dmydz2_b, dmzdz2_b;

	if (pt_glb->if_prescribe_Hext == false) {
#pragma acc serial default(present) async(4)
		{
			update_external_Hfield_half();
		}
	}

#pragma acc parallel default(present) async(4)
	{
#pragma acc loop gang vector private(Hx_anis_crt, Hy_anis_crt, Hz_anis_crt,\
	Hx_anis_glb, Hy_anis_glb, Hz_anis_glb,\
	Hx_elas_crt, Hy_elas_crt, Hz_elas_crt,\
	Hx_elas_glb, Hy_elas_glb, Hz_elas_glb,\
	Hx_exch_glb, Hy_exch_glb, Hz_exch_glb,\
	Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb,\
	Hx_eff_glb, Hy_eff_glb, Hz_eff_glb,\
	precess_termx, precess_termy, precess_termz,\
	damping_termx, damping_termy, damping_termz,\
	x,y,z,\
	mat_type,\
	Ms, K1, K2, B1, B2, \
	Aexi, iDMIi, \
	alpha, gyro, precess,\
	mx_glb_local, my_glb_local, mz_glb_local,\
	mx, my, mz,\
	exx, eyy, ezz, eyz, exz, exy, mat,\
	dmxdx, dmydy, dmzdx, dmzdy, \
	dmxdx2_f, dmydx2_f, dmzdx2_f,\
	dmxdx2_b, dmydx2_b, dmzdx2_b,\
	dmxdy2_f, dmydy2_f, dmzdy2_f,\
	dmxdy2_b, dmydy2_b, dmzdy2_b,\
	dmxdz2_f, dmydz2_f, dmzdz2_f,\
	dmxdz2_b, dmydz2_b, dmzdz2_b)
		for (long int i = 0; i < n; i++) {
			mat_type = pt_glb->material_cell(i);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_FM == true) {
					x = i / (ny * nz);
					y = (i - x * (ny * nz)) / nz;
					z = i - x * (ny * nz) - y * nz;

					mx_glb_local = mx_glb_store(i);
					my_glb_local = my_glb_store(i);
					mz_glb_local = mz_glb_store(i);

					pt_math->transform_vector_glb2crt(mx_glb_local, my_glb_local, mz_glb_local, \
						mx, my, mz);

					Ms = mat->Ms;
					K1 = mat->K1 / mu0 / Ms; K2 = mat->K2 / mu0 / Ms;
					B1 = mat->B1 / mu0 / Ms; B2 = mat->B2 / mu0 / Ms;
					Aexi = mat->Aex; iDMIi = mat->iDMI;

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

					get_m_spatial_derivative_local \
						(x, y, z, \
							Aexi, iDMIi, \
							dmxdx, dmydy, dmzdx, dmzdy, \
							dmxdx2_f, dmydx2_f, dmzdx2_f, \
							dmxdx2_b, dmydx2_b, dmzdx2_b, \
							dmxdy2_f, dmydy2_f, dmzdy2_f, \
							dmxdy2_b, dmydy2_b, dmzdy2_b, \
							dmxdz2_f, dmydz2_f, dmzdz2_f, \
							dmxdz2_b, dmydz2_b, dmzdz2_b);

					//get_laplacian_m_local \
						//	(x, y, z, \
						//		d2mxdx2, d2mxdy2, d2mxdz2, \
						//		d2mydx2, d2mydy2, d2mydz2, \
						//		d2mzdx2, d2mzdy2, d2mzdz2);

						//------------Magnetocrystalline anisotropy field--------------//
					if (mat->anisotropy_type == 1) {
						Hx_anis_crt = -2. * \
							(K1 * (my * my + mz * mz) \
								+ K2 * my * my * mz * mz) * mx;
						Hy_anis_crt = -2. * \
							(K1 * (mx * mx + mz * mz) \
								+ K2 * mx * mx * mz * mz) * my;
						Hz_anis_crt = -2. * \
							(K1 * (mx * mx + my * my) \
								+ K2 * mx * mx * my * my) * mz;
					}
					else if (mat->anisotropy_type == 2) {
						Hx_anis_crt = 0.;
						Hy_anis_crt = 0.;
						Hz_anis_crt = 2. * K1 * mz;
					}
					pt_math->transform_vector_crt2glb(Hx_anis_crt, Hy_anis_crt, Hz_anis_crt, \
						Hx_anis_glb, Hy_anis_glb, Hz_anis_glb);

					//Hx_anis(i) = Hx_anis_glb;
					//Hy_anis(i) = Hy_anis_glb;
					//Hz_anis(i) = Hz_anis_glb;

					//------------------Magnetoelastic effective field---------------//
					Hx_elas_crt = -2. * (B1 * (exx)*mx + \
						B2 * ((exy)*my + (exz)*mz));
					Hy_elas_crt = -2. * (B1 * (eyy)*my + \
						B2 * ((exy)*mx + (eyz)*mz));
					Hz_elas_crt = -2. * (B1 * (ezz)*mz + \
						B2 * ((exz)*mx + (eyz)*my));

					pt_math->transform_vector_crt2glb(Hx_elas_crt, Hy_elas_crt, Hz_elas_crt, \
						Hx_elas_glb, Hy_elas_glb, Hz_elas_glb);

					//Hx_elas(i) = Hx_elas_glb;
					//Hy_elas(i) = Hy_elas_glb;
					//Hz_elas(i) = Hz_elas_glb;

					//-----------------Magnetic exchange field----------------------//
					//----------Averaged parameter for inter-region exchange - Approach 1---------//
					Hx_exch_glb = 2. / mu0 / Ms * \
						(Aij_xf(i) * dmxdx2_f + Aij_xb(i) * dmxdx2_b + Aij_yf(i) * dmxdy2_f + Aij_yb(i) * dmxdy2_b + Aij_zf(i) * dmxdz2_f + Aij_zb(i) * dmxdz2_b);
					Hy_exch_glb = 2. / mu0 / Ms * \
						(Aij_xf(i) * dmydx2_f + Aij_xb(i) * dmydx2_b + Aij_yf(i) * dmydy2_f + Aij_yb(i) * dmydy2_b + Aij_zf(i) * dmydz2_f + Aij_zb(i) * dmydz2_b);
					Hz_exch_glb = 2. / mu0 / Ms * \
						(Aij_xf(i) * dmzdx2_f + Aij_xb(i) * dmzdx2_b + Aij_yf(i) * dmzdy2_f + Aij_yb(i) * dmzdy2_b + Aij_zf(i) * dmzdz2_f + Aij_zb(i) * dmzdz2_b);
					//----------Averaged parameter for inter-region exchange - Approach 2---------//
					//Hx_exch_glb = 2. / mu0 * \
					//	AM_exchange(i) * (dmxdx2_f + dmxdx2_b + dmxdy2_f + dmxdy2_b + dmxdz2_f + dmxdz2_b);
					//Hy_exch_glb = 2. / mu0 * \
					//	AM_exchange(i) * (dmydx2_f + dmydx2_b + dmydy2_f + dmydy2_b + dmydz2_f + dmydz2_b);
					//Hz_exch_glb = 2. / mu0 * \
					//	AM_exchange(i) * (dmzdx2_f + dmzdx2_b + dmzdy2_f + dmzdy2_b + dmzdz2_f + dmzdz2_b);
					// 
					//Hx_exch(i) = 2. * Aex * (d2mxdx2 + d2mxdy2 + d2mxdz2);
					//Hy_exch(i) = 2. * Aex * (d2mydx2 + d2mydy2 + d2mydz2);
					//Hz_exch(i) = 2. * Aex * (d2mzdx2 + d2mzdy2 + d2mzdz2);

					//------------------Interfacial DMI effective field-------------------//
					Hx_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdx;
					Hy_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdy;
					Hz_iDMI_glb = -2. * iDMIi / mu0 / Ms * (dmxdx + dmydy);

					//-------------------LLG Equation------------------------------//
					Hx_eff_glb = Hx_stat(i) + \
						Hx_anis_glb + Hx_elas_glb + \
						Hx_exch_glb + Hx_iDMI_glb + pt_glb->Hext[0] + pt_glb->Hextx_nonunif(i);

					Hy_eff_glb = Hy_stat(i) + \
						Hy_anis_glb + Hy_elas_glb + \
						Hy_exch_glb + Hy_iDMI_glb + pt_glb->Hext[1] + pt_glb->Hexty_nonunif(i);

					Hz_eff_glb = Hz_stat(i) + \
						Hz_anis_glb + Hz_elas_glb + \
						Hz_exch_glb + Hz_iDMI_glb + pt_glb->Hext[2] + pt_glb->Hextz_nonunif(i);

					if (pt_glb->if_EM_backaction_on_M == true) {
						Hx_eff_glb = Hx_eff_glb + pt_EM->DHx_em_cell(i);
						Hy_eff_glb = Hy_eff_glb + pt_EM->DHy_em_cell(i);
						Hz_eff_glb = Hz_eff_glb + pt_EM->DHz_em_cell(i);
					}

					alpha = mat->FM_damping;
					gyro = mat->gyro;
					precess = -gyro / (1 + alpha * alpha) * pt_glb->dt;

					precess_termx = my_glb_local * Hz_eff_glb - mz_glb_local * Hy_eff_glb;
					precess_termy = mz_glb_local * Hx_eff_glb - mx_glb_local * Hz_eff_glb;
					precess_termz = mx_glb_local * Hy_eff_glb - my_glb_local * Hx_eff_glb;

					damping_termx = my_glb_local * precess_termz - mz_glb_local * precess_termy;
					damping_termy = mz_glb_local * precess_termx - mx_glb_local * precess_termz;
					damping_termz = mx_glb_local * precess_termy - my_glb_local * precess_termx;

					dmx_glb_rk3(i) = precess * precess_termx + alpha * precess * damping_termx;
					dmy_glb_rk3(i) = precess * precess_termy + alpha * precess * damping_termy;
					dmz_glb_rk3(i) = precess * precess_termz + alpha * precess * damping_termz;
				}
			}
		}
	}
	if (pt_glb->if_spin_pumping == true) {
		get_J_ISHE_RK3();
	}
}

void magnetic_system::get_J_ISHE_RK3() {
	long i, j, sz;
	unsigned int mat_type;
	material* mat;
	double coeff_Js, coeff_L;
	double distance;

#pragma acc parallel default(present) async(4)
	{
#pragma acc loop gang vector private(i,j,mat_type, mat, coeff_Js,coeff_L, sz, distance)
		for (long k = 0; k < nz; k++) {

			if (pt_glb->if_J_ISHE[k] == true) {
				mat_type = pt_glb->material_cell(0, 0, k);
				mat = &(pt_glb->material_parameters[(mat_type)-1]);

				sz = pt_glb->source_z[k];

				coeff_Js = mat->spin_hall_angle * unit_e / 2. / PI * mat->spin_mix_cond;

				distance = (abs(static_cast<double>((k - sz))) - 0.5) * dz;

				if (k > sz) {
					coeff_L = sinh((pt_glb->thickness_pumping_layer[k] - distance) / mat->spin_diffus_length) \
						/ sinh(pt_glb->thickness_pumping_layer[k] / mat->spin_diffus_length);
				}
				else {
					coeff_L = -1. * sinh((pt_glb->thickness_pumping_layer[k] - distance) / mat->spin_diffus_length) \
						/ sinh(pt_glb->thickness_pumping_layer[k] / mat->spin_diffus_length);
				}

#pragma acc loop seq
				for (i = 0; i < nx; i++) {
#pragma acc loop seq
					for (j = 0; j < ny; j++) {
						Jx_ISHE(i, j, k) = coeff_Js * coeff_L * \
							(mx_glb_store(i, j, sz) * dmz_glb_rk3(i, j, sz) - mz_glb_store(i, j, sz) * dmx_glb_rk3(i, j, sz)\
								+ mx_AFM1_glb_store(i, j, sz) * dmz_AFM1_glb_rk3(i, j, sz) - mz_AFM1_glb_store(i, j, sz) * dmx_AFM1_glb_rk3(i, j, sz)\
								+ mx_AFM2_glb_store(i, j, sz) * dmz_AFM2_glb_rk3(i, j, sz) - mz_AFM2_glb_store(i, j, sz) * dmx_AFM2_glb_rk3(i, j, sz)) \
							/ pt_glb->dt;

						Jy_ISHE(i, j, k) = coeff_Js * coeff_L * \
							(my_glb_store(i, j, sz) * dmz_glb_rk3(i, j, sz) - mz_glb_store(i, j, sz) * dmy_glb_rk3(i, j, sz)\
								+ my_AFM1_glb_store(i, j, sz) * dmz_AFM1_glb_rk3(i, j, sz) - mz_AFM1_glb_store(i, j, sz) * dmy_AFM1_glb_rk3(i, j, sz)\
								+ my_AFM2_glb_store(i, j, sz) * dmz_AFM2_glb_rk3(i, j, sz) - mz_AFM2_glb_store(i, j, sz) * dmy_AFM2_glb_rk3(i, j, sz)) \
							/ pt_glb->dt;
					}
				}
			}

		}
	}
}

void magnetic_system::get_dm_RK4() {
	double Hx_anis_crt, Hy_anis_crt, Hz_anis_crt;
	double Hx_anis_glb, Hy_anis_glb, Hz_anis_glb;
	double Hx_elas_crt, Hy_elas_crt, Hz_elas_crt;
	double Hx_elas_glb, Hy_elas_glb, Hz_elas_glb;
	double Hx_exch_glb, Hy_exch_glb, Hz_exch_glb;
	double Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb;

	double Hx_eff_glb, Hy_eff_glb, Hz_eff_glb;
	double precess_termx, precess_termy, precess_termz;
	double damping_termx, damping_termy, damping_termz;

	long int x, y, z;
	unsigned int mat_type;
	material* mat;
	double Ms, K1, K2, B1, B2;
	double Aexi, iDMIi;
	double alpha, gyro, precess;
	double mx_glb_local, my_glb_local, mz_glb_local;
	double mx, my, mz;
	double exx, eyy, ezz, eyz, exz, exy;
	//double d2mxdx2, d2mxdy2, d2mxdz2;
	//double d2mydx2, d2mydy2, d2mydz2;
	//double d2mzdx2, d2mzdy2, d2mzdz2;
	double dmxdx, dmydy, dmzdx, dmzdy;
	double dmxdx2_f, dmydx2_f, dmzdx2_f;
	double dmxdx2_b, dmydx2_b, dmzdx2_b;
	double dmxdy2_f, dmydy2_f, dmzdy2_f;
	double dmxdy2_b, dmydy2_b, dmzdy2_b;
	double dmxdz2_f, dmydz2_f, dmzdz2_f;
	double dmxdz2_b, dmydz2_b, dmzdz2_b;

	if (pt_glb->if_prescribe_Hext == false) {
#pragma acc serial default(present) async(4)
		{
			update_external_Hfield_full();
		}
	}

#pragma acc parallel default(present) async(4)
	{
#pragma acc loop gang vector private(Hx_anis_crt, Hy_anis_crt, Hz_anis_crt,\
	Hx_anis_glb, Hy_anis_glb, Hz_anis_glb,\
	Hx_elas_crt, Hy_elas_crt, Hz_elas_crt,\
	Hx_elas_glb, Hy_elas_glb, Hz_elas_glb,\
	Hx_exch_glb, Hy_exch_glb, Hz_exch_glb,\
	Hx_iDMI_glb, Hy_iDMI_glb, Hz_iDMI_glb,\
	Hx_eff_glb, Hy_eff_glb, Hz_eff_glb,\
	precess_termx, precess_termy, precess_termz,\
	damping_termx, damping_termy, damping_termz,\
	x,y,z,\
	mat_type,\
	Ms, K1, K2, B1, B2, \
	Aexi, iDMIi, \
	alpha, gyro, precess,\
	mx_glb_local, my_glb_local, mz_glb_local,\
	mx, my, mz,\
	exx, eyy, ezz, eyz, exz, exy, mat,\
	dmxdx, dmydy, dmzdx, dmzdy, \
	dmxdx2_f, dmydx2_f, dmzdx2_f,\
	dmxdx2_b, dmydx2_b, dmzdx2_b,\
	dmxdy2_f, dmydy2_f, dmzdy2_f,\
	dmxdy2_b, dmydy2_b, dmzdy2_b,\
	dmxdz2_f, dmydz2_f, dmzdz2_f,\
	dmxdz2_b, dmydz2_b, dmzdz2_b)
		for (long int i = 0; i < n; i++) {
			mat_type = pt_glb->material_cell(i);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_FM == true) {
					x = i / (ny * nz);
					y = (i - x * (ny * nz)) / nz;
					z = i - x * (ny * nz) - y * nz;

					mx_glb_local = mx_glb_store(i);
					my_glb_local = my_glb_store(i);
					mz_glb_local = mz_glb_store(i);

					pt_math->transform_vector_glb2crt(mx_glb_local, my_glb_local, mz_glb_local, \
						mx, my, mz);

					Ms = mat->Ms;
					K1 = mat->K1 / mu0 / Ms; K2 = mat->K2 / mu0 / Ms;
					B1 = mat->B1 / mu0 / Ms; B2 = mat->B2 / mu0 / Ms;
					Aexi = mat->Aex; iDMIi = mat->iDMI;

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

					get_m_spatial_derivative_local \
						(x, y, z, \
							Aexi, iDMIi, \
							dmxdx, dmydy, dmzdx, dmzdy, \
							dmxdx2_f, dmydx2_f, dmzdx2_f, \
							dmxdx2_b, dmydx2_b, dmzdx2_b, \
							dmxdy2_f, dmydy2_f, dmzdy2_f, \
							dmxdy2_b, dmydy2_b, dmzdy2_b, \
							dmxdz2_f, dmydz2_f, dmzdz2_f, \
							dmxdz2_b, dmydz2_b, dmzdz2_b);

					//get_laplacian_m_local \
						//	(x, y, z, \
						//		d2mxdx2, d2mxdy2, d2mxdz2, \
						//		d2mydx2, d2mydy2, d2mydz2, \
						//		d2mzdx2, d2mzdy2, d2mzdz2);

						//------------Magnetocrystalline anisotropy field--------------//
					if (mat->anisotropy_type == 1) {
						Hx_anis_crt = -2. * \
							(K1 * (my * my + mz * mz) \
								+ K2 * my * my * mz * mz) * mx;
						Hy_anis_crt = -2. * \
							(K1 * (mx * mx + mz * mz) \
								+ K2 * mx * mx * mz * mz) * my;
						Hz_anis_crt = -2. * \
							(K1 * (mx * mx + my * my) \
								+ K2 * mx * mx * my * my) * mz;
					}
					else if (mat->anisotropy_type == 2) {
						Hx_anis_crt = 0.;
						Hy_anis_crt = 0.;
						Hz_anis_crt = 2. * K1 * mz;
					}
					pt_math->transform_vector_crt2glb(Hx_anis_crt, Hy_anis_crt, Hz_anis_crt, \
						Hx_anis_glb, Hy_anis_glb, Hz_anis_glb);

					//Hx_anis(i) = Hx_anis_glb;
					//Hy_anis(i) = Hy_anis_glb;
					//Hz_anis(i) = Hz_anis_glb;

					//------------------Magnetoelastic effective field---------------//
					Hx_elas_crt = -2. * (B1 * (exx)*mx + \
						B2 * ((exy)*my + (exz)*mz));
					Hy_elas_crt = -2. * (B1 * (eyy)*my + \
						B2 * ((exy)*mx + (eyz)*mz));
					Hz_elas_crt = -2. * (B1 * (ezz)*mz + \
						B2 * ((exz)*mx + (eyz)*my));

					pt_math->transform_vector_crt2glb(Hx_elas_crt, Hy_elas_crt, Hz_elas_crt, \
						Hx_elas_glb, Hy_elas_glb, Hz_elas_glb);

					//Hx_elas(i) = Hx_elas_glb;
					//Hy_elas(i) = Hy_elas_glb;
					//Hz_elas(i) = Hz_elas_glb;

					//-----------------Magnetic exchange field----------------------//
					//----------Averaged parameter for inter-region exchange - Approach 1---------//
					Hx_exch_glb = 2. / mu0 / Ms * \
						(Aij_xf(i) * dmxdx2_f + Aij_xb(i) * dmxdx2_b + Aij_yf(i) * dmxdy2_f + Aij_yb(i) * dmxdy2_b + Aij_zf(i) * dmxdz2_f + Aij_zb(i) * dmxdz2_b);
					Hy_exch_glb = 2. / mu0 / Ms * \
						(Aij_xf(i) * dmydx2_f + Aij_xb(i) * dmydx2_b + Aij_yf(i) * dmydy2_f + Aij_yb(i) * dmydy2_b + Aij_zf(i) * dmydz2_f + Aij_zb(i) * dmydz2_b);
					Hz_exch_glb = 2. / mu0 / Ms * \
						(Aij_xf(i) * dmzdx2_f + Aij_xb(i) * dmzdx2_b + Aij_yf(i) * dmzdy2_f + Aij_yb(i) * dmzdy2_b + Aij_zf(i) * dmzdz2_f + Aij_zb(i) * dmzdz2_b);
					//----------Averaged parameter for inter-region exchange - Approach 2---------//
					//Hx_exch_glb = 2. / mu0 * \
					//	AM_exchange(i) * (dmxdx2_f + dmxdx2_b + dmxdy2_f + dmxdy2_b + dmxdz2_f + dmxdz2_b);
					//Hy_exch_glb = 2. / mu0 * \
					//	AM_exchange(i) * (dmydx2_f + dmydx2_b + dmydy2_f + dmydy2_b + dmydz2_f + dmydz2_b);
					//Hz_exch_glb = 2. / mu0 * \
					//	AM_exchange(i) * (dmzdx2_f + dmzdx2_b + dmzdy2_f + dmzdy2_b + dmzdz2_f + dmzdz2_b);

					//Hx_exch(i) = 2. * Aex * (d2mxdx2 + d2mxdy2 + d2mxdz2);
					//Hy_exch(i) = 2. * Aex * (d2mydx2 + d2mydy2 + d2mydz2);
					//Hz_exch(i) = 2. * Aex * (d2mzdx2 + d2mzdy2 + d2mzdz2);

					//------------------Interfacial DMI effective field-------------------//
					Hx_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdx;
					Hy_iDMI_glb = 2. * iDMIi / mu0 / Ms * dmzdy;
					Hz_iDMI_glb = -2. * iDMIi / mu0 / Ms * (dmxdx + dmydy);

					//-------------------LLG Equation------------------------------//
					Hx_eff_glb = Hx_stat(i) + \
						Hx_anis_glb + Hx_elas_glb + \
						Hx_exch_glb + Hx_iDMI_glb + pt_glb->Hext[0] + pt_glb->Hextx_nonunif(i);

					Hy_eff_glb = Hy_stat(i) + \
						Hy_anis_glb + Hy_elas_glb + \
						Hy_exch_glb + Hy_iDMI_glb + pt_glb->Hext[1] + pt_glb->Hexty_nonunif(i);

					Hz_eff_glb = Hz_stat(i) + \
						Hz_anis_glb + Hz_elas_glb + \
						Hz_exch_glb + Hz_iDMI_glb + pt_glb->Hext[2] + pt_glb->Hextz_nonunif(i);

					if (pt_glb->if_EM_backaction_on_M == true) {
						Hx_eff_glb = Hx_eff_glb + pt_EM->DHx_em_cell(i);
						Hy_eff_glb = Hy_eff_glb + pt_EM->DHy_em_cell(i);
						Hz_eff_glb = Hz_eff_glb + pt_EM->DHz_em_cell(i);
					}

					alpha = mat->FM_damping;
					gyro = mat->gyro;
					precess = -gyro / (1 + alpha * alpha) * pt_glb->dt;

					precess_termx = my_glb_local * Hz_eff_glb - mz_glb_local * Hy_eff_glb;
					precess_termy = mz_glb_local * Hx_eff_glb - mx_glb_local * Hz_eff_glb;
					precess_termz = mx_glb_local * Hy_eff_glb - my_glb_local * Hx_eff_glb;

					damping_termx = my_glb_local * precess_termz - mz_glb_local * precess_termy;
					damping_termy = mz_glb_local * precess_termx - mx_glb_local * precess_termz;
					damping_termz = mx_glb_local * precess_termy - my_glb_local * precess_termx;

					dmx_glb_rk4(i) = precess * precess_termx + alpha * precess * damping_termx;
					dmy_glb_rk4(i) = precess * precess_termy + alpha * precess * damping_termy;
					dmz_glb_rk4(i) = precess * precess_termz + alpha * precess * damping_termz;
				}
			}
		}
	}

	if (pt_glb->if_spin_pumping == true) {
		get_J_ISHE_RK4();
	}
}

void magnetic_system::get_J_ISHE_RK4() {
	long i, j, sz;
	unsigned int mat_type;
	material* mat;
	double coeff_Js, coeff_L;
	double distance;

#pragma acc parallel default(present) async(4)
	{
#pragma acc loop gang vector private(i,j,mat_type, mat, coeff_Js,coeff_L, sz, distance)
		for (long k = 0; k < nz; k++) {

			if (pt_glb->if_J_ISHE[k] == true) {
				mat_type = pt_glb->material_cell(0, 0, k);
				mat = &(pt_glb->material_parameters[(mat_type)-1]);

				sz = pt_glb->source_z[k];

				coeff_Js = mat->spin_hall_angle * unit_e / 2. / PI * mat->spin_mix_cond;

				distance = (abs(static_cast<double>((k - sz))) - 0.5) * dz;

				if (k > sz) {
					coeff_L = sinh((pt_glb->thickness_pumping_layer[k] - distance) / mat->spin_diffus_length) \
						/ sinh(pt_glb->thickness_pumping_layer[k] / mat->spin_diffus_length);
				}
				else {
					coeff_L = -1. * sinh((pt_glb->thickness_pumping_layer[k] - distance) / mat->spin_diffus_length) \
						/ sinh(pt_glb->thickness_pumping_layer[k] / mat->spin_diffus_length);
				}

#pragma acc loop seq
				for (i = 0; i < nx; i++) {
#pragma acc loop seq
					for (j = 0; j < ny; j++) {
						Jx_ISHE(i, j, k) = coeff_Js * coeff_L * \
							(mx_glb_store(i, j, sz) * dmz_glb_rk4(i, j, sz) - mz_glb_store(i, j, sz) * dmx_glb_rk4(i, j, sz)\
								+ mx_AFM1_glb_store(i, j, sz) * dmz_AFM1_glb_rk4(i, j, sz) - mz_AFM1_glb_store(i, j, sz) * dmx_AFM1_glb_rk4(i, j, sz)\
								+ mx_AFM2_glb_store(i, j, sz) * dmz_AFM2_glb_rk4(i, j, sz) - mz_AFM2_glb_store(i, j, sz) * dmx_AFM2_glb_rk4(i, j, sz)) \
							/ pt_glb->dt;

						Jy_ISHE(i, j, k) = coeff_Js * coeff_L * \
							(my_glb_store(i, j, sz) * dmz_glb_rk4(i, j, sz) - mz_glb_store(i, j, sz) * dmy_glb_rk4(i, j, sz)\
								+ my_AFM1_glb_store(i, j, sz) * dmz_AFM1_glb_rk4(i, j, sz) - mz_AFM1_glb_store(i, j, sz) * dmy_AFM1_glb_rk4(i, j, sz)\
								+ my_AFM2_glb_store(i, j, sz) * dmz_AFM2_glb_rk4(i, j, sz) - mz_AFM2_glb_store(i, j, sz) * dmy_AFM2_glb_rk4(i, j, sz)) \
							/ pt_glb->dt;
					}
				}
			}

		}
	}
}

void magnetic_system::get_dm() {
	//	double Hx_anis_crt, Hy_anis_crt, Hz_anis_crt;
	//	double Hx_anis_glb, Hy_anis_glb, Hz_anis_glb;
	//	double Hx_elas_crt, Hy_elas_crt, Hz_elas_crt;
	//	double Hx_elas_glb, Hy_elas_glb, Hz_elas_glb;
	//	double Hx_exch_glb, Hy_exch_glb, Hz_exch_glb;
	//
	//	double Hx_eff_glb, Hy_eff_glb, Hz_eff_glb;
	//	double precess_termx, precess_termy, precess_termz;
	//	double damping_termx, damping_termy, damping_termz;
	//
	//	long int x, y, z;
	//	unsigned int mat_type, mat_typen;
	//	material* mat, * matn;
	//	double Ms, K1, K2, B1, B2;
	//	double Aex1, Aex2;
	//	double Ms1, Ms2;
	//	double AM_xf, AM_yf, AM_zf;
	//	double AM_xb, AM_yb, AM_zb;
	//	double alpha, gyro, precess;
	//	double mx_glb_local, my_glb_local, mz_glb_local;
	//	double mx, my, mz;
	//	double exx, eyy, ezz, eyz, exz, exy;
	//	//double d2mxdx2, d2mxdy2, d2mxdz2;
	//	//double d2mydx2, d2mydy2, d2mydz2;
	//	//double d2mzdx2, d2mzdy2, d2mzdz2;
	//	double dmxdx2_f, dmydx2_f, dmzdx2_f;
	//	double dmxdx2_b, dmydx2_b, dmzdx2_b;
	//	double dmxdy2_f, dmydy2_f, dmzdy2_f;
	//	double dmxdy2_b, dmydy2_b, dmzdy2_b;
	//	double dmxdz2_f, dmydz2_f, dmzdz2_f;
	//	double dmxdz2_b, dmydz2_b, dmzdz2_b;
	//
	//	//#pragma acc declare device_resident(Hx_anis_crt, Hy_anis_crt, Hz_anis_crt,\
	//	//	Hx_anis_glb, Hy_anis_glb, Hz_anis_glb,\
	//	//	Hx_elas_crt, Hy_elas_crt, Hz_elas_crt,\
	//	//	Hx_elas_glb, Hy_elas_glb, Hz_elas_glb,\
	//	//	Hx_eff_glb, Hy_eff_glb, Hz_eff_glb,\
	//	//	precess_termx, precess_termy, precess_termz,\
	//	//	damping_termx, damping_termy, damping_termz,\
	//	//	mat_type,\
	//	//	Ms, K1, K2, B1, B2, Aex,\
	//	//	alpha, gyro, precess,\
	//	//	mx_glb_local, my_glb_local, mz_glb_local,\
	//	//	mx, my, mz,\
	//	//	exx, eyy, ezz, eyz, exz, exy,\
	//	//	d2mxdx2, d2mxdy2, d2mxdz2,\
	//	//	d2mydx2, d2mydy2, d2mydz2,\
	//	//	d2mzdx2, d2mzdy2, d2mzdz2) deviceptr(mat)
	//#pragma acc parallel default(present)
	//	{
	//#pragma acc loop gang vector private(Hx_anis_crt, Hy_anis_crt, Hz_anis_crt,\
	//	Hx_anis_glb, Hy_anis_glb, Hz_anis_glb,\
	//	Hx_elas_crt, Hy_elas_crt, Hz_elas_crt,\
	//	Hx_elas_glb, Hy_elas_glb, Hz_elas_glb,\
	//	Hx_exch_glb, Hy_exch_glb, Hz_exch_glb,\
	//	Hx_eff_glb, Hy_eff_glb, Hz_eff_glb,\
	//	precess_termx, precess_termy, precess_termz,\
	//	damping_termx, damping_termy, damping_termz,\
	//	x,y,z,\
	//	mat_type,mat_typen,\
	//	Ms, K1, K2, B1, B2, \
	//	Aex1, Aex2,\
	//	Ms1, Ms2,\
	//	AM_xf, AM_yf, AM_zf,\
	//	AM_xb, AM_yb, AM_zb,\
	//	alpha, gyro, precess,\
	//	mx_glb_local, my_glb_local, mz_glb_local,\
	//	mx, my, mz,\
	//	exx, eyy, ezz, eyz, exz, exy, mat, matn,\
	//	dmxdx2_f, dmydx2_f, dmzdx2_f,\
	//	dmxdx2_b, dmydx2_b, dmzdx2_b,\
	//	dmxdy2_f, dmydy2_f, dmzdy2_f,\
	//	dmxdy2_b, dmydy2_b, dmzdy2_b,\
	//	dmxdz2_f, dmydz2_f, dmzdz2_f,\
	//	dmxdz2_b, dmydz2_b, dmzdz2_b)
	//		for (long int i = 0; i < n; i++) {
	//			mat_type = pt_glb->material_cell(i);
	//			if (mat_type != 0) {
	//				mat = &(pt_glb->material_parameters[(mat_type)-1]);
	//				if (mat->if_FM == true) {
	//					x = i / (ny * nz);
	//					y = (i - x * (ny * nz)) / nz;
	//					z = i - x * (ny * nz) - y * nz;
	//
	//					mx_glb_local = mx_glb(i);
	//					my_glb_local = my_glb(i);
	//					mz_glb_local = mz_glb(i);
	//
	//					pt_math->transform_vector_glb2crt(mx_glb_local, my_glb_local, mz_glb_local, \
	//						mx, my, mz);
	//
	//					Ms = mat->Ms;
	//					K1 = mat->K1 / mu0 / Ms; K2 = mat->K2 / mu0 / Ms;
	//					B1 = mat->B1 / mu0 / Ms; B2 = mat->B2 / mu0 / Ms;
	//
	//					if (FM_surfYZ(x + 1, y, z) == true) {
	//						AM_xf = 0.;
	//					}
	//					else {
	//						if (x == nx - 1) {
	//							mat_typen = pt_glb->material_cell(0, y, z);
	//							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
	//							Aex1 = mat->Aex; Aex2 = matn->Aex;
	//							Ms1 = mat->Ms; Ms2 = matn->Ms;
	//							AM_xf = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
	//						}
	//						else {
	//							mat_typen = pt_glb->material_cell(x + 1, y, z);
	//							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
	//							Aex1 = mat->Aex; Aex2 = matn->Aex;
	//							Ms1 = mat->Ms; Ms2 = matn->Ms;
	//							AM_xf = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
	//						}
	//					}
	//
	//					if (FM_surfYZ(x, y, z) == true) {
	//						AM_xb = 0.;
	//					}
	//					else {
	//						if (x == 0) {
	//							mat_typen = pt_glb->material_cell(nx - 1, y, z);
	//							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
	//							Aex1 = mat->Aex; Aex2 = matn->Aex;
	//							Ms1 = mat->Ms; Ms2 = matn->Ms;
	//							AM_xb = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
	//						}
	//						else {
	//							mat_typen = pt_glb->material_cell(x - 1, y, z);
	//							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
	//							Aex1 = mat->Aex; Aex2 = matn->Aex;
	//							Ms1 = mat->Ms; Ms2 = matn->Ms;
	//							AM_xb = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
	//						}
	//					}
	//
	//					if (FM_surfXZ(x, y + 1, z) == true) {
	//						AM_yf = 0.;
	//					}
	//					else {
	//						if (y == ny - 1) {
	//							mat_typen = pt_glb->material_cell(x, 0, z);
	//							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
	//							Aex1 = mat->Aex; Aex2 = matn->Aex;
	//							Ms1 = mat->Ms; Ms2 = matn->Ms;
	//							AM_yf = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
	//						}
	//						else {
	//							mat_typen = pt_glb->material_cell(x, y + 1, z);
	//							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
	//							Aex1 = mat->Aex; Aex2 = matn->Aex;
	//							Ms1 = mat->Ms; Ms2 = matn->Ms;
	//							AM_yf = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
	//						}
	//					}
	//
	//					if (FM_surfXZ(x, y, z) == true) {
	//						AM_yb = 0.;
	//					}
	//					else {
	//						if (y == 0) {
	//							mat_typen = pt_glb->material_cell(x, ny - 1, z);
	//							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
	//							Aex1 = mat->Aex; Aex2 = matn->Aex;
	//							Ms1 = mat->Ms; Ms2 = matn->Ms;
	//							AM_yb = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
	//						}
	//						else {
	//							mat_typen = pt_glb->material_cell(x, y - 1, z);
	//							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
	//							Aex1 = mat->Aex; Aex2 = matn->Aex;
	//							Ms1 = mat->Ms; Ms2 = matn->Ms;
	//							AM_yb = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
	//						}
	//					}
	//
	//					if (FM_surfXY(x, y, z + 1) == true) {
	//						AM_zf = 0.;
	//					}
	//					else {
	//						if (z == nz - 1) {
	//							mat_typen = pt_glb->material_cell(x, y, 0);
	//							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
	//							Aex1 = mat->Aex; Aex2 = matn->Aex;
	//							Ms1 = mat->Ms; Ms2 = matn->Ms;
	//							AM_zf = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
	//						}
	//						else {
	//							mat_typen = pt_glb->material_cell(x, y, z + 1);
	//							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
	//							Aex1 = mat->Aex; Aex2 = matn->Aex;
	//							Ms1 = mat->Ms; Ms2 = matn->Ms;
	//							AM_zf = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
	//						}
	//					}
	//
	//					if (FM_surfXY(x, y, z) == true) {
	//						AM_zb = 0.;
	//					}
	//					else {
	//						if (z == 0) {
	//							mat_typen = pt_glb->material_cell(x, y, nz - 1);
	//							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
	//							Aex1 = mat->Aex; Aex2 = matn->Aex;
	//							Ms1 = mat->Ms; Ms2 = matn->Ms;
	//							AM_zb = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
	//						}
	//						else {
	//							mat_typen = pt_glb->material_cell(x, y, z - 1);
	//							matn = &(pt_glb->material_parameters[(mat_typen)-1]);
	//							Aex1 = mat->Aex; Aex2 = matn->Aex;
	//							Ms1 = mat->Ms; Ms2 = matn->Ms;
	//							AM_zb = 2. * (Aex1 / Ms1) * (Aex2 / Ms2) / (Aex1 / Ms1 + Aex2 / Ms2);
	//						}
	//					}
	//
	//					exx = pt_elas->exxt0_crt(i) + pt_elas->Dexx_crt(i);
	//					eyy = pt_elas->eyyt0_crt(i) + pt_elas->Deyy_crt(i);
	//					ezz = pt_elas->ezzt0_crt(i) + pt_elas->Dezz_crt(i);
	//					eyz = pt_elas->eyzt0_crt(i) + pt_elas->Deyz_crt(i);
	//					exz = pt_elas->exzt0_crt(i) + pt_elas->Dexz_crt(i);
	//					exy = pt_elas->exyt0_crt(i) + pt_elas->Dexy_crt(i);
	//
	//					get_laplacian_m_local \
	//						(x, y, z, \
	//							dmxdx2_f, dmydx2_f, dmzdx2_f, \
	//							dmxdx2_b, dmydx2_b, dmzdx2_b, \
	//							dmxdy2_f, dmydy2_f, dmzdy2_f, \
	//							dmxdy2_b, dmydy2_b, dmzdy2_b, \
	//							dmxdz2_f, dmydz2_f, dmzdz2_f, \
	//							dmxdz2_b, dmydz2_b, dmzdz2_b);
	//
	//					//get_laplacian_m_local \
	//						//	(x, y, z, \
	//						//		d2mxdx2, d2mxdy2, d2mxdz2, \
	//						//		d2mydx2, d2mydy2, d2mydz2, \
	//						//		d2mzdx2, d2mzdy2, d2mzdz2);
	//
	//						//------------Magnetocrystalline anisotropy field--------------//
	//					if (mat->anisotropy_type == 1) {
	//						Hx_anis_crt = -2. * \
	//							(K1 * (my * my + mz * mz) \
	//								+ K2 * my * my * mz * mz) * mx;
	//						Hy_anis_crt = -2. * \
	//							(K1 * (mx * mx + mz * mz) \
	//								+ K2 * mx * mx * mz * mz) * my;
	//						Hz_anis_crt = -2. * \
	//							(K1 * (mx * mx + my * my) \
	//								+ K2 * mx * mx * my * my) * mz;
	//					}
	//					else if (mat->anisotropy_type == 2) {
	//						Hx_anis_crt = 0.;
	//						Hy_anis_crt = 0.;
	//						Hz_anis_crt = 2. * K1 * mz;
	//					}
	//					pt_math->transform_vector_crt2glb(Hx_anis_crt, Hy_anis_crt, Hz_anis_crt, \
	//						Hx_anis_glb, Hy_anis_glb, Hz_anis_glb);
	//
	//					//Hx_anis(i) = Hx_anis_glb;
	//					//Hy_anis(i) = Hy_anis_glb;
	//					//Hz_anis(i) = Hz_anis_glb;
	//
	//					//------------------Magnetoelastic effective field---------------//
	//					Hx_elas_crt = -2. * (B1 * (exx)*mx + \
	//						B2 * ((exy)*my + (exz)*mz));
	//					Hy_elas_crt = -2. * (B1 * (eyy)*my + \
	//						B2 * ((exy)*mx + (eyz)*mz));
	//					Hz_elas_crt = -2. * (B1 * (ezz)*mz + \
	//						B2 * ((exz)*mx + (eyz)*my));
	//
	//					pt_math->transform_vector_crt2glb(Hx_elas_crt, Hy_elas_crt, Hz_elas_crt, \
	//						Hx_elas_glb, Hy_elas_glb, Hz_elas_glb);
	//
	//					//Hx_elas(i) = Hx_elas_glb;
	//					//Hy_elas(i) = Hy_elas_glb;
	//					//Hz_elas(i) = Hz_elas_glb;
	//
	//					//-----------------Magnetic exchange field----------------------//
	//					Hx_exch_glb = 2. / mu0 * (AM_xf * dmxdx2_f + AM_xb * dmxdx2_b + AM_yf * dmxdy2_f + AM_yb * dmxdy2_b + AM_zf * dmxdz2_f + AM_zb * dmxdz2_b);
	//					Hy_exch_glb = 2. / mu0 * (AM_xf * dmydx2_f + AM_xb * dmydx2_b + AM_yf * dmydy2_f + AM_yb * dmydy2_b + AM_zf * dmydz2_f + AM_zb * dmydz2_b);
	//					Hz_exch_glb = 2. / mu0 * (AM_xf * dmzdx2_f + AM_xb * dmzdx2_b + AM_yf * dmzdy2_f + AM_yb * dmzdy2_b + AM_zf * dmzdz2_f + AM_zb * dmzdz2_b);
	//					//Hx_exch(i) = 2. * Aex * (d2mxdx2 + d2mxdy2 + d2mxdz2);
	//					//Hy_exch(i) = 2. * Aex * (d2mydx2 + d2mydy2 + d2mydz2);
	//					//Hz_exch(i) = 2. * Aex * (d2mzdx2 + d2mzdy2 + d2mzdz2);
	//
	//					//-------------------LLG Equation------------------------------//
	//					Hx_eff_glb = Hx_stat(i) + \
	//						Hx_anis_glb + Hx_elas_glb + \
	//						Hx_exch_glb + pt_EM->DHx_em_cell(i) + pt_glb->Hext[0];
	//
	//					Hy_eff_glb = Hy_stat(i) + \
	//						Hy_anis_glb + Hy_elas_glb + \
	//						Hy_exch_glb + pt_EM->DHy_em_cell(i) + pt_glb->Hext[1];
	//
	//					Hz_eff_glb = Hz_stat(i) + \
	//						Hz_anis_glb + Hz_elas_glb + \
	//						Hz_exch_glb + pt_EM->DHz_em_cell(i) + pt_glb->Hext[2];
	//
	//					alpha = mat->FM_damping;
	//					gyro = mat->gyro;
	//					precess = -gyro / (1 + alpha * alpha) * pt_glb->dt;
	//
	//					precess_termx = my_glb_local * Hz_eff_glb - mz_glb_local * Hy_eff_glb;
	//					precess_termy = mz_glb_local * Hx_eff_glb - mx_glb_local * Hz_eff_glb;
	//					precess_termz = mx_glb_local * Hy_eff_glb - my_glb_local * Hx_eff_glb;
	//
	//					damping_termx = my_glb_local * precess_termz - mz_glb_local * precess_termy;
	//					damping_termy = mz_glb_local * precess_termx - mx_glb_local * precess_termz;
	//					damping_termz = mx_glb_local * precess_termy - my_glb_local * precess_termx;
	//
	//					dmx_glb(i) = precess * precess_termx + alpha * precess * damping_termx;
	//					dmy_glb(i) = precess * precess_termy + alpha * precess * damping_termy;
	//					dmz_glb(i) = precess * precess_termz + alpha * precess * damping_termz;
	//				}
	//			}
	//		}
	//	}

#pragma acc parallel default(present) async(4)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			dmx_glb(id) = dmx_glb_rk1(id) / 6. + dmx_glb_rk2(id) / 3. + dmx_glb_rk3(id) / 3. + dmx_glb_rk4(id) / 6.;
			dmy_glb(id) = dmy_glb_rk1(id) / 6. + dmy_glb_rk2(id) / 3. + dmy_glb_rk3(id) / 3. + dmy_glb_rk4(id) / 6.;
			dmz_glb(id) = dmz_glb_rk1(id) / 6. + dmz_glb_rk2(id) / 3. + dmz_glb_rk3(id) / 3. + dmz_glb_rk4(id) / 6.;
		}
	}
}

void magnetic_system::get_J_ISHE() {
	long i, j, sz;
	unsigned int mat_type;
	material* mat;
	double coeff_Js, coeff_L;
	double distance;

#pragma acc parallel default(present) async(3)
	{
#pragma acc loop gang vector private(i,j,mat_type, mat, coeff_Js,coeff_L, sz, distance)
		for (long k = 0; k < nz; k++) {

			if (pt_glb->if_J_ISHE[k] == true) {
				mat_type = pt_glb->material_cell(0, 0, k);
				mat = &(pt_glb->material_parameters[(mat_type)-1]);

				sz = pt_glb->source_z[k];

				coeff_Js = mat->spin_hall_angle * unit_e / 2. / PI * mat->spin_mix_cond;

				distance = (abs(static_cast<double>((k - sz))) - 0.5) * dz;

				if (k > sz) {
					coeff_L = sinh((pt_glb->thickness_pumping_layer[k] - distance) / mat->spin_diffus_length) \
						/ sinh(pt_glb->thickness_pumping_layer[k] / mat->spin_diffus_length);
				}
				else {
					coeff_L = -1. * sinh((pt_glb->thickness_pumping_layer[k] - distance) / mat->spin_diffus_length) \
						/ sinh(pt_glb->thickness_pumping_layer[k] / mat->spin_diffus_length);
				}

#pragma acc loop seq
				for (i = 0; i < nx; i++) {
#pragma acc loop seq
					for (j = 0; j < ny; j++) {
						Jx_ISHE(i, j, k) = coeff_Js * coeff_L * \
							(mx_glb(i, j, sz) * dmz_glb(i, j, sz) - mz_glb(i, j, sz) * dmx_glb(i, j, sz)\
								+ mx_AFM1_glb(i, j, sz) * dmz_AFM1_glb(i, j, sz) - mz_AFM1_glb(i, j, sz) * dmx_AFM1_glb(i, j, sz)\
								+ mx_AFM2_glb(i, j, sz) * dmz_AFM2_glb(i, j, sz) - mz_AFM2_glb(i, j, sz) * dmx_AFM2_glb(i, j, sz)) \
							/ pt_glb->dt;

						Jy_ISHE(i, j, k) = coeff_Js * coeff_L * \
							(my_glb(i, j, sz) * dmz_glb(i, j, sz) - mz_glb(i, j, sz) * dmy_glb(i, j, sz)\
								+ my_AFM1_glb(i, j, sz) * dmz_AFM1_glb(i, j, sz) - mz_AFM1_glb(i, j, sz) * dmy_AFM1_glb(i, j, sz)\
								+ my_AFM2_glb(i, j, sz) * dmz_AFM2_glb(i, j, sz) - mz_AFM2_glb(i, j, sz) * dmy_AFM2_glb(i, j, sz)) \
							/ pt_glb->dt;
					}
				}
			}

		}
	}
}


void magnetic_system::update_m_RK1() {
#pragma acc parallel default(present) async(3)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			mx_glb_store(id) = mx_glb(id) + dmx_glb_rk1(id) * 0.5;
			my_glb_store(id) = my_glb(id) + dmy_glb_rk1(id) * 0.5;
			mz_glb_store(id) = mz_glb(id) + dmz_glb_rk1(id) * 0.5;
		}
	}

	if (pt_glb->if_prescribed_m == true) {
#pragma acc serial default(present) async(3)
		{
			update_prescribed_m_half();
		}
	}

	if (pt_glb->if_EMdynamic == false && pt_glb->if_mag_1Dmodel == true && pt_glb->if_magnetostatics == true) {
#pragma acc parallel default(present) async(3)
		{
			get_H_static_1D();
		}
	}
}

void magnetic_system::update_m_RK2() {
#pragma acc parallel default(present) async(3)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			mx_glb_store(id) = mx_glb(id) + dmx_glb_rk2(id) * 0.5;
			my_glb_store(id) = my_glb(id) + dmy_glb_rk2(id) * 0.5;
			mz_glb_store(id) = mz_glb(id) + dmz_glb_rk2(id) * 0.5;
		}
	}

	if (pt_glb->if_prescribed_m == true) {
#pragma acc serial default(present) async(3)
		{
			update_prescribed_m_half();
		}
	}

	if (pt_glb->if_EMdynamic == false && pt_glb->if_mag_1Dmodel == true && pt_glb->if_magnetostatics == true) {
#pragma acc parallel default(present) async(3)
		{
			get_H_static_1D();
		}
	}
}

void magnetic_system::update_m_RK3() {
#pragma acc parallel default(present) async(3)
	{
#pragma acc loop gang vector
		for (long id = 0; id < n; id++) {
			mx_glb_store(id) = mx_glb(id) + dmx_glb_rk3(id);
			my_glb_store(id) = my_glb(id) + dmy_glb_rk3(id);
			mz_glb_store(id) = mz_glb(id) + dmz_glb_rk3(id);
		}
	}

	if (pt_glb->if_prescribed_m == true) {
#pragma acc serial default(present) async(3)
		{
			update_prescribed_m_full();
		}
	}

	if (pt_glb->if_EMdynamic == false && pt_glb->if_mag_1Dmodel == true && pt_glb->if_magnetostatics == true) {
#pragma acc parallel default(present) async(3)
		{
			get_H_static_1D();
		}
	}
}

void magnetic_system::update_m() {

	if (pt_glb->if_spin_pumping == true) {
		get_J_ISHE();
	}

#pragma acc parallel default(present) async(3)
	{
		mx_glb += dmx_glb;
		my_glb += dmy_glb;
		mz_glb += dmz_glb;
	}

	if (pt_glb->if_prescribed_m == true) {
#pragma acc serial default(present) async(3)
		{
			update_prescribed_m();
		}
	}

#pragma acc parallel default(present) async(3)
	{
		mx_glb_store = mx_glb;
		my_glb_store = my_glb;
		mz_glb_store = mz_glb;
	}

	if (pt_glb->if_EMdynamic == false && pt_glb->if_mag_1Dmodel == true && pt_glb->if_magnetostatics == true) {
#pragma acc parallel default(present) async(3)
		{
			get_H_static_1D();
		}
	}

	//for (long int i = 0; i < nx; i++) {
	//	for (long int j = 0; j < ny; j++) {
	//		for (long int k = 0; k < nz; k++) {
	//			mx_glb.pt_matrix[i][j][k] = mx_glb.pt_matrix[i][j][k] + dmx_glb.pt_matrix[i][j][k];
	//			my_glb.pt_matrix[i][j][k] = my_glb.pt_matrix[i][j][k] + dmy_glb.pt_matrix[i][j][k];
	//			mz_glb.pt_matrix[i][j][k] = mz_glb.pt_matrix[i][j][k] + dmz_glb.pt_matrix[i][j][k];
	//		}
	//	}
	//}
}

#pragma acc routine seq// nohost
void magnetic_system::get_m_spatial_derivative_local \
(long int& i, long int& j, long int& k, \
	double& A, double& D, \
	double& dmxdx, double& dmydy, double& dmzdx, double& dmzdy, \
	double& dmxdx2_f, double& dmydx2_f, double& dmzdx2_f, \
	double& dmxdx2_b, double& dmydx2_b, double& dmzdx2_b, \
	double& dmxdy2_f, double& dmydy2_f, double& dmzdy2_f, \
	double& dmxdy2_b, double& dmydy2_b, double& dmzdy2_b, \
	double& dmxdz2_f, double& dmydz2_f, double& dmzdz2_f, \
	double& dmxdz2_b, double& dmydz2_b, double& dmzdz2_b)
{
	double mx_fwd, mx_center, mx_bwd;
	double my_fwd, my_center, my_bwd;
	double mz_fwd, mz_center, mz_bwd;

	double dmxdx_fwd, dmydx_fwd, dmzdx_fwd;
	double dmxdx_bwd, dmydx_bwd, dmzdx_bwd;
	double dmxdy_fwd, dmydy_fwd, dmzdy_fwd;
	double dmxdy_bwd, dmydy_bwd, dmzdy_bwd;
	double dmxdz_fwd, dmydz_fwd, dmzdz_fwd;
	double dmxdz_bwd, dmydz_bwd, dmzdz_bwd;

	mx_center = mx_glb_store(i, j, k);
	my_center = my_glb_store(i, j, k);
	mz_center = mz_glb_store(i, j, k);

	//--------------X direction (YZ surface)------------//
	//------------- negative X direction---------------//
	if (FM_surfYZ(i, j, k) == true) {
		//mx_bwd = mx_center; my_bwd = my_center; mz_bwd = mz_center;
		dmxdx_bwd = -D / 2. / A * mz_center;
		dmydx_bwd = 0.;
		dmzdx_bwd = D / 2. / A * mx_center;
	}
	else if (FM_surfYZ(i, j, k) == false) {
		if (i == 0) {
			mx_bwd = mx_glb_store(nx - 1, j, k);
			my_bwd = my_glb_store(nx - 1, j, k);
			mz_bwd = mz_glb_store(nx - 1, j, k);
		}
		else {
			mx_bwd = mx_glb_store(i - 1, j, k);
			my_bwd = my_glb_store(i - 1, j, k);
			mz_bwd = mz_glb_store(i - 1, j, k);
		}

		dmxdx_bwd = (mx_center - mx_bwd) / dx;
		dmydx_bwd = (my_center - my_bwd) / dx;
		dmzdx_bwd = (mz_center - mz_bwd) / dx;
	}
	//------------- positive X direction---------------//
	if (FM_surfYZ(i + 1, j, k) == true) {
		//mx_fwd = mx_center; my_fwd = my_center; mz_fwd = mz_center;
		dmxdx_fwd = -D / 2. / A * mz_center;
		dmydx_fwd = 0.;
		dmzdx_fwd = D / 2. / A * mx_center;
	}
	else if (FM_surfYZ(i + 1, j, k) == false) {
		if (i == nx - 1) {
			mx_fwd = mx_glb_store(0, j, k);
			my_fwd = my_glb_store(0, j, k);
			mz_fwd = mz_glb_store(0, j, k);
		}
		else {
			mx_fwd = mx_glb_store(i + 1, j, k);
			my_fwd = my_glb_store(i + 1, j, k);
			mz_fwd = mz_glb_store(i + 1, j, k);
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
	if (FM_surfXZ(i, j, k) == true) {
		//mx_bwd = mx_center; my_bwd = my_center; mz_bwd = mz_center;
		dmxdy_bwd = 0.;
		dmydy_bwd = -D / 2. / A * mz_center;
		dmzdy_bwd = D / 2. / A * my_center;
	}
	else if (FM_surfXZ(i, j, k) == false) {
		if (j == 0) {
			mx_bwd = mx_glb_store(i, ny - 1, k);
			my_bwd = my_glb_store(i, ny - 1, k);
			mz_bwd = mz_glb_store(i, ny - 1, k);
		}
		else {
			mx_bwd = mx_glb_store(i, j - 1, k);
			my_bwd = my_glb_store(i, j - 1, k);
			mz_bwd = mz_glb_store(i, j - 1, k);
		}

		dmxdy_bwd = (mx_center - mx_bwd) / dy;
		dmydy_bwd = (my_center - my_bwd) / dy;
		dmzdy_bwd = (mz_center - mz_bwd) / dy;
	}
	//------------- positive Y direction---------------//
	if (FM_surfXZ(i, j + 1, k) == true) {
		//mx_fwd = mx_center; my_fwd = my_center; mz_fwd = mz_center;
		dmxdy_fwd = 0.;
		dmydy_fwd = -D / 2. / A * mz_center;
		dmzdy_fwd = D / 2. / A * my_center;
	}
	else if (FM_surfXZ(i, j + 1, k) == false) {
		if (j == ny - 1) {
			mx_fwd = mx_glb_store(i, 0, k);
			my_fwd = my_glb_store(i, 0, k);
			mz_fwd = mz_glb_store(i, 0, k);
		}
		else {
			mx_fwd = mx_glb_store(i, j + 1, k);
			my_fwd = my_glb_store(i, j + 1, k);
			mz_fwd = mz_glb_store(i, j + 1, k);
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
	if (FM_surfXY(i, j, k) == true) {
		//mx_bwd = mx_center; my_bwd = my_center; mz_bwd = mz_center;
		dmxdz_bwd = 0.; dmydz_bwd = 0.; dmzdz_bwd = 0.;
	}
	else if (FM_surfXY(i, j, k) == false) {
		if (k == 0) {
			mx_bwd = mx_glb_store(i, j, nz - 1);
			my_bwd = my_glb_store(i, j, nz - 1);
			mz_bwd = mz_glb_store(i, j, nz - 1);
		}
		else {
			mx_bwd = mx_glb_store(i, j, k - 1);
			my_bwd = my_glb_store(i, j, k - 1);
			mz_bwd = mz_glb_store(i, j, k - 1);
		}

		dmxdz_bwd = (mx_center - mx_bwd) / dz;
		dmydz_bwd = (my_center - my_bwd) / dz;
		dmzdz_bwd = (mz_center - mz_bwd) / dz;
	}
	//------------- positive Z direction---------------//
	if (FM_surfXY(i, j, k + 1) == true) {
		//mx_fwd = mx_center; my_fwd = my_center; mz_fwd = mz_center;
		dmxdz_fwd = 0.; dmydz_fwd = 0.; dmzdz_fwd = 0.;
	}
	else if (FM_surfXY(i, j, k + 1) == false) {
		if (k == nz - 1) {
			mx_fwd = mx_glb_store(i, j, 0);
			my_fwd = my_glb_store(i, j, 0);
			mz_fwd = mz_glb_store(i, j, 0);
		}
		else {
			mx_fwd = mx_glb_store(i, j, k + 1);
			my_fwd = my_glb_store(i, j, k + 1);
			mz_fwd = mz_glb_store(i, j, k + 1);
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

#pragma acc routine seq
void magnetic_system::update_external_Hfield() {
	double local_time;

	local_time = pt_glb->time_device - pt_glb->dt;
	pt_glb->Hext[0] = pt_glb->Hext_stat[0] + pt_glb->Hext_altr[0] * sin(2. * PI * pt_glb->H_altr_freq[0] * local_time);
	pt_glb->Hext[1] = pt_glb->Hext_stat[1] + pt_glb->Hext_altr[1] * sin(2. * PI * pt_glb->H_altr_freq[1] * local_time);
	pt_glb->Hext[2] = pt_glb->Hext_stat[2] + pt_glb->Hext_altr[2] * sin(2. * PI * pt_glb->H_altr_freq[2] * local_time);
}

#pragma acc routine seq
void magnetic_system::update_external_Hfield_half() {
	double local_time;

	local_time = pt_glb->time_device - pt_glb->dt * 0.5;
	pt_glb->Hext[0] = pt_glb->Hext_stat[0] + pt_glb->Hext_altr[0] * sin(2. * PI * pt_glb->H_altr_freq[0] * local_time);
	pt_glb->Hext[1] = pt_glb->Hext_stat[1] + pt_glb->Hext_altr[1] * sin(2. * PI * pt_glb->H_altr_freq[1] * local_time);
	pt_glb->Hext[2] = pt_glb->Hext_stat[2] + pt_glb->Hext_altr[2] * sin(2. * PI * pt_glb->H_altr_freq[2] * local_time);
}

#pragma acc routine seq
void magnetic_system::update_external_Hfield_full() {
	double local_time;

	local_time = pt_glb->time_device;
	pt_glb->Hext[0] = pt_glb->Hext_stat[0] + pt_glb->Hext_altr[0] * sin(2. * PI * pt_glb->H_altr_freq[0] * local_time);
	pt_glb->Hext[1] = pt_glb->Hext_stat[1] + pt_glb->Hext_altr[1] * sin(2. * PI * pt_glb->H_altr_freq[1] * local_time);
	pt_glb->Hext[2] = pt_glb->Hext_stat[2] + pt_glb->Hext_altr[2] * sin(2. * PI * pt_glb->H_altr_freq[2] * local_time);
}

void magnetic_system::get_averagem() {
	double ax, ay, az;
	ax = 0.; ay = 0.; az = 0.;

#pragma acc parallel default(present)
	{
#pragma acc loop gang vector reduction(+:ax)
		for (long int id = 0; id < n; id++) {
			ax = ax + mx_glb(id);
		}

#pragma acc loop gang vector reduction(+:ay)
		for (long int id = 0; id < n; id++) {
			ay = ay + my_glb(id);
		}

#pragma acc loop gang vector reduction(+:az)
		for (long int id = 0; id < n; id++) {
			az = az + mz_glb(id);
		}
	}

	mx_ave = ax / static_cast<double>(pt_glb->NFM);
	my_ave = ay / static_cast<double>(pt_glb->NFM);
	mz_ave = az / static_cast<double>(pt_glb->NFM);
}

#pragma acc routine seq
void magnetic_system::update_prescribed_m_half() {
	unsigned int mat_type;
	material* mat;
	double m0, omega;

	m0 = sin(pt_glb->precess_angle / 180. * PI);
	omega = 2. * PI * pt_glb->precess_frequency;
	for (long i = 0; i < nx; i++) {
		for (long j = 0; j < ny; j++) {
			mat_type = pt_glb->material_cell(i, j, 0);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_FM == true) {
					mx_glb_store(i, j, 0) = sqrt((1. - m0 * m0));
					my_glb_store(i, j, 0) = m0 * sin(omega * (pt_glb->time_device - pt_glb->dt * 0.5));
					mz_glb_store(i, j, 0) = m0 * cos(omega * (pt_glb->time_device - pt_glb->dt * 0.5));
				}
			}
		}
	}
}

#pragma acc routine seq
void magnetic_system::update_prescribed_m_full() {
	unsigned int mat_type;
	material* mat;
	double m0, omega;

	m0 = sin(pt_glb->precess_angle / 180. * PI);
	omega = 2. * PI * pt_glb->precess_frequency;
	for (long i = 0; i < nx; i++) {
		for (long j = 0; j < ny; j++) {
			mat_type = pt_glb->material_cell(i, j, 0);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_FM == true) {
					mx_glb_store(i, j, 0) = sqrt((1. - m0 * m0));
					my_glb_store(i, j, 0) = m0 * sin(omega * pt_glb->time_device);
					mz_glb_store(i, j, 0) = m0 * cos(omega * pt_glb->time_device);
				}
			}
		}
	}
}

#pragma acc routine seq
void magnetic_system::update_prescribed_m() {
	unsigned int mat_type;
	material* mat;
	double m0, omega;

	m0 = sin(pt_glb->precess_angle / 180. * PI);
	omega = 2. * PI * pt_glb->precess_frequency;
	for (long i = 0; i < nx; i++) {
		for (long j = 0; j < ny; j++) {
			mat_type = pt_glb->material_cell(i, j, 0);
			if (mat_type != 0) {
				mat = &(pt_glb->material_parameters[(mat_type)-1]);
				if (mat->if_FM == true) {
					mx_glb(i, j, 0) = sqrt((1. - m0 * m0));
					my_glb(i, j, 0) = m0 * sin(omega * pt_glb->time_device);
					mz_glb(i, j, 0) = m0 * cos(omega * pt_glb->time_device);
				}
			}
		}
	}
}
