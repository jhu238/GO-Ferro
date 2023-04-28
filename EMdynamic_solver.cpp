#include "EMdynamic_system.h"

//void EMdynamic_system::get_dE() {
//	double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
//	double dP;
//
//	long int i, j, k;
//	unsigned int mat_type;
//	material* mat;
//	long int idx1, idy1, idz1;
//	long int idx2, idy2, idz2;
//	long int idx3, idy3, idz3;
//	long int idx4, idy4, idz4;
//
//	double er;
//	double cond;
//	double Jf, Jp, Jishe;
//
//	if (pt_glb->if_Jf_test == true) {
//		update_Jf_test();
//	}
//#pragma acc wait
//
//	//-----------X component-----------//
//#pragma acc parallel loop gang vector \
//private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
//er,cond, Jf,Jp, Jishe, dP,mat_type,mat,dHzdy,dHydz,i,j,k) default(present) async(1)
//	for (long int id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
//		i = id / ((ny + 1) * (nz + 1));
//		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
//		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);
//
//		idx1 = i; idx2 = i; idx3 = i; idx4 = i;
//		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
//			continue; // should update Ex using Liao's absorbing boundary surface
//		}
//		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
//			idy1 = ny - 1; idy2 = ny - 1; idy3 = 0; idy4 = 0;
//		}
//		else {
//			idy1 = j - 1; idy2 = j - 1; idy3 = j; idy4 = j;
//		}
//		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
//			continue; // should update Ex using Liao's absorbing boundary surface
//		}
//		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
//			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
//		}
//		else {
//			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
//		}
//
//		dHzdy = (DHz_em(idx1, idy3, idz4) - DHz_em(idx1, idy1, idz4)) / dy;
//		dHydz = (DHy_em(idx1, idy3, idz2) - DHy_em(idx1, idy3, idz1)) / dz;
//
//		Jf = 0.; Jp = 0.; Jishe = 0.;
//		dP = 0.; er = 0.; cond = 0.;
//		//----------Polarization 1------//
//		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP + pt_fe->dpx_glb(idx1, idy1, idz1);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfx(idx1, idy1, idz1);
//		Jp = Jp + Jpx_n1(idx1, idy1, idz1) + Jpx_n2(idx1, idy1, idz1) + Jpx_n3(idx1, idy1, idz1);
//		Jishe = Jishe + pt_mag->Jx_ISHE(idx1, idy1, idz1);
//
//		//----------Polarization 2------//
//		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+ pt_fe->dpx_glb(idx2, idy2, idz2);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfx(idx2, idy2, idz2);
//		Jp = Jp + Jpx_n1(idx2, idy2, idz2) + Jpx_n2(idx2, idy2, idz2) + Jpx_n3(idx2, idy2, idz2);
//		Jishe = Jishe + pt_mag->Jx_ISHE(idx2, idy2, idz2);
//
//		//----------Polarization 3------//
//		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+ pt_fe->dpx_glb(idx3, idy3, idz3);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfx(idx3, idy3, idz3);
//		Jp = Jp + Jpx_n1(idx3, idy3, idz3) + Jpx_n2(idx3, idy3, idz3) + Jpx_n3(idx3, idy3, idz3);
//		Jishe = Jishe + pt_mag->Jx_ISHE(idx3, idy3, idz3);
//
//		//----------Polarization 4------//
//		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+ pt_fe->dpx_glb(idx4, idy4, idz4);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfx(idx4, idy4, idz4);
//		Jp = Jp + Jpx_n1(idx4, idy4, idz4) + Jpx_n2(idx4, idy4, idz4) + Jpx_n3(idx4, idy4, idz4);
//		Jishe = Jishe + pt_mag->Jx_ISHE(idx4, idy4, idz4);
//
//		Jf = Jf / 4.;
//		Jp = Jp / 4.;
//		Jishe = Jishe / 4.;
//		dP = dP / 4.; er = er / 4.; cond = cond / 4.;
//
//		dDEx_em(i, j, k) = pt_glb->dt / e0 / er * \
//			(dHzdy - dHydz - Jf - Jp - Jishe - cond * DEx_em(i, j, k) - dP/ pt_glb->dt);
//	}
//	//-----------END X component-----------//
//
//	//-----------Y component-----------//
//#pragma acc parallel loop gang vector \
//private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
//Jf,Jp,Jishe,er,cond,dP,mat_type,mat,dHxdz,dHzdx,i,j,k) default(present) async(2)
//	for (long int id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
//		i = id / ((ny) * (nz + 1));
//		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
//		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);
//
//		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
//			continue; // should update Ex using Liao's absorbing boundary surface
//		}
//		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
//			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
//		}
//		else {
//			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
//		}
//
//		idy1 = j; idy2 = j; idy3 = j; idy4 = j;
//
//		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
//			continue; // should update Ex using Liao's absorbing boundary surface
//		}
//		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
//			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
//		}
//		else {
//			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
//		}
//
//		dHxdz = (DHx_em(idx3, idy1, idz4) - DHx_em(idx3, idy1, idz3)) / dz;
//		dHzdx = (DHz_em(idx3, idy1, idz2) - DHz_em(idx1, idy1, idz2)) / dx;
//
//		Jf = 0.; Jp = 0.; Jishe = 0.;
//		dP = 0.; er = 0.; cond = 0.;
//		//----------Polarization 1------//
//		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP =dP+ pt_fe->dpy_glb(idx1, idy1, idz1);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfy(idx1, idy1, idz1);
//		Jp = Jp + Jpy_n1(idx1, idy1, idz1) + Jpy_n2(idx1, idy1, idz1) + Jpy_n3(idx1, idy1, idz1);
//		Jishe = Jishe + pt_mag->Jy_ISHE(idx1, idy1, idz1);
//
//		//----------Polarization 2------//
//		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+ pt_fe->dpy_glb(idx2, idy2, idz2);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfy(idx2, idy2, idz2);
//		Jp = Jp + Jpy_n1(idx2, idy2, idz2) + Jpy_n2(idx2, idy2, idz2) + Jpy_n3(idx2, idy2, idz2);
//		Jishe = Jishe + pt_mag->Jy_ISHE(idx2, idy2, idz2);
//
//		//----------Polarization 3------//
//		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+ pt_fe->dpy_glb(idx3, idy3, idz3);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfy(idx3, idy3, idz3);
//		Jp = Jp + Jpy_n1(idx3, idy3, idz3) + Jpy_n2(idx3, idy3, idz3) + Jpy_n3(idx3, idy3, idz3);
//		Jishe = Jishe + pt_mag->Jy_ISHE(idx3, idy3, idz3);
//
//		//----------Polarization 4------//
//		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+ pt_fe->dpy_glb(idx4, idy4, idz4);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfy(idx4, idy4, idz4);
//		Jp = Jp + Jpy_n1(idx4, idy4, idz4) + Jpy_n2(idx4, idy4, idz4) + Jpy_n3(idx4, idy4, idz4);
//		Jishe = Jishe + pt_mag->Jy_ISHE(idx4, idy4, idz4);
//
//		Jf = Jf / 4.;
//		Jp = Jp / 4.;
//		Jishe = Jishe / 4.;
//		dP = dP / 4.; er = er / 4.; cond = cond / 4.;
//
//		dDEy_em(i, j, k) = pt_glb->dt / e0 / er * \
//			(dHxdz - dHzdx - Jf - Jp - Jishe -cond * DEy_em(i, j, k) - dP/ pt_glb->dt);
//	}
//
//	//-----------END Y component-----------//
//
//	//-----------Z component-----------//
//#pragma acc parallel loop gang vector \
//private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
//er,cond, dP, Jf, Jp, mat_type,mat,dHydx,dHxdy,i,j,k) default(present) async(3)
//	for (long int id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
//		i = id / ((ny + 1) * (nz));
//		j = (id - i * ((ny + 1) * (nz))) / (nz);
//		k = id - i * ((ny + 1) * (nz)) - j * (nz);
//
//		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
//			continue; // should update Ex using Liao's absorbing boundary surface
//		}
//		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
//			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
//		}
//		else {
//			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
//		}
//
//		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
//			continue; // should update Ex using Liao's absorbing boundary surface
//		}
//		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
//			idy1 = ny - 1; idy2 = 0; idy3 = ny - 1; idy4 = 0;
//		}
//		else {
//			idy1 = j - 1; idy2 = j; idy3 = j - 1; idy4 = j;
//		}
//
//		idz1 = k; idz2 = k; idz3 = k; idz4 = k;
//
//		dHydx = (DHy_em(idx3, idy2, idz1) - DHy_em(idx1, idy2, idz1)) / dx;
//		dHxdy = (DHx_em(idx3, idy2, idz1) - DHx_em(idx3, idy1, idz1)) / dy;
//
//		Jf = 0.; Jp = 0.;
//		dP = 0.; er = 0.; cond = 0.;
//		//----------Polarization 1------//
//		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+pt_fe->dpz_glb(idx1, idy1, idz1);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfz(idx1, idy1, idz1);
//		Jp = Jp + Jpz_n1(idx1, idy1, idz1) + Jpz_n2(idx1, idy1, idz1) + Jpz_n3(idx1, idy1, idz1);
//
//		//----------Polarization 2------//
//		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP = dP+ pt_fe->dpz_glb(idx2, idy2, idz2);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfz(idx2, idy2, idz2);
//		Jp = Jp + Jpz_n1(idx2, idy2, idz2) + Jpz_n2(idx2, idy2, idz2) + Jpz_n3(idx2, idy2, idz2);
//
//		//----------Polarization 3------//
//		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP =dP+ pt_fe->dpz_glb(idx3, idy3, idz3);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfz(idx3, idy3, idz3);
//		Jp = Jp + Jpz_n1(idx3, idy3, idz3) + Jpz_n2(idx3, idy3, idz3) + Jpz_n3(idx3, idy3, idz3);
//
//		//----------Polarization 4------//
//		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
//		if (mat_type == 0) {
//			er = er + 1.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			dP =dP+ pt_fe->dpz_glb(idx4, idy4, idz4);
//			er = er + mat->r_permittivity;
//			cond = cond + mat->conductivity;
//		}
//		Jf = Jf + DJfz(idx4, idy4, idz4);
//		Jp = Jp + Jpz_n1(idx4, idy4, idz4) + Jpz_n2(idx4, idy4, idz4) + Jpz_n3(idx4, idy4, idz4);
//
//		Jf = Jf / 4.;
//		Jp = Jp / 4.;
//		dP = dP / 4.; er = er / 4.; cond = cond / 4.;
//
//		dDEz_em(i, j, k) = pt_glb->dt / e0 / er * \
//			(dHydx - dHxdy - Jf - Jp - cond * DEz_em(i, j, k) -dP/ pt_glb->dt);
//	}
//	//-----------END Z component-----------//
//}

void EMdynamic_system::get_dE_RK1() {
	double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
	double dP;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idy1, idz1;
	long int idx2, idy2, idz2;
	long int idx3, idy3, idz3;
	long int idx4, idy4, idz4;

	double er;
	double cond;
	double Jf, Jp, Jishe;
	//double Jf_count, Jp_count, Ji_count;
	double P_count;

	if (pt_glb->if_Jf_input == true) {
		update_Jf_input(); //RK
	}

#pragma acc parallel default(present) async(8)
	{
		get_dJp_RK1();
	}

#pragma acc wait
	//-----------X component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
er,cond, Jf,Jp, Jishe, dP,mat_type,mat,dHzdy,dHydz,i,j,k, P_count) default(present) async(1)
	for (long int id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
		i = id / ((ny + 1) * (nz + 1));
		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);

		idx1 = i; idx2 = i; idx3 = i; idx4 = i;
		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			idy1 = ny - 1; idy2 = ny - 1; idy3 = 0; idy4 = 0;
		}
		else {
			idy1 = j - 1; idy2 = j - 1; idy3 = j; idy4 = j;
		}
		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
		}
		else {
			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
		}

		dHzdy = (DHz_em_store(idx1, idy3, idz4) - DHz_em_store(idx1, idy1, idz4)) / dy;
		dHydz = (DHy_em_store(idx1, idy3, idz2) - DHy_em_store(idx1, idy3, idz1)) / dz;

		Jf = 0.; Jp = 0.; Jishe = 0.;
		dP = 0.; er = 0.; cond = 0.;
		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		P_count = 0.;
		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk1(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			//if (mat->if_spin_pump == true) {
			//	Ji_count = Ji_count + 1.;
			//}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx1, idy1, idz1);
		Jp = Jp + Jpx_n1_store(idx1, idy1, idz1) + Jpx_n2_store(idx1, idy1, idz1) + Jpx_n3_store(idx1, idy1, idz1);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx1, idy1, idz1);

		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk1(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx2, idy2, idz2);
		Jp = Jp + Jpx_n1_store(idx2, idy2, idz2) + Jpx_n2_store(idx2, idy2, idz2) + Jpx_n3_store(idx2, idy2, idz2);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx2, idy2, idz2);

		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk1(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx3, idy3, idz3);
		Jp = Jp + Jpx_n1_store(idx3, idy3, idz3) + Jpx_n2_store(idx3, idy3, idz3) + Jpx_n3_store(idx3, idy3, idz3);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx3, idy3, idz3);

		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk1(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx4, idy4, idz4);
		Jp = Jp + Jpx_n1_store(idx4, idy4, idz4) + Jpx_n2_store(idx4, idy4, idz4) + Jpx_n3_store(idx4, idy4, idz4);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx4, idy4, idz4);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jf = Jf / 4.;
		Jp = Jp / 4.;
		Jishe = Jishe / 4.;
		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

		dDEx_em_rk1(i, j, k) = pt_glb->dt / e0 / er * \
			(dHzdy - dHydz - Jf - Jp - Jishe - cond * DEx_em_store(i, j, k) - dP / pt_glb->dt); //RK
	}
	//-----------END X component-----------//

	//-----------Y component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
Jf,Jp,Jishe,er,cond,dP,mat_type,mat,dHxdz,dHzdx,i,j,k,P_count) default(present) async(2)
	for (long int id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
		}
		else {
			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
		}

		idy1 = j; idy2 = j; idy3 = j; idy4 = j;

		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
		}
		else {
			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
		}

		dHxdz = (DHx_em_store(idx3, idy1, idz4) - DHx_em_store(idx3, idy1, idz3)) / dz;
		dHzdx = (DHz_em_store(idx3, idy1, idz2) - DHz_em_store(idx1, idy1, idz2)) / dx;

		Jf = 0.; Jp = 0.; Jishe = 0.;
		dP = 0.; er = 0.; cond = 0.;
		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		P_count = 0.;
		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk1(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx1, idy1, idz1);
		Jp = Jp + Jpy_n1_store(idx1, idy1, idz1) + Jpy_n2_store(idx1, idy1, idz1) + Jpy_n3_store(idx1, idy1, idz1);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx1, idy1, idz1);

		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk1(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx2, idy2, idz2);
		Jp = Jp + Jpy_n1_store(idx2, idy2, idz2) + Jpy_n2_store(idx2, idy2, idz2) + Jpy_n3_store(idx2, idy2, idz2);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx2, idy2, idz2);

		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk1(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx3, idy3, idz3);
		Jp = Jp + Jpy_n1_store(idx3, idy3, idz3) + Jpy_n2_store(idx3, idy3, idz3) + Jpy_n3_store(idx3, idy3, idz3);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx3, idy3, idz3);

		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk1(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx4, idy4, idz4);
		Jp = Jp + Jpy_n1_store(idx4, idy4, idz4) + Jpy_n2_store(idx4, idy4, idz4) + Jpy_n3_store(idx4, idy4, idz4);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx4, idy4, idz4);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jf = Jf / 4.;
		Jp = Jp / 4.;
		Jishe = Jishe / 4.;
		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

		dDEy_em_rk1(i, j, k) = pt_glb->dt / e0 / er * \
			(dHxdz - dHzdx - Jf - Jp - Jishe - cond * DEy_em_store(i, j, k) - dP / pt_glb->dt); //RK
	}

	//-----------END Y component-----------//

	//-----------Z component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
er,cond, dP, Jf, Jp, mat_type,mat,dHydx,dHxdy,i,j,k,P_count) default(present) async(3)
	for (long int id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
		}
		else {
			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
		}

		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			idy1 = ny - 1; idy2 = 0; idy3 = ny - 1; idy4 = 0;
		}
		else {
			idy1 = j - 1; idy2 = j; idy3 = j - 1; idy4 = j;
		}

		idz1 = k; idz2 = k; idz3 = k; idz4 = k;

		dHydx = (DHy_em_store(idx3, idy2, idz1) - DHy_em_store(idx1, idy2, idz1)) / dx;
		dHxdy = (DHx_em_store(idx3, idy2, idz1) - DHx_em_store(idx3, idy1, idz1)) / dy;

		Jf = 0.; Jp = 0.;
		dP = 0.; er = 0.; cond = 0.;
		//Jf_count = 0.; Jp_count = 0.; 
		P_count = 0.;
		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk1(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx1, idy1, idz1);
		Jp = Jp + Jpz_n1_store(idx1, idy1, idz1) + Jpz_n2_store(idx1, idy1, idz1) + Jpz_n3_store(idx1, idy1, idz1);

		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk1(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx2, idy2, idz2);
		Jp = Jp + Jpz_n1_store(idx2, idy2, idz2) + Jpz_n2_store(idx2, idy2, idz2) + Jpz_n3_store(idx2, idy2, idz2);

		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk1(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx3, idy3, idz3);
		Jp = Jp + Jpz_n1_store(idx3, idy3, idz3) + Jpz_n2_store(idx3, idy3, idz3) + Jpz_n3_store(idx3, idy3, idz3);

		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk1(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx4, idy4, idz4);
		Jp = Jp + Jpz_n1_store(idx4, idy4, idz4) + Jpz_n2_store(idx4, idy4, idz4) + Jpz_n3_store(idx4, idy4, idz4);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jf = Jf / 4.;
		Jp = Jp / 4.;
		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

		dDEz_em_rk1(i, j, k) = pt_glb->dt / e0 / er * \
			(dHydx - dHxdy - Jf - Jp - cond * DEz_em_store(i, j, k) - dP / pt_glb->dt); //RK
	}
	//-----------END Z component-----------//
}

void EMdynamic_system::get_dE_RK2() {
	double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
	double dP;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idy1, idz1;
	long int idx2, idy2, idz2;
	long int idx3, idy3, idz3;
	long int idx4, idy4, idz4;

	double er;
	double cond;
	double Jf, Jp, Jishe;
	//double Jf_count, Jp_count, Ji_count;
	double P_count;

	if (pt_glb->if_Jf_input == true) {
		update_Jf_input_half(); //RK
	}

#pragma acc parallel default(present) async(8)
	{
		get_dJp_RK2();
	}

#pragma acc wait
	//-----------X component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
er,cond, Jf,Jp, Jishe, dP,mat_type,mat,dHzdy,dHydz,i,j,k,P_count) default(present) async(1)
	for (long int id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
		i = id / ((ny + 1) * (nz + 1));
		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);

		idx1 = i; idx2 = i; idx3 = i; idx4 = i;
		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			idy1 = ny - 1; idy2 = ny - 1; idy3 = 0; idy4 = 0;
		}
		else {
			idy1 = j - 1; idy2 = j - 1; idy3 = j; idy4 = j;
		}
		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
		}
		else {
			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
		}

		dHzdy = (DHz_em_store(idx1, idy3, idz4) - DHz_em_store(idx1, idy1, idz4)) / dy;
		dHydz = (DHy_em_store(idx1, idy3, idz2) - DHy_em_store(idx1, idy3, idz1)) / dz;

		Jf = 0.; Jp = 0.; Jishe = 0.;
		dP = 0.; er = 0.; cond = 0.;
		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		P_count = 0.;
		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk2(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx1, idy1, idz1);
		Jp = Jp + Jpx_n1_store(idx1, idy1, idz1) + Jpx_n2_store(idx1, idy1, idz1) + Jpx_n3_store(idx1, idy1, idz1);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx1, idy1, idz1);

		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk2(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx2, idy2, idz2);
		Jp = Jp + Jpx_n1_store(idx2, idy2, idz2) + Jpx_n2_store(idx2, idy2, idz2) + Jpx_n3_store(idx2, idy2, idz2);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx2, idy2, idz2);

		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk2(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx3, idy3, idz3);
		Jp = Jp + Jpx_n1_store(idx3, idy3, idz3) + Jpx_n2_store(idx3, idy3, idz3) + Jpx_n3_store(idx3, idy3, idz3);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx3, idy3, idz3);

		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk2(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx4, idy4, idz4);
		Jp = Jp + Jpx_n1_store(idx4, idy4, idz4) + Jpx_n2_store(idx4, idy4, idz4) + Jpx_n3_store(idx4, idy4, idz4);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx4, idy4, idz4);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jf = Jf / 4.;
		Jp = Jp / 4.;
		Jishe = Jishe / 4.;
		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

		dDEx_em_rk2(i, j, k) = pt_glb->dt / e0 / er * \
			(dHzdy - dHydz - Jf - Jp - Jishe - cond * DEx_em_store(i, j, k) - dP / pt_glb->dt); //RK
	}
	//-----------END X component-----------//

	//-----------Y component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
Jf,Jp,Jishe,er,cond,dP,mat_type,mat,dHxdz,dHzdx,i,j,k,P_count) default(present) async(2)
	for (long int id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
		}
		else {
			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
		}

		idy1 = j; idy2 = j; idy3 = j; idy4 = j;

		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
		}
		else {
			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
		}

		dHxdz = (DHx_em_store(idx3, idy1, idz4) - DHx_em_store(idx3, idy1, idz3)) / dz;
		dHzdx = (DHz_em_store(idx3, idy1, idz2) - DHz_em_store(idx1, idy1, idz2)) / dx;

		Jf = 0.; Jp = 0.; Jishe = 0.;
		dP = 0.; er = 0.; cond = 0.;
		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		P_count = 0.;
		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk2(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx1, idy1, idz1);
		Jp = Jp + Jpy_n1_store(idx1, idy1, idz1) + Jpy_n2_store(idx1, idy1, idz1) + Jpy_n3_store(idx1, idy1, idz1);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx1, idy1, idz1);

		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk2(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx2, idy2, idz2);
		Jp = Jp + Jpy_n1_store(idx2, idy2, idz2) + Jpy_n2_store(idx2, idy2, idz2) + Jpy_n3_store(idx2, idy2, idz2);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx2, idy2, idz2);

		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk2(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx3, idy3, idz3);
		Jp = Jp + Jpy_n1_store(idx3, idy3, idz3) + Jpy_n2_store(idx3, idy3, idz3) + Jpy_n3_store(idx3, idy3, idz3);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx3, idy3, idz3);

		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk2(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx4, idy4, idz4);
		Jp = Jp + Jpy_n1_store(idx4, idy4, idz4) + Jpy_n2_store(idx4, idy4, idz4) + Jpy_n3_store(idx4, idy4, idz4);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx4, idy4, idz4);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jf = Jf / 4.;
		Jp = Jp / 4.;
		Jishe = Jishe / 4.;
		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

		dDEy_em_rk2(i, j, k) = pt_glb->dt / e0 / er * \
			(dHxdz - dHzdx - Jf - Jp - Jishe - cond * DEy_em_store(i, j, k) - dP / pt_glb->dt); //RK
	}

	//-----------END Y component-----------//

	//-----------Z component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
er,cond, dP, Jf, Jp, mat_type,mat,dHydx,dHxdy,i,j,k,P_count) default(present) async(3)
	for (long int id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
		}
		else {
			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
		}

		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			idy1 = ny - 1; idy2 = 0; idy3 = ny - 1; idy4 = 0;
		}
		else {
			idy1 = j - 1; idy2 = j; idy3 = j - 1; idy4 = j;
		}

		idz1 = k; idz2 = k; idz3 = k; idz4 = k;

		dHydx = (DHy_em_store(idx3, idy2, idz1) - DHy_em_store(idx1, idy2, idz1)) / dx;
		dHxdy = (DHx_em_store(idx3, idy2, idz1) - DHx_em_store(idx3, idy1, idz1)) / dy;

		Jf = 0.; Jp = 0.;
		dP = 0.; er = 0.; cond = 0.;
		//Jf_count = 0.; Jp_count = 0.; 
		P_count = 0.;
		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk2(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx1, idy1, idz1);
		Jp = Jp + Jpz_n1_store(idx1, idy1, idz1) + Jpz_n2_store(idx1, idy1, idz1) + Jpz_n3_store(idx1, idy1, idz1);

		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk2(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx2, idy2, idz2);
		Jp = Jp + Jpz_n1_store(idx2, idy2, idz2) + Jpz_n2_store(idx2, idy2, idz2) + Jpz_n3_store(idx2, idy2, idz2);

		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk2(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx3, idy3, idz3);
		Jp = Jp + Jpz_n1_store(idx3, idy3, idz3) + Jpz_n2_store(idx3, idy3, idz3) + Jpz_n3_store(idx3, idy3, idz3);

		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk2(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx4, idy4, idz4);
		Jp = Jp + Jpz_n1_store(idx4, idy4, idz4) + Jpz_n2_store(idx4, idy4, idz4) + Jpz_n3_store(idx4, idy4, idz4);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jf = Jf / 4.;
		Jp = Jp / 4.;
		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

		dDEz_em_rk2(i, j, k) = pt_glb->dt / e0 / er * \
			(dHydx - dHxdy - Jf - Jp - cond * DEz_em_store(i, j, k) - dP / pt_glb->dt); //RK
	}
	//-----------END Z component-----------//
}

void EMdynamic_system::get_dE_RK3() {
	double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
	double dP;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idy1, idz1;
	long int idx2, idy2, idz2;
	long int idx3, idy3, idz3;
	long int idx4, idy4, idz4;

	double er;
	double cond;
	double Jf, Jp, Jishe;
	//double Jf_count, Jp_count, Ji_count;
	double P_count;

	if (pt_glb->if_Jf_input == true) {
		update_Jf_input_half(); //RK
	}

#pragma acc parallel default(present) async(8)
	{
		get_dJp_RK3();
	}

#pragma acc wait
	//-----------X component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
er,cond, Jf,Jp, Jishe, dP,mat_type,mat,dHzdy,dHydz,i,j,k,P_count) default(present) async(1)
	for (long int id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
		i = id / ((ny + 1) * (nz + 1));
		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);

		idx1 = i; idx2 = i; idx3 = i; idx4 = i;
		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			idy1 = ny - 1; idy2 = ny - 1; idy3 = 0; idy4 = 0;
		}
		else {
			idy1 = j - 1; idy2 = j - 1; idy3 = j; idy4 = j;
		}
		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
		}
		else {
			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
		}

		dHzdy = (DHz_em_store(idx1, idy3, idz4) - DHz_em_store(idx1, idy1, idz4)) / dy;
		dHydz = (DHy_em_store(idx1, idy3, idz2) - DHy_em_store(idx1, idy3, idz1)) / dz;

		Jf = 0.; Jp = 0.; Jishe = 0.;
		dP = 0.; er = 0.; cond = 0.;
		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		P_count = 0.;
		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk3(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx1, idy1, idz1);
		Jp = Jp + Jpx_n1_store(idx1, idy1, idz1) + Jpx_n2_store(idx1, idy1, idz1) + Jpx_n3_store(idx1, idy1, idz1);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx1, idy1, idz1);

		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk3(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx2, idy2, idz2);
		Jp = Jp + Jpx_n1_store(idx2, idy2, idz2) + Jpx_n2_store(idx2, idy2, idz2) + Jpx_n3_store(idx2, idy2, idz2);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx2, idy2, idz2);

		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk3(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx3, idy3, idz3);
		Jp = Jp + Jpx_n1_store(idx3, idy3, idz3) + Jpx_n2_store(idx3, idy3, idz3) + Jpx_n3_store(idx3, idy3, idz3);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx3, idy3, idz3);

		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk3(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx4, idy4, idz4);
		Jp = Jp + Jpx_n1_store(idx4, idy4, idz4) + Jpx_n2_store(idx4, idy4, idz4) + Jpx_n3_store(idx4, idy4, idz4);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx4, idy4, idz4);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jf = Jf / 4.;
		Jp = Jp / 4.;
		Jishe = Jishe / 4.;
		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

		dDEx_em_rk3(i, j, k) = pt_glb->dt / e0 / er * \
			(dHzdy - dHydz - Jf - Jp - Jishe - cond * DEx_em_store(i, j, k) - dP / pt_glb->dt); //RK
	}
	//-----------END X component-----------//

	//-----------Y component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
Jf,Jp,Jishe,er,cond,dP,mat_type,mat,dHxdz,dHzdx,i,j,k,P_count) default(present) async(2)
	for (long int id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
		}
		else {
			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
		}

		idy1 = j; idy2 = j; idy3 = j; idy4 = j;

		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
		}
		else {
			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
		}

		dHxdz = (DHx_em_store(idx3, idy1, idz4) - DHx_em_store(idx3, idy1, idz3)) / dz;
		dHzdx = (DHz_em_store(idx3, idy1, idz2) - DHz_em_store(idx1, idy1, idz2)) / dx;

		Jf = 0.; Jp = 0.; Jishe = 0.;
		dP = 0.; er = 0.; cond = 0.;
		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		P_count = 0.;
		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk3(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx1, idy1, idz1);
		Jp = Jp + Jpy_n1_store(idx1, idy1, idz1) + Jpy_n2_store(idx1, idy1, idz1) + Jpy_n3_store(idx1, idy1, idz1);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx1, idy1, idz1);

		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk3(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx2, idy2, idz2);
		Jp = Jp + Jpy_n1_store(idx2, idy2, idz2) + Jpy_n2_store(idx2, idy2, idz2) + Jpy_n3_store(idx2, idy2, idz2);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx2, idy2, idz2);

		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk3(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx3, idy3, idz3);
		Jp = Jp + Jpy_n1_store(idx3, idy3, idz3) + Jpy_n2_store(idx3, idy3, idz3) + Jpy_n3_store(idx3, idy3, idz3);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx3, idy3, idz3);

		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk3(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx4, idy4, idz4);
		Jp = Jp + Jpy_n1_store(idx4, idy4, idz4) + Jpy_n2_store(idx4, idy4, idz4) + Jpy_n3_store(idx4, idy4, idz4);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx4, idy4, idz4);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jf = Jf / 4.;
		Jp = Jp / 4.;
		Jishe = Jishe / 4.;
		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

		dDEy_em_rk3(i, j, k) = pt_glb->dt / e0 / er * \
			(dHxdz - dHzdx - Jf - Jp - Jishe - cond * DEy_em_store(i, j, k) - dP / pt_glb->dt); //RK
	}

	//-----------END Y component-----------//

	//-----------Z component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
er,cond, dP, Jf, Jp, mat_type,mat,dHydx,dHxdy,i,j,k,P_count) default(present) async(3)
	for (long int id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
		}
		else {
			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
		}

		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			idy1 = ny - 1; idy2 = 0; idy3 = ny - 1; idy4 = 0;
		}
		else {
			idy1 = j - 1; idy2 = j; idy3 = j - 1; idy4 = j;
		}

		idz1 = k; idz2 = k; idz3 = k; idz4 = k;

		dHydx = (DHy_em_store(idx3, idy2, idz1) - DHy_em_store(idx1, idy2, idz1)) / dx;
		dHxdy = (DHx_em_store(idx3, idy2, idz1) - DHx_em_store(idx3, idy1, idz1)) / dy;

		Jf = 0.; Jp = 0.;
		dP = 0.; er = 0.; cond = 0.;
		//Jf_count = 0.; Jp_count = 0.; 
		P_count = 0.;
		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk3(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx1, idy1, idz1);
		Jp = Jp + Jpz_n1_store(idx1, idy1, idz1) + Jpz_n2_store(idx1, idy1, idz1) + Jpz_n3_store(idx1, idy1, idz1);

		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk3(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx2, idy2, idz2);
		Jp = Jp + Jpz_n1_store(idx2, idy2, idz2) + Jpz_n2_store(idx2, idy2, idz2) + Jpz_n3_store(idx2, idy2, idz2);

		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk3(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx3, idy3, idz3);
		Jp = Jp + Jpz_n1_store(idx3, idy3, idz3) + Jpz_n2_store(idx3, idy3, idz3) + Jpz_n3_store(idx3, idy3, idz3);

		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk3(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx4, idy4, idz4);
		Jp = Jp + Jpz_n1_store(idx4, idy4, idz4) + Jpz_n2_store(idx4, idy4, idz4) + Jpz_n3_store(idx4, idy4, idz4);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jf = Jf / 4.;
		Jp = Jp / 4.;
		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

		dDEz_em_rk3(i, j, k) = pt_glb->dt / e0 / er * \
			(dHydx - dHxdy - Jf - Jp - cond * DEz_em_store(i, j, k) - dP / pt_glb->dt); //RK
	}
	//-----------END Z component-----------//
}

void EMdynamic_system::get_dE_RK4() {
	double dHxdy, dHxdz, dHydx, dHydz, dHzdx, dHzdy;
	double dP;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idy1, idz1;
	long int idx2, idy2, idz2;
	long int idx3, idy3, idz3;
	long int idx4, idy4, idz4;

	double er;
	double cond;
	double Jf, Jp, Jishe;
	//double Jf_count, Jp_count, Ji_count;
	double P_count;

	if (pt_glb->if_Jf_input == true) {
		update_Jf_input_full(); //RK
	}
#pragma acc parallel default(present) async(8)
	{
		get_dJp_RK4();
	}

#pragma acc wait
	//-----------X component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
er,cond, Jf,Jp, Jishe, dP,mat_type,mat,dHzdy,dHydz,i,j,k,P_count) default(present) async(1)
	for (long int id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
		i = id / ((ny + 1) * (nz + 1));
		j = (id - i * ((ny + 1) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny + 1) * (nz + 1)) - j * (nz + 1);

		idx1 = i; idx2 = i; idx3 = i; idx4 = i;
		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			idy1 = ny - 1; idy2 = ny - 1; idy3 = 0; idy4 = 0;
		}
		else {
			idy1 = j - 1; idy2 = j - 1; idy3 = j; idy4 = j;
		}
		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
		}
		else {
			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
		}

		dHzdy = (DHz_em_store(idx1, idy3, idz4) - DHz_em_store(idx1, idy1, idz4)) / dy;
		dHydz = (DHy_em_store(idx1, idy3, idz2) - DHy_em_store(idx1, idy3, idz1)) / dz;

		Jf = 0.; Jp = 0.; Jishe = 0.;
		dP = 0.; er = 0.; cond = 0.;
		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		P_count = 0.;
		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk4(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx1, idy1, idz1);
		Jp = Jp + Jpx_n1_store(idx1, idy1, idz1) + Jpx_n2_store(idx1, idy1, idz1) + Jpx_n3_store(idx1, idy1, idz1);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx1, idy1, idz1);

		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk4(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx2, idy2, idz2);
		Jp = Jp + Jpx_n1_store(idx2, idy2, idz2) + Jpx_n2_store(idx2, idy2, idz2) + Jpx_n3_store(idx2, idy2, idz2);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx2, idy2, idz2);

		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk4(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx3, idy3, idz3);
		Jp = Jp + Jpx_n1_store(idx3, idy3, idz3) + Jpx_n2_store(idx3, idy3, idz3) + Jpx_n3_store(idx3, idy3, idz3);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx3, idy3, idz3);

		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpx_glb_rk4(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfx(idx4, idy4, idz4);
		Jp = Jp + Jpx_n1_store(idx4, idy4, idz4) + Jpx_n2_store(idx4, idy4, idz4) + Jpx_n3_store(idx4, idy4, idz4);
		Jishe = Jishe + pt_mag->Jx_ISHE(idx4, idy4, idz4);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jf = Jf / 4.;
		Jp = Jp / 4.;
		Jishe = Jishe / 4.;
		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

		dDEx_em_rk4(i, j, k) = pt_glb->dt / e0 / er * \
			(dHzdy - dHydz - Jf - Jp - Jishe - cond * DEx_em_store(i, j, k) - dP / pt_glb->dt); //RK
	}
	//-----------END X component-----------//

	//-----------Y component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
Jf,Jp,Jishe,er,cond,dP,mat_type,mat,dHxdz,dHzdx,i,j,k,P_count) default(present) async(2)
	for (long int id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
		}
		else {
			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
		}

		idy1 = j; idy2 = j; idy3 = j; idy4 = j;

		if ((k == 0 || k == nz) && pt_geo->periodicZ == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((k == 0 || k == nz) && pt_geo->periodicZ == true) {
			idz1 = nz - 1; idz2 = 0; idz3 = nz - 1; idz4 = 0;
		}
		else {
			idz1 = k - 1; idz2 = k; idz3 = k - 1; idz4 = k;
		}

		dHxdz = (DHx_em_store(idx3, idy1, idz4) - DHx_em_store(idx3, idy1, idz3)) / dz;
		dHzdx = (DHz_em_store(idx3, idy1, idz2) - DHz_em_store(idx1, idy1, idz2)) / dx;

		Jf = 0.; Jp = 0.; Jishe = 0.;
		dP = 0.; er = 0.; cond = 0.;
		//Jf_count = 0.; Jp_count = 0.; Ji_count = 0.; 
		P_count = 0.;
		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk4(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx1, idy1, idz1);
		Jp = Jp + Jpy_n1_store(idx1, idy1, idz1) + Jpy_n2_store(idx1, idy1, idz1) + Jpy_n3_store(idx1, idy1, idz1);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx1, idy1, idz1);

		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk4(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx2, idy2, idz2);
		Jp = Jp + Jpy_n1_store(idx2, idy2, idz2) + Jpy_n2_store(idx2, idy2, idz2) + Jpy_n3_store(idx2, idy2, idz2);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx2, idy2, idz2);

		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk4(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx3, idy3, idz3);
		Jp = Jp + Jpy_n1_store(idx3, idy3, idz3) + Jpy_n2_store(idx3, idy3, idz3) + Jpy_n3_store(idx3, idy3, idz3);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx3, idy3, idz3);

		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpy_glb_rk4(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfy(idx4, idy4, idz4);
		Jp = Jp + Jpy_n1_store(idx4, idy4, idz4) + Jpy_n2_store(idx4, idy4, idz4) + Jpy_n3_store(idx4, idy4, idz4);
		Jishe = Jishe + pt_mag->Jy_ISHE(idx4, idy4, idz4);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jf = Jf / 4.;
		Jp = Jp / 4.;
		Jishe = Jishe / 4.;
		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

		dDEy_em_rk4(i, j, k) = pt_glb->dt / e0 / er * \
			(dHxdz - dHzdx - Jf - Jp - Jishe - cond * DEy_em_store(i, j, k) - dP / pt_glb->dt); //RK
	}

	//-----------END Y component-----------//

	//-----------Z component-----------//
#pragma acc parallel loop gang vector \
private(idx1, idy1, idz1,idx2, idy2, idz2,idx3, idy3, idz3,idx4, idy4, idz4,\
er,cond, dP, Jf, Jp, mat_type,mat,dHydx,dHxdy,i,j,k,P_count) default(present) async(3)
	for (long int id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		if ((i == 0 || i == nx) && pt_geo->periodicX == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((i == 0 || i == nx) && pt_geo->periodicX == true) {
			idx1 = nx - 1; idx2 = nx - 1; idx3 = 0; idx4 = 0;
		}
		else {
			idx1 = i - 1; idx2 = i - 1; idx3 = i; idx4 = i;
		}

		if ((j == 0 || j == ny) && pt_geo->periodicY == false) {
			continue; // should update Ex using Liao's absorbing boundary surface
		}
		else if ((j == 0 || j == ny) && pt_geo->periodicY == true) {
			idy1 = ny - 1; idy2 = 0; idy3 = ny - 1; idy4 = 0;
		}
		else {
			idy1 = j - 1; idy2 = j; idy3 = j - 1; idy4 = j;
		}

		idz1 = k; idz2 = k; idz3 = k; idz4 = k;

		dHydx = (DHy_em_store(idx3, idy2, idz1) - DHy_em_store(idx1, idy2, idz1)) / dx;
		dHxdy = (DHx_em_store(idx3, idy2, idz1) - DHx_em_store(idx3, idy1, idz1)) / dy;

		Jf = 0.; Jp = 0.;
		dP = 0.; er = 0.; cond = 0.;
		//Jf_count = 0.; Jp_count = 0.; 
		P_count = 0.;
		//----------Polarization 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk4(idx1, idy1, idz1); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx1, idy1, idz1);
		Jp = Jp + Jpz_n1_store(idx1, idy1, idz1) + Jpz_n2_store(idx1, idy1, idz1) + Jpz_n3_store(idx1, idy1, idz1);

		//----------Polarization 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk4(idx2, idy2, idz2); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx2, idy2, idz2);
		Jp = Jp + Jpz_n1_store(idx2, idy2, idz2) + Jpz_n2_store(idx2, idy2, idz2) + Jpz_n3_store(idx2, idy2, idz2);

		//----------Polarization 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk4(idx3, idy3, idz3); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx3, idy3, idz3);
		Jp = Jp + Jpz_n1_store(idx3, idy3, idz3) + Jpz_n2_store(idx3, idy3, idz3) + Jpz_n3_store(idx3, idy3, idz3);

		//----------Polarization 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		if (mat_type == 0) {
			er = er + 1.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (pt_glb->if_EM_fromP == true) {
				dP = dP + pt_fe->dpz_glb_rk4(idx4, idy4, idz4); //RK
			}
			if (mat->if_FE == true) {
				P_count = P_count + 1.;
			}
			er = er + mat->r_permittivity;
			cond = cond + mat->conductivity;
		}
		Jf = Jf + DJfz(idx4, idy4, idz4);
		Jp = Jp + Jpz_n1_store(idx4, idy4, idz4) + Jpz_n2_store(idx4, idy4, idz4) + Jpz_n3_store(idx4, idy4, idz4);

		if (P_count < 0.5) {
			P_count = 1.;
		}

		Jf = Jf / 4.;
		Jp = Jp / 4.;
		dP = dP / P_count; er = er / 4.; cond = cond / 4.;

		dDEz_em_rk4(i, j, k) = pt_glb->dt / e0 / er * \
			(dHydx - dHxdy - Jf - Jp - cond * DEz_em_store(i, j, k) - dP / pt_glb->dt); //RK
	}
	//-----------END Z component-----------//
}

//void EMdynamic_system::get_dH() {
//	double dExdz, dExdy, dEydx, dEydz, dEzdx, dEzdy;
//	double dMx1, dMy1, dMz1;
//	double dMx2, dMy2, dMz2;
//	double Ms1, Ms2;
//
//	long int i, j, k;
//	unsigned int mat_type;
//	material* mat;
//	long int idx1, idy1, idz1;
//	long int idx2, idy2, idz2;
//
//	//------------------X component-------------------//
//#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
//Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dEzdy,dEydz,i,j,k) default(present) async(1)
//	for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
//		i = id / ((ny) * (nz));
//		j = (id - i * ((ny) * (nz))) / (nz);
//		k = id - i * ((ny) * (nz)) - j * (nz);
//
//		if (i == 0) {
//			if (pt_geo->periodicX == true) {
//				idx1 = nx - 1;
//			}
//			else {
//				idx1 = 0;
//			}
//		}
//		else {
//			idx1 = i - 1;
//		}
//		if (i == nx) {
//			if (pt_geo->periodicX == true) {
//				idx2 = 0;
//			}
//			else {
//				idx2 = nx - 1;
//			}
//		}
//		else {
//			idx2 = i;
//		}
//
//		idy1 = j; idy2 = j;
//		idz1 = k; idz2 = k;
//		dEzdy = (DEz_em(i, j + 1, k) - DEz_em(i, j, k)) / dy;
//		dEydz = (DEy_em(i, j, k + 1) - DEy_em(i, j, k)) / dz;
//
//		//----------Magnetization 1------//
//
//		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
//
//		if (mat_type == 0) {
//			Ms1 = 0.;
//			dMx1 = 0.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			Ms1 = mat->Ms;
//			dMx1 = pt_mag->dmx_glb(idx1, idy1, idz1) * Ms1;
//		}
//		//----------Magnetization 2------//
//
//		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
//
//		if (mat_type == 0) {
//			Ms2 = 0.;
//			dMx2 = 0.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			Ms2 = mat->Ms;
//			dMx2 = pt_mag->dmx_glb(idx2, idy2, idz2) * Ms2;
//		}
//		//---------END calculating temporal change of M--------//
//
//		dDHx_em(i, j, k) = -pt_glb->dt / mu0 * (dEzdy - dEydz) - (dMx1 + dMx2) * 0.5;
//		//if (abs(dDHx_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
//		//	dDHx_em(i, j, k) = 0.;
//		//}
//	}
//
//	//------------------END X component-------------------//
//
//	//------------------Y component-------------------//
//#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
//Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdz,dEzdx,i,j,k) default(present) async(2)
//	for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
//		i = id / ((ny + 1) * (nz));
//		j = (id - i * ((ny + 1) * (nz))) / (nz);
//		k = id - i * ((ny + 1) * (nz)) - j * (nz);
//
//		idx1 = i; idx2 = i;
//		if (j == 0) {
//			if (pt_geo->periodicY == true) {
//				idy1 = ny - 1;
//			}
//			else {
//				idy1 = 0;
//			}
//		}
//		else {
//			idy1 = j - 1;
//		}
//		if (j == ny) {
//			if (pt_geo->periodicY == true) {
//				idy2 = 0;
//			}
//			else {
//				idy2 = ny - 1;
//			}
//		}
//		else {
//			idy2 = j;
//		}
//		idz1 = k; idz2 = k;
//
//		dExdz = (DEx_em(i, j, k + 1) - DEx_em(i, j, k)) / dz;
//		dEzdx = (DEz_em(i + 1, j, k) - DEz_em(i, j, k)) / dx;
//
//		//----------Magnetization 1------//
//
//		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
//
//		if (mat_type == 0) {
//			dMy1 = 0.;
//			Ms1 = 0.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			Ms1 = mat->Ms;
//			dMy1 = pt_mag->dmy_glb(idx1, idy1, idz1) * Ms1;
//		}
//		//----------Magnetization 2------//
//
//		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
//
//		if (mat_type == 0) {
//			dMy2 = 0.;
//			Ms2 = 0.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			Ms2 = mat->Ms;
//			dMy2 = pt_mag->dmy_glb(idx2, idy2, idz2) * Ms2;
//		}
//		//---------END calculating temporal change of M--------//
//
//		dDHy_em(i, j, k) = -pt_glb->dt / mu0 * (dExdz - dEzdx) - (dMy1 + dMy2) * 0.5;
//		//if (abs(dDHy_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
//		//	dDHy_em(i, j, k) = 0.;
//		//}
//	}
//	//------------------END Y component-------------------//
//
//	//------------------Z component-------------------//
//#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
//Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdy,dEydx,i,j,k) default(present) async(3)
//	for (long int id = 0; id < nx * ny * (nz + 1); id++) {
//		i = id / ((ny) * (nz + 1));
//		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
//		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);
//
//		dEydx = (DEy_em(i + 1, j, k) - DEy_em(i, j, k)) / dx;
//		dExdy = (DEx_em(i, j + 1, k) - DEx_em(i, j, k)) / dy;
//
//		idx1 = i; idx2 = i;
//		idy1 = j; idy2 = j;
//		//----------Magnetization 1------//
//		if (k == 0) {
//			if (pt_geo->periodicZ == true) {
//				idz1 = nz - 1;
//			}
//			else {
//				idz1 = 0;
//			}
//		}
//		else {
//			idz1 = k - 1;
//		}
//		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
//
//		if (mat_type == 0) {
//			dMz1 = 0.;
//			Ms1 = 0.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			Ms1 = mat->Ms;
//			dMz1 = pt_mag->dmz_glb(idx1, idy1, idz1) * Ms1;
//		}
//		//----------Magnetization 2------//
//		if (k == nz) {
//			if (pt_geo->periodicZ == true) {
//				idz2 = 0;
//			}
//			else {
//				idz2 = nz - 1;
//			}
//		}
//		else {
//			idz2 = k;
//		}
//		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
//
//		if (mat_type == 0) {
//			dMz2 = 0.;
//			Ms2 = 0.;
//		}
//		else {
//			mat = &(pt_glb->material_parameters[(mat_type)-1]);
//			Ms2 = mat->Ms;
//			dMz2 = pt_mag->dmz_glb(idx2, idy2, idz2) * Ms2;
//		}
//		//---------END calculating temporal change of M--------//
//
//		dDHz_em(i, j, k) = -pt_glb->dt / mu0 * (dEydx - dExdy) - (dMz1 + dMz2) * 0.5;
//		//if (abs(dDHz_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
//		//	dDHz_em(i, j, k) = 0.;
//		//}
//	}
//	//------------------END Z component-------------------//
//}

void EMdynamic_system::get_dH_RK1() {
	double dExdz, dExdy, dEydx, dEydz, dEzdx, dEzdy;
	double dMx1, dMy1, dMz1;
	double dMx2, dMy2, dMz2;
	double Ms1, Ms2;
	double Ms1_AFM1, Ms2_AFM1;
	double Ms1_AFM2, Ms2_AFM2;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idy1, idz1;
	long int idx2, idy2, idz2;
	double M_count;

	//------------------X component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dEzdy,dEydz,i,j,k,M_count) default(present) async(4)
	for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
		i = id / ((ny) * (nz));
		j = (id - i * ((ny) * (nz))) / (nz);
		k = id - i * ((ny) * (nz)) - j * (nz);

		if (i == 0) {
			if (pt_geo->periodicX == true) {
				idx1 = nx - 1;
			}
			else {
				idx1 = 0;
			}
		}
		else {
			idx1 = i - 1;
		}
		if (i == nx) {
			if (pt_geo->periodicX == true) {
				idx2 = 0;
			}
			else {
				idx2 = nx - 1;
			}
		}
		else {
			idx2 = i;
		}

		idy1 = j; idy2 = j;
		idz1 = k; idz2 = k;
		dEzdy = (DEz_em_store(i, j + 1, k) - DEz_em_store(i, j, k)) / dy;
		dEydz = (DEy_em_store(i, j, k + 1) - DEy_em_store(i, j, k)) / dz;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
			dMx1 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx1 = pt_mag->dmx_glb_rk1(idx1, idy1, idz1) * Ms1 \
					+ pt_mag->dmx_AFM1_glb_rk1(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk1(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMx1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
			dMx2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx2 = pt_mag->dmx_glb_rk1(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmx_AFM1_glb_rk1(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk1(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMx2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}
		dDHx_em_rk1(i, j, k) = -pt_glb->dt / mu0 * (dEzdy - dEydz) - (dMx1 + dMx2) / M_count; //RK
		//if (abs(dDHx_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHx_em(i, j, k) = 0.;
		//}
	}

	//------------------END X component-------------------//

	//------------------Y component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdz,dEzdx,i,j,k,M_count) default(present) async(5)
	for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		idx1 = i; idx2 = i;
		if (j == 0) {
			if (pt_geo->periodicY == true) {
				idy1 = ny - 1;
			}
			else {
				idy1 = 0;
			}
		}
		else {
			idy1 = j - 1;
		}
		if (j == ny) {
			if (pt_geo->periodicY == true) {
				idy2 = 0;
			}
			else {
				idy2 = ny - 1;
			}
		}
		else {
			idy2 = j;
		}
		idz1 = k; idz2 = k;

		dExdz = (DEx_em_store(i, j, k + 1) - DEx_em_store(i, j, k)) / dz;
		dEzdx = (DEz_em_store(i + 1, j, k) - DEz_em_store(i, j, k)) / dx;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMy1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy1 = pt_mag->dmy_glb_rk1(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmy_AFM1_glb_rk1(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk1(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMy1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMy2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy2 = pt_mag->dmy_glb_rk1(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmy_AFM1_glb_rk1(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk1(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMy2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}
		dDHy_em_rk1(i, j, k) = -pt_glb->dt / mu0 * (dExdz - dEzdx) - (dMy1 + dMy2) / M_count; //RK
		//if (abs(dDHy_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHy_em(i, j, k) = 0.;
		//}
	}
	//------------------END Y component-------------------//

	//------------------Z component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdy,dEydx,i,j,k,M_count) default(present) async(6)
	for (long int id = 0; id < nx * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		dEydx = (DEy_em_store(i + 1, j, k) - DEy_em_store(i, j, k)) / dx;
		dExdy = (DEx_em_store(i, j + 1, k) - DEx_em_store(i, j, k)) / dy;

		idx1 = i; idx2 = i;
		idy1 = j; idy2 = j;

		M_count = 0.;
		//----------Magnetization 1------//
		if (k == 0) {
			if (pt_geo->periodicZ == true) {
				idz1 = nz - 1;
			}
			else {
				idz1 = 0;
			}
		}
		else {
			idz1 = k - 1;
		}

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMz1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz1 = pt_mag->dmz_glb_rk1(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmz_AFM1_glb_rk1(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk1(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMz1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//
		if (k == nz) {
			if (pt_geo->periodicZ == true) {
				idz2 = 0;
			}
			else {
				idz2 = nz - 1;
			}
		}
		else {
			idz2 = k;
		}
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMz2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz2 = pt_mag->dmz_glb_rk1(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmz_AFM1_glb_rk1(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk1(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMz2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}
		dDHz_em_rk1(i, j, k) = -pt_glb->dt / mu0 * (dEydx - dExdy) - (dMz1 + dMz2) / M_count; //RK
		//if (abs(dDHz_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHz_em(i, j, k) = 0.;
		//}
	}
	//------------------END Z component-------------------//
}

void EMdynamic_system::get_dH_RK2() {
	double dExdz, dExdy, dEydx, dEydz, dEzdx, dEzdy;
	double dMx1, dMy1, dMz1;
	double dMx2, dMy2, dMz2;
	double Ms1, Ms2;
	double Ms1_AFM1, Ms2_AFM1;
	double Ms1_AFM2, Ms2_AFM2;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idy1, idz1;
	long int idx2, idy2, idz2;
	double M_count;

	//------------------X component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dEzdy,dEydz,i,j,k,M_count) default(present) async(4)
	for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
		i = id / ((ny) * (nz));
		j = (id - i * ((ny) * (nz))) / (nz);
		k = id - i * ((ny) * (nz)) - j * (nz);

		if (i == 0) {
			if (pt_geo->periodicX == true) {
				idx1 = nx - 1;
			}
			else {
				idx1 = 0;
			}
		}
		else {
			idx1 = i - 1;
		}
		if (i == nx) {
			if (pt_geo->periodicX == true) {
				idx2 = 0;
			}
			else {
				idx2 = nx - 1;
			}
		}
		else {
			idx2 = i;
		}

		idy1 = j; idy2 = j;
		idz1 = k; idz2 = k;
		dEzdy = (DEz_em_store(i, j + 1, k) - DEz_em_store(i, j, k)) / dy;
		dEydz = (DEy_em_store(i, j, k + 1) - DEy_em_store(i, j, k)) / dz;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
			dMx1 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx1 = pt_mag->dmx_glb_rk2(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmx_AFM1_glb_rk2(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk2(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMx1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
			dMx2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx2 = pt_mag->dmx_glb_rk2(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmx_AFM1_glb_rk2(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk2(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMx2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}
		dDHx_em_rk2(i, j, k) = -pt_glb->dt / mu0 * (dEzdy - dEydz) - (dMx1 + dMx2) / M_count; //RK
		//if (abs(dDHx_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHx_em(i, j, k) = 0.;
		//}
	}

	//------------------END X component-------------------//

	//------------------Y component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdz,dEzdx,i,j,k,M_count) default(present) async(5)
	for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		idx1 = i; idx2 = i;
		if (j == 0) {
			if (pt_geo->periodicY == true) {
				idy1 = ny - 1;
			}
			else {
				idy1 = 0;
			}
		}
		else {
			idy1 = j - 1;
		}
		if (j == ny) {
			if (pt_geo->periodicY == true) {
				idy2 = 0;
			}
			else {
				idy2 = ny - 1;
			}
		}
		else {
			idy2 = j;
		}
		idz1 = k; idz2 = k;

		dExdz = (DEx_em_store(i, j, k + 1) - DEx_em_store(i, j, k)) / dz;
		dEzdx = (DEz_em_store(i + 1, j, k) - DEz_em_store(i, j, k)) / dx;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMy1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy1 = pt_mag->dmy_glb_rk2(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmy_AFM1_glb_rk2(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk2(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMy1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMy2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy2 = pt_mag->dmy_glb_rk2(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmy_AFM1_glb_rk2(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk2(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMy2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}
		dDHy_em_rk2(i, j, k) = -pt_glb->dt / mu0 * (dExdz - dEzdx) - (dMy1 + dMy2) / M_count; //RK
		//if (abs(dDHy_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHy_em(i, j, k) = 0.;
		//}
	}
	//------------------END Y component-------------------//

	//------------------Z component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdy,dEydx,i,j,k,M_count) default(present) async(6)
	for (long int id = 0; id < nx * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		dEydx = (DEy_em_store(i + 1, j, k) - DEy_em_store(i, j, k)) / dx;
		dExdy = (DEx_em_store(i, j + 1, k) - DEx_em_store(i, j, k)) / dy;

		idx1 = i; idx2 = i;
		idy1 = j; idy2 = j;

		M_count = 0.;
		//----------Magnetization 1------//
		if (k == 0) {
			if (pt_geo->periodicZ == true) {
				idz1 = nz - 1;
			}
			else {
				idz1 = 0;
			}
		}
		else {
			idz1 = k - 1;
		}
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMz1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz1 = pt_mag->dmz_glb_rk2(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmz_AFM1_glb_rk2(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk2(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMz1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//
		if (k == nz) {
			if (pt_geo->periodicZ == true) {
				idz2 = 0;
			}
			else {
				idz2 = nz - 1;
			}
		}
		else {
			idz2 = k;
		}
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMz2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz2 = pt_mag->dmz_glb_rk2(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmz_AFM1_glb_rk2(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk2(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMz2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}
		dDHz_em_rk2(i, j, k) = -pt_glb->dt / mu0 * (dEydx - dExdy) - (dMz1 + dMz2) / M_count; //RK
		//if (abs(dDHz_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHz_em(i, j, k) = 0.;
		//}
	}
	//------------------END Z component-------------------//
}

void EMdynamic_system::get_dH_RK3() {
	double dExdz, dExdy, dEydx, dEydz, dEzdx, dEzdy;
	double dMx1, dMy1, dMz1;
	double dMx2, dMy2, dMz2;
	double Ms1, Ms2;
	double Ms1_AFM1, Ms2_AFM1;
	double Ms1_AFM2, Ms2_AFM2;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idy1, idz1;
	long int idx2, idy2, idz2;
	double M_count;

	//------------------X component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dEzdy,dEydz,i,j,k,M_count) default(present) async(4)
	for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
		i = id / ((ny) * (nz));
		j = (id - i * ((ny) * (nz))) / (nz);
		k = id - i * ((ny) * (nz)) - j * (nz);

		if (i == 0) {
			if (pt_geo->periodicX == true) {
				idx1 = nx - 1;
			}
			else {
				idx1 = 0;
			}
		}
		else {
			idx1 = i - 1;
		}
		if (i == nx) {
			if (pt_geo->periodicX == true) {
				idx2 = 0;
			}
			else {
				idx2 = nx - 1;
			}
		}
		else {
			idx2 = i;
		}

		idy1 = j; idy2 = j;
		idz1 = k; idz2 = k;
		dEzdy = (DEz_em_store(i, j + 1, k) - DEz_em_store(i, j, k)) / dy;
		dEydz = (DEy_em_store(i, j, k + 1) - DEy_em_store(i, j, k)) / dz;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
			dMx1 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx1 = pt_mag->dmx_glb_rk3(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmx_AFM1_glb_rk3(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk3(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMx1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
			dMx2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx2 = pt_mag->dmx_glb_rk3(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmx_AFM1_glb_rk3(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk3(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMx2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}
		dDHx_em_rk3(i, j, k) = -pt_glb->dt / mu0 * (dEzdy - dEydz) - (dMx1 + dMx2) / M_count; //RK
		//if (abs(dDHx_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHx_em(i, j, k) = 0.;
		//}
	}

	//------------------END X component-------------------//

	//------------------Y component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdz,dEzdx,i,j,k,M_count) default(present) async(5)
	for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		idx1 = i; idx2 = i;
		if (j == 0) {
			if (pt_geo->periodicY == true) {
				idy1 = ny - 1;
			}
			else {
				idy1 = 0;
			}
		}
		else {
			idy1 = j - 1;
		}
		if (j == ny) {
			if (pt_geo->periodicY == true) {
				idy2 = 0;
			}
			else {
				idy2 = ny - 1;
			}
		}
		else {
			idy2 = j;
		}
		idz1 = k; idz2 = k;

		dExdz = (DEx_em_store(i, j, k + 1) - DEx_em_store(i, j, k)) / dz;
		dEzdx = (DEz_em_store(i + 1, j, k) - DEz_em_store(i, j, k)) / dx;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMy1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy1 = pt_mag->dmy_glb_rk3(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmy_AFM1_glb_rk3(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk3(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMy1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMy2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy2 = pt_mag->dmy_glb_rk3(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmy_AFM1_glb_rk3(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk3(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMy2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}
		dDHy_em_rk3(i, j, k) = -pt_glb->dt / mu0 * (dExdz - dEzdx) - (dMy1 + dMy2) / M_count; //RK
		//if (abs(dDHy_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHy_em(i, j, k) = 0.;
		//}
	}
	//------------------END Y component-------------------//

	//------------------Z component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdy,dEydx,i,j,k,M_count) default(present) async(6)
	for (long int id = 0; id < nx * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		dEydx = (DEy_em_store(i + 1, j, k) - DEy_em_store(i, j, k)) / dx;
		dExdy = (DEx_em_store(i, j + 1, k) - DEx_em_store(i, j, k)) / dy;

		idx1 = i; idx2 = i;
		idy1 = j; idy2 = j;

		M_count = 0.;
		//----------Magnetization 1------//
		if (k == 0) {
			if (pt_geo->periodicZ == true) {
				idz1 = nz - 1;
			}
			else {
				idz1 = 0;
			}
		}
		else {
			idz1 = k - 1;
		}
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMz1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz1 = pt_mag->dmz_glb_rk3(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmz_AFM1_glb_rk3(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk3(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMz1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//
		if (k == nz) {
			if (pt_geo->periodicZ == true) {
				idz2 = 0;
			}
			else {
				idz2 = nz - 1;
			}
		}
		else {
			idz2 = k;
		}
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMz2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz2 = pt_mag->dmz_glb_rk3(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmz_AFM1_glb_rk3(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk3(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMz2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}
		dDHz_em_rk3(i, j, k) = -pt_glb->dt / mu0 * (dEydx - dExdy) - (dMz1 + dMz2) / M_count; //RK
		//if (abs(dDHz_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHz_em(i, j, k) = 0.;
		//}
	}
	//------------------END Z component-------------------//
}

void EMdynamic_system::get_dH_RK4() {
	double dExdz, dExdy, dEydx, dEydz, dEzdx, dEzdy;
	double dMx1, dMy1, dMz1;
	double dMx2, dMy2, dMz2;
	double Ms1, Ms2;
	double Ms1_AFM1, Ms2_AFM1;
	double Ms1_AFM2, Ms2_AFM2;

	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idy1, idz1;
	long int idx2, idy2, idz2;
	double M_count;

	//------------------X component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dEzdy,dEydz,i,j,k,M_count) default(present) async(4)
	for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
		i = id / ((ny) * (nz));
		j = (id - i * ((ny) * (nz))) / (nz);
		k = id - i * ((ny) * (nz)) - j * (nz);

		if (i == 0) {
			if (pt_geo->periodicX == true) {
				idx1 = nx - 1;
			}
			else {
				idx1 = 0;
			}
		}
		else {
			idx1 = i - 1;
		}
		if (i == nx) {
			if (pt_geo->periodicX == true) {
				idx2 = 0;
			}
			else {
				idx2 = nx - 1;
			}
		}
		else {
			idx2 = i;
		}

		idy1 = j; idy2 = j;
		idz1 = k; idz2 = k;
		dEzdy = (DEz_em_store(i, j + 1, k) - DEz_em_store(i, j, k)) / dy;
		dEydz = (DEy_em_store(i, j, k + 1) - DEy_em_store(i, j, k)) / dz;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
			dMx1 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx1 = pt_mag->dmx_glb_rk4(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmx_AFM1_glb_rk4(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk4(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMx1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
			dMx2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMx2 = pt_mag->dmx_glb_rk4(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmx_AFM1_glb_rk4(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmx_AFM2_glb_rk4(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMx2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}
		dDHx_em_rk4(i, j, k) = -pt_glb->dt / mu0 * (dEzdy - dEydz) - (dMx1 + dMx2) / M_count; //RK
		//if (abs(dDHx_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHx_em(i, j, k) = 0.;
		//}
	}

	//------------------END X component-------------------//

	//------------------Y component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdz,dEzdx,i,j,k,M_count) default(present) async(5)
	for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
		i = id / ((ny + 1) * (nz));
		j = (id - i * ((ny + 1) * (nz))) / (nz);
		k = id - i * ((ny + 1) * (nz)) - j * (nz);

		idx1 = i; idx2 = i;
		if (j == 0) {
			if (pt_geo->periodicY == true) {
				idy1 = ny - 1;
			}
			else {
				idy1 = 0;
			}
		}
		else {
			idy1 = j - 1;
		}
		if (j == ny) {
			if (pt_geo->periodicY == true) {
				idy2 = 0;
			}
			else {
				idy2 = ny - 1;
			}
		}
		else {
			idy2 = j;
		}
		idz1 = k; idz2 = k;

		dExdz = (DEx_em_store(i, j, k + 1) - DEx_em_store(i, j, k)) / dz;
		dEzdx = (DEz_em_store(i + 1, j, k) - DEz_em_store(i, j, k)) / dx;

		M_count = 0.;
		//----------Magnetization 1------//

		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMy1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy1 = pt_mag->dmy_glb_rk4(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmy_AFM1_glb_rk4(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk4(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMy1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//

		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMy2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMy2 = pt_mag->dmy_glb_rk4(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmy_AFM1_glb_rk4(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmy_AFM2_glb_rk4(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMy2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}
		dDHy_em_rk4(i, j, k) = -pt_glb->dt / mu0 * (dExdz - dEzdx) - (dMy1 + dMy2) / M_count; //RK
		//if (abs(dDHy_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHy_em(i, j, k) = 0.;
		//}
	}
	//------------------END Y component-------------------//

	//------------------Z component-------------------//
#pragma acc parallel loop gang vector private(dMx1, dMy1, dMz1,dMx2, dMy2, dMz2,\
Ms1_AFM1, Ms2_AFM1,Ms1_AFM2, Ms2_AFM2,\
Ms1, Ms2,mat_type,idx1, idy1, idz1,idx2, idy2, idz2,mat,dExdy,dEydx,i,j,k,M_count) default(present) async(6)
	for (long int id = 0; id < nx * ny * (nz + 1); id++) {
		i = id / ((ny) * (nz + 1));
		j = (id - i * ((ny) * (nz + 1))) / (nz + 1);
		k = id - i * ((ny) * (nz + 1)) - j * (nz + 1);

		dEydx = (DEy_em_store(i + 1, j, k) - DEy_em_store(i, j, k)) / dx;
		dExdy = (DEx_em_store(i, j + 1, k) - DEx_em_store(i, j, k)) / dy;

		idx1 = i; idx2 = i;
		idy1 = j; idy2 = j;

		M_count = 0.;
		//----------Magnetization 1------//
		if (k == 0) {
			if (pt_geo->periodicZ == true) {
				idz1 = nz - 1;
			}
			else {
				idz1 = 0;
			}
		}
		else {
			idz1 = k - 1;
		}
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);

		if (mat_type == 0) {
			dMz1 = 0.;
			Ms1 = 0.; Ms1_AFM1 = 0.; Ms1_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms1 = mat->Ms; Ms1_AFM1 = mat->Ms_AFM1; Ms1_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz1 = pt_mag->dmz_glb_rk4(idx1, idy1, idz1) * Ms1\
					+ pt_mag->dmz_AFM1_glb_rk4(idx1, idy1, idz1) * Ms1_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk4(idx1, idy1, idz1) * Ms1_AFM2; //RK
			}
			else {
				dMz1 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//----------Magnetization 2------//
		if (k == nz) {
			if (pt_geo->periodicZ == true) {
				idz2 = 0;
			}
			else {
				idz2 = nz - 1;
			}
		}
		else {
			idz2 = k;
		}
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);

		if (mat_type == 0) {
			dMz2 = 0.;
			Ms2 = 0.; Ms2_AFM1 = 0.; Ms2_AFM2 = 0.;
		}
		else {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			Ms2 = mat->Ms; Ms2_AFM1 = mat->Ms_AFM1; Ms2_AFM2 = mat->Ms_AFM2;
			if (pt_glb->if_EM_fromM == true) {
				dMz2 = pt_mag->dmz_glb_rk4(idx2, idy2, idz2) * Ms2\
					+ pt_mag->dmz_AFM1_glb_rk4(idx2, idy2, idz2) * Ms2_AFM1\
					+ pt_mag->dmz_AFM2_glb_rk4(idx2, idy2, idz2) * Ms2_AFM2; //RK
			}
			else {
				dMz2 = 0.;
			}

			if (mat->if_FM == true || mat->if_AFM == true) {
				M_count = M_count + 1.;
			}
		}
		//---------END calculating temporal change of M--------//

		if (M_count < 0.5) {
			M_count = 1.;
		}
		dDHz_em_rk4(i, j, k) = -pt_glb->dt / mu0 * (dEydx - dExdy) - (dMz1 + dMz2) / M_count; //RK
		//if (abs(dDHz_em(i, j, k)) < 1.e-6 * (Ms1 + Ms2) * 0.5 * pt_glb->dt) {
		//	dDHz_em(i, j, k) = 0.;
		//}
	}
	//------------------END Z component-------------------//
}


void EMdynamic_system::update_DH_RK1()
{
#pragma acc parallel default(present) async(10)
	{
#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * nz; id++) {
			DHx_em_store(id) = DHx_em(id) + dDHx_em_rk1(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * nz; id++) {
			DHy_em_store(id) = DHy_em(id) + dDHy_em_rk1(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * ny * (nz + 1); id++) {
			DHz_em_store(id) = DHz_em(id) + dDHz_em_rk1(id) * 0.5;
		}
	}

#pragma acc wait(3,10) async(11)
#pragma acc parallel default(present) async(11)
	{
		update_DH_cell();
	}
}

void EMdynamic_system::update_DH_RK2()
{
#pragma acc parallel default(present) async(10)
	{
#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * nz; id++) {
			DHx_em_store(id) = DHx_em(id) + dDHx_em_rk2(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * nz; id++) {
			DHy_em_store(id) = DHy_em(id) + dDHy_em_rk2(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * ny * (nz + 1); id++) {
			DHz_em_store(id) = DHz_em(id) + dDHz_em_rk2(id) * 0.5;
		}
	}

#pragma acc wait(3,10) async(11)
#pragma acc parallel default(present) async(11)
	{
		update_DH_cell();
	}
}

void EMdynamic_system::update_DH_RK3()
{
#pragma acc parallel default(present) async(10)
	{
#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * nz; id++) {
			DHx_em_store(id) = DHx_em(id) + dDHx_em_rk3(id);
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * nz; id++) {
			DHy_em_store(id) = DHy_em(id) + dDHy_em_rk3(id);
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * ny * (nz + 1); id++) {
			DHz_em_store(id) = DHz_em(id) + dDHz_em_rk3(id);
		}
	}

#pragma acc wait(3,10) async(11)
#pragma acc parallel default(present) async(11)
	{
		update_DH_cell();
	}
}

void EMdynamic_system::update_DH()
{
#pragma acc parallel default(present) async(10)
	{
#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * nz; id++) {
			DHx_em(id) = DHx_em(id) + dDHx_em_rk1(id) / 6. + dDHx_em_rk2(id) / 3. + dDHx_em_rk3(id) / 3. + dDHx_em_rk4(id) / 6.;
			DHx_em_store(id) = DHx_em(id);
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * nz; id++) {
			DHy_em(id) = DHy_em(id) + dDHy_em_rk1(id) / 6. + dDHy_em_rk2(id) / 3. + dDHy_em_rk3(id) / 3. + dDHy_em_rk4(id) / 6.;
			DHy_em_store(id) = DHy_em(id);
		}

#pragma acc loop gang vector
		for (long id = 0; id < nx * ny * (nz + 1); id++) {
			DHz_em(id) = DHz_em(id) + dDHz_em_rk1(id) / 6. + dDHz_em_rk2(id) / 3. + dDHz_em_rk3(id) / 3. + dDHz_em_rk4(id) / 6.;
			DHz_em_store(id) = DHz_em(id);
		}
	}

#pragma acc wait(3,10) async(11)
#pragma acc parallel default(present) async(11)
	{
		update_DH_cell();
	}
}

void EMdynamic_system::update_DE_RK1()
{
#pragma acc parallel default(present) async(8)
	{
#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
			DEx_em_store(id) = DEx_em(id) + dDEx_em_rk1(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
			DEy_em_store(id) = DEy_em(id) + dDEy_em_rk1(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
			DEz_em_store(id) = DEz_em(id) + dDEz_em_rk1(id) * 0.5;
		}
	}

	if (pt_glb->if_input_planeEM_E == true) {
#pragma acc parallel default(present) async(8)
		{
			update_planeEM_half();
		}
	}

	if (pt_glb->if_periodic_allsurface == false) {
		update_DE_Boundary_half();
	}
	if (pt_glb->if_PEC_all == true) {
		update_DE_PEC();
	}

#pragma acc parallel default(present) async(8)
	{
		update_DE_cell();
	}

#pragma acc parallel default(present) async(9)
	{
		update_Jp_RK1();
	}
}

void EMdynamic_system::update_DE_RK2()
{
#pragma acc parallel default(present) async(8)
	{
#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
			DEx_em_store(id) = DEx_em(id) + dDEx_em_rk2(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
			DEy_em_store(id) = DEy_em(id) + dDEy_em_rk2(id) * 0.5;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
			DEz_em_store(id) = DEz_em(id) + dDEz_em_rk2(id) * 0.5;
		}
	}

	if (pt_glb->if_input_planeEM_E == true) {
#pragma acc parallel default(present) async(8)
		{
			update_planeEM_half();
		}
	}

	if (pt_glb->if_periodic_allsurface == false) {
		update_DE_Boundary_half();
	}
	if (pt_glb->if_PEC_all == true) {
		update_DE_PEC();
	}

#pragma acc parallel default(present) async(8)
	{
		update_DE_cell();
	}

#pragma acc parallel default(present) async(9)
	{
		update_Jp_RK2();
	}
}

void EMdynamic_system::update_DE_RK3()
{
#pragma acc parallel default(present) async(8)
	{
#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
			DEx_em_store(id) = DEx_em(id) + dDEx_em_rk3(id);
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
			DEy_em_store(id) = DEy_em(id) + dDEy_em_rk3(id);
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
			DEz_em_store(id) = DEz_em(id) + dDEz_em_rk3(id);
		}
	}

	if (pt_glb->if_input_planeEM_E == true) {
#pragma acc parallel default(present) async(8)
		{
			update_planeEM_full();
		}
	}

	if (pt_glb->if_periodic_allsurface == false) {
		update_DE_Boundary_full();
	}
	if (pt_glb->if_PEC_all == true) {
		update_DE_PEC();
	}

#pragma acc parallel default(present) async(8)
	{
		update_DE_cell();
	}

#pragma acc parallel default(present) async(9)
	{
		update_Jp_RK3();
	}
}

void EMdynamic_system::update_DE()
{
	if (pt_glb->if_periodic_allsurface == false) {
		transfer_pointer();
	}

#pragma acc parallel default(present) async(8)
	{
#pragma acc loop gang vector
		for (long id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
			DEx_em(id) = DEx_em(id) + dDEx_em_rk1(id) / 6. + dDEx_em_rk2(id) / 3. + dDEx_em_rk3(id) / 3. + dDEx_em_rk4(id) / 6.;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
			DEy_em(id) = DEy_em(id) + dDEy_em_rk1(id) / 6. + dDEy_em_rk2(id) / 3. + dDEy_em_rk3(id) / 3. + dDEy_em_rk4(id) / 6.;
		}

#pragma acc loop gang vector
		for (long id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
			DEz_em(id) = DEz_em(id) + dDEz_em_rk1(id) / 6. + dDEz_em_rk2(id) / 3. + dDEz_em_rk3(id) / 3. + dDEz_em_rk4(id) / 6.;
		}
	}

	if (pt_glb->if_input_planeEM_E == true) {
#pragma acc parallel default(present) async(8)
		{
			update_planeEM();
		}
	}

	if (pt_glb->if_periodic_allsurface == false) {
		update_DE_Boundary();
	}
	if (pt_glb->if_PEC_all == true) {
		update_DE_PEC();
	}

#pragma acc parallel default(present) async(8)
	{
		DEx_em_store = DEx_em;
		DEy_em_store = DEy_em;
		DEz_em_store = DEz_em;
	}

#pragma acc parallel default(present) async(8)
	{
		update_DE_cell();
	}

#pragma acc parallel default(present) async(9)
	{
		update_Jp();
	}
}

#pragma acc routine gang// nohost
void EMdynamic_system::update_DH_cell() {
	long int i, j, k;
	unsigned int mat_type;
	material* mat;

#pragma acc loop gang vector private(i,j,k,mat_type,mat)
	for (long int id = 0; id < n; id++) {
		i = id / (ny * nz);
		j = (id - i * (ny * nz)) / nz;
		k = id - i * (ny * nz) - j * nz;

		DHx_em_cell(i, j, k) = (DHx_em_store(i, j, k) + DHx_em_store(i + 1, j, k)) / 2.;
		DHy_em_cell(i, j, k) = (DHy_em_store(i, j, k) + DHy_em_store(i, j + 1, k)) / 2.;
		if (pt_glb->if_mag_1Dmodel == true && pt_glb->if_EM_fromM == true) {
			mat_type = pt_glb->material_cell(id);
			if (mat_type == 0) {
				DHz_em_cell(id) = 0.;
			}
			else {
				mat = &(pt_glb->material_parameters[mat_type - 1]);
				DHz_em_cell(id) = -1. * mat->Ms * pt_mag->mz_glb_store(id) \
					- 1. * mat->Ms_AFM1 * pt_mag->mz_AFM1_glb_store(id) \
					- 1. * mat->Ms_AFM2 * pt_mag->mz_AFM2_glb_store(id) \
					- pt_mag->Hz_stat(id);
			}
		}
		else {
			DHz_em_cell(i, j, k) = (DHz_em_store(i, j, k) + DHz_em_store(i, j, k + 1)) / 2.;
		}
	}
}

#pragma acc routine gang// nohost
void EMdynamic_system::update_DE_cell() {
	long int i, j, k;

#pragma acc loop gang vector private(i,j,k)
	for (long int id = 0; id < n; id++) {
		i = id / (ny * nz);
		j = (id - i * (ny * nz)) / nz;
		k = id - i * (ny * nz) - j * nz;

		DEx_em_cell(i, j, k) = (DEx_em_store(i, j, k) + DEx_em_store(i, j + 1, k) + DEx_em_store(i, j, k + 1) + DEx_em_store(i, j + 1, k + 1)) / 4.;
		DEy_em_cell(i, j, k) = (DEy_em_store(i, j, k) + DEy_em_store(i + 1, j, k) + DEy_em_store(i, j, k + 1) + DEy_em_store(i + 1, j, k + 1)) / 4.;
		DEz_em_cell(i, j, k) = (DEz_em_store(i, j, k) + DEz_em_store(i + 1, j, k) + DEz_em_store(i, j + 1, k) + DEz_em_store(i + 1, j + 1, k)) / 4.;
	}
}

#pragma acc routine gang// nohost
void EMdynamic_system::update_planeEM_half() {
	unsigned long i, j;

#pragma acc loop gang vector private(i,j)
	for (unsigned long id = 0; id < nx * ny; id++) {
		i = id / ny;
		j = id - i * ny;
		DEx_em_store(i, j, planeEM_source_z) = pt_glb->planeEM_E[0] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[0] * (pt_glb->time_device - 0.5 * pt_glb->dt));

		DEy_em_store(i, j, planeEM_source_z) = pt_glb->planeEM_E[1] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[1] * (pt_glb->time_device - 0.5 * pt_glb->dt));

		DEz_em_store(i, j, planeEM_source_z) = pt_glb->planeEM_E[2] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[2] * (pt_glb->time_device - 0.5 * pt_glb->dt));
	}
}

#pragma acc routine gang// nohost
void EMdynamic_system::update_planeEM_full() {
	unsigned long i, j;

#pragma acc loop gang vector private(i,j)
	for (unsigned long id = 0; id < nx * ny; id++) {
		i = id / ny;
		j = id - i * ny;
		DEx_em_store(i, j, planeEM_source_z) = pt_glb->planeEM_E[0] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[0] * pt_glb->time_device);

		DEy_em_store(i, j, planeEM_source_z) = pt_glb->planeEM_E[1] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[1] * pt_glb->time_device);

		DEz_em_store(i, j, planeEM_source_z) = pt_glb->planeEM_E[2] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[2] * pt_glb->time_device);
	}
}

#pragma acc routine gang// nohost
void EMdynamic_system::update_planeEM() {
	unsigned long i, j;

#pragma acc loop gang vector private(i,j)
	for (unsigned long id = 0; id < nx * ny; id++) {
		i = id / ny;
		j = id - i * ny;
		DEx_em(i, j, planeEM_source_z) = pt_glb->planeEM_E[0] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[0] * pt_glb->time_device);

		DEy_em(i, j, planeEM_source_z) = pt_glb->planeEM_E[1] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[1] * pt_glb->time_device);

		DEz_em(i, j, planeEM_source_z) = pt_glb->planeEM_E[2] * \
			sin(2. * PI * pt_glb->planeEM_E_freq[2] * pt_glb->time_device);
	}
}

void EMdynamic_system::update_Jf_input() {
	unsigned long int i, j, k;
	double temporal_var, spatial_var;

#pragma acc parallel default(present) async(7)
	{
#pragma acc loop gang vector private(i,j,k, temporal_var, spatial_var)
		for (unsigned long id = 0; id < Jf_n; id++) {
			i = id / (Jf_ny * Jf_nz);
			j = (id - i * (Jf_ny * Jf_nz)) / Jf_nz;
			k = id - i * (Jf_ny * Jf_nz) - j * Jf_nz;

			i = i + pt_glb->Jfin_xi;
			j = j + pt_glb->Jfin_yi;
			k = k + pt_glb->Jfin_zi;

			if (pt_glb->Jf_input_type < 3) {
				temporal_var = sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device - pt_glb->dt));
				spatial_var = PI / (static_cast<double>(Jf_nz) * dz) * (static_cast<double>(k - pt_glb->Jfin_zi) + 0.5) * dz;
				if (pt_glb->Jf_input_type == 0) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
				}
				else if (pt_glb->Jf_input_type == 1) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
				}
				else if (pt_glb->Jf_input_type == 2) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
				}
			}
			else if (pt_glb->Jf_input_type == 3) {
				if ((pt_glb->time_device - pt_glb->dt) <= static_cast<double>(pt_glb->num_Jf_cycles) / pt_glb->Jf_input_freq) {
					temporal_var = sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device - pt_glb->dt));
				}
				else {
					temporal_var = 0.;
				}

				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
			}
			else if (pt_glb->Jf_input_type == 4) {
				temporal_var = pt_glb->time_device - pt_glb->dt - 5. * pt_glb->Jf_input_freq;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
			}
			else if (pt_glb->Jf_input_type == 5) {
				temporal_var = pt_glb->time_device - pt_glb->dt - 5. * pt_glb->Jf_input_freq;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
			}
			else if (pt_glb->Jf_input_type == 6) {
				temporal_var = pt_glb->time_device - pt_glb->dt - 5. * pt_glb->Jf_input_sigma;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_sigma, 2.))) \
						* sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device - pt_glb->dt));
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_sigma, 2.))) \
						* sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device - pt_glb->dt));
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_sigma, 2.))) \
						* sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device - pt_glb->dt));
				}
			}

		}
	}
	//	{
	//#pragma acc loop gang vector private(i,j,temporal_var, spatial_var)
	//		for (long int k = Jf_z1; k < Jf_z2; k++) {
	//			temporal_var = sin(2. * PI * pt_glb->Jf_test_freq * (pt_glb->time_device - pt_glb->dt));
	//			spatial_var = PI / Jf_thickness * (static_cast<double>(k - Jf_z1) + 0.5) * dz;
	//#pragma acc loop seq
	//			for (i = 0; i < nx; i++) {
	//#pragma acc loop seq
	//				for (j = 0; j < ny; j++) {
	//					if (pt_glb->material_cell(i, j, k) != 0) {
	//						if (pt_glb->Jf_test_type == 0) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var;
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var;
	//						}
	//						else if (pt_glb->Jf_test_type == 1) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var * sin(spatial_var);
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var * sin(spatial_var);
	//						}
	//						else if (pt_glb->Jf_test_type == 2) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var * cos(spatial_var);
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var * cos(spatial_var);
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
}

void EMdynamic_system::update_Jf_input_half() {
	unsigned long int i, j, k;
	double temporal_var, spatial_var;

#pragma acc parallel default(present) async(7)
	{
#pragma acc loop gang vector private(i,j,k, temporal_var, spatial_var)
		for (unsigned long id = 0; id < Jf_n; id++) {
			i = id / (Jf_ny * Jf_nz);
			j = (id - i * (Jf_ny * Jf_nz)) / Jf_nz;
			k = id - i * (Jf_ny * Jf_nz) - j * Jf_nz;

			i = i + pt_glb->Jfin_xi;
			j = j + pt_glb->Jfin_yi;
			k = k + pt_glb->Jfin_zi;

			if (pt_glb->Jf_input_type < 3) {
				temporal_var = sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device - 0.5 * pt_glb->dt));
				spatial_var = PI / (static_cast<double>(Jf_nz) * dz) * (static_cast<double>(k - pt_glb->Jfin_zi) + 0.5) * dz;
				if (pt_glb->Jf_input_type == 0) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
				}
				else if (pt_glb->Jf_input_type == 1) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
				}
				else if (pt_glb->Jf_input_type == 2) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
				}
			}
			else if (pt_glb->Jf_input_type == 3) {
				if ((pt_glb->time_device - 0.5 * pt_glb->dt) <= static_cast<double>(pt_glb->num_Jf_cycles) / pt_glb->Jf_input_freq) {
					temporal_var = sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device - 0.5 * pt_glb->dt));
				}
				else {
					temporal_var = 0.;
				}

				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
			}
			else if (pt_glb->Jf_input_type == 4) {
				temporal_var = pt_glb->time_device - 0.5 * pt_glb->dt - 5. * pt_glb->Jf_input_freq;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
			}
			else if (pt_glb->Jf_input_type == 5) {
				temporal_var = pt_glb->time_device - 0.5 * pt_glb->dt - 5. * pt_glb->Jf_input_freq;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
			}
			else if (pt_glb->Jf_input_type == 6) {
				temporal_var = pt_glb->time_device - 0.5 * pt_glb->dt - 5. * pt_glb->Jf_input_sigma;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_sigma, 2.))) \
						* sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device - 0.5 * pt_glb->dt));
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_sigma, 2.))) \
						* sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device - 0.5 * pt_glb->dt));
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_sigma, 2.))) \
						* sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device - 0.5 * pt_glb->dt));
				}
			}

		}
	}

	//#pragma acc parallel default(present) async(7)
	//	{
	//#pragma acc loop gang vector private(i,j,temporal_var, spatial_var)
	//		for (long int k = Jf_z1; k < Jf_z2; k++) {
	//			temporal_var = sin(2. * PI * pt_glb->Jf_test_freq * (pt_glb->time_device - 0.5 * pt_glb->dt));
	//			spatial_var = PI / Jf_thickness * (static_cast<double>(k - Jf_z1) + 0.5) * dz;
	//#pragma acc loop seq
	//			for (i = 0; i < nx; i++) {
	//#pragma acc loop seq
	//				for (j = 0; j < ny; j++) {
	//					if (pt_glb->material_cell(i, j, k) != 0) {
	//						if (pt_glb->Jf_test_type == 0) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var;
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var;
	//						}
	//						else if (pt_glb->Jf_test_type == 1) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var * sin(spatial_var);
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var * sin(spatial_var);
	//						}
	//						else if (pt_glb->Jf_test_type == 2) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var * cos(spatial_var);
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var * cos(spatial_var);
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
}

void EMdynamic_system::update_Jf_input_full() {
	unsigned long int i, j, k;
	double temporal_var, spatial_var;

#pragma acc parallel default(present) async(7)
	{
#pragma acc loop gang vector private(i,j,k, temporal_var, spatial_var)
		for (unsigned long id = 0; id < Jf_n; id++) {
			i = id / (Jf_ny * Jf_nz);
			j = (id - i * (Jf_ny * Jf_nz)) / Jf_nz;
			k = id - i * (Jf_ny * Jf_nz) - j * Jf_nz;

			i = i + pt_glb->Jfin_xi;
			j = j + pt_glb->Jfin_yi;
			k = k + pt_glb->Jfin_zi;

			if (pt_glb->Jf_input_type < 3) {
				temporal_var = sin(2. * PI * pt_glb->Jf_input_freq * pt_glb->time_device);
				spatial_var = PI / (static_cast<double>(Jf_nz) * dz) * (static_cast<double>(k - pt_glb->Jfin_zi) + 0.5) * dz;
				if (pt_glb->Jf_input_type == 0) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
					}
				}
				else if (pt_glb->Jf_input_type == 1) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var * sin(spatial_var);
					}
				}
				else if (pt_glb->Jf_input_type == 2) {
					if (pt_glb->Jf_input_component == 'x') {
						DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'y') {
						DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
					else if (pt_glb->Jf_input_component == 'z') {
						DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var * cos(spatial_var);
					}
				}
			}
			else if (pt_glb->Jf_input_type == 3) {
				if ((pt_glb->time_device) <= static_cast<double>(pt_glb->num_Jf_cycles) / pt_glb->Jf_input_freq) {
					temporal_var = sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device));
				}
				else {
					temporal_var = 0.;
				}

				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * temporal_var;
				}
			}
			else if (pt_glb->Jf_input_type == 4) {
				temporal_var = pt_glb->time_device - 5. * pt_glb->Jf_input_freq;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
			}
			else if (pt_glb->Jf_input_type == 5) {
				temporal_var = pt_glb->time_device - 5. * pt_glb->Jf_input_freq;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * sqrt(exp(1.)) * temporal_var / pt_glb->Jf_input_freq * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_freq, 2.)));
				}
			}
			else if (pt_glb->Jf_input_type == 6) {
				temporal_var = pt_glb->time_device - 5. * pt_glb->Jf_input_sigma;
				if (pt_glb->Jf_input_component == 'x') {
					DJfx(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_sigma, 2.))) \
						* sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device));
				}
				else if (pt_glb->Jf_input_component == 'y') {
					DJfy(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_sigma, 2.))) \
						* sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device));
				}
				else if (pt_glb->Jf_input_component == 'z') {
					DJfz(i, j, k) = pt_glb->Jf_input_amp * exp(-1. * pow(temporal_var, 2.) / (2. * pow(pt_glb->Jf_input_sigma, 2.))) \
						* sin(2. * PI * pt_glb->Jf_input_freq * (pt_glb->time_device));
				}
			}

		}
	}
	//	unsigned int i, j;
	//	double temporal_var, spatial_var;
	//
	//#pragma acc parallel default(present) async(7)
	//	{
	//#pragma acc loop gang vector private(i,j,temporal_var, spatial_var)
	//		for (long int k = Jf_z1; k < Jf_z2; k++) {
	//			temporal_var = sin(2. * PI * pt_glb->Jf_test_freq * pt_glb->time_device);
	//			spatial_var = PI / Jf_thickness * (static_cast<double>(k - Jf_z1) + 0.5) * dz;
	//#pragma acc loop seq
	//			for (i = 0; i < nx; i++) {
	//#pragma acc loop seq
	//				for (j = 0; j < ny; j++) {
	//					if (pt_glb->material_cell(i, j, k) != 0) {
	//						if (pt_glb->Jf_test_type == 0) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var;
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var;
	//						}
	//						else if (pt_glb->Jf_test_type == 1) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var * sin(spatial_var);
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var * sin(spatial_var);
	//						}
	//						else if (pt_glb->Jf_test_type == 2) {
	//							DJfx(i, j, k) = pt_glb->Jf_test_amp * temporal_var * cos(spatial_var);
	//							DJfy(i, j, k) = pt_glb->Jf_test_amp * temporal_var * cos(spatial_var);
	//						}
	//					}
	//				}
	//			}
	//		}
	//	}
}

#pragma acc routine gang
void EMdynamic_system::get_dJp_RK1() {
	unsigned int mat_type;
	material* mat;
	double wp1, wp2, wp3;
	double re1, re2, re3;
	double cn1, cn2, cn3;

#pragma acc loop gang vector private(mat_type, mat, cn1, cn2, cn3, wp1, wp2, wp3, re1, re2, re3)
	for (long int id = 0; id < n; id++) {
		mat_type = pt_glb->material_cell(id);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);

			cn1 = mat->comp_n1; cn2 = mat->comp_n2; cn3 = mat->comp_n3;
			wp1 = mat->omega_plasma_n1; wp2 = mat->omega_plasma_n2; wp3 = mat->omega_plasma_n3;
			re1 = 1. / mat->tao_e_n1; re2 = 1. / mat->tao_e_n2; re3 = 1. / mat->tao_e_n3;

			dJpx_n1_rk1(id) = (e0 * wp1 * wp1 * DEx_em_cell(id) * cn1 - re1 * Jpx_n1_store(id)) * pt_glb->dt;
			dJpy_n1_rk1(id) = (e0 * wp1 * wp1 * DEy_em_cell(id) * cn1 - re1 * Jpy_n1_store(id)) * pt_glb->dt;
			dJpz_n1_rk1(id) = (e0 * wp1 * wp1 * DEz_em_cell(id) * cn1 - re1 * Jpz_n1_store(id)) * pt_glb->dt;

			dJpx_n2_rk1(id) = (e0 * wp2 * wp2 * DEx_em_cell(id) * cn2 - re2 * Jpx_n2_store(id)) * pt_glb->dt;
			dJpy_n2_rk1(id) = (e0 * wp2 * wp2 * DEy_em_cell(id) * cn2 - re2 * Jpy_n2_store(id)) * pt_glb->dt;
			dJpz_n2_rk1(id) = (e0 * wp2 * wp2 * DEz_em_cell(id) * cn2 - re2 * Jpz_n2_store(id)) * pt_glb->dt;

			dJpx_n3_rk1(id) = (e0 * wp3 * wp3 * DEx_em_cell(id) * cn3 - re3 * Jpx_n3_store(id)) * pt_glb->dt;
			dJpy_n3_rk1(id) = (e0 * wp3 * wp3 * DEy_em_cell(id) * cn3 - re3 * Jpy_n3_store(id)) * pt_glb->dt;
			dJpz_n3_rk1(id) = (e0 * wp3 * wp3 * DEz_em_cell(id) * cn3 - re3 * Jpz_n3_store(id)) * pt_glb->dt;
		}
	}

}

#pragma acc routine gang
void EMdynamic_system::get_dJp_RK2() {
	unsigned int mat_type;
	material* mat;
	double wp1, wp2, wp3;
	double re1, re2, re3;
	double cn1, cn2, cn3;

#pragma acc loop gang vector private(mat_type, mat, cn1, cn2, cn3, wp1, wp2, wp3, re1, re2, re3)
	for (long int id = 0; id < n; id++) {
		mat_type = pt_glb->material_cell(id);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);

			cn1 = mat->comp_n1; cn2 = mat->comp_n2; cn3 = mat->comp_n3;
			wp1 = mat->omega_plasma_n1; wp2 = mat->omega_plasma_n2; wp3 = mat->omega_plasma_n3;
			re1 = 1. / mat->tao_e_n1; re2 = 1. / mat->tao_e_n2; re3 = 1. / mat->tao_e_n3;

			dJpx_n1_rk2(id) = (e0 * wp1 * wp1 * DEx_em_cell(id) * cn1 - re1 * Jpx_n1_store(id)) * pt_glb->dt;
			dJpy_n1_rk2(id) = (e0 * wp1 * wp1 * DEy_em_cell(id) * cn1 - re1 * Jpy_n1_store(id)) * pt_glb->dt;
			dJpz_n1_rk2(id) = (e0 * wp1 * wp1 * DEz_em_cell(id) * cn1 - re1 * Jpz_n1_store(id)) * pt_glb->dt;

			dJpx_n2_rk2(id) = (e0 * wp2 * wp2 * DEx_em_cell(id) * cn2 - re2 * Jpx_n2_store(id)) * pt_glb->dt;
			dJpy_n2_rk2(id) = (e0 * wp2 * wp2 * DEy_em_cell(id) * cn2 - re2 * Jpy_n2_store(id)) * pt_glb->dt;
			dJpz_n2_rk2(id) = (e0 * wp2 * wp2 * DEz_em_cell(id) * cn2 - re2 * Jpz_n2_store(id)) * pt_glb->dt;

			dJpx_n3_rk2(id) = (e0 * wp3 * wp3 * DEx_em_cell(id) * cn3 - re3 * Jpx_n3_store(id)) * pt_glb->dt;
			dJpy_n3_rk2(id) = (e0 * wp3 * wp3 * DEy_em_cell(id) * cn3 - re3 * Jpy_n3_store(id)) * pt_glb->dt;
			dJpz_n3_rk2(id) = (e0 * wp3 * wp3 * DEz_em_cell(id) * cn3 - re3 * Jpz_n3_store(id)) * pt_glb->dt;
		}
	}

}

#pragma acc routine gang
void EMdynamic_system::get_dJp_RK3() {
	unsigned int mat_type;
	material* mat;
	double wp1, wp2, wp3;
	double re1, re2, re3;
	double cn1, cn2, cn3;

#pragma acc loop gang vector private(mat_type, mat, cn1, cn2, cn3, wp1, wp2, wp3, re1, re2, re3)
	for (long int id = 0; id < n; id++) {
		mat_type = pt_glb->material_cell(id);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);

			cn1 = mat->comp_n1; cn2 = mat->comp_n2; cn3 = mat->comp_n3;
			wp1 = mat->omega_plasma_n1; wp2 = mat->omega_plasma_n2; wp3 = mat->omega_plasma_n3;
			re1 = 1. / mat->tao_e_n1; re2 = 1. / mat->tao_e_n2; re3 = 1. / mat->tao_e_n3;

			dJpx_n1_rk3(id) = (e0 * wp1 * wp1 * DEx_em_cell(id) * cn1 - re1 * Jpx_n1_store(id)) * pt_glb->dt;
			dJpy_n1_rk3(id) = (e0 * wp1 * wp1 * DEy_em_cell(id) * cn1 - re1 * Jpy_n1_store(id)) * pt_glb->dt;
			dJpz_n1_rk3(id) = (e0 * wp1 * wp1 * DEz_em_cell(id) * cn1 - re1 * Jpz_n1_store(id)) * pt_glb->dt;

			dJpx_n2_rk3(id) = (e0 * wp2 * wp2 * DEx_em_cell(id) * cn2 - re2 * Jpx_n2_store(id)) * pt_glb->dt;
			dJpy_n2_rk3(id) = (e0 * wp2 * wp2 * DEy_em_cell(id) * cn2 - re2 * Jpy_n2_store(id)) * pt_glb->dt;
			dJpz_n2_rk3(id) = (e0 * wp2 * wp2 * DEz_em_cell(id) * cn2 - re2 * Jpz_n2_store(id)) * pt_glb->dt;

			dJpx_n3_rk3(id) = (e0 * wp3 * wp3 * DEx_em_cell(id) * cn3 - re3 * Jpx_n3_store(id)) * pt_glb->dt;
			dJpy_n3_rk3(id) = (e0 * wp3 * wp3 * DEy_em_cell(id) * cn3 - re3 * Jpy_n3_store(id)) * pt_glb->dt;
			dJpz_n3_rk3(id) = (e0 * wp3 * wp3 * DEz_em_cell(id) * cn3 - re3 * Jpz_n3_store(id)) * pt_glb->dt;
		}
	}

}

#pragma acc routine gang
void EMdynamic_system::get_dJp_RK4() {
	unsigned int mat_type;
	material* mat;
	double wp1, wp2, wp3;
	double re1, re2, re3;
	double cn1, cn2, cn3;

#pragma acc loop gang vector private(mat_type, mat, cn1, cn2, cn3, wp1, wp2, wp3, re1, re2, re3)
	for (long int id = 0; id < n; id++) {
		mat_type = pt_glb->material_cell(id);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);

			cn1 = mat->comp_n1; cn2 = mat->comp_n2; cn3 = mat->comp_n3;
			wp1 = mat->omega_plasma_n1; wp2 = mat->omega_plasma_n2; wp3 = mat->omega_plasma_n3;
			re1 = 1. / mat->tao_e_n1; re2 = 1. / mat->tao_e_n2; re3 = 1. / mat->tao_e_n3;

			dJpx_n1_rk4(id) = (e0 * wp1 * wp1 * DEx_em_cell(id) * cn1 - re1 * Jpx_n1_store(id)) * pt_glb->dt;
			dJpy_n1_rk4(id) = (e0 * wp1 * wp1 * DEy_em_cell(id) * cn1 - re1 * Jpy_n1_store(id)) * pt_glb->dt;
			dJpz_n1_rk4(id) = (e0 * wp1 * wp1 * DEz_em_cell(id) * cn1 - re1 * Jpz_n1_store(id)) * pt_glb->dt;

			dJpx_n2_rk4(id) = (e0 * wp2 * wp2 * DEx_em_cell(id) * cn2 - re2 * Jpx_n2_store(id)) * pt_glb->dt;
			dJpy_n2_rk4(id) = (e0 * wp2 * wp2 * DEy_em_cell(id) * cn2 - re2 * Jpy_n2_store(id)) * pt_glb->dt;
			dJpz_n2_rk4(id) = (e0 * wp2 * wp2 * DEz_em_cell(id) * cn2 - re2 * Jpz_n2_store(id)) * pt_glb->dt;

			dJpx_n3_rk4(id) = (e0 * wp3 * wp3 * DEx_em_cell(id) * cn3 - re3 * Jpx_n3_store(id)) * pt_glb->dt;
			dJpy_n3_rk4(id) = (e0 * wp3 * wp3 * DEy_em_cell(id) * cn3 - re3 * Jpy_n3_store(id)) * pt_glb->dt;
			dJpz_n3_rk4(id) = (e0 * wp3 * wp3 * DEz_em_cell(id) * cn3 - re3 * Jpz_n3_store(id)) * pt_glb->dt;
		}
	}

}

#pragma acc routine gang
void EMdynamic_system::update_Jp_RK1() {
#pragma acc loop gang vector
	for (long id = 0; id < n; id++) {
		Jpx_n1_store(id) = Jpx_n1(id) + dJpx_n1_rk1(id) * 0.5;
		Jpy_n1_store(id) = Jpy_n1(id) + dJpy_n1_rk1(id) * 0.5;
		Jpz_n1_store(id) = Jpz_n1(id) + dJpz_n1_rk1(id) * 0.5;
		Jpx_n2_store(id) = Jpx_n2(id) + dJpx_n2_rk1(id) * 0.5;
		Jpy_n2_store(id) = Jpy_n2(id) + dJpy_n2_rk1(id) * 0.5;
		Jpz_n2_store(id) = Jpz_n2(id) + dJpz_n2_rk1(id) * 0.5;
		Jpx_n3_store(id) = Jpx_n3(id) + dJpx_n3_rk1(id) * 0.5;
		Jpy_n3_store(id) = Jpy_n3(id) + dJpy_n3_rk1(id) * 0.5;
		Jpz_n3_store(id) = Jpz_n3(id) + dJpz_n3_rk1(id) * 0.5;
	}
}

#pragma acc routine gang
void EMdynamic_system::update_Jp_RK2() {
#pragma acc loop gang vector
	for (long id = 0; id < n; id++) {
		Jpx_n1_store(id) = Jpx_n1(id) + dJpx_n1_rk2(id) * 0.5;
		Jpy_n1_store(id) = Jpy_n1(id) + dJpy_n1_rk2(id) * 0.5;
		Jpz_n1_store(id) = Jpz_n1(id) + dJpz_n1_rk2(id) * 0.5;
		Jpx_n2_store(id) = Jpx_n2(id) + dJpx_n2_rk2(id) * 0.5;
		Jpy_n2_store(id) = Jpy_n2(id) + dJpy_n2_rk2(id) * 0.5;
		Jpz_n2_store(id) = Jpz_n2(id) + dJpz_n2_rk2(id) * 0.5;
		Jpx_n3_store(id) = Jpx_n3(id) + dJpx_n3_rk2(id) * 0.5;
		Jpy_n3_store(id) = Jpy_n3(id) + dJpy_n3_rk2(id) * 0.5;
		Jpz_n3_store(id) = Jpz_n3(id) + dJpz_n3_rk2(id) * 0.5;
	}
}

#pragma acc routine gang
void EMdynamic_system::update_Jp_RK3() {
#pragma acc loop gang vector
	for (long id = 0; id < n; id++) {
		Jpx_n1_store(id) = Jpx_n1(id) + dJpx_n1_rk3(id);
		Jpy_n1_store(id) = Jpy_n1(id) + dJpy_n1_rk3(id);
		Jpz_n1_store(id) = Jpz_n1(id) + dJpz_n1_rk3(id);
		Jpx_n2_store(id) = Jpx_n2(id) + dJpx_n2_rk3(id);
		Jpy_n2_store(id) = Jpy_n2(id) + dJpy_n2_rk3(id);
		Jpz_n2_store(id) = Jpz_n2(id) + dJpz_n2_rk3(id);
		Jpx_n3_store(id) = Jpx_n3(id) + dJpx_n3_rk3(id);
		Jpy_n3_store(id) = Jpy_n3(id) + dJpy_n3_rk3(id);
		Jpz_n3_store(id) = Jpz_n3(id) + dJpz_n3_rk3(id);
	}
}

#pragma acc routine gang
void EMdynamic_system::update_Jp() {
#pragma acc loop gang vector
	for (long int id = 0; id < n; id++) {
		Jpx_n1(id) = Jpx_n1(id) + dJpx_n1_rk1(id) / 6. + dJpx_n1_rk2(id) / 3. + dJpx_n1_rk3(id) / 3. + dJpx_n1_rk4(id) / 6.;
		Jpx_n1_store(id) = Jpx_n1(id);
		Jpy_n1(id) = Jpy_n1(id) + dJpy_n1_rk1(id) / 6. + dJpy_n1_rk2(id) / 3. + dJpy_n1_rk3(id) / 3. + dJpy_n1_rk4(id) / 6.;
		Jpy_n1_store(id) = Jpy_n1(id);
		Jpz_n1(id) = Jpz_n1(id) + dJpz_n1_rk1(id) / 6. + dJpz_n1_rk2(id) / 3. + dJpz_n1_rk3(id) / 3. + dJpz_n1_rk4(id) / 6.;
		Jpz_n1_store(id) = Jpz_n1(id);
		Jpx_n2(id) = Jpx_n2(id) + dJpx_n2_rk1(id) / 6. + dJpx_n2_rk2(id) / 3. + dJpx_n2_rk3(id) / 3. + dJpx_n2_rk4(id) / 6.;
		Jpx_n2_store(id) = Jpx_n2(id);
		Jpy_n2(id) = Jpy_n2(id) + dJpy_n2_rk1(id) / 6. + dJpy_n2_rk2(id) / 3. + dJpy_n2_rk3(id) / 3. + dJpy_n2_rk4(id) / 6.;
		Jpy_n2_store(id) = Jpy_n2(id);
		Jpz_n2(id) = Jpz_n2(id) + dJpz_n2_rk1(id) / 6. + dJpz_n2_rk2(id) / 3. + dJpz_n2_rk3(id) / 3. + dJpz_n2_rk4(id) / 6.;
		Jpz_n2_store(id) = Jpz_n2(id);
		Jpx_n3(id) = Jpx_n3(id) + dJpx_n3_rk1(id) / 6. + dJpx_n3_rk2(id) / 3. + dJpx_n3_rk3(id) / 3. + dJpx_n3_rk4(id) / 6.;
		Jpx_n3_store(id) = Jpx_n3(id);
		Jpy_n3(id) = Jpy_n3(id) + dJpy_n3_rk1(id) / 6. + dJpy_n3_rk2(id) / 3. + dJpy_n3_rk3(id) / 3. + dJpy_n3_rk4(id) / 6.;
		Jpy_n3_store(id) = Jpy_n3(id);
		Jpz_n3(id) = Jpz_n3(id) + dJpz_n3_rk1(id) / 6. + dJpz_n3_rk2(id) / 3. + dJpz_n3_rk3(id) / 3. + dJpz_n3_rk4(id) / 6.;
		Jpz_n3_store(id) = Jpz_n3(id);
	}
}

void EMdynamic_system::update_DE_PEC() {
#pragma acc parallel default(present) async(8)
	{
		//---------------------X component-----------------------//
#pragma acc loop gang vector 
		for (long int idx = 0; idx < nx * (ny + 1) * (nz + 1); idx++) {
			if (DEx_ifPEC(idx) == true) {
				DEx_em(idx) = 0.;
				DEx_em_store(idx) = 0.;
			}
		}

		//---------------------Y component-----------------------//
#pragma acc loop gang vector 
		for (long int idy = 0; idy < (nx + 1) * ny * (nz + 1); idy++) {
			if (DEy_ifPEC(idy) == true) {
				DEy_em(idy) = 0.;
				DEy_em_store(idy) = 0.;
			}
		}

		//---------------------Z component-----------------------//
#pragma acc loop gang vector 
		for (long int idz = 0; idz < (nx + 1) * (ny + 1) * nz; idz++) {
			if (DEz_ifPEC(idz) == true) {
				DEz_em(idz) = 0.;
				DEz_em_store(idz) = 0.;
			}
		}
	}
}
