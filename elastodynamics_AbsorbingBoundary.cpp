#include "elastic_system.h"

void elastic_system::transfer_pointer() {
	//#pragma acc loop seq
	//	for (long int i = 0; i < 1; i++) {
	//		Dux_glb_t4.nullify_device();
	//		Duy_glb_t4.nullify_device();
	//		Duz_glb_t4.nullify_device();
	//
	//		Dux_glb_t4.get_pointer() = Dux_glb_t3.get_pointer();
	//		Duy_glb_t4.get_pointer() = Duy_glb_t3.get_pointer();
	//		Duz_glb_t4.get_pointer() = Duz_glb_t3.get_pointer();
	//
	//		Dux_glb_t3.get_pointer() = Dux_glb_t2.get_pointer();
	//		Duy_glb_t3.get_pointer() = Duy_glb_t2.get_pointer();
	//		Duz_glb_t3.get_pointer() = Duz_glb_t2.get_pointer();
	//
	//		Dux_glb_t2.get_pointer() = Dux_glb_t1.get_pointer();
	//		Duy_glb_t2.get_pointer() = Duy_glb_t1.get_pointer();
	//		Duz_glb_t2.get_pointer() = Duz_glb_t1.get_pointer();
	//
	//		Dux_glb_t1.get_pointer() = NULL;
	//		Duy_glb_t1.get_pointer() = NULL;
	//		Duz_glb_t1.get_pointer() = NULL;
	//
	//		Dux_glb_t1.initialize_device(nx, ny, 5);
	//		Duy_glb_t1.initialize_device(nx, ny, 5);
	//		Duz_glb_t1.initialize_device(nx, ny, 5);
	//	}
	//#pragma acc wait
	long int i, j, k;
	if (pt_glb->if_1D_ABC == false) {
#pragma acc parallel default(present) async(2)
		{
#pragma acc loop gang vector private(i,j,k)
			for (long int id = 0; id < nx * ny * nz; id++) {
				i = id / (ny * nz);
				j = (id - i * (ny * nz)) / nz;
				k = id - i * (ny * nz) - j * nz;

				Dux_glb_t4(i, j, k) = Dux_glb_t3(i, j, k);
				Duy_glb_t4(i, j, k) = Duy_glb_t3(i, j, k);
				Duz_glb_t4(i, j, k) = Duz_glb_t3(i, j, k);
			}
		}

#pragma acc parallel default(present) async(2)
		{
#pragma acc loop gang vector private(i,j,k)
			for (long int id = 0; id < nx * ny * nz; id++) {
				i = id / (ny * nz);
				j = (id - i * (ny * nz)) / nz;
				k = id - i * (ny * nz) - j * nz;

				Dux_glb_t3(i, j, k) = Dux_glb_t2(i, j, k);
				Duy_glb_t3(i, j, k) = Duy_glb_t2(i, j, k);
				Duz_glb_t3(i, j, k) = Duz_glb_t2(i, j, k);
			}
		}

#pragma acc parallel default(present) async(2)
		{
#pragma acc loop gang vector private(i,j,k)
			for (long int id = 0; id < nx * ny * nz; id++) {
				i = id / (ny * nz);
				j = (id - i * (ny * nz)) / nz;
				k = id - i * (ny * nz) - j * nz;

				Dux_glb_t2(i, j, k) = Dux_glb_t1(i, j, k);
				Duy_glb_t2(i, j, k) = Duy_glb_t1(i, j, k);
				Duz_glb_t2(i, j, k) = Duz_glb_t1(i, j, k);
			}
		}
	}

#pragma acc parallel default(present) async(2)
	{
#pragma acc loop gang vector private(i,j,k)
		for (long int id = 0; id < nx * ny * nz; id++) {
			i = id / (ny * nz);
			j = (id - i * (ny * nz)) / nz;
			k = id - i * (ny * nz) - j * nz;

			Dux_glb_t1(i, j, k) = Dux_glb(i, j, k);
			Duy_glb_t1(i, j, k) = Duy_glb(i, j, k);
			Duz_glb_t1(i, j, k) = Duz_glb(i, j, k);
		}
	}
}

#pragma acc routine gang// nohost
void elastic_system::update_Du_Boundary_half() {
	unsigned int i, j;
	double w, coeff_long_1st, coeff_tran_1st, coeff_long_last, coeff_tran_last;

#pragma acc loop gang vector private(w,coeff_long_1st,coeff_tran_1st,coeff_long_last,coeff_tran_last, i,j)
	for (long int id = 0; id < nx * ny; id++) {
		i = id / ny;
		j = id - i * ny;

		if (pt_glb->if_1D_ABC == false) {
			w = pt_glb->weighting_factor_fourth_Liao;
			Dux_glb(i, j, 0) = (4. * w + 2. * (1. - w)) * Dux_glb_t1(i, j, 1) - \
				(6. * w + (1. - w)) * Dux_glb_t2(i, j, 2) + \
				(4. * w) * Dux_glb_t3(i, j, 3) - \
				w * Dux_glb_t4(i, j, 4);

			Duy_glb(i, j, 0) = (4. * w + 2. * (1. - w)) * Duy_glb_t1(i, j, 1) - \
				(6. * w + (1. - w)) * Duy_glb_t2(i, j, 2) + \
				(4. * w) * Duy_glb_t3(i, j, 3) - \
				w * Duy_glb_t4(i, j, 4);

			Duz_glb(i, j, 0) = (4. * w + 2. * (1. - w)) * Duz_glb_t1(i, j, 1) - \
				(6. * w + (1. - w)) * Duz_glb_t2(i, j, 2) + \
				(4. * w) * Duz_glb_t3(i, j, 3) - \
				w * Duz_glb_t4(i, j, 4);
		}
		else {
			if (pt_glb->if_1D_ABC_ELAST_onlytop == false) {
				coeff_long_1st = (2. * dz - vs_long_1st * pt_glb->dt) / (2. * dz + vs_long_1st * pt_glb->dt);
				coeff_tran_1st = (2. * dz - vs_tran_1st * pt_glb->dt) / (2. * dz + vs_tran_1st * pt_glb->dt);

				Dux_glb_store(i, j, 0) = Dux_glb(i, j, 1) + coeff_tran_1st * (Dux_glb(i, j, 0) - Dux_glb_store(i, j, 1));
				Duy_glb_store(i, j, 0) = Duy_glb(i, j, 1) + coeff_tran_1st * (Duy_glb(i, j, 0) - Duy_glb_store(i, j, 1));
				Duz_glb_store(i, j, 0) = Duz_glb(i, j, 1) + coeff_long_1st * (Duz_glb(i, j, 0) - Duz_glb_store(i, j, 1));
			}

			coeff_long_last = (2. * dz - vs_long_last * pt_glb->dt) / (2. * dz + vs_long_last * pt_glb->dt);
			coeff_tran_last = (2. * dz - vs_tran_last * pt_glb->dt) / (2. * dz + vs_tran_last * pt_glb->dt);

			Dux_glb_store(i, j, nz - 1) = Dux_glb(i, j, nz - 2) + coeff_tran_last * (Dux_glb(i, j, nz - 1) - Dux_glb_store(i, j, nz - 2));
			Duy_glb_store(i, j, nz - 1) = Duy_glb(i, j, nz - 2) + coeff_tran_last * (Duy_glb(i, j, nz - 1) - Duy_glb_store(i, j, nz - 2));
			Duz_glb_store(i, j, nz - 1) = Duz_glb(i, j, nz - 2) + coeff_long_last * (Duz_glb(i, j, nz - 1) - Duz_glb_store(i, j, nz - 2));
		}
	}
}

#pragma acc routine gang// nohost
void elastic_system::update_Du_Boundary_full() {
	unsigned int i, j;
	double w, coeff_long_1st, coeff_tran_1st, coeff_long_last, coeff_tran_last;

#pragma acc loop gang vector private(w,coeff_long_1st,coeff_tran_1st,coeff_long_last,coeff_tran_last, i,j)
	for (long int id = 0; id < nx * ny; id++) {
		i = id / ny;
		j = id - i * ny;

		if (pt_glb->if_1D_ABC == false) {
			w = pt_glb->weighting_factor_fourth_Liao;
			Dux_glb(i, j, 0) = (4. * w + 2. * (1. - w)) * Dux_glb_t1(i, j, 1) - \
				(6. * w + (1. - w)) * Dux_glb_t2(i, j, 2) + \
				(4. * w) * Dux_glb_t3(i, j, 3) - \
				w * Dux_glb_t4(i, j, 4);

			Duy_glb(i, j, 0) = (4. * w + 2. * (1. - w)) * Duy_glb_t1(i, j, 1) - \
				(6. * w + (1. - w)) * Duy_glb_t2(i, j, 2) + \
				(4. * w) * Duy_glb_t3(i, j, 3) - \
				w * Duy_glb_t4(i, j, 4);

			Duz_glb(i, j, 0) = (4. * w + 2. * (1. - w)) * Duz_glb_t1(i, j, 1) - \
				(6. * w + (1. - w)) * Duz_glb_t2(i, j, 2) + \
				(4. * w) * Duz_glb_t3(i, j, 3) - \
				w * Duz_glb_t4(i, j, 4);
		}
		else {
			if (pt_glb->if_1D_ABC_ELAST_onlytop == false) {
				coeff_long_1st = (dz - vs_long_1st * pt_glb->dt) / (dz + vs_long_1st * pt_glb->dt);
				coeff_tran_1st = (dz - vs_tran_1st * pt_glb->dt) / (dz + vs_tran_1st * pt_glb->dt);

				Dux_glb_store(i, j, 0) = Dux_glb(i, j, 1) + coeff_tran_1st * (Dux_glb(i, j, 0) - Dux_glb_store(i, j, 1));
				Duy_glb_store(i, j, 0) = Duy_glb(i, j, 1) + coeff_tran_1st * (Duy_glb(i, j, 0) - Duy_glb_store(i, j, 1));
				Duz_glb_store(i, j, 0) = Duz_glb(i, j, 1) + coeff_long_1st * (Duz_glb(i, j, 0) - Duz_glb_store(i, j, 1));
			}

			coeff_long_last = (dz - vs_long_last * pt_glb->dt) / (dz + vs_long_last * pt_glb->dt);
			coeff_tran_last = (dz - vs_tran_last * pt_glb->dt) / (dz + vs_tran_last * pt_glb->dt);

			Dux_glb_store(i, j, nz - 1) = Dux_glb(i, j, nz - 2) + coeff_tran_last * (Dux_glb(i, j, nz - 1) - Dux_glb_store(i, j, nz - 2));
			Duy_glb_store(i, j, nz - 1) = Duy_glb(i, j, nz - 2) + coeff_tran_last * (Duy_glb(i, j, nz - 1) - Duy_glb_store(i, j, nz - 2));
			Duz_glb_store(i, j, nz - 1) = Duz_glb(i, j, nz - 2) + coeff_long_last * (Duz_glb(i, j, nz - 1) - Duz_glb_store(i, j, nz - 2));
		}
	}
}

#pragma acc routine gang// nohost
void elastic_system::update_Du_Boundary() {
	unsigned int i, j;
	double w, coeff_long_1st, coeff_tran_1st, coeff_long_last, coeff_tran_last;

#pragma acc loop gang vector private(w,coeff_long_1st,coeff_tran_1st,coeff_long_last,coeff_tran_last, i,j)
	for (long int id = 0; id < nx * ny; id++) {
		i = id / ny;
		j = id - i * ny;

		if (pt_glb->if_1D_ABC == false) {
			w = pt_glb->weighting_factor_fourth_Liao;
			Dux_glb(i, j, 0) = (4. * w + 2. * (1. - w)) * Dux_glb_t1(i, j, 1) - \
				(6. * w + (1. - w)) * Dux_glb_t2(i, j, 2) + \
				(4. * w) * Dux_glb_t3(i, j, 3) - \
				w * Dux_glb_t4(i, j, 4);

			Duy_glb(i, j, 0) = (4. * w + 2. * (1. - w)) * Duy_glb_t1(i, j, 1) - \
				(6. * w + (1. - w)) * Duy_glb_t2(i, j, 2) + \
				(4. * w) * Duy_glb_t3(i, j, 3) - \
				w * Duy_glb_t4(i, j, 4);

			Duz_glb(i, j, 0) = (4. * w + 2. * (1. - w)) * Duz_glb_t1(i, j, 1) - \
				(6. * w + (1. - w)) * Duz_glb_t2(i, j, 2) + \
				(4. * w) * Duz_glb_t3(i, j, 3) - \
				w * Duz_glb_t4(i, j, 4);
		}
		else {
			if (pt_glb->if_1D_ABC_ELAST_onlytop == false) {
				coeff_long_1st = (dz - vs_long_1st * pt_glb->dt) / (dz + vs_long_1st * pt_glb->dt);
				coeff_tran_1st = (dz - vs_tran_1st * pt_glb->dt) / (dz + vs_tran_1st * pt_glb->dt);

				Dux_glb(i, j, 0) = Dux_glb_t1(i, j, 1) + coeff_tran_1st * (Dux_glb_t1(i, j, 0) - Dux_glb(i, j, 1));
				Duy_glb(i, j, 0) = Duy_glb_t1(i, j, 1) + coeff_tran_1st * (Duy_glb_t1(i, j, 0) - Duy_glb(i, j, 1));
				Duz_glb(i, j, 0) = Duz_glb_t1(i, j, 1) + coeff_long_1st * (Duz_glb_t1(i, j, 0) - Duz_glb(i, j, 1));
			}

			coeff_long_last = (dz - vs_long_last * pt_glb->dt) / (dz + vs_long_last * pt_glb->dt);
			coeff_tran_last = (dz - vs_tran_last * pt_glb->dt) / (dz + vs_tran_last * pt_glb->dt);

			Dux_glb(i, j, nz - 1) = Dux_glb_t1(i, j, nz - 2) + coeff_tran_last * (Dux_glb_t1(i, j, nz - 1) - Dux_glb(i, j, nz - 2));
			Duy_glb(i, j, nz - 1) = Duy_glb_t1(i, j, nz - 2) + coeff_tran_last * (Duy_glb_t1(i, j, nz - 1) - Duy_glb(i, j, nz - 2));
			Duz_glb(i, j, nz - 1) = Duz_glb_t1(i, j, nz - 2) + coeff_long_last * (Duz_glb_t1(i, j, nz - 1) - Duz_glb(i, j, nz - 2));
		}
	}
}
