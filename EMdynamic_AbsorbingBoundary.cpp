#include "EMdynamic_system.h"

void EMdynamic_system::transfer_pointer() {
	//	//For Liao ABC
	//	if (pt_glb->if_1D_ABC == false) {
	//#pragma acc parallel default(present) async(8)
	//		{
	//			DEx_em_t4 = DEx_em_t3;
	//			DEy_em_t4 = DEy_em_t3;
	//			DEz_em_t4 = DEz_em_t3;
	//		}
	//
	//#pragma acc parallel default(present) async(8)
	//		{
	//			DEx_em_t3 = DEx_em_t2;
	//			DEy_em_t3 = DEy_em_t2;
	//			DEz_em_t3 = DEz_em_t2;
	//		}
	//
	//#pragma acc parallel default(present) async(8)
	//		{
	//			DEx_em_t2 = DEx_em_t1;
	//			DEy_em_t2 = DEy_em_t1;
	//			DEz_em_t2 = DEz_em_t1;
	//		}
	//	}

#pragma acc parallel default(present) async(8)
	{
		DEx_em_t1 = DEx_em;
		DEy_em_t1 = DEy_em;
		DEz_em_t1 = DEz_em;
	}
	//#pragma acc loop gang vector 
	//	for (long int i = 0; i < 1; i++) {
			//DEx_em_t4.nullify_device();
			//DEy_em_t4.nullify_device();
			//DEz_em_t4.nullify_device();
			//DEx_em_t4.get_pointer() = DEx_em_t3.get_pointer();
			//DEy_em_t4.get_pointer() = DEy_em_t3.get_pointer();
			//DEz_em_t4.get_pointer() = DEz_em_t3.get_pointer();
			//DEx_em_t3.get_pointer() = DEx_em_t2.get_pointer();
			//DEy_em_t3.get_pointer() = DEy_em_t2.get_pointer();
			//DEz_em_t3.get_pointer() = DEz_em_t2.get_pointer();
			//DEx_em_t2.get_pointer() = DEx_em_t1.get_pointer();
			//DEy_em_t2.get_pointer() = DEy_em_t1.get_pointer();
			//DEz_em_t2.get_pointer() = DEz_em_t1.get_pointer();
			//DEx_em_t1.get_pointer() = NULL;
			//DEy_em_t1.get_pointer() = NULL;
			//DEz_em_t1.get_pointer() = NULL;
			//DEx_em_t1.initialize_device(nx, ny + 1, nz + 1);
			//DEy_em_t1.initialize_device(nx + 1, ny, nz + 1);
			//DEz_em_t1.initialize_device(nx + 1, ny + 1, nz);
	//	}
	//#pragma acc wait
}

void EMdynamic_system::update_DE_Boundary_half() {
	//double w;
	//double coeff_bot, coeff_top;
	//#pragma acc declare device_resident(w)
	//#pragma acc loop seq
	//	for (long int i = 0; i < 1; i++) {
	//		w = pt_glb->weighting_factor_fourth_Liao;
	//	}
	//#pragma acc wait
	long int i, j, k;
	//--------YZ surface -------------//
	if (pt_geo->periodicX == false) {
#pragma acc parallel default(present) async(8)
		{
			//---------Y component--------//
#pragma acc loop gang vector private(j,k)
			for (long int id = 0; id < ny * (nz + 1); id++) {
				j = id / (nz + 1);
				k = id - j * (nz + 1);

				if (pt_glb->if_PEC_XY == false) {
					//{ //Liao ABC
					//	if (abs(DEy_em_t1(1, j, k)) > 1.e-1) {
					//		DEy_em(0, j, k) = (4. * w + 2. * (1. - w)) * DEy_em_t1(1, j, k) - \
							//			(6. * w + (1. - w)) * DEy_em_t2(2, j, k) + \
							//			(4. * w) * DEy_em_t3(3, j, k) - \
							//			w * DEy_em_t4(4, j, k);
							//	}
							//	else {
							//		DEy_em(0, j, k) = DEy_em_t1(1, j, k);
							//	}
							//	if (abs(DEy_em_t1(nx - 1, j, k)) > 1.e-1) {
							//		DEy_em(nx, j, k) = (4. * w + 2. * (1. - w)) * DEy_em_t1(nx - 1, j, k) - \
							//			(6. * w + (1. - w)) * DEy_em_t2(nx - 2, j, k) + \
							//			(4. * w) * DEy_em_t3(nx - 3, j, k) - \
							//			w * DEy_em_t4(nx - 4, j, k);
							//	}
							//	else {
							//		DEy_em(nx, j, k) = DEy_em_t1(nx - 1, j, k);
							//	}
							//}

					//First order Mur
					DEy_em_store(0, j, k) = DEy_em(1, j, k) + coeff_x_half * (DEy_em(0, j, k) - DEy_em_store(1, j, k));
					DEy_em_store(nx, j, k) = DEy_em(nx - 1, j, k) + coeff_x_half * (DEy_em(nx, j, k) - DEy_em_store(nx - 1, j, k));
				}
				else {
					//PEC
					DEy_em_store(0, j, k) = 0.;
					DEy_em_store(nx, j, k) = 0.;
				}
			}

			//----------Z component--------//
#pragma acc loop gang vector private(j,k)
			for (long int id = 0; id < (ny + 1) * nz; id++) {
				j = id / nz;
				k = id - j * nz;

				if (pt_glb->if_PEC_XY == false) {
					////Liao ABC
					//if (abs(DEz_em_t1(1, j, k)) > 1.e-1) {
					//	DEz_em(0, j, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(1, j, k) - \
							//		(6. * w + (1. - w)) * DEz_em_t2(2, j, k) + \
							//		(4. * w) * DEz_em_t3(3, j, k) - \
							//		w * DEz_em_t4(4, j, k);
							//}
							//else {
							//	DEz_em(0, j, k) = DEz_em_t1(1, j, k);
							//}
							//if (abs(DEz_em_t1(nx - 1, j, k)) > 1.e-1) {
							//	DEz_em(nx, j, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(nx - 1, j, k) - \
							//		(6. * w + (1. - w)) * DEz_em_t2(nx - 2, j, k) + \
							//		(4. * w) * DEz_em_t3(nx - 3, j, k) - \
							//		w * DEz_em_t4(nx - 4, j, k);
							//}
							//else {
							//	DEz_em(nx, j, k) = DEz_em_t1(nx - 1, j, k);
							//}

					//First order Mur ABC
					DEz_em_store(0, j, k) = DEz_em(1, j, k) + coeff_x_half * (DEz_em(0, j, k) - DEz_em_store(1, j, k));
					DEz_em_store(nx, j, k) = DEz_em(nx - 1, j, k) + coeff_x_half * (DEz_em(nx, j, k) - DEz_em_store(nx - 1, j, k));
				}
				else {
					DEz_em_store(0, j, k) = 0.;
					DEz_em_store(nx, j, k) = 0.;
				}
			}
		}
	}

	//--------XZ surface -------------//
	if (pt_geo->periodicY == false) {
#pragma acc parallel default(present) async(8)
		{
			//---------X component--------//
#pragma acc loop gang vector private(i,k)
			for (long int id = 0; id < nx * (nz + 1); id++) {
				i = id / (nz + 1);
				k = id - i * (nz + 1);

				if (pt_glb->if_PEC_XY == false) {
					////Liao ABC
					//w = pt_glb->weighting_factor_fourth_Liao;
					//if (abs(DEx_em_t1(i, 1, k)) > 1.e-1) {
					//	DEx_em(i, 0, k) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, 1, k) - \
					//		(6. * w + (1. - w)) * DEx_em_t2(i, 2, k) + \
					//		(4. * w) * DEx_em_t3(i, 3, k) - \
					//		w * DEx_em_t4(i, 4, k);
					//}
					//else {
					//	DEx_em(i, 0, k) = DEx_em_t1(i, 1, k);
					//}
					//if (abs(DEx_em_t1(i, ny - 1, k)) > 1.e-1) {
					//	DEx_em(i, ny, k) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, ny - 1, k) - \
					//		(6. * w + (1. - w)) * DEx_em_t2(i, ny - 2, k) + \
					//		(4. * w) * DEx_em_t3(i, ny - 3, k) - \
					//		w * DEx_em_t4(i, ny - 4, k);
					//}
					//else {
					//	DEx_em(i, ny, k) = DEx_em_t1(i, ny - 1, k);
					//}

					//First order Mur ABC
					DEx_em_store(i, 0, k) = DEx_em(i, 1, k) + coeff_y_half * (DEx_em(i, 0, k) - DEx_em_store(i, 1, k));
					DEx_em_store(i, ny, k) = DEx_em(i, ny - 1, k) + coeff_y_half * (DEx_em(i, ny, k) - DEx_em_store(i, ny - 1, k));
				}
				else {
					DEx_em_store(i, 0, k) = 0.;
					DEx_em_store(i, ny, k) = 0.;
				}
			}

			//----------Z component--------//
#pragma acc loop gang vector private(i,k)
			for (long int id = 0; id < (nx + 1) * nz; id++) {
				i = id / nz;
				k = id - i * nz;

				if (pt_glb->if_PEC_XY == false) {
					////Liao ABC
					//if (abs(DEz_em_t1(i, 1, k)) > 1.e-1) {
					//	DEz_em(i, 0, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(i, 1, k) - \
					//		(6. * w + (1. - w)) * DEz_em_t2(i, 2, k) + \
					//		(4. * w) * DEz_em_t3(i, 3, k) - \
					//		w * DEz_em_t4(i, 4, k);
					//}
					//else {
					//	DEz_em(i, 0, k) = DEz_em_t1(i, 1, k);
					//}
					//if (abs(DEz_em_t1(i, ny - 1, k)) > 1.e-1) {
					//	DEz_em(i, ny, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(i, ny - 1, k) - \
					//		(6. * w + (1. - w)) * DEz_em_t2(i, ny - 2, k) + \
					//		(4. * w) * DEz_em_t3(i, ny - 3, k) - \
					//		w * DEz_em_t4(i, ny - 4, k);
					//}
					//else {
					//	DEz_em(i, ny, k) = DEz_em_t1(i, ny - 1, k);
					//}

					//First order Mur ABC
					DEz_em_store(i, 0, k) = DEz_em(i, 1, k) + coeff_y_half * (DEz_em(i, 0, k) - DEz_em_store(i, 1, k));
					DEz_em_store(i, ny, k) = DEz_em(i, ny - 1, k) + coeff_y_half * (DEz_em(i, ny, k) - DEz_em_store(i, ny - 1, k));
				}
				else {
					DEz_em_store(i, 0, k) = 0.;
					DEz_em_store(i, ny, k) = 0.;
				}
			}
		}
	}

	//--------XY surface -------------//
	if (pt_geo->periodicZ == false) {
#pragma acc parallel default(present) async(8)
		{
			//---------X component--------//
#pragma acc loop gang vector private(i,j)
			for (long int id = 0; id < nx * (ny + 1); id++) {
				i = id / (ny + 1);
				j = id - i * (ny + 1);

				if (pt_glb->if_PEC_Z == false) {
					if (pt_glb->if_1D_ABC == false) {
						////Liao ABC
						//if (abs(DEx_em_t1(i, j, 1)) > 1.e-1) {
						//	DEx_em(i, j, 0) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, j, 1) - \
						//		(6. * w + (1. - w)) * DEx_em_t2(i, j, 2) + \
						//		(4. * w) * DEx_em_t3(i, j, 3) - \
						//		w * DEx_em_t4(i, j, 4);
						//}
						//else {
						//	DEx_em(i, j, 0) = DEx_em_t1(i, j, 1);
						//}
						//if (abs(DEx_em_t1(i, j, nz - 1)) > 1.e-1) {
						//	DEx_em(i, j, nz) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, j, nz - 1) - \
						//		(6. * w + (1. - w)) * DEx_em_t2(i, j, nz - 2) + \
						//		(4. * w) * DEx_em_t3(i, j, nz - 3) - \
						//		w * DEx_em_t4(i, j, nz - 4);
						//}
						//else {
						//	DEx_em(i, j, nz) = DEx_em_t1(i, j, nz - 1);
						//}

						//First order Mur
						DEx_em_store(i, j, nz) = DEx_em(i, j, nz - 1) + coeff_top_half * (DEx_em(i, j, nz) - DEx_em_store(i, j, nz - 1));
						DEx_em_store(i, j, 0) = DEx_em(i, j, 1) + coeff_bot_half * (DEx_em(i, j, 0) - DEx_em_store(i, j, 1));
					}
					else {
						DEx_em_store(i, j, nz) = DEx_em(i, j, nz - 1) + coeff_top_half * (DEx_em(i, j, nz) - DEx_em_store(i, j, nz - 1));
						if (pt_glb->if_1D_ABC_EM_onlytop == false) {
							DEx_em_store(i, j, 0) = DEx_em(i, j, 1) + coeff_bot_half * (DEx_em(i, j, 0) - DEx_em_store(i, j, 1));
						}
						else {
							DEx_em_store(i, j, 0) = 0.;
						}
					}
				}
				else {
					DEx_em_store(i, j, 0) = 0.;
					DEx_em_store(i, j, nz) = 0.;
				}
			}
			//----------Y component--------//
#pragma acc loop gang vector private(i,j)
			for (long int id = 0; id < (nx + 1) * ny; id++) {
				i = id / ny;
				j = id - i * ny;

				if (pt_glb->if_PEC_Z == false) {
					if (pt_glb->if_1D_ABC == false) {
						////Liao ABC
						//if (abs(DEy_em_t1(i, j, 1)) > 1.e-1) {
						//	DEy_em(i, j, 0) = (4. * w + 2. * (1. - w)) * DEy_em_t1(i, j, 1) - \
						//		(6. * w + (1. - w)) * DEy_em_t2(i, j, 2) + \
						//		(4. * w) * DEy_em_t3(i, j, 3) - \
						//		w * DEy_em_t4(i, j, 4);
						//}
						//else {
						//	DEy_em(i, j, 0) = DEy_em_t1(i, j, 1);
						//}
						//if (abs(DEy_em_t1(i, j, nz - 1)) > 1.e-1) {
						//	DEy_em(i, j, nz) = (4. * w + 2. * (1. - w)) * DEy_em_t1(i, j, nz - 1) - \
						//		(6. * w + (1. - w)) * DEy_em_t2(i, j, nz - 2) + \
						//		(4. * w) * DEy_em_t3(i, j, nz - 3) - \
						//		w * DEy_em_t4(i, j, nz - 4);
						//}
						//else {
						//	DEy_em(i, j, nz) = DEy_em_t1(i, j, nz - 1);
						//}

						//First order Mur ABC
						DEy_em_store(i, j, nz) = DEy_em(i, j, nz - 1) + coeff_top_half * (DEy_em(i, j, nz) - DEy_em_store(i, j, nz - 1));
						DEy_em_store(i, j, 0) = DEy_em(i, j, 1) + coeff_bot_half * (DEy_em(i, j, 0) - DEy_em_store(i, j, 1));
					}
					else {
						DEy_em_store(i, j, nz) = DEy_em(i, j, nz - 1) + coeff_top_half * (DEy_em(i, j, nz) - DEy_em_store(i, j, nz - 1));
						if (pt_glb->if_1D_ABC_EM_onlytop == false) {
							DEy_em_store(i, j, 0) = DEy_em(i, j, 1) + coeff_bot_half * (DEy_em(i, j, 0) - DEy_em_store(i, j, 1));
						}
						else {
							DEy_em_store(i, j, 0) = 0.;
						}
					}
				}
				else {
					DEy_em_store(i, j, 0) = 0.;
					DEy_em_store(i, j, nz) = 0.;
				}
			}
		}
	}
}

void EMdynamic_system::update_DE_Boundary_full() {
	//double w;
	//double coeff_bot, coeff_top;
	//#pragma acc declare device_resident(w)
	//#pragma acc loop seq
	//	for (long int i = 0; i < 1; i++) {
	//		w = pt_glb->weighting_factor_fourth_Liao;
	//	}
	//#pragma acc wait
	long int i, j, k;

	//--------YZ surface -------------//
	if (pt_geo->periodicX == false) {
#pragma acc parallel default(present) async(8)
		{
			//---------Y component--------//
#pragma acc loop gang vector private(j,k)
			for (long int id = 0; id < ny * (nz + 1); id++) {
				j = id / (nz + 1);
				k = id - j * (nz + 1);
				if (pt_glb->if_PEC_XY == false) {
					////Liao ABC
					//if (abs(DEy_em_t1(1, j, k)) > 1.e-1) {
					//	DEy_em(0, j, k) = (4. * w + 2. * (1. - w)) * DEy_em_t1(1, j, k) - \
						//		(6. * w + (1. - w)) * DEy_em_t2(2, j, k) + \
						//		(4. * w) * DEy_em_t3(3, j, k) - \
						//		w * DEy_em_t4(4, j, k);
						//}
						//else {
						//	DEy_em(0, j, k) = DEy_em_t1(1, j, k);
						//}
						//if (abs(DEy_em_t1(nx - 1, j, k)) > 1.e-1) {
						//	DEy_em(nx, j, k) = (4. * w + 2. * (1. - w)) * DEy_em_t1(nx - 1, j, k) - \
						//		(6. * w + (1. - w)) * DEy_em_t2(nx - 2, j, k) + \
						//		(4. * w) * DEy_em_t3(nx - 3, j, k) - \
						//		w * DEy_em_t4(nx - 4, j, k);
						//}
						//else {
						//	DEy_em(nx, j, k) = DEy_em_t1(nx - 1, j, k);
						//}

					//First order Mur ABC
					DEy_em_store(0, j, k) = DEy_em(1, j, k) + coeff_x * (DEy_em(0, j, k) - DEy_em_store(1, j, k));
					DEy_em_store(nx, j, k) = DEy_em(nx - 1, j, k) + coeff_x * (DEy_em(nx, j, k) - DEy_em_store(nx - 1, j, k));
				}
				else {
					DEy_em_store(0, j, k) = 0.;
					DEy_em_store(nx, j, k) = 0.;
				}
			}

			//----------Z component--------//
#pragma acc loop gang vector private(j,k)
			for (long int id = 0; id < (ny + 1) * nz; id++) {
				j = id / nz;
				k = id - j * nz;

				if (pt_glb->if_PEC_XY == false) {
					////Liao ABC
					//if (abs(DEz_em_t1(1, j, k)) > 1.e-1) {
					//	DEz_em(0, j, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(1, j, k) - \
						//		(6. * w + (1. - w)) * DEz_em_t2(2, j, k) + \
						//		(4. * w) * DEz_em_t3(3, j, k) - \
						//		w * DEz_em_t4(4, j, k);
						//}
						//else {
						//	DEz_em(0, j, k) = DEz_em_t1(1, j, k);
						//}
						//if (abs(DEz_em_t1(nx - 1, j, k)) > 1.e-1) {
						//	DEz_em(nx, j, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(nx - 1, j, k) - \
						//		(6. * w + (1. - w)) * DEz_em_t2(nx - 2, j, k) + \
						//		(4. * w) * DEz_em_t3(nx - 3, j, k) - \
						//		w * DEz_em_t4(nx - 4, j, k);
						//}
						//else {
						//	DEz_em(nx, j, k) = DEz_em_t1(nx - 1, j, k);
						//}

					//First order Mur ABC
					DEz_em_store(0, j, k) = DEz_em(1, j, k) + coeff_x * (DEz_em(0, j, k) - DEz_em_store(1, j, k));
					DEz_em_store(nx, j, k) = DEz_em(nx - 1, j, k) + coeff_x * (DEz_em(nx, j, k) - DEz_em_store(nx - 1, j, k));
				}
				else {
					DEz_em_store(0, j, k) = 0.;
					DEz_em_store(nx, j, k) = 0.;
				}
			}
		}
	}

	//--------XZ surface -------------//
	if (pt_geo->periodicY == false) {
#pragma acc parallel default(present) async(8)
		{
			//---------X component--------//
#pragma acc loop gang vector private(i,k)
			for (long int id = 0; id < nx * (nz + 1); id++) {
				i = id / (nz + 1);
				k = id - i * (nz + 1);

				if (pt_glb->if_PEC_XY == false) {
					////Liao ABC
					//if (abs(DEx_em_t1(i, 1, k)) > 1.e-1) {
					//	DEx_em(i, 0, k) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, 1, k) - \
						//		(6. * w + (1. - w)) * DEx_em_t2(i, 2, k) + \
						//		(4. * w) * DEx_em_t3(i, 3, k) - \
						//		w * DEx_em_t4(i, 4, k);
						//}
						//else {
						//	DEx_em(i, 0, k) = DEx_em_t1(i, 1, k);
						//}
						//if (abs(DEx_em_t1(i, ny - 1, k)) > 1.e-1) {
						//	DEx_em(i, ny, k) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, ny - 1, k) - \
						//		(6. * w + (1. - w)) * DEx_em_t2(i, ny - 2, k) + \
						//		(4. * w) * DEx_em_t3(i, ny - 3, k) - \
						//		w * DEx_em_t4(i, ny - 4, k);
						//}
						//else {
						//	DEx_em(i, ny, k) = DEx_em_t1(i, ny - 1, k);
						//}

					//First order Mur ABC
					DEx_em_store(i, 0, k) = DEx_em(i, 1, k) + coeff_y * (DEx_em(i, 0, k) - DEx_em_store(i, 1, k));
					DEx_em_store(i, ny, k) = DEx_em(i, ny - 1, k) + coeff_y * (DEx_em(i, ny, k) - DEx_em_store(i, ny - 1, k));
				}
				else {
					DEx_em_store(i, 0, k) = 0.;
					DEx_em_store(i, ny, k) = 0.;
				}
			}
			//----------Z component--------//
#pragma acc loop gang vector private(i,k)
			for (long int id = 0; id < (nx + 1) * nz; id++) {
				i = id / nz;
				k = id - i * nz;

				if (pt_glb->if_PEC_XY == false) {
					////Liao ABC
					//if (abs(DEz_em_t1(i, 1, k)) > 1.e-1) {
					//	DEz_em(i, 0, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(i, 1, k) - \
						//		(6. * w + (1. - w)) * DEz_em_t2(i, 2, k) + \
						//		(4. * w) * DEz_em_t3(i, 3, k) - \
						//		w * DEz_em_t4(i, 4, k);
						//}
						//else {
						//	DEz_em(i, 0, k) = DEz_em_t1(i, 1, k);
						//}
						//if (abs(DEz_em_t1(i, ny - 1, k)) > 1.e-1) {
						//	DEz_em(i, ny, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(i, ny - 1, k) - \
						//		(6. * w + (1. - w)) * DEz_em_t2(i, ny - 2, k) + \
						//		(4. * w) * DEz_em_t3(i, ny - 3, k) - \
						//		w * DEz_em_t4(i, ny - 4, k);
						//}
						//else {
						//	DEz_em(i, ny, k) = DEz_em_t1(i, ny - 1, k);
						//}

					//First order Mur ABC
					DEz_em_store(i, 0, k) = DEz_em(i, 1, k) + coeff_y * (DEz_em(i, 0, k) - DEz_em_store(i, 1, k));
					DEz_em_store(i, ny, k) = DEz_em(i, ny - 1, k) + coeff_y * (DEz_em(i, ny, k) - DEz_em_store(i, ny - 1, k));
				}
				else {
					DEz_em_store(i, 0, k) = 0.;
					DEz_em_store(i, ny, k) = 0.;
				}
			}
		}
	}

	//--------XY surface -------------//
	if (pt_geo->periodicZ == false) {
#pragma acc parallel default(present) async(8)
		{
			//---------X component--------//
#pragma acc loop gang vector private(i,j)
			for (long int id = 0; id < nx * (ny + 1); id++) {
				i = id / (ny + 1);
				j = id - i * (ny + 1);

				if (pt_glb->if_PEC_Z == false) {
					if (pt_glb->if_1D_ABC == false) {
						////Liao ABC
						//if (abs(DEx_em_t1(i, j, 1)) > 1.e-1) {
						//	DEx_em(i, j, 0) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, j, 1) - \
							//		(6. * w + (1. - w)) * DEx_em_t2(i, j, 2) + \
							//		(4. * w) * DEx_em_t3(i, j, 3) - \
							//		w * DEx_em_t4(i, j, 4);
							//}
							//else {
							//	DEx_em(i, j, 0) = DEx_em_t1(i, j, 1);
							//}
							//if (abs(DEx_em_t1(i, j, nz - 1)) > 1.e-1) {
							//	DEx_em(i, j, nz) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, j, nz - 1) - \
							//		(6. * w + (1. - w)) * DEx_em_t2(i, j, nz - 2) + \
							//		(4. * w) * DEx_em_t3(i, j, nz - 3) - \
							//		w * DEx_em_t4(i, j, nz - 4);
							//}
							//else {
							//	DEx_em(i, j, nz) = DEx_em_t1(i, j, nz - 1);
							//}

						//First order Mur ABC
						DEx_em_store(i, j, nz) = DEx_em(i, j, nz - 1) + coeff_top * (DEx_em(i, j, nz) - DEx_em_store(i, j, nz - 1));
						DEx_em_store(i, j, 0) = DEx_em(i, j, 1) + coeff_bot * (DEx_em(i, j, 0) - DEx_em_store(i, j, 1));
					}
					else {
						DEx_em_store(i, j, nz) = DEx_em(i, j, nz - 1) + coeff_top * (DEx_em(i, j, nz) - DEx_em_store(i, j, nz - 1));
						if (pt_glb->if_1D_ABC_EM_onlytop == false) {
							DEx_em_store(i, j, 0) = DEx_em(i, j, 1) + coeff_bot * (DEx_em(i, j, 0) - DEx_em_store(i, j, 1));
						}
						else {
							DEx_em_store(i, j, 0) = 0.;
						}
					}
				}
				else {
					DEx_em_store(i, j, 0) = 0.;
					DEx_em_store(i, j, nz) = 0.;
				}
			}
			//----------Y component--------//
#pragma acc loop gang vector private(i,j)
			for (long int id = 0; id < (nx + 1) * ny; id++) {
				i = id / ny;
				j = id - i * ny;

				if (pt_glb->if_PEC_Z == false) {
					if (pt_glb->if_1D_ABC == false) {
						////Liao ABC
						//if (abs(DEy_em_t1(i, j, 1)) > 1.e-1) {
						//	DEy_em(i, j, 0) = (4. * w + 2. * (1. - w)) * DEy_em_t1(i, j, 1) - \
							//		(6. * w + (1. - w)) * DEy_em_t2(i, j, 2) + \
							//		(4. * w) * DEy_em_t3(i, j, 3) - \
							//		w * DEy_em_t4(i, j, 4);
							//}
							//else {
							//	DEy_em(i, j, 0) = DEy_em_t1(i, j, 1);
							//}
							//if (abs(DEy_em_t1(i, j, nz - 1)) > 1.e-1) {
							//	DEy_em(i, j, nz) = (4. * w + 2. * (1. - w)) * DEy_em_t1(i, j, nz - 1) - \
							//		(6. * w + (1. - w)) * DEy_em_t2(i, j, nz - 2) + \
							//		(4. * w) * DEy_em_t3(i, j, nz - 3) - \
							//		w * DEy_em_t4(i, j, nz - 4);
							//}
							//else {
							//	DEy_em(i, j, nz) = DEy_em_t1(i, j, nz - 1);
							//}

						//First order Mur ABC
						DEy_em_store(i, j, nz) = DEy_em(i, j, nz - 1) + coeff_top * (DEy_em(i, j, nz) - DEy_em_store(i, j, nz - 1));
						DEy_em_store(i, j, 0) = DEy_em(i, j, 1) + coeff_bot * (DEy_em(i, j, 0) - DEy_em_store(i, j, 1));
					}
					else {
						DEy_em_store(i, j, nz) = DEy_em(i, j, nz - 1) + coeff_top * (DEy_em(i, j, nz) - DEy_em_store(i, j, nz - 1));
						if (pt_glb->if_1D_ABC_EM_onlytop == false) {
							DEy_em_store(i, j, 0) = DEy_em(i, j, 1) + coeff_bot * (DEy_em(i, j, 0) - DEy_em_store(i, j, 1));
						}
						else {
							DEy_em_store(i, j, 0) = 0.;
						}
					}
				}
				else {
					DEy_em_store(i, j, 0) = 0.;
					DEy_em_store(i, j, nz) = 0.;
				}
			}
		}
	}
}

void EMdynamic_system::update_DE_Boundary() {
	//double w;
	//double coeff_bot, coeff_top;
	//#pragma acc declare device_resident(w)
	//#pragma acc loop seq
	//	for (long int i = 0; i < 1; i++) {
	//		w = pt_glb->weighting_factor_fourth_Liao;
	//	}
	//#pragma acc wait
	long int i, j, k;

	//--------YZ surface -------------//
	if (pt_geo->periodicX == false) {
#pragma acc parallel default(present) async(8)
		{
			//---------Y component--------//
#pragma acc loop gang vector private(j,k)
			for (long int id = 0; id < ny * (nz + 1); id++) {
				j = id / (nz + 1);
				k = id - j * (nz + 1);

				if (pt_glb->if_PEC_XY == false) {
					////Liao ABC
					//if (abs(DEy_em_t1(1, j, k)) > 1.e-1) {
					//	DEy_em(0, j, k) = (4. * w + 2. * (1. - w)) * DEy_em_t1(1, j, k) - \
						//		(6. * w + (1. - w)) * DEy_em_t2(2, j, k) + \
						//		(4. * w) * DEy_em_t3(3, j, k) - \
						//		w * DEy_em_t4(4, j, k);
						//}
						//else {
						//	DEy_em(0, j, k) = DEy_em_t1(1, j, k);
						//}
						//if (abs(DEy_em_t1(nx - 1, j, k)) > 1.e-1) {
						//	DEy_em(nx, j, k) = (4. * w + 2. * (1. - w)) * DEy_em_t1(nx - 1, j, k) - \
						//		(6. * w + (1. - w)) * DEy_em_t2(nx - 2, j, k) + \
						//		(4. * w) * DEy_em_t3(nx - 3, j, k) - \
						//		w * DEy_em_t4(nx - 4, j, k);
						//}
						//else {
						//	DEy_em(nx, j, k) = DEy_em_t1(nx - 1, j, k);
						//}

					//First order Mur ABC
					DEy_em(0, j, k) = DEy_em_t1(1, j, k) + coeff_x * (DEy_em_t1(0, j, k) - DEy_em(1, j, k));
					DEy_em(nx, j, k) = DEy_em_t1(nx - 1, j, k) + coeff_x * (DEy_em_t1(nx, j, k) - DEy_em(nx - 1, j, k));
				}
				else {
					DEy_em(0, j, k) = 0.;
					DEy_em(nx, j, k) = 0.;
				}
			}
			//----------Z component--------//
#pragma acc loop gang vector private(j,k)
			for (long int id = 0; id < (ny + 1) * nz; id++) {
				j = id / nz;
				k = id - j * nz;

				if (pt_glb->if_PEC_XY == false) {
					////Liao ABC
					//if (abs(DEz_em_t1(1, j, k)) > 1.e-1) {
					//	DEz_em(0, j, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(1, j, k) - \
						//		(6. * w + (1. - w)) * DEz_em_t2(2, j, k) + \
						//		(4. * w) * DEz_em_t3(3, j, k) - \
						//		w * DEz_em_t4(4, j, k);
						//}
						//else {
						//	DEz_em(0, j, k) = DEz_em_t1(1, j, k);
						//}
						//if (abs(DEz_em_t1(nx - 1, j, k)) > 1.e-1) {
						//	DEz_em(nx, j, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(nx - 1, j, k) - \
						//		(6. * w + (1. - w)) * DEz_em_t2(nx - 2, j, k) + \
						//		(4. * w) * DEz_em_t3(nx - 3, j, k) - \
						//		w * DEz_em_t4(nx - 4, j, k);
						//}
						//else {
						//	DEz_em(nx, j, k) = DEz_em_t1(nx - 1, j, k);
						//}

					//First order Mur ABC
					DEz_em(0, j, k) = DEz_em_t1(1, j, k) + coeff_x * (DEz_em_t1(0, j, k) - DEz_em(1, j, k));
					DEz_em(nx, j, k) = DEz_em_t1(nx - 1, j, k) + coeff_x * (DEz_em_t1(nx, j, k) - DEz_em(nx - 1, j, k));
				}
				else {
					DEz_em(0, j, k) = 0.;
					DEz_em(nx, j, k) = 0.;
				}
			}
		}
	}

	//--------XZ surface -------------//
	if (pt_geo->periodicY == false) {
#pragma acc parallel default(present) async(8)
		{
			//---------X component--------//
#pragma acc loop gang vector private(i,k)
			for (long int id = 0; id < nx * (nz + 1); id++) {
				i = id / (nz + 1);
				k = id - i * (nz + 1);

				if (pt_glb->if_PEC_XY == false) {
					////Liao ABC
					//if (abs(DEx_em_t1(i, 1, k)) > 1.e-1) {
					//	DEx_em(i, 0, k) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, 1, k) - \
						//		(6. * w + (1. - w)) * DEx_em_t2(i, 2, k) + \
						//		(4. * w) * DEx_em_t3(i, 3, k) - \
						//		w * DEx_em_t4(i, 4, k);
						//}
						//else {
						//	DEx_em(i, 0, k) = DEx_em_t1(i, 1, k);
						//}
						//if (abs(DEx_em_t1(i, ny - 1, k)) > 1.e-1) {
						//	DEx_em(i, ny, k) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, ny - 1, k) - \
						//		(6. * w + (1. - w)) * DEx_em_t2(i, ny - 2, k) + \
						//		(4. * w) * DEx_em_t3(i, ny - 3, k) - \
						//		w * DEx_em_t4(i, ny - 4, k);
						//}
						//else {
						//	DEx_em(i, ny, k) = DEx_em_t1(i, ny - 1, k);
						//}

					//First order Mur ABC
					DEx_em(i, 0, k) = DEx_em_t1(i, 1, k) + coeff_y * (DEx_em_t1(i, 0, k) - DEx_em(i, 1, k));
					DEx_em(i, ny, k) = DEx_em_t1(i, ny - 1, k) + coeff_y * (DEx_em_t1(i, ny, k) - DEx_em(i, ny - 1, k));
				}
				else {
					DEx_em(i, 0, k) = 0.;
					DEx_em(i, ny, k) = 0.;
				}
			}
			//----------Z component--------//
#pragma acc loop gang vector private(i,k)
			for (long int id = 0; id < (nx + 1) * nz; id++) {
				i = id / nz;
				k = id - i * nz;
				if (pt_glb->if_PEC_XY == false) {
					////Liao ABC
					//if (abs(DEz_em_t1(i, 1, k)) > 1.e-1) {
					//	DEz_em(i, 0, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(i, 1, k) - \
						//		(6. * w + (1. - w)) * DEz_em_t2(i, 2, k) + \
						//		(4. * w) * DEz_em_t3(i, 3, k) - \
						//		w * DEz_em_t4(i, 4, k);
						//}
						//else {
						//	DEz_em(i, 0, k) = DEz_em_t1(i, 1, k);
						//}
						//if (abs(DEz_em_t1(i, ny - 1, k)) > 1.e-1) {
						//	DEz_em(i, ny, k) = (4. * w + 2. * (1. - w)) * DEz_em_t1(i, ny - 1, k) - \
						//		(6. * w + (1. - w)) * DEz_em_t2(i, ny - 2, k) + \
						//		(4. * w) * DEz_em_t3(i, ny - 3, k) - \
						//		w * DEz_em_t4(i, ny - 4, k);
						//}
						//else {
						//	DEz_em(i, ny, k) = DEz_em_t1(i, ny - 1, k);
						//}

					//First order Mur ABC
					DEz_em(i, 0, k) = DEz_em_t1(i, 1, k) + coeff_y * (DEz_em_t1(i, 0, k) - DEz_em(i, 1, k));
					DEz_em(i, ny, k) = DEz_em_t1(i, ny - 1, k) + coeff_y * (DEz_em_t1(i, ny, k) - DEz_em(i, ny - 1, k));
				}
				else {
					DEz_em(i, 0, k) = 0.;
					DEz_em(i, ny, k) = 0.;
				}
			}
		}
	}

	//--------XY surface -------------//
	if (pt_geo->periodicZ == false) {
#pragma acc parallel default(present) async(8)
		{
			//---------X component--------//
#pragma acc loop gang vector private(i,j)
			for (long int id = 0; id < nx * (ny + 1); id++) {
				i = id / (ny + 1);
				j = id - i * (ny + 1);

				if (pt_glb->if_PEC_Z == false) {
					if (pt_glb->if_1D_ABC == false) {
						////Liao ABC
						//if (abs(DEx_em_t1(i, j, 1)) > 1.e-1) {
						//	DEx_em(i, j, 0) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, j, 1) - \
							//		(6. * w + (1. - w)) * DEx_em_t2(i, j, 2) + \
							//		(4. * w) * DEx_em_t3(i, j, 3) - \
							//		w * DEx_em_t4(i, j, 4);
							//}
							//else {
							//	DEx_em(i, j, 0) = DEx_em_t1(i, j, 1);
							//}
							//if (abs(DEx_em_t1(i, j, nz - 1)) > 1.e-1) {
							//	DEx_em(i, j, nz) = (4. * w + 2. * (1. - w)) * DEx_em_t1(i, j, nz - 1) - \
							//		(6. * w + (1. - w)) * DEx_em_t2(i, j, nz - 2) + \
							//		(4. * w) * DEx_em_t3(i, j, nz - 3) - \
							//		w * DEx_em_t4(i, j, nz - 4);
							//}
							//else {
							//	DEx_em(i, j, nz) = DEx_em_t1(i, j, nz - 1);
							//}

						//First order Mur ABC
						DEx_em(i, j, nz) = DEx_em_t1(i, j, nz - 1) + coeff_top * (DEx_em_t1(i, j, nz) - DEx_em(i, j, nz - 1));
						DEx_em(i, j, 0) = DEx_em_t1(i, j, 1) + coeff_bot * (DEx_em_t1(i, j, 0) - DEx_em(i, j, 1));
					}
					else {
						DEx_em(i, j, nz) = DEx_em_t1(i, j, nz - 1) + coeff_top * (DEx_em_t1(i, j, nz) - DEx_em(i, j, nz - 1));
						if (pt_glb->if_1D_ABC_EM_onlytop == false) {
							DEx_em(i, j, 0) = DEx_em_t1(i, j, 1) + coeff_bot * (DEx_em_t1(i, j, 0) - DEx_em(i, j, 1));
						}
						else {
							DEx_em(i, j, 0) = 0.;
						}
					}
				}
				else {
					DEx_em(i, j, 0) = 0.;
					DEx_em(i, j, nz) = 0.;
				}
			}

			//----------Y component--------//
#pragma acc loop gang vector private(i,j)
			for (long int id = 0; id < (nx + 1) * ny; id++) {
				i = id / ny;
				j = id - i * ny;

				if (pt_glb->if_PEC_Z == false) {
					if (pt_glb->if_1D_ABC == false) {
						////Liao ABC
						//if (abs(DEy_em_t1(i, j, 1)) > 1.e-1) {
						//	DEy_em(i, j, 0) = (4. * w + 2. * (1. - w)) * DEy_em_t1(i, j, 1) - \
						//		(6. * w + (1. - w)) * DEy_em_t2(i, j, 2) + \
						//		(4. * w) * DEy_em_t3(i, j, 3) - \
						//		w * DEy_em_t4(i, j, 4);
						//}
						//else {
						//	DEy_em(i, j, 0) = DEy_em_t1(i, j, 1);
						//}
						//if (abs(DEy_em_t1(i, j, nz - 1)) > 1.e-1) {
						//	DEy_em(i, j, nz) = (4. * w + 2. * (1. - w)) * DEy_em_t1(i, j, nz - 1) - \
						//		(6. * w + (1. - w)) * DEy_em_t2(i, j, nz - 2) + \
						//		(4. * w) * DEy_em_t3(i, j, nz - 3) - \
						//		w * DEy_em_t4(i, j, nz - 4);
						//}
						//else {
						//	DEy_em(i, j, nz) = DEy_em_t1(i, j, nz - 1);
						//}

						//First order Mur ABC
						DEy_em(i, j, nz) = DEy_em_t1(i, j, nz - 1) + coeff_top * (DEy_em_t1(i, j, nz) - DEy_em(i, j, nz - 1));
						DEy_em(i, j, 0) = DEy_em_t1(i, j, 1) + coeff_bot * (DEy_em_t1(i, j, 0) - DEy_em(i, j, 1));
					}
					else {
						DEy_em(i, j, nz) = DEy_em_t1(i, j, nz - 1) + coeff_top * (DEy_em_t1(i, j, nz) - DEy_em(i, j, nz - 1));
						if (pt_glb->if_1D_ABC_EM_onlytop == false) {
							DEy_em(i, j, 0) = DEy_em_t1(i, j, 1) + coeff_bot * (DEy_em_t1(i, j, 0) - DEy_em(i, j, 1));
						}
						else {
							DEy_em(i, j, 0) = 0.;
						}
					}
				}
				else {
					DEy_em(i, j, 0) = 0.;
					DEy_em(i, j, nz) = 0.;
				}
			}
		}
	}
}

