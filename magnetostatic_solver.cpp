#include "magnetic_system.h"

//void magnetic_system::get_H_static() {
//	unsigned int mat_type;
//	material* mat;
//	double Ms;
//
//	//#pragma acc declare device_resident(Ms,mat_type) deviceptr(mat)
//
//	if (pt_glb->if_mag_1Dmodel == true) {
//		//#pragma acc loop gang 
//		//		for (long int i = 0; i < nx; i++) {
//		//#pragma acc loop worker
//		//			for (long int j = 0; j < ny; j++) {
//		//#pragma acc loop vector private(mat_type,mat,Ms)
//		//				for (long int k = 0; k < nz; k++) {
//
//#pragma acc parallel default(present)
//		{
//#pragma acc loop gang vector private(mat_type,mat,Ms)
//			for (long int i = 0; i < n; i++) {
//				mat_type = pt_glb->material_cell(i);
//				if (mat_type != 0) {
//					mat = &(pt_glb->material_parameters[static_cast<long int>(mat_type) - 1]);
//					if (mat->if_FM == true) {
//						Ms = mat->Ms;
//						//--Calculate stray field in 1D model: Hz = -Msmz--//
//						Hx_stat(i) = 0.;
//						Hy_stat(i) = 0.;
//						Hz_stat(i) = -Ms * mz_glb(i);
//					}
//				}
//			}
//		}
//		//	}
//		//}
//	}
//	//	else {
//	//#pragma acc loop gang
//	//		for (long int i = 0; i < nx; i++) {
//	//#pragma acc loop worker
//	//			for (long int j = 0; j < ny; j++) {
//	//#pragma acc loop vector private(mat_type,mat,Ms)
//	//				for (long int k = 0; k < nz; k++) {
//	//					mat_type = (pt_glb->material_cell(i, j, k));
//	//					if (mat_type != 0) {
//	//						mat = &(pt_glb->material_parameters[static_cast<long int>(mat_type) - 1]);
//	//						if (mat->if_FM == true) {
//	//							Ms = mat->Ms;
//	//							//--Calculate non-normalized M for H_stat--//
//	//							Mx_glb(i, j, k) = Ms * mx_glb(i, j, k);
//	//							My_glb(i, j, k) = Ms * my_glb(i, j, k);
//	//							Mz_glb(i, j, k) = Ms * mz_glb(i, j, k);
//	//						}
//	//					}
//	//				}
//	//			}
//	//		}
//	//#pragma acc wait
//	//
//	////#pragma acc routine(cufftPlan3d)(cufftExecD2Z)(cufftDestroy) gang nohost	
//	////		cufftPlan3d(&plan, nx, ny, nz, CUFFT_D2Z);
//	////		cufftExecD2Z(plan, Mx_glb.matrix, (cufftDoubleComplex*)(Mx_k3D.matrix));
//	////		cufftExecD2Z(plan, My_glb.matrix, (cufftDoubleComplex*)(My_k3D.matrix));
//	////		cufftExecD2Z(plan, Mz_glb.matrix, (cufftDoubleComplex*)(Mz_k3D.matrix));
//	////		cufftDestroy(plan);
//	////
//	////		//pt_math->fourier_transform3D(Mx_glb, Mx_k3D);// mx_k3D.multiply(1. / static_cast<double>(nx * ny * nz));
//	////		//pt_math->fourier_transform3D(My_glb, My_k3D);// my_k3D.multiply(1. / static_cast<double>(nx * ny * nz));
//	////		//pt_math->fourier_transform3D(Mz_glb, Mz_k3D);// mz_k3D.multiply(1. / static_cast<double>(nx * ny * nz));
//	////#pragma acc wait
//	//
//	//#pragma acc loop gang
//	//		for (long int i = 0; i < nx; i++) {
//	//#pragma acc loop worker
//	//			for (long int j = 0; j < ny; j++) {
//	//#pragma acc loop vector
//	//				for (long int k = 0; k < nz / 2 + 1; k++) {
//	//					Hx_stat_k(i, j, k) = \
//	//						(Mx_k3D(i, j, k) * Axyzk(i, j, k) + \
//	//							My_k3D(i, j, k) * Bxyzk(i, j, k) + \
//	//							Mz_k3D(i, j, k) * Bxzyk(i, j, k)) / (-dx * dy * dz);
//	//
//	//					Hy_stat_k(i, j, k) = \
//	//						(My_k3D(i, j, k) * Ayzxk(i, j, k) + \
//	//							Mz_k3D(i, j, k) * Byzxk(i, j, k) + \
//	//							Mx_k3D(i, j, k) * Byxzk(i, j, k)) / (-dx * dy * dz);
//	//
//	//
//	//					Hz_stat_k(i, j, k) = \
//	//						(Mz_k3D(i, j, k) * Azxyk(i, j, k) + \
//	//							Mx_k3D(i, j, k) * Bzxyk(i, j, k) + \
//	//							My_k3D(i, j, k) * Bzyxk(i, j, k)) / (-dx * dy * dz);
//	//				}
//	//			}
//	//		}
//	//#pragma acc wait
//	//
//	//		
//	////		cufftPlan3d(&plan, nx, ny, nz, CUFFT_Z2D);
//	////		cufftExecZ2D(plan, (cufftDoubleComplex*)(Hx_stat_k.matrix), Hx_stat.matrix);
//	////		cufftExecZ2D(plan, (cufftDoubleComplex*)(Hy_stat_k.matrix), Hy_stat.matrix);
//	////		cufftExecZ2D(plan, (cufftDoubleComplex*)(Hz_stat_k.matrix), Hz_stat.matrix);
//	////		cufftDestroy(plan);
//	////
//	////		//pt_math->fourier_transform3D(Hx_stat_k, Hx_stat);
//	////		//pt_math->fourier_transform3D(Hy_stat_k, Hy_stat);
//	////		//pt_math->fourier_transform3D(Hz_stat_k, Hz_stat);
//	////#pragma acc wait
//	//	}
//	//#pragma acc wait
//}

void magnetic_system::get_H_static() {
	unsigned int mat_type;
	material* mat;
	double Ms;
	double Ms1, Ms2;

	if (pt_glb->if_demag_factor == false) {
#pragma acc parallel default(present)
		{
#pragma acc loop gang vector private(mat_type,mat,Ms,Ms1,Ms2)
			for (long int i = 0; i < n; i++) {
				mat_type = pt_glb->material_cell(i);
				if (mat_type != 0) {
					mat = &(pt_glb->material_parameters[mat_type - 1]);
					if (mat->if_FM == true) {
						Ms = mat->Ms;
						//--Calculate non-normalized M for H_stat--//
						Mx_glb(i) = Ms * mx_glb_store(i);
						My_glb(i) = Ms * my_glb_store(i);
						Mz_glb(i) = Ms * mz_glb_store(i);
					}
					else if (mat->if_AFM == true) {
						Ms1 = mat->Ms_AFM1;
						Ms2 = mat->Ms_AFM2;

						Mx_glb(i) = Ms1 * mx_AFM1_glb_store(i) + Ms2 * mx_AFM2_glb_store(i);
						My_glb(i) = Ms1 * my_AFM1_glb_store(i) + Ms2 * my_AFM2_glb_store(i);
						Mz_glb(i) = Ms1 * mz_AFM1_glb_store(i) + Ms2 * mz_AFM2_glb_store(i);
					}
				}
			}
		}

#pragma acc host_data use_device(Mx_glb.matrix, My_glb.matrix, Mz_glb.matrix, Mx_k3D.matrix, My_k3D.matrix, Mz_k3D.matrix)
		{
			cufftExecD2Z(plan_magneto_D2Z, Mx_glb.matrix, (cufftDoubleComplex*)(Mx_k3D.matrix));
			cufftExecD2Z(plan_magneto_D2Z, My_glb.matrix, (cufftDoubleComplex*)(My_k3D.matrix));
			cufftExecD2Z(plan_magneto_D2Z, Mz_glb.matrix, (cufftDoubleComplex*)(Mz_k3D.matrix));
			//
			//		//pt_math->fourier_transform3D(Mx_glb, Mx_k3D);// mx_k3D.multiply(1. / static_cast<double>(nx * ny * nz));
			//		//pt_math->fourier_transform3D(My_glb, My_k3D);// my_k3D.multiply(1. / static_cast<double>(nx * ny * nz));
			//		//pt_math->fourier_transform3D(Mz_glb, Mz_k3D);// mz_k3D.multiply(1. / static_cast<double>(nx * ny * nz));
		}
#pragma acc wait

#pragma acc parallel default(present)
		{
#pragma acc loop gang vector
			for (long int i = 0; i < nx * ny * (nz / 2 + 1); i++) {
				Hx_stat_k(i) = \
					(Mx_k3D(i) * Axyzk(i) + \
						My_k3D(i) * Bxyzk(i) + \
						Mz_k3D(i) * Bxzyk(i)) / (-dx * dy * dz);

				Hy_stat_k(i) = \
					(My_k3D(i) * Ayzxk(i) + \
						Mz_k3D(i) * Byzxk(i) + \
						Mx_k3D(i) * Byxzk(i)) / (-dx * dy * dz);


				Hz_stat_k(i) = \
					(Mz_k3D(i) * Azxyk(i) + \
						Mx_k3D(i) * Bzxyk(i) + \
						My_k3D(i) * Bzyxk(i)) / (-dx * dy * dz);
			}
		}

#pragma acc host_data use_device(Hx_stat_k.matrix, Hy_stat_k.matrix, Hz_stat_k.matrix, Hx_stat.matrix, Hy_stat.matrix, Hz_stat.matrix)
		{
			cufftExecZ2D(plan_magneto_Z2D, (cufftDoubleComplex*)(Hx_stat_k.matrix), Hx_stat.matrix);
			cufftExecZ2D(plan_magneto_Z2D, (cufftDoubleComplex*)(Hy_stat_k.matrix), Hy_stat.matrix);
			cufftExecZ2D(plan_magneto_Z2D, (cufftDoubleComplex*)(Hz_stat_k.matrix), Hz_stat.matrix);
		}
#pragma acc wait
		//pt_math->fourier_transform3D(Hx_stat_k, Hx_stat);
		//pt_math->fourier_transform3D(Hy_stat_k, Hy_stat);
		//pt_math->fourier_transform3D(Hz_stat_k, Hz_stat);
	}
	else if (pt_glb->if_demag_factor == true) {
#pragma acc parallel default(present)
		{
#pragma acc loop gang vector private(mat_type,mat,Ms,Ms1,Ms2)
			for (long int i = 0; i < n; i++) {
				mat_type = pt_glb->material_cell(i);
				if (mat_type != 0) {
					mat = &(pt_glb->material_parameters[mat_type - 1]);
					if (mat->if_FM == true) {
						Ms = mat->Ms;
						//--Calculate non-normalized M for H_stat--//
						Hx_stat(i) = Ms * mx_glb_store(i) * pt_glb->demag_fac_x;
						Hy_stat(i) = Ms * my_glb_store(i) * pt_glb->demag_fac_y;
						Hz_stat(i) = Ms * mz_glb_store(i) * pt_glb->demag_fac_z;
					}
					else if (mat->if_AFM == true) {
						Ms1 = mat->Ms_AFM1;
						Ms2 = mat->Ms_AFM2;

						Hx_stat(i) = (Ms1 * mx_AFM1_glb_store(i) + Ms2 * mx_AFM2_glb_store(i)) * pt_glb->demag_fac_x;
						Hy_stat(i) = (Ms1 * my_AFM1_glb_store(i) + Ms2 * my_AFM2_glb_store(i)) * pt_glb->demag_fac_y;
						Hz_stat(i) = (Ms1 * mz_AFM1_glb_store(i) + Ms2 * mz_AFM2_glb_store(i)) * pt_glb->demag_fac_z;
					}
				}
			}
		}
	}
}

#pragma acc routine gang
void magnetic_system::get_H_static_1D() {
	unsigned int mat_type;
	material* mat;
	double Ms;
	double Ms1, Ms2;
#pragma acc loop gang vector private(mat_type,mat,Ms,Ms1,Ms2)
	for (long int i = 0; i < n; i++) {
		mat_type = pt_glb->material_cell(i);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[mat_type - 1]);
			if (mat->if_FM == true) {
				Ms = mat->Ms;
				//--Calculate stray field in 1D model: Hz = -Msmz--//
				Hx_stat(i) = 0.;
				Hy_stat(i) = 0.;
				Hz_stat(i) = -Ms * mz_glb_store(i);
			}
			else if (mat->if_AFM == true) {
				Ms1 = mat->Ms_AFM1;
				Ms2 = mat->Ms_AFM2;

				Hx_stat(i) = 0.;
				Hy_stat(i) = 0.;
				Hz_stat(i) = -Ms1 * mz_AFM1_glb_store(i) - Ms2 * mz_AFM2_glb_store(i);
			}
		}
	}
}

