#include "elastic_system.h"

//void elastic_system::initialize(magnetic_system* pt_mag_in, ferroelectric_system* pt_ele_in)
void elastic_system::initialize_host() {

	//pt_mag = pt_mag_in;
	//pt_ele = pt_ele_in;

	nx = pt_geo->nx_system;
	ny = pt_geo->ny_system;
	nz = pt_geo->nz_system;
	nz21 = nz / 2 + 1;

	n = nx * ny * nz;
	scalenn = 1. / static_cast<double>(n);

	dx = pt_geo->dx;
	dy = pt_geo->dy;
	dz = pt_geo->dz;

	//		Calculate homogeneous elastic stiffness for static solver//
	c11_homo = pt_glb->material_parameters[0].c11;
	c12_homo = pt_glb->material_parameters[0].c12;
	c44_homo = pt_glb->material_parameters[0].c44;

	for (long int i = 1; i < pt_glb->num_materials; i++) {
		if (pt_glb->material_parameters[i].c11 > c11_homo) {
			c11_homo = pt_glb->material_parameters[i].c11;
		}
		if (pt_glb->material_parameters[i].c12 > c12_homo) {
			c12_homo = pt_glb->material_parameters[i].c12;
		}
		if (pt_glb->material_parameters[i].c44 > c44_homo) {
			c44_homo = pt_glb->material_parameters[i].c44;
		}
	}

	c11_homo = c11_homo * pt_glb->scale_elasto;
	c12_homo = c12_homo * pt_glb->scale_elasto;
	c44_homo = c44_homo * pt_glb->scale_elasto;
	//		END Calculate homogeneous elastic stiffness for static solver//

	// Calculate longitudinal and transverse sound speed in first and last layers for ABC
	int i, j;
	unsigned int mat_type;
	material* mat;

	c11_1st = 0.; c44_1st = 0.; density_1st = 1.e-8;
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (pt_glb->material_cell(i, j, 0) != 0) {
				mat_type = pt_glb->material_cell(i, j, 0);
				mat = &(pt_glb->material_parameters[mat_type - 1]);
				c11_1st = mat->c11;
				c44_1st = mat->c44;
				density_1st = mat->density;
				j = ny; i = nx;
			}
		}
	}

	c11_last = 0.; c44_last = 0.; density_last = 1.e-8;
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (pt_glb->material_cell(i, j, nz-1) != 0) {
				mat_type = pt_glb->material_cell(i, j, nz - 1);
				mat = &(pt_glb->material_parameters[mat_type - 1]);
				c11_last = mat->c11;
				c44_last = mat->c44;
				density_last = mat->density;
				j = ny; i = nx;
			}
		}
	}

	vs_long_1st = sqrt((c11_1st / density_1st));
	vs_tran_1st = sqrt((c44_1st / density_1st));

	vs_long_last = sqrt((c11_last / density_last));
	vs_tran_last = sqrt((c44_last / density_last));


	source_z = pt_glb->input_strain_sourcez;
	//if (pt_glb->if_strain_input_bottom == true) {
	//	source_z = 0;
	//}
	//else {
	//	source_z = pt_geo->pt_nz_layer[0] / 5;
	//}
	//END Calculating sound speed

	//Initialization
	Gik_11.initialize(nx, ny, nz21); Gik_12.initialize(nx, ny, nz21); Gik_13.initialize(nx, ny, nz21);
	Gik_21.initialize(nx, ny, nz21); Gik_22.initialize(nx, ny, nz21); Gik_23.initialize(nx, ny, nz21);
	Gik_31.initialize(nx, ny, nz21); Gik_32.initialize(nx, ny, nz21); Gik_33.initialize(nx, ny, nz21);

	duxdx_glb_surfYZ.initialize(nx + 1, ny, nz); duydx_glb_surfYZ.initialize(nx + 1, ny, nz); duzdx_glb_surfYZ.initialize(nx + 1, ny, nz);
	duxdy_glb_surfXZ.initialize(nx, ny + 1, nz); duydy_glb_surfXZ.initialize(nx, ny + 1, nz); duzdy_glb_surfXZ.initialize(nx, ny + 1, nz);
	duxdz_glb_surfXY.initialize(nx, ny, nz + 1); duydz_glb_surfXY.initialize(nx, ny, nz + 1); duzdz_glb_surfXY.initialize(nx, ny, nz + 1);

	Dux_glb.initialize(nx, ny, nz); Duy_glb.initialize(nx, ny, nz); Duz_glb.initialize(nx, ny, nz);
	Dux_glb_store.initialize(nx, ny, nz); Duy_glb_store.initialize(nx, ny, nz); Duz_glb_store.initialize(nx, ny, nz);
	//dDux_glb.initialize(nx, ny, nz); dDuy_glb.initialize(nx, ny, nz); dDuz_glb.initialize(nx, ny, nz);
	dDux_glb_rk1.initialize(nx, ny, nz); dDuy_glb_rk1.initialize(nx, ny, nz); dDuz_glb_rk1.initialize(nx, ny, nz);
	dDux_glb_rk2.initialize(nx, ny, nz); dDuy_glb_rk2.initialize(nx, ny, nz); dDuz_glb_rk2.initialize(nx, ny, nz);
	dDux_glb_rk3.initialize(nx, ny, nz); dDuy_glb_rk3.initialize(nx, ny, nz); dDuz_glb_rk3.initialize(nx, ny, nz);
	dDux_glb_rk4.initialize(nx, ny, nz); dDuy_glb_rk4.initialize(nx, ny, nz); dDuz_glb_rk4.initialize(nx, ny, nz);

	vx_glb.initialize(nx, ny, nz); vy_glb.initialize(nx, ny, nz); vz_glb.initialize(nx, ny, nz);
	vx_glb_store.initialize(nx, ny, nz); vy_glb_store.initialize(nx, ny, nz); vz_glb_store.initialize(nx, ny, nz);
	//dvx_glb.initialize(nx, ny, nz); dvy_glb.initialize(nx, ny, nz); dvz_glb.initialize(nx, ny, nz);

	force_x.initialize(nx, ny, nz);			force_y.initialize(nx, ny, nz);			force_z.initialize(nx, ny, nz);
	force_x_store.initialize(nx, ny, nz);	force_y_store.initialize(nx, ny, nz);	force_z_store.initialize(nx, ny, nz);

	dvx_glb_rk1.initialize(nx, ny, nz); dvy_glb_rk1.initialize(nx, ny, nz); dvz_glb_rk1.initialize(nx, ny, nz);
	dvx_glb_rk2.initialize(nx, ny, nz); dvy_glb_rk2.initialize(nx, ny, nz); dvz_glb_rk2.initialize(nx, ny, nz);
	dvx_glb_rk3.initialize(nx, ny, nz); dvy_glb_rk3.initialize(nx, ny, nz); dvz_glb_rk3.initialize(nx, ny, nz);
	dvx_glb_rk4.initialize(nx, ny, nz); dvy_glb_rk4.initialize(nx, ny, nz); dvz_glb_rk4.initialize(nx, ny, nz);

	Dexx_glb.initialize(nx, ny, nz); Deyy_glb.initialize(nx, ny, nz); Dezz_glb.initialize(nx, ny, nz);
	Deyz_glb.initialize(nx, ny, nz); Dexz_glb.initialize(nx, ny, nz); Dexy_glb.initialize(nx, ny, nz);

	Dexx_crt.initialize(nx, ny, nz); Deyy_crt.initialize(nx, ny, nz); Dezz_crt.initialize(nx, ny, nz);
	Deyz_crt.initialize(nx, ny, nz); Dexz_crt.initialize(nx, ny, nz); Dexy_crt.initialize(nx, ny, nz);

	exxt0_glb.initialize(nx, ny, nz); eyyt0_glb.initialize(nx, ny, nz); ezzt0_glb.initialize(nx, ny, nz);
	eyzt0_glb.initialize(nx, ny, nz); exzt0_glb.initialize(nx, ny, nz); exyt0_glb.initialize(nx, ny, nz);

	exxt0_crt.initialize(nx, ny, nz); eyyt0_crt.initialize(nx, ny, nz); ezzt0_crt.initialize(nx, ny, nz);
	eyzt0_crt.initialize(nx, ny, nz); exzt0_crt.initialize(nx, ny, nz); exyt0_crt.initialize(nx, ny, nz);

	Dexx0_glb.initialize(nx, ny, nz); Deyy0_glb.initialize(nx, ny, nz); Dezz0_glb.initialize(nx, ny, nz);
	Deyz0_glb.initialize(nx, ny, nz); Dexz0_glb.initialize(nx, ny, nz); Dexy0_glb.initialize(nx, ny, nz);

	Dexx0_crt.initialize(nx, ny, nz); Deyy0_crt.initialize(nx, ny, nz); Dezz0_crt.initialize(nx, ny, nz);
	Deyz0_crt.initialize(nx, ny, nz); Dexz0_crt.initialize(nx, ny, nz); Dexy0_crt.initialize(nx, ny, nz);

	exx0t0_crt.initialize(nx, ny, nz); eyy0t0_crt.initialize(nx, ny, nz); ezz0t0_crt.initialize(nx, ny, nz);
	eyz0t0_crt.initialize(nx, ny, nz); exz0t0_crt.initialize(nx, ny, nz); exy0t0_crt.initialize(nx, ny, nz);

	Dux_glb_t1.initialize(nx, ny, nz); Duy_glb_t1.initialize(nx, ny, nz); Duz_glb_t1.initialize(nx, ny, nz);
	if (pt_glb->if_1D_ABC == false) {
		Dux_glb_t2.initialize(nx, ny, nz); Duy_glb_t2.initialize(nx, ny, nz); Duz_glb_t2.initialize(nx, ny, nz);
		Dux_glb_t3.initialize(nx, ny, nz); Duy_glb_t3.initialize(nx, ny, nz); Duz_glb_t3.initialize(nx, ny, nz);
		Dux_glb_t4.initialize(nx, ny, nz); Duy_glb_t4.initialize(nx, ny, nz); Duz_glb_t4.initialize(nx, ny, nz);
	}

	stress_free_surfYZ.initialize(nx + 1, ny, nz);
	stress_free_surfXZ.initialize(nx, ny + 1, nz);
	stress_free_surfXY.initialize(nx, ny, nz + 1);

	//--------Initialize intermediate variables for flexoelectrics--------//
	if (pt_glb->if_FE_all == true && pt_glb->if_flexo == true) {
		dexxt0dx_glb.initialize(nx, ny, nz); dexxt0dy_glb.initialize(nx, ny, nz); dexxt0dz_glb.initialize(nx, ny, nz);
		deyyt0dx_glb.initialize(nx, ny, nz); deyyt0dy_glb.initialize(nx, ny, nz); deyyt0dz_glb.initialize(nx, ny, nz);
		dezzt0dx_glb.initialize(nx, ny, nz); dezzt0dy_glb.initialize(nx, ny, nz); dezzt0dz_glb.initialize(nx, ny, nz);
		deyzt0dx_glb.initialize(nx, ny, nz); deyzt0dy_glb.initialize(nx, ny, nz); deyzt0dz_glb.initialize(nx, ny, nz);
		dexzt0dx_glb.initialize(nx, ny, nz); dexzt0dy_glb.initialize(nx, ny, nz); dexzt0dz_glb.initialize(nx, ny, nz);
		dexyt0dx_glb.initialize(nx, ny, nz); dexyt0dy_glb.initialize(nx, ny, nz); dexyt0dz_glb.initialize(nx, ny, nz);

		dDexxdx_glb.initialize(nx, ny, nz); dDexxdy_glb.initialize(nx, ny, nz); dDexxdz_glb.initialize(nx, ny, nz);
		dDeyydx_glb.initialize(nx, ny, nz); dDeyydy_glb.initialize(nx, ny, nz); dDeyydz_glb.initialize(nx, ny, nz);
		dDezzdx_glb.initialize(nx, ny, nz); dDezzdy_glb.initialize(nx, ny, nz); dDezzdz_glb.initialize(nx, ny, nz);
		dDeyzdx_glb.initialize(nx, ny, nz); dDeyzdy_glb.initialize(nx, ny, nz); dDeyzdz_glb.initialize(nx, ny, nz);
		dDexzdx_glb.initialize(nx, ny, nz); dDexzdy_glb.initialize(nx, ny, nz); dDexzdz_glb.initialize(nx, ny, nz);
		dDexydx_glb.initialize(nx, ny, nz); dDexydy_glb.initialize(nx, ny, nz); dDexydz_glb.initialize(nx, ny, nz);
	}

	//--------Initialize intermediate variables for elasto-static solver--------//	
	exx_ext_crt.initialize(nx, ny, nz); eyy_ext_crt.initialize(nx, ny, nz); ezz_ext_crt.initialize(nx, ny, nz);
	exy_ext_crt.initialize(nx, ny, nz); exz_ext_crt.initialize(nx, ny, nz); eyz_ext_crt.initialize(nx, ny, nz);

	if (pt_glb->if_read_strain_ext == true) {
		read_external_strain();
	}

	stress0homo_xx.initialize(nx, ny, nz); stress0homo_yy.initialize(nx, ny, nz); stress0homo_zz.initialize(nx, ny, nz);
	stress0homo_xy.initialize(nx, ny, nz); stress0homo_xz.initialize(nx, ny, nz); stress0homo_yz.initialize(nx, ny, nz);

	stress0_xx.initialize(nx, ny, nz); stress0_yy.initialize(nx, ny, nz); stress0_zz.initialize(nx, ny, nz);
	stress0_xy.initialize(nx, ny, nz); stress0_xz.initialize(nx, ny, nz); stress0_yz.initialize(nx, ny, nz);

	stress_xx.initialize(nx, ny, nz); stress_yy.initialize(nx, ny, nz); stress_zz.initialize(nx, ny, nz);
	stress_xy.initialize(nx, ny, nz); stress_xz.initialize(nx, ny, nz); stress_yz.initialize(nx, ny, nz);

	duxdx.initialize(nx, ny, nz); duxdy.initialize(nx, ny, nz); duxdz.initialize(nx, ny, nz);
	duydx.initialize(nx, ny, nz); duydy.initialize(nx, ny, nz); duydz.initialize(nx, ny, nz);
	duzdx.initialize(nx, ny, nz); duzdy.initialize(nx, ny, nz); duzdz.initialize(nx, ny, nz);

	stress0homok_xx.initialize(nx, ny, nz21); stress0homok_yy.initialize(nx, ny, nz21); stress0homok_zz.initialize(nx, ny, nz21);
	stress0homok_xy.initialize(nx, ny, nz21); stress0homok_xz.initialize(nx, ny, nz21); stress0homok_yz.initialize(nx, ny, nz21);

	stressk_xx.initialize(nx, ny, nz21); stressk_yy.initialize(nx, ny, nz21); stressk_zz.initialize(nx, ny, nz21);
	stressk_xy.initialize(nx, ny, nz21); stressk_xz.initialize(nx, ny, nz21); stressk_yz.initialize(nx, ny, nz21);

	duxdxk.initialize(nx, ny, nz21); duxdyk.initialize(nx, ny, nz21); duxdzk.initialize(nx, ny, nz21);
	duydxk.initialize(nx, ny, nz21); duydyk.initialize(nx, ny, nz21); duydzk.initialize(nx, ny, nz21);
	duzdxk.initialize(nx, ny, nz21); duzdyk.initialize(nx, ny, nz21); duzdzk.initialize(nx, ny, nz21);

	if (pt_glb->if_elastostatic_film == true) {
		exx0_ave_film = new double[nz];
		eyy0_ave_film = new double[nz];
		ezz0_ave_film = new double[nz];

		eyz0_ave_film = new double[nz];
		exz0_ave_film = new double[nz];
		exy0_ave_film = new double[nz];

		for (unsigned long i = 0; i < nz; i++) {
			exx0_ave_film[i] = 0.;
			eyy0_ave_film[i] = 0.;
			ezz0_ave_film[i] = 0.;

			eyz0_ave_film[i] = 0.;
			exz0_ave_film[i] = 0.;
			exy0_ave_film[i] = 0.;
		}
	}

	e_sum = 0.;
	e_elas[0] = 0.; e_elas[1] = 0.; e_elas[2] = 0.; e_elas[3] = 0.;
	tolerance = 0.;
	cufftPlan3d(&plan_elasto_D2Z, nx, ny, nz, CUFFT_D2Z);
	cufftPlan3d(&plan_elasto_Z2D, nx, ny, nz, CUFFT_Z2D);
}

void elastic_system::initialize_device() {
	double chi;
	double kqx, kqy, kqz, knorm;
	double d0;

	long i, j, k;

	double strain_in[3][3];
	double strain_out[3][3];
	unsigned int mat_type;
	material* mat;
	double c11, c12, c44;

	//------------------------------------------//
	//		Determine if stress-free surfaces	//
	//------------------------------------------//

#pragma acc parallel default(present) async(5)
	{
		//------------X Direction (YZ plane)----------------//
#pragma acc loop gang vector private(i,j,k)
		for (long int id = 0; id < (nx + 1) * ny * nz; id++) {
			i = id / (ny * nz);
			j = (id - i * ny * nz) / nz;
			k = id - i * ny * nz - j * nz;

			if (i == 0 || i == nx) {
				if (pt_geo->periodicX == true) {
					if (pt_glb->material_cell(0, j, k) != 0 && pt_glb->material_cell(nx - 1, j, k) != 0) {
						stress_free_surfYZ(i, j, k) = false;
					}
					else {
						stress_free_surfYZ(i, j, k) = true;
					}
				}
				else if (pt_geo->periodicX == false) {
					stress_free_surfYZ(i, j, k) = true;
				}
			}
			else {
				if (pt_glb->material_cell(i - 1, j, k) != 0 && pt_glb->material_cell(i, j, k) != 0) {
					stress_free_surfYZ(i, j, k) = false;
				}
				else {
					stress_free_surfYZ(i, j, k) = true;
				}
			}
		}

		//------------Y Direction (XZ plane)----------------//
#pragma acc loop gang vector private(i,j,k)
		for (long int id = 0; id < nx * (ny + 1) * nz; id++) {
			i = id / ((ny + 1) * nz);
			j = (id - i * (ny + 1) * nz) / nz;
			k = id - i * (ny + 1) * nz - j * nz;

			if (j == 0 || j == ny) {
				if (pt_geo->periodicY == true) {
					if (pt_glb->material_cell(i, 0, k) != 0 && pt_glb->material_cell(i, ny - 1, k) != 0) {
						stress_free_surfXZ(i, j, k) = false;
					}
					else {
						stress_free_surfXZ(i, j, k) = true;
					}
				}
				else if (pt_geo->periodicY == false) {
					stress_free_surfXZ(i, j, k) = true;
				}
			}
			else {
				if (pt_glb->material_cell(i, j - 1, k) != 0 && pt_glb->material_cell(i, j, k) != 0) {
					stress_free_surfXZ(i, j, k) = false;
				}
				else {
					stress_free_surfXZ(i, j, k) = true;
				}
			}
		}

		//------------Z Direction (XY plane)----------------//
#pragma acc loop gang vector private(i,j,k)
		for (long int id = 0; id < nx * ny * (nz + 1); id++) {
			i = id / (ny * (nz + 1));
			j = (id - i * ny * (nz + 1)) / (nz + 1);
			k = id - i * ny * (nz + 1) - j * (nz + 1);

			if (k == 0 || k == nz) {
				if (pt_geo->periodicZ == true) {
					if (pt_glb->material_cell(i, j, 0) != 0 && pt_glb->material_cell(i, j, nz - 1) != 0) {
						stress_free_surfXY(i, j, k) = false;
					}
					else {
						stress_free_surfXY(i, j, k) = true;
					}
				}
				else if (pt_geo->periodicZ == false) {
					stress_free_surfXY(i, j, k) = true;
				}
			}
			else {
				if (pt_glb->material_cell(i, j, k - 1) != 0 && pt_glb->material_cell(i, j, k) != 0) {
					stress_free_surfXY(i, j, k) = false;
				}
				else {
					stress_free_surfXY(i, j, k) = true;
				}
			}
		}
		//		END Determine if stress-free surfaces	//


		//----Initialize intermediate variables for elastostatic solver------//
#pragma acc loop gang vector private(chi,kqx, kqy, kqz, knorm,d0,i,j,k)
		for (long int id = 0; id < nx * ny * nz21; id++) {
			i = id / (ny * nz21);
			j = (id - i * ny * nz21) / nz21;
			k = id - i * ny * nz21 - j * nz21;

			chi = (c11_homo - c12_homo - 2. * c44_homo) / c44_homo;
			pt_math->transform_vector_glb2crt(pt_glb->kx[i], pt_glb->ky[j], pt_glb->kz[k], \
				kqx, kqy, kqz);
			knorm = pt_glb->k_norm(i, j, k);

			if (knorm < 1.e-8) {
				d0 = 1;
			}
			else {
				d0 = c11_homo * pow(knorm, 3.) + chi * (c11_homo + c12_homo) * knorm * \
					(kqx * kqx * kqy * kqy \
						+ kqy * kqy * kqz * kqz \
						+ kqz * kqz * kqx * kqx) \
					+ chi * chi * (c11_homo + 2 * c12_homo + c44_homo) * \
					kqx * kqx * kqy * kqy * kqz * kqz;
			}

			Gik_11(i, j, k) = (c44_homo * pow(knorm, 2.) + \
				(c11_homo + c12_homo) * chi * kqy * kqy * kqz * kqz \
				+ (c11_homo - c44_homo) * knorm * (kqy * kqy + kqz * kqz)) \
				/ (c44_homo * d0);

			Gik_22(i, j, k) = (c44_homo * pow(knorm, 2.) + \
				(c11_homo + c12_homo) * chi * kqz * kqz * kqx * kqx \
				+ (c11_homo - c44_homo) * knorm * (kqz * kqz + kqx * kqx)) \
				/ (c44_homo * d0);

			Gik_33(i, j, k) = (c44_homo * pow(knorm, 2.) + \
				(c11_homo + c12_homo) * chi * kqx * kqx * kqy * kqy \
				+ (c11_homo - c44_homo) * knorm * (kqx * kqx + kqy * kqy)) \
				/ (c44_homo * d0);

			Gik_12(i, j, k) = -(c12_homo + c44_homo) * (knorm + chi * kqz * kqz) \
				* kqx * kqy / (c44_homo * d0);
			Gik_21(i, j, k) = Gik_12(i, j, k);

			Gik_23(i, j, k) = -(c12_homo + c44_homo) * (knorm + chi * kqx * kqx) \
				* kqy * kqz / (c44_homo * d0);
			Gik_32(i, j, k) = Gik_23(i, j, k);

			Gik_13(i, j, k) = -(c12_homo + c44_homo) * (knorm + chi * kqy * kqy) \
				* kqx * kqz / (c44_homo * d0);
			Gik_31(i, j, k) = Gik_13(i, j, k);
		}
	}

	if (pt_glb->if_input_uandv == true) {
#pragma acc parallel default(present)  async(6)
		{
#pragma acc loop gang vector
			for (long id = 0; id < n; id++) {
				vx_glb_store(id) = vx_glb(id);
				vy_glb_store(id) = vy_glb(id);
				vz_glb_store(id) = vz_glb(id);
				Dux_glb_store(id) = Dux_glb(id);
				Duy_glb_store(id) = Duy_glb(id);
				Duz_glb_store(id) = Duz_glb(id);
			}
		}

#pragma acc wait(5) async(6)
#pragma acc parallel default(present) async(6) 
		{
			get_dudxyz_glb();
		}
#pragma acc parallel default(present) async(6) 
		{
			get_Dstrain();
		}

		if (pt_glb->if_FE_all == true && pt_glb->if_flexo == true) {
#pragma acc parallel default(present) async(6) 
			{
				get_Dstrain_gradient();
			}
		}
	}

	if (pt_glb->if_input_straint0 == true) {
#pragma acc parallel default(present) async(7)
		{
#pragma acc loop gang vector private(strain_in, strain_out, mat_type)
			for (unsigned long id = 0; id < n; id++) {
				mat_type = pt_glb->material_cell(id);
				if (mat_type == 0) {
					exxt0_crt(id) = 0.;
					eyyt0_crt(id) = 0.;
					ezzt0_crt(id) = 0.;
					eyzt0_crt(id) = 0.;
					exzt0_crt(id) = 0.;
					exyt0_crt(id) = 0.;
				}
				else {
					strain_in[0][0] = exxt0_glb(id);
					strain_in[1][1] = eyyt0_glb(id);
					strain_in[2][2] = ezzt0_glb(id);
					strain_in[2][1] = eyzt0_glb(id);
					strain_in[1][2] = strain_in[2][1];
					strain_in[2][0] = exzt0_glb(id);
					strain_in[0][2] = strain_in[2][0];
					strain_in[1][0] = exyt0_glb(id);
					strain_in[0][1] = strain_in[1][0];

					pt_math->transform_matrix_glb2crt(strain_in, strain_out);

					exxt0_crt(id) = strain_out[0][0];
					eyyt0_crt(id) = strain_out[1][1];
					ezzt0_crt(id) = strain_out[2][2];
					eyzt0_crt(id) = strain_out[1][2];
					exzt0_crt(id) = strain_out[0][2];
					exyt0_crt(id) = strain_out[0][1];
				}
			}
		}

		if (pt_glb->if_FE_all == true && pt_glb->if_flexo == true) {
#pragma acc wait(5) async(7)
#pragma acc parallel default(present) async(7)
			{
				get_strain_static_glb_gradient();
			}
		}
	}

#pragma acc parallel default(present) async(8)
	{
#pragma acc loop gang vector private(strain_in, strain_out,i,j,k, c11, c12, c44, mat_type, mat)
		for (unsigned long id = 0; id < n; id++) {
			if (pt_glb->if_read_strain_ext == true) {
				strain_in[0][0] = exx_ext_crt(id);
				strain_in[1][1] = eyy_ext_crt(id);
				strain_in[2][2] = ezz_ext_crt(id);
				strain_in[2][1] = eyz_ext_crt(id);
				strain_in[1][2] = strain_in[2][1];
				strain_in[2][0] = exz_ext_crt(id);
				strain_in[0][2] = strain_in[2][0];
				strain_in[1][0] = exy_ext_crt(id);
				strain_in[0][1] = strain_in[1][0];
			}
			else if (pt_glb->if_read_strain_ext == false) {
				i = id / (ny * nz);
				j = (id - i * (ny * nz)) / nz;
				k = id - i * (ny * nz) - j * nz;

				mat_type = pt_glb->material_cell(id);
				if (mat_type != 0) {
					mat = &(pt_glb->material_parameters[mat_type - 1]);
					c11 = mat->c11; c12 = mat->c12; c44 = mat->c44;
				}

				if (k > pt_glb->strain_ext_nzi - 2 && k < pt_glb->strain_ext_nzf) {
					strain_in[0][0] = pt_glb->exx_ext_glb;
					strain_in[1][1] = pt_glb->eyy_ext_glb;
					strain_in[2][2] = -c12 / c11 * (pt_glb->exx_ext_glb + pt_glb->eyy_ext_glb);
					strain_in[2][1] = 0.;
					strain_in[1][2] = strain_in[2][1];
					strain_in[2][0] = 0.;
					strain_in[0][2] = strain_in[2][0];
					strain_in[1][0] = pt_glb->exy_ext_glb;
					strain_in[0][1] = strain_in[1][0];
				}
				else {
					strain_in[0][0] = 0.;
					strain_in[1][1] = 0.;
					strain_in[2][2] = 0.;
					strain_in[2][1] = 0.;
					strain_in[1][2] = 0.;
					strain_in[2][0] = 0.;
					strain_in[0][2] = 0.;
					strain_in[1][0] = 0.;
					strain_in[0][1] = 0.;
				}
			}
			pt_math->transform_matrix_glb2crt(strain_in, strain_out);

			exx_ext_crt(id) = strain_out[0][0];
			eyy_ext_crt(id) = strain_out[1][1];
			ezz_ext_crt(id) = strain_out[2][2];

			eyz_ext_crt(id) = strain_out[1][2];
			exz_ext_crt(id) = strain_out[0][2];
			exy_ext_crt(id) = strain_out[0][1];
		}
	}
}

void elastic_system::copy_to_device() {
	long int n21 = nx * ny * (nz / 2 + 1);

#pragma acc enter data copyin(this)

#pragma acc serial default(present)
	{
		this->pt_geo = &(geometry_parameters::geo);
		this->pt_glb = &(global_parameters::glb);
		this->pt_math = &(mathlib::mlb);
	}

#pragma acc enter data copyin(this->stress_free_surfYZ)
#pragma acc enter data copyin(this->stress_free_surfXZ)
#pragma acc enter data copyin(this->stress_free_surfXY)

#pragma acc enter data copyin(this->stress_free_surfYZ.matrix[0:(nx + 1)*ny*nz])
#pragma acc enter data copyin(this->stress_free_surfXZ.matrix[0:nx*(ny + 1)*nz])
#pragma acc enter data copyin(this->stress_free_surfXY.matrix[0:nx*ny*(nz + 1)])


#pragma acc enter data copyin(this->Gik_11,this->Gik_12,this->Gik_13)
#pragma acc enter data copyin(this->Gik_21,this->Gik_22,this->Gik_23)
#pragma acc enter data copyin(this->Gik_31,this->Gik_32,this->Gik_33)

#pragma acc enter data copyin(this->Gik_11.matrix[0:n21],this->Gik_12.matrix[0:n21],this->Gik_13.matrix[0:n21])
#pragma acc enter data copyin(this->Gik_21.matrix[0:n21],this->Gik_22.matrix[0:n21],this->Gik_23.matrix[0:n21])
#pragma acc enter data copyin(this->Gik_31.matrix[0:n21],this->Gik_32.matrix[0:n21],this->Gik_33.matrix[0:n21])


#pragma acc enter data copyin(this->duxdx_glb_surfYZ)
#pragma acc enter data copyin(this->duydx_glb_surfYZ)
#pragma acc enter data copyin(this->duzdx_glb_surfYZ)

#pragma acc enter data copyin(this->duxdy_glb_surfXZ)
#pragma acc enter data copyin(this->duydy_glb_surfXZ)
#pragma acc enter data copyin(this->duzdy_glb_surfXZ)

#pragma acc enter data copyin(this->duxdz_glb_surfXY)
#pragma acc enter data copyin(this->duydz_glb_surfXY)
#pragma acc enter data copyin(this->duzdz_glb_surfXY)

#pragma acc enter data copyin(this->duxdx_glb_surfYZ.matrix[0:(nx + 1)*ny*nz])
#pragma acc enter data copyin(this->duydx_glb_surfYZ.matrix[0:(nx + 1)*ny*nz])
#pragma acc enter data copyin(this->duzdx_glb_surfYZ.matrix[0:(nx + 1)*ny*nz])

#pragma acc enter data copyin(this->duxdy_glb_surfXZ.matrix[0:nx*(ny + 1)*nz])
#pragma acc enter data copyin(this->duydy_glb_surfXZ.matrix[0:nx*(ny + 1)*nz])
#pragma acc enter data copyin(this->duzdy_glb_surfXZ.matrix[0:nx*(ny + 1)*nz])

#pragma acc enter data copyin(this->duxdz_glb_surfXY.matrix[0:nx*ny*(nz + 1)])
#pragma acc enter data copyin(this->duydz_glb_surfXY.matrix[0:nx*ny*(nz + 1)])
#pragma acc enter data copyin(this->duzdz_glb_surfXY.matrix[0:nx*ny*(nz + 1)])

#pragma acc enter data copyin(this->Dux_glb,this->Duy_glb, this->Duz_glb)
#pragma acc enter data copyin(this->Dux_glb_store,this->Duy_glb_store, this->Duz_glb_store)
	//#pragma acc enter data copyin(this->dDux_glb,this->dDuy_glb, this->dDuz_glb)
#pragma acc enter data copyin(this->dDux_glb_rk1,this->dDuy_glb_rk1, this->dDuz_glb_rk1)
#pragma acc enter data copyin(this->dDux_glb_rk2,this->dDuy_glb_rk2, this->dDuz_glb_rk2)
#pragma acc enter data copyin(this->dDux_glb_rk3,this->dDuy_glb_rk3, this->dDuz_glb_rk3)
#pragma acc enter data copyin(this->dDux_glb_rk4,this->dDuy_glb_rk4, this->dDuz_glb_rk4)

#pragma acc enter data copyin(this->Dux_glb.matrix[0:n],this->Duy_glb.matrix[0:n], this->Duz_glb.matrix[0:n])
#pragma acc enter data copyin(this->Dux_glb_store.matrix[0:n],this->Duy_glb_store.matrix[0:n], this->Duz_glb_store.matrix[0:n])
//#pragma acc enter data copyin(this->dDux_glb.matrix[0:n],this->dDuy_glb.matrix[0:n], this->dDuz_glb.matrix[0:n])
#pragma acc enter data copyin(this->dDux_glb_rk1.matrix[0:n],this->dDuy_glb_rk1.matrix[0:n], this->dDuz_glb_rk1.matrix[0:n])
#pragma acc enter data copyin(this->dDux_glb_rk2.matrix[0:n],this->dDuy_glb_rk2.matrix[0:n], this->dDuz_glb_rk2.matrix[0:n])
#pragma acc enter data copyin(this->dDux_glb_rk3.matrix[0:n],this->dDuy_glb_rk3.matrix[0:n], this->dDuz_glb_rk3.matrix[0:n])
#pragma acc enter data copyin(this->dDux_glb_rk4.matrix[0:n],this->dDuy_glb_rk4.matrix[0:n], this->dDuz_glb_rk4.matrix[0:n])

#pragma acc enter data copyin(this->force_x,		this->force_y,			this->force_z)
#pragma acc enter data copyin(this->force_x_store,	this->force_y_store,	this->force_z_store)
#pragma acc enter data copyin(this->force_x.matrix[0:n],		this->force_y.matrix[0:n],			this->force_z.matrix[0:n])
#pragma acc enter data copyin(this->force_x_store.matrix[0:n],	this->force_y_store.matrix[0:n],	this->force_z_store.matrix[0:n])

#pragma acc enter data copyin(this->vx_glb,this->vy_glb, this->vz_glb)
#pragma acc enter data copyin(this->vx_glb_store,this->vy_glb_store, this->vz_glb_store)
//#pragma acc enter data copyin(this->dvx_glb,this->dvy_glb, this->dvz_glb)	
#pragma acc enter data copyin(this->dvx_glb_rk1,this->dvy_glb_rk1, this->dvz_glb_rk1)
#pragma acc enter data copyin(this->dvx_glb_rk2,this->dvy_glb_rk2, this->dvz_glb_rk2)
#pragma acc enter data copyin(this->dvx_glb_rk3,this->dvy_glb_rk3, this->dvz_glb_rk3)
#pragma acc enter data copyin(this->dvx_glb_rk4,this->dvy_glb_rk4, this->dvz_glb_rk4)

#pragma acc enter data copyin(this->vx_glb.matrix[0:n],this->vy_glb.matrix[0:n], this->vz_glb.matrix[0:n])
#pragma acc enter data copyin(this->vx_glb_store.matrix[0:n],this->vy_glb_store.matrix[0:n], this->vz_glb_store.matrix[0:n])
//#pragma acc enter data copyin(this->dvx_glb.matrix[0:n],this->dvy_glb.matrix[0:n], this->dvz_glb.matrix[0:n])	
#pragma acc enter data copyin(this->dvx_glb_rk1.matrix[0:n],this->dvy_glb_rk1.matrix[0:n], this->dvz_glb_rk1.matrix[0:n])
#pragma acc enter data copyin(this->dvx_glb_rk2.matrix[0:n],this->dvy_glb_rk2.matrix[0:n], this->dvz_glb_rk2.matrix[0:n])
#pragma acc enter data copyin(this->dvx_glb_rk3.matrix[0:n],this->dvy_glb_rk3.matrix[0:n], this->dvz_glb_rk3.matrix[0:n])
#pragma acc enter data copyin(this->dvx_glb_rk4.matrix[0:n],this->dvy_glb_rk4.matrix[0:n], this->dvz_glb_rk4.matrix[0:n])

#pragma acc enter data copyin(this->Dexx_glb)
#pragma acc enter data copyin(this->Deyy_glb)
#pragma acc enter data copyin(this->Dezz_glb)
#pragma acc enter data copyin(this->Deyz_glb)
#pragma acc enter data copyin(this->Dexz_glb)
#pragma acc enter data copyin(this->Dexy_glb)

#pragma acc enter data copyin(this->exxt0_glb)
#pragma acc enter data copyin(this->eyyt0_glb)
#pragma acc enter data copyin(this->ezzt0_glb)
#pragma acc enter data copyin(this->eyzt0_glb)
#pragma acc enter data copyin(this->exzt0_glb)
#pragma acc enter data copyin(this->exyt0_glb)

#pragma acc enter data copyin(this->Dexx_crt)
#pragma acc enter data copyin(this->Deyy_crt)
#pragma acc enter data copyin(this->Dezz_crt)
#pragma acc enter data copyin(this->Deyz_crt)
#pragma acc enter data copyin(this->Dexz_crt)
#pragma acc enter data copyin(this->Dexy_crt)

#pragma acc enter data copyin(this->exxt0_crt)
#pragma acc enter data copyin(this->eyyt0_crt)
#pragma acc enter data copyin(this->ezzt0_crt)
#pragma acc enter data copyin(this->eyzt0_crt)
#pragma acc enter data copyin(this->exzt0_crt)
#pragma acc enter data copyin(this->exyt0_crt)

#pragma acc enter data copyin(this->Dexx0_glb)
#pragma acc enter data copyin(this->Deyy0_glb)
#pragma acc enter data copyin(this->Dezz0_glb)
#pragma acc enter data copyin(this->Deyz0_glb)
#pragma acc enter data copyin(this->Dexz0_glb)
#pragma acc enter data copyin(this->Dexy0_glb)

#pragma acc enter data copyin(this->Dexx0_crt)
#pragma acc enter data copyin(this->Deyy0_crt)
#pragma acc enter data copyin(this->Dezz0_crt)
#pragma acc enter data copyin(this->Deyz0_crt)
#pragma acc enter data copyin(this->Dexz0_crt)
#pragma acc enter data copyin(this->Dexy0_crt)

#pragma acc enter data copyin(this->exx0t0_crt)
#pragma acc enter data copyin(this->eyy0t0_crt)
#pragma acc enter data copyin(this->ezz0t0_crt)
#pragma acc enter data copyin(this->eyz0t0_crt)
#pragma acc enter data copyin(this->exz0t0_crt)
#pragma acc enter data copyin(this->exy0t0_crt)

#pragma acc enter data copyin(this->Dexx_glb.matrix[0:n])
#pragma acc enter data copyin(this->Deyy_glb.matrix[0:n])
#pragma acc enter data copyin(this->Dezz_glb.matrix[0:n])
#pragma acc enter data copyin(this->Deyz_glb.matrix[0:n])
#pragma acc enter data copyin(this->Dexz_glb.matrix[0:n])
#pragma acc enter data copyin(this->Dexy_glb.matrix[0:n])

#pragma acc enter data copyin(this->exxt0_glb.matrix[0:n])
#pragma acc enter data copyin(this->eyyt0_glb.matrix[0:n])
#pragma acc enter data copyin(this->ezzt0_glb.matrix[0:n])
#pragma acc enter data copyin(this->eyzt0_glb.matrix[0:n])
#pragma acc enter data copyin(this->exzt0_glb.matrix[0:n])
#pragma acc enter data copyin(this->exyt0_glb.matrix[0:n])

#pragma acc enter data copyin(this->Dexx_crt.matrix[0:n])
#pragma acc enter data copyin(this->Deyy_crt.matrix[0:n])
#pragma acc enter data copyin(this->Dezz_crt.matrix[0:n])
#pragma acc enter data copyin(this->Deyz_crt.matrix[0:n])
#pragma acc enter data copyin(this->Dexz_crt.matrix[0:n])
#pragma acc enter data copyin(this->Dexy_crt.matrix[0:n])

#pragma acc enter data copyin(this->exxt0_crt.matrix[0:n])
#pragma acc enter data copyin(this->eyyt0_crt.matrix[0:n])
#pragma acc enter data copyin(this->ezzt0_crt.matrix[0:n])
#pragma acc enter data copyin(this->eyzt0_crt.matrix[0:n])
#pragma acc enter data copyin(this->exzt0_crt.matrix[0:n])
#pragma acc enter data copyin(this->exyt0_crt.matrix[0:n])

#pragma acc enter data copyin(this->Dexx0_glb.matrix[0:n])
#pragma acc enter data copyin(this->Deyy0_glb.matrix[0:n])
#pragma acc enter data copyin(this->Dezz0_glb.matrix[0:n])
#pragma acc enter data copyin(this->Deyz0_glb.matrix[0:n])
#pragma acc enter data copyin(this->Dexz0_glb.matrix[0:n])
#pragma acc enter data copyin(this->Dexy0_glb.matrix[0:n])

#pragma acc enter data copyin(this->Dexx0_crt.matrix[0:n])
#pragma acc enter data copyin(this->Deyy0_crt.matrix[0:n])
#pragma acc enter data copyin(this->Dezz0_crt.matrix[0:n])
#pragma acc enter data copyin(this->Deyz0_crt.matrix[0:n])
#pragma acc enter data copyin(this->Dexz0_crt.matrix[0:n])
#pragma acc enter data copyin(this->Dexy0_crt.matrix[0:n])

#pragma acc enter data copyin(this->exx0t0_crt.matrix[0:n])
#pragma acc enter data copyin(this->eyy0t0_crt.matrix[0:n])
#pragma acc enter data copyin(this->ezz0t0_crt.matrix[0:n])
#pragma acc enter data copyin(this->eyz0t0_crt.matrix[0:n])
#pragma acc enter data copyin(this->exz0t0_crt.matrix[0:n])
#pragma acc enter data copyin(this->exy0t0_crt.matrix[0:n])

#pragma acc enter data copyin(this->Dux_glb_t1,this->Duy_glb_t1,this->Duz_glb_t1)
#pragma acc enter data copyin(this->Dux_glb_t1.matrix[0:nx*ny*nz],this->Duy_glb_t1.matrix[0:nx*ny*nz],this->Duz_glb_t1.matrix[0:nx*ny*nz])

	if (pt_glb->if_1D_ABC == false) {
#pragma acc enter data copyin(this->Dux_glb_t2,this->Duy_glb_t2,this->Duz_glb_t2)
#pragma acc enter data copyin(this->Dux_glb_t3,this->Duy_glb_t3,this->Duz_glb_t3)
#pragma acc enter data copyin(this->Dux_glb_t4,this->Duy_glb_t4,this->Duz_glb_t4)

#pragma acc enter data copyin(this->Dux_glb_t2.matrix[0:nx*ny*nz],this->Duy_glb_t2.matrix[0:nx*ny*nz],this->Duz_glb_t2.matrix[0:nx*ny*nz])
#pragma acc enter data copyin(this->Dux_glb_t3.matrix[0:nx*ny*nz],this->Duy_glb_t3.matrix[0:nx*ny*nz],this->Duz_glb_t3.matrix[0:nx*ny*nz])
#pragma acc enter data copyin(this->Dux_glb_t4.matrix[0:nx*ny*nz],this->Duy_glb_t4.matrix[0:nx*ny*nz],this->Duz_glb_t4.matrix[0:nx*ny*nz])
	}

#pragma acc enter data copyin(exx_ext_crt, eyy_ext_crt, ezz_ext_crt)
#pragma acc enter data copyin(exy_ext_crt, exz_ext_crt, eyz_ext_crt)

#pragma acc enter data copyin(exx_ext_crt.matrix[0:n], eyy_ext_crt.matrix[0:n], ezz_ext_crt.matrix[0:n])
#pragma acc enter data copyin(exy_ext_crt.matrix[0:n], exz_ext_crt.matrix[0:n], eyz_ext_crt.matrix[0:n])

#pragma acc enter data copyin(stress0homo_xx , stress0homo_yy , stress0homo_zz )
#pragma acc enter data copyin(stress0homo_xy , stress0homo_xz , stress0homo_yz )

#pragma acc enter data copyin(stress0_xx , stress0_yy , stress0_zz )
#pragma acc enter data copyin(stress0_xy , stress0_xz , stress0_yz )

#pragma acc enter data copyin(stress_xx , stress_yy , stress_zz )
#pragma acc enter data copyin(stress_xy , stress_xz , stress_yz )

#pragma acc enter data copyin(duxdx , duxdy , duxdz )
#pragma acc enter data copyin(duydx , duydy , duydz )
#pragma acc enter data copyin(duzdx , duzdy , duzdz )

#pragma acc enter data copyin(stress0homok_xx , stress0homok_yy , stress0homok_zz )
#pragma acc enter data copyin(stress0homok_xy , stress0homok_xz , stress0homok_yz )

#pragma acc enter data copyin(stressk_xx , stressk_yy , stressk_zz )
#pragma acc enter data copyin(stressk_xy , stressk_xz , stressk_yz )

#pragma acc enter data copyin(duxdxk , duxdyk , duxdzk )
#pragma acc enter data copyin(duydxk , duydyk , duydzk )
#pragma acc enter data copyin(duzdxk , duzdyk , duzdzk )

#pragma acc enter data copyin(stress0homo_xx.matrix[0:n], stress0homo_yy.matrix[0:n], stress0homo_zz.matrix[0:n])
#pragma acc enter data copyin(stress0homo_xy.matrix[0:n], stress0homo_xz.matrix[0:n], stress0homo_yz.matrix[0:n])

#pragma acc enter data copyin(stress0_xx.matrix[0:n], stress0_yy.matrix[0:n], stress0_zz.matrix[0:n])
#pragma acc enter data copyin(stress0_xy.matrix[0:n], stress0_xz.matrix[0:n], stress0_yz.matrix[0:n])

#pragma acc enter data copyin(stress_xx.matrix[0:n], stress_yy.matrix[0:n], stress_zz.matrix[0:n])
#pragma acc enter data copyin(stress_xy.matrix[0:n], stress_xz.matrix[0:n], stress_yz.matrix[0:n])

#pragma acc enter data copyin(duxdx.matrix[0:n], duxdy.matrix[0:n], duxdz.matrix[0:n])
#pragma acc enter data copyin(duydx.matrix[0:n], duydy.matrix[0:n], duydz.matrix[0:n])
#pragma acc enter data copyin(duzdx.matrix[0:n], duzdy.matrix[0:n], duzdz.matrix[0:n])

#pragma acc enter data copyin(stress0homok_xx.matrix[0:n21], stress0homok_yy.matrix[0:n21], stress0homok_zz.matrix[0:n21])
#pragma acc enter data copyin(stress0homok_xy.matrix[0:n21], stress0homok_xz.matrix[0:n21], stress0homok_yz.matrix[0:n21])

#pragma acc enter data copyin(stressk_xx.matrix[0:n21], stressk_yy.matrix[0:n21], stressk_zz.matrix[0:n21])
#pragma acc enter data copyin(stressk_xy.matrix[0:n21], stressk_xz.matrix[0:n21], stressk_yz.matrix[0:n21])

#pragma acc enter data copyin(duxdxk.matrix[0:n21], duxdyk.matrix[0:n21], duxdzk.matrix[0:n21])
#pragma acc enter data copyin(duydxk.matrix[0:n21], duydyk.matrix[0:n21], duydzk.matrix[0:n21])
#pragma acc enter data copyin(duzdxk.matrix[0:n21], duzdyk.matrix[0:n21], duzdzk.matrix[0:n21])

	if (pt_glb->if_elastostatic_film == true) {
#pragma acc enter data copyin(exx0_ave_film[0:nz],eyy0_ave_film[0:nz],ezz0_ave_film[0:nz])
#pragma acc enter data copyin(eyz0_ave_film[0:nz],exz0_ave_film[0:nz],exy0_ave_film[0:nz])	
	}

	if (pt_glb->if_FE_all == true && pt_glb->if_flexo == true) {
#pragma acc enter data copyin(dexxt0dx_glb, dexxt0dy_glb, dexxt0dz_glb)
#pragma acc enter data copyin(deyyt0dx_glb, deyyt0dy_glb, deyyt0dz_glb)
#pragma acc enter data copyin(dezzt0dx_glb, dezzt0dy_glb, dezzt0dz_glb)
#pragma acc enter data copyin(deyzt0dx_glb, deyzt0dy_glb, deyzt0dz_glb)
#pragma acc enter data copyin(dexzt0dx_glb, dexzt0dy_glb, dexzt0dz_glb)
#pragma acc enter data copyin(dexyt0dx_glb, dexyt0dy_glb, dexyt0dz_glb)

#pragma acc enter data copyin(dDexxdx_glb, dDexxdy_glb, dDexxdz_glb)
#pragma acc enter data copyin(dDeyydx_glb, dDeyydy_glb, dDeyydz_glb)
#pragma acc enter data copyin(dDezzdx_glb, dDezzdy_glb, dDezzdz_glb)
#pragma acc enter data copyin(dDeyzdx_glb, dDeyzdy_glb, dDeyzdz_glb)
#pragma acc enter data copyin(dDexzdx_glb, dDexzdy_glb, dDexzdz_glb)
#pragma acc enter data copyin(dDexydx_glb, dDexydy_glb, dDexydz_glb)

#pragma acc enter data copyin(dexxt0dx_glb.matrix[0:n], dexxt0dy_glb.matrix[0:n], dexxt0dz_glb.matrix[0:n])
#pragma acc enter data copyin(deyyt0dx_glb.matrix[0:n], deyyt0dy_glb.matrix[0:n], deyyt0dz_glb.matrix[0:n])
#pragma acc enter data copyin(dezzt0dx_glb.matrix[0:n], dezzt0dy_glb.matrix[0:n], dezzt0dz_glb.matrix[0:n])
#pragma acc enter data copyin(deyzt0dx_glb.matrix[0:n], deyzt0dy_glb.matrix[0:n], deyzt0dz_glb.matrix[0:n])
#pragma acc enter data copyin(dexzt0dx_glb.matrix[0:n], dexzt0dy_glb.matrix[0:n], dexzt0dz_glb.matrix[0:n])
#pragma acc enter data copyin(dexyt0dx_glb.matrix[0:n], dexyt0dy_glb.matrix[0:n], dexyt0dz_glb.matrix[0:n])

#pragma acc enter data copyin(dDexxdx_glb.matrix[0:n], dDexxdy_glb.matrix[0:n], dDexxdz_glb.matrix[0:n])
#pragma acc enter data copyin(dDeyydx_glb.matrix[0:n], dDeyydy_glb.matrix[0:n], dDeyydz_glb.matrix[0:n])
#pragma acc enter data copyin(dDezzdx_glb.matrix[0:n], dDezzdy_glb.matrix[0:n], dDezzdz_glb.matrix[0:n])
#pragma acc enter data copyin(dDeyzdx_glb.matrix[0:n], dDeyzdy_glb.matrix[0:n], dDeyzdz_glb.matrix[0:n])
#pragma acc enter data copyin(dDexzdx_glb.matrix[0:n], dDexzdy_glb.matrix[0:n], dDexzdz_glb.matrix[0:n])
#pragma acc enter data copyin(dDexydx_glb.matrix[0:n], dDexydy_glb.matrix[0:n], dDexydz_glb.matrix[0:n])
	}

}

void elastic_system::copy_Dstrain_from_device() {
#pragma acc update host(this->Dexx_glb.matrix[0:n])
#pragma acc update host(this->Deyy_glb.matrix[0:n])
#pragma acc update host(this->Dezz_glb.matrix[0:n])
#pragma acc update host(this->Deyz_glb.matrix[0:n])
#pragma acc update host(this->Dexz_glb.matrix[0:n])
#pragma acc update host(this->Dexy_glb.matrix[0:n])
}

void elastic_system::copy_straint0_from_device() {
#pragma acc update host(this->exxt0_glb.matrix[0:n])
#pragma acc update host(this->eyyt0_glb.matrix[0:n])
#pragma acc update host(this->ezzt0_glb.matrix[0:n])
#pragma acc update host(this->eyzt0_glb.matrix[0:n])
#pragma acc update host(this->exzt0_glb.matrix[0:n])
#pragma acc update host(this->exyt0_glb.matrix[0:n])
}

void elastic_system::copy_uandv_from_device() {
#pragma acc update host(this->Dux_glb.matrix[0:n],this->Duy_glb.matrix[0:n], this->Duz_glb.matrix[0:n])
#pragma acc update host(this->vx_glb.matrix[0:n],this->vy_glb.matrix[0:n], this->vz_glb.matrix[0:n])
}

void elastic_system::copy_elastoforce_from_device() {
#pragma acc update host(this->force_x_store.matrix[0:n],this->force_y_store.matrix[0:n], this->force_z_store.matrix[0:n])
}

void elastic_system::copy_eigenstraint0_from_device() {
#pragma acc update host(this->exx0t0_crt.matrix[0:n])
#pragma acc update host(this->eyy0t0_crt.matrix[0:n])
#pragma acc update host(this->ezz0t0_crt.matrix[0:n])
#pragma acc update host(this->eyz0t0_crt.matrix[0:n])
#pragma acc update host(this->exz0t0_crt.matrix[0:n])
#pragma acc update host(this->exy0t0_crt.matrix[0:n])
}

void elastic_system::read_external_strain() {
	long x = 0;
	long y = 0;
	long z = 0;
	double inxx, inyy, inzz, inyz, inxz, inxy;

	std::ifstream file("strain_ext.in");
	if (file.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			file >> x >> y >> z >> inxx >> inyy >> inzz >> inyz >> inxz >> inxy;
			exx_ext_crt(x - 1, y - 1, z - 1) = inxx;
			eyy_ext_crt(x - 1, y - 1, z - 1) = inyy;
			ezz_ext_crt(x - 1, y - 1, z - 1) = inzz;
			eyz_ext_crt(x - 1, y - 1, z - 1) = inyz;
			exz_ext_crt(x - 1, y - 1, z - 1) = inxz;
			exy_ext_crt(x - 1, y - 1, z - 1) = inxy;
		}
	}
	file.close();
}

