#include "EMdynamic_system.h"

void EMdynamic_system::initialize_host(){

//	pt_mag=pt_mag_in;
//	pt_ele=pt_ele_in;

	nx = pt_geo->nx_system;
	ny = pt_geo->ny_system;
	nz = pt_geo->nz_system;

	n = nx * ny * nz;

	dx = pt_geo->dx;
	dy = pt_geo->dy;
	dz = pt_geo->dz;

	planeEM_source_z = pt_geo->pt_nz_layer[0] / 5;

	int i, j;
	unsigned int mat_type;
	material* mat;

	sqrt_er_bot = 1.;
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (pt_glb->material_cell(i, j, 0) != 0) {
				mat_type = pt_glb->material_cell(i, j, 0);
				mat = &(pt_glb->material_parameters[mat_type - 1]);
				sqrt_er_bot = sqrt(mat->r_permittivity);
				j = ny; i = nx;
			}
		}
	}

	sqrt_er_top = 1.;
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (pt_glb->material_cell(i, j, nz - 1) != 0) {
				mat_type = pt_glb->material_cell(i, j, nz - 1);
				mat = &(pt_glb->material_parameters[mat_type - 1]);
				sqrt_er_top = sqrt(mat->r_permittivity);
				j = ny; i = nx;
			}
		}
	}

	coeff_top_half = (2. * dz - c / sqrt_er_top * pt_glb->dt) / (2. * dz + c / sqrt_er_top * pt_glb->dt);
	coeff_bot_half = (2. * dz - c / sqrt_er_bot * pt_glb->dt) / (2. * dz + c / sqrt_er_bot * pt_glb->dt);
	coeff_x_half = (2. * dx - c * pt_glb->dt) / (2. * dx + c * pt_glb->dt);
	coeff_y_half = (2. * dy - c * pt_glb->dt) / (2. * dy + c * pt_glb->dt);

	coeff_top = (dz - c / sqrt_er_top * pt_glb->dt) / (dz + c / sqrt_er_top * pt_glb->dt);
	coeff_bot = (dz - c / sqrt_er_bot * pt_glb->dt) / (dz + c / sqrt_er_bot * pt_glb->dt);
	coeff_x = (dx - c * pt_glb->dt) / (dx + c * pt_glb->dt);
	coeff_y = (dy - c * pt_glb->dt) / (dy + c * pt_glb->dt);

	w = pt_glb->weighting_factor_fourth_Liao;

	if (pt_glb->if_Jf_input == true) {
		Jf_nx = (pt_glb->Jfin_xf - pt_glb->Jfin_xi + 1);
		Jf_ny = (pt_glb->Jfin_yf - pt_glb->Jfin_yi + 1);
		Jf_nz = (pt_glb->Jfin_zf - pt_glb->Jfin_zi + 1);
		Jf_n = Jf_nx * Jf_ny * Jf_nz;
		//Jf_z1 = pt_geo->pt_nz_layer[0];
		//Jf_z2 = pt_geo->pt_nz_layer[0]+ pt_geo->pt_nz_layer[1];
		//Jf_thickness = static_cast<double>(pt_geo->pt_nz_layer[1]) * dz;
	}

	//Initialization
	DHx_em.initialize(nx + 1, ny, nz); DHx_em_store.initialize(nx + 1, ny, nz);
	DHy_em.initialize(nx, ny + 1, nz); DHy_em_store.initialize(nx, ny + 1, nz);
	DHz_em.initialize(nx, ny, nz + 1); DHz_em_store.initialize(nx, ny, nz + 1);

	dDHx_em_rk1.initialize(nx + 1, ny, nz);
	dDHy_em_rk1.initialize(nx, ny + 1, nz);
	dDHz_em_rk1.initialize(nx, ny, nz + 1);

	dDHx_em_rk2.initialize(nx + 1, ny, nz);
	dDHy_em_rk2.initialize(nx, ny + 1, nz);
	dDHz_em_rk2.initialize(nx, ny, nz + 1);

	dDHx_em_rk3.initialize(nx + 1, ny, nz);
	dDHy_em_rk3.initialize(nx, ny + 1, nz);
	dDHz_em_rk3.initialize(nx, ny, nz + 1);

	dDHx_em_rk4.initialize(nx + 1, ny, nz);
	dDHy_em_rk4.initialize(nx, ny + 1, nz);
	dDHz_em_rk4.initialize(nx, ny, nz + 1);

	DEx_em.initialize(nx, ny + 1, nz + 1); DEx_em_store.initialize(nx, ny + 1, nz + 1);
	DEy_em.initialize(nx + 1, ny, nz + 1); DEy_em_store.initialize(nx + 1, ny, nz + 1);
	DEz_em.initialize(nx + 1, ny + 1, nz); DEz_em_store.initialize(nx + 1, ny + 1, nz);

	dDEx_em_rk1.initialize(nx, ny + 1, nz + 1);
	dDEy_em_rk1.initialize(nx + 1, ny, nz + 1);
	dDEz_em_rk1.initialize(nx + 1, ny + 1, nz);

	dDEx_em_rk2.initialize(nx, ny + 1, nz + 1);
	dDEy_em_rk2.initialize(nx + 1, ny, nz + 1);
	dDEz_em_rk2.initialize(nx + 1, ny + 1, nz);

	dDEx_em_rk3.initialize(nx, ny + 1, nz + 1);
	dDEy_em_rk3.initialize(nx + 1, ny, nz + 1);
	dDEz_em_rk3.initialize(nx + 1, ny + 1, nz);

	dDEx_em_rk4.initialize(nx, ny + 1, nz + 1);
	dDEy_em_rk4.initialize(nx + 1, ny, nz + 1);
	dDEz_em_rk4.initialize(nx + 1, ny + 1, nz);

	DHx_em_cell.initialize(nx, ny, nz); DHy_em_cell.initialize(nx, ny, nz); DHz_em_cell.initialize(nx, ny, nz);
	DEx_em_cell.initialize(nx, ny, nz); DEy_em_cell.initialize(nx, ny, nz); DEz_em_cell.initialize(nx, ny, nz);

	DJfx.initialize(nx, ny, nz); DJfy.initialize(nx, ny, nz); DJfz.initialize(nx, ny, nz);

	Jpx_n1.initialize(nx, ny, nz); Jpy_n1.initialize(nx, ny, nz); Jpz_n1.initialize(nx, ny, nz);
	Jpx_n2.initialize(nx, ny, nz); Jpy_n2.initialize(nx, ny, nz); Jpz_n2.initialize(nx, ny, nz);
	Jpx_n3.initialize(nx, ny, nz); Jpy_n3.initialize(nx, ny, nz); Jpz_n3.initialize(nx, ny, nz);

	Jpx_n1_store.initialize(nx, ny, nz); Jpy_n1_store.initialize(nx, ny, nz); Jpz_n1_store.initialize(nx, ny, nz);
	Jpx_n2_store.initialize(nx, ny, nz); Jpy_n2_store.initialize(nx, ny, nz); Jpz_n2_store.initialize(nx, ny, nz);
	Jpx_n3_store.initialize(nx, ny, nz); Jpy_n3_store.initialize(nx, ny, nz); Jpz_n3_store.initialize(nx, ny, nz);

	dJpx_n1_rk1.initialize(nx, ny, nz); dJpy_n1_rk1.initialize(nx, ny, nz); dJpz_n1_rk1.initialize(nx, ny, nz);
	dJpx_n2_rk1.initialize(nx, ny, nz); dJpy_n2_rk1.initialize(nx, ny, nz); dJpz_n2_rk1.initialize(nx, ny, nz);
	dJpx_n3_rk1.initialize(nx, ny, nz); dJpy_n3_rk1.initialize(nx, ny, nz); dJpz_n3_rk1.initialize(nx, ny, nz);

	dJpx_n1_rk2.initialize(nx, ny, nz); dJpy_n1_rk2.initialize(nx, ny, nz); dJpz_n1_rk2.initialize(nx, ny, nz);
	dJpx_n2_rk2.initialize(nx, ny, nz); dJpy_n2_rk2.initialize(nx, ny, nz); dJpz_n2_rk2.initialize(nx, ny, nz);
	dJpx_n3_rk2.initialize(nx, ny, nz); dJpy_n3_rk2.initialize(nx, ny, nz); dJpz_n3_rk2.initialize(nx, ny, nz);

	dJpx_n1_rk3.initialize(nx, ny, nz); dJpy_n1_rk3.initialize(nx, ny, nz); dJpz_n1_rk3.initialize(nx, ny, nz);
	dJpx_n2_rk3.initialize(nx, ny, nz); dJpy_n2_rk3.initialize(nx, ny, nz); dJpz_n2_rk3.initialize(nx, ny, nz);
	dJpx_n3_rk3.initialize(nx, ny, nz); dJpy_n3_rk3.initialize(nx, ny, nz); dJpz_n3_rk3.initialize(nx, ny, nz);

	dJpx_n1_rk4.initialize(nx, ny, nz); dJpy_n1_rk4.initialize(nx, ny, nz); dJpz_n1_rk4.initialize(nx, ny, nz);
	dJpx_n2_rk4.initialize(nx, ny, nz); dJpy_n2_rk4.initialize(nx, ny, nz); dJpz_n2_rk4.initialize(nx, ny, nz);
	dJpx_n3_rk4.initialize(nx, ny, nz); dJpy_n3_rk4.initialize(nx, ny, nz); dJpz_n3_rk4.initialize(nx, ny, nz);

	if (pt_glb->if_periodic_allsurface == false) {
		DEx_em_t1.initialize(nx, ny + 1, nz + 1); DEy_em_t1.initialize(nx + 1, ny, nz + 1); DEz_em_t1.initialize(nx + 1, ny + 1, nz);
		if (pt_glb->if_1D_ABC == false) {
			DEx_em_t2.initialize(nx, ny + 1, nz + 1); DEy_em_t2.initialize(nx + 1, ny, nz + 1); DEz_em_t2.initialize(nx + 1, ny + 1, nz);
			DEx_em_t3.initialize(nx, ny + 1, nz + 1); DEy_em_t3.initialize(nx + 1, ny, nz + 1); DEz_em_t3.initialize(nx + 1, ny + 1, nz);
			DEx_em_t4.initialize(nx, ny + 1, nz + 1); DEy_em_t4.initialize(nx + 1, ny, nz + 1); DEz_em_t4.initialize(nx + 1, ny + 1, nz);
		}
	}

	if (pt_glb->if_PEC_all == true) {
		DEx_ifPEC.initialize(nx, ny + 1, nz + 1);
		DEy_ifPEC.initialize(nx + 1, ny, nz + 1);
		DEz_ifPEC.initialize(nx + 1, ny + 1, nz);

		initialize_PEC_edges();
	}
}

void EMdynamic_system::initialize_PEC_edges() {
	long int i, j, k;
	unsigned int mat_type;
	material* mat;
	long int idx1, idy1, idz1;
	long int idx2, idy2, idz2;
	long int idx3, idy3, idz3;
	long int idx4, idy4, idz4;

	//-----------------X component-------------------//
	for (long int id = 0; id < nx * (ny + 1) * (nz + 1); id++) {
		DEx_ifPEC(id) = false;

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

		//----------Cell 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (mat->if_PEC == true) {
				DEx_ifPEC(id) = true;
			}
		}

		//----------Cell 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (mat->if_PEC == true) {
				DEx_ifPEC(id) = true;
			}
		}

		//----------Cell 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (mat->if_PEC == true) {
				DEx_ifPEC(id) = true;
			}
		}

		//----------Cell 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (mat->if_PEC == true) {
				DEx_ifPEC(id) = true;
			}
		}
	}

	//-----------------Y component-------------------//
	for (long int id = 0; id < (nx + 1) * ny * (nz + 1); id++) {
		DEy_ifPEC(id) = false;

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

		//----------Cell 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (mat->if_PEC == true) {
				DEy_ifPEC(id) = true;
			}
		}

		//----------Cell 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (mat->if_PEC == true) {
				DEy_ifPEC(id) = true;
			}
		}

		//----------Cell 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (mat->if_PEC == true) {
				DEy_ifPEC(id) = true;
			}
		}

		//----------Cell 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);	
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (mat->if_PEC == true) {
				DEy_ifPEC(id) = true;
			}
		}
	}

	//-----------------Z component-------------------//
	for (long int id = 0; id < (nx + 1) * (ny + 1) * nz; id++) {
		DEz_ifPEC(id) = false;

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

		//----------Cell 1------//
		mat_type = pt_glb->material_cell(idx1, idy1, idz1);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (mat->if_PEC == true) {
				DEz_ifPEC(id) = true;
			}
		}

		//----------Cell 2------//
		mat_type = pt_glb->material_cell(idx2, idy2, idz2);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (mat->if_PEC == true) {
				DEz_ifPEC(id) = true;
			}
		}

		//----------Cell 3------//
		mat_type = pt_glb->material_cell(idx3, idy3, idz3);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (mat->if_PEC == true) {
				DEz_ifPEC(id) = true;
			}
		}

		//----------Cell 4------//
		mat_type = pt_glb->material_cell(idx4, idy4, idz4);
		if (mat_type != 0) {
			mat = &(pt_glb->material_parameters[(mat_type)-1]);
			if (mat->if_PEC == true) {
				DEz_ifPEC(id) = true;
			}
		}
	}
}

void EMdynamic_system::initialize_device() {
	if (pt_glb->if_input_em == true) {

#pragma acc parallel default(present) async(3)
		{
			DHx_em_store = DHx_em;
			DHy_em_store = DHy_em;
			DHz_em_store = DHz_em;
			DEx_em_store = DEx_em;
			DEy_em_store = DEy_em;
			DEz_em_store = DEz_em;
		}

#pragma acc parallel default(present) async(3)
		{
			update_DH_cell();
			update_DE_cell();
		}
	}

	if (pt_glb->if_input_Jp == true) {
#pragma acc parallel default(present) async(4)
		{
			Jpx_n1_store = Jpx_n1;
			Jpy_n1_store = Jpy_n1;
			Jpz_n1_store = Jpz_n1;
			Jpx_n2_store = Jpx_n2;
			Jpy_n2_store = Jpy_n2;
			Jpz_n2_store = Jpz_n2;
			Jpx_n3_store = Jpx_n3;
			Jpy_n3_store = Jpy_n3;
			Jpz_n3_store = Jpz_n3;
		}
	}
}

void EMdynamic_system::copy_to_device() {
#pragma acc enter data copyin(this)

#pragma acc serial default(present)
	{
		this->pt_geo = &(geometry_parameters::geo);
		this->pt_glb = &(global_parameters::glb);
		this->pt_math = &(mathlib::mlb);
	}

#pragma acc enter data copyin(this->DHx_em, this->DHx_em_store)
#pragma acc enter data copyin(this->DHy_em, this->DHy_em_store)
#pragma acc enter data copyin(this->DHz_em, this->DHz_em_store)

#pragma acc enter data copyin(this->DHx_em.matrix[0:(nx+1)*ny*nz], this->DHx_em_store.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->DHy_em.matrix[0:nx*(ny+1)*nz], this->DHy_em_store.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->DHz_em.matrix[0:nx*ny*(nz+1)], this->DHz_em_store.matrix[0:nx*ny*(nz+1)])

#pragma acc enter data copyin(this->DEx_em, this->DEx_em_store)
#pragma acc enter data copyin(this->DEy_em, this->DEy_em_store)
#pragma acc enter data copyin(this->DEz_em, this->DEz_em_store)

#pragma acc enter data copyin(this->DEx_em.matrix[0:nx*(ny+1)*(nz+1)], this->DEx_em_store.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->DEy_em.matrix[0:(nx+1)*ny*(nz+1)], this->DEy_em_store.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->DEz_em.matrix[0:(nx+1)*(ny+1)*nz], this->DEz_em_store.matrix[0:(nx+1)*(ny+1)*nz])

	if (pt_glb->if_PEC_all == true) {
#pragma acc enter data copyin(this->DEx_ifPEC)
#pragma acc enter data copyin(this->DEy_ifPEC)
#pragma acc enter data copyin(this->DEz_ifPEC)

#pragma acc enter data copyin(this->DEx_ifPEC.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->DEy_ifPEC.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->DEz_ifPEC.matrix[0:(nx+1)*(ny+1)*nz])
	}

#pragma acc enter data copyin(this->dDHx_em_rk1)
#pragma acc enter data copyin(this->dDHy_em_rk1)
#pragma acc enter data copyin(this->dDHz_em_rk1)
#pragma acc enter data copyin(this->dDHx_em_rk1.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->dDHy_em_rk1.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->dDHz_em_rk1.matrix[0:nx*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEx_em_rk1)
#pragma acc enter data copyin(this->dDEy_em_rk1)
#pragma acc enter data copyin(this->dDEz_em_rk1)
#pragma acc enter data copyin(this->dDEx_em_rk1.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->dDEy_em_rk1.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEz_em_rk1.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->dDHx_em_rk2)
#pragma acc enter data copyin(this->dDHy_em_rk2)
#pragma acc enter data copyin(this->dDHz_em_rk2)
#pragma acc enter data copyin(this->dDHx_em_rk2.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->dDHy_em_rk2.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->dDHz_em_rk2.matrix[0:nx*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEx_em_rk2)
#pragma acc enter data copyin(this->dDEy_em_rk2)
#pragma acc enter data copyin(this->dDEz_em_rk2)
#pragma acc enter data copyin(this->dDEx_em_rk2.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->dDEy_em_rk2.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEz_em_rk2.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->dDHx_em_rk3)
#pragma acc enter data copyin(this->dDHy_em_rk3)
#pragma acc enter data copyin(this->dDHz_em_rk3)
#pragma acc enter data copyin(this->dDHx_em_rk3.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->dDHy_em_rk3.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->dDHz_em_rk3.matrix[0:nx*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEx_em_rk3)
#pragma acc enter data copyin(this->dDEy_em_rk3)
#pragma acc enter data copyin(this->dDEz_em_rk3)
#pragma acc enter data copyin(this->dDEx_em_rk3.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->dDEy_em_rk3.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEz_em_rk3.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->dDHx_em_rk4)
#pragma acc enter data copyin(this->dDHy_em_rk4)
#pragma acc enter data copyin(this->dDHz_em_rk4)
#pragma acc enter data copyin(this->dDHx_em_rk4.matrix[0:(nx+1)*ny*nz])
#pragma acc enter data copyin(this->dDHy_em_rk4.matrix[0:nx*(ny+1)*nz])
#pragma acc enter data copyin(this->dDHz_em_rk4.matrix[0:nx*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEx_em_rk4)
#pragma acc enter data copyin(this->dDEy_em_rk4)
#pragma acc enter data copyin(this->dDEz_em_rk4)
#pragma acc enter data copyin(this->dDEx_em_rk4.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->dDEy_em_rk4.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->dDEz_em_rk4.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->DHx_em_cell,this->DHy_em_cell,this->DHz_em_cell)
#pragma acc enter data copyin(this->DEx_em_cell,this->DEy_em_cell,this->DEz_em_cell)

#pragma acc enter data copyin(this->DHx_em_cell.matrix[0:n],this->DHy_em_cell.matrix[0:n],this->DHz_em_cell.matrix[0:n])
#pragma acc enter data copyin(this->DEx_em_cell.matrix[0:n],this->DEy_em_cell.matrix[0:n],this->DEz_em_cell.matrix[0:n])

#pragma acc enter data copyin(this->DJfx,this->DJfy,this->DJfz)
#pragma acc enter data copyin(this->DJfx.matrix[0:n],this->DJfy.matrix[0:n],this->DJfz.matrix[0:n])

#pragma acc enter data copyin(this->Jpx_n1,this->Jpy_n1,this->Jpz_n1)
#pragma acc enter data copyin(this->Jpx_n2,this->Jpy_n2,this->Jpz_n2)
#pragma acc enter data copyin(this->Jpx_n3,this->Jpy_n3,this->Jpz_n3)
#pragma acc enter data copyin(this->Jpx_n1_store,this->Jpy_n1_store,this->Jpz_n1_store)
#pragma acc enter data copyin(this->Jpx_n2_store,this->Jpy_n2_store,this->Jpz_n2_store)
#pragma acc enter data copyin(this->Jpx_n3_store,this->Jpy_n3_store,this->Jpz_n3_store)

#pragma acc enter data copyin(this->Jpx_n1.matrix[0:n],this->Jpy_n1.matrix[0:n],this->Jpz_n1.matrix[0:n])
#pragma acc enter data copyin(this->Jpx_n2.matrix[0:n],this->Jpy_n2.matrix[0:n],this->Jpz_n2.matrix[0:n])
#pragma acc enter data copyin(this->Jpx_n3.matrix[0:n],this->Jpy_n3.matrix[0:n],this->Jpz_n3.matrix[0:n])
#pragma acc enter data copyin(this->Jpx_n1_store.matrix[0:n],this->Jpy_n1_store.matrix[0:n],this->Jpz_n1_store.matrix[0:n])
#pragma acc enter data copyin(this->Jpx_n2_store.matrix[0:n],this->Jpy_n2_store.matrix[0:n],this->Jpz_n2_store.matrix[0:n])
#pragma acc enter data copyin(this->Jpx_n3_store.matrix[0:n],this->Jpy_n3_store.matrix[0:n],this->Jpz_n3_store.matrix[0:n])

#pragma acc enter data copyin(this->dJpx_n1_rk1,			 this->dJpy_n1_rk1,				this->dJpz_n1_rk1)
#pragma acc enter data copyin(this->dJpx_n2_rk1,			 this->dJpy_n2_rk1,				this->dJpz_n2_rk1)
#pragma acc enter data copyin(this->dJpx_n3_rk1,			 this->dJpy_n3_rk1,				this->dJpz_n3_rk1)
#pragma acc enter data copyin(this->dJpx_n1_rk1.matrix[0:n], this->dJpy_n1_rk1.matrix[0:n], this->dJpz_n1_rk1.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n2_rk1.matrix[0:n], this->dJpy_n2_rk1.matrix[0:n], this->dJpz_n2_rk1.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n3_rk1.matrix[0:n], this->dJpy_n3_rk1.matrix[0:n], this->dJpz_n3_rk1.matrix[0:n])

#pragma acc enter data copyin(this->dJpx_n1_rk2,			 this->dJpy_n1_rk2,				this->dJpz_n1_rk2)
#pragma acc enter data copyin(this->dJpx_n2_rk2,			 this->dJpy_n2_rk2,				this->dJpz_n2_rk2)
#pragma acc enter data copyin(this->dJpx_n3_rk2,			 this->dJpy_n3_rk2,				this->dJpz_n3_rk2)
#pragma acc enter data copyin(this->dJpx_n1_rk2.matrix[0:n], this->dJpy_n1_rk2.matrix[0:n], this->dJpz_n1_rk2.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n2_rk2.matrix[0:n], this->dJpy_n2_rk2.matrix[0:n], this->dJpz_n2_rk2.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n3_rk2.matrix[0:n], this->dJpy_n3_rk2.matrix[0:n], this->dJpz_n3_rk2.matrix[0:n])

#pragma acc enter data copyin(this->dJpx_n1_rk3,			 this->dJpy_n1_rk3,				this->dJpz_n1_rk3)
#pragma acc enter data copyin(this->dJpx_n2_rk3,			 this->dJpy_n2_rk3,				this->dJpz_n2_rk3)
#pragma acc enter data copyin(this->dJpx_n3_rk3,			 this->dJpy_n3_rk3,				this->dJpz_n3_rk3)
#pragma acc enter data copyin(this->dJpx_n1_rk3.matrix[0:n], this->dJpy_n1_rk3.matrix[0:n], this->dJpz_n1_rk3.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n2_rk3.matrix[0:n], this->dJpy_n2_rk3.matrix[0:n], this->dJpz_n2_rk3.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n3_rk3.matrix[0:n], this->dJpy_n3_rk3.matrix[0:n], this->dJpz_n3_rk3.matrix[0:n])

#pragma acc enter data copyin(this->dJpx_n1_rk4,			 this->dJpy_n1_rk4,				this->dJpz_n1_rk4)
#pragma acc enter data copyin(this->dJpx_n2_rk4,			 this->dJpy_n2_rk4,				this->dJpz_n2_rk4)
#pragma acc enter data copyin(this->dJpx_n3_rk4,			 this->dJpy_n3_rk4,				this->dJpz_n3_rk4)
#pragma acc enter data copyin(this->dJpx_n1_rk4.matrix[0:n], this->dJpy_n1_rk4.matrix[0:n], this->dJpz_n1_rk4.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n2_rk4.matrix[0:n], this->dJpy_n2_rk4.matrix[0:n], this->dJpz_n2_rk4.matrix[0:n])
#pragma acc enter data copyin(this->dJpx_n3_rk4.matrix[0:n], this->dJpy_n3_rk4.matrix[0:n], this->dJpz_n3_rk4.matrix[0:n])

	if (pt_glb->if_periodic_allsurface == false) {
#pragma acc enter data copyin(this->DEx_em_t1)
#pragma acc enter data copyin(this->DEy_em_t1)
#pragma acc enter data copyin(this->DEz_em_t1)
		if (pt_glb->if_1D_ABC == false) {
#pragma acc enter data copyin(this->DEx_em_t2)
#pragma acc enter data copyin(this->DEy_em_t2)
#pragma acc enter data copyin(this->DEz_em_t2)

#pragma acc enter data copyin(this->DEx_em_t3)
#pragma acc enter data copyin(this->DEy_em_t3)
#pragma acc enter data copyin(this->DEz_em_t3)

#pragma acc enter data copyin(this->DEx_em_t4)
#pragma acc enter data copyin(this->DEy_em_t4)
#pragma acc enter data copyin(this->DEz_em_t4)
		}

#pragma acc enter data copyin(this->DEx_em_t1.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->DEy_em_t1.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->DEz_em_t1.matrix[0:(nx+1)*(ny+1)*nz])
		if (pt_glb->if_1D_ABC == false) {
#pragma acc enter data copyin(this->DEx_em_t2.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->DEy_em_t2.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->DEz_em_t2.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->DEx_em_t3.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->DEy_em_t3.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->DEz_em_t3.matrix[0:(nx+1)*(ny+1)*nz])

#pragma acc enter data copyin(this->DEx_em_t4.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc enter data copyin(this->DEy_em_t4.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc enter data copyin(this->DEz_em_t4.matrix[0:(nx+1)*(ny+1)*nz])
		}
	}
}

void EMdynamic_system::copy_from_device() {
#pragma acc update host(DHx_em_cell.matrix[0:n],DHy_em_cell.matrix[0:n],DHz_em_cell.matrix[0:n])
#pragma acc update host(DEx_em_cell.matrix[0:n],DEy_em_cell.matrix[0:n],DEz_em_cell.matrix[0:n])
}

void EMdynamic_system::copyYee_from_device() {
#pragma acc update host(DHx_em.matrix[0:(nx+1)*ny*nz])
#pragma acc update host(DHy_em.matrix[0:nx*(ny+1)*nz])
#pragma acc update host(DHz_em.matrix[0:nx*ny*(nz+1)])

#pragma acc update host(DEx_em.matrix[0:nx*(ny+1)*(nz+1)])
#pragma acc update host(DEy_em.matrix[0:(nx+1)*ny*(nz+1)])
#pragma acc update host(DEz_em.matrix[0:(nx+1)*(ny+1)*nz])
}

void EMdynamic_system::copyJp_from_device() {
#pragma acc update host(Jpx_n1.matrix[0:n],Jpy_n1.matrix[0:n],Jpz_n1.matrix[0:n])
#pragma acc update host(Jpx_n2.matrix[0:n],Jpy_n2.matrix[0:n],Jpz_n2.matrix[0:n])
#pragma acc update host(Jpx_n3.matrix[0:n],Jpy_n3.matrix[0:n],Jpz_n3.matrix[0:n])
}
