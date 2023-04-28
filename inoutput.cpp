#include "inoutput.h"
#include <string>
#include <fstream>
#include<stdio.h>
#include"mathlib.h"

//---------------------------------//
//				INPUT				//
//---------------------------------//

void inoutput::input_m(magnetic_system* mag) {
	long x = 0;
	long y = 0;
	long z = 0;
	double inx, iny, inz;
	std::ifstream file("m.in");
	if (file.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			file >> x >> y >> z >> inx >> iny >> inz;
			mag->mx_glb(x - 1, y - 1, z - 1) = inx;
			mag->my_glb(x - 1, y - 1, z - 1) = iny;
			mag->mz_glb(x - 1, y - 1, z - 1) = inz;
		}
	}
	file.close();
}

void inoutput::input_AFMm(magnetic_system* mag) {
	long x = 0;
	long y = 0;
	long z = 0;
	double inx1, iny1, inz1;
	double inx2, iny2, inz2;
	std::ifstream file("m_AFM.in");
	if (file.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			file >> x >> y >> z >> inx1 >> iny1 >> inz1 >> inx2 >> iny2 >> inz2;
			mag->mx_AFM1_glb(x - 1, y - 1, z - 1) = inx1;
			mag->my_AFM1_glb(x - 1, y - 1, z - 1) = iny1;
			mag->mz_AFM1_glb(x - 1, y - 1, z - 1) = inz1;

			mag->mx_AFM2_glb(x - 1, y - 1, z - 1) = inx2;
			mag->my_AFM2_glb(x - 1, y - 1, z - 1) = iny2;
			mag->mz_AFM2_glb(x - 1, y - 1, z - 1) = inz2;
		}
	}
	file.close();
}

void inoutput::input_Hstat(magnetic_system* mag) {
	long x = 0;
	long y = 0;
	long z = 0;
	double inx, iny, inz;

	std::ifstream file("Hstat.in");
	if (file.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			file >> x >> y >> z >> inx >> iny >> inz;
			mag->Hx_stat(x - 1, y - 1, z - 1) = inx;
			mag->Hy_stat(x - 1, y - 1, z - 1) = iny;
			mag->Hz_stat(x - 1, y - 1, z - 1) = inz;
		}
	}
	file.close();
}

void inoutput::input_pandq(ferroelectric_system* fe) {
	long x = 0;
	long y = 0;
	long z = 0;
	double inx, iny, inz;

	std::ifstream filep("p.in");
	if (filep.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			filep >> x >> y >> z >> inx >> iny >> inz;
			fe->px_glb(x - 1, y - 1, z - 1) = inx;
			fe->py_glb(x - 1, y - 1, z - 1) = iny;
			fe->pz_glb(x - 1, y - 1, z - 1) = inz;
		}
	}
	filep.close();

	std::ifstream fileq("q.in");
	if (fileq.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			fileq >> x >> y >> z >> inx >> iny >> inz;
			fe->qx_glb(x - 1, y - 1, z - 1) = inx;
			fe->qy_glb(x - 1, y - 1, z - 1) = iny;
			fe->qz_glb(x - 1, y - 1, z - 1) = inz;
		}
	}
	fileq.close();
}

void inoutput::input_Estat(ferroelectric_system* fe) {
	long x = 0;
	long y = 0;
	long z = 0;
	double inx, iny, inz;

	std::ifstream file("Estat.in");
	if (file.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			file >> x >> y >> z >> inx >> iny >> inz;
			fe->Ex_stat(x - 1, y - 1, z - 1) = inx;
			fe->Ey_stat(x - 1, y - 1, z - 1) = iny;
			fe->Ez_stat(x - 1, y - 1, z - 1) = inz;
		}
	}
	file.close();
}


void inoutput::input_em_Yee(EMdynamic_system* em) {
	long x = 0;
	long y = 0;
	long z = 0;
	double inn;

	std::ifstream fileEx("Ex_Yee.in");
	if (fileEx.is_open()) {
		for (unsigned long i = 0; i < nx * (ny + 1) * (nz + 1); i++)
		{
			fileEx >> x >> y >> z >> inn;
			em->DEx_em(x - 1, y - 1, z - 1) = inn;
		}
	}
	fileEx.close();

	std::ifstream fileEy("Ey_Yee.in");
	if (fileEy.is_open()) {
		for (unsigned long i = 0; i < (nx + 1) * ny * (nz + 1); i++)
		{
			fileEy >> x >> y >> z >> inn;
			em->DEy_em(x - 1, y - 1, z - 1) = inn;
		}
	}
	fileEy.close();

	std::ifstream fileEz("Ez_Yee.in");
	if (fileEz.is_open()) {
		for (unsigned long i = 0; i < (nx + 1) * (ny + 1) * nz; i++)
		{
			fileEz >> x >> y >> z >> inn;
			em->DEz_em(x - 1, y - 1, z - 1) = inn;
		}
	}
	fileEz.close();

	std::ifstream fileHx("Hx_Yee.in");
	if (fileHx.is_open()) {
		for (unsigned long i = 0; i < (nx + 1) * ny * nz; i++)
		{
			fileHx >> x >> y >> z >> inn;
			em->DHx_em(x - 1, y - 1, z - 1) = inn;
		}
	}
	fileHx.close();

	std::ifstream fileHy("Hy_Yee.in");
	if (fileHy.is_open()) {
		for (unsigned long i = 0; i < nx * (ny + 1) * nz; i++)
		{
			fileHy >> x >> y >> z >> inn;
			em->DHy_em(x - 1, y - 1, z - 1) = inn;
		}
	}
	fileHy.close();

	std::ifstream fileHz("Hz_Yee.in");
	if (fileHz.is_open()) {
		for (unsigned long i = 0; i < nx * ny * (nz + 1); i++)
		{
			fileHz >> x >> y >> z >> inn;
			em->DHz_em(x - 1, y - 1, z - 1) = inn;
		}
	}
	fileHz.close();
}

void inoutput::input_uandv(elastic_system* elasto) {
	long x = 0;
	long y = 0;
	long z = 0;
	double inx, iny, inz;

	std::ifstream fileu("u.in");
	if (fileu.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			fileu >> x >> y >> z >> inx >> iny >> inz;
			elasto->Dux_glb(x - 1, y - 1, z - 1) = inx;
			elasto->Duy_glb(x - 1, y - 1, z - 1) = iny;
			elasto->Duz_glb(x - 1, y - 1, z - 1) = inz;
		}
	}
	fileu.close();

	std::ifstream filev("v.in");
	if (filev.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			filev >> x >> y >> z >> inx >> iny >> inz;
			elasto->vx_glb(x - 1, y - 1, z - 1) = inx;
			elasto->vy_glb(x - 1, y - 1, z - 1) = iny;
			elasto->vz_glb(x - 1, y - 1, z - 1) = inz;
		}
	}
	filev.close();
}

void inoutput::input_elastoforce(elastic_system* elasto) {
	long x = 0;
	long y = 0;
	long z = 0;
	double inx, iny, inz;

	std::ifstream filef("elastoforce.in");
	if (filef.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			filef >> x >> y >> z >> inx >> iny >> inz;
			elasto->force_x_store(x - 1, y - 1, z - 1) = inx;
			elasto->force_y_store(x - 1, y - 1, z - 1) = iny;
			elasto->force_z_store(x - 1, y - 1, z - 1) = inz;
		}
	}
	filef.close();
}

void inoutput::input_straint0(elastic_system* elasto) {
	long x = 0;
	long y = 0;
	long z = 0;
	double inxx, inyy, inzz, inyz, inxz, inxy;

	std::ifstream file("strain_stat.in");
	if (file.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			file >> x >> y >> z >> inxx >> inyy >> inzz >> inyz >> inxz >> inxy;
			elasto->exxt0_glb(x - 1, y - 1, z - 1) = inxx;
			elasto->eyyt0_glb(x - 1, y - 1, z - 1) = inyy;
			elasto->ezzt0_glb(x - 1, y - 1, z - 1) = inzz;
			elasto->eyzt0_glb(x - 1, y - 1, z - 1) = inyz;
			elasto->exzt0_glb(x - 1, y - 1, z - 1) = inxz;
			elasto->exyt0_glb(x - 1, y - 1, z - 1) = inxy;
		}
	}
	file.close();
}

void inoutput::input_eigenstraint0(elastic_system* elasto) {
	long x = 0;
	long y = 0;
	long z = 0;
	double inxx, inyy, inzz, inyz, inxz, inxy;

	std::ifstream file("eigenstraint0_crt.in");
	if (file.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			file >> x >> y >> z >> inxx >> inyy >> inzz >> inyz >> inxz >> inxy;
			elasto->exx0t0_crt(x - 1, y - 1, z - 1) = inxx;
			elasto->eyy0t0_crt(x - 1, y - 1, z - 1) = inyy;
			elasto->ezz0t0_crt(x - 1, y - 1, z - 1) = inzz;
			elasto->eyz0t0_crt(x - 1, y - 1, z - 1) = inyz;
			elasto->exz0t0_crt(x - 1, y - 1, z - 1) = inxz;
			elasto->exy0t0_crt(x - 1, y - 1, z - 1) = inxy;
		}
	}
	file.close();
}

void inoutput::input_Jp(EMdynamic_system* em) {
	long x = 0;
	long y = 0;
	long z = 0;
	double inx, iny, inz;

	std::ifstream file1("Jp1.in");
	if (file1.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			file1 >> x >> y >> z >> inx >> iny >> inz;
			em->Jpx_n1(x - 1, y - 1, z - 1) = inx;
			em->Jpy_n1(x - 1, y - 1, z - 1) = iny;
			em->Jpz_n1(x - 1, y - 1, z - 1) = inz;
		}
	}
	file1.close();

	std::ifstream file2("Jp2.in");
	if (file2.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			file2 >> x >> y >> z >> inx >> iny >> inz;
			em->Jpx_n2(x - 1, y - 1, z - 1) = inx;
			em->Jpy_n2(x - 1, y - 1, z - 1) = iny;
			em->Jpz_n2(x - 1, y - 1, z - 1) = inz;
		}
	}
	file2.close();

	std::ifstream file3("Jp3.in");
	if (file3.is_open()) {
		for (unsigned long i = 0; i < nx * ny * nz; i++)
		{
			file3 >> x >> y >> z >> inx >> iny >> inz;
			em->Jpx_n3(x - 1, y - 1, z - 1) = inx;
			em->Jpy_n3(x - 1, y - 1, z - 1) = iny;
			em->Jpz_n3(x - 1, y - 1, z - 1) = inz;
		}
	}
	file3.close();
}

//---------------------------------//
//				OUTPUT			   //
//---------------------------------//

void inoutput::output_m(unsigned long long int& nstep, magnetic_system* pt_mag) {
	std::string filename = "m." + std::to_string(nstep) + ".dat";
	//std::ofstream file(filename.c_str());
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_file;
	pt_file = fopen(filename.c_str(), "w");

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				fprintf(pt_file, "%ld  %ld  %ld  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, pt_mag->mx_glb(i, j, k), pt_mag->my_glb(i, j, k), pt_mag->mz_glb(i, j, k));
			}
		}
	}

	fclose(pt_file);
}

void inoutput::output_magcell(unsigned long long int& nstep, magnetic_system* pt_mag) {
	unsigned long long index, i, j, k;

	std::string filename = "m_magcell." + std::to_string(nstep) + ".dat";
	FILE* pt_file;
	pt_file = fopen(filename.c_str(), "w");

	for (unsigned long id = 0; id < pt_glb->NFM; id++) {
		index = pt_glb->magcell_index[id];
		i = index / (ny * nz);
		j = (index - i * (ny * nz)) / nz;
		k = index - i * (ny * nz) - j * nz;
		fprintf(pt_file, "%ld  %ld  %ld  %.7e  %.7e  %.7e\n", \
			i + 1, j + 1, k + 1, pt_mag->mx_glb(index), pt_mag->my_glb(index), pt_mag->mz_glb(index));
	}

	fclose(pt_file);
}

void inoutput::output_AFMm(unsigned long long int& nstep, magnetic_system* pt_mag) {
	std::string filename = "m_AFM." + std::to_string(nstep) + ".dat";
	//std::ofstream file(filename.c_str());
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_file;
	pt_file = fopen(filename.c_str(), "w");

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				fprintf(pt_file, "%ld  %ld  %ld  %.7e  %.7e  %.7e  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, \
					pt_mag->mx_AFM1_glb(i, j, k), pt_mag->my_AFM1_glb(i, j, k), pt_mag->mz_AFM1_glb(i, j, k), \
					pt_mag->mx_AFM2_glb(i, j, k), pt_mag->my_AFM2_glb(i, j, k), pt_mag->mz_AFM2_glb(i, j, k));
			}
		}
	}

	fclose(pt_file);
}

void inoutput::output_p(unsigned long long int& nstep, ferroelectric_system* pt_fe) {
	std::string filename = "p." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_file;
	pt_file = fopen(filename.c_str(), "w");

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				fprintf(pt_file, "%ld  %ld  %ld  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, pt_fe->px_glb(i, j, k), pt_fe->py_glb(i, j, k), pt_fe->pz_glb(i, j, k));
			}
		}
	}

	fclose(pt_file);
}

void inoutput::output_q(unsigned long long int& nstep, ferroelectric_system* pt_fe) {
	std::string filename = "q." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_file;
	pt_file = fopen(filename.c_str(), "w");

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				fprintf(pt_file, "%ld  %ld  %ld  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, pt_fe->qx_glb(i, j, k), pt_fe->qy_glb(i, j, k), pt_fe->qz_glb(i, j, k));
			}
		}
	}

	fclose(pt_file);
}

void inoutput::output_Estat(unsigned long long int& nstep, ferroelectric_system* pt_fe) {
	std::string filename = "Estat." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_file;
	pt_file = fopen(filename.c_str(), "w");

	double Ex, Ey, Ez;

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				Ex = pt_fe->Ex_stat(i, j, k);
				Ey = pt_fe->Ey_stat(i, j, k);
				Ez = pt_fe->Ez_stat(i, j, k);
				fprintf(pt_file, "%ld  %ld  %ld  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, Ex, Ey, Ez);
			}
		}
	}

	fclose(pt_file);
}

void inoutput::output_Eem(unsigned long long int& nstep, EMdynamic_system* pt_em) {
	std::string filename = "Eem." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_file;
	pt_file = fopen(filename.c_str(), "w");

	double Ex, Ey, Ez;

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				Ex = pt_em->DEx_em_cell(i, j, k);
				Ey = pt_em->DEy_em_cell(i, j, k);
				Ez = pt_em->DEz_em_cell(i, j, k);
				fprintf(pt_file, "%ld  %ld  %ld  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, Ex, Ey, Ez);
			}
		}
	}

	fclose(pt_file);
}

void inoutput::output_EemYee(unsigned long long int& nstep, EMdynamic_system* pt_em) {
	double Ex, Ey, Ez;

	std::string filenamex = "Ex_Yee." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_filex;
	pt_filex = fopen(filenamex.c_str(), "w");

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny + 1; j++) {
			for (long int k = 0; k < nz + 1; k++) {
				Ex = pt_em->DEx_em(i, j, k);

				fprintf(pt_filex, "%ld  %ld  %ld  %.7e\n", \
					i + 1, j + 1, k + 1, Ex);
			}
		}
	}
	fclose(pt_filex);

	std::string filenamey = "Ey_Yee." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_filey;
	pt_filey = fopen(filenamey.c_str(), "w");

	for (long int i = 0; i < nx + 1; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz + 1; k++) {
				Ey = pt_em->DEy_em(i, j, k);

				fprintf(pt_filey, "%ld  %ld  %ld  %.7e\n", \
					i + 1, j + 1, k + 1, Ey);
			}
		}
	}
	fclose(pt_filey);

	std::string filenamez = "Ez_Yee." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_filez;
	pt_filez = fopen(filenamez.c_str(), "w");

	for (long int i = 0; i < nx + 1; i++) {
		for (long int j = 0; j < ny + 1; j++) {
			for (long int k = 0; k < nz; k++) {
				Ez = pt_em->DEz_em(i, j, k);

				fprintf(pt_filez, "%ld  %ld  %ld  %.7e\n", \
					i + 1, j + 1, k + 1, Ez);
			}
		}
	}
	fclose(pt_filez);
}

void inoutput::output_Hstat(unsigned long long int& nstep, magnetic_system* pt_mag) {
	std::string filename = "Hstat." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_file;
	pt_file = fopen(filename.c_str(), "w");

	double Hx, Hy, Hz;

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				Hx = pt_mag->Hx_stat(i, j, k);
				Hy = pt_mag->Hy_stat(i, j, k);
				Hz = pt_mag->Hz_stat(i, j, k);
				fprintf(pt_file, "%ld  %ld  %ld  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, Hx, Hy, Hz);
			}
		}
	}

	fclose(pt_file);
}

void inoutput::output_Hem(unsigned long long int& nstep, EMdynamic_system* pt_em) {
	std::string filename = "Hem." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_file;
	pt_file = fopen(filename.c_str(), "w");

	double Hx, Hy, Hz;

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				Hx = pt_em->DHx_em_cell(i, j, k);
				Hy = pt_em->DHy_em_cell(i, j, k);
				Hz = pt_em->DHz_em_cell(i, j, k);
				fprintf(pt_file, "%ld  %ld  %ld  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, Hx, Hy, Hz);
			}
		}
	}

	fclose(pt_file);
}

void inoutput::output_HemYee(unsigned long long int& nstep, EMdynamic_system* pt_em) {
	double Hx, Hy, Hz;

	std::string filenamex = "Hx_Yee." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_filex;
	pt_filex = fopen(filenamex.c_str(), "w");

	for (long int i = 0; i < nx + 1; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				Hx = pt_em->DHx_em(i, j, k);

				fprintf(pt_filex, "%ld  %ld  %ld  %.7e\n", \
					i + 1, j + 1, k + 1, Hx);
			}
		}
	}
	fclose(pt_filex);

	std::string filenamey = "Hy_Yee." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_filey;
	pt_filey = fopen(filenamey.c_str(), "w");

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny + 1; j++) {
			for (long int k = 0; k < nz; k++) {
				Hy = pt_em->DHy_em(i, j, k);

				fprintf(pt_filey, "%ld  %ld  %ld  %.7e\n", \
					i + 1, j + 1, k + 1, Hy);
			}
		}
	}
	fclose(pt_filey);

	std::string filenamez = "Hz_Yee." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_filez;
	pt_filez = fopen(filenamez.c_str(), "w");

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz + 1; k++) {
				Hz = pt_em->DHz_em(i, j, k);

				fprintf(pt_filez, "%ld  %ld  %ld  %.7e\n", \
					i + 1, j + 1, k + 1, Hz);
			}
		}
	}
	fclose(pt_filez);
}

void inoutput::output_em_onecell(unsigned long long int& nstep, EMdynamic_system* pt_em) {
	FILE* pt_file = fopen("EM_onecell.dat", "a");

	double Hx, Hy, Hz;
	double Ex, Ey, Ez;
	unsigned long long index = pt_glb->em_onecell_index - 1;

	Hx = pt_em->DHx_em_cell(index);
	Hy = pt_em->DHy_em_cell(index);
	Hz = pt_em->DHz_em_cell(index);
	
	Ex = pt_em->DEx_em_cell(index);
	Ey = pt_em->DEy_em_cell(index);
	Ez = pt_em->DEz_em_cell(index);
	
	fprintf(pt_file, "%lld %.7e  %.7e  %.7e %.7e  %.7e  %.7e\n", \
		nstep, Ex, Ey, Ez, Hx, Hy, Hz);

	fclose(pt_file);
}

void inoutput::output_strain(unsigned long long int& nstep, elastic_system* pt_elasto) {
	std::string filename = "strain." + std::to_string(nstep) + ".dat";

	FILE* pt_file;
	pt_file = fopen(filename.c_str(), "w");

	double exx, eyy, ezz, eyz, exz, exy;

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {

				exx = pt_elasto->exxt0_glb(i, j, k) + pt_elasto->Dexx_glb(i, j, k);
				eyy = pt_elasto->eyyt0_glb(i, j, k) + pt_elasto->Deyy_glb(i, j, k);
				ezz = pt_elasto->ezzt0_glb(i, j, k) + pt_elasto->Dezz_glb(i, j, k);
				eyz = pt_elasto->eyzt0_glb(i, j, k) + pt_elasto->Deyz_glb(i, j, k);
				exz = pt_elasto->exzt0_glb(i, j, k) + pt_elasto->Dexz_glb(i, j, k);
				exy = pt_elasto->exyt0_glb(i, j, k) + pt_elasto->Dexy_glb(i, j, k);

				fprintf(pt_file, "%ld  %ld  %ld  %.7e  %.7e  %.7e  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, exx, eyy, ezz, eyz, exz, exy);
			}
		}
	}

	fclose(pt_file);
}

void inoutput::output_straint0(unsigned long long int& nstep, elastic_system* pt_elasto) {
	std::string filename = "strain_stat." + std::to_string(nstep) + ".dat";

	FILE* pt_file;
	pt_file = fopen(filename.c_str(), "w");

	double exx, eyy, ezz, eyz, exz, exy;

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {

				exx = pt_elasto->exxt0_glb(i, j, k);
				eyy = pt_elasto->eyyt0_glb(i, j, k);
				ezz = pt_elasto->ezzt0_glb(i, j, k);
				eyz = pt_elasto->eyzt0_glb(i, j, k);
				exz = pt_elasto->exzt0_glb(i, j, k);
				exy = pt_elasto->exyt0_glb(i, j, k);

				fprintf(pt_file, "%ld  %ld  %ld  %.7e  %.7e  %.7e  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, exx, eyy, ezz, eyz, exz, exy);
			}
		}
	}

	fclose(pt_file);
}

void inoutput::output_uandv(unsigned long long int& nstep, elastic_system* pt_elasto) {
	double ux, uy, uz, vx, vy, vz;

	std::string filenameu = "u." + std::to_string(nstep) + ".dat";

	FILE* pt_fileu;
	pt_fileu = fopen(filenameu.c_str(), "w");
	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				ux = pt_elasto->Dux_glb(i, j, k);
				uy = pt_elasto->Duy_glb(i, j, k);
				uz = pt_elasto->Duz_glb(i, j, k);
				fprintf(pt_fileu, "%ld  %ld  %ld  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, ux, uy, uz);
			}
		}
	}
	fclose(pt_fileu);

	std::string filenamev = "v." + std::to_string(nstep) + ".dat";

	FILE* pt_filev;
	pt_filev = fopen(filenamev.c_str(), "w");
	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				vx = pt_elasto->vx_glb(i, j, k);
				vy = pt_elasto->vy_glb(i, j, k);
				vz = pt_elasto->vz_glb(i, j, k);
				fprintf(pt_filev, "%ld  %ld  %ld  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, vx, vy, vz);
			}
		}
	}
	fclose(pt_filev);
}

void inoutput::output_elastoforce(unsigned long long int& nstep, elastic_system* pt_elasto) {
	double fx, fy, fz;

	std::string filenamef = "elastoforce." + std::to_string(nstep) + ".dat";

	FILE* pt_filef;
	pt_filef = fopen(filenamef.c_str(), "w");
	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				fx = pt_elasto->force_x_store(i, j, k);
				fy = pt_elasto->force_y_store(i, j, k);
				fz = pt_elasto->force_z_store(i, j, k);
				fprintf(pt_filef, "%ld  %ld  %ld  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, fx, fy, fz);
			}
		}
	}
	fclose(pt_filef);
}

void inoutput::output_eigenstraint0_crt(unsigned long long int& nstep, elastic_system* pt_elasto) {
	std::string filename = "eigenstraint0_crt." + std::to_string(nstep) + ".dat";

	FILE* pt_file;
	pt_file = fopen(filename.c_str(), "w");

	double exx0, eyy0, ezz0, eyz0, exz0, exy0;

	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {

				exx0 = pt_elasto->exx0t0_crt(i, j, k);
				eyy0 = pt_elasto->eyy0t0_crt(i, j, k);
				ezz0 = pt_elasto->ezz0t0_crt(i, j, k);
				eyz0 = pt_elasto->eyz0t0_crt(i, j, k);
				exz0 = pt_elasto->exz0t0_crt(i, j, k);
				exy0 = pt_elasto->exy0t0_crt(i, j, k);

				fprintf(pt_file, "%ld  %ld  %ld  %.7e  %.7e  %.7e  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, exx0, eyy0, ezz0, eyz0, exz0, exy0);
			}
		}
	}

	fclose(pt_file);
}

void inoutput::output_Jp(unsigned long long int& nstep, EMdynamic_system* pt_em) {
	double Jx, Jy, Jz;

	std::string filename1 = "Jp1." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_file1;
	pt_file1 = fopen(filename1.c_str(), "w");
	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				Jx = pt_em->Jpx_n1(i, j, k);
				Jy = pt_em->Jpy_n1(i, j, k);
				Jz = pt_em->Jpz_n1(i, j, k);
				fprintf(pt_file1, "%ld  %ld  %ld  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, Jx, Jy, Jz);
			}
		}
	}
	fclose(pt_file1);

	std::string filename2 = "Jp2." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_file2;
	pt_file2 = fopen(filename2.c_str(), "w");
	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				Jx = pt_em->Jpx_n2(i, j, k);
				Jy = pt_em->Jpy_n2(i, j, k);
				Jz = pt_em->Jpz_n2(i, j, k);
				fprintf(pt_file2, "%ld  %ld  %ld  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, Jx, Jy, Jz);
			}
		}
	}
	fclose(pt_file2);

	std::string filename3 = "Jp3." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_file3;
	pt_file3 = fopen(filename3.c_str(), "w");
	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				Jx = pt_em->Jpx_n3(i, j, k);
				Jy = pt_em->Jpy_n3(i, j, k);
				Jz = pt_em->Jpz_n3(i, j, k);
				fprintf(pt_file3, "%ld  %ld  %ld  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, Jx, Jy, Jz);
			}
		}
	}
	fclose(pt_file3);
}

void inoutput::output_Jishe(unsigned long long int& nstep, magnetic_system* pt_mag) {
	double Jx, Jy;

	std::string filename = "J_ISHE." + std::to_string(nstep) + ".dat";
	//FILE* pt_file = fopen(filename.c_str(), "w");
	FILE* pt_file;
	pt_file = fopen(filename.c_str(), "w");
	for (long int i = 0; i < nx; i++) {
		for (long int j = 0; j < ny; j++) {
			for (long int k = 0; k < nz; k++) {
				Jx = pt_mag->Jx_ISHE(i, j, k);
				Jy = pt_mag->Jy_ISHE(i, j, k);
				fprintf(pt_file, "%ld  %ld  %ld  %.7e  %.7e  %.7e\n", \
					i + 1, j + 1, k + 1, Jx, Jy, 0.);
			}
		}
	}
	fclose(pt_file);
}

void inoutput::output_averagem(unsigned long long int& nstep, magnetic_system* pt_mag) {
	FILE* pt_file = fopen("average_m.dat", "a");
	fprintf(pt_file, "%lld %.7e  %.7e  %.7e\n", \
		nstep, pt_mag->mx_ave, pt_mag->my_ave, pt_mag->mz_ave);
	fclose(pt_file);
}

void inoutput::output_averageAFMm(unsigned long long int& nstep, magnetic_system* pt_mag) {
	FILE* pt_file = fopen("average_AFMm.dat", "a");
	fprintf(pt_file, "%lld %.7e  %.7e  %.7e  %.7e  %.7e  %.7e\n", \
		nstep, \
		pt_mag->mx_AFM1_ave, pt_mag->my_AFM1_ave, pt_mag->mz_AFM1_ave,\
		pt_mag->mx_AFM2_ave, pt_mag->my_AFM2_ave, pt_mag->mz_AFM2_ave);
	fclose(pt_file);
}

void inoutput::output_averagep(unsigned long long int& nstep, ferroelectric_system* pt_fe) {
	FILE* pt_file = fopen("average_p.dat", "a");
	fprintf(pt_file, "%lld %.7e  %.7e  %.7e\n", \
		nstep, pt_fe->px_ave, pt_fe->py_ave, pt_fe->pz_ave);
	fclose(pt_file);
}

void inoutput::get_dimension(long int& nx_in, long int& ny_in, long int& nz_in, global_parameters* in_glb) {
	nx = nx_in;
	ny = ny_in;
	nz = nz_in;
	pt_glb = in_glb;
}

inoutput inoutput::io;
