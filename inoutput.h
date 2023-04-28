#pragma once
#include "global.h"
#include "magnetic_system.h"
#include "ferroelectric_system.h"
#include "EMdynamic_system.h"
#include "elastic_system.h"

class inoutput {
private:
	long int nx, ny, nz;
	class global_parameters* pt_glb;
public:
	void input_m(magnetic_system*);
	void input_AFMm(magnetic_system*);
	void input_Hstat(magnetic_system*);
	void input_pandq(ferroelectric_system*);
	void input_Estat(ferroelectric_system*);
	void input_em_Yee(EMdynamic_system*);
	void input_uandv(elastic_system*);
	void input_elastoforce(elastic_system*);
	void input_straint0(elastic_system*);
	void input_eigenstraint0(elastic_system*);
	void input_Jp(EMdynamic_system*);

	void output_m(unsigned long long int&, magnetic_system*);
	void output_magcell(unsigned long long int&, magnetic_system*);
	void output_AFMm(unsigned long long int&, magnetic_system*);

	void output_p(unsigned long long int&, ferroelectric_system*);
	void output_q(unsigned long long int&, ferroelectric_system*);

	void output_Estat(unsigned long long int&, ferroelectric_system*);
	void output_Eem(unsigned long long int&, EMdynamic_system*);
	void output_EemYee(unsigned long long int&, EMdynamic_system*);

	void output_Hstat(unsigned long long int&, magnetic_system*);
	void output_Hem(unsigned long long int&, EMdynamic_system*);
	void output_HemYee(unsigned long long int&, EMdynamic_system*);

	void output_em_onecell(unsigned long long int&, EMdynamic_system*);

	void output_strain(unsigned long long int&, elastic_system*);
	void output_straint0(unsigned long long int&, elastic_system*);
	void output_uandv(unsigned long long int&, elastic_system*);
	void output_elastoforce(unsigned long long int&, elastic_system*);
	void output_eigenstraint0_crt(unsigned long long int&, elastic_system*);

	void output_Jp(unsigned long long int&, EMdynamic_system*);
	void output_Jishe(unsigned long long int&, magnetic_system*);

	void output_averagem(unsigned long long int&, magnetic_system*);
	void output_averageAFMm(unsigned long long int&, magnetic_system*);
	void output_averagep(unsigned long long int&, ferroelectric_system*);

public:
	void get_dimension(long int&, long int&, long int&, global_parameters*);
	static inoutput io;
};
