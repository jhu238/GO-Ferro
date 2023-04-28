#pragma once
//#include "mathlib.h"
//#include "EMdynamic_system.h"
//class EMdynamic_system;
#include "openacc.h"

template <typename T>
class matrix3d {
public:
	long int dimx, dimy, dimz;
	long int dim;
	friend class mathlib;
	T* matrix;
	//public:
public:
	//matrix3d();
	//matrix3d(const long int& dimx_in, const long int& dimy_in, const long int& dimz_in);
	void initialize(const long int& dimx_in, const long int& dimy_in, const long int& dimz_in) {
		dimx = dimx_in; dimy = dimy_in; dimz = dimz_in;
		dim = dimx * dimy * dimz;
		matrix = new T[dim];
		for (long int i = 0; i < dim; i++) {
			matrix[i] = static_cast<T>(0);
		}
	}

	
	void initialize_device(const long int& dimx_in, const long int& dimy_in, const long int& dimz_in) {
//#pragma acc loop seq
//		for (long int i = 0; i < 1; i++) {
		dimx = dimx_in; dimy = dimy_in; dimz = dimz_in;
		dim = dimx * dimy * dimz;
		matrix = (T*)(acc_malloc(dim * sizeof(T)));
		//}
		/*matrix = new T[n];*/
		//#pragma acc loop independent
		//for (long int i = 0; i < n; i++) {
		//	matrix[i] = static_cast<T>(0);
		//}
	}


	//	void transformTo1D();
	//	void transformTo3D();

	#pragma acc routine seq
	T& operator()(const long int& idx, const long int& idy, const long int& idz) {
		//unsigned int id;
		///*#pragma acc declare device_resident(id)*/
		//id = idx * (dimz * dimy) + idy * (dimz)+idz;
		return matrix[idx * (dimz * dimy) + idy * (dimz)+idz];
	}

#pragma acc routine seq// nohost
	T& operator()(const long int& id) {
		//unsigned int id;
		///*#pragma acc declare device_resident(id)*/
		//id = idx * (dimz * dimy) + idy * (dimz)+idz;
		return matrix[id];
	}

	#pragma acc routine gang //nohost
	void operator+=(const matrix3d<T>& matrix_other) {
		#pragma acc loop gang vector
		for (long int i = 0; i < dim; i++) {
			matrix[i] = matrix[i] + matrix_other.matrix[i];
		}
	}

	#pragma acc routine gang //nohost
	void operator*=(const T& multiplier) {
		#pragma acc loop gang vector
		for (long int i = 0; i < dim; i++) {
			matrix[i] = matrix[i] * multiplier;
		}
	}

	#pragma acc routine gang //nohost
	void operator=(const matrix3d<T>& matrix_other) {
		#pragma acc loop gang vector	
		for (long int i = 0; i < dim; i++) {
			matrix[i] = matrix_other.matrix[i];
		}
	}

	#pragma acc routine seq
	T*& get_pointer() {
		return matrix;
	}
	//	void multiply(const T& multiplier);
	//friend void EMdynamic_system::transfer_pointer();

	void nullify() {
		delete matrix;
		matrix = 0;
	}

	void nullify_device() {
//#pragma acc loop seq
//		for (long int i = 0; i < 1; i++) {
			acc_free(matrix);
		//}
	}
};

//template <typename T>
//class matrix2d {
//private:
//	long int dimx, dimy;
//	size_t n;
//	friend class mathlib;
//	T* matrix = 0;
//public:
//	//matrix2d();
//	//matrix2d(const size_t& dimx_in, const size_t& dimy_in);
//	void initialize(const size_t& dimx_in, const size_t& dimy_in);
//	T& operator()(const size_t&, const size_t&);
//	void nullify();
//};
