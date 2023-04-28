#include "matrix.h"

//template<typename T>
//matrix3d<T>::matrix3d() {
//
//}
//
//template<typename T>
//matrix3d<T>::matrix3d(const size_t& dimx_in, const size_t& dimy_in, const size_t& dimz_in)
//	:dimx(dimx_in), dimy(dimy_in), dimz(dimz_in) {
//	n = dimx * dimy * dimz;
//	//pt_matrix = new T * *[dimx];
//	//for (size_t i = 0; i < dimx; i++) {
//	//	pt_matrix[i] = new T * [dimy];
//	//	for (size_t j = 0; j < dimy; j++) {
//	//		pt_matrix[i][j] = new T[dimz];
//	//		for (size_t k = 0; k < dimz; k++) {
//	//			pt_matrix[i][j][k] = static_cast<T>(0);
//	//		}
//	//	}
//	//}
//
//	matrix = new T[n];
//	for (size_t i = 0; i < n; i++) {
//		matrix[i] = static_cast<T>(0);
//	}
//}

//template<typename T>
//void matrix3d<T>::initialize(const size_t& dimx_in, const size_t& dimy_in, const size_t& dimz_in) {
//	dimx = dimx_in; dimy = dimy_in; dimz = dimz_in;
//	n = dimx * dimy * dimz;
//	//pt_matrix = new T * *[dimx];
//	//for (size_t i = 0; i < dimx; i++) {
//	//	pt_matrix[i] = new T * [dimy];
//	//	for (size_t j = 0; j < dimy; j++) {
//	//		pt_matrix[i][j] = new T[dimz];
//	//		for (size_t k = 0; k < dimz; k++) {
//	//			pt_matrix[i][j][k] = static_cast<T>(0);
//	//		}
//	//	}
//	//}
//
//	matrix = new T[n];
//	for (size_t i = 0; i < n; i++) {
//		matrix[i] = static_cast<T>(0);
//	}
//}

//template<typename T>
//T& matrix3d<T>::operator()(const size_t& idx, const size_t& idy, const size_t& idz) {
//	unsigned int id = idx * (dimz * dimy) + idy * (dimz)+idz;
//	return matrix[id];
//}
//
//template<typename T>
//void matrix3d<T>::operator=(const matrix3d<T>& matrix_other) {
//	for (size_t i = 0; i < n; i++) {
//		matrix[i] = matrix_other.matrix[i];
//	}
//}
//
//template<typename T>
//void matrix3d<T>::operator+=(const matrix3d<T>& matrix_other) {
//	for (size_t i = 0; i < n; i++) {
//		matrix[i] = matrix[i] + matrix_other.matrix[i];
//	}
//}
//
//template<typename T>
//void matrix3d<T>::operator*=(const T& multiplier) {
//	for (size_t i = 0; i < n; i++) {
//		matrix[i] = matrix[i] * multiplier;
//	}
//}

// template<typename T>
// T* & matrix3d<T>::get_pointer() {
//	 return matrix;
//}

//template<typename T>
//void matrix3d<T>::nullify() {
//	//for (size_t i = 0; i < dimx; i++) {
//	//	for (size_t j = 0; j < dimy; j++) {
//	//		delete[] pt_matrix[i][j];
//	//		pt_matrix[i][j] = NULL;
//	//	}
//	//	delete[] pt_matrix[i];
//	//	pt_matrix[i] = NULL;
//	//}
//	//delete pt_matrix;
//	//pt_matrix = NULL;
//
//	delete matrix;
//	matrix = 0;
//}


//template<typename T>
//void matrix3d<T>::transformTo1D() {
//	for (size_t k = 0; k < dimz; k++)
//	{
//		for (size_t j = 0; j < dimy; j++)
//		{
//			for (size_t i = 0; i < dimx; i++)
//			{
//				pt_matrix_1D[k * (dimx * dimy) + j * dimx + i] = pt_matrix[i][j][k];
//			}
//		}
//	}
//}
//
//template<typename T>
//void matrix3d<T>::transformTo3D() {
//	for (size_t k = 0; k < dimz; k++)
//	{
//		for (size_t j = 0; j < dimy; j++)
//		{
//			for (size_t i = 0; i < dimx; i++)
//			{
//				pt_matrix[i][j][k] = pt_matrix_1D[k * (dimx * dimy) + j * dimx + i];
//			}
//		}
//	}
//}


//template<typename T>
//void matrix3d<T>::multiply(const T& multiplier) {
//	for (size_t i = 0; i < dimx; i++) {
//		for (size_t j = 0; j < dimy; j++) {
//			for (size_t k = 0; k < dimz; k++) {
//				pt_matrix[i][j][k] = pt_matrix[i][j][k] * multiplier;
//			}
//		}
//	}
//}

//template<typename T>
//matrix2d<T>::matrix2d() {
//
//}
//
//template <typename T>
//matrix2d<T>::matrix2d(const size_t& dimx_in, const size_t& dimy_in)
//	:dimx(dimx_in), dimy(dimy_in) {
//	n = dimx * dimy;
//	matrix = new T[n];
//	for (size_t i = 0; i < n; i++) {
//		matrix[i] = static_cast<T>(0);
//	}
//	//pt_matrix = new T * [dimx];
//	//for (size_t i = 0; i < dimx; i++) {
//	//	pt_matrix[i] = new T[dimy];
//	//	for (size_t j = 0; j < dimy; j++) {
//	//		pt_matrix[i][j] = static_cast<T>(0);
//	//	}
//	//}
//}

//template <typename T>
//void matrix2d<T>::initialize(const size_t& dimx_in, const size_t& dimy_in) {
//	dimx = dimx_in; dimy = dimy_in;
//	n = dimx * dimy;
//	matrix = new T[n];
//	for (size_t i = 0; i < n; i++) {
//		matrix[i] = static_cast<T>(0);
//	}
//
//	//pt_matrix = new T * [dimx];
//	//for (size_t i = 0; i < dimx; i++) {
//	//	pt_matrix[i] = new T[dimy];
//	//	for (size_t j = 0; j < dimy; j++) {
//	//		pt_matrix[i][j] = static_cast<T>(0);
//	//	}
//	//}
//}
//
//template<typename T>
//T& matrix2d<T>::operator()(const size_t& idx, const size_t& idy) {
//	size_t id = idx * dimy + idy;
//	return matrix[id];
//}
//
//
//template<typename T>
//void matrix2d<T>::nullify() {
//	//for (size_t i = 0; i < dimx; i++) {
//	//	delete[] pt_matrix[i];
//	//	pt_matrix[i] = NULL;
//	//}
//	delete matrix;
//	matrix = 0;
//}
