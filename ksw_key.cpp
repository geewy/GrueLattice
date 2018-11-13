/**
  * @file ksw_key.cpp
  * @date may 2018
  * @author geewy
  * @brief Definition of Key-Switching key class
  *
  * See ksw_key.hpp for more details
  */

#include "ksw_key.hpp"

template<class Ring, uint64_t dim, uint64_t dimP, uint64_t Base, uint64_t Ksize>
KSWKey<Ring, dim, dimP, Base, Ksize>::KSWKey() : init_ok(false){}

template<class Ring, uint64_t dim, uint64_t dimP, uint64_t Base, uint64_t Ksize>
KSWKey<Ring, dim, dimP, Base, Ksize>::KSWKey(const Ring secret[dim], const Ring secret_p[dimP], const double variance){
	for(size_t i = 0; i < dim; ++i){
		int64_t pow = 1;
		for(size_t j = 0; j < Ksize; ++j){
			SECRET[i][j] = RLWE<Ring, dimP>();
			const Ring tmp = secret[i] * pow;
			SECRET[i][j].encrypt(secret_p, tmp, 1, 1, variance);
			SECRET[i][j].decomp_fft();
			pow *= Base;
		}
	}
	init_ok = true;
}

template<class Ring1, uint64_t dim1, uint64_t dimP, uint64_t Base, uint64_t Ksize>
void KSWKey<Ring1, dim1, dimP, Base, Ksize>::init(const Ring1 secret[dim1], const Ring1 secret_p[dimP], const double variance){
	for(size_t i = 0; i < dim1; ++i){
		int64_t pow = 1;
		for(size_t j = 0; j < Ksize; ++j){
			SECRET[i][j] = RLWE<Ring1, dimP>();
			const Ring1 tmp = secret[i] * pow;
			SECRET[i][j].encrypt(secret_p, tmp, 1, 1, variance);
			SECRET[i][j].decomp_fft();
			pow *= Base;
		}
	}
	init_ok = true;
}

template class KSWKey<Rp1>;
template class KSWKey<Rp2>;
template class KSWKey<Rp12, 3, 1, B_Qp, K_Qp>;
template class KSWKey<Rz, P1, N, B_Qp_2, K_Qp_2>;
