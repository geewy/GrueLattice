/**
 * @file rgsw.cpp
 * @date may 2018
 * @author geewy
 * @brief Declaration of RGSW encryption
 * 
 * See rgsw.hpp for more details
 */

#include "operations.hpp"
#include "rgsw.hpp"

template<class Ring>
template<class Rt>
void RGSW<Ring>::encrypt(const Ring *const secret, const Rt &message, const double variance){
	for(size_t i = 0; i < 2*K_Q; ++i){
		for(size_t j = 0; j < 2-1; ++j)
			Cipher[i][j] = Ring::sample_a();
		const Ring error = Ring::sample_e(variance);
		Cipher[i][2-1] = dot_product<Ring, 2-1>(secret, Cipher[i]) + error;
	}
	
	Ring mu;
	mu = U * Ring(message);
	for(size_t i = 0; i < K_Q; ++i){
		const Ring tmp = Matrix[i][0] * mu;
		for(size_t j = 0; j < 2; ++j)
			Cipher[i + K_Q * j][j] += tmp;
	}
	for(size_t i = 0; i < 2*K_Q; ++i)
		for(size_t j = 0; j < 2; ++j)
			Cipher[i][j].decomp_fft();
}

template<class Ring>
template<class Rt>
void RGSW<Ring>::decrypt(const Ring *const secret, Rt &dec) const{
	Ring curl_val = Cipher[2 * K_Q - (K_Q-1) - 1][2-1];
	curl_val.recomp_fft(SPLIT_Q);
	for(size_t j = 0; j < 2-1; ++j){
		Ring tmp = secret[j];
		tmp *= Cipher[2 * K_Q - (K_Q-1) - 1][j];
		curl_val -= tmp;
	}
	dec = curl_val.exact_rounding(Q, T);	
}

template class RGSW<Rp1>;
template class RGSW<Rp2>;

template void Rp1GSW::encrypt<CirculantRing<Zt, P1> >(const Rp1 *const secret, const CirculantRing<Zt, P1> &message, const double variance);
template void Rp2GSW::encrypt<CirculantRing<Zt, P2> >(const Rp2 *const secret, const CirculantRing<Zt, P2> &message, const double variance);

/** For tests only **/
template void Rp1GSW::decrypt<CirculantRing<Zt, P1> >(const Rp1 *const secret, CirculantRing<Zt, P1> &dec) const;
