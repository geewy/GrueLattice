/**
  * @file rlwe.cpp
  * @date may 2018
  * @author geewy
  * @brief Declaration of RLWE encryption and other operations
  *
  * See rlwe.hpp for more details
  */

#include <memory>
#include <cmath>
#include "rlwe.hpp"
#include "parameters.hpp"
#include "gadget.hpp"
#include "operations.hpp"
#include "ksw_key.hpp"

template<class Ring, uint64_t dim>
RLWE<Ring, dim>::RLWE(){
    array = new Ring[dim];
}

template<class Ring, uint64_t dim>
RLWE<Ring, dim>::RLWE(const RLWE &src){
    array = new Ring[dim];
    std::copy(src.array, src.array + dim, array);
    elem = src.elem;
}

template<class Ring, uint64_t dim>
RLWE<Ring, dim>::RLWE(Ring *new_array, Ring new_elem){
    array = new Ring[dim];
    
    if(new_array != nullptr)
        std::copy(new_array, new_array + dim, array);
        
    elem = new_elem;
}

template<class Ring, uint64_t dim>
RLWE<Ring, dim>::RLWE(Ring *new_array, bool copy){
    elem = new_array[dim];
    if(copy){
        array = new Ring[dim];
        std::copy(new_array, new_array + dim, array);
    }
    else
		array = std::move(new_array); 
}

template<class Ring, uint64_t dim>
RLWE<Ring, dim>::~RLWE(){
    delete[] array;
}

template<class Ring, uint64_t dim>
template <class Zp>
void RLWE<Ring, dim>:: mod(Zp* array_b) const{
    for(size_t i = 0 ; i < dim ; ++i)
        array_b[i] = Zp((int64_t) array[i].get_data()[0]);
        
    array_b[dim] = Zp((int64_t) elem.get_data()[0]);
}

template<class Ring, uint64_t dim>
template<class Rt>
void RLWE<Ring, dim>::encrypt(const Ring *const secret, const Rt &message, const uint64_t encMod, const uint64_t plainMod, const double variance){
    const Ring error = Ring::sample_e(variance);
    
    for(size_t i = 0 ; i < dim ; ++i)
        array[i] = Ring::sample_a();
        
    size_t size = Rt::get_degree();
    int64_t *elem_data = new int64_t[size];
    
    for(size_t i = 0 ; i < size ; ++i)
        elem_data[i] = floor((double)(message.get_data()[i]) * (double) encMod / (double) plainMod + 0.5);
        
    elem = Ring(elem_data);
    delete[] elem_data;
    elem += dot_product<Ring, dim>(secret, array) + error;
    return;
}

template<class Ring, uint64_t dim>
template<class Rt>
void RLWE<Ring, dim>::decrypt(const Ring *const secret, Rt &dec, const uint64_t encMod, const uint64_t plainMod) const{
    dec = (elem - dot_product<Ring, dim>(secret, array)).exact_rounding(encMod, plainMod);
}

template<class Ring, uint64_t dim>
double RLWE<Ring, dim>::noise(const Ring *const secret, const uint64_t encMod, const uint64_t plainMod) const{
    Ring res = elem - dot_product<Ring, dim>(secret, array);
    Ring signed_noise = res - (encMod/plainMod) * res.exact_rounding(encMod, plainMod);
    return signed_noise.get_norm();
}

template<class Ring, uint64_t dim>
template<class Ring_p, uint64_t degree, uint64_t degree_fft>
RLWE<Ring_p, dim> RLWE<Ring, dim>::trace() const{
    Ring_p new_array[dim + 1];
    
    for(size_t i = 0 ; i < dim ; ++i)
        new_array[i] = array[i].template trace<degree, degree_fft>();
        
    new_array[dim] = elem.template trace<degree, degree_fft>();
    return RLWE<Ring_p, dim>(new_array);
}

template<class Ring, uint64_t dim>
template <class Ring2, uint64_t new_mod, uint64_t old_mod>
RLWE<Ring2, dim> RLWE<Ring, dim>::mod_switch() const{
    Ring2 array_dst[dim], elem_dst;
    
    for(size_t i = 0 ; i < dim ; ++i)
        array_dst[i] = array[i].template rounding<Ring2>(old_mod, new_mod);
        
    elem_dst = elem.template rounding<Ring2>(old_mod, new_mod);
    return RLWE<Ring2, dim>(array_dst, elem_dst);
}

template<class Ring, uint64_t dim>
template <uint64_t dimP, uint64_t Base, uint64_t Ksize>
RLWE<Ring, dimP> RLWE<Ring, dim>::key_switch(const KSWKey<Ring, dim, dimP, Base, Ksize> &ksw) const{
    Ring array_decomp[dim * Ksize];
    
    for(size_t i = 0 ; i < dim ; ++i)
        Gadget<Ring, Base, Ksize>::g_invT(array_decomp + i * Ksize, array[i]);
        
    Ring *res = Ring::keymult(array_decomp, ksw);
    
    for(size_t i = 0 ; i < dimP ; ++i)
        res[i] = - res[i];
        
    res[dimP] = elem - res[dimP];
    return RLWE<Ring, dimP>(std::move(res), false);
}

template<class Ring, uint64_t dim>
RLWE<Ring, dim> RLWE<Ring, dim>::ext_mult(const RGSW<Ring> &rgsw) const{
    RLWE<Ring, dim> result(*this);
    result.ext_mult_inplace(rgsw);
    return result;
}

template<class Ring, uint64_t dim>
void RLWE<Ring, dim>::ext_mult_inplace(const RGSW<Ring> &rgsw){
    Ring vec[(dim + 1) * K_Q];
    
    for(size_t i = 0 ; i < dim ; ++i)
        Gadget<Ring, B_Q, K_Q>::g_invT(vec + i * K_Q, U_inv * array[i]);
        
    Gadget<Ring, B_Q, K_Q>::g_invT(vec + dim * K_Q, U_inv * elem);
    Ring *res = vec * rgsw;
    std::copy(res, res+dim, array);
    elem = res[dim];
    delete[] res;
}

template<class Ring, uint64_t dim>
void RLWE<Ring, dim>::galois_inplace(const uint64_t alpha){
    for(size_t i = 0 ; i < dim ; ++i)
        array[i].galois_inplace(alpha);
        
    elem.galois_inplace(alpha);
}

template<class Ring, uint64_t dim>
RLWE<Ring, dim> RLWE<Ring, dim>::galois(const uint64_t alpha) const{
    RLWE<Ring, dim> result(*this);
    result.galois_inplace(alpha);
    return result;
}

template<class Ring, uint64_t dim>
RLWE<Ring, dim> &RLWE<Ring, dim>::operator=(const RLWE &other){
    if(this != &other){
        std::copy(other.array, other.array + dim, array);
        elem = other.elem;
    }
    return *this;
}

template<class Ring, uint64_t dim>
RLWE<Ring, dim> &RLWE<Ring, dim>::operator=(RLWE && other) noexcept{
    if(this != &other){
        delete[] array;
        Ring *old_value = std::move(other.array);
        other.array = std::forward<Ring *>(nullptr);
        array = old_value;

        elem = other.elem;
    }
    return *this;
}

template<class Ring>
void RLWE1<Ring>::ext_exp_mult_add(const RGSW<Ring> &rgsw, const uint64_t alpha, const KSWKey<Ring, 1, 1, B_Q, K_Q> &key_alpha, const uint64_t beta, const KSWKey<Ring, 1, 1, B_Q, K_Q> &key_beta){
    if(beta != 1){
        galois_inplace(beta);
        key_switch_inplace<B_Q, K_Q>(key_beta);
    }
    ext_mult_inplace(rgsw);
    if(alpha != 1){
        galois_inplace(alpha);
        key_switch_inplace<B_Q, K_Q>(key_alpha);
    }
}
template<class Ring>
template <uint64_t Base, uint64_t Ksize>
void RLWE1<Ring>::key_switch_inplace(const KSWKey<Ring, 1, 1, Base, Ksize> &ksw){
    Ring array_decomp[1 * Ksize];
    for(size_t i = 0 ; i < 1 ; ++i)
        Gadget<Ring, Base, Ksize>::g_invT(array_decomp + i * Ksize, array[i]);
    Ring *res = Ring::keymult(array_decomp, ksw);

    array[0] = -res[0];
    elem -= res[1];
    delete[] res;
}

template class RLWE<Rz, P1>;
template class RLWE<Rz, N>;
template class RLWE<CirculantRing<Zp12, 1, 1>, N>;
template class RLWE<Rp1, 1>;
template class RLWE<Rp2, 1>;
template class RLWE<Rp1_crt, 1>;
template class RLWE<Rp2_crt, 1>;
template class RLWE<Rp12_crt, 3>;
template class RLWE<Rp12, 3>;
template class RLWE<Rp12, 1>;

/** LWE **/
template RLWE<Rz, N> RLWE<Rz, P1>::key_switch<N, B_Qp_2, K_Qp_2>(const KSWKeyLWE &ksw) const;
template RLWE<CirculantRing<Zp12, 1, 1>, N> RLWE<Rz, N>::mod_switch<CirculantRing<Zp12, 1, 1>, P1*P2, Qp>() const;

/** ExtExpInner **/
template void RLWE<CirculantRing<Zp12, 1, 1>, N>::mod<Zp1>(Zp1 *cipher) const;
template void RLWE<CirculantRing<Zp12, 1, 1>, N>::mod<Zp2>(Zp2 *cipher) const;
template void RLWE1<Rp1>::key_switch_inplace<B_Q, K_Q>(const KSWKeyRp1 &ksw);
template void RLWE1<Rp2>::key_switch_inplace<B_Q, K_Q>(const KSWKeyRp2 &ksw);
template void Rp1LWE::ext_exp_mult_add(const Rp1GSW &rgsw, const uint64_t alpha, const KSWKeyRp1 &key_alpha, const uint64_t beta,  const KSWKeyRp1 &key_beta);
template void Rp2LWE::ext_exp_mult_add(const Rp2GSW &rgsw, const uint64_t alpha, const KSWKeyRp2 &key_alpha, const uint64_t beta,  const KSWKeyRp2 &key_beta);

/** ExpCRT **/
template RLWE<Rp1_crt, 1> RLWE<Rp1, 1>::mod_switch<Rp1_crt, Qcrt, Q>() const;
template RLWE<Rp2_crt, 1> RLWE<Rp2, 1>::mod_switch<Rp2_crt, Qcrt, Q>() const;
template Rp12LWE RLWE<Rp12_crt, 3>::mod_switch<Rp12, Qp, Qcrt>() const;

/** FunExpExtract **/
template RLWE<Rp12, 1> Rp12LWE::key_switch<1, B_Qp, K_Qp>(const KSWKeyRp12 &ksw) const;
template RLWE<Rp1_p, 1> RLWE<Rp12 , 1>::trace<Rp1_p, P1, FFT_DIM1>() const;

/** input/output **/
template void LWE::encrypt<CirculantRing<Zt, 1, 1> >(const Rz *secret, const CirculantRing<Zt, 1, 1> &message, const uint64_t encMod, const uint64_t plainMod, const double variance);
template void LWE::decrypt<CirculantRing<Zt, 1, 1> >(const Rz *secret, CirculantRing<Zt, 1, 1> &dec, const uint64_t encMod, const uint64_t plainMod) const;

/** Key-switching keys **/
template void RLWE<Rz, N>::encrypt<Rz>(const Rz *secret, const Rz &message, const uint64_t encMod, const uint64_t plainMod, const double variance);
template void RLWE<Rp1, 1>::encrypt<Rp1>(const Rp1 *secret, const Rp1 &message, const uint64_t encMod, const uint64_t plainMod, const double variance);
template void RLWE<Rp2, 1>::encrypt<Rp2>(const Rp2 *secret, const Rp2 &message, const uint64_t encMod, const uint64_t plainMod, const double variance);
template void RLWE<Rp12, 1>::encrypt<Rp12>(const Rp12 *secret, const Rp12 &message, const uint64_t encMod, const uint64_t plainMod, const double variance);

/** For tests only **/
template void RLWE<Rp1, 1>::encrypt<CirculantRing<Zt, P1> >(const Rp1 *secret, const CirculantRing<Zt, P1> &message, const uint64_t encMod, const uint64_t plainMod, const double variance);
template void RLWE<Rp1, 1>::decrypt<CirculantRing<Zt, P1> >(const Rp1 *secret, CirculantRing<Zt, P1> &dec, const uint64_t encMod, const uint64_t plainMod) const;
template void Rp12LWE::encrypt<CirculantRing<Zt, P1*P2, FFT_DIM2> >(const Rp12 *secret, const CirculantRing<Zt, P1*P2, FFT_DIM2> &message, const uint64_t encMod, const uint64_t plainMod, const double variance);
template void Rp12LWE::decrypt<CirculantRing<Zt, P1*P2, FFT_DIM2> >(const Rp12 *secret, CirculantRing<Zt, P1*P2, FFT_DIM2> &dec, const uint64_t encMod, const uint64_t plainMod) const;
template void RLWE<Rz, N>::encrypt<CirculantRing<Zt, 1, 1> >(const Rz *secret, const CirculantRing<Zt, 1, 1>  &message, const uint64_t encMod, const uint64_t plainMod, const double variance);
template void RLWE<Rp2, 1>::encrypt<CirculantRing<Zt, P2> >(const Rp2 *secret, const CirculantRing<Zt, P2> &message, const uint64_t encMod, const uint64_t plainMod, const double variance);
template void RLWE<CirculantRing<Zp12, 1, 1>, N>::decrypt<CirculantRing<Zt, 1, 1> >(const CirculantRing<Zp12, 1, 1> *secret, CirculantRing<Zt, 1, 1> &dec, const uint64_t encMod, const uint64_t plainMod) const;
template void RLWE<Rp12, 1>::decrypt<CirculantRing<Zt, P1*P2, FFT_DIM2> >(const Rp12 *secret, CirculantRing<Zt, P1*P2, FFT_DIM2> &dec, const uint64_t encMod, const uint64_t plainMod) const;

/** For stats only **/
template void RLWE<Rp1_crt, 1>::encrypt<CirculantRing<Zt, P1> >(const Rp1_crt *secret, const CirculantRing<Zt, P1> &message, const uint64_t encMod, const uint64_t plainMod, const double variance);
template void RLWE<Rp2_crt, 1>::encrypt<CirculantRing<Zt, P2> >(const Rp2_crt *secret, const CirculantRing<Zt, P2> &message, const uint64_t encMod, const uint64_t plainMod, const double variance);
template void RLWE<Rp12_crt, 3>::encrypt<CirculantRing<Zt, P1*P2, FFT_DIM2> >(const Rp12_crt *secret, const CirculantRing<Zt, P1*P2, FFT_DIM2> &message, const uint64_t encMod, const uint64_t plainMod, const double variance);
