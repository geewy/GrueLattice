#ifndef GL_OPERATIONS_H
#define GL_OPERATIONS_H

/**
  * @file operations.hpp
  * @date may 2018
  * @author geewy
  * @brief Declaration of some elementary operations of the scheme
  *
  * Some operations are not directly included in RLWE or RGSW classes. They are regrouped here.
  */

#include "ksw_key.hpp"
#include "rlwe.hpp"
#include "predicate.hpp"

void print_param();

/**
 * @brief Compute the dot product of two vectors
 * @param vec1 the first vector
 * @param vec2 the second vector
 * @tparam Ring : the class representing vector elements
 * @tparam size : the size of the vectors
 * @return the dot product <vec1.vec2> in Ring
 */
template<class Ring_op, uint64_t size>
Ring_op dot_product(const Ring_op *const vec1, const Ring_op *const vec2);
/**
 * @brief ExpCRT operation
 * @param cp : the RLWE encryption in Rp
 * @param cq : the RLWE encryption in Rq
 * @param cpq : the computed RLWE encryption in Rpq
 */
void exp_crt(RLWE<Rp12_crt, 3> &cpq, const RLWE<Rp1_crt, 1> &cp, const RLWE<Rp2_crt, 1> &cq);
/**
 * @brief Computes the key for ciphertexts from ExpCRT
 * @param sp : the key of the RLWE encryption in Rp
 * @param sq : the key of the RLWE encryption in Rq
 * @param spq : the resulting key for the RLWE encryption in Rpq
 */
void crt_key(Rp12 *spq, const Rp1 &sp, const Rp2 &sq);
/**
 * @brief ExtExpInner operation (optimized)
 * @tparam Ring : the class of the ring elements
 * @tparam Yelem:  the class of elements of y
 * @tparam inMod : the modulus of the input elements
 * @param encMod the modulus of the encryptions
 * @param plainMod the plaintext modulus
 * @param length:  the length of the inner product
 * @param Yelem the l plaintext elements
 * @param rgswEnc the l RGSW encryptions of T^{x_i}
 * @param kswArray the array of key-switching keys
 * @return an RLWE encryption of T^{<x, y>}
 */
template<class Ring_op2, class Yelem, uint64_t inMod>
RLWE1<Ring_op2> ext_exp_inner(const uint64_t encMod, const uint64_t plainMod, const size_t length, const Yelem *plainElem, const RGSW<Ring_op2> *rgswEnc, const KSWKey<Ring_op2, 1, 1, B_Q, K_Q> *kswArray);
/**
 * @brief fun_exp_extract
 * @param func : the function to evaluate on m
 * @param c_Zm : a RLWE encryption in Rpq of Z^m
 * @param KSW : a Key-Switching key from s_pq in Rpq to s
 * @return an LWE encryption of dimension n of F(m)
 */
LWE fun_exp_extract(const fftw_complex *func, const Rp12LWE &c_Zm, const KSWKeyRp12 &KSW);
/**
 * @brief Accumulation, part of the preparation phase
 * @tparam Ring the ring of the key elements
 * @tparam Elem the class of the intermediary elements
 * @tparam R_E the ring of the resulting elements
 * @tparam mod the modulus
 * @param enc encryption of m mod pq
 * @param BK : Bootstrapping keys (RGSW encryptions of s_i)
 * @param KS : KeySwitching Keys phi_alpha(s)-> s for alpha in Zp1
 * @return an RLWE encryption of Z^{-<a.s> mod p}
 */
template<class Ring_op, class Elem, class R_E, uint64_t mod>
RLWE<R_E, 1> accumulation(const RLWE<CirculantRing<Zp12, 1, 1>, N> &enc, const RGSW<Ring_op> BK[N], const KSWKey<Ring_op> KS[mod]);
/**
 * @brief Preparation phase, from LWE(m) to Rp12LWE(Z^m)
 * @param enc : LWE encryption of m
 * @param S_LWE : KeySwitching Key for an LWE encryption
 * @param Xsi : Rp1GSW encryptions of s_i
 * @param KSp1 : KeySwitching Keys phi_alpha(s)-> s for alpha in Zp1
 * @param Ysi : Rp2GSW encryptions of s_i
 * @param KSp2 : KeySwitching Keys phi_alpha(s)-> s for alpha in Zp2
 * @return an Rp12LWE encryption of Z^m
 */
Rp12LWE preparation(const LWE &enc, const KSWKeyLWE &S_LWE, const Rp1GSW Xsi[N], const KSWKeyRp1 KSp1[P1], const Rp2GSW Ysi[N], const KSWKeyRp2 KSp2[P2]);
/**
 * @brief Combination phase, from LWE(m_k) to LWE(m) where m = coef_0*m_0 + coef_1*m_1 + coef_{k-1}*m_{k-1}
 * @param num : the number of element to combine
 * @param coefs : the k coefs for the combination
 * @param c_i : the elements
 * @return the combination of the inputs
 */
LWE combination(const size_t num, const int64_t *coefs, const LWE *c_i);
/**
 * @brief Compute the function F on the input data
 * @param num : the number of inputs
 * @param coefs : the coefficients for the combination
 * @param c_i t: he inputs (LWE encryptions)
 * @param func : the function to compute
 * @param Xsi : the bootstrapping keys in the first ring
 * @param KSp1 : the Key-switching keys in the first ring
 * @param Ysi : the bootstrapping keys in the second ring
 * @param KSp2 : the Key-switching keys in the second ring
 * @param KSW : the key-switching key for the function extraction
 * @return an LWE encryption of F(m_1, ..., m_k)
 */
LWE gate(const size_t num, const int64_t *coefs, const LWE *c_i, const fftw_complex *func, const KSWKeyLWE &S_LWE, const Rp1GSW Xsi[N],  const KSWKeyRp1 KSp1[P1], const Rp2GSW Ysi[N], const KSWKeyRp2 KSp2[P2], const KSWKeyRp12 &KSW);
/**
 * @brief Generate the bootstrapping keys for ExtExpInner
 * @tparam Ring : the class for the ring elements
 * @tparam size : the length of the LWE key
 * @tparam dim : the dimension of the key
 * @param mod : ciphertext modulus
 * @param plainMod : plaintext modulus
 * @param Tsi : the RGSW encryptions of T^s_i
 * @param key : the LWE key
 * @param key_p : the key for the RGSW encryptions
 * @param variance : the variance of the encryption
 */
template<class Ring_op, uint64_t size, uint64_t dim_op>
void gen_bootstrapping_keys(const uint64_t mod, const uint64_t plainMod, RGSW<Ring_op> Tsi[size], const Rz key[size], const Ring_op &key_p, const double variance);
/**
 * @brief Generate the key-switching keys for ExtExpInner
 * @tparam Ring : the class for the ring elements
 * @tparam dim : the dimension of the key
 * @param keys : the resulting keys
 * @param key_p the secret key for which to compute the KS keys
 */
template<class Ring_op, uint64_t dim_op>
void gen_keyswitching_keys(KSWKey<Ring_op, 1, 1, B_Q, K_Q> keys[dim_op], const Ring_op &key_p, const double variance);
/**
 * @brief Generate the key-switching key for FunExpExtract
 * @param KSW : the key-switching key
 * @param key_pq : the key in R_pq after ExpCRT
 * @param key : the LWE key
 */
void gen_funexpextract_key(KSWKeyRp12 *KSW, const Rp12 key_pq[3], const Rz key[P1], const double variance);

#endif
