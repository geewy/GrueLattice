#ifndef GL_RLWE_H
#define GL_RLWE_H

/**
  * @file rlwe.hpp
  * @date may 2018
  * @author geewy
  * @brief Declaration of RLWE encryption and other operations
  *
  * This class defines how to encrypt and decrypt following the symmetric LWE encryption scheme over ring.
  * It can be used for LWE (ring of degree 1), Ring-LWE (ring of arbitrary degree) and n-dimensionnal variant of Ring-LWE 
  * (array of ring elements)
  * The object holds the ciphertext c=(a, b) as an array of ring elements for a (of size 1 for LWE/Ring-LWE), 
  * plus a ring element for b.
  * It is templated with the ring, the ciphetext modulus and the dimension
  */

#include <iostream>
#include "circulant_ring.hpp"
#include "rgsw.hpp"

template<class Ring1, uint64_t dim1, uint64_t dimP, uint64_t Base, uint64_t Ksize>
class KKSWKey{};

/**
 * @class RLWE
 * @brief Class for encryption/decryption and other associated operations
 * @tparam Ring : the class for the elements (=Z for LWE)
 * @tparam dim : the ciphertext dimension
 */
template<class Ring, uint64_t dim>
class RLWE{
	protected:
		Ring *array;	/*! < array of ring elements of size dim */
		Ring elem;	/*! < one ring element */
	public:
		/**
		 * @brief Default constructor
		 */
		RLWE();
		/**
		 * @brief Copy constructor
		 * @param src : RLWE object to duplicate
		 */
		RLWE(const RLWE &src);
		/**
		 * @brief Constructs an RLWE encryption from a representation
		 * @param new_array : an array of size dim
		 * @param new_elem : a ring element
		 */
		RLWE(Ring *new_array, Ring new_elem);
		/**
		 * @brief Constructs an RLWE encryption from a representation
		 * @param new_array : an array of size dim+1, elem is in the last coefficient
		 */
		RLWE(Ring *new_array, bool copy=true);
		 /**
		 * @brief Default destructor
		 */
		~RLWE();
		/**
		 * @brief Accessor to the internal representation
		 * @return a pointer to the internal
		 */
		inline Ring *get_array() const {return array;}
		/**
		 * @brief Accessor to the internal representation
		 * @return the ring element elem
		 */
		inline Ring get_elem() const {return elem;}
		/**
		 * @brief Encryption function. It stores the ciphertext in the object.
		 * @tparam Rt : the plaintext set
		 * @param secret : the secret to use for the encryption
		 * @param message : the plaintext data to encrypt
		 * @param encMod : ciphertext modulus
		 * @param plainMod : plaintext modulus
		 * @param variance : variance for the noise
		 */
		template<class Rt>
		void encrypt(const Ring *const secret, const Rt &message, const uint64_t encMod, const uint64_t plainMod, const double variance);
		/**
		 * @brief Decryption function. It decrypts the internal ciphertext
		 * @tparam Rt : the plaintext set
		 * @param secret : the secret to use for the decryption
		 * @param dec : receive the decryption of the internal ciphertext
		 * @param encMod : ciphertext modulus
		 * @param plainMod : plaintext modulus
		 */
		template <class Rt>
		void decrypt(const Ring *const secret, Rt &dec, const uint64_t encMod, const uint64_t plainMod) const;
		/**
		 * @brief Return the noise of the ciphertext
		 * @param secret : the secret key
		 * @param encMod : ciphertext modulus
		 * @param plainMod : plaintext modulus
		 * @return the noise magnitude
		 */
		double noise(const Ring *const secret, const uint64_t encMod, const uint64_t plainMod) const;
		/**
		 * @brief Computes the Trace of each internal ring elements over another ring of degree "degree" (templated)
		 * @tparam Ring_p : the other ring
		 * @tparam degree : the degree of the other ring
		 * @tparam degree_fft : ftt equivalent of the ring's degree
		 * @return a RLWE encryption over the ring Ring_p
		 */
		template<class Ring_p, uint64_t degree, uint64_t degree_fft>
		RLWE<Ring_p, dim> trace() const;
		/**
		 * @brief Compute the modulo on each internal elements, makes sense for LWE only
		 * @tparam Zp : the class for the resulting elements
		 * @param array_b : an array of size dim+1 with ring elements modulo rhs
		 */
		template <class Zp>
		void mod(Zp* array_b) const;
		/**
		 * @brief ModSwitch operation (out-of-place)
		 * @tparam Ring2 : the new ciphertext set
		 * @tparam new_mod : the new ciphertext modulus
		 * @tparam old_mod : the old ciphertext modulus
		 * @return an RLWE encryption of the same plaintext under the new modulus
		 */
		template<class Ring2, uint64_t new_mod, uint64_t old_mod>
		RLWE<Ring2, dim> mod_switch() const;
		/**
		 * @brief KeySwitch operation (out-of-place)
		 * @tparam dimP : the dimension of the destination modulus
		 * @tparam Base : the basis for the decomposition
		 * @tparam Ksize : the size of the decompositions
		 * @param ksw : the Key-Switching key
		 * @return an RLWE encryption of the same plaintext under the new key
		 */
		template <uint64_t dimP, uint64_t Base, uint64_t Ksize>
		RLWE<Ring, dimP> key_switch(const KSWKey <Ring, dim, dimP, Base, Ksize> &ksw) const;
		/**
		 * @brief ExtMult operation (out-of-place)
		 * @param rgsw : the RGSW ciphertext to multiply with
		 * @return the RLWE encryption result of the external multiplication
		 */
		RLWE ext_mult(const RGSW<Ring> &rgsw) const;
		/**
		 * @brief ExtMult operation (in-place)
		 * @param rgsw : the RGSW ciphertext to multiply with
		 */
		void ext_mult_inplace(const RGSW<Ring> &rgsw);
		/**
		 * @brief Applies the Galois morphism to each internal ring elements (out-of-place)
		 * @param alpha : the morphism is psi_alpha : T -> T^alpha
		 * @return the resulting RLWE encryption
		 */
		RLWE galois(const uint64_t alpha) const;
		/**
		 * @brief Applies the Galois morphism to each internal ring elements (in-place)
		 * @param alpha : the morphism is psi_alpha : T -> T^alpha
		 */
		void galois_inplace(const uint64_t alpha);
		/**
		 * @brief Fallback generic multiplication function for cases where ring elements are small nor in FFT
		 * @param rhs : the ring element to multiply with
		 */
		inline void mult(const fftw_complex *rhs){
			for(size_t i = 0; i < dim; ++i)
				array[i].mult(rhs);
			elem.mult(rhs);
		}
		/**
		 * @brief Fallback generic multiplication function for cases where ring elements are small nor in FFT
		 * @param rhs the ring element to multiply with
		 */
		inline void mult(const Ring &rhs){
			for(size_t i = 0 ; i < dim; ++i)
				array[i].mult(rhs);
			elem.mult(rhs);
		}
		 /**
		 * @brief Compute FFT from data for each element in internal representation
		 */
		inline void decomp_fft(){
			for(size_t i = 0; i < dim; ++i)
				array[i].decomp_fft();
			elem.decomp_fft();
		}
		/**
		 * @brief Compute data from fft for each element in internal representation
		 */
		inline void recomp_fft(){
			for(size_t i = 0; i < dim; ++i)
				array[i].recomp_fft(SPLIT_Q);
			elem.recomp_fft(SPLIT_Q);
		}
		inline RLWE &operator += (const RLWE &rhs){
			for(size_t i = 0; i < dim; ++i)
				array[i] += rhs.array[i];
			elem += rhs.elem;
			return *this;
		}
		inline RLWE &operator -= (const RLWE &rhs){
			for(size_t i = 0; i < dim; ++i)
				array[i] -= rhs.array[i];
			elem -= rhs.elem;
			return *this;
		}
		inline RLWE &operator *= (const Ring &rhs){
			for(size_t i = 0; i < dim; ++i)
				array[i] *= rhs;
			elem *= rhs;
			return *this;
		}
		inline RLWE &operator *= (const int64_t &rhs){
			for(size_t i = 0; i < dim; ++i)
				array[i] *= rhs;
			elem *= rhs;
			return *this;
		}
		
		RLWE &operator = (const RLWE &other);
		RLWE &operator = (RLWE && other) noexcept;
		friend inline RLWE operator + (RLWE lhs, const RLWE &rhs){lhs += rhs; return lhs;}
		friend inline RLWE operator - (RLWE lhs, const RLWE &rhs){lhs -= rhs; return lhs;}
		friend inline RLWE operator * (RLWE lhs, const int64_t &rhs){lhs *= rhs; return lhs;}
		friend inline RLWE operator * (RLWE lhs, const Ring &rhs){lhs *= rhs; return lhs;}
		friend inline RLWE operator * (const Ring &lhs, RLWE rhs){rhs *= lhs; return rhs;}
		friend std::ostream &operator << (std::ostream &os, const RLWE &obj){
			os << "(";
			for(size_t i = 0; i < dim; ++i)
				os << obj.array[i] << ", ";
			os << obj.elem << ")";
			return os;
		}		
};


template<class Ring2>
class RLWE1 : public RLWE<Ring2, 1>{
	private:
		using RLWE<Ring2, 1>::array;
		using RLWE<Ring2, 1>::elem;
	public:
		using RLWE<Ring2, 1>::ext_mult_inplace;
		using RLWE<Ring2, 1>::galois_inplace;
		using RLWE<Ring2, 1>::RLWE;
		/**
		 * @brief KeySwitch operation (in-place). The new and old keys have the same dimension
		 * @tparam Base : the basis for the decomposition
		 * @tparam Ksize : the size of the decompositions
		 * @param ksw : the Key-Switching key
		 */
		template <uint64_t Base, uint64_t Ksize>
		void key_switch_inplace(const KSWKey <Ring2, 1, 1, Base, Ksize> &ksw);
		 /**
		 * @brief ExtExpMultAdd operation (in-place)
		 * @param rgsw : the RGSW encryption to multiply and add
		 * @param alpha : the coefficient for the multiplication
		 * @param key_alpha : a Key-Switching key from psi_alpha(s) to s
		 * @param beta : inverse of alpha
		 * @param key_beta : a Key-Switching key from psi_beta(s) to s
		 */
		void ext_exp_mult_add(const RGSW<Ring2> &rgsw, const uint64_t alpha, const KSWKey<Ring2, 1, 1, B_Q, K_Q> &key_alpha, const uint64_t beta, const KSWKey<Ring2, 1, 1, B_Q, K_Q> &key_beta);
};


typedef RLWE<Rz, P1> LWE;
typedef RLWE1<Rp1> Rp1LWE;
typedef RLWE1<Rp2> Rp2LWE;
typedef RLWE<Rp12, 3> Rp12LWE;

#endif
