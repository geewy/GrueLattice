#ifndef GL_RGSW_H
#define GL_RGSW_H

/**
 * @file rgsw.hpp
 * @date may 2018
 * @author geewy
 * @brief Declaration of RGSW encryption
 * 
 * This class defines how to encrypt and decrypt following the symetric GSW encryption scheme over ring.
 * The object holds the ciphertext as a matrix.
 * For each template setting, the class also includes the gadget matrix.
 */

#include <iostream>
#include "gadget.hpp"
#include "parameters.hpp"

template<class Ring, uint64_t dim>
class RLWE;

/**
 * @class RGSW
 * @brief Store Ring-GSW ciphertexts and provide encryption/decryption routines
 * @tparam Ring : the class for representing its elements
 */
template<class Ring>
class RGSW{
	private:
		static int64_t Matrix[2*K_Q][2];	/*! < the gadget matrix for fixed class */
		Ring Cipher[2*K_Q][2];	/*! < the ciphertext */
	public:
		/**
		 * @brief Initialize the gadget matrix for the class in use. To be called once per setting by the user
		 */
		static void Gadget_init(){Gadget<Ring, B_Q, K_Q>::template gadget<2>(Matrix);}
		/**
		 * @brief Encryption function. It stores the ciphertext in the object
		 * @tparam Rt : the plaintext set
		 * @param secret : the secret to use for the encryption
		 * @param message : the plaintext data to encrypt
		 */
		template<class Rt>
		void encrypt(const Ring *const secret, const Rt &message, const double variance);
		/**
		 * @brief Decryption function. It decrypts the internal ciphertext
		 * @tparam Rt : the plaintext set
		 * @param secret : the secret to use for the decryption
		 * @param dec : receives the decryption of the internal ciphertext
		 */
		template<class Rt>
		void decrypt(const Ring *const secret, Rt &dec) const;
		inline const Ring &operator () (const size_t ite1, const size_t ite2) const{return Cipher[ite1][ite2];}
		friend std::ostream &operator << (std::ostream &os, const RGSW &obj){
			for(size_t i = 0; i < 2*K_Q; ++i){
				for(size_t j = 0; j < 2; ++j)
					os << obj.A[i][j] << " ";
				os << std::endl;
			}
			return os;
		}
};

typedef RGSW<Rp1> Rp1GSW;
typedef RGSW<Rp2> Rp2GSW;

#endif
