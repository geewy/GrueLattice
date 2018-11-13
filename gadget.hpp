#ifndef GL_GADGET_H
#define GL_GADGET_H

/**
  * @file gadget.hpp
  * @date may 2018
  * @author geewy
  * @brief Definition of Gadget tool
  *
  * The Gadget allows to compute decomposition of integer elements into a basis "Base".
  * It has no internal data, just class methods.
  * The class is templated on the ring it operates on, the basis that is used for the decomposition and the 
  * size of these decompositions.
  */

#include <cmath>
#include <cstring>

/**
 * @class Gadget
 * @brief The tool for decomposition  of elements in Ring smaller than Base^{Ksize-1} into basis Base
 * @tparam Ring : the class representing elements
 * @tparam Base : the integer basis for the decomposition
 * @tparam Ksize :  the size of the decomposition
 */
template<class Ring, uint64_t Base, uint64_t Ksize>
class Gadget{
	public:
		/**
		 * @brief Computes the gadget matrix of size nSize : G_nSize = I_{nSize+1}(x)g
		 * @tparam nSize : the size of the gadget
		 * @return the gadget matrix
		 */
		template<uint64_t nSize>
		static void gadget(int64_t Gadget_G[nSize*Ksize][nSize]){
			for(size_t i = 0; i < nSize * Ksize; ++i)
				for(size_t j = 0; j < nSize; ++j)
					Gadget_G[i][j] = 0;
			int64_t Bexpo = 1;
			for(size_t i = 0; i < Ksize; ++i){
				for(size_t j = 0; j < nSize; ++j)
					Gadget_G[i + Ksize*j][j] = Bexpo;
				Bexpo *= Base;
			}
		}
		/**
		 * @brief Compute the decomposition
		 * @param decom : the returned decomposition
		 * @param val : the element to decompose
		 */
		static void g_invT(Ring *decomp, const Ring &val){
			Ring temp(val);
			for(size_t expo = 0; expo < Ksize; ++expo){
				decomp[expo] = temp % Base;
				decomp[expo].balance(Base);
				temp -= decomp[expo];
				temp = temp / Base;
			}
		}
		 /**
		 * @brief Compute the decomposition and outputs the results in FFT
		 * @param dec_fft the returned decomposition in FFT
		 * @param val the element to decompose
		 */
		template <class Ring2>
		static void g_invT_fft(fftw_complex *dec_fft[Ksize], const Ring2 &val){
			Ring2 decomp[Ksize];
			g_invT(decomp, val);
			for (size_t i = 0 ; i < Ksize ; ++i)
				decomp[i].compute_fft(dec_fft[i]);
		}
		/**
		 * @brief Accessor to the basis for the decomposition
		 * @return Base
		 */
		static inline uint64_t get_basis(){return Base;}
		/**
		 * @brief Accessor to the size of the decomposition
		 * @return K
		 */
		static inline uint64_t get_size(){return Ksize;}
};
#endif

