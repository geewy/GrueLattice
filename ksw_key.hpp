#ifndef GL_KSW_KEY_H
#define GL_KSW_KEY_H

/**
  * @file ksw_key.hpp
  * @date may 2018
  * @author geewy
  * @brief Definition of Key-Switching key class
  *
  * The KSWkey object is a convenient container for a Key-Switching for secret key s to secret key s' (sp)
  * It contains only basic initialization functions and accessors and is extensively templated for performance.
  */

#include <assert.h>
#include "rlwe.hpp"

/**
 * @class KSWKey
 * @brief The key-Switching key from s to sp, aka RLWE encryptions under sp of scaled s.
 * @tparam Ring : the class representing ring elements
 * @tparam dim : the dimension of the source key s
 * @tparam dimP : the dimension of the destination key sp
 * @tparam Base : the integer basis for the decomposition
 * @tparam Ksize : the size of the decomposition
 */
template<class Ring1, uint64_t dim1 = 1, uint64_t dimP = 1, uint64_t Base = B_Q, uint64_t Ksize = K_Q>
class KSWKey{
	private:
		RLWE<Ring1, dimP>SECRET[dim1][Ksize];
		bool init_ok;
	public:
		KSWKey();
		KSWKey(const Ring1 secret[dim1], const Ring1 secret_p[dimP], const double variance);
		void init(const Ring1 secret[dim1], const Ring1 secret_p[dimP], const double variance);
		inline const RLWE<Ring1, dimP> &operator () (const size_t idx_n, const size_t idx_K) const{
			assert(init_ok);
			return SECRET[idx_n][idx_K];
		}
		inline uint64_t get_basis() const{return Base;}
		inline uint64_t get_size() const{return Ksize;}
};

typedef KSWKey<Rp1> KSWKeyRp1;
typedef KSWKey<Rp2> KSWKeyRp2;
typedef KSWKey<Rp12, 3, 1, B_Qp, K_Qp> KSWKeyRp12;
typedef KSWKey<Rz, P1, N, B_Qp_2, K_Qp_2> KSWKeyLWE;

#endif
