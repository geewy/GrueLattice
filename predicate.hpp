#ifndef GL_PREDICATE_H
#define GL_PREDICATE_H

/**
  * @file predicate.hpp
  * @date may 2018
  * @author geewy
  * @brief Declaration of the class Fun
  *
  * The Eva object stores everything needed to compute the function. The most important aspect is its operator ()
  */

#include <iostream>
#include "parameters.hpp"
#include "integer_mod.hpp"

/**
 * @class Eva
 * @brief The container with its call operator
 */
class Eva{
	public:
		/**
		 * @brief Allows to evaluate the function represented bu the object
		 * @param point : the point where to evaluate the function
		 * @return F(point)
		 */
		inline Zt operator()(const Zp12 point) const{
			return eval((int64_t) std::floor((double) point * (double) T / (double) (P1*P2) + 0.5));
		}
		inline Zt eval(const Zt i) const{
			#if INPUT_BIT > 10
			/** 63 bits threshold **/
				int64_t input = (int64_t) i;
				int64_t result = (input < 0) ? 1 : 0;
			#else
			/** 6 bits parity **/
				int64_t input  = (int64_t) i;
				input = input & 0x3F;
				int64_t result = __builtin_popcount(input) % 2;
			#endif
				Zt res(result);
				return res;
		}
};

#endif
