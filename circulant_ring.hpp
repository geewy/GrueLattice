#ifndef GL_CIRCULANT_RING_H
#define GL_CIRCULANT_RING_H

/**
 * @file circulant_ring.hpp
 * @date may 2018
 * @author geewy
 * @brief Declaration of the CirculantRing class
 * 
 * This class provides everything needed to represent and operate with elements in a cirulant ring of degree d, 
 * such as R = Z[T]/(T^d-1).
 * The degree is templated for better performance. Several constructors and operators ared defined for easy use.
 * For each degree an FTT class is associated for efficient multiplication.
 */

#include <iostream>
#include <utility>
#include <memory>
#include <cstring>
#include <assert.h>
#include "fft.hpp"
#include "integer_mod.hpp"
#include "predicate.hpp"

template<class Ring>
class RGSW;

template<class Ring, uint64_t dim, uint64_t dimP, uint64_t Base, uint64_t Ksize>
class KSWKey;

/**
 * @class CirculantRing
 * @brief Class representing elements in a circulant ring of degree "degree"
 * @tparam degree : the dregree of the ring
 * @tparam degree_fft : ftt equivalent of the ring's degree
 */
template<class Int, uint64_t degree, uint64_t degree_fft = FFT_DIM1>
class CirculantRing{
	private:
		static FFT<degree_fft> fft;	/*! < FFT object associated to the CirculantRing class, for each degree */
		Int *elem;	/*! < Internal representation of the ring element */
		fftw_complex *msb_fft;
		fftw_complex *lsb_fft;
		bool fft_status;
		static Int tmp[degree];
		static fftw_complex fft_cache1[degree_fft/2+1] __attribute__ ((aligned(32)));
		static fftw_complex fft_cache2[degree_fft/2+1] __attribute__ ((aligned(32)));
		static fftw_complex fft_cache3[degree_fft/2+1] __attribute__ ((aligned(32)));
	public:
		/**
		 * @brief Default constructor, doing only memory allocation
		 */
		CirculantRing();
		/**
		 * @brief Copy constructor 
		 * @param src : ring element to duplicate
		 */
		CirculantRing(const CirculantRing &src);
		/**
		 * @brief Copy constructor where src has another element representation
		 * @param src : ring element to reinterpret and duplicate
		 */
		template<class Int2>
		CirculantRing(const CirculantRing<Int2, degree, degree_fft> &src);
		/**
		 * @brief Copy constructor where src has another element representation
		 * @param src : ring element to reinterpret and duplicate
		 */
		template<uint64_t degree_fft2>
		CirculantRing(const CirculantRing<Int, degree, degree_fft2> &src);
		/**
		 * @brief Constructs a ring element from a reprensation
		 * @param values : an array corresponding to the internal representation
		 * @param copy : indicates if the representation is copied to a new internal array or just references
		 */
		CirculantRing(Int *values, bool copy=true);
		/**
		 * @brief Constructs a ring element from a reprensation
		 * @param msb_fft_values : the FFT representation of the most significant part
		 * @param lsb_fft_values : the FFT representation of the least significant part
		 * @param copy : indicates if the representation is copied to a new internal array or just references
		 */
		CirculantRing(fftw_complex *msb_fft_values, fftw_complex *lsb_fft_values, bool copy=true);
		/**
		 * @brief Constructs a ring element from a constant term
		 * @param val : the constant term in the internal representation
		 */
		CirculantRing(const Int val);
		CirculantRing(const int64_t *values);
		CirculantRing(const Eva &E);
		/**
		 * @brief Defaul destructor
		 */
		~CirculantRing();
		/**
		 * @brief Accessor to the FFT dimension
		 * @return the FFT dimension
		 */
		static inline size_t get_degree_fft(){return degree_fft;}
		static inline size_t get_degree(){return degree;}
		/**
		 * @brief Accessor to the internal representation
		 * @return a pointer to the internal representation
		 */
		inline Int *get_data() const{
			assert(fft_status == false);
			return elem;
		}
		/**
		 * @brief Accessor to the msb FFT internal representation
		 * @return a pointer to the msb FFT internal representation
		 */
		inline fftw_complex *get_msb_fft() const{
			assert(fft_status);
			return msb_fft;
		}
		/**
		 * @brief Accessor to the lsb FFT internal representation
		 * @return a pointer to the lsb FFT internal representation
		 */
		inline fftw_complex *get_lsb_fft() const{
			assert(fft_status);
			return lsb_fft;
		}
		/**
		 * @brief Accessor to the internal flag fft_status
		 * @return its value
		 */
		inline bool get_fft_status() const{return fft_status;}
		/**
		 * @brief Compute the 2-norm of the coefficient representation
		 * @return the norm
		 */
		inline double get_norm() const{
			double res = 0;
			for(size_t i = 0; i < degree; ++i)
				res += (double) elem[i] * (double) elem[i];
			return res;
		}
		inline bool is_zero() const{
			bool result = true;
			for(size_t i = 0; i < degree && result; ++i)
				result = result && (elem[i] == 0);
			return result;
		}
		static inline FFT<degree_fft> *get_fft_obj(){return &fft;}
		/**
		 * @brief Class method for Gaussian sampling of ring element
		 * @param stddev : the standard deviation
		 * @return a ring element whose coefficients are drawn from a Gaussian distribution of standard deviation stddev
		 */
		static CirculantRing gaussian_sample(const double stddev);
		/**
		 * @brief Class method uniformly sampling a ring element
		 * @return a ring element whose coefficients are drawn uniformly
		 */
		static CirculantRing uniform_sample();
		 /**
		 * @brief Class method sampling a binary ring element
		 * @return a ring element with coefficients in {0, 1}
		 */
		static CirculantRing binary_sample();
		/**
		 * @brief Class method uniformly sampling a ring element a with contraint a(1) = 0
		 * @return a ring element whose d-1 coefficients are drawn uniformly and the last meets the contraint
		 */
		static CirculantRing sample_a();
		/**
		 * @brief Class method sampling a ring elements with d/3 1, d/3 -1 and 0 elsewhere with operator norm < 0.8 sqrt(d ln d)
		 * @return a ring element satisfying the contraints
		 */
		static CirculantRing sample_s(const double density = 1/3);
		 /**
		 * @brief Class method sampling a ring elements serving as error term in encryptions
		 * @return a ring element satisfying the contraints for errors
		 */
		static CirculantRing sample_e(const double variance = 8);
		/**
		 * @brief Compute the fft of the current
		 * @param elem_fft : the FFT coefficients
		 */
		void compute_fft(fftw_complex *elem_fft) const;
		 /**
		 * @brief Update elem from elem_fft
		 * @param elem_fft : the FFT coefficients
		 */
		inline void compute_fft_inv(fftw_complex *elem_fft){compute_fft_inv(elem_fft, elem);}
		/**
		 * @brief Compute the inverse FFT of elem_fft and store to values
		 * @param elem_fft : the FFT coefficients
		 * @param values the result of the inverse FFT
		 */
		void compute_fft_inv(fftw_complex *elem_fft, Int *values) const;
		/**
		 * @brief Split the polynomial into the most and least significant parts and compute the FFT of each
		 */
		void decomp_fft();
		/**
		 * @brief Regroup the polynomial into an internal element elem
		 * @param cons : the constant used to regroup
		 */
		void recomp_fft(const int64_t cons);
		/**
		 * @brief Split the ring element into two halves, one with the most significant part and the other with 
		 * the least significant part
		 * @param msb : the most significant part
		 * @param lsb : the least signigicant part
		 * @return the : threshold value between the part
		 */
		int64_t split(CirculantRing &msb, CirculantRing &lsb) const;
		/**
		 * @brief Recompose the ring element from two halves, one with the most significant part and the other with 
		 * the least significant part
		 * @param msb : the most significant part
		 * @param lsb : the least signigicant part
		 * @return the threshold value between the part
		 */
		void unsplit(CirculantRing &msb, CirculantRing &lsb);
		/**
		 * @brief Shift the values by a factor shifter from least to most significant
		 * @param shifter : the shifting factor
		 * @return the shifted values
		 */
		inline void shift(const int64_t shifter){
			for(size_t i = 0; i < degree; ++i)
				elem[i] = ((int64_t) elem[i] % shifter) * shifter;
		}
		/**
		 * @brief Scales the ring element by a factor new_modulus/old_modulus and round (out-of-place)
		 * @param old_modulus : denominator of the ratio
		 * @param new_modulus : numerator of the ratio
		 * @return the scaled ring element
		 */
		CirculantRing exact_rounding(const uint64_t old_modulus, const uint64_t new_modulus) const;
		/**
		 * @brief Scales the ring element by a factor new_modulus/old_modulus and probabilistically round it (out-of-place)
		 * @param old_modulus : denominator of the ratio
		 * @param new_modulus : numerator of the ratio
		 * @return the scaled ring element
		 */
		template<class Ring2>
		Ring2 rounding(const uint64_t old_modulus, const uint64_t new_modulus) const;
		/**
		 * @brief Computes the Trace of the ring element over another circulant ring of degree degreeP (templated)
		 * @tparam degreeP :  the degree of the other circulant ring
		 * @tparam degree2_fft : the fft dimension of the result ring element
		 * @return an element in CirculantRing<Int, degreeP>
		 */
		template<uint64_t degreeP, uint64_t degree2_fft>
		CirculantRing<Int, degreeP, degree2_fft> trace() const;
		/**
		 * @brief Applies the Galois morphism to the ring element (out-of-place)
		 * @param alpha : the morphism is psi_alpha : T -> T^alpha
		 * @return the resulting ring element
		 */
		CirculantRing galois(const uint64_t alpha) const;
		/**
		 * @brief Applies the Galois morphism to the ring element (in-place)
		 * @param alpha : the morphism is psi_alpha : T -> T^alpha
		 */
		void galois_inplace(const uint64_t alpha);
		/**
		 * @brief Computes the tensor product between the ring element of degree d and another of degree d2
		 * @tparam degree2 : the degree of the other ring
		 * @tparam degree2_fft : the fft dimension of the other ring
		 * @tparam degree3_fft : the fft dimension of the result ring
		 * @param rhs : a ring element in a circulant ring of degree degree2
		 * @return the tensor product in the circulant ring of degree degree x degree2
		 */
		template<uint64_t degree2, uint64_t degree2_fft, uint64_t degree3_fft>
		CirculantRing<Int, degree*degree2, degree3_fft> tensor_product(const CirculantRing<Int, degree2, degree2_fft> &rhs) const;
		/**
		 * @brief Balances the coefficients between -B/2 and B/2
		 * @param range : the range
		 */
		inline void balance(const uint64_t range){
			for(size_t i = 0; i < degree; ++i)
				elem[i].balance(range);
		}
		 /**
		 * @brief Fallback generic multiplication function for cases where ring elements are in one FFT
		 * @param f : the ring element to multiply with
		 */
		void mult(const fftw_complex *f);
		/**
		 * @brief Fallback generic multiplication function for cases where ring elements are small nor in FFT
		 * @param rhs : the ring element to multiply with
		 */
		void mult(const CirculantRing &rhs);

		inline CirculantRing &operator += (const CirculantRing &rhs){
			for(size_t i = 0; i < degree; ++i)
				elem[i] += rhs.elem[i];
			return *this;
		}
		inline CirculantRing &operator -= (const CirculantRing &rhs){
			for(size_t i = 0; i < degree; ++i)
				elem[i] -= rhs.elem[i];
			return *this;
		}
		inline CirculantRing &operator *= (const uint64_t &rhs){
			const Int tmp(rhs);
			for(size_t i = 0; i < degree; ++i)
				elem[i] *= tmp;
			return *this;
		}
		CirculantRing &operator = (const CirculantRing &other);
		CirculantRing &operator = (CirculantRing &&other) noexcept;
		CirculantRing &operator = (const int &other);
		CirculantRing &operator *= (const CirculantRing &rhs);
		friend inline CirculantRing operator + (CirculantRing lhs, const CirculantRing &rhs) { lhs += rhs; return lhs; }
		friend inline CirculantRing operator - (CirculantRing lhs, const CirculantRing &rhs) { lhs -= rhs; return lhs; }
		friend inline CirculantRing operator * (CirculantRing lhs, const CirculantRing &rhs) { lhs *= rhs; return lhs; }
		friend inline bool operator != (const CirculantRing &lhs, const CirculantRing &rhs) { return !(lhs == rhs); }
		friend inline CirculantRing operator - (CirculantRing lhs) { return Int(0) - lhs; }
		friend inline CirculantRing operator * (const int64_t &lhs, const CirculantRing &rhs) { return rhs * lhs; }
		friend CirculantRing operator * (const CirculantRing &lhs, const int64_t &rhs){
			Int *res = new Int[degree];
			for(size_t i = 0 ; i < degree ; ++i)
				res[i] = lhs.elem[i] * rhs;
			return CirculantRing(res, false);
		}

		friend CirculantRing operator / (const CirculantRing &lhs, const int64_t &rhs){
			Int *res = new Int[degree];
			for(size_t i = 0 ; i < degree ; ++i)
				res[i] = lhs.elem[i] / rhs;
			return CirculantRing(res, false);
		}
		friend inline CirculantRing operator % (CirculantRing lhs, const int64_t &rhs){
			for(size_t i = 0 ; i < degree ; ++i)
				lhs.elem[i] = lhs.elem[i] % rhs;  
			return lhs;
		}
		friend inline bool operator == (const CirculantRing &lhs, const CirculantRing &rhs){
			bool res = true;
			for(size_t i = 0 ; i < degree  && res ; ++i)
				res = res && (lhs.elem[i] == rhs.elem[i]);
			return res;
		}
		friend std::ostream &operator << (std::ostream &os, const CirculantRing<Int, degree, degree_fft> &obj){
			os << "(";
			for(size_t i = 0; i < degree-1; ++i)
				os << obj.elem[i] << ",";
			os << obj.elem[degree-1] << ")" << std::flush;
			return os;
		}
		/**
		 * @brief Vector-Matrix Multiplication
		 * @param val : the row vector to multiply to the ciphertext matrix
		 * @param rhs : the ciphertext
		 * @return a vector row vector equal to val * rhs
		 */
		friend CirculantRing *operator * (CirculantRing *val, const RGSW<CirculantRing> &rhs){
			static const int64_t thresh = Int::get_split_threshold();
			CirculantRing *res = new CirculantRing[2];
			CirculantRing res_lsb[2];
			fftw_complex val_fft[2* K_Q][degree_fft/2+1];
			
			for(size_t j = 0 ; j < 2*K_Q ; ++j)
				val[j].compute_fft(val_fft[j]);

			for(size_t i = 0 ; i < 2 ; ++i){
				std::memset(fft_cache1, 0, (degree_fft/2+1) * 2 * sizeof(double));
				std::memset(fft_cache2, 0, (degree_fft/2+1) * 2 * sizeof(double));
				for(size_t j = 0 ; j < 2*K_Q ; ++j){
					fft::product<degree_fft>(fft_cache3, val_fft[j], rhs(j, i).get_msb_fft());
					for(size_t k = 0 ; k < degree_fft/2+1 ; ++k){
						fft_cache1[k][0] += fft_cache3[k][0];
						fft_cache1[k][1] += fft_cache3[k][1];
					}
					fft::product<degree_fft>(fft_cache3, val_fft[j], rhs(j, i).get_lsb_fft());
					for(size_t k = 0 ; k < degree_fft/2+1 ; ++k){
						fft_cache2[k][0] += fft_cache3[k][0];
						fft_cache2[k][1] += fft_cache3[k][1];
					}
				}
				res[i].compute_fft_inv(fft_cache1);
				res[i].shift(thresh);

				res_lsb[i].compute_fft_inv(fft_cache2);
				res[i] += res_lsb[i];
			}
			return res;
		}
		
		/**
		 * @brief Vector-Matrix Multiplication
		 * @param val : the row vector to multiply to the ciphertext matrix
		 * @param rhs : the ciphertext
		 * @return a vector row vector equal to val * rhs
		 */
		template <uint64_t dim, uint64_t dimP, uint64_t Base, uint64_t Ksize>
		static CirculantRing *keymult(CirculantRing *val, const KSWKey<CirculantRing, dim, dimP, Base, Ksize> &rhs){
			static const int64_t thresh = Int::get_split_threshold();
			CirculantRing *res = new CirculantRing[dimP + 1];
			CirculantRing res_lsb;
			fftw_complex *val_fft[dim * Ksize];
			
			for(size_t j = 0 ; j < dim * Ksize ; ++j){
				val_fft[j] = fftw_alloc_complex(degree_fft/2+1);
				val[j].compute_fft(val_fft[j]);
			}

			for(size_t idx = 0 ; idx < dimP ; ++idx){
				std::memset(fft_cache1, 0, (degree_fft/2+1) * 2 * sizeof(double));
				std::memset(fft_cache2, 0, (degree_fft/2+1) * 2 * sizeof(double));
				for(size_t i = 0 ; i < dim ; ++i){
					for(size_t j = 0 ; j < Ksize ; ++j){
						fft::product<degree_fft>(fft_cache3, val_fft[i * Ksize + j], rhs(i, j).get_array()[idx].get_msb_fft());
						for(size_t k = 0 ; k < degree_fft/2+1 ; ++k){
							fft_cache1[k][0] += fft_cache3[k][0];
							fft_cache1[k][1] += fft_cache3[k][1];
						}
						fft::product<degree_fft>(fft_cache3, val_fft[i * Ksize + j], rhs(i, j).get_array()[idx].get_lsb_fft());
						for(size_t k = 0 ; k < degree_fft/2+1 ; ++k){
							fft_cache2[k][0] += fft_cache3[k][0];
							fft_cache2[k][1] += fft_cache3[k][1];
						}
					}
				}
				res[idx].compute_fft_inv(fft_cache1);
				res[idx].shift(thresh);

				res_lsb.compute_fft_inv(fft_cache2);
				res[idx] += res_lsb;
			}
			{
				std::memset(fft_cache1, 0, (degree_fft/2+1) * 2 * sizeof(double));
				std::memset(fft_cache2, 0, (degree_fft/2+1) * 2 * sizeof(double));
				for(size_t i = 0 ; i < dim ; ++i){
					for (size_t j = 0 ; j < Ksize ; ++j){
						fft::product<degree_fft>(fft_cache3, val_fft[i * Ksize + j], rhs(i, j).get_elem().get_msb_fft());
						for(size_t k = 0 ; k < degree_fft/2+1 ; ++k){
							fft_cache1[k][0] += fft_cache3[k][0];
							fft_cache1[k][1] += fft_cache3[k][1];
						}
						fft::product<degree_fft>(fft_cache3, val_fft[i * Ksize + j], rhs(i, j).get_elem().get_lsb_fft());
						for(size_t k = 0 ; k < degree_fft/2+1 ; ++k){
							fft_cache2[k][0] += fft_cache3[k][0];
							fft_cache2[k][1] += fft_cache3[k][1];
						}
					}
				}
				res[dimP].compute_fft_inv(fft_cache1);
				res[dimP].shift(thresh);

				res_lsb.compute_fft_inv(fft_cache2);
				res[dimP] += res_lsb;
			}
			for(size_t j = 0 ; j < dim * Ksize ; ++j)
				fftw_free(val_fft[j]);
			return res;
		}
};

template<class Int, uint64_t degree, uint64_t degree_fft>
FFT<degree_fft> CirculantRing<Int, degree, degree_fft>::fft;
template<class Int, uint64_t degree, uint64_t degree_fft>
Int CirculantRing<Int, degree, degree_fft>::tmp[degree];
template<class Int, uint64_t degree, uint64_t degree_fft>
fftw_complex CirculantRing<Int, degree, degree_fft>::fft_cache1[degree_fft/2+1] __attribute__ ((aligned (32)));
template<class Int, uint64_t degree, uint64_t degree_fft>
fftw_complex CirculantRing<Int, degree, degree_fft>::fft_cache2[degree_fft/2+1] __attribute__ ((aligned (32)));
template<class Int, uint64_t degree, uint64_t degree_fft>
fftw_complex CirculantRing<Int, degree, degree_fft>::fft_cache3[degree_fft/2+1] __attribute__ ((aligned (32)));

typedef CirculantRing<Zqp, 1, 1> Rz;
typedef CirculantRing<Zq, P1> Rp1;
typedef CirculantRing<Zq, P2> Rp2;
typedef CirculantRing<Zqp, P1> Rp1_p;
typedef CirculantRing<Zqp, P2> Rp2_p;
typedef CirculantRing<Zqcrt, P1> Rp1_crt;
typedef CirculantRing<Zqcrt, P2> Rp2_crt;
typedef CirculantRing<Zqp, P1*P2, FFT_DIM2> Rp12;
typedef CirculantRing<Zqcrt, P1*P2, FFT_DIM2> Rp12_crt;

#endif
