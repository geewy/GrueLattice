/**
 * @file circulant_ring.cpp
 * @date may 2018
 * @author geewy
 * @brief Declaration of the CirculantRing class
 * 
 * See circulant_ring.hpp for more details
 */

#include <algorithm>
#include <cmath>
#include <random>
#include <limits>
#include "circulant_ring.hpp"
#include "operations.hpp"

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft>::CirculantRing() : msb_fft(nullptr), lsb_fft(nullptr), fft_status(false){
	elem = new Int[degree];
}

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft>::CirculantRing(const CirculantRing &src){
	fft_status = src.fft_status;
	if(fft_status){
		elem = nullptr;
		msb_fft = fftw_alloc_complex(degree_fft/2+1);
		lsb_fft = fftw_alloc_complex(degree_fft/2+1);
		if(degree < 2000){
			std::memcpy(msb_fft, src.msb_fft, (degree_fft/2+1) * 2 * sizeof(double));
			std::memcpy(lsb_fft, src.lsb_fft, (degree_fft/2+1) * 2 * sizeof(double));
		}
		else{
			for(size_t i = 0; i < degree_fft/2+1; ++i){
				msb_fft[i][0] = src.msb_fft[i][0];
                msb_fft[i][1] = src.msb_fft[i][1];
                lsb_fft[i][0] = src.lsb_fft[i][0];
				lsb_fft[i][1] = src.lsb_fft[i][1];
			}
		}
	}
	else{
		elem = new Int[degree];
		msb_fft = nullptr;
		lsb_fft = nullptr;
		std::copy(src.elem, src.elem + degree, elem);
	}
}

template<class Int, uint64_t degree, uint64_t degree_fft>
template<class Int2>
CirculantRing<Int, degree, degree_fft>::CirculantRing(const CirculantRing<Int2, degree, degree_fft> &src){
	assert(! src.get_fft_status());
	elem = new Int[degree];
	Int2 *src_elem = src.get_data();
	for(size_t i = 0; i < degree; ++i)
		elem[i] = src_elem[i];
	fft_status = false;
	msb_fft = nullptr;
	lsb_fft = nullptr;	
}

template<class Int, uint64_t degree, uint64_t degree_fft>
template<uint64_t degree_fft2>
CirculantRing<Int, degree, degree_fft>::CirculantRing(const CirculantRing<Int, degree, degree_fft2> &src){
	assert(! src.get_fft_status());
	const Int *src_elem = src.get_data();
	elem = new Int[degree];
	std::copy(src_elem, src_elem + degree, elem);
	fft_status = false;
	msb_fft = nullptr;
	lsb_fft = nullptr;
}

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft>::CirculantRing(Int *values, bool copy){
	if(copy){
		elem = new Int[degree];
		std::copy(values, values + degree, elem);	
	}
	else
		elem = std::move(values);
	
	fft_status = false;
	msb_fft = nullptr;
	lsb_fft = nullptr;
}

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft>::CirculantRing (fftw_complex *msb_fft_values, fftw_complex *lsb_fft_values, bool copy){
    if(copy){
        msb_fft = fftw_alloc_complex(degree_fft/2+1);
        lsb_fft = fftw_alloc_complex(degree_fft/2+1);
        std::memcpy(msb_fft, msb_fft_values, (degree_fft/2+1) * 2 * sizeof(double));
        std::memcpy(lsb_fft, lsb_fft_values, (degree_fft/2+1) * 2 * sizeof(double));
    }
    else{
        msb_fft = std::move(msb_fft_values);
        lsb_fft = std::move(lsb_fft_values);
    }
    fft_status = true;
    elem = nullptr;
}

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft>::CirculantRing (const Int val){
    elem = new Int[degree];
    elem[0] = val;
    std::fill(elem + 1, elem + degree, 0);
    fft_status = false;
    msb_fft = nullptr;
    lsb_fft = nullptr;
}
template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft>::CirculantRing (const int64_t *values){
    elem = new Int[degree];
    for(size_t i = 0 ; i < degree ; ++i)
        elem[i] = values[i];
        
    fft_status = false;
    msb_fft = nullptr;
    lsb_fft = nullptr;
}
template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft>::CirculantRing (const Eva &E){
    elem = new Int[degree];
    elem[0] = E(0);
    for(size_t i = 1 ; i < degree ; ++i)
        elem[degree - i] = E(i);
        
    fft_status = false;
    msb_fft = nullptr;
    lsb_fft = nullptr;
}
template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft>::~CirculantRing(){
    if(fft_status){
        fftw_free(msb_fft);
        fftw_free(lsb_fft);
        msb_fft = nullptr;
        lsb_fft = nullptr;
    }
    else{
        delete[] elem;
        elem = nullptr;
    }
}

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft> &CirculantRing<Int, degree, degree_fft>::operator=(const int &other){
    assert(! fft_status);
    elem[0] = Int(other);
    std::fill(elem + 1, elem + degree, 0);
    fft_status = false;
    msb_fft = nullptr;
    lsb_fft = nullptr;
    return *this;
}

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft> &CirculantRing<Int, degree, degree_fft>::operator=(const CirculantRing &other){
    if(this != &other){
        if(other.fft_status){
            if(! fft_status){
                delete[] elem;
                msb_fft = fftw_alloc_complex(degree_fft/2+1);
                lsb_fft = fftw_alloc_complex(degree_fft/2+1);
                fft_status = true;
            }
            std::memcpy(msb_fft, other.msb_fft, (degree_fft/2+1) * 2 * sizeof(double));
            std::memcpy(lsb_fft, other.lsb_fft, (degree_fft/2+1) * 2 * sizeof(double));
        }
        else{
            if(fft_status){
                fftw_free(msb_fft);
                fftw_free(lsb_fft);
                elem = new Int[degree];
                fft_status = false;
            }
            std::copy(other.elem, other.elem + degree, elem);
        }
    }
    return *this;
}
template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft> &CirculantRing<Int, degree, degree_fft>::operator=(CirculantRing &&other) noexcept{
    if(this != &other){
        if(fft_status){
            fftw_free(msb_fft);
            fftw_free(lsb_fft);
        }
        else
            delete[] elem;

        if(other.fft_status){
            fftw_complex *old_msb = std::move(other.msb_fft);
            other.msb_fft = std::forward<fftw_complex *>(nullptr);
            msb_fft = old_msb;
            fftw_complex *old_lsb = std::move(other.lsb_fft);
            other.lsb_fft = std::forward<fftw_complex *>(nullptr);
            lsb_fft = old_lsb;
        }
        else{
            Int *old_value = std::move(other.elem);
            other.elem = std::forward<Int *>(nullptr);
            elem = old_value;
        }
        fft_status = other.fft_status;
    }
    return *this;
}

template<class Int, uint64_t degree, uint64_t degree_fft>
void CirculantRing<Int, degree, degree_fft>::mult(const fftw_complex *f){
    const int64_t thresh = Int::get_split_threshold();
    CirculantRing lsb, msb;
    split(msb, lsb);

    msb.compute_fft(fft_cache1);
    fft::product<degree_fft>(fft_cache2, fft_cache1, f);

    compute_fft_inv(fft_cache2);
    for(size_t i = 0 ; i < degree ; ++i)
        elem[i] = ((int64_t) elem[i] % thresh) * thresh;

    lsb.compute_fft(fft_cache1);
    fft::product<degree_fft>(fft_cache2, fft_cache1, f);
    compute_fft_inv(fft_cache2, tmp);
    for(size_t i = 0 ; i < degree ; ++i)
        elem[i] += tmp[i];
}

template<class Int, uint64_t degree, uint64_t degree_fft>
void CirculantRing<Int, degree, degree_fft>::mult(const CirculantRing &rhs){
    rhs.compute_fft(fft_cache3);

    const int64_t spliter = Int::get_split_threshold();
    CirculantRing lsb, msb;
    split(msb, lsb);

    msb.compute_fft(fft_cache1);
    fft::product<degree_fft>(fft_cache2, fft_cache1, fft_cache3);

    compute_fft_inv(fft_cache2);
    shift(spliter);

    lsb.compute_fft(fft_cache1);
    fft::product<degree_fft>(fft_cache2, fft_cache1, fft_cache3);
    compute_fft_inv(fft_cache2, tmp);
    for(size_t i = 0 ; i < degree ; ++i)
        elem[i] += tmp[i];
}

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft> &CirculantRing<Int, degree, degree_fft>::operator*=(const CirculantRing &rhs){
    assert (rhs.fft_status);
    const int64_t thresh = Int::get_split_threshold();
    compute_fft(fft_cache1);

    fft::product<degree_fft>(fft_cache2, fft_cache1, rhs.get_msb_fft());

    compute_fft_inv(fft_cache2);
    shift(thresh);

    fft::product<degree_fft>(fft_cache2, fft_cache1, rhs.get_lsb_fft());
    compute_fft_inv(fft_cache2, tmp);
    for(size_t i = 0 ; i < degree ; ++i)
        elem[i] += tmp[i];

    return *this;
}

template<class Int, uint64_t degree, uint64_t degree_fft>
void CirculantRing<Int, degree, degree_fft>::compute_fft(fftw_complex *elem_fft) const{
    double *in = fft.get_in();
    for(size_t i = 0 ; i < degree ; ++i)
        in[i] = (double) elem[i];
    std::fill(in + degree, in + degree_fft, 0.0);
    fft.forward(elem_fft);
}

template<class Int, uint64_t degree, uint64_t degree_fft>
void CirculantRing<Int, degree, degree_fft>::compute_fft_inv(fftw_complex *elem_fft, Int *values) const{
    double *in = fft.get_in();
    fft.backward(elem_fft);

    if(degree > 1)
        for (size_t i = 0 ; i < degree ; ++i)
			values[i] = (int64_t)std::floor(in[i] + 0.5) + (int64_t)std::floor(in[i+degree]+0.5);   
 
    else
		values[0] = (int64_t)std::floor(in[0] + 0.5);
}

template<class Int, uint64_t degree, uint64_t degree_fft>
void CirculantRing<Int, degree, degree_fft>::decomp_fft(){
    assert(! fft_status);
    CirculantRing msb, lsb;
    split(msb, lsb);

    msb_fft = fftw_alloc_complex(degree_fft/2+1);
    lsb_fft = fftw_alloc_complex(degree_fft/2+1);
    msb.compute_fft(msb_fft);
    lsb.compute_fft(lsb_fft);
    delete[] elem;
    fft_status = true;
}

template<class Int, uint64_t degree, uint64_t degree_fft>
void CirculantRing<Int, degree, degree_fft>::recomp_fft(const int64_t cons){
    assert(fft_status);
    Int *values = new Int[degree];
    compute_fft_inv(msb_fft, values);
    elem = new Int[degree];
    for(size_t i = 0 ; i < degree ; ++i)
        elem[i] = (values[i] % cons) * cons;

    compute_fft_inv(lsb_fft, values);
    for(size_t i = 0 ; i < degree ; ++i)
        elem[i] += values[i];
    
    fftw_free(msb_fft);
    fftw_free(lsb_fft);
    delete[] values;
    fft_status = false;
}

template<class Int, uint64_t degree, uint64_t degree_fft>
int64_t CirculantRing<Int, degree, degree_fft>::split(CirculantRing &msb, CirculantRing &lsb) const{
    Int *most_significant = new Int[degree];
    Int *least_significant = new Int[degree];
    int64_t used_threshold = elem[0].split(most_significant[0], least_significant[0]);
    
    for(size_t i = 1 ; i < degree ; ++i)
        elem[i].split(most_significant[i], least_significant[i]);
        
    msb = CirculantRing(most_significant, false);
    lsb = CirculantRing(least_significant, false);

    return used_threshold;
}

template<class Int, uint64_t degree, uint64_t degree_fft>
void CirculantRing<Int, degree, degree_fft>::unsplit(CirculantRing &msb, CirculantRing &lsb){
    for(size_t i = 0 ; i < degree ; ++i)
        elem[i].unsplit(msb.elem[i], lsb.elem[i]);
}

template<class Int, uint64_t degree, uint64_t degree_fft>
template<class Ring2>
Ring2 CirculantRing<Int, degree, degree_fft>::rounding(const uint64_t old_modulus, const uint64_t new_modulus) const{
    int64_t *res = new int64_t[degree];
    for(size_t i = 0 ; i < degree ; ++i)
        res[i] = (int64_t)std::floor((double) elem[i] * (double) new_modulus / (double) old_modulus + 0.5);
        
    Ring2 result(res);
    delete[] res;
    return result;
}

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft> CirculantRing<Int, degree, degree_fft>::exact_rounding(const uint64_t old_modulus, const uint64_t new_modulus) const{
    Int *res = new Int[degree];
    for(size_t i = 0 ; i < degree ; ++i)
        res[i] = (int64_t)std::floor((double) elem[i] * (double) new_modulus / (double) old_modulus + 0.5);
    return CirculantRing<Int, degree, degree_fft>(res, false);
}


template<class Int>
static Int sample_long(const double variance){
    double rander = 0;
    for(size_t i = 0 ; i < 12 ; ++i)
        rander += static_cast<double> (rand()) / static_cast<double> (RAND_MAX) - 0.5;
    return (int64_t)floor(rander * std::sqrt(variance) + 0.5);
}

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft> CirculantRing<Int, degree, degree_fft>::gaussian_sample(const double variance){
    Int *res = new Int[degree];
    for(size_t i = 0 ; i < degree ; ++i)
        res[i] = sample_long<Int>(variance);
    return CirculantRing<Int, degree, degree_fft>(res, false);
}

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft> CirculantRing<Int, degree, degree_fft>::uniform_sample(){
    Int *res = new Int[degree];
    for(size_t i = 0 ; i < degree ; ++i)
        res[i] = (((int64_t)rand())<<31) + (int64_t)rand();
    return CirculantRing<Int, degree, degree_fft>(res, false);
}

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft> CirculantRing<Int, degree, degree_fft>::binary_sample(){
    Int *res = new Int[degree];
    for(size_t i = 0 ; i < degree ; ++i)
        res[i] = rand() % 2;
    return CirculantRing<Int, degree, degree_fft>(res, false);
}

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft> CirculantRing<Int, degree, degree_fft>::sample_a(){
    Int *res = new Int[degree];
    if(degree == 1)
        res[0] = (((int64_t)rand())<<31) + (int64_t)rand();

    else{
        Int sum;
        sum = 0;
        for(size_t i = 0 ; i < degree-1 ; ++i){
            res[i] = (((int64_t)rand())<<31) + (int64_t)rand();
            sum += res[i];
        }
        res[degree-1] = -sum;
    }
    return CirculantRing<Int, degree, degree_fft>(res, false);
}

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft> CirculantRing<Int, degree, degree_fft>::sample_s(const double density){
    Int *res = new Int[degree];
    if(degree == 1)
        res[0] = rand() % 3 - 1;
    
    else{
        std::fill(res, res + degree, 0);

        size_t position, defined;
        defined = 0;
        while(defined < degree*density){
            position = rand() % degree;
            if(res[position] == 0){
                res[position] = 1;
                ++defined;
            }
        }
        defined = 0;
        while(defined < degree*density){
            position = rand() % degree;
            if(res[position] == 0){
                res[position] = -1;
                ++defined;
            }
        }
    }
    return CirculantRing<Int, degree, degree_fft>(res, false);
}

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft> CirculantRing<Int, degree, degree_fft>::sample_e(const double variance){
    if(variance > 1e6)
        return gaussian_sample(variance);
    else{
        Int *res = new Int[degree];
        if(degree == 1)
            res[0] = rand() % 3 - 1;
        
        else{
            std::fill(res, res + degree, 0);

            for(size_t i = 0 ; i < variance*degree/2 ; ++i){
                res[rand() % degree] += 1;
                res[rand() % degree] -= 1;
            }
        }
        return CirculantRing<Int, degree, degree_fft>(res, false);
    }
}

template<class Int, uint64_t degree, uint64_t degree_fft>
CirculantRing<Int, degree, degree_fft> CirculantRing<Int, degree, degree_fft>::galois(const uint64_t alpha) const{
    CirculantRing<Int, degree, degree_fft> result(*this);
    result.galois_inplace(alpha);
    return result;
}

template<class Int, uint64_t degree, uint64_t degree_fft>
void CirculantRing<Int, degree, degree_fft>::galois_inplace(const uint64_t alpha){
    Int *tmp = new Int[degree];
    size_t current = 0;
    for(size_t i = 0 ; i < degree ; ++i){
        tmp[current] = elem[i];
        current += alpha;
        current %= degree;
    }
    std::copy(tmp, tmp + degree, elem);
    delete[] tmp;
}

template<class Int, uint64_t degree, uint64_t degree_fft>
template<uint64_t degree2, uint64_t degree2_fft, uint64_t degree3_fft>
CirculantRing<Int, degree*degree2, degree3_fft> CirculantRing<Int, degree, degree_fft>::tensor_product(const CirculantRing<Int, degree2, degree2_fft>& rhs) const{
    const uint64_t new_dim = degree*degree2;
    Int *result = new Int[degree*degree2];
    Int *rhs_x = rhs.get_data();
    size_t current = 0;
    size_t current_i = 0;
    for(size_t i = 0 ; i < degree ; ++i){
        current = current_i;
        for (size_t j = 0 ; j < degree2 ; ++j){
            result[current] = elem[i] * rhs_x[j];
            current += degree;
            current %= new_dim;
        }
        current_i += degree2;
        current_i %= new_dim;
    }
    return CirculantRing<Int, degree*degree2, degree3_fft>(result, false);
}

template<class Int, uint64_t pq, uint64_t degree_fft>
template <uint64_t degreeP, uint64_t degree2_fft>
CirculantRing<Int, degreeP, degree2_fft> CirculantRing<Int, pq, degree_fft>::trace() const{
    const uint64_t mod = pq/degreeP;
    Int *elem_res = new Int[degreeP];

    size_t idx = 0;
    for(size_t i = 0 ; i < degreeP ; ++i){
        elem_res[i] = elem[idx];
        idx += mod;
    }
    return CirculantRing<Int, degreeP, degree2_fft>(elem_res, false);
}

template class CirculantRing<Zt, 1, 1>;
template class CirculantRing<Zt, P1>;
template class CirculantRing<Zt, P2>;
template class CirculantRing<Zqp, 1, 1>;
template class CirculantRing<Zp12, 1, 1>;
template class CirculantRing<Zq, P1>;
template class CirculantRing<Zq, P2>;
template class CirculantRing<Zqcrt, P1>;
template class CirculantRing<Zqcrt, P2>;
template class CirculantRing<Zqcrt, P1*P2, FFT_DIM2>;
template class CirculantRing<Zqp, P1*P2, FFT_DIM2>;

template Rp1::CirculantRing<Zt>(const CirculantRing<Zt, P1> &src);
template Rp2::CirculantRing<Zt>(const CirculantRing<Zt, P2> &src);

template CirculantRing<Zt, P1>::CirculantRing<Zq>(const Rp1 &src);
template CirculantRing<Zt, P2>::CirculantRing<Zq>(const Rp2 &src);
template CirculantRing<Zt, 1, 1>::CirculantRing<Zqp>(const Rz &src);

template Rp12_crt Rp1_crt::tensor_product<P2, FFT_DIM1, FFT_DIM2>(const Rp2_crt &rhs) const;
template Rp12 Rp1_p::tensor_product<P2, FFT_DIM1, FFT_DIM2>(const Rp2_p &rhs) const;

template CirculantRing<Zp12, 1, 1> Rz::rounding<CirculantRing<Zp12, 1, 1> >(const uint64_t old_modulus, const uint64_t new_modulus) const;
template Rp1_crt Rp1::rounding<Rp1_crt>(const uint64_t old_modulus, const uint64_t new_modulus) const;
template Rp2_crt Rp2::rounding<Rp2_crt>(const uint64_t old_modulus, const uint64_t new_modulus) const;
template Rp12 Rp12_crt::rounding<Rp12>(const uint64_t old_modulus, const uint64_t new_modulus) const;

template Rp1_p Rp12::trace<P1, FFT_DIM1>() const;

/** For tests only **/
template class CirculantRing<Zt, P1*P2, FFT_DIM2>;
template CirculantRing<Zt, 1, 1>::CirculantRing<Zp12>(const CirculantRing<Zp12, 1, 1> &src);
template CirculantRing<Zt, P1*P2, FFT_DIM2>::CirculantRing<Zqp>(const Rp12 &src);
template CirculantRing<Zp12, 1, 1>::CirculantRing<Zqp>(const Rz &src);

/** For stats only **/
template Rp1_crt::CirculantRing<Zq>(const Rp1 &src);
template Rp12_crt::CirculantRing<Zqp>(const Rp12 &src);
template Rp1::CirculantRing<Zqcrt>(const Rp1_crt &src);
template Rp2::CirculantRing<Zqcrt>(const Rp2_crt &src);
template Rp12::CirculantRing<Zqcrt>(const Rp12_crt &src);
