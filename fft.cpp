/**
 * @file fft.cpp
 * @date may 2018
 * @author geewy
 * @brief Declaration of FFT utility
 * 
 * See fft.hpp for more details.
 */

#include <iostream>
#include <string>
#include <cstring>
#include "fft.hpp"
#include "parameters.hpp"

template<size_t n>
FFT<n>::FFT(){
	in = fftw_alloc_real(n);
	out = fftw_alloc_complex(n/2+1);
	
	std::string file_name = "wisdom" + std::to_string(n);
	int existing_wisdom = fftw_import_wisdom_from_filename(file_name.c_str());
	
	if(existing_wisdom == 0)
		std::cerr << "Generating FFTW plan..." << std::flush;
	
	forward_plan = fftw_plan_dft_r2c_1d(n, in, out, PLAN_MODE);
	backward_plan = fftw_plan_dft_c2r_1d(n, out, in, PLAN_MODE);
	
	if(existing_wisdom == 0){
		std::cerr << "Done." << std::endl;
		fftw_export_wisdom_to_filename(file_name.c_str());
	}
}

template<size_t n>
FFT<n>::~FFT(){
	fftw_destroy_plan(forward_plan);
	fftw_destroy_plan(backward_plan);
	fftw_free(in);
	fftw_free(out);
}

template<size_t n>
void FFT<n>::forward(fftw_complex *res){
	fftw_execute(forward_plan);
	std::memcpy(res, out, (n/2+1) * 2 * sizeof(double));
}

template<size_t n>
void FFT<n>::backward(const fftw_complex *data){
	for(size_t i = 0; i < n/2+1; ++i){
		out[i][0] = data[i][0] / (n);
		out[i][1] = data[i][1] / (n);
	}
	fftw_execute(backward_plan);	
}

namespace fft{
	template<size_t n>
	void product(fftw_complex *res_fft, const fftw_complex *lhs_fft, const fftw_complex *rhs_fft){
		for(size_t i = 0; i < n/2+1; ++i){
			const double lhs_r = lhs_fft[i][0];
			const double lhs_i = lhs_fft[i][1];
			const double rhs_r = rhs_fft[i][0];
			const double rhs_i = rhs_fft[i][1];
			const double x0y0 = lhs_r * rhs_r;
			const double x1y1 = lhs_i * rhs_i;
			const double x0y1 = lhs_r * rhs_i;
			const double x1y0 = lhs_i * rhs_r;
			res_fft[i][0] = x0y0 - x1y1;
			res_fft[i][1] = x0y1 + x1y0;
		}
	}
	
	template<size_t p, size_t q>
	void tensor_product(fftw_complex *res_fft, const fftw_complex * lhs_fft, const fftw_complex* rhs_fft){
		fftw_complex *lhs_big = fftw_alloc_complex(p*q/2+1);
		fftw_complex *rhs_big = fftw_alloc_complex(p*q/2+1);
		
		for(size_t i = 0; i < p/2+1; ++i){
			for(size_t j = 0; j < q/2; ++j){
				lhs_big[j * p + i][0] = lhs_fft[i][0];
				lhs_big[j * p + i][1] = lhs_fft[i][1];
			}
		}	
		for(size_t i = 1; i < p/2+1; ++i){
			for(size_t j = 0; j < q/2; ++j){
				lhs_big[j * p + p - i][0] = lhs_fft[i][0];
				lhs_big[j * p + p - i][1] = lhs_fft[i][1];
			}
		}
		for(size_t j = 0; j < q/2+1; ++j){
			for(size_t i = 0; i <= p/2; ++i){
				rhs_big[i * q + j][0] = rhs_fft[j][0];
				rhs_big[i * q + j][1] = rhs_fft[j][1];
			}
		}
		for(size_t j = 1; j < q/2+1; ++j){
			for(size_t i = 0; i <= p/2; ++i){
				rhs_big[i * q + q - j][0] = rhs_fft[j][0];
				rhs_big[i * q + q - j][1] = rhs_fft[j][1];
			}
		}
	}
}

template class FFT<1>;
template class FFT<FFT_DIM1>;
template class FFT<FFT_DIM2>;

template void fft::product<1>(fftw_complex *res_fft, const fftw_complex * lhs_fft, const fftw_complex* rhs_fft);
template void fft::product<FFT_DIM1>(fftw_complex *res_fft, const fftw_complex * lhs_fft, const fftw_complex* rhs_fft);
template void fft::product<FFT_DIM2>(fftw_complex *res_fft, const fftw_complex * lhs_fft, const fftw_complex* rhs_fft);

template void fft::tensor_product<P1, P2>(fftw_complex *res_fft, const fftw_complex * lhs_fft, const fftw_complex* rhs_fft);
