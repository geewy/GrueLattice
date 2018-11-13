#include <chrono>
#include <cmath>
#include <ctime>
#include "circulant_ring.hpp"
#include "operations.hpp"
#include "rgsw.hpp"

void print_param(){
	std::cout << "T = " << T << ",";
	std::cout << "Q = " << Q << ",";
	std::cout << "Qp = " << Qp << ",";
	std::cout << "P1 = " << P1 << ",";
	std::cout << "P2 = " << P2 << ",";
	std::cout << "N = " << N << ",";
	std::cout << "B_Q = " << B_Q << ",";
	std::cout << "K_Q = " << K_Q << ",";
	std::cout << "B_Qp = " << B_Qp << ",";
	std::cout << "K_Qp = " << K_Qp << std::endl;
}

template<class Ring_op, uint64_t size>
Ring_op dot_product(const Ring_op *const vec1, const Ring_op *const vec2){
	Ring_op res;
	res = 0;
	Ring_op *vec2_rw = new Ring_op[size];
	for(size_t i = 0; i < size; ++i){
		vec2_rw[i] = vec2[i];
		vec2_rw[i].decomp_fft();
	}
	for(size_t i = 0; i < size; ++i)
		res += vec1[i] * vec2_rw[i];
		
	delete[] vec2_rw;
	return res;
}

void exp_crt(RLWE<Rp12_crt, 3> &cpq, const RLWE<Rp1_crt, 1> &cp, const RLWE<Rp2_crt, 1> &cq){
	const RLWE<Rp1_crt, 1> c_p = cp.galois(P2inv_modP1);
    const Rp1_crt array_p = c_p.get_array()[0];
	const Rp1_crt elem_p = c_p.get_elem();
	
	const RLWE<Rp2_crt, 1> c_q = cq.galois(P1inv_modP2);
    const Rp2_crt array_q = c_q.get_array()[0];
	const Rp2_crt elem_q = c_q.get_elem();
	
	Rp12_crt array_crt[3];
	array_crt[0] = QpoverT_inv * array_p.tensor_product<P2, FFT_DIM1, FFT_DIM2>(array_q);
	array_crt[1] = QpoverT_inv * array_p.tensor_product<P2, FFT_DIM1, FFT_DIM2>(elem_q);
	array_crt[2] = QpoverT_inv * elem_p.tensor_product<P2, FFT_DIM1, FFT_DIM2>(array_q);
	Rp12_crt elem_crt = QpoverT_inv * elem_p.tensor_product<P2, FFT_DIM1, FFT_DIM2>(elem_q);
	
	cpq = RLWE<Rp12_crt, 3>(array_crt, elem_crt);
}

void crt_key(Rp12 *spq, const Rp1 &sp, const Rp2 &sq){
	const Rp1_p s_p = sp.galois(P2inv_modP1);
	const Rp2_p s_q = sq.galois(P1inv_modP2);
	
	spq[0] = - s_p.tensor_product<P2, FFT_DIM1, FFT_DIM2>(s_q);
	spq[1] = s_p.tensor_product<P2, FFT_DIM1, FFT_DIM2>(Rp2_p(1));
	spq[2] = Rp1_p(1).tensor_product<P2, FFT_DIM1, FFT_DIM2>(s_q);
}

template<class Ring_op2, class Yelem, uint64_t inMod>
RLWE1<Ring_op2> ext_exp_inner(const uint64_t encMod, const uint64_t plainMod, const size_t length, const Yelem *plainElem, const RGSW<Ring_op2> *rgswEnc, const KSWKey<Ring_op2, 1, 1, B_Q, K_Q> *kswArray){
	Ring_op2 zeros[1];
	zeros[0] = 0;
	Ring_op2 bit((int64_t) (encMod / plainMod));
	RLWE1<Ring_op2> result(zeros, bit);
	
	Yelem current;
	size_t one[length], other[length], count_one = 0, count_other = 0;
	for(size_t i = 0; i < length; ++i){
		if(plainElem[i] == 0) /*! < nothing to add to ext_exp_inner */
			continue;
		if(plainElem[i] == 1) /*! < no galois + ksw requires */
			one[count_one++] = i;
		else
			other[count_other++] = i;
	}
	
	if(count_other > 0){
		for(size_t i = 0; i < count_other-1; ++i){
			current = plainElem[other[i]] * plainElem[other[i+1]].inv();
			result.ext_mult_inplace(rgswEnc[other[i]]);
			
			if(current != 1){
				result.galois_inplace((size_t) current);
				result.key_switch_inplace(kswArray[other[i]]);
			}
		}
		current = plainElem[other[count_other - 1]];
		result.ext_mult_inplace(rgswEnc[other[count_other - 1]]);
		result.galois_inplace((size_t) current);
		result.key_switch_inplace(kswArray[(size_t) current]);
	}
	
	for(size_t i = 0; i < count_one; ++i)
		result.ext_mult_inplace(rgswEnc[one[i]]);
	
	return result;		
}

LWE fun_exp_extract(const fftw_complex *func, const Rp12LWE &c_Zm, const KSWKeyRp12 &KSW){
	std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> ksw_time, mult_time, trace_time;

    start = std::chrono::system_clock::now();
    RLWE<Rp12 , 1> c_1 = c_Zm.key_switch(KSW);
    end = std::chrono::system_clock::now();
    ksw_time = end - start;

    start = std::chrono::system_clock::now();
    c_1.mult(func);
    end = std::chrono::system_clock::now();
    mult_time = end - start;

    start = std::chrono::system_clock::now();
    RLWE<Rp1_p, 1> c_2 = c_1.template trace<Rp1_p, P1, FFT_DIM1>();
    Zqp *array = c_2.get_array()->get_data();
    Rz cipher[P1 + 1];
    cipher[0] = array[0];
    
    for (size_t i = 1 ; i < P1 ; ++i)
        cipher[P1 - i] = array[i];
    
    cipher[P1] = c_2.get_elem().get_data()[0];
    end = std::chrono::system_clock::now();
    trace_time = end - start;
    std::cerr << "In N_f -> KeySwitch: " << ksw_time.count() << "s, Mult: " << mult_time.count() << "s, Tr*:" << trace_time.count() << "s" << std::endl;
	
	return LWE(cipher);
}

template<class Ring_op, class Elem, class R_E, uint64_t mod>
RLWE<R_E, 1> accumulation(const RLWE<CirculantRing<Zp12, 1, 1>, N> &enc, const RGSW<Ring_op> BK[N], const KSWKey<Ring_op> KS[mod]){
	Elem *cipher_p = new Elem[N+1];
	enc.mod<Elem>(cipher_p);
	for(size_t i = 0; i < N; ++i)
		cipher_p[i] = - cipher_p[i];
		
	Zq *Tb_coefs = new Zq[mod];
	std::fill(Tb_coefs, Tb_coefs + mod, 0);
	Tb_coefs[(size_t) cipher_p[N]] = 1;
	Ring_op Tn(Tb_coefs, false);
	
	RLWE<Ring_op, 1> acc_p = ext_exp_inner<Ring_op, Elem, mod>(Q, T, N, cipher_p, BK, KS);
	acc_p.mult(Tn);
	return acc_p.template mod_switch<R_E, Qcrt, Q>();
}

Rp12LWE preparation(const LWE &enc, const KSWKeyLWE &S_LWE, const Rp1GSW Xsi[N], const KSWKeyRp1 KSp1[P1], const Rp2GSW Ysi[N], const KSWKeyRp2 KSp2[P2]){
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> ksw_time, msw_time, acc1_time, acc2_time, crt_time;
	
	start = std::chrono::system_clock::now();
	RLWE<Rz, N> enc_short = enc.key_switch(S_LWE);
	end = std::chrono::system_clock::now();
	ksw_time = end - start;
	
	start = std::chrono::system_clock::now();
    RLWE<CirculantRing<Zp12, 1, 1>, N> cipher_short = enc_short.mod_switch<CirculantRing<Zp12, 1, 1>, P1*P2, Qp>();
    end = std::chrono::system_clock::now();
    msw_time = end - start;

    start = std::chrono::system_clock::now();
    RLWE<Rp1_crt, 1> crt_in_p = accumulation<Rp1, Zp1, Rp1_crt, P1>(cipher_short, Xsi, KSp1);
    end = std::chrono::system_clock::now();
    acc1_time = end - start;

    start = std::chrono::system_clock::now();
    RLWE<Rp2_crt, 1> crt_in_q = accumulation<Rp2, Zp2, Rp2_crt, P2>(cipher_short, Ysi, KSp2);
    end = std::chrono::system_clock::now();
    acc2_time = end - start;

    RLWE<Rp12_crt, 3> result;
    start = std::chrono::system_clock::now();
    exp_crt(result, crt_in_p, crt_in_q);
    end = std::chrono::system_clock::now();
	crt_time = end - start;
	
	std::cerr << "In L_c -> KeySwitchLWE: " << ksw_time.count() << "s, ModSwitch: " << msw_time.count() << "s, ExtExpInner 1: " << acc1_time.count() << "s, ExtExpInner 2: " << acc2_time.count() << "s, ExpCRT: " << crt_time.count() << "s" << std::endl;
	return result.template mod_switch<Rp12, Qp, Qcrt>();
}

LWE combination(const size_t num, const int64_t *coefs, const LWE *c_i){
	LWE result = c_i[0] * coefs[0];
	for(size_t i = 0; i < num; ++i)
		result += c_i[i] * coefs[i];
	return result;
}

LWE gate(const size_t num, const int64_t *coefs, const LWE *c_i, const fftw_complex *func, const KSWKeyLWE &S_LWE, const Rp1GSW Xsi[N],  const KSWKeyRp1 KSp1[P1], const Rp2GSW Ysi[N], const KSWKeyRp2 KSp2[P2], const KSWKeyRp12 &KSW){
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> comb_time, prep_time, func_time;
	
	
    start = std::chrono::system_clock::now();
    LWE enc = combination(num, coefs, c_i);
    end = std::chrono::system_clock::now();
    comb_time = end - start;

    start = std::chrono::system_clock::now();
    Rp12LWE cipher_Zm = preparation(enc, S_LWE, Xsi, KSp1, Ysi, KSp2);
    end = std::chrono::system_clock::now();
    prep_time = end - start;

    start = std::chrono::system_clock::now();
    LWE Gate = fun_exp_extract(func, cipher_Zm, KSW);
    end = std::chrono::system_clock::now();
    func_time = end - start;

    std::cerr << "Combination: " << comb_time.count() << "s, L_c: " << prep_time.count() << "s, N_f: " << func_time.count() << "s" << std::endl;
	return Gate;
}

template<class Ring_op, uint64_t size, uint64_t dim_op>
void gen_bootstrapping_keys(const uint64_t mod, const uint64_t plainMod, RGSW<Ring_op> Tsi[size], const Rz key[size], const Ring_op &key_p, const double variance){
	Zt *coefs = new Zt[dim_op];
    for (size_t j = 0 ; j < dim_op ; ++j)
        coefs[j] = 0;
    for (size_t i = 0 ; i < size ; ++i){
        const size_t position = (dim_op + (int64_t)(key[i].get_data()[0])) % dim_op;
        coefs[position] = 1;
        CirculantRing<Zt, dim_op> T_si(coefs); 
        Tsi[i].encrypt(&key_p, T_si, variance);
        coefs[position] = 0;
    }
	delete[] coefs;
}

template<class Ring_op, uint64_t dim_op>
void gen_keyswitching_keys(KSWKey<Ring_op, 1, 1, B_Q, K_Q> keys[dim_op], const Ring_op &key_p, const double variance){
	Ring_op key_alpha[1];
    for (size_t i = 1 ; i < dim_op ; ++i){
        key_alpha[0] = key_p.galois(i);
        keys[i].init(key_alpha, &key_p, variance);
	}
}

void gen_funexpextract_key(KSWKeyRp12 *KSW, const Rp12 key_pq[3], const Rz key[P1], const double variance){
	Zqp *coefs = new Zqp[P1*P2];
    for (size_t i = 0 ; i < P1*P2 ; ++i)
        coefs[i] = 0;
    for (size_t i = 0 ; i < P1 ; ++i)
        coefs[i*P2] = key[i].get_data()[0];
    Rp12 key2(coefs, false);
	KSW->init(key_pq, &key2, variance);
}

template void gen_bootstrapping_keys<Rp1, N, P1>(const uint64_t encMod, const uint64_t plainMod, Rp1GSW Tsi[N], const Rz key[N], const Rp1 &key_p, const double variance);
template void gen_bootstrapping_keys<Rp2, N, P2>(const uint64_t encMod, const uint64_t plainMod, Rp2GSW Tsi[N], const Rz key[N], const Rp2 &key_p, const double variance);

template void gen_keyswitching_keys<Rp1, P1>(KSWKeyRp1 keys[P1], const Rp1 &key_p, const double variance);
template void gen_keyswitching_keys<Rp2, P2>(KSWKeyRp2 keys[P2], const Rp2 &key_p, const double variance);

template RLWE1<Rp1> ext_exp_inner<Rp1, Zp1, P1>(const uint64_t encMod, const uint64_t plainMod, const size_t length, const Zp1 *Yelem, const Rp1GSW *rgswElem, const KSWKeyRp1 *KSW);
template RLWE1<Rp2> ext_exp_inner<Rp2, Zp2, P2>(const uint64_t encMod, const uint64_t plainMod, const size_t length, const Zp2 *Yelem, const Rp2GSW *rgswElem, const KSWKeyRp2 *KSW);

template RLWE<Rp1_crt, 1> accumulation<Rp1, Zp1, Rp1_crt, P1>(const RLWE<CirculantRing<Zp12, 1, 1>, N> &cipher, const Rp1GSW BK[N], const KSWKeyRp1 KS[P1]);
template RLWE<Rp2_crt, 1> accumulation<Rp2, Zp2, Rp2_crt, P2>(const RLWE<CirculantRing<Zp12, 1, 1>, N> &cipher, const Rp2GSW BK[N], const KSWKeyRp2 KS[P2]);

template Rz dot_product<Rz, P1>(const Rz *const vec1, const Rz *const vec2);
template Rz dot_product<Rz, N>(const Rz *const vec1, const Rz *const vec2);
template Rp1 dot_product<Rp1, 1>(const Rp1 *const vec1, const Rp1 *const vec2);
template Rp2 dot_product<Rp2, 1>(const Rp2 *const vec1, const Rp2 *const vec2);
template Rp12 dot_product<Rp12, 1>(const Rp12 *const vec1, const Rp12 *const vec2);
template Rp12 dot_product<Rp12, 3>(const Rp12 *const vec1, const Rp12 *const vec2);
template Rp1_crt dot_product<Rp1_crt, 1>(const Rp1_crt *const vec1, const Rp1_crt *const vec2);
template Rp2_crt dot_product<Rp2_crt, 1>(const Rp2_crt * const vec1, const Rp2_crt *const vec2);
template Rp12_crt dot_product<Rp12_crt, 3>(const Rp12_crt * const v1, const Rp12_crt *const vec2);
template CirculantRing<Zp12, 1, 1> dot_product<CirculantRing<Zp12, 1, 1>, N>(const CirculantRing<Zp12, 1, 1> *const vec1, const CirculantRing<Zp12, 1, 1> *const vec2);
