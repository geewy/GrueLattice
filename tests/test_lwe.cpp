#include <iostream>
#include <ctime>
#include <cmath>
#include "../ksw_key.hpp"
#include "../gadget.hpp"
#include "../circulant_ring.hpp"
#include "../predicate.hpp"
#include "../operations.hpp"

int main()
{
    int seed = time(NULL);
    srand(seed);
    Rz secret[P1];

    for(size_t i = 0 ; i < P1 ; ++i)
        secret[i] = Rz::sample_s();

    for(unsigned int i = 0 ; i < 3 ; ++i){
        /** Testing basic encryption/decryption **/
        CirculantRing<Zt, 1, 1> message = CirculantRing<Zt, 1, 1>::uniform_sample();
        LWE lwe;
        lwe.encrypt(secret, message, Qp, T, VARIANCE_INPUT);
        CirculantRing<Zt, 1, 1> decrypt;
        lwe.decrypt(secret, decrypt, Qp, T);
        if(! (message - decrypt).is_zero()){
            std::cerr << "Encrypt/decrypt failed." << std::endl;
            return 1;
        }

        /** Testing addition **/
        CirculantRing<Zt, 1, 1> message1 = CirculantRing<Zt, 1, 1>::uniform_sample();
        LWE lwe1;
        lwe1.encrypt(secret, message1, Qp, T, VARIANCE_INPUT);
        lwe.decrypt(secret, decrypt, Qp, T);
        lwe1 += lwe;
        lwe1.decrypt(secret, decrypt, Qp, T);
        if(!(decrypt - (message + message1)%T).is_zero()){
            std::cerr << "Addition failed." << std::endl;
            return 1;
        }

        /** Testing modulus switching **/
        Rz secret_ms2[N];
        CirculantRing<Zp12, 1, 1> secret_ms[N];
        for(size_t i = 0 ; i < N ; ++i){
            secret_ms2[i] = Rz::sample_s();
            secret_ms[i] = secret_ms2[i];
        }
        RLWE<Rz, N> lwe2;
        RLWE<CirculantRing<Zp12, 1, 1>, N> lwe_ms;
        lwe2.encrypt(secret_ms2, message, Qp, T, VARIANCE_INPUT);
        lwe_ms = lwe2.mod_switch<CirculantRing<Zp12, 1, 1>, P1*P2, Qp>();
        lwe_ms.decrypt(secret_ms, decrypt, P1*P2, T);
        if(! (message - decrypt).is_zero()){
            std::cerr << "ModSwitch failed." << std::endl;
            return 1;
        }

        /** Testing key switching **/
        Rp12LWE lwe3;
        RLWE<Rp12, 1> lwe4;
        Rp12 secret1[3];
        for(size_t j=0 ; j<3 ; ++j)
            secret1[j] = Rp12::sample_s();
        Rp12 secret2[1];
        secret2[0] = Rp12::sample_s();
        CirculantRing<Zt, P1*P2, FFT_DIM2>  message_ksw = CirculantRing<Zt, P1*P2, FFT_DIM2>::uniform_sample();
        lwe3.encrypt(secret1, message_ksw, Qp, T, VARIANCE_INPUT);
        KSWKeyRp12 *ksw = new KSWKeyRp12(secret1, secret2, VARIANCE_FUNEXPEXTRACT);
        lwe4 = lwe3.key_switch(*ksw);
        delete ksw;
        CirculantRing<Zt, P1*P2, FFT_DIM2> decrypt_ksw;
        lwe4.decrypt(secret2, decrypt_ksw, Qp, T);
        if (! (message_ksw - decrypt_ksw).is_zero()){
            std::cerr << "KeySwitch failed." << std::endl;
            return 1;
        }

        /** Testing Galois **/
        uint64_t alpha = 2;
        CirculantRing<Zt, P1> message_galois = CirculantRing<Zt, P1>::uniform_sample();
        Rp1 secret_g[1];
        secret_g[0] = Rp1::sample_s();
        Rp1 secret_galois[1] = {secret_g[0].galois(alpha)};
        Rp1LWE lwe_galois;
        lwe_galois.encrypt(secret_g, message_galois, Q, T, VARIANCE_INPUT);
        lwe_galois.galois_inplace(alpha);
        message_galois.galois_inplace(alpha);
        CirculantRing<Zt, P1> decrypt_galois;
        lwe_galois.decrypt(secret_galois, decrypt_galois, Q, T);
        if(! (decrypt_galois - message_galois).is_zero()){
            std::cerr << "Galois failed." << std::endl;
            return 1;
        }

        size_t message4 = rand() % T;
        size_t message3 = (double) message4 * (double) (P1*P2) / (double) T;
        Zt *Zm_coefs = new Zt[P1*P2];
        for(size_t j = 0 ; j < P1*P2 ; ++j) 
			Zm_coefs[j] = 0;
        Zm_coefs[message3 % (P1*P2)] = 1;
        CirculantRing<Zt, P1*P2, FFT_DIM2> Zm(Zm_coefs, false);
        
        /** Testing ExpCRT **/
        Zt Xmp_coefs[P1];
        for(size_t j = 0 ; j < P1 ; ++j) 
			Xmp_coefs[j] = 0;
        Xmp_coefs[message3 % P1] = 1;
        CirculantRing<Zt, P1> Xmp(Xmp_coefs);
        Zt Ymq_coefs[P2];
        for(size_t j = 0 ; j < P2 ; ++j) 
			Ymq_coefs[j] = 0;
        Ymq_coefs[message3 % P2] = 1;
        CirculantRing<Zt, P2> Ymq(Ymq_coefs);
        Rp1 secret_p = Rp1::sample_s();
        Rp2 secret_q = Rp2::sample_s();
        Rp12 secret_pq[3];
        crt_key(secret_pq, secret_p, secret_q);
        Rp1LWE cipher_p;
        cipher_p.encrypt(&secret_p, Xmp, Q, T, VARIANCE_INPUT);
        RLWE<Rp1_crt, 1> cipher_p_crt = cipher_p.template mod_switch<Rp1_crt, Qcrt, Q>();
        Rp2LWE cipher_q;
        cipher_q.encrypt(&secret_q, Ymq, Q, T, VARIANCE_INPUT);
        RLWE<Rp2_crt, 1> cipher_q_crt = cipher_q.template mod_switch<Rp2_crt, Qcrt, Q>();
        RLWE<Rp12_crt, 3> result;
        exp_crt(result, cipher_p_crt, cipher_q_crt);
        Rp12LWE cipher_pq = result.template mod_switch<Rp12, Qp, Qcrt>();
        CirculantRing<Zt, P1*P2, FFT_DIM2> decrypt;
        cipher_pq.decrypt(secret_pq, decrypt, Qp, T);
        if (! (Zm - decrypt).is_zero()){
            std::cerr << "ExpCRT failed." << std::endl;
            return 1;
        }

        /** Testing function extraction **/
        Fun eval;
        Rp12 Rf(eval);
        fftw_complex *func = fftw_alloc_complex(FFT_DIM2);
        Rf.compute_fft(func);
        Rp12LWE rlwe_eval;
        Rp12 secret_eval[3];
        for(size_t j=0 ; j<3 ; ++j)
            secret_eval[j] = Rp12::sample_s();
        rlwe_eval.encrypt(secret_eval, Zm, Qp, T, VARIANCE_INPUT);
        Rz secret_p_eval[P1];
        for(size_t j = 0 ; j < P1 ; ++j)
			secret_p_eval[j] = Rz::sample_s();
        Rp12 seccret_pp_eval[1];
        Zqp *coefs = new Zqp[P1*P2];
        for(size_t k = 0 ; k < P1*P2 ; ++k)
            coefs[k] = 0;
        for(size_t k = 0 ; k < P1 ; ++k)
            coefs[k * P2] = secret_p_eval[k].get_data()[0];
        secret_pp_eval[0] = Rp12(coefs, false);
        KSWKeyRp12 *Secret = new KSWKeyRp12();
        gen_funexpextract_key(Secret, secret_pq, secret_p_eval, VARIANCE_FUNEXPEXTRACT);
        LWE lwe_eval = fun_exp_extract(func, cipher_pq, *Secret);
        delete Secret;
		CirculantRing<Zt, 1, 1> dec_eval;
        lwe_eval.decrypt(secret_p_eval, dec_eval, Qp, T);
		if(! (dec_eval - CirculantRing<Zt, 1, 1>(eval.eval(message4))).is_zero()){
            std::cerr << "FunExpExtract failed." << std::endl;
            return 1;
        }
    }
    return 0;
}
