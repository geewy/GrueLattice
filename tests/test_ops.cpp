#include <iostream>
#include <ctime>
#include <cstring>

#include "../operations.hpp"

typedef CirculantRing<Zt, P1> MyRpt;

int main(){
    int seed = time(NULL);
    srand(seed);
    Rp1GSW::Gadget_init();

    Rp1 secret[1];
    secret[0] = Rp1::sample_s();

    for (unsigned int idx = 0 ; idx < 5 ; ++idx){
        /** Testing basic encryption and decryption **/
        size_t message = rand() % (P1);
        Zt Xm_coefs[P1];
        for(size_t j = 0 ; j < P1 ; ++j) 
			Xm_coefs[j] = 0;
        Xm_coefs[message] = 1;
        MyRpt Xm(Xm_coefs);
        Rp1LWE rlwe;
        rlwe.encrypt(secret, Xm, Q, T, VARIANCE_INPUT);
        MyRpt decrypt;
        rlwe.decrypt(secret, decrypt, Q, T);
        if(! (Xm - derypt).is_zero()){
            std::cerr << "RLWE Encrypt/decrypt failed." << std::endl;
            return 1;
        }

        size_t message_p = rand() % (P1);
        Zt Xmp_coefs[P1];
        for(size_t j = 0 ; j < P1 ; ++j) 
			Xmp_coefs[j] = 0;
        Xmp_coefs[message_p] = 1;
        MyRpt Xmp(Xmp_coefs);
        Rp1GSW rgsw;
        rgsw.encrypt(secret, Xmp, VARIANCE_INPUT);
        rgsw.decrypt(secret, decrypt);
        if(! (Xmp - decrypt).is_zero()){
            std::cerr << "RGSW Encrypt/decrypt failed." << std::endl;
            return 1;
        }

        /** Testing ExtMult **/
        RLWE<Rp1, 1> rlwe2 = rlwe.ext_mult(rgsw);
        Zt Xm_mp_coefs[P1];
        for(size_t j = 0 ; j < P1 ; ++j) 
			Xm_mp_coefs[j] = 0;
        Xm_mp_coefs[(message + message_p) % P1] = 1;
        MyRpt Xm_mp(Xm_mp_coefs);
        rlwe2.decrypt(secret, decrypt, Q, T);
        if(! (Xm_mp - decrypt).is_zero()){
            std::cerr << "ExtMult failed." << std::endl;
            return 1;
        }

        /** Testing ExtExpMultAdd **/
        uint64_t alpha = 2;
        uint64_t beta = (size_t)Zp1(alpha).inv();
        Rp1 secret_alpha[1], secret_beta[1];
        secret_alpha[0] = secret[0].galois(alpha);
        secret_beta[0]  = secret[0].galois(beta);
        KSWKeyRp1 Secret_alpha(secret_alpha, secret, VARIANCE_ACC);
        KSWKeyRp1 Secret_beta(secret_beta, secret, VARIANCE_ACC);
        rlwe.ext_exp_mult_add(rgsw, alpha, S_alpha, beta, S_beta);
        Zt Xam_mp_coefs[P1];
        for(size_t j = 0 ; j < P1 ; ++j) 
			Xam_mp_coefs[j] = 0;
        Xam_mp_coefs[(alpha * message_p + message) % P1] = 1;
        MyRpt Xam_mp(Xam_mp_coefs);
        rlwe.decrypt(secret, decrypt, Q, T);
        if(! (Xam_mp - decrypt).is_zero()){
            std::cerr << "ExtExpMultAdd failed." << std::endl;
            return 1;
        }

        /** Testing ExtExpInner **/
        KSWKeyRp1 Secret[P1];
        gen_keyswitching_keys<Rp1, P1>(Secret, secret[0], VARIANCE_ACC);
        Rp1GSW Cipher[N];
        Zqp x1[N];
        Rz secret_x[N];
        Zp1 x2[N];
        for(size_t i = 0 ; i < N ; ++i){
            x1[i] = rand();
            x2[i] = x1[i] % P1;
            secret_x[i] = x1[i];
        }
        gen_bootstrapping_keys<Rp1, N, P1>(Q, T, C, secret_x, secret[0], VARIANCE_ACC);
        Zp1 y[N];
        for(size_t i = 0 ; i < N ; ++i)
            y[i] = rand();
		Zp1 xy = 0;
        for(size_t i = 0 ; i < N ; ++i)
            xy += x2[i] * y[i];
        Zt sol_coefs[P1];
        for(size_t j = 0 ; j < P1 ; ++j) 
			sol_coefs[j] = 0;
        sol_coefs[(size_t)xy] = 1;
        MyRpt sol(sol_coefs);
        Rp1LWE cipher = ext_exp_inner<Rp1, Zp1, P1>(Q, T, N, y, Cipher, Secret);
        cipher.decrypt(secret, decrypt, Q, T);
        if(! (sol - decrypt).is_zero()){
            std::cerr << "ExtExpInner failed." << std::endl;
            return 1;
        }
    }
    return 0;
}
