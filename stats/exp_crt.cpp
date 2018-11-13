#include "../rlwe.hpp"
#include "../operations.hpp"
#include "../parameters.hpp"

#define NB_TESTS 100

int main(){
    int seed = time(NULL);
    srand(seed);
    print_param();

    size_t message;
    Zt Xmp_coefs[P1], Ymq_coefs[P2];
    Rp1 secret_p_tmp[1];
    Rp2 secret_q_tmp[1];
    Rp12 secret_pq_tmp[3];
    Rp1_crt secret_p[1];
    Rp2_crt secret_q[1];
    Rp12_crt secret_pq[3];
    RLWE<Rp1_crt, 1> cipher_p;
    RLWE<Rp2_crt, 1> cipher_q;
    RLWE<Rp12_crt, 3> cipher_pq;

    double noise = 0;
    double cumul_noise = 0;

    for (size_t idx = 0 ; idx < NB_TESTS ; ++idx){
        message = rand() % (P1*P2);

        for(size_t j = 0 ; j < P1 ; ++j) 
			Xmp_coefs[j] = 0;
        Xmp_coefs[message % P1] = 1;
        CirculantRing<Zt, P1> Xmp(Xmp_coefs);

        for(size_t j = 0 ; j < P2 ; ++j) 
			Ymq_coefs[j] = 0;
        Ymq_coefs[message % P2] = 1;
        CirculantRing<Zt, P2> Ymq(Ymq_coefs);

        secret_p[0] = Rp1_crt::sample_s(DENSITY_KEY);
        secret_q[0] = Rp2_crt::sample_s(DENSITY_KEY);
        secret_p_tmp[0] = secret_p[0];
        secret_q_tmp[0] = secret_q[0];

        cipher_p.encrypt(secret_p, Xmp, Qcrt, T, std::pow(2, 6.41));
        cipher_q.encrypt(secret_q, Ymq, Qcrt, T, std::pow(2, 6.41));

        exp_crt(cipher_pq, cipher_p, cipher_q);
        crt_key(secret_pq_tmp, secret_p_tmp[0], secret_q_tmp[0]);
        secret_pq[0] = secret_pq_tmp[0];
        secret_pq[1] = secret_pq_tmp[1];
        secret_pq[2] = secret_pq_tmp[2];

        noise = cipher_pq.noise(secret_pq, Qcrt, T);
        cumul_noise += noise;
    }
    std::cout << cumul_noise << " " << cumul_noise / NB_TESTS / (P1*P2) << " = 2^" << std::log2(cumul_noise / NB_TESTS / (P1*P2)) << std::endl;

    return 0;
}
