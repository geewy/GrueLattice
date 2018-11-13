#include "../fft.hpp"
#include "../operations.hpp"
#include "../rlwe.hpp"
#include "../parameters.hpp"

#define NB_TESTS 1000

int main(){
    int seed = time(NULL);
    srand(seed);
    print_param();

    LWE cipher_i[6], res;
    CirculantRing<Zt, 1, 1> message_i[6];
    Rz secret[P1];

    const size_t nb_bit = INPUT_BIT;
    int64_t coefs[k];
    if(nb_bit > 10){
        for(size_t i = 0 ; i < nb_bit ; ++i)
            coefs[i] = 1;
    }
    else{
        coefs[0] = 1;
        for(size_t i = 1 ; i < nb_bit ; ++i)
            coefs[i] = 2 * coefs[i-1];
    }

    double noise = 0;
    double cumul_noise = 0;

    for(size_t idx = 0 ; idx < NB_TESTS ; ++idx){
        for (size_t i = 0 ; i < P1 ; ++i)
            secret[i] = Rz::sample_s(DENSITY_KEY);
        for(size_t i = 0 ; i < 6 ; ++i){
            message_i[i] = rand() % 2;
            cipher_i[i].encrypt(secret, message_i[i], Qp, T, std::pow(2, 80.81));
        }
		res = combination(6, coefs, c_i);
        noise = res.noise(secret, Qp, T);
        cumul_noise += noise;
    }
    std::cout << cumul_noise << " " << cumul_noise / NB_TESTS << " = 2^" << std::log2(cumul_noise / NB_TESTS) << std::endl;

    return 0;
}
