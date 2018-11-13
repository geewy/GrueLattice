#include <iostream>
#include <ctime>

#include "../circulant_ring.hpp"
#include "../parameters.hpp"
#include "../gadget.hpp"
#include "../rgsw.hpp"

typedef CirculantRing<Zt, P1> MyRpt;

int main(){
    int seed = time(NULL);
    srand(seed);
    Rp1GSW::Gadget_init();

    Rp1 secret[N];
    for(size_t i = 0 ; i < N ; ++i)
        secret[i] = Rp1::sample_s();

    for(unsigned int i = 0 ; i < 10 ; ++i){
        /** Testing basic encryption/decryption **/
        MyRpt message = MyRpt::uniform_sample();
        Rp1GSW rgsw;
        rgsw.encrypt(secret, message, VARIANCE_INPUT);
        MyRpt decrypt;
        rgsw.decrypt(secret, decrypt);

        if(! (message - decrypt).is_zero()){
            std::cerr << "Encrypt/decrypt failed." << std::endl;
            return 1;
        }
    }
    return 0;
}
