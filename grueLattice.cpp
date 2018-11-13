#include <chrono>
#include <ctime>
#include <iostream>
#include <valgrind/callgrind.h>
#include "operations.hpp"
#include "ksw_key.hpp"

int maint(int argc, char *argv[]){
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> keyGen_time, gate_time;
	
	std::cerr << "Generating the gadget..." << std::flush;
	int seed = time(NULL);
	srand(seed);
	Rp1GSW::Gadget_init();
	Rp2GSW::Gadget_init();
	std::cerr << "Done." << std::endl;
	
	uint64_t encMod = Qp;
	uint64_t plainMod = T;
	size_t nb_gate  = INPUT_BIT;
	size_t nb_iter = (argc > 1) ? atoi(argv[1]) : 5;
	
	/** Input data **/
	LWE cipher_i[nb_gate];
	Eva eval;
	Rp12 Ring_f(eval);
	fftw_complex *func = fftw_alloc_complex(FFT_DIM2);
	Ring_f.compute_fft(func);
	
	int64_t coefs[nb_gate];
	if(nb_gate > 10){
		std::cerr << "Computing 63 bits threshold gate." << std::endl;
		for(size_t i = 0; i < nb_gate; ++i)
			coefs[i] = 1;
	}
	else{
		std::cerr << "Computing 6 bits parity gate." << std::endl;
		coefs[0] = 1;
		for(size_t i = 1; i < nb_gate; ++i)
			coefs[i] = 2 * coefs[i-1];
	}
	
	start = std::chrono::system_clock::now();
	/** Generating secret keys **/
	std::cerr << "Generating secret keys..." << std::flush;
	Rz *secret = new Rz[P1];
	for(size_t i = 0; i < P1; ++i)
		secret[i] = Rz::sample_s(DENSITY_KEY);
	Rz *secret2 = new Rz[N];
	for(size_t i = 0; i < N; ++i)
		secret2[i] = Rz::sample_s(DENSITY_KEY_SMALL);
	Rp1 secret_p = Rp1::sample_s(DENSITY_KEY_ACC);
	Rp2 secret_p2 = Rp2::sample_s(DENSITY_KEY_ACC);
	Rp12 secret_pq[3];
	crt_key(secret_pq, secret_p, secret_p2);
	std::cerr << "Done." << std::endl;
	
	/** Derive key material **/
	std::cerr << "Generating LWE-keyswitching key..." << std::flush;
	KSWKeyLWE *secret_lwe = new KSWKeyLWE(secret, secret2, VARIANCE_KSWLWE);
	std::cerr << "Done." << std::endl;
	
	std::cerr << "Generating bootstrapping keys..." << std::flush;
	Rp1GSW *Xsi = new Rp1GSW[P1];
	gen_bootstrapping_keys<Rp1, N, P1>(encMod, plainMod, Xsi, secret2, secret_p, VARIANCE_ACC);
	std::cerr << "1/2" << std::flush;
	Rp2GSW *Ysi = new Rp2GSW[N];
	gen_bootstrapping_keys<Rp2, N, P2>(encMod, plainMod, Ysi, secret2, secret_p2, VARIANCE_ACC);
	std::cerr << "Done." << std::endl;
	
	std::cerr << "Generating keyswitching keys..." << std::flush;
	KSWKeyRp1 *KSW1 = new KSWKeyRp1[P1];
	gen_keyswitching_keys<Rp1, P1>(KSW1, secret_p, VARIANCE_ACC);
	KSWKeyRp2 *KSW2 = new KSWKeyRp2[P2];
	gen_keyswitching_keys<Rp2, P2>(KSW2, secret_p2, VARIANCE_ACC);
	std::cerr << "Done." << std::endl;
	
	std::cerr << "Generating function extraction key..." << std::flush;
	KSWKeyRp12 *Secret = new KSWKeyRp12();
	gen_funexpextract_key(Secret, secret_pq, secret, VARIANCE_FUNEXPEXTRACT);
	std::cerr << "Done." << std::endl;
	
	end = std::chrono::system_clock::now();
	keyGen_time = end - start;
	std::cerr << "key generation time : " << keyGen_time.count() << std::endl;
	
	size_t nb_success = 0;
	std::chrono::duration<double> total_time = std::chrono::duration<double>::zero();
	for(size_t iter = 0; iter < nb_iter; ++iter){
		std::cerr << " ====== Test " << iter << "/" << nb_iter << " ======" << std::endl;
		/** Setup values **/
		std::cerr << "Generating values..." <<std::flush;
		CirculantRing<Zt, 1, 1> message_i[nb_gate];
		for(size_t i = 0; i < nb_gate; ++i){
			message_i[i] = rand() % 2;
			cipher_i[i].encrypt(secret, message_i[i], encMod, plainMod, VARIANCE_INPUT);
		}
		int64_t message = 0;
		for(size_t i = 0; i < nb_gate; ++i)
			message += coefs[i] * message_i[i].get_data()[0].get_value();
		std::cerr << "Done." << std::endl;
		std::cerr << "Computing eval(" << message << ") ->" << eval.eval(message) << std::endl;
		std::cerr << "Timings are " << std::endl;
		
		/** Launching boolean circuits **/
		CALLGRIND_START_INSTRUMENTATION;
		CALLGRIND_TOGGLE_COLLECT;
		start = std::chrono::system_clock::now();
		LWE result = gate(nb_gate, coefs, cipher_i, func, *secret_lwe, Xsi, KSW1, Ysi, KSW2, *Secret);
		end = std::chrono::system_clock::now();
		gate_time = end - start;
		CALLGRIND_TOGGLE_COLLECT;
		CALLGRIND_STOP_INSTRUMENTATION;
		
		CirculantRing<Zt, 1, 1> decrypt;
        result.decrypt(secret, decrypt, encMod, plainMod);
        
        total_time += gate_time;
        std::cerr << "Total gate time " << gate_time.count() << "s" << std::endl;
        std::cout << "Actual result eval(" << message << ") = " << decrypt.get_data()[0] << " (expected " << eval.eval(message) << ")" << std::endl << std::endl;
        if((decrypt - CirculantRing<Zt, 1, 1>(eval.eval(message))).is_zero())
            ++nb_success;
	}
	std::cout << std::endl << "Success: " << nb_success << "/" << nb_iter << " (" << nb_success * 100 / nb_iter << " %), average gate time: " << total_time.count()/nb_iter << "s" <<  std::endl << std::endl;


    fftw_free(func);
    delete[] secret;
    delete Secret;
    delete[] Xsi;
    delete[] Ysi;
    delete[] KSW1;
    delete[] KSW2;
    delete secret_lwe;

    return 0;
}
