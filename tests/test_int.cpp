#include <iostream>
#include "../integer_mod.hpp"

int main (){
    Zq int_a(10);
    Zq int_b(3);
    Zq int_c(0);
    Zq int_d(32);
    Zq int_e(-32);
    Zq int_eh, int_el;
    int mod = 11;

    std::cout << "int_a " << int_a << std::endl;
    std::cout << "int_b " << int_b << std::endl;
    std::cout << "int_c " << int_c << std::endl;
    std::cout << "int_d " << int_d << std::endl;
    std::cout << "int_e " << int_e << std::endl;
    std::cout << "int_a - int_b " << int_a - int_b << std::endl;
    std::cout << "int_b - int_a " << int_b - int_a << std::endl;
    std::cout << "int_a * int_b " << int_a * int_b << std::endl;
    std::cout << "int_a % mod " << int_a % mod << std::endl;
    std::cout << "int_b % mod " << int_b % mod << std::endl;
    Zq res = int_a * 2;
    Zq res0 = int_b * 2;
    std::cout << "int_a * 2 " << int_a * 2 << std::endl;
    std::cout << "int_b * 2 " << int_b * 2 << std::endl;
    std::cout << "res " << res << std::endl;
    std::cout << "res0 " << res0 << std::endl;
    Zq res1 = int_a / 2;
    Zq res2 = int_b / 2;
    std::cout << "int_a / 2 " << int_a / 2 << std::endl;
    std::cout << "int_b / 2 " << int_b / 2 << std::endl;
    std::cout << "res1 " << res1 << std::endl;
    std::cout << "res2 " << res2 << std::endl;
    Zq res3 = (int_a/2) % 5;
    Zq res4 = (int_b/2) % 5;
    std::cout << "(int_a / 2) % 5 " << (int_a/2)%5 << std::endl;
    std::cout << "(int_b / 2) % 5 " << (int_b/2)%5 << std::endl;
    std::cout << "res3 " << res3 << std::endl;
    std::cout << "res4 " << res4 << std::endl;

    int64_t constant = int_e.split(int_eh, int_el);
    std::cout << int_e << " = " << int_eh << " * " << constant << " + " << int_el << std::endl;

    return 0;
}
