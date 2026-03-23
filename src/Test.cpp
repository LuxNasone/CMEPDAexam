#include <iostream>

#include "main.h"

void Minv_test(){

    Float_t pt[100];
    pt[0] = 6;
    pt[1] = 5;

    Float_t eta[100];
    eta[0] = 0.5;
    eta[1] = 1.5;

    Float_t mass[100];
    mass[0] = 0.105;
    mass[1] = 0.105;

    std::cout << M_inv(pt, eta, mass) << std::endl;
}