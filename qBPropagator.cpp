//
// Created by callum on 20/01/2020.
//

#include <iostream>
#include <cmath>
#include "qBPropagator.h"

qBPropagator::qBPropagator(int M_x, int N_s, double L, double N, double R_g, double f, double flory_higgs) : Propagator(M_x, N_s, L, N, R_g, f, flory_higgs) {

}

inline double sech(double x) {
    return 1.0 / cosh(x);
}

double qBPropagator::w(double x, int n) {
    /*int index = x / this->GetDeltaX();
    if(n < this->GetF() * this->GetN_s()) {
        return w_B[index] ;
    } else {
        return w_A[index];
    }*/
    return (1.0  - 2.0 * pow(sech(3.0 * (x - this->GetL() / 2.0)/ (2 * this->GetRg())),2.0)) / this->GetN();
}

void qBPropagator::Set_Fields(double *w_A, double *w_B) {
    this->w_A = w_A;
    this->w_B = w_B;
}
