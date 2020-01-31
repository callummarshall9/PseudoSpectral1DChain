//
// Created by callum on 20/01/2020.
//

#include <iostream>
#include <cmath>
#include "qPropagator.h"

qPropagator::qPropagator(int M_x, int N_s, double L, double N, double R_g, double f, double flory_higgs, FieldMethod field_method) : Propagator(M_x, N_s, L, N, R_g, f, flory_higgs), field_method(field_method) {

}

inline double sech(double x) {
    return 1.0 / cosh(x);
}

double qPropagator::w(double x, int n) {
    int index = x / this->GetDeltaX();
    if(field_method == SimpleMixing) {
        if(n < (this->GetN_s() - this->GetF() * this->GetN_s())) {
            return w_A[index];
        } else {
            return w_B[index];
        }
    } else if(field_method == Langevindgn) {
        if(n < this->GetF() * this->GetN_s()) {
            return (-w_A[index] - w_B[index])  / this->GetN();
        } else {
            return (-w_A[index] + w_B[index])  / this->GetN();
        }
    } else {
        std::cout << "Error: SCFT method not specified." << '\n';
        return 0;
    }
}

void qPropagator::Set_Fields(double *w_A, double *w_B) {
    this->w_A = w_A;
    this->w_B = w_B;
}
