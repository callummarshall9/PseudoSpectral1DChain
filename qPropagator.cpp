//
// Created by callum on 20/01/2020.
//

#include <iostream>
#include <cmath>
#include "mathutils.h"
#include "qPropagator.h"

qPropagator::qPropagator(int M_x, int N_s, double L, double N, double R_g, double f, double flory_higgs, FieldMethod field_method) : Propagator(M_x, N_s, L, N, R_g, f, flory_higgs, field_method) {

}



double qPropagator::w(double x, int n) {
    int index = x / delta_x;
    //Scale the fields by N.
    //Between [0,fN] field = w_A
    //Between [fN,N] field = w_B
    //q starts from s=0
    if(s < f * N) {
        if(field_method == SimpleMixing) {
            return w_A[index] / N;
        } else if(field_method == Langevindgn) {
            return (-w_A[index] - w_B[index]) / N;
        } else {
            std::cout << "Error: SCFT method not specified." << '\n';
            return 0;
        }
    } else {
        if(field_method == SimpleMixing) {
            return w_B[index] / N;
        } else if(field_method == Langevindgn) {
            return (-w_A[index] + w_B[index]) / N;
        } else {
            std::cout << "Error: SCFT method not specified." << '\n';
            return 0;
        }
    }
}

