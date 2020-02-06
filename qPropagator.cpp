//
// Created by callum on 20/01/2020.
//

#include <iostream>
#include <cmath>
#include "mathutils.h"
#include "qPropagator.h"

qPropagator::qPropagator(int M_x, int N_s, double L, double N, double R_g, double f, double flory_higgs) : Propagator(M_x, N_s, L, N, R_g, f, flory_higgs) {

}



double qPropagator::w(double x, int n) {
    int index = x / delta_x;
    //Scale the fields by N.
    //Between [0,fN] field = w_A
    //Between [fN,N] field = w_B
    //q starts from s=0
    if(s < f * N) {
        return w_A[index] / N;
    } else {
        return w_B[index] / N;

    }
}

