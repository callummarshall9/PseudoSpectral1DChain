//
// Created by callum on 20/01/2020.
//

#include <iostream>
#include <cmath>
#include "mathutils.h"
#include "qDaggerPropagator.h"

qDaggerPropagator::qDaggerPropagator(int M_x, int N_s, double L, double N, double R_g, double f, double flory_higgs) : Propagator(M_x, N_s, L, N, R_g, f, flory_higgs) {

}


double qDaggerPropagator::w(double x, int n) {
    int index = x / delta_x;
    //Scale the fields by N.
    //Between [0,fN] field = w_A
    //Between [fN,N] field = w_B
    //q dagger starts from s=N
    if(s < (N - f * N)) {//f=0.7, s=0.7N. from other way, s=N-0.7N=0.3N
        return w_B[index] / N;
    } else {
        return w_A[index] / N;

    }
}


