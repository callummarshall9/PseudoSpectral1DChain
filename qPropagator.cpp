//
// Created by callum on 20/01/2020.
//

#include <iostream>
#include <cmath>
#include "mathutils.h"
#include "qPropagator.h"

qPropagator::qPropagator(int M_x, int N_s, double L, double N, double R_g, double f, double chiN) : Propagator(M_x, N_s, L, N, R_g, f, chiN) {

}



double qPropagator::w(double x, int n) {
    double field = 10.0 * sin(2 * M_PI * (x - L / 2.0) / L) / N;
return field;
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

