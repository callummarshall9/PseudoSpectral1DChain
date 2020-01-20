//
// Created by callum on 20/01/2020.
//

#ifndef SINGLE_CHAIN_QAPROPAGATOR_H
#define SINGLE_CHAIN_QAPROPAGATOR_H


#include "Propagator.h"

class qAPropagator : public Propagator {
public:
    qAPropagator(int M_x, int N_s, double L, double N, double R_g, double f, double flory_higgs);
    double w(double x, int n) override ;
    void Set_Fields(double* w_A, double* w_B);
private:
    double *w_A, *w_B;
};


#endif //SINGLE_CHAIN_QAPROPAGATOR_H
