//
// Created by callum on 20/01/2020.
//

#ifndef SINGLE_CHAIN_QPROPAGATOR_H
#define SINGLE_CHAIN_QPROPAGATOR_H

#include "FieldMethod.h"
#include "Propagator.h"

class qPropagator : public Propagator {
public:
    qPropagator(int M_x, int N_s, double L, double N, double R_g, double f, double flory_higgs, FieldMethod field_method);
    double w(double x, int n) override ;
};


#endif //SINGLE_CHAIN_QPROPAGATOR_H
