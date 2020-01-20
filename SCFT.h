//
// Created by callum on 20/01/2020.
//
#include <stdio.h>
#include <iostream>
#include <string>
#include <map>
#include <random>

#ifndef SINGLE_CHAIN_SCFT_H
#define SINGLE_CHAIN_SCFT_H


#include "Propagator.h"
#include "qAPropagator.h"
#include "qBPropagator.h"

class SCFT {
public:
    SCFT(int M_x, int N_s, double L, double N, double R_g, double f, double flory_higgs, double mixing_parameter);
    void Run();
    void Determine_Density_Differences();
    double Determine_Error();
    void Save(std::string file_name);

private:
    qAPropagator q_propagator;
    qBPropagator q_star_propagator;
    int M_x,N_s;
    double L,  N,  R_g,  f,  flory_higgs, mixing_parameter;
    double delta_x, delta_s;
    double *w_A, *w_B, *phi_A, *phi_B;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<> dist;
    void Update_Fields();

};


#endif //SINGLE_CHAIN_SCFT_H
