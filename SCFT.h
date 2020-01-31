//
// Created by callum on 20/01/2020.
//
#include <stdio.h>
#include <iostream>
#include <string>
#include <map>
#include <random>
#include "FieldMethod.h"

#ifndef SINGLE_CHAIN_SCFT_H
#define SINGLE_CHAIN_SCFT_H


#include "Propagator.h"
#include "qPropagator.h"
#include "qDaggerPropagator.h"

class SCFT {
public:
    SCFT(int M_x, int N_s, double L, double N, double R_g, double f, double flory_huggins, double mixing_parameter, FieldMethod field_type);
    void Run();
    void Determine_Density_Differences();
    double Determine_Error();
    void Save(std::string file_name);
    void Cleanup();
    double DetermineVariance();
    double DetermineVarianceTotal();
private:
    FieldMethod field_method;
    qPropagator q_propagator;
    qDaggerPropagator q_star_propagator;
    int M_x,N_s;
    double L,  N,  R_g,  f,  flory_huggins, mixing_parameter;
    double delta_x, delta_s;
    double *w_A, *w_B, *phi_A, *phi_B;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<> dist;
    std::default_random_engine norm_generator;
    std::normal_distribution<double> gauss_noise;
    void Update_Fields();

};


#endif //SINGLE_CHAIN_SCFT_H
