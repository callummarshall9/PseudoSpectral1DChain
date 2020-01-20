//
// Created by callum on 20/01/2020.
//

#ifndef SINGLE_CHAIN_PROPAGATOR_H
#define SINGLE_CHAIN_PROPAGATOR_H


#include <fftw3.h>

class Propagator {
public:
    Propagator(int M_x, int N_s, double L, double N, double R_g, double f, double flory_higgs);

    virtual double w(double x, int n);
    void Save(const char *file_name);

    double GetF();
    double GetDeltaX();
    int GetN_s();
    double** GetPropagator();
    double GetQ();
    double GetL();
    double GetRg();
    double GetN();

    double DetermineQ();

    void Output_Parameters();

    void Propagate();

private:
    void Apply_Q_Operator(double** q, int src_N_s, int dest_N_s);
    void Propagate_n(int n);
    int M_x, N_s;
    double delta_x, delta_s,b, R_g,Q,L,f,N, flory_higgs;

    double** q;
    fftw_complex *a_j, *h_j, *q_complex;
    fftw_plan p1,p2;
};


#endif //SINGLE_CHAIN_PROPAGATOR_H
