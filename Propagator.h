#ifndef SINGLE_CHAIN_PROPAGATOR_H
#define SINGLE_CHAIN_PROPAGATOR_H


#include <fftw3.h>

class Propagator {
public:
    Propagator(int M_x, int N_s, double L, double N,
            double R_g, double f, double chiN);

    virtual double w(double x, int n);
    void Save(const char *file_name);

    double** GetPropagator();

    ~Propagator();

    double DetermineQ();
    double Find_Max();
    double Find_Min();

    void Output_Parameters();
    void Set_Fields(double* w_A, double* w_B);

    void Propagate();
    double Q;

protected:
    double delta_x, delta_s,b, R_g,L,f,N, chiN, s;
    int M_x, N_s;
    double* w_A;
    double* w_B;
private:
    void Apply_W_To_Q(double** q, int src_N_s, int dest_N_s);
    void Propagate_n(int n);
    double** q;
    fftw_complex *a_j, *h_j, *q_complex;
    fftw_plan p1,p2;
};


#endif //SINGLE_CHAIN_PROPAGATOR_H
