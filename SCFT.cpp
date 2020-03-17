//
// Created by callum on 20/01/2020.
//

#include <fstream>
#include "SCFT.h"

SCFT::SCFT(int M_x, int N_s, double L, double N, double R_g, double f, double chiN, double mixing_parameter) : gen(rd()), dist(0.0,1.0),
    M_x(M_x), N_s(N_s), L(L), N(N), R_g(R_g), f(f), chiN(chiN), mixing_parameter(mixing_parameter), q_propagator(M_x, N_s, L, N, R_g, f, chiN),
                                                                                                                      q_star_propagator(M_x, N_s, L, N, R_g, f, chiN),
                                                                                                gauss_noise(0, sqrt(2 * 0.00001)) {//TODO: Add sqrt(2.0.01) as parameter.
    w_A = new double[M_x];
    w_B = new double[M_x];
    delta_s = N / N_s;
    delta_x = L / M_x;
    for(int i = 0; i < M_x; i++) {
        //w_A[i] = dist(gen) * 0.1;
        //w_B[i] = dist(gen) * 0.1;
        const double value = 1.0 * sin(2 * M_PI * (i * delta_x - L / 2.0) / L);
        w_A[i] = -value;
        w_B[i] = value;
    }
    /*w_B[0] = -5;
    int midpoint = M_x / 2;
    w_B[midpoint] = 5;*/

    q_propagator.Set_Fields(w_A, w_B);
    q_star_propagator.Set_Fields(w_A, w_B);
    phi_A = new double[M_x];
    phi_B = new double[M_x];
}

double SCFT::Determine_Error(int index) {
    bool output = false;
    //std::ofstream output_fields("output_fields_" + std::to_string(index) + ".csv");

    double sum = 0.0;
    for(int i = 0; i < M_x; i++) {
        double compressibility_condiiton = 0.5 * (w_A[i] + w_B[i] - chiN);
        double w_a_out = chiN * phi_B[i] + compressibility_condiiton;
        double w_b_out = chiN * phi_A[i] + compressibility_condiiton;
        //output_fields << i << "," << w_A[i] << "," <<  w_a_out << "," << w_B[i] << "," << w_b_out << '\n';
        double difference_in_w_A = fabs(pow(w_a_out - w_A[i], 2.0));
        double difference_in_w_B = fabs(pow(w_b_out - w_B[i], 2.0));
        sum = sum + (difference_in_w_A + difference_in_w_B) * delta_x;
    }
    sum = sqrt(1.0 / (M_x * delta_x) * sum) / std::max(1.0, chiN);
    //output_fields.close();
    return sum;
}

double SCFT::DetermineVariance() {
    double mean_squares = 0.0;
    for(int i = 0; i < M_x; i++) {
        mean_squares += pow(phi_A[i],2.0);
    }
    mean_squares /= M_x;
    double mean = 0.0;
    for(int i = 0; i < M_x; i++) {
        mean += phi_A[i];
    }
    mean /= M_x;
    return mean_squares - mean * mean;
}

double SCFT::DetermineVarianceTotal() {
    double mean_squares = 0.0;
    for(int i = 0; i < M_x; i++) {
        mean_squares += pow(phi_A[i] + phi_B[i],2.0);
    }
    mean_squares /= M_x;
    double mean = 0.0;
    for(int i = 0; i < M_x; i++) {
        mean += phi_A[i] + phi_B[i];
    }
    mean /= M_x;
    return mean_squares - mean * mean;
}

void SCFT::Update_Fields() {
    for(int i = 0; i < M_x; i++) {
        double compressibility_condition = 0.5 * (w_A[i] + w_B[i] - chiN);
        double existing_w_a = w_A[i];
        double existing_w_b = w_B[i];
        double density_B_scalar = phi_B[i];
        double density_A_scalar = phi_A[i];
        double w_a_out = chiN *  phi_B[i] + compressibility_condition;
        double w_b_out = chiN * phi_A[i] + compressibility_condition;
        double w_a = (w_A[i] * (1.0 - mixing_parameter) + mixing_parameter * w_a_out);
        double w_b = (w_B[i] * (1.0 - mixing_parameter) + mixing_parameter * w_b_out);
        w_A[i] = w_a;
        w_B[i] = w_b;
        //std::cout << "bob!" << '\n';
    }

    double average_field = Determine_Average();
    for(int i = 0; i < M_x; i++) {
        w_A[i] -= average_field;
        w_B[i] -= average_field;
    }

}

void SCFT::Run() {

    int index = 1.0;
    double field_error_threshold = pow(10.0,-3.0);
    double variance_threshold = pow(10.0,-5.0);
    double field_error = Determine_Error(index);
    double variance_error = sqrt(DetermineVarianceTotal());
    while(field_error > field_error_threshold || variance_error > variance_threshold) {
        std::cout << index << " - Error: " << Determine_Error(index) << '\n';
        std::cout << index << " - Stdev: " << sqrt(DetermineVarianceTotal()) << '\n';
        std::cout << index << " - Field average: " << Determine_Average() << '\n';
        q_propagator.Propagate();
        q_star_propagator.Propagate();
        Determine_Density_Differences();
        Update_Fields();
        std::string file_name = "hell_" + std::to_string(index) + ".csv";
        if(index % 200 == 0 || index < 20) {
            Save(file_name);
        }
        index++;
        if(index > 2000) {
            break;
        }
        field_error = Determine_Error(index);
        variance_error = sqrt(DetermineVarianceTotal());
    }
    std::cout << index << " - Error: " << field_error << '\n';
    std::cout << index << " - Stdev: " << variance_error << '\n';
    std::cout << index << " - Field average: " << Determine_Average() << '\n';
    std::string file_name = "hell_" + std::to_string(index) + ".csv";
    Save(file_name);
}

void SCFT::Save(std::string file_name) {
    std::ofstream output;
    output.open(file_name);
    output << "x,q,q*,w_A,w_B,phi_A,phi_B,density" << '\n';
    for(int i = 0; i < M_x; i++) {
        double** q = q_propagator.GetPropagator();
        double** q_star = q_propagator.GetPropagator();
        double x = delta_x * i;
        double density = phi_A[i] + phi_B[i];
        output << x << "," << q[i][N_s - 1] << "," << q_star[i][N_s-1] <<","<< w_A[i]  << "," << w_B[i] << "," << phi_A[i] << "," << phi_B[i] << "," << density << '\n';
    }
    output.close();
}

double SCFT::Determine_Average() {
    double w_a_b_average = 0.0;
    for(int i = 0; i < M_x; i++) {
        w_a_b_average += w_A[i] + w_B[i];
    }
    w_a_b_average /= (2.0 * M_x);
    return w_a_b_average;
}

void SCFT::Determine_Density_Differences() {
    double** q = q_propagator.GetPropagator();
    double** q_star = q_star_propagator.GetPropagator();
    double q_Q = q_propagator.Q;
    double q_star_Q = q_star_propagator.Q;
    if(q_Q - q_star_Q > pow(10.0, -4.0)) {
        std::cout << "Warning q_Q != q_star_Q" << '\n';
        std::cout << q_Q << "," << q_star_Q << '\n';
    }
    double Q = q_Q;

    for(int i = 0; i < M_x; i++) {
        double sum = 0.0;
        for(int s = 0; s < f * N_s; s++) {
            sum += q[i][s] * q_star[i][N_s - s] * delta_s;
        }
        phi_A[i] = sum / (Q * N);
        sum = 0.0;
        for(int s = f * N_s; s < N_s; s++) {
            sum += q[i][s] * q_star[i][N_s - s] * delta_s;
        }
        phi_B[i] = sum / (Q * N);
    }

}

void SCFT::output_max() {
    std::cout << "Q propagator max value: " << q_propagator.Find_Max() << '\n';
    std::cout << "Q dagger propagator max value: " << q_propagator.Find_Max() << '\n';
}

void SCFT::output_min() {
    std::cout << "Q propagator min value: " << q_propagator.Find_Min() << '\n';
    std::cout << "Q dagger propagator min value: " << q_propagator.Find_Min() << '\n';
}