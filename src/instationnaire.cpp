#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "instationnaire.hpp"

Instationnaire::Instationnaire(double kappa_prime, double Te_prime, double Phi_p_prime, 
                                double hc_prime, double Lx_prime, double Ly_prime, double Lz_prime, 
                                int M_prime, int Mx_prime, int My_prime, int Mz_prime, 
                                double rho_prime, double Cp_prime, double t_final_prime, 
                                bool flux_constant_prime, int N_prime):
    Stationnaire::Stationnaire(kappa_prime, Te_prime, Phi_p_prime, 
                                hc_prime, Lx_prime, Ly_prime, Lz_prime, 
                                M_prime, Mx_prime, My_prime, Mz_prime){
    rho = rho_prime;
    Cp = Cp_prime;
    t_final = t_final_prime;
    N = N_prime;

    flux_constant = flux_constant_prime;

}

void Instationnaire::solve(){
    double p = 2 * (Ly + Lz); 
    double S = Ly * Lz;
    double h = Lx / M;
    double dt = t_final / N;
    double d = rho * Cp / dt;
    double a_i = -kappa / (h * h);  double a_M = -kappa / h; 
    double b_i = d + 2 * kappa / (h * h) + hc * p / S;  double b_0 = kappa / h; double b_M = kappa / h;
    double c_i = -kappa / (h * h);  double c_0 = - kappa / h;
    
    double f_i = hc * p * Te / S;   double f_0 = Phi_p;

    // Remplisaige de A et F
    A[0][1] = c_0;  A[0][0] = b_0;  F[0] = f_0;
    for (int i = 1; i < M; i++){
        A[i][i+1] = c_i;
        A[i][i] = b_i;
        A[i][i-1] = a_i;

        F[i] = f_i;
    }
    A[M][M] = b_M;  A[M][M-1] = a_M;    F[M] = 0;
    

    // Decomposition LU de A
    std::cout << std::endl << "Facto LU en cours ...";
    factoLU();

    // Verification facto LU
    // verificationFactoLU();
    std::cout << std::endl << "Facto LU terminee";

    /* Le probleme a resoudre de maniere iterative est A*T_k+1 = d*T_k + F */

    // Initialisation
    for (int i = 0; i < M+1; i++){
        T[i] = Te;
    }
    double t = dt;

    // Recherche des intervalles correspondants et ecriture dans le fichier a l'instant 0
    searchInterval();       // La recherche des intervalles se fait une seule fois
    interpolationLineaire();
    bool succes = ecriture3D(0, flux_constant);

    // Iteration et remplissage des differentes solutions
    if (flux_constant){
        std::cout << "\nSCENARIO 1: Le flux Phi_p est constant";
        std::cout << "\nResolution iterative du systeme en cours ...";
        std::cout << "\nEcriture des fichiers '" << "instationnaire_constant.*.vtk" << "' dans le repertoire 'paraview_data' en cours ...";
        std::cout << std::flush;        // Pour des soucis d'affichage
        
        // Les 6 temps auxqels on affiche la solution
        int times[6];
        times[0] = 15;
        times[1] = 30;
        times[2] = 60;
        times[3] = 90;
        times[4] = 150;
        times[5] = 210;
        T_to_plot = new double*[6];
        for (int i = 0; i < 6; i++)
            T_to_plot[i] = new double[M+1];
        
        int k = 0;
        for (int i = 1; i < N+1; i++){
            F[0] = Phi_p;
            // F[M] = 0;
            for (int j = 1; j < M; j++){
                F[j] = f_i + d * T[j];
            }

            // Remplissage des solutions aux differents temps a afficher
            t = i*dt;
            // if (t >= times[k] && t <= times[5]){
            if (k < 6 && t <= times[k] && times[k] < t+dt){
                for (int j = 0; j < M+1; j++)
                    T_to_plot[k][j] = T[j];
                k ++;
            }
            
            // // descneteRemontee n'utlise que les valeurs de T_n+1 en commencant par la derniere. Donc pas besoin de faire T_n = T_n+1!
            // std::cout << "\n" << i << "-eme descente-remontee en cours ...";
            descenteRemontee();
            // verificationSolution();

            // Pour l'ecriture dans le fichier pour la visualisation 
            interpolationLineaire();
            if (succes == false)
                continue;
            else
                succes = ecriture3D(i, true);        // true car le flux est constant
        }        
    }else{
        std::cout << "\nSCENARIO 2: Le flux Phi_p est active et desactive toutes les 30 sec" << "\t";
        std::cout << "\nResolution iterative du systeme en cours ...";
        std::cout << "\nEcriture des fichiers '" << "instationnaire_non_constant.*.vtk" << "' dans le repertoire 'paraview_data' en cours ...";
        std::cout << std::flush;
        
        // Les 3 points en lesquels on visualise la solution
        X_to_plot = new double*[3];
        for (int i = 0; i < 3; i++)
            X_to_plot[i] = new double[N+1];
        X_to_plot[0][0] = T[0];
        X_to_plot[1][0] = T[(int)(M/2)];
        X_to_plot[2][0] = T[M];
        
        int k = 1;
        for (int i = 1; i < N+1; i++){
            // Activation / desactivation du flux
            for (int j = k; j < t_final/30; j++){
                if (i*dt > 30*j){
                    if (j % 2 == 1)
                        F[0] = 0;
                    else
                        F[0] = Phi_p;
                    k ++;
                    break;
                }
            }

            // F[M] = 0;
            // std::cout << "\nt = " << i*dt << "\t\tF[0] = " << F[0];
            for (int j = 1; j < M; j++){
                F[j] = f_i + d * T[j];
            }
            // std::cout << "\n" << i << "-eme descente-remontee en cours ...";
            descenteRemontee();
            // verificationSolution();
            
            // Remplissage des temperatures aux positions qu'on va afficher 1D
            X_to_plot[0][i] = T[0];
            X_to_plot[1][i] = T[(int)(M/2)];
            X_to_plot[2][i] = T[M];

            // Pour l'ecriture dans le fichier pour la visualisation 
            interpolationLineaire();
            if (succes == false)
                continue;
            else
                succes = ecriture3D(i, false);        // false car le flux n'est pas constant
        } 
    }
    if(succes)
        std::cout << "\nEcriture des fichiers vtk OK" << std::flush;
    else
        std::cout << "\nEcriture des fichiers vtk pas OK" << std::flush;
}

void Instationnaire::displayPlot() const{
    // Ecriture dans le fichier pour faire le graphe
    if(flux_constant){
        std::cout << std::endl << "Ecriture du fichier 'data_instationnaire_const.csv' en cours ..." << std::flush;
        std::ofstream write_to_file("data_instationnaire_const.csv");
        if (write_to_file){
            write_to_file << "x," << " T_15, " << " T_30," << " T_60," << " T_90," << " T_150," << " T_210\n";
            for (int i = 0; i < M+1; i++){
                write_to_file << x[i] << ", " << T_to_plot[0][i] << ", " << T_to_plot[1][i] << ", " << T_to_plot[2][i] 
                                << ", " << T_to_plot[3][i] << ", " << T_to_plot[4][i] << ", " << T_to_plot[5][i] << "\n";
            }
            write_to_file.close();

            // tracer du graphique a l'aide de gnuplot
            std::cout << std::endl << "Ecriture du fichier csv OK\nAffichage 1D de la solution instationnaire via 'gnuplot'" << std::flush;
            system("gnuplot gnuplot_command_insta_const.dat");

        } else
            std::cout << "\nErreur lors de l'ecriture du fichier csv, du coup pas de visulation 1D !" << std::flush;
    }else{
        std::cout << std::endl << "Ecriture du fichier 'data_instationnaire_non_const.csv' en cours ..." << std::flush;
        std::ofstream write_to_file("data_instationnaire_non_const.csv");
        write_to_file << "t," << " X_0, " << " X_M/2," << " X_M\n";
        if (write_to_file){
            for (int i = 0; i < N+1; i++){
                write_to_file << i*t_final/N << ", " << X_to_plot[0][i] << ", " << X_to_plot[1][i] << ", " << X_to_plot[2][i] <<"\n";
            }
            write_to_file.close();

            // tracer du graphique a l'aide de gnuplot
            std::cout << std::endl << "Ecriture du fichier csv OK\nAffichage 1D de la solution instationnaire via 'gnuplot'" << std::flush;
            system("gnuplot gnuplot_command_insta_non_const.dat"); 

        } else
            std::cout << std::endl << "Erreur lors de l'ecriture du fichier csv, du coup pas de visulation 1D !" << std::flush;  
    }
}

Instationnaire::~Instationnaire(){
    if (flux_constant){
        for (int i = 0; i < 6; i++){
            delete [] T_to_plot[i];
        }
        delete [] T_to_plot;
    } else{
        for (int i = 0; i < 3; i++){
            delete [] X_to_plot[i];
        }
        delete [] X_to_plot;
    }
}