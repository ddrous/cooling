#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include "stationnaire.hpp"

Stationnaire::Stationnaire(double kappa_prime, double Te_prime, double Phi_p_prime, 
                        double hc_prime, double Lx_prime, double Ly_prime, double Lz_prime, int M_prime,
                        int Mx_prime, int My_prime, int Mz_prime){
    // Parametres physiques du systeme
    kappa = kappa_prime;
    Te = Te_prime;
    Phi_p = Phi_p_prime;
    hc = hc_prime;

    // Parametres geometriques du systeme
    Lx = Lx_prime;
    Ly = Ly_prime;
    Lz = Lz_prime;

    // Parametres de simulation
    M = M_prime;

    // Differents points de discretisation
    double h = Lx / M; 
    x = new double[M+1];
    for (int i = 0; i < M+1; i++)
        x[i] = i * h;
    
    // Initialisation des matrice A, L et U
    A = new double*[M+1];
    for (int i = 0; i < M+1; i++){
        A[i] = new double[M+1];
    }
    for (int i = 0; i < M+1; i++){
        for (int j = 0; j < M+1; j++)
            A[i][j] = 0;            // Pour eviter les surprises!        
    }

    L = new double*[M+1];
    for (int i = 0; i < M+1; i++){
        L[i] = new double[M+1];     // L et U ne sont pas caree, mmalheuresement. Mais necessaire pour la verification!
    }

    U = new double*[M+1];
    for (int i = 0; i < M+1; i++){
        U[i] = new double[M+1];
    }

    // for (int i = 0; i < M+1; i++){
    //     for (int j = 0; j < M+1; j++)
    //         A[i][j] = 0;
    // }
    
    // Initialisation du vecteur T = X (voire enonce)
    T = new double[M+1];
    for (int i = 0; i < M+1; i++){
        T[i] = 50;
    }
    

    // Initialisation du vecteur second membre F
    F = new double[M+1]; 

    // Initialisation du vecteur analytique Texact
    Texact = new double[M+1];

    // Initialisation des parametres de visualisation 3D
    Mx = Mx_prime;     
    My = My_prime;     
    Mz = Mz_prime;
    T_hat = new double[Mx+1];
    x_hat = new double[Mx+1];
    for (int i = 0; i < Mx+1; i++)
        x_hat[i] = i * Lx / Mx;      
    x_hat_indices = new int[Mx+1];
}

// Fonction specifiques pour la facto LU
void Stationnaire::factoLU(){
    double a_i = A[1][0];   double b_i = A[1][1];   double c_i = A[1][2];
    double b_etoile = A[0][0];  double c_etoile = A[0][1] / b_etoile;       // b_etoile et c_etoile juste temporaires (a refaire !)
    L[0][0] = b_etoile;     U[0][1] = c_etoile;
    L[1][0] = A[1][0];          U[0][0] = 1;
    for (int i = 1; i < M; i++){
        b_etoile = b_i - a_i * U[i-1][i];
        c_etoile = c_i / b_etoile;
        L[i][i] = b_etoile;     U[i][i+1] = c_etoile;
        L[i][i-1] = a_i;        U[i][i] = 1;
    }
    b_etoile = A[M][M] - A[M][M-1] * U[M-1][M];
    L[M][M-1] = A[M][M-1];
    L[M][M] = b_etoile;         U[M][M] = 1;

    // // Affichage de L
    // std::cout << "\n\nMatrice L";
    // for (int i = 0; i < M+1; i++){
    //     std::cout << "\n";
    //     for (int j = 0; j < M+1; j++)
    //         std::cout << std::setw(8) << std::left << L[i][j] << "\t";
    // }

    // // Affichage de U
    // std::cout << "\n\nMatrice U";
    // for (int i = 0; i < M+1; i++){
    //     std::cout << "\n";
    //     for (int j = 0; j < M+1; j++)
    //         std::cout << std::setw(8) << std::left << U[i][j] << "\t";
    // }
}

// Fonction specifiques pour la verification de la facto LU
void Stationnaire::verificationFactoLU() const{
    std::cout << "\nFacto LU en cours ...";

    double ** LU = new double *[M+1];
    for (int i = 0; i < M+1; i++)
        LU[i] = new double[M+1];
    for (int i = 0; i < M+1; i++){
        for (int j = 0; j < M+1; j++){
            for (int k = 0; k < M+1; k++)
                LU[i][j] += L[i][k] * U[k][j];
        }
    }
    bool flag = true;
    for (int i = 0; i < M+1; i++){
        for (int j = 0; j < M+1; j++){
            if (A[i][j] != LU[i][j]){
                flag = false;
                break;
            }
        }
        if (flag)
            break;             
    }
    if (flag)
        std::cout << "\nFacto LU OK";
    else
        std::cout << "\nFacto LU pas OK";
    // // Affichage de LU
    // std::cout << "\n\nMatrice LU";
    // for (int i = 0; i < M+1; i++){
    //     std::cout << "\n";
    //     for (int j = 0; j < M+1; j++)
    //         std::cout << std::setw(8) << std::left << LU[i][j] << "\t";
    // }

    // Supression du pointeur inutile
    for (int i = 0; i < M+1; i++){
        delete [] LU[i];
    }
    delete [] LU;
    
}

// Fonction specifiques pour la resolution lineaire tridiagonale
void Stationnaire::descenteRemontee(){
    double *Y = new double[M+1];                                // Ne pas oublier de retirer la declaration de Y (ineficace pour instationnaire)
    // Resolution de LY = F (Descente)
    Y[0] = F[0] / L[0][0];
    for (int i = 1; i < M+1; i++)
        Y[i] = (F[i] - L[i][i-1] * Y[i-1]) / L[i][i];

    // Resolution de UT = Y (Remontee)
    T[M] = Y[M];
    for (int i = M-1; i > -1; i--)
        T[i] = Y[i] - U[i][i+1] * T[i+1];

    // On supprime Y
    delete [] Y; 
}

// Fonction qui recherche les intervalles pour l'interpolation
void Stationnaire::searchInterval(){
    /* Recherche de j tel que x_i00 = x_hat[i] appatient a [x[j]; x[j+1]
        et le resultat est stocke dans x_hat_indices*/
    int j;
    int k = 0;
    for (int i = 0; i < Mx; i++){
        // std::cout << x_hat[i] << "\n";
        // for (j = k; j < M; j++){
        //     if (x[j] <= x_hat[i] && x_hat[i] < x[j+1]){
        //         // std::cout << "\n" <<  x[j] << ", " << x_hat[i] << ", " << x[j+1];
        //         // k++;     // On incremente k dans les cas ou Mx > M afin d'accelerer la recherche des intervalles
        //         break;
        //     }
        // }
        // x_hat_indices[i] = j;

        j = 1;
        while (j < M && x_hat[i] > x[j]){
            j += 1;
        }
        // std::cout << "\n" <<  x[j-1] << ", " << x_hat[i] << ", " << x[j];
        x_hat_indices[i] = j-1;
    }
    x_hat_indices[Mx] = M;
    // std::cout << "\n" << x_hat[Mx] << ", " << x[M];
}

// Fonction qui va calculer l'interpolant lineaire T_hat
void Stationnaire::interpolationLineaire(){
    int j;
    double a, b;
    for (int i = 0; i < Mx; i++){
        j = x_hat_indices[i];
        a = (T[j+1] - T[j])/(x[j+1] - x[j]);
        b = (T[j] * x[j+1] - T[j+1] * x[j])/(x[j+1] - x[j]);
        
        T_hat[i] = a * x_hat[i] + b;
        // std::cout << T[j] << ", " << T_hat[i] << ", " << T[j+1] << "\n";        
    }
    // Ne pas oublier les dernier points x_hat[Mx] et x[M] qui coincident
    T_hat[Mx] = T[M];
    // std::cout << T_hat[Mx] << ", " << T[M] << "\n";
}

// Fonction pour ecrire les resultats 3d
bool Stationnaire::ecriture3D(int indice_iteration, bool flux_de_chaleur_constant){   
    // Ecriture dans le fichier pour faire la visualisation avec paraview

    // Variable a retourner
    bool succes = true;

    int nb_points = (Mx)*(My)*(Mz);     // en realite, nb_points = (Mx+1)*(My+1)*(Mz+1)
    
    std::string nom_du_fichier = "";
    std::string nom_de_solution = "";
    if (indice_iteration == -1){
        nom_du_fichier = "stationnaire.vtk";
        nom_de_solution = "Temperature_Stationnaire";
    }
    else{
        if(flux_de_chaleur_constant){
            nom_du_fichier = "instationnaire_constant." + std::to_string(indice_iteration) + ".vtk";
            nom_de_solution = "Temperature_Instationnaire_Scenario_1";
        }
        else{
            nom_du_fichier = "instationnaire_non_constant." + std::to_string(indice_iteration) + ".vtk";
            nom_de_solution = "Temperature_Instationnaire_Scenario_2";
        }
    }

    // std::cout << "\nEcriture dans le fichier '" << nom_du_fichier << "' en cours ...";
    std::ofstream for_paraview("paraview_data/" + nom_du_fichier);
    if(for_paraview){
        for_paraview << "# vtk DataFile Version 2.0\nvtk output\nASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS " ;
        for_paraview << Mx << " " << My << " " << Mz << "\n";
        for_paraview << "POINTS " << nb_points << " float\n";
        for (int k = 0; k < Mz; k++){
            for (int j = 0; j < My; j++){
                for (int i = 0; i < Mx; i++)
                    for_paraview << i << " " << j << " " << k << "\n";
            } 
        }
        for_paraview << "POINT_DATA " << nb_points << "\nFIELD FieldData 1\n" + nom_de_solution << " 1 " << nb_points << " float \n";
        for (int k = 0; k < Mz; k++){
            for (int j = 0; j < My; j++){
                for (int i = 0; i < Mx; i++)
                            for_paraview << T_hat[i] << "\n";
            } 
        }
        // std::cout << "\nEcriture OK";
        succes = true;
    }else{
        // std::cout << "\nErreur lors de l'ecriture";
        succes = false;
    }
    for_paraview.close();

    return succes;
};

// Fonction specifiques pour la verification de la solution
void Stationnaire::verificationSolution() const{
    // double ** LU = new double *[M+1];
    // for (int i = 0; i < M+1; i++)
    //     LU[i] = new double[M+1];
   
    // for (int i = 0; i < M+1; i++){
    //     for (int j = 0; j < M+1; j++){
    //         for (int k = 0; k < M+1; k++)
    //             LU[i][j] += L[i][k] * U[k][j];
    //     }
    // }

    double * F_prime = new double[M+1];
    for (int i = 0; i < M+1; i++){
        F_prime[i] = 0;
        for (int j = 0; j < M+1; j++)
            F_prime[i] += A[i][j] * T[j];
    }
    
    double flag = true;
    for (int i = 0; i < M+1; i++){
        if (fabs((F_prime[i] - F[i])) > 1e-3){
            flag = false;
            break;
        }
    }
    if (flag)
        std::cout << "\nSolution OK";
    else
        std::cout << "\nSolution pas OK";
    // // Affichage de L*U*T
    // std::cout << "\n\nL*U*T\n";
    // for (int i = 0; i < M+1; i++)
    //     std::cout << std::setw(8) << F_prime[i] << "\t";
    // // Affichage de F
    // std::cout << "\n\nSecond membre F\n";
    // for (int i = 0; i < M+1; i++)
    //     std::cout << std::setw(8) << F[i] << "\t";

    // On suprrime F
    delete [] F_prime;
}

void Stationnaire::solve(){
    double p = 2 * (Ly + Lz);       // et non p = 2 * (Lx + Ly);
    double S = Ly * Lz;
    double h = Lx / M;
    double a_i = -kappa / (h * h);  double a_M = -kappa / h; 
    double b_i = 2 * kappa / (h * h) + hc * p / S;  double b_0 = kappa / h; double b_M = kappa / h;
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

    // Resolution de LY = F
    std::cout << "\nResolution du systeme lineaire en cours ...";
    descenteRemontee();

    // Verification AT = F
    // verificationSolution();

    // Recherche de l'intervalle correspondat a chaque point x_i00
    searchInterval();

    // Calcul de l'interpolnat lineaire
    interpolationLineaire();

    std::cout << "\nEcriture du fichier '" << "stationnaire.vtk" << "' dans le repertoire 'paraview_data' en cours ...";
    // Ecriture 3D de la solution
    bool succes = ecriture3D(-1, true);     // -1 correspond a la solution stationnaire, solution a l'infini, et TRUE puisque le flux de chaleur ne varie pas
    if(succes)
        std::cout << "\nEcriture OK" << std::flush;
    else
        std::cout << "\nEcriture PAS OK" << std::flush;
}

void Stationnaire::solveAnalytiqueStationnaire(){
    // Calcul de la solution analytique
    double p = 2 * (Ly + Lz);       // et non p = 2 * (Lx + Ly);
    double S = Ly * Lz;
    double a = hc * p / (kappa * S);
    for (int i = 0; i < M+1; i++){
        Texact[i] = Te + (Phi_p * cosh(sqrt(a)*Lx) * cosh(sqrt(a)*(Lx-x[i]))) / 
                            (kappa * sqrt(a)*sinh(sqrt(a)*Lx) * cosh(sqrt(a)*Lx));      
    }
}

void Stationnaire::displayA() const{
    std::cout << "\n--> Matrice A <--";
    for (int i = 0; i < M+1; i++){
        std::cout << "\n";
        for (int j = 0; j < M+1; j++)
            std::cout << std::setw(8) << std::left << A[i][j] << "\t";
    }   
}

void Stationnaire::displayTexact() const{
    std::cout << "\n--> Solution analytique Texact pour le modele stationnaire <--\n";
    for (int i = 0; i < M+1; i++)
        std::cout << std::setw(8) << Texact[i] << "\t";
}

void Stationnaire::displayT() const{
    std::cout << "\n--> Solution numerique T <--\n";
    for (int i = 0; i < M+1; i++)
        std::cout << std::setw(8) << T[i] << "\t";
}

void Stationnaire::displayF() const{
    std::cout << "\n--> Second membre F <--\n";
    for (int i = 0; i < M+1; i++)
        std::cout << std::setw(8) << F[i] << "\t";
}

void Stationnaire::displayPlot() const{
    // Ecriture dans le fichier pour faire le graphe
    std::cout << std::endl << "Ecriture du fichier 'data_stationnaire.csv' en cours ...";
    std::ofstream write_to_file("data_stationnaire.csv");
    if(write_to_file){
        write_to_file << "x" << ", " << "T" << ", " << "Texact\n" ;
        for (int i = 0; i < M+1; i++){
            write_to_file << x[i] << ", " << Texact[i] << ", " << T[i] << "\n";
        }
        write_to_file.close();

        // On trace le graphique a l'aide de gnuplot
        std::cout << std::endl << "Ecriture du fichier csv OK\nAffichage 1D de la solution stationnaire via 'gnuplot'" << std::flush;
        system("gnuplot gnuplot_command_stat.dat");

    } else
        std::cout << std::endl << "Erreur lors de l'ecriture du fichier csv, du coup pas de visulation 1D !" << std::flush;
}

Stationnaire::~Stationnaire(){
    for (int i = 0; i < M+1; i++){
        delete [] A[i];
        delete [] L[i];
        delete [] U[i];
    }
    delete [] A;
    delete [] L;
    delete [] U;
    
    delete [] x;
    delete [] T;
    delete [] F;
    delete [] Texact;
    delete [] x_hat;
    delete [] x_hat_indices;
    delete [] T_hat;
}