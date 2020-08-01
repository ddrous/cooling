#include <iostream>
#include <fstream>
#include <map>
#include "stationnaire.hpp"
#include "instationnaire.hpp"

int main(int argc, char * argv[]){

    std::cout << "\n==> Simulation thermique d'une ailette d'un dissipateur de chaleur <==" << std::endl;

    std::map<std::string, int> noms_des_variables;
    std::map<std::string, int> :: iterator it;

    noms_des_variables["Lx"] = 0;
    noms_des_variables["Ly"] = 1;
    noms_des_variables["Lz"] = 2;
    noms_des_variables["M"] = 3;
    noms_des_variables["Phi"] = 4;
    noms_des_variables["hc"] = 5;
    noms_des_variables["Te"] = 6;
    noms_des_variables["stationary"] = 7;
    noms_des_variables["TFinal"] = 8;
    noms_des_variables["N"] = 9;
    noms_des_variables["Mx"] = 10;
    noms_des_variables["My"] = 11;
    noms_des_variables["Mz"] = 12;
    noms_des_variables["rho"] = 13;
    noms_des_variables["Cp"] = 14;
    noms_des_variables["kappa"] = 15;
    noms_des_variables["flux_constant"] = 16;

    double * valeurs_des_variables = new double[17];

    if (argc == 2){
        std::string nom_de_la_config = argv[1];
        std::ifstream config_file("config_files/" + nom_de_la_config);

        bool read_success = true;
        std::string nom_de_la_variable;
        if(config_file){
            // std::getline(config_file >> std::ws, nom_de_la_variable, '\n');
            std::cout <<nom_de_la_variable << "\nLecture des parametres du systeme en cours ... ";
            for (int i = 0; i < 17; i++){
                config_file >> nom_de_la_variable;
                it = noms_des_variables.find(nom_de_la_variable);
                if (it != noms_des_variables.end()){
                    config_file >> valeurs_des_variables[it->second];
                    std::cout << "\n" << nom_de_la_variable << " = " << valeurs_des_variables[it->second] << "";
                }
                else
                    read_success = false;
            }
            
            config_file.close();
        
        }else{
            std::cout << "\nErreur lors de l'ouverture du fichier de configuration !\n" << std::endl;
            exit(-1);
        }

        if (read_success == true){
            std::cout << "\nLecture des 17 parametres physiques et geometriques OK!";
        } else{
            std::cout << "\nERREUR: Fichier config mal ecrit, mauvaise lecture !\n" << std::endl;
            exit(-1);        
        }

    }else{
        std::cout << "\nERREUR: Fournissez un (et un seul) fichier de configuration !\n" << std::endl;
        exit(-1);
    }

    // Definition des variables du probleme dans leur ordre 
    // Specifique au le probleme stationnaire
    double kappa = valeurs_des_variables[15];
    double Te = valeurs_des_variables[6];
    double Phi_p = valeurs_des_variables[4];
    double hc = valeurs_des_variables[5];            // Ventilateur etteint ou allume
    double Lx = valeurs_des_variables[0];           // Les longuers sont lues en m
    double Ly = valeurs_des_variables[1];
    double Lz = valeurs_des_variables[2];
    int M = (int)valeurs_des_variables[3];
    // Necessaires pour le probleme instationnaire
    double rho = valeurs_des_variables[13];
    double Cp = valeurs_des_variables[14];
    double t_final = valeurs_des_variables[8];
    bool flux_constant = (bool)valeurs_des_variables[16];
    bool stationary = (bool)valeurs_des_variables[7];
    int N = (int)valeurs_des_variables[9];
    // Insdispensables pour la simualtion avec paraview
    int Mx = (int)valeurs_des_variables[10];
    int My = (int)valeurs_des_variables[11];
    int Mz = (int)valeurs_des_variables[12];

    // Partie de resolution du probleme
    if (stationary != 0){        
        Stationnaire * probleme = new Stationnaire(kappa, Te, Phi_p, hc, Lx, Ly, Lz, M, Mx, My, Mz);
        probleme->solve();
        // probleme->displayA();
        probleme->solveAnalytiqueStationnaire();
        // probleme->displayTexact();
        // probleme->displayT();
        probleme->displayPlot();
        delete probleme;
    }

    else{
        Stationnaire * probleme = new Instationnaire(kappa, Te, Phi_p, hc, Lx, Ly, Lz, M, Mx, My, Mz, rho, Cp, t_final, flux_constant, N);
        probleme->solve();
        probleme->solveAnalytiqueStationnaire();
        probleme->displayPlot();
        delete probleme;
    }

    delete [] valeurs_des_variables;

    std::cout << std::endl << "==> Roussel D. Nzoyem N. <==" << std::endl << std::endl;
    return 0;
}