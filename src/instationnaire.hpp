#ifndef _INSTAT
#define _INSTAT

#include "stationnaire.hpp"

// Classe traitaint du probleme instationnaire
class Instationnaire: public Stationnaire{
    public:
        // Les constructeurs
        Instationnaire(){}

        Instationnaire(double kappa, double Te, double Phi_p, double hc, double Lx, 
                        double Ly, double Lz, int M, int Mx, int My, int Mz, 
                        double rho, double Cp, double t_final, bool flux_constant, int N);

        // Fonction de resolution specifique au probleme instationnaire
        void solve() override;

        // Fonction d'affichage du probleme en 1D pour le modele instationnaire
        void displayPlot () const override;

        // Destructeur
        ~ Instationnaire() override;

    protected: 
        // Parametres physiques supplementaires du systeme
        double rho;
        double Cp; 

        // Parametres geometriques du systeme

        // Parametres supplementaires de simulation 
        double t_final;
        int N;

        // Pour le deuxieme cas de solve (Phi_p variant)
        bool flux_constant;

        // Differents points en temps(15, 30, etc) ou la solution sera trace
        double ** T_to_plot;
        
        // Differents points en espace(0, M/2, etc) ou la solution sera trace
        double ** X_to_plot;

};



#endif