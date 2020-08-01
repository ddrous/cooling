#ifndef STAT
#define STAT

// // Nouvelle classe matrice pour traiter les matrices qui apparaissent
// class Matrice{
//     public:
//         Matrice(){
//             MAT = new double *[M+1];
//             for (int i = 0; i < M+1; i++){
//                 MAT[i] = new double [M+1];
//             }
//         }
//         int M;
//         double ** MAT;        
// };

// Classe traitaint du probleme stationnaire
class Stationnaire{
    public:
        // Constructeur vide
        Stationnaire(){}

        // Constructeur qui initialise le probleme stationnnaire
        Stationnaire(double kappa, double Te, double Phi_p, 
                            double hc, double Lx, double Ly, double Lz, int M, 
                            int Mx, int My, int Mz);

        // Fonction qui va resoudre le probleme
        virtual void solve();

        // Fonction specifiques pour la facto LU
        void factoLU();

        // Fonction specifiques pour la verification de la facto LU
        void verificationFactoLU() const;

        // Fonction specifiques pour la resolution lineaire 
        void descenteRemontee();

        // Fonction qui recherche les intervalles pour l'interpolation
        void searchInterval();

        // Fonction qui va calculer l'interpolant lineaire T_hat
        void interpolationLineaire();

        // Fonction pour ecrire les resultats 3d
        bool ecriture3D(int indice_temps, bool flux_de_chaleur_constant);      // indice temps correspond au temps auquel on ecrit

        // Fonction specifiques pour la verification de la solution
        void verificationSolution() const;

        // Calcule la solution analytique T_exact dans le cas du modele stataionnaire
        void solveAnalytiqueStationnaire();

        // Fonction d'affichage du probleme en 1d
        virtual void displayPlot() const;

        // Fonction pour afficher A
        void displayA () const;

        // Fonction pour afficher T_exact
        void displayTexact() const;

        // Fonction pour afficher la solution numerique T
        void displayT() const;
        
        // Fonction pour afficher le second membre F
        void displayF() const;

        // Destructeur
        virtual ~ Stationnaire();

    protected:      
        // Parametres physiques du systeme
        double kappa;
        double Te;
        double Phi_p;
        double hc;

        // Parametres geometriques du systeme
        double Lx, Ly, Lz;

        // Parametres de simulation
        int M;
        double * x;    // Points de discretisation en espace

        double ** A;
        double ** L;
        double ** U; 
        double * T;     // T = X pour le modele stationnaire
        double * F;

        // Solution analytique
        double * Texact;

        // Parametres de visualisation 3D
        int Mx, My, Mz;
        // x_hat contient contient deux lignes de coordonees, 
        double * x_hat;
        // x_hat_indices; contient les indices des differents x_hat
        int * x_hat_indices;
        // Temperature 3D T_hat
        double * T_hat;

};

#endif