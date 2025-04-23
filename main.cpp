#include <filesystem>
#include <chrono>
#include <fstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <ctime>
#include <nlohmann/json.hpp>
#include <set>
#include <iomanip>
#include <numeric>
#include <cmath>
using namespace std;
using json = nlohmann::json;

using Chromosome = vector<int>;// Un chromosome est un vecteur d'entiers représentant les positions des capteurs
using Population = vector<Chromosome>;// Une population est un vecteur de chromosomes

int N_CAPTEURS = 0;
int N_EMPLACEMENTS = 0;
int TAILLE_POPULATION = 0;
int TAILLE_CHROMOSOME = N_CAPTEURS * N_EMPLACEMENTS ;
int MAX_GENERATIONS = 10;
// Structures de données pour les emplacements et points d'intérêt
vector<pair<double, double>> emplacements; // Liste des positions possibles (x,y)
vector<pair<double, double>> points_interet;// Liste des points à couvrir (x,y)
vector<vector<double>> matrice_distance; // Matrice des distances entre emplacements et POI
vector<double> rayons_capteurs; // Rayons de couverture de chaque capteur

void lireDonneesDepuisJSON(const string& nomFichier) {    // Lecture du fichier JSON contenant la configuration et les données
    ifstream fichier(nomFichier);
    if (!fichier) {
        cerr << "Erreur ouverture du fichier JSON.\n";
        exit(1);
    }
    json data;
    fichier >> data;// Lecture des données JSON
        // Initialisation des constantes globales
    N_CAPTEURS = data["capteurs"].size();
    N_EMPLACEMENTS = data["emplacements"].size();
    TAILLE_POPULATION = data["genetic_params"]["population_size"];
    TAILLE_CHROMOSOME = N_CAPTEURS * N_EMPLACEMENTS;
    // Chargement des emplacements et points d'intérêt
    for (const auto& e : data["emplacements"])
        emplacements.emplace_back(e["x"], e["y"]);

    for (const auto& p : data["points_interet"])
        points_interet.emplace_back(p["x"], p["y"]);
    // Chargement des rayons des capteursx
    for (const auto& capteur : data["capteurs"])
        rayons_capteurs.push_back(capteur["rayon"]);
}

double calculerDistance(const pair<double, double>& a, const pair<double, double>& b) {
    return sqrt(pow(a.first - b.first, 2) + pow(a.second - b.second, 2));    // Calcul de la distance euclidienne entre deux points

}
//matrice de distance 
void calculerMatriceDistance() {    // Création de la matrice des distances entre tous les emplacements et POI
    matrice_distance.resize(N_EMPLACEMENTS, vector<double>(points_interet.size()));
    for (int i = 0; i < N_EMPLACEMENTS; ++i) {
        for (int j = 0; j < points_interet.size(); ++j)
            matrice_distance[i][j] = calculerDistance(emplacements[i], points_interet[j]);
    }
}

double fonctionFitness(const Chromosome& individu) {    // Évaluation de la qualité d'une solution
    vector<bool> couverts(points_interet.size(), false);// Tableau pour suivre les POI couverts
    int capteurs_utiles = 0;

    for (int i = 0; i < N_CAPTEURS; ++i) {
        if (individu[i] == 0) continue; // Capteur non utilisé   |// Si le capteur n'est pas utilisé

        int e = individu[i] - i * N_EMPLACEMENTS - 1;// Calculer l'indice de l'emplacement
        double rayon = rayons_capteurs[i];// Rayon du capteur courant
        bool a_couvert = false;// Indique si le capteur couvre au moins un POI
     // Vérification de la couverture des points d'intérêt
        for (int j = 0; j < points_interet.size(); ++j) {
            if (!couverts[j] && matrice_distance[e][j] <= rayon) {
                couverts[j] = true;
                a_couvert = true;
            }
        }

        if (a_couvert) capteurs_utiles++;
    }
    // Calcul du score final
    int nb_couverts = count(couverts.begin(), couverts.end(), true);
    double penalite = 0.8 * capteurs_utiles;  // Pénalité pour l'utilisation de capteurss
    return static_cast<double>(nb_couverts) - penalite;
}//Cette fonction de fitness équilibre deux objectifs contradictoires : maximiser la couverture des points d'intérêt tout en minimisant le nombre de capteurs utilisés

void sauvegarderMeilleur(const Chromosome& individu, double fitness) {
    json j;
    j["fitness"] = fitness;
    j["capteurs"] = individu;
    ofstream out("best_result.json");
    out << setw(4) << j << endl;
}

Chromosome tournoiSelection(const Population& population, int k = 5) {
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distrib(0, population.size() - 1);
    Chromosome meilleur_individu;
    double meilleure_fitness = -1e9;

    for (int i = 0; i < k; ++i) {
        int index = distrib(gen);
        const Chromosome& candidat = population[index];
        double f = fonctionFitness(candidat);
        if (f > meilleure_fitness) {
            meilleure_fitness = f;
            meilleur_individu = candidat;
        }
    }
    return meilleur_individu;
}

pair<Chromosome, Chromosome> croisement(const Chromosome& parent1, const Chromosome& parent2) {
    int taille = parent1.size();
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(1, taille - 2);
    int point = dist(gen);

    Chromosome enfant1 = parent1;
    Chromosome enfant2 = parent2;

    for (int i = point; i < taille; ++i)
        swap(enfant1[i], enfant2[i]);

    return {enfant1, enfant2};
}

void mutation(Chromosome& individu, double taux_mutation = 0.1) {
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> proba(0.0, 1.0);

    for (int i = 0; i < individu.size(); ++i) {
        if (proba(gen) < taux_mutation) {
            int base = i * N_EMPLACEMENTS;
            uniform_int_distribution<> dist_emplacement(base + 1, base + N_EMPLACEMENTS); // Correction ici
            individu[i] = dist_emplacement(gen);
        }
    }
}
Population initialiserPopulation() {
    Population population;
    random_device rd;
    mt19937 gen(rd()); // Déclaration ici

    for (int i = 0; i < TAILLE_POPULATION; ++i) {
        Chromosome individu;
        for (int j = 0; j < N_CAPTEURS; ++j) {
            uniform_real_distribution<> chance(0.0, 1.0); // Pour décider de l'activation

            if (chance(gen) < 0.3) {
                individu.push_back(0); // Capteur non utilisé
            } else {
                int base = j * N_EMPLACEMENTS;
                uniform_int_distribution<> distrib(base + 1, base + N_EMPLACEMENTS);
                individu.push_back(distrib(gen));
            }
        }
        population.push_back(individu);
    }

    return population; // Ajouté pour corriger le warning
}
// ... (ton code précédent : fonctions inchangées jusqu'au main)

// CSV : écrire l'entête
void initialiserFichierResultats(const string& nomFichierCSV) {
    ofstream fichier(nomFichierCSV);
    fichier << "Fichier,POI,Capteurs,Faisable,Capteurs_utilisés,Temps_ms,Generations\n";
    fichier.close();
}

// CSV : enregistrer les résultats d’une expérience
void enregistrerResultat(const string& fichier_json, int nb_poi, int nb_capteurs, bool faisable,
                         int capteurs_utilises, long long temps_ms, int generation, const string& fichier_csv) {
    ofstream fichier(fichier_csv, ios::app);
    fichier << fichier_json << "," << nb_poi << "," << nb_capteurs << ","
            << (faisable ? "Oui" : "Non") << "," << capteurs_utilises << ","
            << temps_ms << "," << generation << "\n";
    fichier.close();
}

// Point d'entrée principal
int main() {
    const string dossier = "experiences/";
    const string fichier_csv = "resultats1.csv";
if (!filesystem::exists(fichier_csv)) {
    initialiserFichierResultats(fichier_csv);
}


    for (const auto& entry : filesystem::directory_iterator(dossier)) {
        string nom_fichier = entry.path().string();

        // Début chrono
        auto debut = chrono::high_resolution_clock::now();

        // Lecture et init
        emplacements.clear();
        points_interet.clear();
        rayons_capteurs.clear();
        matrice_distance.clear();

        lireDonneesDepuisJSON(nom_fichier);
        
        calculerMatriceDistance();
        Population population = initialiserPopulation();

        Chromosome meilleur_global;
        double meilleur_fitness_global = -1e9;
        int generation_finale = -1;

        bool faisable = false;

        for (int generation = 1; generation <= MAX_GENERATIONS; ++generation) {
            double meilleure_fitness = -1e9;
            Chromosome meilleur_local;

            for (const Chromosome& ind : population) {
                vector<bool> couverts(points_interet.size(), false);
                for (int k = 0; k < N_CAPTEURS; ++k) {
                    if (ind[k] == 0) continue;
                    int e = ind[k] - k * N_EMPLACEMENTS - 1;
                    if (e < 0 || e >= matrice_distance.size()) continue;
                    double rayon = rayons_capteurs[k];
                    for (int j = 0; j < points_interet.size(); ++j)
                        if (!couverts[j] && matrice_distance[e][j] <= rayon)
                            couverts[j] = true;
                }
                int nb_couverts = count(couverts.begin(), couverts.end(), true);
                bool solution_faisable = (nb_couverts == points_interet.size());

                double fitness = fonctionFitness(ind);
                if (solution_faisable && fitness > meilleure_fitness) {
                    meilleure_fitness = fitness;
                    meilleur_local = ind;
                    generation_finale = generation;
                    faisable = true;
                }
            }

            if (faisable) {
                meilleur_global = meilleur_local;
                break;
            }

            // Sélection, croisement, mutation
            Population nouvelle_gen;
            for (int i = 0; i < TAILLE_POPULATION; ++i)
                nouvelle_gen.push_back(tournoiSelection(population));

            Population croisee;
            for (int i = 0; i < TAILLE_POPULATION; i += 2) {
                auto [e1, e2] = croisement(nouvelle_gen[i], nouvelle_gen[(i + 1) % TAILLE_POPULATION]);
                croisee.push_back(e1);
                croisee.push_back(e2);
            }
            for (Chromosome& ind : croisee)
                mutation(ind);

            population = croisee;
        }

        // Fin chrono
        auto fin = chrono::high_resolution_clock::now();
        auto temps_ms = chrono::duration_cast<chrono::milliseconds>(fin - debut).count();

        int capteurs_utiles = 0;
        if (faisable) {
            for (int i = 0; i < N_CAPTEURS; ++i) {
                if (meilleur_global[i] != 0) capteurs_utiles++;
            }
        }

        enregistrerResultat(
            nom_fichier.substr(nom_fichier.find_last_of("/\\") + 1),
            points_interet.size(),
            N_CAPTEURS,
            faisable,
            capteurs_utiles,
            temps_ms,
            generation_finale == -1 ? MAX_GENERATIONS : generation_finale,
            fichier_csv
        );
    }

    cout << "\n✅ Toutes les expériences ont été traitées.\n";
    return 0;
}
