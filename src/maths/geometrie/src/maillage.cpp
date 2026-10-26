#include "../maillage.hpp"

void Maillage::ajoute_cercle(const Real pas_maillage, const Cercle &Ob)
{
    // Verification pas_maillage > 0
    if (pas_maillage <= 0)
    {
        cout << "ERROR: Le pas du maillage doit-être strictement positif." << endl;
        exit(-1);
    }

    // Construction du cercle
    const Index nb_segments = static_cast<Index>(std::round(2.0 * pi * Ob.rayon / pas_maillage));
    const Real pas_cercle = 2 * pi / nb_segments;

    for (Index i = 0; i < nb_segments; i++)
    {
        Point A(Ob.centre.x + Ob.rayon * cos(pas_cercle * i), Ob.centre.y + Ob.rayon * sin(pas_cercle * i));
        Point B(Ob.centre.x + Ob.rayon * cos(pas_cercle * (i + 1)), Ob.centre.y + Ob.rayon * sin(pas_cercle * (i + 1)));
        this->push_back(Segment(A, B));
    }

    return;
}

void Maillage::export_maillage(const string &filename) const
{
    ofstream f(filename);

    if (!f.is_open())
    {
        cout << "ERROR: Le fichier " << filename << " n'a pas pu être ouvert" << endl;
        exit(-1);
    }

    for (Index i = 0; i < this->size(); i++)
    {
        const Segment &s = this->operator[](i);
        f << s.P1.x << " " << s.P1.y << " ";
        f << s.P2.x << " " << s.P2.y << endl;
    }

    f.close();

    return;
}

vector<Point> Maillage::PointsMaillage() const
{ // only is useful for connex meshings
    vector<Point> points(this->size());

    for (Index i = 0; i < this->size(); i++)
    {
        points[i] = this->operator[](i).P1;
    }

    cout << "out of PointsMaillage 0 " << endl;

    return points;
}

Maillage::Maillage(const vector<Cercle> &cercles, const Real pas_maillage)
{
    // Verification du strict positivité du pas du maillage fait dans ajoute cercle
    for (Index i = 0; i < cercles.size(); i++)
    {
        ajoute_cercle(pas_maillage, cercles[i]);
    }
}

Maillage::Maillage(const vector<Point> &points)
{
    // suppose que aucun segment ne s'interecte
    for (Index i = 0; i < points.size() - 1; i++)
    {
        this->push_back(Segment(points[i], points[i + 1]));
    }

    this->push_back(Segment(points[points.size()], points[0]));
}

/* Fonctions associées à la class */

std::vector<Point> construire_maillage_solution(const Maillage &maillage)
{
    std::vector<Point> pointsSolution;

    const Real xmin = -L / 2.0;
    const Real xmax = L / 2.0;
    const Real ymin = -L / 2.0;
    const Real ymax = L / 2.0;

    // Cercle de reference
    // (utilise uniquement pour la solution analytique du cercle)
    const Real distanceMin = rayon + delta;

    // -- Boite englobante de la frontiere (pour accelerer les tests) --

    Real bxmin = std::numeric_limits<Real>::max();
    Real bxmax = -bxmin;
    Real bymin = bxmin;
    Real bymax = -bxmin;

    for (Index i = 0; i < maillage.size(); i++)
    {
        bxmin = std::min<Real>(bxmin, maillage[i].P1.x);
        bxmax = std::max<Real>(bxmax, maillage[i].P1.x);
        bymin = std::min<Real>(bymin, maillage[i].P1.y);
        bymax = std::max<Real>(bymax, maillage[i].P1.y);
    }

    // -- Point exterieur a l'obstacle et a au moins delta de la frontiere --

    auto exterieurValide = [&](const Point &M) -> bool
    {
        // Point suffisamment loin de la boite englobante :
        // il est nécessairement extérieur à l'obstacle.
        if (M.x < bxmin - delta || M.x > bxmax + delta || M.y < bymin - delta || M.y > bymax + delta)
        {
            return true;
        }

        bool dedans = false;
        Real d2min = std::numeric_limits<Real>::max();

        // Décalage du rayon de test pour éviter les sommets du polygone.
        const Real epsilon = 1e-10;
        const Real yTest = M.y + epsilon;

        for (Index i = 0; i < maillage.size(); i++)
        {
            const Point &A = maillage[i].P1;
            const Point &B = maillage[i].P2;

            // Test d'intersection du rayon horizontal avec le segment.
            if ((A.y > yTest) != (B.y > yTest))
            {
                const Real xIntersection = A.x + (yTest - A.y) * (B.x - A.x) / (B.y - A.y);

                if (M.x < xIntersection)
                    dedans = !dedans;
            }

            // Distance au segment.
            const Real ex = B.x - A.x;
            const Real ey = B.y - A.y;
            const Real len2 = ex * ex + ey * ey;

            Real t = 0.0;

            if (len2 > 0.0)
            {
                t = std::max<Real>(0.0, std::min<Real>(1.0, ((M.x - A.x) * ex + (M.y - A.y) * ey) / len2));
            }

            const Real dx = M.x - (A.x + t * ex);
            const Real dy = M.y - (A.y + t * ey);

            d2min = std::min(d2min, dx * dx + dy * dy);
        }

        return !dedans && d2min >= delta * delta;
    };
    // -- 1) Grille reguliere, filtree par la vraie forme de l'obstacle --

    for (Real x = xmin; x <= xmax; x += pasSolution)
    {
        for (Real y = ymin; y <= ymax; y += pasSolution)
        {
            const Point M(x, y);

            if (exterieurValide(M))
                pointsSolution.emplace_back(M);
        }
    }

    // -- 2) Couches de points le long de la frontiere --

    // Tous les pasProche (abscisse curviligne), on place des points
    // sur la normale sortante a des distances multiples de pasMaillage.
    const Real pasProche = 0.02;

    const std::vector<Real> multiplesProches = {2.0, 4.0, 8.0, 16.0, 32.0};

    Real arc = pasProche;

    for (Index i = 0; i < maillage.size(); i++)
    {
        const Real ex = maillage[i].P2.x - maillage[i].P1.x;
        const Real ey = maillage[i].P2.y - maillage[i].P1.y;
        const Real len = std::sqrt(ex * ex + ey * ey);

        arc += len;

        if (arc < pasProche || len <= 0.0)
            continue;

        arc = 0.0;

        // Normale sortante (contour anti-horaire).
        const Real nx = ey / len;
        const Real ny = -ex / len;

        const Point &C = maillage[i].milieu;

        for (const Real m : multiplesProches)
        {
            const Real d = std::max<Real>(delta, m * pasMaillage);
            const Point M(C.x + d * nx, C.y + d * ny);

            if (M.x < xmin || M.x > xmax || M.y < ymin || M.y > ymax)
            {
                continue;
            }

            if (exterieurValide(M))
                pointsSolution.emplace_back(M);
        }
    }
    return pointsSolution;
}

Arbre construire_arbre_maillage(const Maillage &maillage, const Index N_leafs)
{
    Arbre arbre = std::make_unique<Noeud>();
    // -- Initialisation du 1er noeud --

    arbre->niveau = 0;

    const std::vector<Point> &points = maillage.PointsMaillage();
    arbre->indices.resize(points.size());

    arbre->boite.x_min = std::numeric_limits<double>::infinity();
    arbre->boite.x_max = -std::numeric_limits<double>::infinity();

    arbre->boite.y_min = std::numeric_limits<double>::infinity();
    arbre->boite.y_max = -std::numeric_limits<double>::infinity();

    for (Index i = 0; i < points.size(); ++i)
    {
        arbre->indices[i] = i;

        const Point &point = points[i];

        arbre->boite.x_min = std::min(arbre->boite.x_min, point.x);
        arbre->boite.x_max = std::max(arbre->boite.x_max, point.x);

        arbre->boite.y_min = std::min(arbre->boite.y_min, point.y);
        arbre->boite.y_max = std::max(arbre->boite.y_max, point.y);
    }

    arbre->boite.calculer_diametre();

    // -- Construction récursive de l'arbre --
    return arbre;
}


void subdivision_arbre_recursive( Noeud &noeud, const Index N_leafs, std::vector<Point> &points){

    if(noeud.indices.size() > N_leafs){
        const double largeur = (noeud.boite.x_max - noeud.boite.x_min);
        const double hauteur = (noeud.boite.y_max - noeud.boite.y_min);

        //separation geometrique en 2 groupe divisée selon la plus longue direction
        //trouver quelle est la direction choisie (x ou y) puis definir nouvelles boites englobantes definies comme les moities de la boite englobante du parent
        const bool couper_x  = (largeur >= hauteur);
        const double milieu = couper_x ? (noeud.boite.x_min +noeud.boite.x_max)/2. : (noeud.boite.y_min +noeud.boite.y_max)/2.;
        noeud.fils_droit = std::make_unique<Noeud>();
        noeud.fils_gauche = std::make_unique<Noeud>();
        noeud.fils_gauche->niveau = noeud.niveau +1;
        noeud.fils_droit->niveau = noeud.niveau +1;

        //fils droit
        noeud.fils_droit->boite.x_min = noeud.boite.x_min;
        noeud.fils_droit->boite.y_min = noeud.boite.y_min;
        noeud.fils_droit->boite.x_max = couper_x ? milieu : noeud.boite.x_max;
        noeud.fils_droit->boite.y_max = couper_x ? noeud.boite.y_max : milieu;

        //fils gauche
        noeud.fils_gauche->boite.x_min = couper_x ? milieu : noeud.boite.x_min;
        noeud.fils_gauche->boite.y_min = couper_x ? noeud.boite.y_min : milieu;
        noeud.fils_gauche->boite.x_max = noeud.boite.x_max;
        noeud.fils_gauche->boite.x_max = noeud.boite.x_max;

        //iteration sur les points indexés dans le parent pour séparation sur les 2 fils
        for(Index i=0; i<noeud.indices.size(); i++){
            const double coord = couper_x ? points[noeud.indices[i]].x : points[noeud.indices[i]].y;
            if(coord <= milieu){
                noeud.fils_droit->indices.push_back(noeud.indices[i]);
            }else{
                noeud.fils_gauche->indices.push_back(noeud.indices[i]);
            }
        }
        

    }else{
        return;
    }
    subdivision_arbre_recursive(*noeud.fils_gauche, N_leafs, points);
    subdivision_arbre_recursive(*noeud.fils_droit, N_leafs, points);

}

double distance(const BoiteEnglobante &boite_a, const BoiteEnglobante &boite_b)
{
    double distance_x = std::max({0.0, boite_a.x_min - boite_b.x_max, boite_b.x_min - boite_b.x_max});
    double distance_y = std::max({0.0, boite_a.y_min - boite_b.y_max, boite_b.y_min - boite_b.y_max});
    return std::sqrt(distance_x * distance_x + distance_y * distance_y);
}

bool sont_admissibles(const BoiteEnglobante &boite_a, const BoiteEnglobante &boite_b, const double eta)
{
    return std::min(boite_a.diametre, boite_b.diametre) <= eta * distance(boite_a, boite_b);
}

ArbreBlocs construire_bloc_recursif(const Noeud &tau, const Noeud &sigma, const double eta)
{
    // Création du bloc associé aux deux clusters.
    auto bloc = std::make_unique<Bloc>();

    bloc->cluster_lignes = &tau;
    bloc->cluster_colonnes = &sigma;

    // 1. Le bloc est admissible : c'est une feuille.
    if (sont_admissibles(tau.boite, sigma.boite, eta))
    {
        bloc->type = TypeBloc::Admissible;
        return bloc;
    }

    // 2. Le bloc n'est pas admissible et ne peut plus
    // être subdivisé si l'un des clusters est une feuille.
    if (tau.est_feuille() || sigma.est_feuille())
    {
        bloc->type = TypeBloc::Dense;
        return bloc;
    }

    // 3. Le bloc n'est pas admissible :
    // on le subdivise en quatre blocs fils.
    bloc->type = TypeBloc::Interne;

    bloc->fils[0][0] = construire_bloc_recursif(*tau.fils_gauche, *sigma.fils_gauche, eta);
    bloc->fils[0][1] = construire_bloc_recursif(*tau.fils_gauche, *sigma.fils_droit, eta);
    bloc->fils[1][0] = construire_bloc_recursif(*tau.fils_droit, *sigma.fils_gauche, eta);
    bloc->fils[1][1] = construire_bloc_recursif(*tau.fils_droit, *sigma.fils_droit, eta);

    return bloc;
}

ArbreBlocs construire_arbre_blocs(const Arbre &arbre_clusters, const double eta)
{
    if (!arbre_clusters)
    {
        return nullptr;
    }

    return construire_bloc_recursif(*arbre_clusters, *arbre_clusters, eta);
}
