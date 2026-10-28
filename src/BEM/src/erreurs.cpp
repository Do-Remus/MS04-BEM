#include "../erreurs.hpp"

ErreursSolution calculer_erreurs_solution(const std::vector<Point> &pointsSolution, const Vecteur &u_cached)
{
    // Répartition des points entre les processus
    ErreursSolution result;

    const Index nbPoints = pointsSolution.size();
    const Index nbPointsBase = nbPoints / mpi_size;
    const Index reste = nbPoints % mpi_size;
    const Index nbPointsLocaux = nbPointsBase + (mpi_rank < static_cast<int>(reste) ? 1 : 0);
    const Index indiceDebut = mpi_rank * nbPointsBase + std::min(static_cast<Index>(mpi_rank), reste);

    std::vector<Complex> u_analytique_local(nbPointsLocaux, Complex(0.0));
    std::vector<Real> erreurs_local(nbPointsLocaux, 0.0);

    Real sommeErreurLocale = 0.0;
    Real sommeAnalytiqueLocale = 0.0;
    Real erreurMaxLocale = 0.0;

    // Calcul des solutions analytiques et des erreurs locales
    for (Index i = 0; i < nbPointsLocaux; i++)
    {
        const Index j = indiceDebut + i;
        const Point &Pj = pointsSolution[j];
        const Complex uExact = u_N_plus_analytique(Pj, rayon, idxTroncature);

        const Complex difference = u_cached[j] - uExact;

        const Real erreur = erreur_relative(u_cached[j], uExact);

        u_analytique_local[i] = uExact;
        erreurs_local[i] = erreur;

        sommeErreurLocale += std::norm(difference);
        sommeAnalytiqueLocale += std::norm(uExact);

        erreurMaxLocale = std::max(erreurMaxLocale, erreur);
    }

    // Réduction des erreurs L2 et maximale sur le processus 0
    Real sommeErreur = 0.0;
    Real sommeAnalytique = 0.0;
    Real erreurMax = 0.0;

    MPI_Reduce(&sommeErreurLocale, &sommeErreur, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sommeAnalytiqueLocale, &sommeAnalytique, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&erreurMaxLocale, &erreurMax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // Préparation du rassemblement des résultats sur le processus 0
    std::vector<int> nombresLocaux;
    std::vector<int> deplacements;

    if (mpi_rank == 0)
    {
        nombresLocaux.resize(mpi_size);
        deplacements.resize(mpi_size);

        for (int rank = 0; rank < mpi_size; rank++)
        {
            nombresLocaux[rank] = static_cast<int>(nbPointsBase + (rank < static_cast<int>(reste) ? 1 : 0));
            deplacements[rank] = rank * static_cast<int>(nbPointsBase) + std::min(rank, static_cast<int>(reste));
        }

        result.u_analytique.resize(nbPoints, Complex(0.0));
        result.erreurs.resize(nbPoints, 0.0);
    }

    // Rassemblement de u_analytique et des erreurs sur le processus 0
    MPI_Gatherv(
        u_analytique_local.data(),
        static_cast<int>(nbPointsLocaux),
        mpi_complex_type(),
        mpi_rank == 0 ? result.u_analytique.data() : nullptr,
        mpi_rank == 0 ? nombresLocaux.data() : nullptr,
        mpi_rank == 0 ? deplacements.data() : nullptr,
        mpi_complex_type(),
        0,
        MPI_COMM_WORLD);

    MPI_Gatherv(
        erreurs_local.data(),
        static_cast<int>(nbPointsLocaux),
        MPI_DOUBLE,
        mpi_rank == 0 ? result.erreurs.data() : nullptr,
        mpi_rank == 0 ? nombresLocaux.data() : nullptr,
        mpi_rank == 0 ? deplacements.data() : nullptr,
        MPI_DOUBLE,
        0,
        MPI_COMM_WORLD);

    // Calcul et stockage des erreurs globales sur le processus 0
    if (mpi_rank == 0)
    {
        if (sommeAnalytique > 0.0)
        {
            result.erreurL2 = std::sqrt(sommeErreur / sommeAnalytique);
        }
        else
        {
            result.erreurL2 = std::numeric_limits<Real>::quiet_NaN();
        }

        result.erreurMax = erreurMax;
    }

    return result;
}

Real calculer_erreur_relative_p(const Maillage &maillage, const Vecteur &vect_p_cached)
{
    const Index nbSegments = maillage.size();

    Real sommeErreurLocale = 0.0;
    Real sommeAnalytiqueLocale = 0.0;

    for (Index i = mpi_rank; i < nbSegments; i += mpi_size)
    {
        const Complex pExact = p_analytique(maillage[i].milieu, idxTroncature);
        const Complex difference = vect_p_cached[i] - pExact;

        sommeErreurLocale += std::norm(difference);
        sommeAnalytiqueLocale += std::norm(pExact);
    }

    Real sommeErreur = 0.0;
    Real sommeAnalytique = 0.0;

    MPI_Allreduce(&sommeErreurLocale, &sommeErreur, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&sommeAnalytiqueLocale, &sommeAnalytique, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (sommeAnalytique <= 0.0)
        return std::numeric_limits<Real>::quiet_NaN();

    return std::sqrt(sommeErreur / sommeAnalytique);
}

Real calculer_erreur_relative_matrice(const Matrice &A, const Vecteur &vecteur_x, const Vecteur &vecteur_Hx)
{
    const Index N = static_cast<Index>(A.rows());

    const Vecteur vecteur_Ax = A * vecteur_x;
    const Vecteur difference = vecteur_Ax - vecteur_Hx;

    Real sommeErreurLocale = 0.0;
    Real sommeReferenceLocale = 0.0;

    for (Index i = mpi_rank; i < N; i += mpi_size)
    {
        sommeErreurLocale += std::norm(difference[i]);
        sommeReferenceLocale += std::norm(vecteur_Ax[i]);
    }

    Real sommeErreur = 0.0;
    Real sommeReference = 0.0;

    MPI_Allreduce(&sommeErreurLocale, &sommeErreur, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&sommeReferenceLocale, &sommeReference, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return std::sqrt(sommeErreur / std::max(sommeReference, 1e-30));
}
