# Documentation du projet

## Classes

### Point (point.hpp)

Classe Point contenant des coordonnées, x et y, sous forme de double.

Constructeurs:
- Point() : créer le point (0,0)
- Point(double x, double y) : créer le point (x,y)

Public:
- x : double
- y : double

## Cercle (cercle.hpp)

Classe Cercle contenant un centre (Point) et un rayon (double).

Constructeurs:
- Cercle() : créer le cercle de centre (0,0) et rayon 1
- Cercle(Point A, rayon r) : créer le cercle de centre A et de rayon r (une erreur sera levé si ce dernier n'est pas strictement positif)

Public:
- centre : Point
- rayon : double

## Segment (segment.hpp)

Classe Segment contenant les deux bouts du segment sous forme de points.

Constructeurs:
- Segment() : créer le segment du point nul au point nul
- Segment(Point A, Point B) : créer le segment du point A au point B

Public:
- P1 : Point
- P2 : Point

## Maillage (maillage.hpp)

Classe Maillage liste de segments représentant un maillage.

Constructeurs:
- Maillage() : créer un maillage contenant 0 Segments
- Maillage(const vector<Cercle> cercles, const double pas_maillage) : créer le maillage approximant les cercles d'une liste de cercles donnée pour un pas de maillage donné

Public:
- cf. class vector<Segment>
- ajoute_cercle(const double pas_mailage, const Cercle &Ob) : ajoute les segments approximant un cercle pour un pas de maillage donné
