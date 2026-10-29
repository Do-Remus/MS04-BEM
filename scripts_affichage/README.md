# Visualisation MS04 / BEM

## Fichiers

- `main.py` : propose les calculs disponibles, puis les visualisations à générer.
- `io_utils.py` : lecture des fichiers et alias de nomenclature.
- `visualiser_solution.py` : champs exacts/calculés, erreur et animation GIF.
- `visualiser_clusters.py` : animation GIF et image PNG de l'arbre de clusters.
- `visualiser_blocs.py` : animation GIF et image PNG de l'arbre des blocs.

## Installation

Depuis le dossier du projet :

```bash
python -m pip install -r scripts_affichage/requirements.txt
```

Adaptez le chemin si vous avez placé les scripts ailleurs.

## Lancement

Depuis la racine du projet :

```bash
python scripts_affichage/main.py
```

Le script cherche par défaut dans `outputs/` des sous-dossiers qui contiennent
`solution.txt` ou un fichier `u_*.txt`. Il conserve également une compatibilité
de base avec l'ancien format où les fichiers se trouvaient directement dans
`outputs/`.

Les figures sont enregistrées dans :

```text
outputs/<dossier_du_calcul>/visualisations/
```

## Nomenclature

Les noms reconnus sont définis dans `io_utils.py`, dans `ALIAS_FICHIERS`.
Pour ajouter un nouveau nom, ajoutez-le à la liste correspondant à la clé.
Exemple :

```python
"solution": ["solution.txt", "u_*.txt", "solution_finale.txt"]
```

Les fichiers de solution doivent contenir au moins sept colonnes :

`x y uExact_re uExact_im uCached_re uCached_im erreur`

Les fichiers de maillage de frontière contiennent quatre colonnes :

`x1 y1 x2 y2`

Le maillage des points de solution contient `x y`, avec éventuellement le nombre
de points sur la première ligne.

## À vérifier côté C++

Pour que le script détecte plusieurs calculs, le code C++ doit créer un sous-dossier
par calcul, par exemple :

```text
outputs/cercle_k10_R1_.../
```

et exporter `solution.txt`, `maillage.txt` et les CSV dans ce même dossier.
Si `output_folder` ne contient que `cercle_k10_...`, ajoutez le préfixe `outputs/`
dans le chemin de sortie ou construisez le chemin avec `std::filesystem::path`.

Les noms de colonnes des CSV doivent correspondre aux noms utilisés par les scripts.
