% Charger les données depuis le fichier (modifier 'data.txt' avec le vrai nom)
data = load('outputs/output.txt');

% Extraire les colonnes
x = data(:, 1);
y = data(:, 2);


% Extraction de la partie réelle de P depuis la 3ème colonne
% On suppose que la colonne 3 est sous la forme (Re, Im), donc on garde uniquement Re
 P_real = sscanf(data(:, 3), '(%f,%f)');  % Si stocké sous format texte "(Re, Im)"

% Créer une grille régulière pour interpolation
xq = linspace(min(x), max(x), 100);  % 100 points en X
yq = linspace(min(y), max(y), 100);  % 100 points en Y
[Xq, Yq] = meshgrid(xq, yq);

% Interpolation des valeurs de P_real sur la grille
Pq = griddata(x, y, P_real, Xq, Yq, 'cubic');

% Affichage de la carte de chaleur
figure;
imagesc(xq, yq, Pq);
colorbar;
xlabel('X');
ylabel('Y');
title('Partie réelle de la pression P');
set(gca, 'YDir', 'normal'); % Correction pour l'orientation de l'image

