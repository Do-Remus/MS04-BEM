% Charger les données depuis le fichier (modifier 'data.txt' avec le vrai nom)
data = load('u_N_plusResu_k=100.000000_a=1.000000_N=150.txt');

% Extraire les colonnes
x = data(:, 1);
y = data(:, 2);
P_real= data(:,3);
P_imag = data(:,4);

% Extraction de la partie réelle de P depuis la 3ème colonne
% On suppose que la colonne 3 est sous la forme (Re, Im), donc on garde uniquement Re
%P_real = sscanf(data(:, 3), '(%f,%f)');  % Si stocké sous format texte "(Re, Im)"

% Créer une grille régulière pour interpolation
xq = linspace(min(x), max(x), 100);  % 100 points en X
yq = linspace(min(y), max(y), 100);  % 100 points en Y
[Xq, Yq] = meshgrid(xq, yq);

% Interpolation des valeurs de P_real sur la grille
Pq = griddata(x, y, P_real, Xq, Yq, 'linear');



% Affichage de la carte de chaleur
figure;
imagesc(xq, yq, Pq);
colorbar;
xlabel('X');
ylabel('Y');
axis equal;
title('Partie réelle de la pression P');
%set(gca, 'YDir', 'normal'); % Correction pour l'orientation de l'image
% Modifier la taille des axes et de la colorbar
set(gca, 'FontSize', 40);
set(colorbar, 'FontSize', 40);

hold on; % Superposition des cercles noirs

% Charger les cercles depuis maillage.txt
maillage_data = load('outputs/obstacles.txt');

% Déterminer si maillage.txt contient les centres et rayons des cercles
if size(maillage_data, 2) == 3
    % Format : (xc, yc, r) avec xc, yc centre et r rayon
    for i = 1:size(maillage_data, 1)
        xc = maillage_data(i, 1);
        yc = maillage_data(i, 2);
        r = maillage_data(i, 3);

        % Tracer et remplir le cercle en noir
        theta = linspace(0, 2*pi, 100); % 100 points pour un cercle lisse
        x_circle = xc + r * cos(theta);
        y_circle = yc + r * sin(theta);
        fill(x_circle, y_circle, 'k', 'EdgeColor', 'none'); % Cercle noir sans bord
    end
else
    disp("Le fichier maillage.txt ne contient pas des cercles sous la forme (xc, yc, r)");
endif




hold off; % Arrêter la superposition



% Interpolation des valeurs de P_imag sur la grille
Pq = griddata(x, y, P_imag, Xq, Yq, 'linear');



% Affichage de la carte de chaleur
figure;
imagesc(xq, yq, Pq);
colorbar;
xlabel('X');
ylabel('Y');
axis equal;
title('Partie imaginaire de la pression P');
%set(gca, 'YDir', 'normal'); % Correction pour l'orientation de l'image
% Modifier la taille des axes et de la colorbar
set(gca, 'FontSize', 40);
set(colorbar, 'FontSize', 40);

hold on; % Superposition des cercles noirs

% Charger les cercles depuis maillage.txt
maillage_data = load('outputs/obstacles.txt');

% Déterminer si maillage.txt contient les centres et rayons des cercles
if size(maillage_data, 2) == 3
    % Format : (xc, yc, r) avec xc, yc centre et r rayon
    for i = 1:size(maillage_data, 1)
        xc = maillage_data(i, 1);
        yc = maillage_data(i, 2);
        r = maillage_data(i, 3);

        % Tracer et remplir le cercle en noir
        theta = linspace(0, 2*pi, 100); % 100 points pour un cercle lisse
        x_circle = xc + r * cos(theta);
        y_circle = yc + r * sin(theta);
        fill(x_circle, y_circle, 'k', 'EdgeColor', 'none'); % Cercle noir sans bord
    end
else
    disp("Le fichier maillage.txt ne contient pas des cercles sous la forme (xc, yc, r)");
endif




hold off; % Arrêter la superposition


% Interpolation des valeurs de norme de P sur la grille
P_norm = sqrt(P_real.^2 + P_imag.^2);
Pq = griddata(x, y, P_norm, Xq, Yq, 'linear');



% Affichage de la carte de chaleur
figure;
imagesc(xq, yq, Pq);
colorbar;
xlabel('X');
ylabel('Y');
axis equal;
title('Norme de la pression complexe P');
%set(gca, 'YDir', 'normal'); % Correction pour l'orientation de l'image
% Modifier la taille des axes et de la colorbar
set(gca, 'FontSize', 40);
set(colorbar, 'FontSize', 40);

hold on; % Superposition des cercles noirs

% Charger les cercles depuis maillage.txt
maillage_data = load('outputs/obstacles.txt');

% Déterminer si maillage.txt contient les centres et rayons des cercles
if size(maillage_data, 2) == 3
    % Format : (xc, yc, r) avec xc, yc centre et r rayon
    for i = 1:size(maillage_data, 1)
        xc = maillage_data(i, 1);
        yc = maillage_data(i, 2);
        r = maillage_data(i, 3);

        % Tracer et remplir le cercle en noir
        theta = linspace(0, 2*pi, 100); % 100 points pour un cercle lisse
        x_circle = xc + r * cos(theta);
        y_circle = yc + r * sin(theta);
        fill(x_circle, y_circle, 'k', 'EdgeColor', 'none'); % Cercle noir sans bord
    end
else
    disp("Le fichier maillage.txt ne contient pas des cercles sous la forme (xc, yc, r)");
endif




hold off; % Arrêter la superposition















