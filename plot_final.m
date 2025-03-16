% Charger les données de output.txt
output_data = load('outputs/output.txt');

% Extraire les colonnes
x = output_data(:, 1);
y = output_data(:, 2);
P_real = output_data(:, 3);  % Partie réelle de P

% Création d'une grille régulière
xq = linspace(min(x), max(x), 100);
yq = linspace(min(y), max(y), 100);
[Xq, Yq] = meshgrid(xq, yq);

% Interpolation des valeurs de P_real sur la grille
Pq = griddata(x, y, P_real, Xq, Yq, 'cubic');

% Affichage de la heatmap
figure;
imagesc(xq, yq, Pq);
set(gca, 'YDir', 'normal'); % Correction orientation
colormap(jet); % Palette de couleurs
colorbar;
xlabel('X');
ylabel('Y');
title('Partie réelle de P(x, y)');

hold on; % Superposition des cercles noirs

% Charger les segments de maillage.txt
maillage_data = load('outputs/maillage.txt');

% Tracer les segments en noir
for i = 1:size(maillage_data, 1)
    x1 = maillage_data(i, 1);
    y1 = maillage_data(i, 2);
    x2 = maillage_data(i, 3);
    y2 = maillage_data(i, 4);
    plot([x1, x2], [y1, y2], 'k', 'LineWidth', 2); % Segment noir épais
end

hold off; % Arrêter la superposition

