% Load the data file (change 'data.txt' to your actual filename)
data = load('outputs/test.txt');

% Extract columns
x1 = data(:,1);
y1 = data(:,2);
x2 = data(:,3);
y2 = data(:,4);



% Plot segments
hold on; % Keep all segments on the same figure
for i = 1:length(x1)
    plot([x1(i), x2(i)], [y1(i), y2(i)], '-b', 'LineWidth', 2); % Blue lines
end
for(i =1: length(x1)/100)
    plot([x1(i), y2(i)-y1(i)], [y1(1), x1(i)-x2(i)], '-r', 'LineWidth', 2)
end
%Point normale(AB.P2.y - AB.P1.y, AB.P1.x - AB.P2.x);
hold off;

% Formatting
xlabel("X-axis");
ylabel("Y-axis");
title("Line segments between points");
grid on;
axis equal; % Ensure correct aspect ratio





