%% 1. Load and get the bands near Fermi surface
% Define the data path for loading energy bands
path = "replydata/Vmax_29/Lm_114/";

% Load energy band data from the specified path
load(path + "enk.mat");

% Extract dimensions of the Hamiltonian
dimH = size(Enk, 3);
valband = ceil(dimH / 2); % Index of the valence band

% Calculate Fermi energy (Ef) based on the middle bands
Ef = (min(Enk(:, :, valband), [], 'all') + max(Enk(:, :, valband + 1), [], 'all')) / 2;

% Select 20 bands near the charge neutrality point (CNP) and shift by Ef
Enk1 = Enk(:, :, valband - 1 - 9:valband - 1 + 10) - Ef;

%% 2. Define the Brillouin Zone (BZ)
% Set the center of the Brillouin Zone
center = [0, 0];

% Create a hexagonal BZ
BZ = nsidedpoly(6, 'Center', center, 'SideLength', norm(Gm1) / sqrt(3));
BZ = rotate(BZ, 30); % Rotate BZ by 30 degrees
BZ = translate(BZ, 0, norm(Gm1) * 2 / sqrt(3)); % Translate BZ vertically

% Get boundary vertices of the BZ
[verticesX, verticesY] = boundary(BZ);

%% 3. Delete the data outside the BZ
% Remove energy data points outside the BZ
knum = size(Kx, 1);
for i = 1:knum
    for j = 1:knum
        inside = inpolygon(Kx(i, j), Ky(i, j), verticesX, verticesY); % Check if inside BZ
        if inside == 0
            Enk1(i, j, :) = NaN; % Assign NaN to data outside BZ
        end
    end
end

%% 4. Plot the energy bands
clear gca;
figure('Color', 'white');

% Select bands to plot
band1 = 9; 
band2 = 14;
part = 1:1:knum;

% Plot selected bands in 3D
for j = band1:band2
    s = surf(Kx(part, part), Ky(part, part), 1000 .* Enk1(part, part, j), 'FaceAlpha', 1);
    s.EdgeColor = 'none'; 
    hold on;
end

% Determine energy range for visualization
Emin = min(Enk(:, :, band1), [], 'all');
Emax = max(Enk(:, :, band2), [], 'all');
dz1 = 0; dz2 = 0; Lz = 200;

% Plot the Brillouin Zone edges
hexplot(BZ, -Lz - dz1, Lz - dz2);
hold on;

% Adjust view and lighting for better visualization
daspect([1, 1, 4000]);
lightangle(90, 5);
material dull;
camlight('headlight');
view([90, 5]);

% Customize axis and labels
ax = gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
grid off;
set(ax, 'FontSize', 18, 'FontName', 'Arial', 'LineWidth', 1);

% Annotate vertices of the Brillouin Zone
shift1 = 0.003; shift2 = 0.002;
text(verticesX(1) + shift1, verticesY(1) + shift2, -Lz + 8, "Y", 'FontSize', 18, 'FontName', 'Arial');
% (Repeat for other vertices as needed...)

% Save the plot as a PDF
% print('3dbands_114', '-dpdf', '-r1000');

%% Hexplot Function: Plot the hexagonal prism edges
function hexplot(hexagon, height1, height2)
    % Get vertices of the hexagon
    vertices = hexagon.Vertices;

    % Create bottom and top vertices for the prism
    verticesBottom = [vertices, height1 * ones(size(vertices, 1), 1)];
    verticesTop = [vertices, height2 * ones(size(vertices, 1), 1)];

    % Plot the edges of the prism
    hold on;
    for i = 1:size(vertices, 1)
        j = mod(i, size(vertices, 1)) + 1; % Wrap around to the first vertex
        plot3([verticesBottom(i, 1), verticesTop(i, 1)], ...
              [verticesBottom(i, 2), verticesTop(i, 2)], ...
              [verticesBottom(i, 3), verticesTop(i, 3)], 'k', 'LineWidth', 1, 'LineStyle', '--');
    end

    % Plot top and bottom faces of the prism
    fill3(verticesBottom(:, 1), verticesBottom(:, 2), verticesBottom(:, 3), 'w', 'FaceAlpha', 0, 'LineWidth', 1, 'LineStyle', '--');
    fill3(verticesTop(:, 1), verticesTop(:, 2), verticesTop(:, 3), 'w', 'FaceAlpha', 0, 'LineWidth', 1, 'LineStyle', '--');
end