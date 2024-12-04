%% 1. Initialize Parameters
c = 0.001; % Broadening factor
Nband = 20; % Number of bands near CNP
Enum = 1000; % Number of energy points for DOS
dos = zeros(Nband, Enum); % Initialize DOS array
Nmax = zeros(1, Nband); % Initialize carrier density array
Vlist = 0.0:0.005:0.1; % List of potential values

% Create a figure for DOS plots
figure('Color', 'white');

%%2. Loop over Vmax values
for i = 1:21
    % Define the data path for the current Vmax
    path = "data/vhs185/Vmax_" + num2str(i) + "/Lm185";
    
    % Load energy band data
    load(path + "/enk.mat");
    
    % Extract dimensions and calculate Fermi energy
    dimH = size(Enk, 3);
    valband = ceil(dimH / 2);
    Ef = (min(Enk(:, :, valband), [], 'all') + max(Enk(:, :, valband + 1), [], 'all')) / 2;
    
    % Select 20 bands near CNP and shift by Ef
    Enk1 = Enk(:, :, valband - 1 - 9:valband - 1 + 10) - Ef;
    
    % Calculate DOS
    Emin = min(Enk1, [], 'all');
    Emax = max(Enk1, [], 'all');
    [Eaxis, Dos, dE] = common.Others.calc_dos(Enk1, c, Enum, Emin, Emax, Nband, 0);
    
    % Calculate carrier density
    Lm = 18.5; % Moiré superlattice size (nm)
    a = 2.46; % Lattice constant (Angstrom)
    [~, Efindex] = min(abs(Eaxis)); % Index of Fermi energy
    cellsize = Lm * 10 / a; % Cell size in lattice units
    n = cumtrapz(Eaxis, Dos) ./ (cellsize^2) ./ (sqrt(3) / 2 * a^2) * 10^16; % Carrier density
    nf = n(Efindex); % Carrier density at Fermi energy
    n = n - nf; % Shift carrier density to normalize at Ef
    
    % Store DOS data
    dos(i, :) = Dos;
    n = 2 .* n; % Scale density
    
    % Plot DOS for the current Vmax
    plot(n, Dos, 'Linestyle', '-', 'Color', 'k', 'LineWidth', 1.5);
    hold on;
end

% Customize plot appearance
xlabel("$n(cm^{-2})$", 'Interpreter', 'latex');
ylabel('DOS');
xlim([-4.5, 4.5] .* 10^12);
set(gca, 'Fontsize', 20, 'FontName', 'Times New Roman', 'LineWidth', 0.8);

%%3. Plot Density vs Potential
% Load precomputed Nmax values
A=load("data/density_vs_potential");
Vlist=A.Vlist;Nmax=A.Nmax;
% Calculate normalized carrier density
n0 = 1 / (cellsize^2) / (sqrt(3) / 2 * a^2) * 10^16; % Normalization factor
n1 = Nmax(2:end) .* 10^12 ./ 2 ./ n0; % Normalized carrier density
V = Vlist(2:end);

% Plot normalized density vs potential
figure('Color', 'white');
plot(V .* 1000, n1, '.', 'MarkerSize', 20); % Scatter plot
xlabel('$V_{max}$(meV)', 'Interpreter', 'latex');
ylabel("$n/n_0$", 'Interpreter', 'latex');
set(gca, 'Fontsize', 20, 'FontName', 'Arial', 'LineWidth', 0.8);

% Fit a linear model to the data
p = polyfit(n1, V, 1);
yfit = polyval(p, n1);
hold on;
plot(1000 .* yfit, n1, '-', 'LineWidth', 2); % Plot the fitted line

%% 4. Save the plot
% print('./fig/vhspos', '-dpdf');

%% 5. Calculate Effective Potential
nx = -3; % Example density value
yfit = polyval(p, nx); % Interpolate potential at nx
Veff = 2 * yfit / 3 / sqrt(3); % Calculate effective potential

% Output effective potential
disp("Effective Potential (Veff): " + Veff);