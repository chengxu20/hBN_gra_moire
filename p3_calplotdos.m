% %% 1. Load and get the bands near Fermi surface
% % Define maximum potential and data path
% Vmax = 11;
% path = "data/dos/Vmax_" + num2str(Vmax) + "/";
% 
% % Load energy band data
% load(path + "enk.mat");
% 
% % Determine dimensions and locate the valence band
% dimH = size(Enk, 3);
% valband = ceil(dimH / 2);
% 
% % Calculate the Fermi energy (Ef) as the average of valence and conduction band edges
% Ef = (min(Enk(:, :, valband), [], 'all') + max(Enk(:, :, valband + 1), [], 'all')) / 2;
% 
% % Extract 20 bands near the charge neutrality point (CNP) and shift by Ef
% Enk1 = Enk(:, :, valband - 1 - 9:valband - 1 + 10) - Ef;
% 
% %% 2. Calculate the density of states (DOS)
% dimH = size(Enk, 3); % Redundant but kept for consistency
% c = 0.0005; % Broadening factor
% Nband = 20; % Number of bands considered for DOS calculation
% Enum = 1000; % Number of energy points in the DOS axis
% Emin = min(Enk1, [], 'all'); % Minimum energy
% Emax = max(Enk1, [], 'all'); % Maximum energy
% 
% % Compute DOS using a Gaussian broadening method
% [Eaxis, Dos, dE] = common.Others.calc_dos(Enk1, c, Enum, Emin, Emax, Nband, 1);
% 
% % Set the x-axis limits for visualization
% xlim([Emin, Emax]);
% 
% %% 3. Calculate the carrier density
% Lm = 20; % Moiré superlattice size
% a = 2.46; % Lattice constant in angstroms
% 
% % Find the index of the Fermi energy in the DOS axis
% [~, Efindex] = min(abs(Eaxis));
% 
% % Calculate cell size in terms of lattice constant
% cellsize = Lm * 10 / a;
% 
% % Compute cumulative carrier density (integrated DOS)
% n = cumtrapz(Eaxis, Dos) ./ (cellsize^2) ./ (sqrt(3) / 2 * a^2) * 10^16;
% 
% % Shift carrier density to normalize at Ef
% nf = n(Efindex);
% n = n - nf;

% Optional visualization of carrier density and DOS
% Uncomment to plot carrier density vs DOS
% figure('Color', 'white');
% plot(2 .* n, Dos, 'LineStyle', '-', 'Color', 'k', 'LineWidth', 1.5);
% xlabel("$n(cm^{-2})$", 'Interpreter', 'latex');
% ylabel('DOS');
% xlim([-4, 4] .* 10^12);
% set(gca, 'FontSize', 20, 'FontName', 'Times New Roman', 'LineWidth', 0.8);

% Save data if needed
% save("data/dos/dos" + num2str(Vmax), 'Eaxis', 'n', 'Dos');

%% 4. Plot DOS for different Vmax together
% Load and plot DOS for Vmax = 29 meV
load('data/dos/dos29.mat');
figure('Color', 'white');
plot(2 .* n, Dos - 100, 'LineStyle', '-', 'Color', 'b', 'LineWidth', 1.5);
hold on;

% Load and plot DOS for Vmax = 11 meV
load('data/dos/dos11.mat');
plot(2 .* n, Dos, 'LineStyle', '-', 'Color', 'r', 'LineWidth', 1.5);

% Customize plot appearance
xlabel("$n(cm^{-2})$", 'Interpreter', 'latex'); % Carrier density axis
ylabel('DOS'); % DOS axis
legend('V_{max}=29 meV', 'V_{max}=11 meV'); % Add legend
legend('boxoff'); % Remove legend box
xlim([-4, 4] .* 10^12); % Set x-axis limits
set(gca, 'FontSize', 20, 'FontName', 'Arial', 'LineWidth', 0.8);

% Save the figure as a PDF
% Uncomment to save the plot
% exportgraphics(gcf, 'dostogether34.eps', 'ContentType', 'vector');
% print('./fig/dostogether', '-dpdf');