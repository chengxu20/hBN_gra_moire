%% 1. Initialize the system
% Define system parameters
name = "NlayerGra/hBN"; % System name
phase = 90; % Rotation angle in degrees
Vmax = 0.029; % Maximum potential strength
Vmoire = 2 * Vmax / 3 / sqrt(3); % Moiré potential strength (scaled)
V = [0.0, Vmoire]; % Potential array
Nlayer = 2; % Number of layers
q_cut = 5; % Momentum cutoff
lm = 114; % Moiré superlattice size (length scale)
valley = 2; % Valley index (±1 for K/K' valleys)
delta = 0.9 * 0.33 / 5 * 0; % Energy gap parameter (set to 0 in this case)

% Tight-binding hopping parameters
t0 = -3.16; t1 = 0.381; t2 = 0; t3 = 0.38; t4 = 0.14;
hop = [t0, t1, t2, t3, t4]; % Hopping integrals for nearest and next-nearest neighbors

% Other system parameters
align = 1; % Alignment parameter
bfield = [0, 0]; % Magnetic field (set to zero)

% Initialize the system object for twisted graphene on hBN
NLG = system.NGra_twistedhBN(name, phase, V, q_cut, lm, Nlayer, valley, align, delta, hop, bfield);
NLG.Alpha = [1, 1, 0, 0, 0]; % Alpha parameters for the system

%% 2. Define momentum space (Q points)
[Q1, ~] = continuum.Others.Q_position(NLG.Lm, NLG.Q_cut, 1); % Compute initial Q positions
t = -pi / 6; % Rotation angle (30 degrees)
C30 = [cos(t), sin(t); -sin(t), cos(t)]; % 2D rotation matrix
Q = (C30 * Q1')'; % Rotate Q points by 30 degrees
NLG.Q0 = Q; % Assign rotated Q points to the system

% Compute Hamiltonian dimension
dimH = length(Q) * Nlayer * 2 * abs(valley);

% Compute reciprocal lattice vectors
[Gm1, Gm2] = NLG.reciprocal_vectors();

%% 3. Solve and plot the band structure
band = Band(NLG); % Call the Band function to compute band structure
ylim([-200, 200]); % Set y-axis limits for the band structure plot

% Extract key energy levels from the band structure
E1 = max(band(:, dimH / 2 + 2), [], 'all'); % Maximum energy of conduction band
E2 = min(band(:, dimH / 2 + 3), [], 'all'); % Minimum energy of next conduction band
E3 = min(band(:, dimH / 2 - 1), [], 'all'); % Minimum energy of valence band
E4 = max(band(:, dimH / 2 - 2), [], 'all'); % Maximum energy of previous valence band

% Set plot labels
ylabel('E (meV)'); % y-axis label
set(gca, 'FontName', 'Arial', 'FontSize', 20); % Customize axis appearance

% Zoom into specific energy range and highlight energy levels
ylim([-140, -100]); % Adjust y-axis limits for zoomed view
yline(1000 * E3, '--', 'LineWidth', 2, 'Color', 'r'); % Highlight E3
yline(1000 * E4, '--', 'LineWidth', 2, 'Color', 'r'); % Highlight E4

%% 4. Auxiliary functions
% Band function to calculate and plot band structure
function band = Band(NLG)
    % Extract parameters
    Nlayer = NLG.Nlayer;
    [Q, ~] = continuum.Others.Q_position(NLG.Lm, NLG.Q_cut, 1);
    dimH = length(Q) * Nlayer * 2 * abs(NLG.Valley);
    [Gm1, Gm2] = NLG.reciprocal_vectors();

    % Define high-symmetry points in momentum space
    K1 = (Gm1 + 2 * Gm2) / 3; % K-point
    K2 = (2 * Gm1 + Gm2) / 3; % K'-point
    KM = (Gm1 + Gm2) / 2; % M-point
    Gamma = [0, 0]; % Gamma point

    % Define k-path for band structure calculation
    knum = 20; % Number of points per segment
    kpath = {Gamma; K2; KM; Gamma}; % Path through high-symmetry points
    step = norm(kpath{2} - kpath{1}) / knum; % Step size along k-path
    [Kpath, Kindex] = common.BZ.make_path(kpath, step); % Generate k-path

    % Solve for band structure along the k-path
    band = common.Solve.band_solve(NLG, Kpath, dimH);

    % Shift energy levels relative to the Fermi energy
    Ef = max(band(:, dimH / 2)); % Define Fermi energy
    band = band - Ef; % Shift bands to set Ef = 0

    % Plot the band structure
    point_labels = {'\Gamma', 'K', 'M', '\Gamma'}; % Labels for high-symmetry points
    common.Others.bands_plot(1000 * band, {"line", "k", 2}, point_labels, Kindex); % Plot bands
    % ylim([-0.03, 0.03]); % Optional zoom into specific energy range
end