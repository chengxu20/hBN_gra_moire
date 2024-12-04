% Calculate the density of states (DOS) at a serial wavelength
% Uncomment the following lines to use a parallel pool (optional)
% parpool('local', 96, 'IdleTimeout', 2400); % Initialize a parallel pool with 96 workers
% distcomp.feature('LocalUseMpiexec', false); % Avoid using mpiexec for local computations

% Parameters for the calculation
Vlist = [0.0452, 0.0134]; % Maximum potential values
alphalist = [0, 0.3]; % Layer-specific alpha values
knum = 400; % Number of k-points in the k-space mesh

% Nested loops to iterate over potential values and alpha values
for i = 1:2
    for j = 1:2

        % System setup parameters
        name = "NlayerGra/hBN"; % System name
        phase = 90; % Twist angle in degrees
        Vmax = Vlist(i); % Maximum potential for this iteration
        Vm = 2 * Vmax / 3 / sqrt(3); % Scaled potential value
        V = [0.0, 2 * Vmax / 3 / sqrt(3)]; % Potential vector
        Nlayer = 2; % Number of layers
        q_cut = 5; % Momentum cutoff
        lm = 200; % Length scale
        valley = 1; % Valley index
        delta = 0.0; % Energy offset
        % Tight-binding hopping parameters
        % t0: nearest neighbor, t1-t4: further neighbor hoppings
        t0 = -3.16; t1 = 0.381; t2 = 0; t3 = 0.38; t4 = 0.14;
        hop = [t0, t1, t2, t3, t4]; % Tight-binding parameters
        align = 1; % Alignment parameter
        bfield = [0, 0]; % Magnetic field vector

        % Initialize the system object
        NLG = system.NGra_twistedhBN(name, phase, V, q_cut, lm, Nlayer, valley, align, delta, hop, bfield);
        NLG.Alpha = [1, alphalist(j), 0, 0, 0]; % Set alpha parameters for the system

        % Compute the Q positions for the momentum cutoff
        [Q1, ~] = continuum.Others.Q_position(NLG.Lm, NLG.Q_cut, 1);
        t = -pi/6; % 30-degree rotation
        C30 = [cos(t), sin(t); -sin(t), cos(t)]; % Rotation matrix
        Q = (C30 * Q1')'; % Rotated Q points
        NLG.Q0 = Q; % Assign rotated Q points to the system

        % Compute reciprocal lattice vectors
        [Gm1, Gm2] = NLG.reciprocal_vectors();

        % Compute the Hamiltonian dimension
        dimH = length(Q) * Nlayer * 2 * abs(valley);

        % Generate a k-space mesh
        [Kx, Ky] = common.BZ.kmesh(knum, Gm1, Gm2);

        % Solve for eigenvalues (energy bands) over the k-space mesh
        [~, Enk] = common.Solve.cal_kmesh(NLG, Kx, Ky, dimH);

        % Create output directory and save results
        outputDir = "data/Vmax_" + num2str(Vm) + "/Vm_layer2_" + num2str(NLG.Alpha(2));
        mkdir(outputDir); % Create directory if it doesn't exist
        save(outputDir + "/enk", 'Enk', 'Kx', 'Ky'); % Save energy bands and k-space mesh

        % Display completion message and time cost
        disp("Finished calculation " + num2str(i) + " and the time cost is:");
        toc(tic);
    end
end
tic; % Start timing for the entire script