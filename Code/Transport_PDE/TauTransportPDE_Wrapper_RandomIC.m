% Wrapper script for generating output files for TauTransportPDE with
% random initialization

% 1. Define directories for saving outputs
curpath = cd;
simpath = [curpath filesep 'SampleFiles'];
if ~isfolder(simpath)
    mkdir(simpath)
end
simstr = 'random_init_ant_test'; % for saving the outputs

% 2a. Define actively tuned parameters 
beta = 1e-6; % default = 1e-6
gamma = 2e-5; % default = 2e-5
delta = 1; % default for anterograde/retrograde/no-bias = 1/0.01/1
epsilon = 0.01; % default for anterograde/retrograde/no-bias = 0.01/1/0.35
frac = 0.92; % default = 0.92
alpha = 0; % Not explored in the model (default = 0)

% 2b. Define other parameters
lambda = 0.01; % default = 0.01
L_int = 1000; % default = 1000; in microns
L1 = 200; % default = 200
L2 = 200; % default = 200
L_ais = 40; % default = 40
L_syn = 40; % default = 40
T = 5e7; % default = 5e7
tsteps = 250; % default = 250

% 2c. Define xmesh resolution
resmesh = 'coarse'; % 'fine' or 'coarse' - use 'coarse' for faster, less precise simulations

% 2d. Define number of random initial conditions to try and set rng seed
numiters = 5; % default = 100 
rng(0);

% 3. Run TauTransportPDE with random intialization and calculate outputs
if strcmp(resmesh,'coarse')
    n_mat = zeros(numiters,tsteps,250);
elseif strcmp(resmesh,'fine')
    n_mat = zeros(numiters,tsteps,1000);
end
m_mat = n_mat;
biasmat = zeros(numiters,tsteps);
for i = 1:numiters
    fprintf('Simulation %d/%d \n',i,numiters)
    
    % Create random IC with same total pathology (184)
    n0 = rand(1,10); m0 = rand(1,10);
    N1_0 = rand; N2_0 = rand;  M1_0 = rand; M2_0 = rand;
    n0 = [zeros(1,L_ais),interp1(linspace(1,920,10),n0,1:920,'makima'),zeros(1,L_syn)]; 
    n0(n0<0) = 0;
    m0 = [zeros(1,L_ais),interp1(linspace(1,920,10),m0,1:920,'makima'),zeros(1,L_syn)]; 
    m0(m0<0) = 0;
    n_init0 = [N1_0*ones(1,L1),n0,N2_0*ones(1,L2)];
    m_init0 = [M1_0*ones(1,L1),m0,M2_0*ones(1,L2)];
    nonnorm_total = trapz(n_init0) + trapz(m_init0);
    rescale_factor = 184/nonnorm_total;
    
    % Run TauTransportPDE
    [n, m, xmesh, trange, jn, jm] = TauTransportPDE('alpha', alpha,...
                                    'beta', beta,...
                                    'gamma', gamma,...
                                    'delta', delta,...
                                    'epsilon', epsilon,...
                                    'lambda', lambda,...
                                    'frac',frac,...
                                    'L_int', L_int,...
                                    'L1', L1,...
                                    'L2', L2,...
                                    'n0', n0*rescale_factor,...
                                    'm0', m0*rescale_factor,...
                                    'N1_0', N1_0*rescale_factor,...
                                    'N2_0', N2_0*rescale_factor,...
                                    'M1_0', M1_0*rescale_factor,...
                                    'M2_0', M2_0*rescale_factor,...
                                    'T', T,...
                                    'tsteps', tsteps,...
                                    'resmesh', resmesh);

    % Calculate N1(t), N2(t), and n(x,t) (and M analogs)
    GM1_bound = find(xmesh == L1);
    GM2_bound = find(xmesh == (L1+L_int));
    ais_bound = find(xmesh == (L1+40));
    syn_bound = find(xmesh == (L1+L_int-40));
    N1 = zeros(1,tsteps);
    M1 = zeros(1,tsteps);
    N2 = zeros(1,tsteps);
    M2 = zeros(1,tsteps);
    n_mean = zeros(1,tsteps);
    m_mean = zeros(1,tsteps);
    GM1_inds = 1:GM1_bound;
    GM2_inds = GM2_bound+1:length(xmesh);
    WM_inds = GM1_bound+1:GM2_bound;
    for k = 1:tsteps
        N1(k) = trapz(xmesh(GM1_inds),n(k,GM1_inds))/L1;
        M1(k) = trapz(xmesh(GM1_inds),m(k,GM1_inds))/L1;
        N2(k) = trapz(xmesh(GM2_inds),n(k,GM2_inds))/L2;
        M2(k) = trapz(xmesh(GM2_inds),m(k,GM2_inds))/L2;
        n_mean(k) = trapz(xmesh(WM_inds),n(k,WM_inds))/L_int; 
        m_mean(k) = trapz(xmesh(WM_inds),m(k,WM_inds))/L_int;
        total1 = N1(k) + M1(k); total2 = N2(k) + M2(k);
        biasmat(i,k) = (total2 - total1)/(total2 + total1);
    end
    n_mat(i,:,:) = n;
    m_mat(i,:,:) = m;
end

% 4. Save output file
save([simpath filesep simstr '.mat']) % saves all outputs, including parameters
clear