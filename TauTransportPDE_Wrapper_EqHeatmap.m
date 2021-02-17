% Wrapper script for generating output files for TauTransportPDE for
% generating equilibrium heatmaps with respect to delta and epsilon

% 1. Define directories for saving outputs
curpath = cd;
simpath = '/Users/justintorok/Documents/MATLAB/TAUPDE/SimResults/';
simstr = 'heatmap_plot_default'; % for saving the outputs

% 2a. Define actively tuned parameters 
beta = 1e-6; % default = 1e-6
gamma = 2e-5; % default = 2e-5
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

% 2c. Define ICs
n0 = [zeros(1,L_ais), 20*ones(1,(L_int-L_ais-L_syn)), zeros(1,L_syn)]; % default = [zeros(1,L_ais), 20*ones(1,(L_int-L_ais-L_syn)), zeros(1,L_syn)]
m0 = zeros(1,L_int); % default = zeros(1,L_int)
N1_0 = 0; % default = 0
N2_0 = 0; % default = 0
M1_0 = 0; % default = 0
M2_0 = 0; % default = 0
n_init0 = [N1_0*ones(1,L1),n0,N2_0*ones(1,L2)];
m_init0 = [M1_0*ones(1,L1),m0,M2_0*ones(1,L2)];

total_mass = 184; % default = 184, ensures all simulations have the same mass
nonnorm_total = trapz(n_init0) + trapz(m_init0);
rescale_factor = total_mass/nonnorm_total; 

% 2d. Define xmesh resolution
resmesh = 'coarse'; % 'fine' or 'coarse' - use 'coarse' for faster, less precise simulations

% 2e. Define ranges for delta and epsilon, instantiate matrices for storing n
% and m
deltalist = 0:0.05:1;
epsilonlist = 0:0.05:1;
if strcmp(resmesh,'coarse')
    n_mat = zeros(length(deltalist),length(epsilonlist),tsteps,250);
elseif strcmp(resmesh,'fine')
    n_mat = zeros(length(deltalist),length(epsilonlist),tsteps,1000);
end
m_mat = n_mat;
biasmat = zeros(length(deltalist),length(epsilonlist),tsteps);

% 3. Run TauTransportPDE
for i = 1:length(deltalist)
    for j = 1:length(epsilonlist)
        fprintf('Simulation %d/%d \n',j+(i-1)*length(epsilonlist),length(deltalist)*length(epsilonlist))
        delta = deltalist(i);
        epsilon = epsilonlist(j);
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
            biasmat(i,j,k) = (total2 - total1)/(total2 + total1);
        end
        n_mat(i,j,:,:) = n;
        m_mat(i,j,:,:) = m;
    end
end

% 4. Save output file
cd(simpath)
save([simstr '.mat'])
cd(curpath)
clear