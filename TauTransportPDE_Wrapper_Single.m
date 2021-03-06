% Wrapper script for generating output files for TauTransportPDE for single
% parameter sets

% 1. Define directories for saving outputs
curpath = '~/Documents/MATLAB';
simpath = [curpath filesep 'SampleFiles'];
if ~isfolder(simpath)
    mkdir(simpath)
end
simstr = 'test'; % for saving the outputs

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

% 3. Run TauTransportPDE
[n, m, xmesh, trange, jn, jm] = TauTransportPDE('alpha', alpha,...
                                'beta', beta,...
                                'gamma', gamma,...
                                'delta', delta,...
                                'epsilon', epsilon,...
                                'lambda', lambda,...
                                'L_int', L_int,...
                                'L1', L1,...
                                'L2', L2,...
                                'n0', n0*rescale_factor,...
                                'm0', m0*rescale_factor,...
                                'N1_0', N1_0*rescale_factor,...
                                'N2_0', N2_0*rescale_factor,...
                                'M1_0', M1_0*rescale_factor,...
                                'M2_0', M2_0*rescale_factor,...
                                'frac', frac,...
                                'T', T,...
                                'tsteps', tsteps,...
                                'resmesh', resmesh,...
                                'L_ais', L_ais,...
                                'L_syn', L_syn);

% 4. Calculate N1(t), N2(t), and n(x,t) (and M analogs)
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
for i = 1:tsteps
    N1(i) = trapz(xmesh(GM1_inds),n(i,GM1_inds))/L1;
    M1(i) = trapz(xmesh(GM1_inds),m(i,GM1_inds))/L1;
    N2(i) = trapz(xmesh(GM2_inds),n(i,GM2_inds))/L2;
    M2(i) = trapz(xmesh(GM2_inds),m(i,GM2_inds))/L2;
    n_mean(i) = trapz(xmesh(WM_inds),n(i,WM_inds))/L_int; 
    m_mean(i) = trapz(xmesh(WM_inds),m(i,WM_inds))/L_int;
end        

% 5. Save output file
save([simpath filesep simstr '.mat']) % saves all outputs, including parameters
clear