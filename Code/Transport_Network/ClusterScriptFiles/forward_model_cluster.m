function []=forward_model_cluster(beta_in, gamma_in, lambda_in, delta_in, epsilon_in, iteration)

%% Define Directories and Filepaths

curpath = '/Users/nbarron/Desktop/Tau_Transport';
p = genpath(curpath);
addpath(p);

%% 2. Define Model Input Parameters

beta = str2double(beta_in);
gamma1 = str2double(gamma_in);
gamma2 = 0; %str2double(gamma_in);
lambda1 = str2double(lambda_in);
lambda2 = str2double(lambda_in);
delta = str2double(delta_in);
epsilon = str2double(epsilon_in);

frac = 0.7;

study = 'Hurtado';
connectome_subset = 'Hippocampus+PC+RSP';

%% Other Model Settings

L_int = 1000;
L1 = 200;
L2 = 200;
L_ais = 40;
L_syn = 40;
T = 1;
dt = [];
trange = [0:0.0025:0.1, 0.1050:0.005:0.3, 0.31:0.01:1];
resmesh = 'coarse';                             
plotting = 0;
reltol = 1e-4;
abstol = 1e-4;
fsolvetol = 1e-6;
%len_scale = 1e-3;
init_rescale = 0.0020;
conn_thresh = 0.8;
use_sr_flux = 0;
use_sr_w1 = 0;

sim_no = str2double(iteration);

%% Load Connectome

%matdir = [cd filesep 'MatFiles']; %% update this
matdir = '/Users/nbarron/Desktop/Tau_Transport/MatFiles';

load([matdir filesep 'Connectomes.mat'],'Connectomes'); % more updated version of connectome, should be minor
Conn = Connectomes.default;
Conn = Conn - diag(diag(Conn)); % remove the diagonal

thresh_C = 0.8 * mean(nonzeros(Conn(:)));

Conn(Conn < thresh_C) = 0;
Adj = logical(Conn);

switch connectome_subset
    case 'Hippocampus'
        Adj = Adj([27:37 (27+213):(37+213)], [27:37 (27+213):(37+213)]);
    case 'Hippocampus+PC+RSP'
        adjinds = [27:37,78:80,147];
        adjinds = [adjinds,adjinds+213];
        Adj = Adj(adjinds,adjinds);
    case 'RH'
        Adj = Adj(1:213,1:213);
    case 'LH'
        Adj = Adj(214:end,214:end);
    case 'Single'
        Adj = 1;
end

[edge_rows, edge_cols] = find(Adj);

edge_count = length(edge_rows);

%% Import Neural Net

net_f_dir = '/Users/nbarron/Desktop/NTM Edge Transmission/nn_model_f_wide.pt';
net_w_dir = '/Users/nbarron/Desktop/NTM Edge Transmission/nn_model_w_wide.pt';

input_size = [edge_count, 6];

net_f = get_nn_model_pytorch(net_f_dir, input_size);
net_w = get_nn_model_pytorch(net_w_dir, input_size);

%% Run NTM Model

output_struct = struct;

parpool(3)

parfor i = 1:3

    init_path = {};
    
    switch i
        case 1
            init_path = {'Entorhinal area, lateral part_L'};
        case 2
            init_path = {'Field CA1_L'};
        case 3
            init_path = {'Field CA1_L'; 'Entorhinal area, lateral part_L'};

    end

    model_output = NetworkTransportModel(matdir,'beta',beta,...
                                    'gamma1',gamma1,...
                                    'gamma2',gamma2,...
                                    'frac',frac,...
                                    'lambda1',lambda1,...
                                    'lambda2',lambda2,...
                                    'delta',delta,...
                                    'epsilon',epsilon,...
                                    'L_int',L_int,...
                                    'L1',L1,...
                                    'L2',L2,...
                                    'L_ais',L_ais,...
                                    'L_syn',L_syn,...,
                                    'T',T,...
                                    'dt',dt,...
                                    'trange',trange,...
                                    'resmesh',resmesh,...
                                    'plotting',plotting,...
                                    'reltol',reltol,...
                                    'abstol',abstol,...
                                    'fsolvetol',fsolvetol,...
                                    'init_rescale',init_rescale,...
                                    'init_path',init_path,...
                                    'study',study,...
                                    'connectome_subset',connectome_subset,...
                                    'sim_no',sim_no,...
                                    'conn_thresh',conn_thresh,...
                                    'use_sr_flux',use_sr_flux,...
                                    'use_sr_w1',use_sr_w1,...
                                    'use_nn',0,...
                                    'net_f',[],...
                                    'net_w',[]);

    output_struct(i).Model_Outputs = model_output;
    output_struct(i).Init_path = init_path
end

save_data = struct;

save_data.output = output_struct;

save_data.params.beta = beta;
save_data.params.gamma1 = gamma1;
save_data.params.gamma2 = gamma2;
save_data.params.lambda1 = lambda1;
save_data.params.lambda2 = lambda2;
save_data.params.delta = delta;
save_data.params.epsilon = epsilon;
save_data.params.frac = frac;
save_data.params.study = study;
save_data.params.conn_subset = connectome_subset;

%% Save Model Outputs 

savepath = '/Users/nbarron/Desktop';
savename = 'test_outputs';
savename = strcat(savename, '_', iteration);

save([savepath filesep savename '.mat'],'save_data');
