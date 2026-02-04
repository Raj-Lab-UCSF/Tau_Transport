%% Load Data

matdir = '/Users/nbarron/Desktop/NTM_Test';
load([matdir filesep 'hippocampome_final_round2_v2.mat'])
%load('/Users/nbarron/Desktop/frac07_sims/output_data/sim_outputs_1.mat')

%% Grab Parameter Values
index = 1;
params = output_struct.Parameter_Grid(index,:);
%params = save_data.params;

%% Get Params

% beta = params(1);
% gamma1 = params(2);
% gamma2 = params(3);
% %frac = params(4);
% frac = 0.7;
% lambda1 = params(5);
% lambda2 = params(6);
% delta = params(7);
% epsilon = params(8);

% beta = 1e-6;
% gamma1 = 0.004;
% gamma2 = params(3);
% %frac = params(4);
% frac = 0.7;
% lambda1 = 0.025;
% lambda2 = 0.025;
% delta = 10;
% epsilon = 100;

beta = 1e-6;
gamma1 = 5e-5;
gamma2 = params(3);
%frac = params(4);
frac = 0.7;
lambda1 = 0.01;
lambda2 = 0.01;
delta = 10;
epsilon = 50;

%% Grab other parameter values

other_data = output_struct.Simulations(1,index).Model_Outputs;

sim_settings = other_data.Sim;


%% Next section

L_int = sim_settings.L_int;
L1 = sim_settings.L1;
L2 = sim_settings.L2;
L_ais = sim_settings.L_ais;
L_syn = sim_settings.L_syn;
T = sim_settings.T;
dt = sim_settings.dt;
trange = sim_settings.trange;
resmesh = sim_settings.resmesh;                             
plotting = 0;
reltol = sim_settings.rel_tol;
abstol = sim_settings.abs_tol;
fsolvetol = sim_settings.fsolve_tol;
init_rescale = 0.02; %sim_settings.init_rescale;
%init_path = sim_settings.init_path;
init_path = {'Entorhinal area, lateral part_L'};
study = sim_settings.study;
connectome_subset = sim_settings.connectome_subset;
sim_no = 1;
conn_thresh = 0.8;
use_sr_flux = 0;
use_sr_w1 = 0;

%% Load Connectome

matdir = [cd filesep 'MatFiles'];

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

%% Load Neural Net

%net_f_dir = '/Users/nbarron/Desktop/Tau_Transport/Code/Transport_Network/NeuralNet/nn_traced_f.pt';
%net_w_dir = '/Users/nbarron/Desktop/Tau_Transport/Code/Transport_Network/NeuralNet/nn_traced_w.pt';

net_f_dir = '/Users/nbarron/Desktop/NTM Edge Transmission/nn_model_f_e5.pt';
net_w_dir = '/Users/nbarron/Desktop/NTM Edge Transmission/nn_model_w_e5.pt';

input_size = [edge_count, 6];

net_f = get_nn_model_pytorch(net_f_dir, input_size);
net_w = get_nn_model_pytorch(net_w_dir, input_size);
%net_f = get_nn_model_onnx(net_f_dir, input_size);
%net_w = get_nn_model_onnx(net_w_dir, input_size);

%% Run NTM Model

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
                                'net_f',net_f,...
                                'net_w',net_w);




