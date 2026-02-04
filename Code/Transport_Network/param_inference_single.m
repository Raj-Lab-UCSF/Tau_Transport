function [] = param_inference_single(sim_num, rep_num)

%% Load Simulated Data Using Numeric Solver

%sim_num = 2;
%rep_num = 2;

sim_num_str = num2str(sim_num);
rep_num_str = num2str(rep_num);

sim_data_path = strcat('/Users/nbarron/Desktop/Tau_Transport/frac07_sims/output_data/sim_outputs_',sim_num_str,'.mat');

load(sim_data_path);

index = 1;

N_emp_sim = save_data.output(index).Model_Outputs.Predicted.N;
init_path_emp_sim = save_data.output(index).Init_path;

gamma_true = save_data.params.gamma1;
lambda_true = save_data.params.lambda1;
delta_true = save_data.params.delta;
epsilon_true = save_data.params.epsilon;

%% Load Connectome

matdir = '/Users/nbarron/Desktop/Tau_Transport/MatFiles';

connectome_subset = 'Hippocampus+PC+RSP';

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

%net_f_dir = '/Users/nbarron/Desktop/Tau_Transport/neural_networks/nn_model_f_e5.pt';
%net_w_dir = '/Users/nbarron/Desktop/Tau_Transport/neural_networks/nn_model_w_e5.pt';

net_f_dir = '/Users/nbarron/Desktop/NTM Edge Transmission/NN-GS-Figure/output_data_e6_bias/model_e6__FValue_4_80_10_200_0.01.pt';
net_w_dir = '/Users/nbarron/Desktop/NTM Edge Transmission/NN-GS-Figure/output_data_e6_bias/model_e6__W1_12_40_3_200_0.01.pt';

input_size = [edge_count, 6];

net_f = get_nn_model_pytorch(net_f_dir, input_size);
net_w = get_nn_model_pytorch(net_w_dir, input_size);

%% Define Loss Functions

loss_fun = @(x)NTM_Loss_single(matdir, net_f, net_w, init_path_emp_sim, N_emp_sim, x(1), x(2), x(3), x(4));

loss_fun_gamma = @(x)NTM_Loss(matdir, net_f, net_w, init_path_emp_sim, N_emp_sim, x(1), lambda_true, delta_true, epsilon_true);
loss_fun_lambda = @(x)NTM_Loss(matdir, net_f, net_w, init_path_emp_sim, N_emp_sim, gamma_true, x(1), delta_true, epsilon_true);
loss_fun_delta = @(x)NTM_Loss(matdir, net_f, net_w, init_path_emp_sim, N_emp_sim, gamma_true, lambda_true, x(1), epsilon_true);
loss_fun_epsilon = @(x)NTM_Loss(matdir, net_f, net_w, init_path_emp_sim, N_emp_sim, gamma_true, lambda_true, delta_true, x(1));

%% Test Loss Function

gamma_init = (rand() * (5e-2 - 5e-5)) + 5e-5;
lambda_init = (rand() * (0.1 - 0.01)) + 0.01;
delta_init = (rand() * 90) + 10;
epsilon_init = (rand() * 90) + 10;

%x0 = [gamma_init, lambda_init, delta_init, epsilon_init];

x0 = [gamma_init];

%x_min = [5e-5,0.01,10,10];
%x_max = [5e-2,0.1,100,100]

x_min = [5e-5];
x_max = [5e-2];


%initial_loss = loss_fun(x0);

%% Minimize Loss Function

xmin = fminsearchbnd(loss_fun_gamma,x0,x_min,x_max);

%% Save Files

%x0 = [2e-4, 0.05, 25, 25];
%xmin = x0;

savepath = strcat('/Users/nbarron/Desktop/inf_results/inf_output_gamma_',sim_num_str,'_',rep_num_str,'.mat');

output_struct = struct;
output_struct.x0 = x0;
output_struct.xmin = xmin;

gamma_true = save_data.params.gamma1;
lambda_true = save_data.params.lambda1;
delta_true = save_data.params.delta;
epsilon_true = save_data.params.epsilon;

output_struct.xtrue = [gamma_true, lambda_true, delta_true, epsilon_true];

save(savepath,'output_struct');