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

%net_f_dir = '/Users/nbarron/Desktop/Tau_Transport/neural_networks/nn_model_f_e5.pt';
%net_w_dir = '/Users/nbarron/Desktop/Tau_Transport/neural_networks/nn_model_w_e5.pt';

%net_f_dir = '/Users/nbarron/Desktop/NTM Edge Transmission/nn_model_f_e6.pt';
%net_w_dir = '/Users/nbarron/Desktop/NTM Edge Transmission/nn_model_w_e6.pt';

net_f_dir = '/Users/nbarron/Desktop/NTM Edge Transmission/NN-GS-Figure/output_data_e6_bias/model_e6__FValue_4_80_10_200_0.01.pt';
net_w_dir = '/Users/nbarron/Desktop/NTM Edge Transmission/NN-GS-Figure/output_data_e6_bias/model_e6__W1_12_40_3_200_0.01.pt';

input_size = [edge_count, 6];

net_f = get_nn_model_pytorch(net_f_dir, input_size);
net_w = get_nn_model_pytorch(net_w_dir, input_size);

for i = 41:42

    filename = ['/Users/nbarron/Desktop/Tau_Transport/frac07_sims/output_data/sim_outputs_' num2str(i) '.mat'];

    if exist(filename, 'file')

        load(filename)

        index = 1;
        init_path = save_data.output(index).Init_path;

        beta = 1e-6;
        gamma1 = save_data.params.gamma1;
        gamma2 = 0;
         
        frac = 0.7;
        lambda1 = save_data.params.lambda1;
        lambda2 = save_data.params.lambda2;
        delta = save_data.params.delta;
        epsilon = save_data.params.epsilon;

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
        init_rescale = 0.0020;
        conn_thresh = 0.8;
        use_sr_flux = 0;
        use_sr_w1 = 0;
            
        study = 'Hurtado';
        connectome_subset = 'Hippocampus+PC+RSP';

        sim_no = 1;

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
                                'use_nn',1,...
                                'net_f',net_f,...
                                'net_w',net_w);

        N_nn = model_output.Predicted.N;
        FVal_nn = model_output.Predicted.F;
        W1_nn = model_output.Predicted.W1;

        output_struct_nn = struct;
        output_struct_nn.N = N_nn;
        output_struct_nn.F = FVal_nn;
        output_struct_nn.W = W1_nn;

        dir_out = '/Users/nbarron/Desktop/frac07_sims/nn_sims_gs/nn_pred_';
        savepath = [dir_out num2str(i) '.mat'];

        save(savepath,'output_struct_nn');

    end

end
