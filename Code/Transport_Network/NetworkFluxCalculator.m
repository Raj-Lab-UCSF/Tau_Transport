function [network_flux,mass_edge] = NetworkFluxCalculator(tau_x0,tau_xL,matdir,varargin)

if nargin < 3
    matdir = [cd filesep 'MatFiles'];
end

% % % 1. Preset values of flexible parameters
beta_ = 1e-06;
gamma1_ = 2e-05;
gamma2_ = 0;
delta_ = 1;
epsilon_ = 0.01;
lambda1_ = 0.01; %0.02  0.01 
lambda2_ = 0.01; %0.04  0.01;
frac_ = 0.7; % Average fraction of n diffusing (Konsack 2007) 0.92 - NOW USING 0.7!
L_int_ = 1000; % in micrometers
L1_ = 200;
L2_ = 200; 
L_ais_ = 40;
L_syn_ = 40;
resmesh_ = 'coarse';
reltol_ = 1e-4;
abstol_ = 1e-4;
fsolvetol_ = 1e-6;
connectome_subset_ = 'Hippocampus';
len_scale_ = 1e-3;
time_scale_ = 1;
conn_thresh_ = 'default';
use_sr_flux_ = 0;
sr_fun_flux_ = [];
sr_fun_em_ = [];

ip = inputParser;
validScalar = @(x) isnumeric(x) && isscalar(x) && (x>=0);
addParameter(ip, 'beta', beta_, validScalar);
addParameter(ip, 'gamma1', gamma1_, validScalar);
addParameter(ip, 'gamma2', gamma2_, validScalar);
addParameter(ip, 'delta', delta_, validScalar);
addParameter(ip, 'epsilon', epsilon_, validScalar);
addParameter(ip, 'frac', frac_, validScalar);
addParameter(ip, 'lambda1', lambda1_, validScalar);
addParameter(ip, 'lambda2', lambda2_, validScalar);
addParameter(ip, 'L_int', L_int_, validScalar);
addParameter(ip, 'L1', L1_, validScalar);
addParameter(ip, 'L2', L2_, validScalar);
addParameter(ip, 'resmesh', resmesh_);
addParameter(ip, 'L_ais', L_ais_);
addParameter(ip, 'L_syn', L_syn_);
addParameter(ip, 'reltol', reltol_, validScalar);
addParameter(ip, 'abstol', abstol_, validScalar);
addParameter(ip, 'fsolvetol', fsolvetol_, validScalar);
addParameter(ip, 'connectome_subset', connectome_subset_);
addParameter(ip, 'len_scale', len_scale_, validScalar);
addParameter(ip, 'time_scale', time_scale_, validScalar);
addParameter(ip, 'conn_thresh', conn_thresh_, validScalar);
addParameter(ip, 'use_sr_flux', use_sr_flux_);
addParameter(ip, 'sr_fun_flux', sr_fun_flux_);
addParameter(ip, 'sr_fun_em', sr_fun_em_);
addParameter(ip, 'net', [])
addParameter(ip, 'use_nn', 0)
%addParameter(ip, 'z_mu')
%addParameter(ip, 'z_sigma')

parse(ip, varargin{:});
beta_new = ip.Results.beta*ip.Results.time_scale;
gamma1_new = ip.Results.gamma1*ip.Results.time_scale;
gamma2_new = ip.Results.gamma2*ip.Results.time_scale;
L1_new = ip.Results.L1 * ip.Results.len_scale;
L2_new = ip.Results.L2 * ip.Results.len_scale;
L_int_new = ip.Results.L_int * ip.Results.len_scale;
L_ais_new = ip.Results.L_ais * ip.Results.len_scale;
L_syn_new = ip.Results.L_syn * ip.Results.len_scale;

use_nn = ip.Results.use_nn;
net = ip.Results.net;
%z_mu = ip.Results.z_mu;
%z_sigma = ip.Results.z_sigma;

load([matdir filesep 'Connectomes.mat'],'Connectomes'); % more updated version of connectome, should be minor
Conn = Connectomes.default;
Conn = Conn - diag(diag(Conn)); % remove the diagonal
if strcmp(ip.Results.conn_thresh,'default')
    thresh_C = 0.8 * mean(nonzeros(Conn(:)));
else
    thresh_C = ip.Results.conn_thresh * mean(nonzeros(Conn(:)));
end
Conn(Conn < thresh_C) = 0;
Adj = logical(Conn);
switch ip.Results.connectome_subset
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

nroi = size(Adj,1);

if use_nn
    
    [edge_rows, edge_cols] = find(Adj);

    edge_count = length(edge_rows);

    N1_list = zeros(edge_count,1);
    N2_list = zeros(edge_count,1);

    for j = 1:edge_count
        row_j = edge_rows(j);
        col_j = edge_cols(j);
        N1_list(j,1) = tau_x0(row_j,col_j);
        %N1_list(j,1) = tau_xL(row_j);
        N2_list(j,1) = tau_xL(col_j);
    end

    gamma_row = repmat(gamma1_new, edge_count, 1);
    lambda_row = repmat(ip.Results.lambda1, edge_count, 1);
    delta_row = repmat(ip.Results.delta, edge_count, 1);
    epsilon_row = repmat(ip.Results.epsilon, edge_count, 1);

    input_feats = [gamma_row, lambda_row, delta_row, epsilon_row, N1_list, N2_list];

    %input_feats

    input_feats_scaled = scale_inputs(input_feats);

    %input_feats_scaled

    input_feats_scaled_dl = dlarray(input_feats_scaled, 'UUU');

    nn_preds_scaled_dl = predict(net, input_feats_scaled_dl);

    nn_preds_scaled = extractdata(nn_preds_scaled_dl);

    nn_preds = unscale_outputs_f(nn_preds_scaled);

    %nn_preds

    network_flux = zeros(nroi);

    for j = 1:edge_count
        row_j = edge_rows(j);
        col_j = edge_cols(j);
        if tau_xL(col_j) > 0 || tau_x0(row_j,col_j) > 0
            network_flux(row_j, col_j) = nn_preds(j);
        end
    end

    mass_edge = zeros(nroi, nroi); % UNUSED BUT THROWS AN ERROR IF NOT SET

elseif ~isempty(ip.Results.sr_fun_flux) && ~isempty(ip.Results.sr_fun_em) && logical(ip.Results.use_sr_flux)
% % % 2a. Use symbolic expression from DSO
    
% theta = {gamma1, lambda, delta, epsilon, N1, N2}
fprintf('Using DSO Expression\n')

% Adj = readmatrix([matdir filesep 'mouse_adj_matrix_19_01.csv']);
% N1_mat = repmat(tau_xL,1,nroi);
% N2_mat = repmat(tau_xL.',nroi,1);
N1_mat = repmat(tau_xL,1,nroi);
N2_mat = repmat(tau_xL.',nroi,1);
network_flux = ip.Results.sr_fun_flux(gamma1_new,...
                                      ip.Results.lambda1,... % only when lambda1 == lambda2
                                      ip.Results.delta,...
                                      ip.Results.epsilon,...
                                      N1_mat,...
                                      N2_mat);

mass_edge = ip.Results.sr_fun_em(gamma1_new,...
                                  ip.Results.lambda1,... % only when lambda1 == lambda2
                                  ip.Results.delta,...
                                  ip.Results.epsilon,...
                                  N1_mat,...
                                  N2_mat);

% N1 = tau_xL;
% network_flux = zeros(nroi); mass_edge = network_flux;
% for i = 1:length(N1)
%     for j = 1:length(N1)
%         network_flux(i,j) = ip.Results.sr_fun_flux(gamma1_new,...
%                                                 ip.Results.lambda1,... % only when lambda1 == lambda2
%                                                 ip.Results.delta,...
%                                                 ip.Results.epsilon,...
%                                                 N1(i),...
%                                                 N1(j));
% 
%         mass_edge(i,j) = ip.Results.sr_fun_em(gamma1_new,...
%                                             ip.Results.lambda1,... % only when lambda1 == lambda2
%                                             ip.Results.delta,...
%                                             ip.Results.epsilon,...
%                                             N1(i),...
%                                             N1(j));
%     end
% end

zeroAdjinds = find(Adj(:) == 0);
network_flux(zeroAdjinds) = 0;
mass_edge(zeroAdjinds) = 0;
% mass_edge(mass_edge < 0) = 0; % bit of a kludge to avoid negatives from DSO

else
fprintf('Using NTM Integration\n')
% % % 2b. Definition of static constants
v_a = 0.7*ip.Results.len_scale * ip.Results.time_scale; % Average velocity (um/s) of anterograde active transpot (Konsack 2007)
v_r = 0.7*ip.Results.len_scale * ip.Results.time_scale; % Average velocity (um/s) of retrograde active transport (Konsack 2007)
diff_n = 12*ip.Results.len_scale^2 * ip.Results.time_scale; % Diffusivity (um^2/s) of n (Konsack 2007)

% % % 3b. Definition of the (inhomogeneous) xmesh
L_total = L1_new + L2_new + L_int_new; % size of the system
if strcmp(ip.Results.resmesh, 'fine')
    num_comp = 1000; % number of xmesh points
    num_ext = 100; % number of GM compartments per GM region
    num_int = num_comp - 2*(num_ext);
    xmesh1 = [linspace(0,L1_new-10*ip.Results.len_scale,num_ext-40),...
        (L1_new-9.75*ip.Results.len_scale):0.25*ip.Results.len_scale:L1_new];
    xmesh2 = [(L1_new+L_int_new+0.25*ip.Results.len_scale):0.25*ip.Results.len_scale:...
        (L1_new+L_int_new+10*ip.results.len_scale),...
        linspace(L1_new+L_int_new+10.25*ip.Results.len_scale,L_total,num_ext-40)];
    xmesh_int = [(L1_new+0.25*ip.Results.len_scale):0.25*ip.Results.len_scale:(L1_new+L_ais_new),...
        linspace(L1_new+L_ais_new+0.25*ip.Results.len_scale,L1_new+L_int_new-L_syn_new,...
        num_int-((L_ais_new+L_syn_new)/(0.25*ip.Results.len_scale))),...
        (L1_new+L_int_new-(L_syn_new-0.25*ip.Results.len_scale)):0.25*ip.Results.len_scale:(L1_new+L_int_new)];
elseif strcmp(ip.Results.resmesh, 'coarse')
    num_comp = 250; % number of xmesh points
    num_ext = 25; % number of compartments per SD region
    num_int = num_comp - 2*(num_ext);
    xmesh1 = [linspace(0,L1_new-10*ip.Results.len_scale,num_ext-5),...
        (L1_new-8*ip.Results.len_scale):2*ip.Results.len_scale:L1_new];
    xmesh2 = [(L1_new+L_int_new+2*ip.Results.len_scale):2*ip.Results.len_scale:...
        (L1_new+L_int_new+10*ip.Results.len_scale),...
        linspace(L1_new+L_int_new+12*ip.Results.len_scale,L_total,num_ext-5)];
    xmesh_int = [(L1_new+2*ip.Results.len_scale):2*ip.Results.len_scale:(L1_new+L_ais_new),...
        linspace(L1_new+L_ais_new+2*ip.Results.len_scale,L1_new+L_int_new-L_syn_new,...
        num_int-((L_ais_new+L_syn_new)/(2*ip.Results.len_scale))),...
        (L1_new+L_int_new-(L_syn_new-2*ip.Results.len_scale)):2*ip.Results.len_scale:(L1_new+L_int_new)];
end
xmesh = [xmesh1, xmesh_int, xmesh2];

% % % 4b. Steady State Calculation
% Follows the derivation of Michiel Bertsh

% % % 4b1. Presynaptic somatodendritic compartment
presyn_mask = spatial_mask('presyn');
xmesh_presyn = xmesh(presyn_mask);
n0 = @(B) B;
options = odeset('RelTol',ip.Results.reltol,'AbsTol',ip.Results.abstol,'NonNegative',1:length(n0));
n_ss_presyn = @(A,B) ode45(@(x,n)ode_ss_n(x,A,n,diff_n),[0,L1_new],n0(B),options);
n_ss_presyn = @(A,B,x) deval(n_ss_presyn(A,B),x);
x1 = xmesh_presyn(end);

% % % 4b2. Axon initial segment
ais_mask = spatial_mask('ais');
xmesh_ais = xmesh(ais_mask);
n_ss_ais = @(A,B) ode45(@(x,n)ode_ss_n(x,A,n,diff_n*ip.Results.lambda1),[L1_new,L1_new+L_ais_new],n_ss_presyn(A,B,x1),options);
n_ss_ais = @(A,B,x) deval(n_ss_ais(A,B),x);
x2 = xmesh_ais(end);

% % % 4b3. Axon
axon_mask = spatial_mask('axon');
xmesh_axon = xmesh(axon_mask);
x3 = xmesh_axon(end);
n_ss_axon = @(A,B) ode45(@(x,n)ode_ss_axon(x,A,n),[L1_new+L_ais_new,...
   L1_new+L_int_new-L_syn_new],n_ss_ais(A,B,x2),options);
n_ss_axon = @(A,B,x) deval(n_ss_axon(A,B),x);
% nx_init = @(A_,B) max((-A_.*x1/diff_n + B - A_.*(x2 - x1)/(ip.Results.lambda1*diff_n)),0);
% options = odeset('RelTol',ip.Results.reltol,'AbsTol',ip.Results.abstol,'NonNegative',1:length(nx_init));
% n_ss_axon = @(A,B) ode45(@(x,n)ode_ss_axon(x,A,n),[L1_new+L_ais_new,...
%    L1_new+L_int_new-L_syn_new],nx_init(A,B),options);
% n_ss_axon = @(A,B,x) deval(n_ss_axon(A,B),x);

% % % 4b4. Synaptic cleft
syncleft_mask = spatial_mask('syncleft');
xmesh_syncleft = xmesh(syncleft_mask);
x4 = xmesh_syncleft(end);
% n_ss_syncleft = @(A,B,x) max((n_ss_axon(A,B,x3) - A.*(x-x3)/(diff_n*ip.Results.lambda2)),0);
n_ss_syncleft = @(A,B) ode45(@(x,n)ode_ss_n(x,A,n,(diff_n*ip.Results.lambda2)),...
    [L1_new+L_int_new-L_syn_new,L1_new+L_int_new],n_ss_axon(A,B,x3),options);
n_ss_syncleft = @(A,B,x) deval(n_ss_syncleft(A,B),x);

% % % 4b5. Postsynaptic somatodendritic compartment
postsyn_mask = spatial_mask('postsyn');
xmesh_postsyn = xmesh(postsyn_mask);
x5 = xmesh_postsyn(end);
n_ss_postsyn = @(A,B,x) (n_ss_syncleft(A,B,x4) - A.*(x-x4)/diff_n); 
f_ss=@(A,B,C)(n_ss_syncleft(A,B,x4) - A.*(x5-x4)/diff_n-C);

% % % 5b. Flux calculation on network 

% Adj = readmatrix([matdir filesep 'mouse_adj_matrix_19_01.csv']);
load([matdir filesep 'Connectomes.mat'],'Connectomes'); % more updated version of connectome, should be minor
Conn = Connectomes.default;
Conn = Conn - diag(diag(Conn)); % remove the diagonal
if strcmp(ip.Results.conn_thresh,'default')
    thresh_C = 0.8 * mean(nonzeros(Conn(:)));
else
    thresh_C = ip.Results.conn_thresh * mean(nonzeros(Conn(:)));
end
Conn(Conn < thresh_C) = 0;
Adj = logical(Conn);
switch ip.Results.connectome_subset
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

nroi = size(Adj,1);
network_flux = zeros(nroi);
for i = 1:nroi
%     fprintf('ROI %d/%d\n',i,nroi)
    Adj_in = logical(Adj(:,i));
    tau_xL_i = tau_xL(i);
    tau_xL_i = repmat(tau_xL_i,length(Adj_in),1);
    i_app = logical(((tau_x0(:,i) > 0) + (tau_xL_i > 0)) .* (Adj_in)); 
    tau_x0_i = tau_x0(i_app,i);
    tau_xL_i = tau_xL_i(i_app);
    if ~isempty(tau_x0_i)
        x0 = zeros(length(tau_x0_i),1);
        options = optimset('TolFun',ip.Results.fsolvetol,'Display','off');
        fun_ss = @(A) f_ss(A,tau_x0_i,tau_xL_i);
        network_flux(i_app,i) = fsolve(fun_ss,x0,options);
    end
end
%%% 5b1.mass edge calculation
network_flux_mass_edge=network_flux(:);
tau_x0=tau_x0(:);
n_ss_presyn_eval= n_ss_presyn(network_flux_mass_edge, tau_x0,xmesh_presyn);
n_ss_ais_eval=n_ss_ais(network_flux_mass_edge,tau_x0,xmesh_ais);
n_ss_axon_eval=n_ss_axon(network_flux_mass_edge,tau_x0, xmesh_axon);
n_ss_syncleft_eval=n_ss_syncleft(network_flux_mass_edge,tau_x0,xmesh_syncleft);
n_ss_postsyn_eval=n_ss_postsyn(network_flux_mass_edge,tau_x0, xmesh_postsyn);
% n_ss_eval=[n_ss_presyn_eval n_ss_ais_eval n_ss_axon_eval n_ss_syncleft_eval n_ss_postsyn_eval];

n_ss_eval_1=[n_ss_presyn_eval n_ss_ais_eval n_ss_axon_eval];
n_m_ss=n_ss_eval_1+(ip.Results.gamma1 * n_ss_eval_1.^2)./(ip.Results.beta-ip.Results.gamma2 *n_ss_eval_1);
n_m_ss_postsyn=n_ss_postsyn_eval+(ip.Results.gamma1 * n_ss_postsyn_eval.^2)./(ip.Results.beta-ip.Results.gamma2 *n_ss_postsyn_eval);
mass_edge=trapz([xmesh_presyn xmesh_ais xmesh_axon], n_m_ss,2)+trapz(xmesh_syncleft, n_ss_syncleft_eval,2)+trapz(xmesh_postsyn, n_m_ss_postsyn,2);
mass_edge=reshape(mass_edge,nroi,nroi);
end

% % % 6b. Functions
    function nprime=ode_ss_n(x,A,n,D) %#ok<INUSD> 
        nprime=1/D*-A;
    end

    function nprime=ode_ss_axon(x,A,n) %#ok<INUSD> 
%         nprime = -A/(ip.Results.frac*diff_n)+(1/diff_n)*((1-ip.Results.frac)./...
%         ip.Results.frac).*n.*((v_a*(1+ip.Results.delta.*n).*(1-((gamma_new...
%         *ip.Results.epsilon.*n.^2)./beta_new))-v_r));
        nprime=1/diff_n*-A/ip.Results.frac+1/diff_n*((1-ip.Results.frac)./...
            ip.Results.frac).*n.*((v_a*(1+ip.Results.delta.*n).*(1-((gamma1_new...
            *ip.Results.epsilon.*n.^2)./(beta_new-gamma2_new.*n)))-v_r)); 
    end

    function [maskvals] = spatial_mask(compartment)
        switch compartment
            case 'presyn'
                maskvals = (xmesh <= L1_new);
            case 'ais'
                maskvals = logical(-1 + (xmesh > L1_new) + ...
                    (xmesh < (L1_new + L_ais_new)));
            case 'axon'
                maskvals = logical(-1 + (xmesh >= L1_new + L_ais_new) + ...
                    (xmesh <= (L1_new + L_int_new - L_syn_new)));
            case 'syncleft'
                maskvals = logical(-1 + (xmesh > L1_new + L_int_new - L_syn_new) + ...
                    (xmesh < (L1_new + L_int_new))); 
            case 'postsyn'
                maskvals = logical(-1 + (xmesh >= L1_new + L_int_new) + ...
                    (xmesh <= (L1_new + L_int_new + L2_new)));
            otherwise
                error('Incorrect compartment specification')
        end
    end

end
