function model_outputs = SingleEdgeFluxCalculator_WrapperFun(matdir,varargin)
% Updated 24/11/14
if nargin < 1
    matdir = [cd filesep 'MatFiles'];
end

beta_ = 1e-04;
gamma1_ = 2e-03;
gamma2_ = 0;
delta_ = 1;
epsilon_ = 0.01;
lambda1_ = 0.01;
lambda2_ = 0.01;
tau_x0_ = 0.01;
tau_xL_ = 0.01;
% dt_ = 0.01;
% T_ = 0.1;
% trange_ = [];
frac_ = 0.92; % Average fraction of n diffusing (Konsack 2007)
L_int_ = 1000; % in micrometers
L1_ = 200;
L2_ = 200; 
L_ais_ = 40;
L_syn_ = 40;
resmesh_ = 'coarse';
% plotting_ = 1;
reltol_ = 1e-4;
abstol_ = 1e-4;
fsolvetol_ = 1e-6;
% connectome_subset_ = 'Hippocampus';
time_scale_ = 1;
len_scale_ = 1e-3;
sim_no_ = 1;
% conn_thresh_ = 'default';
use_sr_flux_ = 0;
sr_fun_flux_ = [];
sr_fun_em_ = [];
use_sr_w1_ = 0;
sr_fun_w1_ = [];

ip = inputParser;
% validChar = @(x) ischar(x);
validScalar = @(x) isnumeric(x) && isscalar(x) && (x>=0);
validLogical = @(x) validScalar(x) && (x == 0 || x == 1);
addParameter(ip, 'beta', beta_, validScalar);
addParameter(ip, 'gamma1', gamma1_, validScalar);
addParameter(ip, 'gamma2', gamma2_, validScalar);
addParameter(ip, 'delta', delta_, validScalar);
addParameter(ip, 'epsilon', epsilon_, validScalar);
addParameter(ip, 'frac', frac_, validScalar);
addParameter(ip, 'lambda1', lambda1_, validScalar);
addParameter(ip, 'lambda2', lambda2_, validScalar);
addParameter(ip, 'tau_x0', tau_x0_, validScalar);
addParameter(ip, 'tau_xL', tau_xL_, validScalar);
addParameter(ip, 'L_int', L_int_, validScalar);
addParameter(ip, 'L1', L1_, validScalar);
addParameter(ip, 'L2', L2_, validScalar);
addParameter(ip, 'resmesh', resmesh_);
addParameter(ip, 'L_ais', L_ais_);
addParameter(ip, 'L_syn', L_syn_);
addParameter(ip, 'reltol', reltol_, validScalar);
addParameter(ip, 'abstol', abstol_, validScalar);
addParameter(ip, 'fsolvetol', fsolvetol_, validScalar);
% addParameter(ip, 'connectome_subset', connectome_subset_);
addParameter(ip, 'len_scale', len_scale_, validScalar);
addParameter(ip, 'time_scale', time_scale_, validScalar);
% addParameter(ip, 'conn_thresh', conn_thresh_);
addParameter(ip, 'use_sr_flux', use_sr_flux_, validLogical);
addParameter(ip, 'sr_fun_flux', sr_fun_flux_);
addParameter(ip, 'sr_fun_em', sr_fun_em_);
addParameter(ip, 'use_sr_w1', use_sr_w1_, validLogical);
addParameter(ip, 'sr_fun_w1', sr_fun_w1_);

% addParameter(ip, 'study', study_);
% addParameter(ip, 'init_rescale', init_rescale_, validScalar);
% addParameter(ip, 'dt', dt_);
% addParameter(ip, 'T', T_);
% addParameter(ip, 'trange', trange_);
% addParameter(ip, 'init_path', init_path_);
% addParameter(ip, 'plotting', plotting_, validLogical);
addParameter(ip, 'sim_no', sim_no_, validScalar);
parse(ip, varargin{:});

% load([matdir filesep 'Mouse_Tauopathy_Data_HigherQ.mat'],'mousedata_struct'); 
% load([matdir filesep 'DefaultAtlas.mat'],'DefaultAtlas'); 
% load([matdir filesep 'CCF_labels.mat'],'CCF_labels');
% load([matdir filesep 'Connectomes.mat'],'Connectomes');
% Conn = Connectomes.default;
% if strcmp(ip.Results.conn_thresh,'default')
%     thresh_C = 0.8 * mean(nonzeros(Conn(:)));
% else
%     thresh_C = ip.Results.conn_thresh * mean(nonzeros(Conn(:)));
% end
% Conn(Conn < thresh_C) = 0;
% Adj = logical(Conn);
% Conn = readmatrix([matdir filesep 'mouse_connectome_19_01.csv']);
% Adj = readmatrix([matdir filesep 'mouse_adj_matrix_19_01.csv']);

% if ~isempty(ip.Results.init_path)
%     init_path = zeros(size(Conn,1),1);
%     for i = 1:length(ip.Results.init_path)
%         reghemstr = ip.Results.init_path{i};
%         reghemcell = split(reghemstr,'_');
%         reglog = ismember(CCF_labels(:,1),reghemcell{1});
%         if strcmp(reghemcell{2},'L')
%             hemlog = ismember(CCF_labels(:,4),'Left Hemisphere'); 
%         elseif strcmp(reghemcell{2},'R')
%             hemlog = ismember(CCF_labels(:,4),'Right Hemisphere'); 
%         else 
%             hemlog = ones(size(Conn,1),1);
%         end
%         init_path((reglog + hemlog) == 2) = 1;
%     end
% elseif isnan(mousedata_struct.(ip.Results.study).seed)
%     init_path = logical(mousedata_struct.(ip.Results.study).data(:,1));
%     init_path = DataToCCF_Transport(init_path,ip.Results.study,matdir);
% else
%     init_path = logical(mousedata_struct.(ip.Results.study).seed);
%     init_path = DataToCCF_Transport(init_path,ip.Results.study,matdir);
% end
% beta_new = ip.Results.beta * ip.Results.time_scale;
% gamma1_new = ip.Results.gamma1 * ip.Results.time_scale;
% gamma2_new = ip.Results.gamma2 * ip.Results.time_scale;
% taufun = @(x) ip.Results.init_rescale - (x +...
%     (gamma1_new * x.^2)./(beta_new - gamma2_new * x));
% options_taufun = optimset('TolFun',ip.Results.fsolvetol,'Display','off');
% init_rescale_n = fsolve(taufun,0,options_taufun);
% init_tau = init_rescale_n * init_path;

% switch ip.Results.connectome_subset
%     case 'Hippocampus'
%         inds = ismember(CCF_labels(:,3),'Hippocampus');
%     case 'Hippocampus+PC+RSP'
%         inds_hipp = ismember(CCF_labels(:,3),'Hippocampus');
%         inds_pc = ismember(CCF_labels(:,1),'Piriform area');
%         inds_rsp = ismember(CCF_labels(:,3),'Retrosplenial Area');
%         inds = logical(inds_hipp + inds_pc + inds_rsp);
%     case 'RH'
%         inds = ismember(CCF_labels(:,4),'Right Hemisphere');
%     case 'LH'
%         inds = ismember(CCF_labels(:,4),'Left Hemisphere');
%     case 'Single'
%         Adj = 1;
%     otherwise
%         inds = logical(ones(size(Conn,1),1)); %#ok<LOGL> 
% end
% Adj = 1; Conn = 1;
% Vol = DefaultAtlas.volumes(inds);
% init_tau = init_tau(inds);
% nroi = size(Adj,1);
% regnamecell = CCF_labels(inds,:);
% regnames = cell(size(regnamecell,1),1);
% for i = 1:length(regnames)
%     regname = regnamecell{i,1};
%     reghem = regnamecell{i,4};
%     if strcmp(reghem,'Right Hemisphere')
%         regnames{i} = [regname ' RH'];
%     else
%         regnames{i} = [regname ' LH'];
%     end
% end
% i_nonzero_init_tau = init_tau > 0;
% i_zero = ~i_nonzero_init_tau;

% if isempty(ip.Results.trange)
%     t = 0:ip.Results.dt:ip.Results.T;
% else
%     t = ip.Results.trange;
% end
% nt = length(t);
% N = zeros(nroi,nt); 
% N(:,1) = init_tau(:);
% netw_flux= zeros([nroi,size(N)]);
% W_1=zeros([nroi,size(N)]);
% W_2=zeros([nroi,size(N)]);
% S_ss=zeros([nroi,size(N)]);
% R_ss=zeros([nroi,size(N)]);
% Mass_edge=zeros([nroi,nroi,nt]);
% N_adj_0 = N(:,1) .* Adj;

[netw_flux,Mass_edge] = NetworkFluxCalculator(ip.Results.tau_x0,...
                                ip.Results.tau_xL,matdir,...
                                'beta',ip.Results.beta,...
                                'gamma1',ip.Results.gamma1,...
                                'gamma2',ip.Results.gamma2,...
                                'delta',ip.Results.delta,...
                                'epsilon',ip.Results.epsilon,...
                                'lambda1',ip.Results.lambda1,...
                                'lambda2',ip.Results.lambda2,...
                                'frac',ip.Results.frac,...
                                'L_int',ip.Results.L_int,...
                                'L1',ip.Results.L1,...
                                'L2',ip.Results.L2,...
                                'L_ais',ip.Results.L_ais,...
                                'L_syn',ip.Results.L_syn,...
                                'resmesh',ip.Results.resmesh,...
                                'reltol',ip.Results.reltol,...
                                'abstol',ip.Results.abstol,...
                                'fsolvetol',ip.Results.fsolvetol,... % compute the steady state network flux at time t0
                                'connectome_subset','Single',...
                                'time_scale',ip.Results.time_scale,...
                                'len_scale',ip.Results.len_scale,...
                                'use_sr_flux',ip.Results.use_sr_flux,...
                                'sr_fun_flux',ip.Results.sr_fun_flux,...
                                'sr_fun_em',ip.Results.sr_fun_em); 

model_outputs = struct;
% model_outputs.Predicted.N = N;
% model_outputs.Predicted.M = M;
model_outputs.Predicted.F = netw_flux;
% model_outputs.Predicted.W1 = W_1;
model_outputs.Predicted.EdgeMass = Mass_edge;
model_outputs.Parameters.beta = ip.Results.beta;
model_outputs.Parameters.gamma1 = ip.Results.gamma1;
model_outputs.Parameters.gamma2 = ip.Results.gamma2;
model_outputs.Parameters.delta = ip.Results.delta;
model_outputs.Parameters.epsilon = ip.Results.epsilon;
model_outputs.Parameters.lambda1 = ip.Results.lambda1;
model_outputs.Parameters.lambda2 = ip.Results.lambda2;
model_outputs.Parameters.tau_x0 = ip.Results.tau_x0;
model_outputs.Parameters.tau_xL = ip.Results.tau_xL;
model_outputs.Parameters.frac = ip.Results.frac;
model_outputs.Sim.L1 = ip.Results.L1;
model_outputs.Sim.L2 = ip.Results.L2;
model_outputs.Sim.L_int = ip.Results.L_int;
model_outputs.Sim.L_ais = ip.Results.L_ais;
model_outputs.Sim.L_syn = ip.Results.L_syn;
% if isempty(ip.Results.trange)
%     model_outputs.Sim.dt = ip.Results.dt;
%     model_outputs.Sim.T = ip.Results.T;
%     model_outputs.Sim.trange = 0:ip.Results.dt:ip.Results.T;
% else
%     model_outputs.Sim.dt = [];
%     model_outputs.Sim.T = ip.Results.trange(end);
%     model_outputs.Sim.trange = ip.Results.trange;
% end
model_outputs.Sim.len_scale = ip.Results.len_scale;
model_outputs.Sim.time_scale = ip.Results.time_scale;
% model_outputs.Sim.connectome_subset = ip.Results.connectome_subset;
% model_outputs.Sim.region_names = regnames;
% model_outputs.Sim.C = Conn;
% model_outputs.Sim.conn_thresh = ip.Results.conn_thresh;
% model_outputs.Sim.study = ip.Results.study;
% model_outputs.Sim.init_rescale = ip.Results.init_rescale;
% model_outputs.Sim.init_path = init_tau;
model_outputs.Sim.resmesh = ip.Results.resmesh;
model_outputs.Sim.rel_tol = ip.Results.reltol;
model_outputs.Sim.abs_tol = ip.Results.abstol;
model_outputs.Sim.fsolve_tol = ip.Results.fsolvetol;
% model_outputs.Sim.conn_thresh = ip.Results.conn_thresh;
model_outputs.Sim.use_sr_flux = ip.Results.use_sr_flux;
model_outputs.Sim.sr_fun_flux = ip.Results.sr_fun_flux;
model_outputs.Sim.sr_fun_em = ip.Results.sr_fun_em;
model_outputs.Sim.use_sr_w1 = ip.Results.use_sr_w1;
model_outputs.Sim.sr_fun_w1 = ip.Results.sr_fun_w1;

end