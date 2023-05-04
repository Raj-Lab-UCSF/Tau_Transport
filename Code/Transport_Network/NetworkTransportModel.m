function model_outputs = NetworkTransportModel(matdir,varargin)
% Translated from ode_modelBFRTT03_09.m
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
study_ = 'Hurtado';
init_path_ = [];
init_rescale_ = 1e-2;
dt_ = 0.01;
T_ = 0.1;
frac_ = 0.92; % Average fraction of n diffusing (Konsack 2007)
L_int_ = 1000; % in micrometers
L1_ = 200;
L2_ = 200; 
L_ais_ = 40;
L_syn_ = 40;
resmesh_ = 'coarse';
plotting_ = 1;
reltol_ = 1e-4;
abstol_ = 1e-4;
fsolvetol_ = 1e-6;
connectome_subset_ = 'Hippocampus';
time_scale_ = 1;
len_scale_ = 1e-3;

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

addParameter(ip, 'study', study_);
addParameter(ip, 'init_rescale', init_rescale_, validScalar);
addParameter(ip, 'dt', dt_, validScalar);
addParameter(ip, 'T', T_, validScalar);
addParameter(ip, 'init_path', init_path_);
addParameter(ip, 'plotting', plotting_, validLogical);
parse(ip, varargin{:});

load([matdir filesep 'Mouse_Tauopathy_Data_HigherQ.mat'],'mousedata_struct'); 
load([matdir filesep 'DefaultAtlas.mat'],'DefaultAtlas'); 
load([matdir filesep 'CCF_labels.mat'],'CCF_labels');
Conn = readmatrix([matdir filesep 'mouse_connectome_19_01.csv']);
Adj = readmatrix([matdir filesep 'mouse_adj_matrix_19_01.csv']);

if ~isempty(ip.Results.init_path)
    init_path = zeros(size(Conn,1),1);
    for i = 1:length(ip.Results.init_path)
        reghemstr = ip.Results.init_path{i};
        reghemcell = split(reghemstr,'_');
        reglog = ismember(CCF_labels(:,1),reghemcell{1});
        if strcmp(reghemcell{2},'L')
            hemlog = ismember(CCF_labels(:,4),'Left Hemisphere'); 
        elseif strcmp(reghemcell{2},'R')
            hemlog = ismember(CCF_labels(:,4),'Right Hemisphere'); 
        else 
            hemlog = ones(size(Conn,1),1);
        end
        init_path((reglog + hemlog) == 2) = 1;
    end
elseif isnan(mousedata_struct.(ip.Results.study).seed)
    init_path = logical(mousedata_struct.(ip.Results.study).data(:,1));
    init_path = DataToCCF(init_path,ip.Results.study,matdir);
else
    init_path = logical(mousedata_struct.(ip.Results.study).seed);
    init_path = DataToCCF(init_path,ip.Results.study,matdir);
end
beta_new = ip.Results.beta * ip.Results.time_scale;
gamma1_new = ip.Results.gamma1 * ip.Results.time_scale;
gamma2_new = ip.Results.gamma2 * ip.Results.time_scale;
taufun = @(x) ip.Results.init_rescale - (x +...
    (gamma1_new * x.^2)./(beta_new - gamma2_new * x));
options_taufun = optimset('TolFun',ip.Results.fsolvetol,'Display','off');
init_rescale_n = fsolve(taufun,0,options_taufun);
init_tau = init_rescale_n * init_path;

switch ip.Results.connectome_subset
    case 'Hippocampus'
        inds = ismember(CCF_labels(:,3),'Hippocampus');
    case 'Hippocampus+PC+RSP'
        inds_hipp = ismember(CCF_labels(:,3),'Hippocampus');
        inds_pc = ismember(CCF_labels(:,1),'Piriform area');
        inds_rsp = ismember(CCF_labels(:,3),'Retrosplenial Area');
        inds = logical(inds_hipp + inds_pc + inds_rsp);
    case 'RH'
        inds = ismember(CCF_labels(:,4),'Right Hemisphere');
    case 'LH'
        inds = ismember(CCF_labels(:,4),'Left Hemisphere');
    otherwise
        inds = logical(ones(size(Conn,1),1)); %#ok<LOGL> 
end
CCF_labels(inds,1)
find(inds)
Adj = Adj(inds,inds);
Conn = Conn(inds,inds);
Vol = DefaultAtlas.volumes(inds);
init_tau = init_tau(inds);
nroi = size(Adj,1);
% i_nonzero_init_tau = init_tau > 0;
% i_zero = ~i_nonzero_init_tau;

t = 0:ip.Results.dt:ip.Results.T;
nt = length(t);
N = zeros(nroi,nt); 
N(:,1) = init_tau(:);
netw_flux= zeros([nroi,size(N)]);
W_1=zeros([nroi,size(N)]);
W_2=zeros([nroi,size(N)]);
S_ss=zeros([nroi,size(N)]);
R_ss=zeros([nroi,size(N)]);
N_adj_0 = N(:,1) .* Adj;
fprintf('Calculating initial flux\n')
netw_flux(:,:,1) = NetworkFluxCalculator(N_adj_0,N(:,1),matdir,...
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
                                'connectome_subset',ip.Results.connectome_subset,...
                                'time_scale',ip.Results.time_scale,...
                                'len_scale',ip.Results.len_scale); 

for h = 1:(nt-1)
    fprintf('Time step %d/%d\n',h,nt-1)
    N_adj_in = N(:,h) .* Adj; %incoming edges for columns
    fprintf('Calculating weights\n')
    [W_1(:,:,h),W_2(:,:,h),R_ss(:,:,h),S_ss(:,:,h)] =...
                                    MassBalanceCalculator(netw_flux(:,:,h),...
                                    N_adj_in,N(:,h),matdir,...
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
                                    'fsolvetol',ip.Results.fsolvetol,...
                                    'connectome_subset',ip.Results.connectome_subset,...
                                    'time_scale',ip.Results.time_scale,...
                                    'len_scale',ip.Results.len_scale); 
    V_ss_1_h = S_ss(:,:,h);
    V_ss_2_h = R_ss(:,:,h);
    v_1 = diag(Conn*V_ss_1_h.'); % the contributions of the mass change's values in V_ss_1 on each node are by rows 
    v_2 = diag(Conn.'*V_ss_2_h); % the contributions of the mass change's values in V_ss_2 on eanch node are by columns
    F_in = 6*30*24*(60)^2 * netw_flux(:,:,h); % convert from seconds to 180 days
    F_out = F_in.';
    m_t= gamma1_new.*N(:,h).*(2*beta_new-gamma2_new.*...
        N(:,h)./(beta_new-gamma2_new.*N(:,h)).^2); % 2.5e-5
    N(:,h+1) = N(:,h)+(1./(Vol.*(1+m_t)+v_1+v_2)).*((diag((Conn.'*F_in)) - diag((Conn*F_out)))*ip.Results.dt); 
    % ((diag((Conn.'*F_in)) - diag((Conn*F_out)))*k + beta*m(:,h)*k-gamma1*(n(:,h).*n(:,h))*k-gamma2*(n(:,h).*m(:,h))*k);
    N_adj_h1 = N(:,h+1).*Adj;

%     Gamma_h = ip.Results.beta * m(:,h) * ip.Results.dt - ...
%         ip.Results.gamma1 * (n(:,h).*n(:,h) + n(:,h).*m(:,h)) * ip.Results.dt;
%     n(:,h+1) = n(:,h) + (diag((Conn.' * F_in)) - diag((Conn * F_out))) * ip.Results.dt...
%         + Gamma_h;
%     m(:,h+1) = m(:,h) - Gamma_h;
    fprintf('Calculating new flux\n')
    netw_flux(:,:,h+1) = NetworkFluxCalculator(N_adj_h1,N(:,h+1),matdir,...
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
                                    'fsolvetol',ip.Results.fsolvetol,...
                                    'connectome_subset',ip.Results.connectome_subset,...
                                    'time_scale',ip.Results.time_scale,...
                                    'len_scale',ip.Results.len_scale);
end
M = (gamma1_new * N.^2)./(beta_new - gamma2_new * N);
% mass_cons_check = 0;
% if mass_cons_check
%     Mass_tot = sum(N,1)+sum(M,1);
%     fprintf('Total Node Mass = %d\n',Mass_tot)
% %     Mass_edge = NaN(length(Mass_tot),1);
% %     for h=1:(nt-1)
% %         Mass_edge(h) = sum(Mass_tot_edge(:,:,h),'all');
% %     end
% %     fprintf('Total Edge Mass = %d\n',Mass_edge)
% end

model_outputs = struct;
model_outputs.Predicted.N = N;
model_outputs.Predicted.M = M;
model_outputs.Predicted.F = netw_flux;
model_outputs.Parameters.beta = ip.Results.beta;
model_outputs.Parameters.gamma1 = ip.Results.gamma1;
model_outputs.Parameters.gamma2 = ip.Results.gamma2;
model_outputs.Parameters.delta = ip.Results.delta;
model_outputs.Parameters.epsilon = ip.Results.epsilon;
model_outputs.Parameters.lambda1 = ip.Results.lambda1;
model_outputs.Parameters.lambda2 = ip.Results.lambda2;
model_outputs.Parameters.frac = ip.Results.frac;
model_outputs.Sim.L1 = ip.Results.L1;
model_outputs.Sim.L2 = ip.Results.L2;
model_outputs.Sim.L_int = ip.Results.L_int;
model_outputs.Sim.L_ais = ip.Results.L_ais;
model_outputs.Sim.L_syn = ip.Results.L_syn;
model_outputs.Sim.dt = ip.Results.dt;
model_outputs.Sim.T = ip.Results.T;
model_outputs.Sim.len_scale = ip.Results.len_scale;
model_outputs.Sim.time_scale = ip.Results.time_scale;
model_outputs.Sim.connectome_subset = ip.Results.connectome_subset;
model_outputs.Sim.study = ip.Results.study;
model_outputs.Sim.init_rescale = ip.Results.init_rescale;
model_outputs.Sim.init_path = init_tau;
model_outputs.Sim.resmesh = ip.Results.resmesh;
model_outputs.Sim.rel_tol = ip.Results.reltol;
model_outputs.Sim.abs_tol = ip.Results.abstol;
model_outputs.Sim.fsolve_tol = ip.Results.fsolvetol;

if ip.Results.plotting
    figure
    subplot(2,1,1)
    plot(t,N);
    xlabel('t');
    ylabel('N(t)');
    title("N,M distributions on the network",'Fontsize',12);
    txt = ['$\mathbf{\lambda_1} = $' num2str(ip.Results.lambda1) ','...
        '$\mathbf{\lambda_2} = $' num2str(ip.Results.lambda2) ',' , ......
        '$\mathbf{\beta} = $' num2str(ip.Results.beta) ',' , ...
        '$\mathbf{\gamma1}=$' num2str(ip.Results.gamma1) ',',...
        '$\mathbf{\epsilon}=$' num2str(ip.Results.epsilon) ',',...
        '$\mathbf{\delta}=$' num2str(ip.Results.delta)];
    subtitle(txt,'Interpreter','latex');
    subplot(2,1,2)
    plot(t,M);
    xlabel('t');
    ylabel('M(t)');

%     figure
%     subplot(2,1,1)
%     plot(t,N(i_nonzero_init_tau,:));
%     subplot(2,1,2)
%     plot(t,M(i_nonzero_init_tau,:));
end

end


% tmax=0.05;
% trange=[0,tmax];
% n_tau0=init_tau.';
% m_tau0=zeros(1,nroi);
% xtau0=[n_tau0(:); m_tau0(:)];
% size(xtau0)
% odeopts = odeset('NonNegative', 1:length(xtau0), 'RelTol',1e-3,'AbsTol',1e-3);
% [t,tau_sol]=ode23s(@ODE_model_BFRTT, trange, xtau0, odeopts);
% nt = size(tau_sol,1)
%      
%   n = tau_sol(:, 1:nroi); 
%    n = reshape(n, [nt, nroi]);
%    n = permute(n, [2,1]);
% m = tau_sol(:, nroi+1:2*nroi); 
%    m = reshape(m, [nt, nroi]);
%    m = permute(m, [2,1]);

% Runge-Kutta method (2th order)
% for h=1: 2%(nt-1)
% 
%  
%  n_adj_in=n(:,h).*Adj; %incoming edges %for columns
% 
% 
%   A_in(:,:,h)=flux_calculator11_15(n_adj_in,n(:,h));
%  
%  F_in=A_in(:,:,h);
% F_out=F_in.';
% 
%  n1=(diag((Conn.'*F_in)) - diag((Conn*F_out))) + beta*m(:,h)-gamma*(n(:,h).*n(:,h)+n(:,h).*m(:,h));
%  m1=-beta*m(:,h)+gamma*(n(:,h).*(n(:,h)+m(:,h)));
%  n_adj_in_1=(n(:,h)+k*n1).*Adj; %incoming edges %for columns
%  A_in_1=flux_calculator11_15(n_adj_in_1,(n(:,h)+k*n1));
%  
%  F_in_1=A_in_1;
%  F_out_1=F_in_1.';
%  n2=(diag((Conn.'*F_in_1)) - diag((Conn*F_out_1)))+ beta*(m(:,h)+k*n1)-gamma*((n(:,h)+k*n1).*(n(:,h)+k*n1)+(n(:,h)+k*n1).*(m(:,h)+k*n1));
%  m2=-beta*(m(:,h)+k*m1)+gamma*((n(:,h)+k*m1).*(n(:,h)+k*m1)+(n(:,h)+k*m1).*(m(:,h)+k*m1));
% 
%  n(:,h+1)=n(:,h)+k*0.5*(n1+n2);
%  m(:,h+1)=m(:,h)+k*0.5*(m1+m2);
% end




%     function tau_prime=ODE_model_BFRTT(t,x)
% 
%            n=x(1:nroi);
% 
%            m=x(nroi+1:2*nroi) ;  
%         
%            %nprime=zeros(size(n));
%           % mprime=zeros(size(m));
%            n_adj_in=n.*Adj; %incoming edges %for columns
% 
%        
% 
%            A_in=flux_calculator11_15(n_adj_in,n);
% 
%             A_out=A_in.';
%            nprime=(diag((Conn.'*A_in)) - diag(Conn*A_out))+ beta*m-gamma*(n.*n+n.*m);
%           mprime=-beta*m+gamma*(n.*(n+m));
% %              nprime=(diag((Conn.'*A_in)) - diag(Conn*A_out'))+ beta*m-gamma*(n.*n+n.*m);
% %              mprime=-beta*m+gamma*(n.*(n+m));
%    tau_prime=[nprime,mprime];
%    tau_prime=tau_prime(:);
%     end




%  
   %discrete model implementation
%      for h=1: (nt-1)
%        
%         for i=1:nroi
%        
%           x_adj_in=Adj(:,i)>0;
%           n_ad_in=n(x_adj_in,h);
%           x_adj_out=Adj(i,:)>0;
%           n_ad_out=n(x_adj_out,h);
%           A_in(x_adj_in,i)=flux_calculator(n_ad_in,n(i,h)); %incoming flux in the vertex i
%           A_out(i,x_adj_out)=flux_calculator(n(i,h),n_ad_out); %outcoming flux in the vertex i
%          
%          
%         end
%          n(:,h+1)=n(:,h)+(diag((Conn.'*A_in)) - diag(Conn*A_out'))*k+ beta*m(:,h)*k-gamma*(n(:,h).*n(:,h)+n(:,h).*m(:,h))*k;
%          m(:,h+1)=m(:,h)-beta*m(:,h)*k+gamma*(n(:,h).*(n(:,h)+m(:,h)))*k;
%       end
                
%ode approach
%     function tau_prime=ODE_model_BFRTT(x,t)
% 
%            n=x(1:nroi);
%            m=x(nroi+1:2*nroi);
%            A_in= zeros(nroi);
%            A_out=zeros(nroi);
%            n_adj_in=Adj.*n;
%            
%            n_adj_out=(Adj.').*n;
%            for i=1:nroi
%                 x_adj_in=Adj(:,i)>0;
%                
%                x_adj_out=Adj(i,:)>0;
%                n_adj=n_adj(n_adj_in(i));
%              A_in(x_adj_in,i)=flux_calculator(n_adj(n_adj>0),n(i)); %incoming flux in the vertex i
%            
%                n_adj_1=n_adj_out(i);
%              A_out(i,x_adj_out)=flux_calculator(n(i),n_adj_1(n_adj_1>0));
%            end
%              nprime=(diag((Conn.'*A_in)) - diag(Conn*A_out'))+ beta*m-gamma*(n.*n+n.*m);
%              mprime=-beta*m+gamma*(n.*(n+m));
%              tau_prime=[nprime,mprime];
%     end


%          x_adj_u=Adj_u(:,i)>0;
%          n_ad_u=n(x_adj_u,1);
         %Conn_ad=Conn(x_adj,i);
%       %Conn_ad=Conn(i,n_ad~=0);
%       n_ad=n(x_adj,h);
%       Conn_ad=Conn(i,x_adj);
%       if n(i,h)==0
%           n_ad=n_ad(n_ad~=0);
%          Conn_ad=Conn(i,n_ad~=0);
%       end
        
       %  A_u(x_adj_u,i)=flux_calculator(n(i,1),n_ad_u.');
         %n(i,h+1)=n(i,h)+A*Conn_ad*k;
        %end
        %N=(diag((Conn.'*Adj_l+Conn.'*Adj_u)) - diag(Conn*Adj_l.'+Conn*Adj_u.'))%*k;

       
       
       
       
       %end

%     flux=flux+Conn(i,j)*A(i,j,h)*k;
%   %for i=1:nroi
%     for j=1:nroi
%       if (C(i,j)~=0) && (abs(n(i,h))>=1e-06 && abs(n(i,h))>=1e-6);
%       
%           A(i,j,h)=flux_calculator(n(i,h),n(j,h));
%           flux=flux+Conn(i,j)*A(i,j,h)*k;
%       end
%     end
%   n(i,h+1)=n(i,h)+flux;
%   %end


%[M,i]=min(abs((n(:,2)-3.744e-05)),[],'all')
% [M,i]=min(n,[],'all')
% N=length(n(n>0))

%end