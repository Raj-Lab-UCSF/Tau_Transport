function model_outputs = NetworkTransportModel(matdir,varargin)
if nargin < 1
    matdir = [cd filesep 'FromVeronica'];
end

load([matdir filesep '220919_MouseDataforVeronica.mat'],'baseline_pathology'); 
Conn = readmatrix([matdir filesep 'mouse_connectome_19_01.csv']);
nroi = size(Conn,1);
Adj = readmatrix([matdir filesep 'mouse_adj_matrix_19_01.csv']);

beta_ = 1e-06;
gamma1_ = 2e-05;
gamma2_ = 2e-06;
delta_ = 1;
epsilon_ = 0.01;
lambda1_ = 0.01;
lambda2_ = 0.01;
study_ = 'Hurtado';
init_path_ = [];
init_rescale_ = 1e-4;
dt_ = 0.005;
T_ = 0.05;
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

ip = inputParser;
% validChar = @(x) ischar(x);
validScalar = @(x) isnumeric(x) && isscalar(x) && (x>=0);
validNonnegative = @(x) isnumeric(x) && (sum(x>=0) == length(x));
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

addParameter(ip, 'study', study_);
addParameter(ip, 'init_rescale', init_rescale_, validScalar);
addParameter(ip, 'dt', dt_, validScalar);
addParameter(ip, 'T', T_, validScalar);
addParameter(ip, 'init_path', init_path_, validNonnegative);
addParameter(ip, 'plotting', plotting_, validLogical);
parse(ip, varargin{:});

if strcmp(ip.Results.study,'custom') && ~isempty(ip.Results.init_path)
    init_tau = ip.Results.init_rescale * ip.Results.init_path;
else
    init_tau = ip.Results.init_rescale * baseline_pathology.(ip.Results.study);
end

i_nonzero_init_tau = init_tau > 0;
t = 0:ip.Results.dt:ip.Results.T;
nt = length(t);
n = zeros(nroi,nt);
m = zeros(nroi,nt);
n(:,1) = init_tau(:);
F_in_mat = zeros([nroi,size(n)]);
F_out_mat = F_in_mat;

for h = 1:(nt-1)
    fprintf('Time step %d/%d\n',h,nt-1)
    tic
    n_adj_in = n(:,h) .* Adj; %incoming edges for columns
    F_in = NetworkFluxCalculator(n_adj_in, n(:,h), matdir,...
                                    'beta',ip.Results.beta,...
                                    'gamma1',ip.Results.gamma1,...
                                    'gamma2',ip.Results.gamma2,...
                                    'delta',ip.Results.delta,...
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
                                    'fsolvetol',ip.Results.fsolvetol);
    F_out = F_in.';
    Gamma_h = ip.Results.beta * m(:,h) * ip.Results.dt - ...
        ip.Results.gamma1 * (n(:,h).*n(:,h) + n(:,h).*m(:,h)) * ip.Results.dt;
    n(:,h+1) = n(:,h) + (diag((Conn.' * F_in)) - diag((Conn * F_out))) * ip.Results.dt...
        + Gamma_h;
    m(:,h+1) = m(:,h) - Gamma_h;
    F_in_mat(:,:,h) = F_in;
    F_out_mat(:,:,h) = F_out;
    toc
end

model_outputs = struct;
model_outputs.Predicted.n = n;
model_outputs.Predicted.m = m;
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
    plot(t,n);
    xlabel('t');
    ylabel('n(t)');
    title("n,m distributions on the network",'Fontsize',12);
    txt = ['$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $'...
     num2str(lambda2_) ',' , '$\mathbf{\epsilon}=$' num2str(epsilon_) ',',...
     '$\mathbf{\delta}=$' num2str(delta_) ];
    subtitle(txt,'Interpreter','latex');
    subplot(2,1,2)
    plot(t,m);
    xlabel('t');
    ylabel('m(t)');

    figure
    subplot(2,1,1)
    plot(t,n(i_nonzero_init_tau,:));
    subplot(2,1,2)
    plot(t,m(i_nonzero_init_tau,:));
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