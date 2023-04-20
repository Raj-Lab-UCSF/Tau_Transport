function [N,M] = ode_modelBFRTT03_09(varargin)

%load("220919_MouseDataforVeronica.mat");
load 220919_MouseDataforVeronica.mat  injection_sites; 
load DefaultAtlas.mat DefaultAtlas;
Conn=readmatrix('mouse_connectome_19_01.csv');
% Row_weigh_degree=mean(Conn,2);
% Col_weigh_degree=mean(Conn,1);
% Conn_norm_row=(1./Row_weigh_degree).*Conn;
% Conn_norm_col=(1./Col_weigh_degree).*Conn;
% Conn_norm_col=Conn_norm_col.';
%Conn=Conn(1:213,1:213);
Conn=Conn([27:37 27+213:37+213], [27:37 27+213:37+213]);
% Conn_norm_row=Conn_norm_row([27:37 27+213:37+213], [27:37 27+213:37+213]);
% Conn_norm_col=Conn_norm_col([27:37 27+213:37+213], [27:37 27+213:37+213]);
%Conn=Conn./max(Conn,[],'all');
nroi=size(Conn,1);
Adj=readmatrix('mouse_adj_matrix_19_01.csv');
Adj=Adj([27:37 27+213:37+213], [27:37 27+213:37+213]);
Vol=DefaultAtlas.volumes([27:37 27+213:37+213]);
%Vol=Vol./mean(Vol,1);

%Vol=volumes
%Adj=Adj(1:213, 1:213);
time_scale=6*30*24*(60)^2;
beta_ =1e-06; %1e-02;
gamma1_= 2e-05;%2e-01;
gamma2_=0; %2e-05*time_scale;
delta_ = 1; %the parameters delta,epsilon,lambda1, lambda2 are rewritten here just for plotting
epsilon_ = 0.01;
lambda1_ = 0.01;
lambda2_ = 0.01;%0.01;
init_rescale_ = 1e-2;
dt_ = 0.01;  %0.5e5;%24*(60)^2;
T_ = 0.1; %5*dt_;
frac_ = 0.92; % Average fraction of n diffusing (Konsack 2007)
L_int_ = 1000; % in micrometers
L1_ = 200;
L2_ = 200; 
L_ais_ = 40;
L_syn_ = 40;
resmesh_ = 'coarse';

reltol_ = 1e-4;
abstol_ = 1e-4;
fsolvetol_ = 1e-6;

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
addParameter(ip, 'init_rescale', init_rescale_, validScalar);
addParameter(ip, 'dt', dt_, validScalar);
addParameter(ip, 'T', T_, validScalar);
parse(ip, varargin{:});






initial_tau=injection_sites.Hurtado;
%initial_tau=initial_tau(1:213);
%initial_tau([36 36+213]);
initial_tau([37, 37+213])=0; % initiation of the pathology in region 36 and 36+213
initial_tau=initial_tau([27:37 27+213:37+213]);
init_tau=ip.Results.init_rescale*initial_tau;
i_nonzero_init_tau=init_tau>0;
i_zero=init_tau==0;

  t = 0:ip.Results.dt:ip.Results.T;
  nt=length(t);
  N=zeros(nroi,nt);
  %m=zeros(nroi,nt);
  N(:,1)=init_tau(:);
 netw_flux= zeros([nroi,size(N)]);
 %netw_flux_1=zeros([nroi,size(N)]);
 W_1=zeros([nroi,size(N)]);
 W_2=zeros([nroi,size(N)]);
 S_ss=zeros([nroi,size(N)]);
 R_ss=zeros([nroi,size(N)]);
 %Mass_tot_edge=zeros([nroi,nroi,nt-1]);
N_adj_0=N(:,1).*Adj;
 netw_flux(:,:,1)=NetworkFluxCalculator_veronica(N_adj_0,N(:,1), 'beta',ip.Results.beta,...
                                    'gamma1',ip.Results.gamma1,...
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
                                    'fsolvetol',ip.Results.fsolvetol); %compute the steady state  network flux at time t0

 
% netw_flux_1(:,:,1)=netw_flux(:,:,1);

for h=1:(nt-1)
  %Ad_in=logical(Adj);
  fprintf('Time step %d/%d\n',h,nt-1) 
 tic
 %incoming edges %for columns
%     %C_1=c(j);
%     C_1=repmat(n(:,h),1,length(Ad_in));
%    C_1=C_1.';
%      i_app=logical(((n_adj_in>0)+(C_1>0)).*(Ad_in));
% m(:,h)=(ip.Results.gamma1 * n(:,h).^2)./(ip.Results.beta-ip.Results.gamma2 * n(:,h));
  N_adj_in=N(:,h).*Adj;

 %[W_1(:,:,h),W_2(:,:,h),R_ss(:,:,h),S_ss(:,:,h),Mass_tot_edge(:,:,h)]=mass_balance_network_03_09(netw_flux(:,:,h),N_adj_in,N(:,h));

    fprintf('Calculating weights\n')
 [W_1(:,:,h),W_2(:,:,h),R_ss(:,:,h),S_ss(:,:,h)]=mass_balance_network_03_09(netw_flux(:,:,h),N_adj_in,N(:,h));
 
 V_ss_1_h=S_ss(:,:,h);
 V_ss_2_h=R_ss(:,:,h);
 %v_1=diag(V_ss_1_h*(Conn.')); %the contributions of the mass change's values in S_ss on each node are by rows weighted by the Conn matrix
 %v_2= diag(Conn.'*V_ss_2_h);%the contributions of the mass change's values in R_ss on eanch node are by columns
 v_1=diag(Conn*V_ss_1_h.');%the contributions of the mass change's values in V_ss_1 on each node are by rows 
 v_2= diag(Conn.'*V_ss_2_h);%the contributions of the mass change's values in V_ss_2 on eanch node are by columns
 %size(v_1)
 %size(v_2)
 m_t= ip.Results.gamma1.*N(:,h).*(2*ip.Results.beta-ip.Results.gamma2.*N(:,h)./(ip.Results.beta-ip.Results.gamma2.*N(:,h)).^2); %2.5e-5
 %v_1=sum(V_ss_1_h,2);
 %v_2=sum(V_ss_2_h,1).';
%  if h==1
%     netw_flux(:,:,h)=time_scale.*netw_flux(:,:,h);
%      netw_flux_1(:,:,1)=netw_flux(:,:,1);
%  end
 F_in=time_scale.*netw_flux(:,:,h);
 F_out=F_in.';
%  k_2=1./(1+v_1+v_2);
%  k_2(i_nonzero_init_tau)
 %Gamma=((diag((Conn.'*F_in)) - diag((Conn*F_out)))*k + beta*m(:,h)*k-gamma1*(n(:,h).*n(:,h))*k-gamma2*(n(:,h).*m(:,h))*k);
 %Gamma_h =ip.Results.beta * m(:,h) * ip.Results.dt -(ip.Results.gamma1 * n(:,h).*n(:,h)+ip.Results.gamma2 * n(:,h).*m(:,h)) * ip.Results.dt;
 %Gamma_h(i_nonzero_init_tau)
 N(:,h+1)=N(:,h)+(1./(Vol.*(1+m_t) +v_1+v_2)).*( (diag((Conn.'*F_in)) - diag((Conn*F_out)) )*ip.Results.dt    )  ; %((diag((Conn.'*F_in)) - diag((Conn*F_out)))*k + beta*m(:,h)*k-gamma1*(n(:,h).*n(:,h))*k-gamma2*(n(:,h).*m(:,h))*k);
 %n(:,h+1)=n(:,h)+(1./(1+v)).*( (diag((Conn_norm_col*F_in)) - diag((Conn_norm_row*F_out)) )*ip.Results.dt    )  ; 
 % m(:,h+1)=m(:,h)-Gamma_h;
  N_adj_h1=N(:,h+1).*Adj;
 netw_flux(:,:,h+1)=NetworkFluxCalculator_veronica(N_adj_h1,N(:,h+1), 'beta',ip.Results.beta,...
                                    'gamma1',ip.Results.gamma1,...
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

% netw_flux_1(:,:,h+1)=time_scale.*netw_flux_1(:,:,h+1);
%  alpha_2=N(:,h+1)-N(:,h);
%  alpha_1=alpha_2.';
%  netw_flux(:,:,h+1)=netw_flux(:,:,h)+1/(time_scale).*(alpha_1.*W_1(:,:,h)+alpha_2.*W_2(:,:,h));
 toc
end
%max(abs(netw_flux-netw_flux_1),[],'all')
 M=(ip.Results.gamma1 * N.^2)./(ip.Results.beta-ip.Results.gamma2 * N);

% Mass_tot=sum(N,1)+sum(M,1)
% for h=1:(nt-1)
% Mass_edge(h)=sum(Mass_tot_edge(:,:,h),'all');
% end
% Mass_edge
%  %  
%  if isempty(n(:,h+1)<0)==0
%       h_1=ip.Results.dt/10;
%       t_1 = t(h):h_1:t(h+1);
%   nt_1=length(t_1);
%   n_1=zeros(nroi,nt_1);
%   n_1(:,1)=n(:,h);
%  netw_flux_app= zeros([nroi,size(n_1)]);
%  netw_flux_app(:,:,1)=netw_flux(:,:,h);
%  W_1_app=zeros([nroi,size(n_1)]);
%  W_2_app=zeros([nroi,size(n_1)]);
%  V_1_app=zeros([nroi,size(n_1)]);
%  V_2_app=zeros([nroi,size(n_1)]);
%  for k=1:nt_1
%        n_adj_in_1=n_1(:,k).*Adj;
%        [W_1_app(:,:,k),W_2_app(:,:,k),V_1_app(:,:,k),V_2_app(:,:,k)]=mass_balance_network_03_09(netw_flux_app(:,:,k),n_adj_in_1,n_1(:,k));
%  
%  
%       V_ss_1_k=V_1_app(:,:,k);
%       V_ss_2_k=V_2_app(:,:,k);
%       
%  v_1_app=diag(Conn_norm_row*V_ss_1_k.');%the contributions of the mass change's values in V_ss_1 on each node are by rows 
%  v_2_app= diag(Conn_norm_col*V_ss_2_k);%the contributions of the mass change's values in V_ss_2 on eanch node are by columns
%  
%  v_app=1e-4*(1./Vol).*(v_1_app+v_2_app); %2.5e-5
%  %v_1=sum(V_ss_1_h,2);
%  %v_2=sum(V_ss_2_h,1).';
%  F_in=netw_flux_app(:,:,k);
%  F_out=F_in.';
% %  k_2=1./(1+v_1+v_2);
% %  k_2(i_nonzero_init_tau)
%  %Gamma=((diag((Conn.'*F_in)) - diag((Conn*F_out)))*k + beta*m(:,h)*k-gamma1*(n(:,h).*n(:,h))*k-gamma2*(n(:,h).*m(:,h))*k);
%  Gamma_h = 0; %ip.Results.beta * m(:,h) * ip.Results.dt -(ip.Results.gamma1 * n(:,h).*n(:,h)+ip.Results.gamma2 * n(:,h).*m(:,h)) * ip.Results.dt;
%  %Gamma_h(i_nonzero_init_tau)
% % n(:,h+1)=n(:,h)+(1./(1+v)).*( (diag((Conn.'*F_in)) - diag((Conn*F_out)) )*ip.Results.dt +Gamma_h    )  ; %((diag((Conn.'*F_in)) - diag((Conn*F_out)))*k + beta*m(:,h)*k-gamma1*(n(:,h).*n(:,h))*k-gamma2*(n(:,h).*m(:,h))*k);
%  n_1(:,k+1)=n_1(:,k)+(1./(1+v_app)).*( (diag((Conn_norm_col*F_in)) - diag((Conn_norm_row*F_out)) )*h_1 +Gamma_h    )  ; 
%   alpha_2_app=n_1(:,k+1)-n_1(:,k);
%  alpha_1_app=alpha_2_app.';
%  netw_flux_app(:,:,k+1)=netw_flux_app(:,:,k)+alpha_1_app.*W_1_app(:,:,k)+alpha_2_app.*W_2_app(:,:,k);
%  end
% 
%  n(:,h+1)=n_1(:,end);
%  netw_flux(:,:,h+1)= netw_flux_app(:,:,end);
% end




 %n(i_nonzero_init_tau,h+1)
 %(1./(1+v_1+v_2)).
 
 %n(:,h+1)=n(:,h)+ beta*m(:,h)*k-gamma1*(n(:,h).*n(:,h))*k-gamma2*(n(:,h).*m(:,h))*k;

%n(:,h+1)=max(n(:,h+1),1e-8*ones(nroi,1));

 %m(:,h+1)=m(:,h)-beta*m(:,h)*k+gamma1*(n(:,h).*n(:,h))*k+gamma2*(n(:,h).*m(:,h))*k;










figure
subplot(2,1,1)
plot(t,N);


xlabel('t');
ylabel('N(t)');
 title("N,M distributions on the network",'Fontsize',12);
%  txt = ['$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' ,'$\mathbf{\gamma}=$' num2str(gamma1) ',', '$\mathbf{\beta}=$' num2str(beta) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
%  subtitle(txt,'Interpreter','latex');
subplot(2,1,2)
plot(t,M);
xlabel('t');
ylabel('M(t)');


%[M,i]=min(abs((n(:,2)-3.744e-05)),[],'all')
% [Min,i]=min(N,[],'all')
% N_nonzero=length(N(N>0))
% 
% figure
% subplot(2,1,1)
% plot(t,N(i_nonzero_init_tau,:));
%  title("N,M distributions on the network",'Fontsize',12);
%  %txt = ['$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' ,'$\mathbf{\gamma}=$' num2str(gamma1) ',', '$\mathbf{\beta}=$' num2str(beta) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
%  %subtitle(txt,'Interpreter','latex');
% subplot(2,1,2)
% plot(t,M(i_nonzero_init_tau,:));
% 
% figure
% subplot(2,1,1)
% plot(t,N(i_zero,:));
%  title("N,M distributions on the network",'Fontsize',12);
%  %txt = ['$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' ,'$\mathbf{\gamma}=$' num2str(gamma1) ',', '$\mathbf{\beta}=$' num2str(beta) ','  '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
% % subtitle(txt,'Interpreter','latex');
% subplot(2,1,2)
% plot(t,M(i_zero,:));

end
