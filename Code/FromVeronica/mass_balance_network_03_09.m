%function A_1=flux_calculator(b,c,flag,varargin)
%b array dim 1*i c matrix dim i*j
%function [W_1_flux,W_2_flux,R_ss,S_ss,mass_tot_edge ]=mass_balance_network_03_09(network_flux,tau_x0,tau_xL,varargin)
function [W_1_flux,W_2_flux,R_ss,S_ss]=mass_balance_network_03_09(network_flux,tau_x0,tau_xL,varargin)


alpha_ = 0; % Not currently explored in the model
len_scale=1e-3; %scale factor from um to mm

beta_ = 1e-06;
gamma1_ =2e-05;
gamma2_=0;
delta_ =1;
epsilon_ =  0.01;
lambda1_ = 0.01;
lambda2_ = 0.01;
frac_ = 0.92; % Average fraction of n diffusing (Konsack 2007)
L_int_ = 1000; % in micrometers
L1_ = 200;
L2_ = 200; 
L_ais_ = 40;
L_syn_ = 40;

resmesh_ = 'coarse';

ip = inputParser;
validScalar = @(x) isnumeric(x) && isscalar(x) && (x>=0);
addParameter(ip, 'alpha', alpha_, validScalar);
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
%addParameter(ip, 'total_mass',total_mass_);
parse(ip, varargin{:});
% % % 2. Definition of constants

v_a = 0.7*len_scale; % Average velocity (mm/s) of anterograde active transpot (Konsack 2007)
v_r = 0.7*len_scale; % Average velocity (mm/s) of retrograde active transport (Konsack 2007)
diff_n =12*len_scale^2; %12; % Diffusivity (mm^2/s) of n (Konsack 2007)
% diff_ratio = 0; % Diffusivity of m relative to n (0.01 is a heuristic)
% % % 3. Definition of the (inhomogeneous) xmesh
L_total = ip.Results.L1 + ip.Results.L2 + ip.Results.L_int; % size of the system
if strcmp(ip.Results.resmesh, 'fine')
    num_comp = 1000; % number of xmesh points
    num_ext = 100; % number of GM compartments per GM region
    num_int = num_comp - 2*(num_ext);
    xmesh1 = [linspace(0,ip.Results.L1-10,num_ext-40),(ip.Results.L1-9.75):0.25:ip.Results.L1];
    xmesh2 = [(ip.Results.L1+ip.Results.L_int+0.25):0.25:(ip.Results.L1+ip.Results.L_int+10),...
        linspace(ip.Results.L1+ip.Results.L_int+10.25,L_total,num_ext-40)];
    xmesh_int = [(ip.Results.L1+0.25):0.25:(ip.Results.L1+ip.Results.L_ais),...
        linspace(ip.Results.L1+ip.Results.L_ais+0.25,ip.Results.L1+ip.Results.L_int-ip.Results.L_syn,num_int-(4*(ip.Results.L_ais+ip.Results.L_syn))),...
        (ip.Results.L1+ip.Results.L_int-(ip.Results.L_syn-0.25)):0.25:(ip.Results.L1+ip.Results.L_int)];
elseif strcmp(ip.Results.resmesh, 'coarse')
    num_comp = 250; % number of xmesh points
    num_ext = 25; % number of compartments per SD region
    num_int = num_comp - 2*(num_ext);
    xmesh1 = [linspace(0,ip.Results.L1-10,num_ext-5),(ip.Results.L1-8):2:ip.Results.L1];
    xmesh2 = [(ip.Results.L1+ip.Results.L_int+2):2:(ip.Results.L1+ip.Results.L_int+10),...
        linspace(ip.Results.L1+ip.Results.L_int+12,L_total,num_ext-5)];
    xmesh_int = [(ip.Results.L1+2):2:(ip.Results.L1+ip.Results.L_ais),...
        linspace(ip.Results.L1+ip.Results.L_ais+2,ip.Results.L1+ip.Results.L_int-ip.Results.L_syn,num_int-(0.5*(ip.Results.L_ais+ip.Results.L_syn))),...
        (ip.Results.L1+ip.Results.L_int-(ip.Results.L_syn-2)):2:(ip.Results.L1+ip.Results.L_int)];
end
xmesh = [xmesh1, xmesh_int, xmesh2];

%Network stuff
Adj=readmatrix('mouse_adj_matrix_19_01.csv');
%Adj=Adj(1:213,1:213);
Adj=Adj([27:37 27+213:37+213], [27:37 27+213:37+213]);
nroi=size(Adj,1);


% Follows the derivation of Michiel Bertsh
 presyn_mask = spatial_mask('presyn');
 xmesh_presyn = len_scale.*xmesh(presyn_mask);
 n0=@(B) B;
 options = odeset('RelTol',1e-4,'AbsTol',1e-4,'NonNegative',1:length(n0));
n_ss_presyn = @(A,B) ode45(@(x,n)ode_ss_n(x,A,n,diff_n),len_scale.*[0,ip.Results.L1],B,options);
n_ss_presyn = @(A,B,x) deval(n_ss_presyn(A,B),x);
 %x0=xmesh_presyn(1);
 %n_ss_presyn = @(A,B,x) max((B - A.*x/diff_n),0); % B = n_ss(0) and is an integration constant
%m_ss_presyn = @(A,B,x) ((ip.Results.gamma*n_ss_presyn(A,B,x).^2)./(ip.Results.beta - n_ss_presyn(A,B,x).*ip.Results.gamma));
%m_ss_presyn = @(A,x) ((ip.Results.gamma*n_ss_presyn(A,x).^2)./(ip.Results.beta));
% % % 5b. Steady state of axon initial segment
 ais_mask = spatial_mask('ais');
 xmesh_ais= len_scale.*xmesh(ais_mask);
 x1 = xmesh_presyn(end);
n_ss_ais = @(A,B) ode45(@(x,n)ode_ss_n(x,A,n,diff_n*ip.Results.lambda1),len_scale.*[ip.Results.L1,ip.Results.L1+ip.Results.L_ais],n_ss_presyn(A,B,x1),options);
n_ss_ais = @(A,B,x) deval(n_ss_ais(A,B),x);
 %n_ss_ais = @(A,B,x)max((B - A.*x1/diff_n - A.*(x-x1)/(diff_n*ip.Results.lambda1)),0);
%n_ss_ais=@(A,B,x)max((n_ss_presyn(A,B,x1)-A*(x-x1)/(diff_n*ip.Results.lambda1)),0);
% m_ss_ais = @(A,B,x) (ip.Results.gamma*n_ss_ais(A,B,x).^2)./(ip.Results.beta - n_ss_ais(A,B,x)*ip.Results.gamma);
%m_ss_ais = @(A,x) ((ip.Results.gamma*n_ss_ais(A,x).^2)./(ip.Results.beta));

% % 5c. Steady state of axon
 axon_mask = spatial_mask('axon');
 xmesh_axon= len_scale.*xmesh(axon_mask);
 x2 = xmesh_ais(end);
%nx_init = @(A_,B)n_ss_ais(A_,B,x2) ;         %max((-A_.*x1/diff_n + B - A_.*(x2 - x1)/(ip.Results.lambda1*diff_n)),0);
%nx_init=@(A_,B)  n_ss_ais(A_,B,x2);
 %options = odeset('RelTol',1e-4,'AbsTol',1e-4,'NonNegative',1:length(nx_init));

 n_ss_axon = @(A,B) ode45(@(x,n)ode_ss_axon(x,A,n),len_scale.*[ip.Results.L1+ip.Results.L_ais,...
           ip.Results.L1+ip.Results.L_int-ip.Results.L_syn],n_ss_ais(A,B,x2),options);

 n_ss_axon= @(A,B,x) deval(n_ss_axon(A,B),x);
%  m_ss_axon =@(A,B,x)(ip.Results.gamma*n_ss_axon(A,B,x).^2)./ ...
%                   (ip.Results.beta - ip.Results.gamma*n_ss_axon(A,B,x));
%     
%m_ss_axon =@(A,x)(ip.Results.gamma*n_ss_axon(A,x).^2)./ ...
               %   (ip.Results.beta);
      

 syncleft_mask = spatial_mask('syncleft');

 xmesh_syncleft = len_scale.*xmesh(syncleft_mask);
 x3 = xmesh_axon(end);
  n_ss_syncleft = @(A,B) ode45(@(x,n)ode_ss_n(x,A,n,(diff_n*ip.Results.lambda2)),len_scale.*[ ip.Results.L1+ip.Results.L_int-ip.Results.L_syn, ip.Results.L1+ip.Results.L_int],n_ss_axon(A,B,x3),options);
n_ss_syncleft = @(A,B,x) deval(n_ss_syncleft(A,B),x);   
 %n_ss_syncleft = @(A,B,x)max((n_ss_axon(A,B,x3) - A.*(x-x3)/(diff_n*ip.Results.lambda2)),0);
    
% m_ss_syncleft =@(A,B,x) M1_0_*x;    
    

% % % 5e. Steady state of postsynaptic somatodendritic compartment
 postsyn_mask = spatial_mask('postsyn');

 xmesh_postsyn = len_scale.*xmesh(postsyn_mask);
 x4 = xmesh_syncleft(end);
% n_ss_postsyn = @(A,B) ode45(@(x,n)ode_ss_n(x,A,n,diff_n),len_scale.*[ip.Results.L1+ip.Results.L_int,L_total],n_ss_syncleft(A,B,x4),options);
%n_ss_postsyn = @(A,B,x) deval(n_ss_postsyn(A,B),x);    
 n_ss_postsyn =@(A,B,x) (n_ss_syncleft(A,B,x4) - A.*x/diff_n);
     
% m_ss_postsyn =@(A,B,x) (ip.Results.gamma*n_ss_postsyn(A,B,x).^2)./ ...
%                 (ip.Results.beta - ip.Results.gamma*n_ss_postsyn(A,B,x));
%m_ss_postsyn =@(A,x) (ip.Results.gamma*n_ss_postsyn(A,x).^2)./ ...
              %  (ip.Results.beta );
       
    
x5=xmesh_postsyn(end);
% %% a.Linearized  equation for v=dn/db or dn/dC in presynaptic Somatodendritic Compartment
 q_in=@(V0) V0;
 options_q = odeset('RelTol',1e-6,'AbsTol',1e-6,'NonNegative',1:length(q_in));
q_ss_presyn = @(W,V0) ode45(@(x,q)ode_ss_n(x,W,q,diff_n),len_scale.*[0,ip.Results.L1],q_in(V0),options_q);
q_ss_presyn = @(W,V0,x) deval(q_ss_presyn(W,V0),x);
 %q_ss_presyn = @(W,V0,x) max((V0- W.*x/diff_n),0);
% % %b.Equation for v=dn/db or dn/dC in axon initial segment
 q_ss_ais = @(W,V0) ode45(@(x,q)ode_ss_n(x,W,q,(diff_n*ip.Results.lambda1)),len_scale.*[ip.Results.L1,ip.Results.L1+ip.Results.L_ais],q_ss_presyn(W,V0,x1),options_q);
q_ss_ais = @(W,V0,x) deval(q_ss_ais(W,V0),x);

%q_ss_ais = @(W,V0,x)max((V0 - W.*x1/diff_n - W.*(x-x1)/(diff_n*ip.Results.lambda1)),0);
%c. Equation for v=dn/db or dn/dC in axon 
% qx_init= @(W_,V0)max((-W_.*x1/diff_n + V0 - W_.*(x2 - x1)/(ip.Results.lambda1*diff_n)),0);
%nx_init=@(A_,B)  n_ss_ais(A_,B,x2);
 %options = odeset('RelTol',1e-4,'AbsTol',1e-4,'NonNegative',1:length(qx_init));

 q_ss_axon = @(W,A,B,V0) ode45(@(x,q)ode_q_axon(x,W,q,A,B),len_scale.*[ip.Results.L1+ip.Results.L_ais,...
           ip.Results.L1+ip.Results.L_int-ip.Results.L_syn],q_ss_ais(W,V0,x2),options_q);

 q_ss_axon= @(W,A,B,V0,x) deval(q_ss_axon(W,A,B,V0),x);
%Equation for v=dn/db or dn/dC in synaptic cleft
 q_ss_syncleft = @(W,A,B,V0) ode45(@(x,q)ode_ss_n(x,W,q,(diff_n*ip.Results.lambda2)),len_scale.*[ ip.Results.L1+ip.Results.L_int-ip.Results.L_syn, ip.Results.L1+ip.Results.L_int],q_ss_axon(W,A,B,V0,x3),options_q);
q_ss_syncleft = @(W,A,B,V0,x) deval(q_ss_syncleft(W,A,B,V0),x); 
%q_ss_syncleft = @(W,A,B,V0,x)max((q_ss_axon(W,A,B,V0,x3) - W.*(x-x3)/(diff_n*ip.Results.lambda2)),0);
%Equation for v=dn/db or dn/dC in post synaptic Somatodendritic Compartment
 q_ss_postsyn =@(W,A,B,V0,x) (q_ss_syncleft(W,A,B,V0,x4) - W.*x/diff_n);
 
 f_q_ss=@(W,A,B,V0,V_L)(q_ss_postsyn(W,A,B,V0,x5)-V_L) ;  %((q_ss_syncleft(W,A,B,V0,x4)- W.*x5/diff_n -V_L));    

W_1_flux=zeros(nroi);
W_2_flux=zeros(nroi);
% V_ss_1=zeros(nroi);
% V_ss_2=zeros(nroi);
A_0=[0;0];
B_0=[0;0];
V_0_0=[1;0];
V_L_0_0=[0;1];
g_0_0=[0;0];
q_flux_0_0=@(W)f_q_ss(W,A_0,B_0,V_0_0,V_L_0_0);
W_1_0=fsolve(q_flux_0_0,g_0_0);

% V_ss_0=[v_ss_presyn(W_1_0,V_0_0,xmesh_presyn) v_ss_ais(W_1_0,V_0_0,xmesh_ais) v_ss_axon(W_1_0,A_0,B_0,V_0_0,xmesh_axon) v_ss_syncleft(W_1_0,A_0,B_0,V_0_0,xmesh_syncleft) v_ss_postsyn(W_1_0,A_0,B_0,V_0_0,xmesh_postsyn)];
% v_ss_int_0=trapz(xmesh,V_ss_0,2);
% W_1_flux= 0.0011*ones(nroi)-diag(0.0011*ones(nroi,1));
% W_2_flux=-0.0011*ones(nroi)-diag(-0.0011*ones(nroi,1));
% V_ss_1=756.4906*ones(nroi) -diag(756.4906*ones(nroi,1));
%  V_ss_2=643.5094*ones(nroi)-diag(643.5094*ones(nroi,1));
 for j=1:nroi
     
     Ad_in=logical(Adj(:,j));
      C_1=tau_xL(j);
    C_1=repmat(C_1,1,length(Ad_in));
   C_1=C_1.';
   i_app_0=logical((tau_x0(:,j)+C_1==0 ).*(Ad_in));
   i_len=ones(nroi,1);
   i_app_00=i_len(i_app_0,1);
   W_1_flux(i_app_0,j)=W_1_0(1)*i_app_00 ;  %ones(length(i_app_00),1);
   W_2_flux(i_app_0,j)=W_1_0(2)*i_app_00;%     %ones(length(i_app_00),1);
%    V_ss_1(i_app_0,j)=v_ss_int_0(1)*i_app_00;%     ones(length(i_app_00),1);
%    V_ss_2(i_app_0,j)=v_ss_int_0(2)*i_app_00;     %ones(length(i_app_00),1);

     i_app=logical(((tau_x0(:,j)>0)+(C_1>0)).*(Ad_in));
    B_1=tau_x0(i_app,j);
   if ~isempty(B_1)>0

      A_1=network_flux(i_app,j);
      % B_1=repmat(B_1,2,1);
% A_1=repmat(A_1,2,1);
% V0_1=[ones(size(B_1)); zeros(size(B_1))];
   V0_1_1=ones(size(B_1));
  V0_1_2=zeros(size(B_1));
% V_L_1=[zeros(size(B_1)); ones(size(B_1))];
  V_L_1_1=zeros(size(B_1));
  V_L_1_2=ones(size(B_1));
 q_flux_1=@(W)f_q_ss(W,A_1,B_1,V0_1_1,V_L_1_1);
q_flux_2=@(W)f_q_ss(W,A_1,B_1,V0_1_2,V_L_1_2);
%       W_1_sample=rand(length(B_1),10); % N random numbers in the interval (a,b) with the formula r = a + (b-a).*rand(N,1)
%       %W_1_sample=[0 W_1_sample];
%       W_2_sample=-1+rand(length(B_1),10); % N random numbers in the interval (a,b) with the formula r = a + (b-a).*rand(N,1)
%       %W_2_sample=[0 W_2_sample];
%       f_test_1=zeros(10,1);
%       f_test_2=zeros(10,1);
%       for i=1:10
%           W_1_test=W_1_sample(:,i);
%          W_2_test=W_2_sample(:,i);
%          f_test_1(i)=norm(f_v_ss(W_1_test,A_1,B_1,V0_1_1,V_L_1_1));
%          f_test_2(i)=norm(f_v_ss(W_2_test,A_1,B_1,V0_1_2,V_L_1_2));
%       end
%       [m_test_1,i_test_1]=min(f_test_1,[],'all');
%        [m_test_2,i_test_2]=min(f_test_2,[],'all');
%       g0_1=W_1_sample(:,i_test_1);
%       g0_2=W_2_sample(:,i_test_2);
%  g0_1=0.1*ones(length(V0_1_1),1);%0.1
%  g0_2=-0.1*ones(length(V0_1_1),1); %-0.1  
 g0_1=zeros(length(V0_1_1),1);
 g0_2=zeros(length(V0_1_1),1);
% v_flux_1=@(W)f_v_ss(W,A_1,B_1,V0_1_1,V_L_1_1);
% v_flux_2=@(W)f_v_ss(W,A_1,B_1,V0_1_2,V_L_1_2);
options=optimset('TolFun',1e-06);
W_1_flux(i_app,j)=fsolve(q_flux_1,g0_1,options);
W_2_flux(i_app,j)=fsolve(q_flux_2,g0_2,options);

%         W_1_flux(i_app,j)=W_1(1:length(B_1));
%         W_2_flux(i_app,j)=W_1(length(B_1)+1 :end);
% W_1=[W_1_flux(i_app,j); W_2_flux(i_app,j)];
% size(W_1)

%Evaluation of v

% 
%  v_ss_presyn_eval= v_ss_presyn(W_1, V0_1,xmesh_presyn);
% 
%  v_ss_ais_eval=v_ss_ais(W_1,V0_1,xmesh_ais);
% 
%   B_1=repmat(B_1,2,1);
%  A_1=repmat(A_1,2,1);
%  v_ss_axon_eval=v_ss_axon(W_1,A_1,B_1,V0_1,xmesh_axon);
% 
%  v_ss_syncleft_eval=v_ss_syncleft(W_1,A_1,B_1,V0_1,xmesh_syncleft);
% 
% 
%  v_ss_postsyn_eval=v_ss_postsyn(W_1,A_1,B_1,V0_1, xmesh_postsyn);
%  V_ss=[v_ss_presyn_eval v_ss_ais_eval v_ss_axon_eval v_ss_syncleft_eval v_ss_postsyn_eval];
%  v_ss_int=trapz(xmesh,V_ss,2);
% % size(v_ss_int)
%   V_ss_1(i_app,j)=v_ss_int(1:(length(B_1)/2)); % the contribution on each node is given by the sum over the rows
%   V_ss_2(i_app,j)=v_ss_int((length(B_1)/2)+1 : end); % the contribution on each node is given by the sum over the columns                

   end
 end
%xmesh= len_scale.*xmesh; 
network_flux=network_flux(:);
tau_x0=tau_x0(:);
 n_ss_presyn_eval= n_ss_presyn(network_flux, tau_x0,xmesh_presyn);

 n_ss_ais_eval=n_ss_ais(network_flux,tau_x0,xmesh_ais);

  
 n_ss_axon_eval=n_ss_axon(network_flux,tau_x0, xmesh_axon);

 n_ss_syncleft_eval=n_ss_syncleft(network_flux,tau_x0,xmesh_syncleft);


 n_ss_postsyn_eval=n_ss_postsyn(network_flux,tau_x0, xmesh_postsyn);
 %n_ss_eval=[n_ss_presyn_eval n_ss_ais_eval n_ss_axon_eval n_ss_syncleft_eval n_ss_postsyn_eval];
 n_ss_eval_1=[n_ss_presyn_eval n_ss_ais_eval n_ss_axon_eval];
%n_m_ss=n_ss_eval_1+(ip.Results.gamma1 * n_ss_eval_1.^2)./(ip.Results.beta-ip.Results.gamma2 *n_ss_eval_1);
% n_m_ss_postsyn=n_ss_postsyn_eval+(ip.Results.gamma1 * n_ss_postsyn_eval.^2)./(ip.Results.beta-ip.Results.gamma2 *n_ss_postsyn_eval);
 
%mass_tot_edge=trapz([xmesh_presyn xmesh_ais xmesh_axon], n_m_ss,2)+trapz(xmesh_syncleft, n_ss_syncleft_eval,2)+trapz(xmesh_postsyn, n_m_ss_postsyn,2);
%mass_tot_edge=reshape(mass_tot_edge,nroi,nroi);

 
 
 n_ss_eval_1=[n_ss_eval_1; n_ss_eval_1];
 %n_ss_eval_2=[n_ss_syncleft_eval; n_ss_syncleft_eval];
 n_ss_eval_2=[n_ss_postsyn_eval; n_ss_postsyn_eval];
 

%  size(n_ss_eval)

W_flux=[W_1_flux(:);W_2_flux(:)];
v0_1=[ones(length(W_flux)/2,1);zeros(length(W_flux)/2,1)];

A_flux=[network_flux(:);network_flux(:)];
b=[tau_x0(:);tau_x0(:)];

 q_ss_presyn_eval_1= q_ss_presyn(W_flux, v0_1,xmesh_presyn);

 q_ss_ais_eval_1=q_ss_ais(W_flux,v0_1,xmesh_ais);

  
 q_ss_axon_eval_1=q_ss_axon(W_flux,A_flux,b,v0_1, xmesh_axon);

 q_ss_syncleft_eval_1=q_ss_syncleft(W_flux,A_flux,b,v0_1,xmesh_syncleft);


 q_ss_postsyn_eval_1=q_ss_postsyn(W_flux,A_flux,b,v0_1, xmesh_postsyn);
 %q_ss=[q_ss_presyn_eval_1 q_ss_ais_eval_1 q_ss_axon_eval_1 q_ss_syncleft_eval_1 q_ss_postsyn_eval_1];
 q_ss_1=[q_ss_presyn_eval_1 q_ss_ais_eval_1 q_ss_axon_eval_1];
 %q_ss_2= q_ss_syncleft_eval_1;
 q_ss_2=q_ss_postsyn_eval_1;
 %size(q_ss)


%v_n_ss=(1+ip.Results.gamma1.*n_ss_eval.*(2*ip.Results.beta-ip.Results.gamma2.*n_ss_eval)./(ip.Results.beta-ip.Results.gamma2.*n_ss_eval).^2).*q_ss;
v_n_ss_1=(1+ip.Results.gamma1.*n_ss_eval_1.*(2*ip.Results.beta-ip.Results.gamma2.*n_ss_eval_1)./(ip.Results.beta-ip.Results.gamma2.*n_ss_eval_1).^2).*q_ss_1;
v_n_ss_2=q_ss_syncleft_eval_1;
v_n_ss_3=(1+ip.Results.gamma1.*n_ss_eval_2.*(2*ip.Results.beta-ip.Results.gamma2.*n_ss_eval_2)./(ip.Results.beta-ip.Results.gamma2.*n_ss_eval_2).^2).*q_ss_2;
%size(v_n_ss)
%Cc_ss=(1+ip.Results.beta.*n_ss_eval.*(2*ip.Results.beta+ip.Results.gamma2.*n_ss_eval)./(ip.Results.beta-ip.Results.gamma2.*n_ss_eval).^2).*v_ss_2;
%V_ss_int=trapz(xmesh,v_n_ss,2);
V_ss_int=trapz([xmesh_presyn xmesh_ais xmesh_axon] ,v_n_ss_1,2)+trapz(xmesh_syncleft ,v_n_ss_2,2)+trapz(xmesh_postsyn ,v_n_ss_3,2);
%C_ss_int=trapz(xmesh,Cc_ss,2);
 S_ss_int=V_ss_int(1:nroi*nroi);
 R_ss_int=V_ss_int(nroi*nroi+1:2*nroi*nroi);
S_ss=reshape(S_ss_int,nroi,nroi);
R_ss=reshape(R_ss_int,nroi,nroi);
 %  v_ss_int=trapz(xmesh,v_ss_1,2);
% v_ss_int_1=v_ss_int(1:nroi*nroi);
% v_ss_int_2=v_ss_int(nroi*nroi+1:2*nroi*nroi);
%   V_ss_1=reshape(v_ss_int_1,nroi,nroi);
%     V_ss_2=reshape(v_ss_int_2,nroi,nroi);




    function qprime=ode_q_axon(x,W,q,A_1,B_1) 
   n=n_ss_axon(A_1,B_1,x);
 qprime=1/diff_n*-W/ip.Results.frac+1/diff_n*((1-ip.Results.frac)./ip.Results.frac).*q.*(v_a+2*v_a*ip.Results.delta.*n-(3*v_a*ip.Results.beta.*ip.Results.gamma1.*ip.Results.epsilon.*n.^2-2*v_a*ip.Results.gamma1.*ip.Results.epsilon.*ip.Results.gamma2.*n.^3)./(ip.Results.beta-ip.Results.gamma2.*n).^2-(4*v_a*ip.Results.delta.*ip.Results.beta.*ip.Results.gamma1.*ip.Results.epsilon.*n.^3-3*v_a*ip.Results.delta.*ip.Results.gamma1.*ip.Results.gamma2.*ip.Results.epsilon.*n.^4)./(ip.Results.beta-ip.Results.gamma2.*n).^2-v_r);
 %qprime=-W/(ip.Results.frac*diff_n)+(1/diff_n)*((1-ip.Results.frac)./ip.Results.frac).*q.*((v_a+2*v_a*ip.Results.delta.*n-3*v_a*ip.Results.gamma1.*ip.Results.epsilon.*n.^2./ip.Results.beta-4*v_a*ip.Results.gamma1.*ip.Results.epsilon.*ip.Results.delta.*n.^3./ip.Results.beta-v_r));
 end
 
 function nprime=ode_ss_n(x,A,n,D)
        nprime=1/D*-A;
 end

  function nprime=ode_ss_axon(x,A,n) 
 nprime=1/diff_n*-A/ip.Results.frac+1/diff_n*((1-ip.Results.frac)./ip.Results.frac).*n.*((v_a*(1+ip.Results.delta.*n).*(1-((ip.Results.gamma1*ip.Results.epsilon.*n.^2)./(ip.Results.beta-ip.Results.gamma2.*n)))-v_r));
%nprime=-A/(ip.Results.frac*diff_n)+(1/diff_n)*((1-ip.Results.frac)./ip.Results.frac).*n.*((v_a*(1+ip.Results.delta.*n).*(1-((ip.Results.gamma1*ip.Results.epsilon.*n.^2)./ip.Results.beta))-v_r));
 end


   function [maskvals] = spatial_mask(compartment)
        switch compartment
            case 'presyn'
                maskvals = (xmesh <= ip.Results.L1);
            case 'ais'
                maskvals = logical(-1 + (xmesh > ip.Results.L1) + ...
                    (xmesh < (ip.Results.L1 + ip.Results.L_ais)));
            case 'axon'
                maskvals = logical(-1 + (xmesh >= ip.Results.L1 + ip.Results.L_ais) + ...
                    (xmesh <= (ip.Results.L1 + ip.Results.L_int - ip.Results.L_syn)));
            case 'syncleft'
                maskvals = logical(-1 + (xmesh > ip.Results.L1 + ip.Results.L_int - ip.Results.L_syn) + ...
                    (xmesh < (ip.Results.L1 + ip.Results.L_int))); 
            case 'postsyn'
                maskvals = logical(-1 + (xmesh >= ip.Results.L1 + ip.Results.L_int) + ...
                    (xmesh <= (ip.Results.L1 + ip.Results.L_int + ip.Results.L2)));
%             otherwise
%                 error('Incorrect compartment specification')
        end
    end

end
