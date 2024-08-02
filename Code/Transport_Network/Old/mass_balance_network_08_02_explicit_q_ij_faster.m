function [W_1_flux,W_2_flux,R_ss,S_ss]=mass_balance_network_08_02_explicit_q_ij_faster(network_flux,tau_x0,tau_xL,matdir,varargin)

if nargin < 4
    matdir = [cd filesep 'MatFiles'];
end

% % % 1. Preset values of flexible parameters
beta_ = 1e-06;
gamma1_ = 2e-05;
gamma2_ = 0;
delta_ = 1;
epsilon_ = 0.01;
lambda1_ = 0.01;%0.02  0.01 
lambda2_ = 0.01; %0.04  0.01;
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
connectome_subset_ = 'Hippocampus';
len_scale_ = 1e-3;
time_scale_ = 1;
conn_thresh_ = 'default';

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
addParameter(ip, 'conn_thresh', conn_thresh_);
parse(ip, varargin{:});
beta_new = ip.Results.beta * ip.Results.time_scale;
gamma1_new = ip.Results.gamma1 * ip.Results.time_scale;
gamma2_new = ip.Results.gamma2 * ip.Results.time_scale;
L1_new = ip.Results.L1 * ip.Results.len_scale;
L2_new = ip.Results.L2 * ip.Results.len_scale;
L_int_new = ip.Results.L_int * ip.Results.len_scale;
L_ais_new = ip.Results.L_ais * ip.Results.len_scale;
L_syn_new = ip.Results.L_syn * ip.Results.len_scale;

% % % 2. Definition of static constants

v_a = 0.7*ip.Results.len_scale * ip.Results.time_scale; % Average velocity (um/s) of anterograde active transpot (Konsack 2007)
v_r = 0.7*ip.Results.len_scale * ip.Results.time_scale; % Average velocity (um/s) of retrograde active transport (Konsack 2007)
diff_n = 12*ip.Results.len_scale^2 * ip.Results.time_scale; % Diffusivity (um^2/s) of n (Konsack 2007)

% % % 3. Definition of the (inhomogeneous) xmesh
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

%Network stuff
Adj = readmatrix([matdir filesep 'mouse_adj_matrix_19_01.csv']);
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
end
nroi = size(Adj,1);

% Follows the derivation of Michiel Bertsh
 presyn_mask = spatial_mask('presyn');
 xmesh_presyn = ip.Results.len_scale.*xmesh(presyn_mask);
  x0=xmesh_presyn(1);
 n0=@(B) B;
 options = odeset('RelTol',1e-4,'AbsTol',1e-4,'NonNegative',1:length(n0));
n_ss_presyn = @(A,B) ode45(@(x,n)ode_ss_n(x,A,n,diff_n),[0,L1_new],B,options);
n_ss_presyn = @(A,B,x) deval(n_ss_presyn(A,B),x);

% % % 5b. Steady state of axon initial segment
 ais_mask = spatial_mask('ais');
 xmesh_ais= ip.Results.len_scale.*xmesh(ais_mask);
 x1 = xmesh_presyn(end);
n_ss_ais = @(A,B) ode45(@(x,n)ode_ss_n(x,A,n,diff_n*ip.Results.lambda1),[L1_new,L1_new+L_ais_new],n_ss_presyn(A,B,x1),options);
n_ss_ais = @(A,B,x) deval(n_ss_ais(A,B),x);
 %n_ss_ais = @(A,B,x)(B - A.*x1/diff_n - A.*(x-x1)/(diff_n*ip.Results.lambda1));

% % 5c. Steady state of axon
 axon_mask = spatial_mask('axon');
 xmesh_axon= ip.Results.len_scale.*xmesh(axon_mask)
 x2 = xmesh_ais(end);

n_ss_axon = @(A,B) ode45(@(x,n)ode_ss_axon(x,A,n),[L1_new+L_ais_new,...
           L1_new+L_int_new-L_syn_new],n_ss_ais(A,B,x2),options);

 n_ss_axon= @(A,B,x) deval(n_ss_axon(A,B),x);
      

 syncleft_mask = spatial_mask('syncleft');

 xmesh_syncleft = ip.Results.len_scale.*xmesh(syncleft_mask);
 x3 = xmesh_axon(end);

n_ss_syncleft = @(A,B) ode45(@(x,n)ode_ss_n(x,A,n,(diff_n*ip.Results.lambda2)),...
    [L1_new+L_int_new-L_syn_new,L1_new+L_int_new],n_ss_axon(A,B,x3),options);
n_ss_syncleft = @(A,B,x) deval(n_ss_syncleft(A,B),x);   
 %n_ss_syncleft = @(A,B,x)(n_ss_axon(A,B,x3) - A.*(x-x3)/(diff_n*ip.Results.lambda2));
    
    

% % % 5e. Steady state of postsynaptic somatodendritic compartment
 postsyn_mask = spatial_mask('postsyn');

 xmesh_postsyn = ip.Results.len_scale.*xmesh(postsyn_mask);
 x4 = xmesh_syncleft(end);
% n_ss_postsyn = @(A,B) ode45(@(x,n)ode_ss_n(x,A,n,diff_n),ip.Results.len_scale.*[ip.Results.L1+ip.Results.L_int,L_total],n_ss_syncleft(A,B,x4),options);
%n_ss_postsyn = @(A,B,x) deval(n_ss_postsyn(A,B),x);
% n_ss_postsyn = @(A,B) ode45(@(x,n)ode_ss_n(x,A,n,diff_n),ip.Results.len_scale.*[ip.Results.L1+ip.Results.L_int, L_total],n_ss_syncleft(A,B,x4),options);
% n_ss_postsyn = @(A,B,x) deval(n_ss_postsyn(A,B),x);  
 n_ss_postsyn =@(A,B,x) (n_ss_syncleft(A,B,x4) - A.*(x-x4)/diff_n);
     

    
x5=xmesh_postsyn(end);

q_in=@(V0) V0;
 options_q = odeset('RelTol',1e-6,'AbsTol',1e-6,'NonNegative',1:length(q_in));
 %options_q = odeset('RelTol',1e-6,'AbsTol',1e-6);

 %q_ss_presyn = @(W,V0) ode45(@(x,q)ode_ss_n(x,W,q,diff_n),ip.Results.len_scale.*[0,ip.Results.L1],q_in(V0),options_q);
%q_ss_presyn = @(W,V0,x) deval(q_ss_presyn(W,V0),x);
q_ss_presyn = @(W,V0,x) (V0- W.*x/diff_n);
% % %b.Equation for v=dn/db or dn/dC in axon initial segment
% q_ss_ais = @(W,V0) ode45(@(x,q)ode_ss_n(x,W,q,(diff_n*ip.Results.lambda1)),ip.Results.len_scale.*[ip.Results.L1,ip.Results.L1+ip.Results.L_ais],q_ss_presyn(W,V0,x1),options_q);
%q_ss_ais = @(W,V0,x) deval(q_ss_ais(W,V0),x);

q_ss_ais = @(W,V0,x)(V0 - W.*x1/diff_n - W.*(x-x1)/(diff_n*ip.Results.lambda1));
%c. Equation for v=dn/db or dn/dC in axon 
% qx_init= @(W_,V0)(-W_.*x1/diff_n + V0 - W_.*(x2 - x1)/(ip.Results.lambda1*diff_n));
%nx_init=@(A_,B)  n_ss_ais(A_,B,x2);
 %options = odeset('RelTol',1e-4,'AbsTol',1e-4,'NonNegative',1:length(qx_init));

 q_ss_axon = @(W,A,B,V0) ode45(@(x,q)ode_q_axon(x,W,q,A,B),ip.Results.len_scale.*[ip.Results.L1+ip.Results.L_ais,...
           ip.Results.L1+ip.Results.L_int-ip.Results.L_syn],q_ss_ais(W,V0,x2),options_q);

 q_ss_axon= @(W,A,B,V0,x) deval(q_ss_axon(W,A,B,V0),x);
 q_ss_syncleft_10=@(S,x)S.*((x4-x3)/(diff_n*ip.Results.lambda2)+(x5-x4)/(diff_n)-(x-x3)./(diff_n*ip.Results.lambda2));
q_ss_syncleft_01=@(S,x)1-S.*((x4-x3)/(diff_n*ip.Results.lambda2)+(x5-x4)/(diff_n)-(x-x3)./(diff_n*ip.Results.lambda2));
q_ss_postsyn_10=@(S,x)S.*((x5-x4)/(diff_n)-(x-x4)./(diff_n));
q_ss_postsyn_01=@(S,x)1-S.*((x5-x4)/(diff_n)-(x-x4)./(diff_n));



W_1_flux=zeros(nroi);

tic
fprintf('W_1 calculation \n')
A_0=[0];
B_0=[0];
[W_1_0]=Q_calculation(A_0,B_0);

 for j=1:nroi
     
     Ad_in=logical(Adj(:,j));
      C_1=tau_xL(j);
    C_1=repmat(C_1,1,length(Ad_in));
   C_1=C_1.';
   i_app_0=logical((tau_x0(:,j)+C_1==0 ).*(Ad_in));
   i_len=ones(nroi,1);
   i_app_00=i_len(i_app_0,1);
    W_1_flux(i_app_0,j)=W_1_0*i_app_00 ; 

     i_app=logical(((tau_x0(:,j)>0)+(C_1>0)).*(Ad_in));
    B_1=tau_x0(i_app,j);
   if ~isempty(B_1)>0
      A_1=network_flux(i_app,j);
    
[W_1_flux(i_app,j)]=Q_calculation(A_1,B_1);

   end
 end
 toc

fprintf('n_ij calculation \n')
tic
network_flux=network_flux(:);
tau_x0=tau_x0(:);
 
      
        
 n_ss_presyn_eval= n_ss_presyn(network_flux, tau_x0,xmesh_presyn);

 n_ss_ais_eval=n_ss_ais(network_flux,tau_x0,xmesh_ais);

  
 n_ss_axon_eval=n_ss_axon(network_flux,tau_x0, xmesh_axon);

 %n_ss_syncleft_eval=n_ss_syncleft(network_flux,tau_x0,xmesh_syncleft);


 n_ss_postsyn_eval=n_ss_postsyn(network_flux,tau_x0, xmesh_postsyn);
toc

fprintf('W_2 calculation \n')
tic
 % W_2 calculation
 F_ax=fun_F_n_axon(xmesh_axon,n_ss_axon_eval);
 int_num_Q01=trapz(xmesh_axon,F_ax,2);
   Q_num01=exp(int_num_Q01);
   Q_num01=reshape(Q_num01,nroi,nroi);
   W_2_flux=-Q_num01.*W_1_flux;
toc
 S_num10=exp(-int_num_Q01);
   S_num10=reshape(S_num10,nroi,nroi);
 S_10=S_num10.*W_1_flux;
 S_10=S_10(:);
 S_01=1.*W_1_flux;
 S_01=S_01(:);

tic
 

 n_ss_eval_1=[n_ss_presyn_eval n_ss_ais_eval n_ss_axon_eval];

 
 n_ss_eval_1=[n_ss_eval_1; n_ss_eval_1];

 n_ss_eval_2=[n_ss_postsyn_eval; n_ss_postsyn_eval];
 

fprintf('q_ij calculation \n')

W_flux=[W_1_flux(:);W_2_flux(:)];

v0_1=[ones(length(W_flux)/2,1);zeros(length(W_flux)/2,1)];

A_flux=[network_flux(:);network_flux(:)];
b=[tau_x0(:);tau_x0(:)];
tic
 q_ss_presyn_eval_1= q_ss_presyn(W_flux, v0_1,xmesh_presyn);

 q_ss_ais_eval_1=q_ss_ais(W_flux,v0_1,xmesh_ais);

 
 q_ss_axon_eval_1=q_ss_axon(W_flux,A_flux,b,v0_1, xmesh_axon);
 

 %q_ss_syncleft_eval_1=q_ss_syncleft(W_flux,A_flux,b,v0_1,xmesh_syncleft);
q_ss_syncleft_eval_10=q_ss_syncleft_10(S_10,xmesh_syncleft);
q_ss_syncleft_eval_01=q_ss_syncleft_01(S_01,xmesh_syncleft);
q_ss_syncleft_eval_1=[q_ss_syncleft_eval_10;q_ss_syncleft_eval_01 ];

 %q_ss_postsyn_eval_1=q_ss_postsyn(W_flux,A_flux,b,v0_1, xmesh_postsyn);
 q_ss_postsyn_eval_10=q_ss_postsyn_10(S_10,xmesh_postsyn);
q_ss_postsyn_eval_01=q_ss_postsyn_01(S_01,xmesh_postsyn);
q_ss_postsyn_eval_1=[q_ss_postsyn_eval_10;q_ss_postsyn_eval_01 ];
toc
 
 %q_ss=[q_ss_presyn_eval_1 q_ss_ais_eval_1 q_ss_axon_eval_1 q_ss_syncleft_eval_1 q_ss_postsyn_eval_1];
 q_ss_1=[q_ss_presyn_eval_1 q_ss_ais_eval_1 q_ss_axon_eval_1];
 %q_ss_2= q_ss_syncleft_eval_1;
 q_ss_2=q_ss_postsyn_eval_1;
 
 
fprintf('C_ij calculation \n')

tic
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
toc
 %  v_ss_int=trapz(xmesh,v_ss_1,2);
% v_ss_int_1=v_ss_int(1:nroi*nroi);
% v_ss_int_2=v_ss_int(nroi*nroi+1:2*nroi*nroi);
%   V_ss_1=reshape(v_ss_int_1,nroi,nroi);
%     V_ss_2=reshape(v_ss_int_2,nroi,nroi);
    function [Q_10]=Q_calculation(A,B)
             % n=@(x)n_ss_axon(A,B,x);
             
            n_mesh=n_ss_axon(A,B,xmesh_axon);
           F=fun_F_n_axon(xmesh_axon,n_mesh);
  
        int_num=trapz(xmesh_axon,F,2);
        Q_num=exp(int_num);
        
  %int=xmesh_axon;
  Q_den_sd=1/diff_n*(x1-x0);
  Q_den_ais=1/(ip.Results.lambda1*diff_n)*(x2-x1);
  int_den_1=zeros(length(B),length(xmesh_axon));
  %tic
  for i=2:length(xmesh_axon)
       n_mesh_1=n_ss_axon(A,B,xmesh_axon(1:i));
      F1=fun_F_n_axon(xmesh_axon(1:i),n_mesh_1);
      int_den_1(:,i)=trapz(xmesh_axon(1:i),F1, 2);

  end
  %toc
   Q_den_axon=trapz(xmesh_axon,exp(int_den_1),2);
  Q_den_axon=1/(ip.Results.frac.*diff_n).*Q_den_axon;
  Q_den_syn_cleft=1/(ip.Results.lambda2*diff_n)*(x4-x3)*Q_num;
  Q_den_sd_postsyn=1/(diff_n)*(x5-x4)*Q_num;
 Q_den=Q_den_sd+Q_den_ais+Q_den_axon+Q_den_syn_cleft+Q_den_sd_postsyn;
% Q_01=-(Q_num./Q_den); 
 Q_10=(1./Q_den);
% Q=[Q_10;Q_01]


end

  function F=fun_F_n_axon(x,n)
     %F=-(1/diff_n)*((1-ip.Results.frac)./ip.Results.frac).*((v_a+2*v_a.*ip.Results.delta.*n-3*v_a.*ip.Results.gamma1.*ip.Results.epsilon.*n.^2./ip.Results.beta-4*v_a.*ip.Results.gamma1.*ip.Results.epsilon.*ip.Results.delta.*n.^3./ip.Results.beta-v_r));
   F=-1/diff_n*((1-ip.Results.frac)./ip.Results.frac).*(v_a+2*v_a*ip.Results.delta.*n-(3*v_a*ip.Results.beta.*ip.Results.gamma1.*ip.Results.epsilon.*n.^2-2*v_a*ip.Results.gamma1.*ip.Results.epsilon.*ip.Results.gamma2.*n.^3)./(ip.Results.beta-ip.Results.gamma2.*n).^2-(4*v_a*ip.Results.delta.*ip.Results.beta.*ip.Results.gamma1.*ip.Results.epsilon.*n.^3-3*v_a*ip.Results.delta.*ip.Results.gamma1.*ip.Results.gamma2.*ip.Results.epsilon.*n.^4)./(ip.Results.beta-ip.Results.gamma2.*n).^2-v_r);
  %this is F_n/a(x)
  end
  



   function qprime=ode_q_axon(x,W,q,A_1,B_1)
   
   n=n_ss_axon(A_1,B_1,x);
  
   F=fun_F_n_axon(x,n);
   
 %qprime=1/diff_n*-W/ip.Results.frac+1/diff_n*((1-ip.Results.frac)./ip.Results.frac).*q.*(v_a+2*v_a*ip.Results.delta.*n-(3*v_a*ip.Results.beta.*ip.Results.gamma1.*ip.Results.epsilon.*n.^2-2*v_a*ip.Results.gamma1.*ip.Results.epsilon.*ip.Results.gamma2.*n.^3)./(ip.Results.beta-ip.Results.gamma2.*n).^2-(4*v_a*ip.Results.delta.*ip.Results.beta.*ip.Results.gamma1.*ip.Results.epsilon.*n.^3-3*v_a*ip.Results.delta.*ip.Results.gamma1.*ip.Results.gamma2.*ip.Results.epsilon.*n.^4)./(ip.Results.beta-ip.Results.gamma2.*n).^2-v_r);
 
 qprime=1/diff_n*-W/ip.Results.frac-F.*q;
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
