function [W_1_flux,W_2_flux,R_ss,S_ss] = MassBalanceCalculator(network_flux,tau_x0,tau_xL,matdir,varargin)
% Translated from mass_balance_network_03_09.m
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

% n_m_ss=n_ss_eval+(gamma1_new * n_ss_eval.^2)./(beta_new-ip.Results.gamma2 *n_ss_eval);
%  mass_tot_edge=trapz(xmesh,n_m_ss,2);
%  mass_tot_edge=reshape(mass_tot_edge,nroi,nroi);
%  size(n_ss_eval)

% % % 4. Steady State Calculation and Linearized Equation Definitions for v = dn/dB or v = dn/dC 
% Follows the derivation of Michiel Bertsh

% % % 4a. Presynaptic somatodendritic compartment
presyn_mask = spatial_mask('presyn');
xmesh_presyn = xmesh(presyn_mask);
x1 = xmesh_presyn(end);
n0=@(B) B;
options = odeset('RelTol',1e-4,'AbsTol',1e-4,'NonNegative',1:length(n0));
% n_ss_presyn = @(A,B,x) max((B - A.*x/diff_n),0); % B = n_ss(0) and is an integration constant
n_ss_presyn = @(A,B) ode45(@(x,n)ode_ss_n(x,A,n,diff_n),[0,L1_new],B,options);
n_ss_presyn = @(A,B,x) deval(n_ss_presyn(A,B),x);

% q_ss_presyn = @(W,V0,x) max((V0 - W.*x/diff_n),0);
q_in=@(V0) V0;
options_q = odeset('RelTol',1e-6,'AbsTol',1e-6,'NonNegative',1:length(q_in));
q_ss_presyn = @(W,V0) ode45(@(x,q)ode_ss_n(x,W,q,diff_n),[0,L1_new],q_in(V0),options_q);
q_ss_presyn = @(W,V0,x) deval(q_ss_presyn(W,V0),x);


% % % 4b. Axon initial segment
ais_mask = spatial_mask('ais');
xmesh_ais = xmesh(ais_mask);
x2 = xmesh_ais(end);
% n_ss_ais = @(A,B,x)max((B - A.*x1/diff_n - A.*(x-x1)/(diff_n*ip.Results.lambda1)),0);
n_ss_ais = @(A,B) ode45(@(x,n)ode_ss_n(x,A,n,diff_n*ip.Results.lambda1),[L1_new,L1_new+L_ais_new],n_ss_presyn(A,B,x1),options);
n_ss_ais = @(A,B,x) deval(n_ss_ais(A,B),x);

% q_ss_ais = @(W,V0,x)max((V0 - W.*x1/diff_n - W.*(x-x1)/(diff_n*ip.Results.lambda1)),0);
q_ss_ais = @(W,V0) ode45(@(x,q)ode_ss_n(x,W,q,(diff_n*ip.Results.lambda1)),[L1_new,L1_new+L_ais_new],q_ss_presyn(W,V0,x1),options_q);
q_ss_ais = @(W,V0,x) deval(q_ss_ais(W,V0),x);

% % % 4c. Axon
axon_mask = spatial_mask('axon');
xmesh_axon = xmesh(axon_mask);
x3 = xmesh_axon(end);
% nx_init = @(A_,B) max((-A_.*x1/diff_n + B - A_.*(x2 - x1)/(ip.Results.lambda1*diff_n)),0);
% options_n = odeset('RelTol',ip.Results.reltol,'AbsTol',ip.Results.abstol,'NonNegative',1:length(nx_init));
% n_ss_axon = @(A,B) ode45(@(x,n)ode_ss_axon(x,A,n),[L1_new+L_ais_new,...
%    L1_new+L_int_new-L_syn_new],nx_init(A,B),options_n); 
n_ss_axon = @(A,B) ode45(@(x,n)ode_ss_axon(x,A,n),[L1_new+L_ais_new,...
           L1_new+L_int_new-L_syn_new],n_ss_ais(A,B,x2),options);
n_ss_axon = @(A,B,x) deval(n_ss_axon(A,B),x);

% qx_init= @(W_,V0)max((-W_.*x1/diff_n + V0 - W_.*(x2 - x1)/(ip.Results.lambda1*diff_n)),0);
% options_q = odeset('RelTol',ip.Results.reltol,'AbsTol',ip.Results.abstol,'NonNegative',1:length(qx_init));
% q_ss_axon = @(W,A,B,V0) ode45(@(x,q)ode_q_axon(x,W,q,A,B),[L1_new+L_ais_new,...
%            L1_new+L_int_new-L_syn_new],qx_init(W,V0),options_q);
q_ss_axon = @(W,A,B,V0) ode45(@(x,q)ode_q_axon(x,W,q,A,B),[L1_new+L_ais_new,...
       L1_new+L_int_new-L_syn_new],q_ss_ais(W,V0,x2),options_q);
q_ss_axon = @(W,A,B,V0,x) deval(q_ss_axon(W,A,B,V0),x);

% % % 4d. Synaptic cleft
syncleft_mask = spatial_mask('syncleft');
xmesh_syncleft = xmesh(syncleft_mask);
x4 = xmesh_syncleft(end);
% n_ss_syncleft = @(A,B,x) max((n_ss_axon(A,B,x3) - A.*(x-x3)/(diff_n*ip.Results.lambda2)),0);
n_ss_syncleft = @(A,B) ode45(@(x,n)ode_ss_n(x,A,n,(diff_n*ip.Results.lambda2)),...
    [L1_new+L_int_new-L_syn_new,L1_new+L_int_new],n_ss_axon(A,B,x3),options);
n_ss_syncleft = @(A,B,x) deval(n_ss_syncleft(A,B),x);   

% q_ss_syncleft = @(W,A,B,V0,x)max((q_ss_axon(W,A,B,V0,x3) - W.*(x-x3)/(diff_n*ip.Results.lambda2)),0);
q_ss_syncleft = @(W,A,B,V0) ode45(@(x,q)ode_ss_n(x,W,q,(diff_n*ip.Results.lambda2)),...
    [L1_new+L_int_new-L_syn_new,L1_new+L_int_new],q_ss_axon(W,A,B,V0,x3),options_q);
q_ss_syncleft = @(W,A,B,V0,x) deval(q_ss_syncleft(W,A,B,V0),x); 

% % % 4e. Postsynaptic somatodendritic compartment
postsyn_mask = spatial_mask('postsyn');
xmesh_postsyn = xmesh(postsyn_mask);
x5 = xmesh_postsyn(end);
% n_ss_postsyn = @(A,B,x) max((n_ss_syncleft(A,B,x4) - A.*x/diff_n),0);
n_ss_postsyn = @(A,B,x) (n_ss_syncleft(A,B,x4) - A.*x/diff_n);
     
% q_ss_postsyn = @(W,A,B,V0,x) max((q_ss_syncleft(W,A,B,V0,x4) - W.*x/diff_n),0);
% f_q_ss = @(W,A,B,V0,V_L)((q_ss_syncleft(W,A,B,V0,x4)- W.*x5/diff_n -V_L));  
q_ss_postsyn =@(W,A,B,V0,x) (q_ss_syncleft(W,A,B,V0,x4) - W.*x/diff_n);
f_q_ss=@(W,A,B,V0,V_L)(q_ss_postsyn(W,A,B,V0,x5)-V_L); 

% % % 5. Perform shooting problem 
Adj = readmatrix([matdir filesep 'mouse_adj_matrix_19_01.csv']);
switch ip.Results.connectome_subset
    case 'Hippocampus'
        Adj = Adj([27:37 (27+213):(37+213)], [27:37 (27+213):(37+213)]);
    case 'Hippocampus+PC+RSP'
        adjinds = [27:37,78:80,147];
        adjinds = [adjinds,adjinds+213]
        Adj = Adj(adjinds,adjinds);
    case 'RH'
        Adj = Adj(1:213,1:213);
    case 'LH'
        Adj = Adj(214:end,214:end);
end
nroi = size(Adj,1);
W_1_flux = zeros(nroi);
W_2_flux = zeros(nroi);
A_0 = [0;0];
B_0 = [0;0];
V_0_0 = [1;0];
V_L_0_0 = [0;1];
g_0_0 = [0;0];
q_flux_0_0 = @(W)f_q_ss(W,A_0,B_0,V_0_0,V_L_0_0);
options = optimset('TolFun',1e-06,'Display','off');
W_1_0 = fsolve(q_flux_0_0,g_0_0,options);

for j = 1:nroi
    Ad_in = logical(Adj(:,j));
    C_1 = tau_xL(j);
    C_1 = repmat(C_1,1,length(Ad_in));
    C_1 = C_1.';
    i_app_0 = logical((tau_x0(:,j)+C_1==0).*(Ad_in));
    i_len = ones(nroi,1);
    i_app_00 = i_len(i_app_0,1);
    W_1_flux(i_app_0,j) = W_1_0(1)*i_app_00 ;  % ones(length(i_app_00),1);
    W_2_flux(i_app_0,j) = W_1_0(2)*i_app_00;  % ones(length(i_app_00),1);
    i_app = logical(((tau_x0(:,j)>0)+(C_1>0)).*(Ad_in));
    B_1 = tau_x0(i_app,j);
    
    if ~isempty(B_1)>0
        A_1 = network_flux(i_app,j);
        V0_1_1 = ones(size(B_1));
        V0_1_2 = zeros(size(B_1));
        V_L_1_1 = zeros(size(B_1));
        V_L_1_2 = ones(size(B_1));
        q_flux_1 = @(W)f_q_ss(W,A_1,B_1,V0_1_1,V_L_1_1);
        q_flux_2 = @(W)f_q_ss(W,A_1,B_1,V0_1_2,V_L_1_2);
        g0_1 = zeros(length(V0_1_1),1);
        g0_2 = zeros(length(V0_1_1),1);
        options = optimset('TolFun',1e-06,'Display','off');
        W_1_flux(i_app,j) = fsolve(q_flux_1,g0_1,options);
        W_2_flux(i_app,j) = fsolve(q_flux_2,g0_2,options);
    end
end
network_flux_vec = network_flux(:);
tau_x0_vec = tau_x0(:);
n_ss_presyn_eval= n_ss_presyn(network_flux_vec,tau_x0_vec,xmesh_presyn);
n_ss_ais_eval = n_ss_ais(network_flux_vec,tau_x0_vec,xmesh_ais);
n_ss_axon_eval = n_ss_axon(network_flux_vec,tau_x0_vec, xmesh_axon);
% n_ss_syncleft_eval = n_ss_syncleft(network_flux_vec,tau_x0_vec,xmesh_syncleft);
n_ss_postsyn_eval = n_ss_postsyn(network_flux_vec,tau_x0_vec, xmesh_postsyn);

% n_ss_eval = [n_ss_presyn_eval n_ss_ais_eval n_ss_axon_eval n_ss_syncleft_eval n_ss_postsyn_eval];
% n_ss_eval = [n_ss_eval; n_ss_eval];
n_ss_eval_1=[n_ss_presyn_eval n_ss_ais_eval n_ss_axon_eval];
n_ss_eval_1 = [n_ss_eval_1; n_ss_eval_1];
 %n_ss_eval_2=[n_ss_syncleft_eval; n_ss_syncleft_eval];
n_ss_eval_2 = [n_ss_postsyn_eval; n_ss_postsyn_eval];

W_flux = [W_1_flux(:); W_2_flux(:)];
v0_1 = [ones(length(W_flux)/2,1); zeros(length(W_flux)/2,1)];
A_flux = [network_flux_vec; network_flux_vec];
b = [tau_x0_vec; tau_x0_vec];

q_ss_presyn_eval_1 = q_ss_presyn(W_flux,v0_1,xmesh_presyn);
q_ss_ais_eval_1 = q_ss_ais(W_flux,v0_1,xmesh_ais);
q_ss_axon_eval_1 = q_ss_axon(W_flux,A_flux,b,v0_1,xmesh_axon);
q_ss_syncleft_eval_1 = q_ss_syncleft(W_flux,A_flux,b,v0_1,xmesh_syncleft);
q_ss_postsyn_eval_1 = q_ss_postsyn(W_flux,A_flux,b,v0_1,xmesh_postsyn);
% q_ss = [q_ss_presyn_eval_1 q_ss_ais_eval_1 q_ss_axon_eval_1... 
%     q_ss_syncleft_eval_1 q_ss_postsyn_eval_1];
q_ss_1=[q_ss_presyn_eval_1 q_ss_ais_eval_1 q_ss_axon_eval_1];
%q_ss_2= q_ss_syncleft_eval_1;
q_ss_2=q_ss_postsyn_eval_1;

% v_n_ss = (1+gamma1_new.*n_ss_eval.*(2*beta_new-ip.Results.gamma2.*n_ss_eval)...
%     ./(beta_new-ip.Results.gamma2.*n_ss_eval).^2).*q_ss;
%v_n_ss=(1+ip.Results.gamma1.*n_ss_eval.*(2*ip.Results.beta-ip.Results.gamma2.*n_ss_eval)./(ip.Results.beta-ip.Results.gamma2.*n_ss_eval).^2).*q_ss;
v_n_ss_1 = (1+gamma1_new.*n_ss_eval_1.*(2*beta_new-gamma2_new.*n_ss_eval_1)...
    ./(beta_new-gamma2_new.*n_ss_eval_1).^2).*q_ss_1;
v_n_ss_2 = q_ss_syncleft_eval_1;
v_n_ss_3 = (1+gamma1_new.*n_ss_eval_2.*(2*beta_new-gamma2_new.*n_ss_eval_2)...
    ./(beta_new-gamma2_new.*n_ss_eval_2).^2).*q_ss_2;
% V_ss_int = trapz(xmesh,v_n_ss,2);
V_ss_int = trapz([xmesh_presyn xmesh_ais xmesh_axon],v_n_ss_1,2)+...
    trapz(xmesh_syncleft,v_n_ss_2,2)+trapz(xmesh_postsyn,v_n_ss_3,2);
S_ss_int = V_ss_int(1:nroi*nroi);
S_ss = reshape(S_ss_int,nroi,nroi);
R_ss_int = V_ss_int(nroi*nroi+1:2*nroi*nroi);
R_ss = reshape(R_ss_int,nroi,nroi);

% % % 6. Functions

    function nprime=ode_ss_n(x,A,n,D) %#ok<INUSD> 
        nprime = 1/D*-A;
    end

    function qprime=ode_q_axon(x,W,q,A_1,B_1) 
        n = n_ss_axon(A_1,B_1,x);
%         qprime = 1/diff_n*-W/ip.Results.frac+1/diff_n*((1-ip.Results.frac)./...
%             ip.Results.frac).*q.*(v_a+2*v_a*ip.Results.delta.*n-...
%             (3*v_a*beta_new.*gamma1_new.*ip.Results.epsilon.*...
%             n.^2-2*v_a*gamma1_new.*ip.Results.epsilon.*ip.Results.gamma2.*n.^3)./...
%             (beta_new-ip.Results.gamma2.*n).^2-(4*v_a*ip.Results.delta.*...
%             beta_new.*gamma1_new.*ip.Results.epsilon.*n.^3-3*...
%             v_a*ip.Results.delta.*gamma1_new.*ip.Results.gamma2.*...
%             ip.Results.epsilon.*n.^4)./(beta_new-ip.Results.gamma2.*n).^2-v_r);
        qprime = 1/diff_n*-W/ip.Results.frac+1/diff_n*((1-ip.Results.frac)...
            ./ip.Results.frac).*q.*(v_a+2*v_a*ip.Results.delta.*n-...
            (3*v_a*beta_new.*gamma1_new.*ip.Results.epsilon.*n.^2-2*v_a*...
            gamma1_new.*ip.Results.epsilon.*gamma2_new.*n.^3)...
            ./(beta_new-gamma2_new.*n).^2-(4*v_a*ip.Results.delta.*beta_new...
            .*gamma1_new.*ip.Results.epsilon.*n.^3-3*v_a*ip.Results.delta.*...
            gamma1_new.*gamma2_new.*ip.Results.epsilon.*n.^4)...
            ./(beta_new-gamma2_new.*n).^2-v_r);
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
