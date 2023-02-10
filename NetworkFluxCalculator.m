function network_flux = NetworkFluxCalculator(tau_x0,tau_xL,matdir,varargin)
if nargin < 3
    matdir = [cd filesep 'FromVeronica'];
end

beta_ = 1e-06;
gamma1_ = 2e-05;
gamma2_ = 2e-06;
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
parse(ip, varargin{:});

% % % 2. Definition of constants
v_a = 0.7; % Average velocity (um/s) of anterograde active transpot (Konsack 2007)
v_r = 0.7; % Average velocity (um/s) of retrograde active transport (Konsack 2007)
diff_n = 12; % Diffusivity (um/s^2) of n (Konsack 2007)
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

% % % 5. Steady State Calculation
% Follows the derivation of Michiel Bertsh

% % % 5a. Presynaptic somatodendritic compartment
presyn_mask = spatial_mask('presyn');
xmesh_presyn = xmesh(presyn_mask);
%n_ss_presyn = @(A,B,x) max((B - A.*x/diff_n),0); 
%m_ss_presyn = @(A,B,x) ((ip.Results.gamma*n_ss_presyn(A,B,x).^2)./(ip.Results.beta - n_ss_presyn(A,B,x)*ip.Results.gamma));
%m_ss_presyn = @(A,x) ((ip.Results.gamma*n_ss_presyn(A,x).^2)./(ip.Results.beta));

% % % 5b. Axon initial segment
ais_mask = spatial_mask('ais');
xmesh_ais= xmesh(ais_mask);
x1 = xmesh_presyn(end);
%n_ss_ais = @(A,B,x) max((B - A.*x1/diff_n - A.*(x-x1)/(diff_n*ip.Results.lambda1)),0);
%n_ss_ais=@(A,B,x)max((n_ss_presyn(A,B,x1)-A*(x-x1)/(diff_n*ip.Results.lambda1)),0);
%m_ss_ais = @(A,B,x) (ip.Results.gamma*n_ss_ais(A,B,x).^2)./(ip.Results.beta - n_ss_ais(A,B,x)*ip.Results.gamma);
%m_ss_ais = @(A,x) (ip.Results.gamma*n_ss_ais(A,x).^2)./(ip.Results.beta);

% % % 5c. Axon
axon_mask = spatial_mask('axon');
xmesh_axon= xmesh(axon_mask);
x2 = xmesh_ais(end);
nx_init = @(A_,B) max((-A_.*x1/diff_n + B - A_.*(x2 - x1)/(ip.Results.lambda1*diff_n)),0);
%nx_init=@(A_,B)  n_ss_ais(A_,B,x2);
options = odeset('RelTol',ip.Results.reltol,'AbsTol',ip.Results.abstol,'NonNegative',1:length(nx_init));
n_ss_axon = @(A,B) ode45(@(x,n)ode_ss_axon(x,A,n),[ip.Results.L1+ip.Results.L_ais,...
   ip.Results.L1+ip.Results.L_int-ip.Results.L_syn],nx_init(A,B),options);
n_ss_axon = @(A,B,x) deval(n_ss_axon(A,B),x);
%       m_ss_axon = @(A,B,x) (ip.Results.gamma*n_ss_axon(A,B,x).^2)./ ...
%                   (ip.Results.beta - ip.Results.gamma*n_ss_axon(A,B,x));
%m_ss_axon = @(A,x) (ip.Results.gamma*n_ss_axon(A,x).^2)./ ...
       % (ip.Results.beta);

% % % 5d. Synaptic cleft
syncleft_mask = spatial_mask('syncleft');
xmesh_syncleft = xmesh(syncleft_mask);
x3 = xmesh_axon(end);
n_ss_syncleft = @(A,B,x) max((n_ss_axon(A,B,x3) - A.*(x-x3)/(diff_n*ip.Results.lambda2)),0);
%m_ss_syncleft = M1_0_*ones(1,length(xmesh_syncleft));

% % % 5e. Postsynaptic somatodendritic compartment
postsyn_mask = spatial_mask('postsyn');
xmesh_postsyn = xmesh(postsyn_mask);
x4 = xmesh_syncleft(end);
%       n_ss_postsyn = @(A,B,x) max((n_ss_syncleft(A,B,x4) - A.*x/diff_n),0); 
% m_ss_postsyn = @(A,B,x) (ip.Results.gamma*n_ss_postsyn(A,B,x).^2)./ ...
       % (ip.Results.beta - ip.Results.gamma*n_ss_postsyn(A,B,x));
% m_ss_postsyn = @(A,x) (ip.Results.gamma*n_ss_postsyn(A,x).^2)./ ...
%          (ip.Results.beta);
x5 = xmesh_postsyn(end);
%f_ss=@(A)(n_ss_axon(A,x3)- (A/(diff_n*ip.Results.lambda2))*((1-ip.Results.lambda2)*x4+ip.Results.lambda2*x5-x3)-C_1);
% f_ss=@(A,B,C)(n_ss_postsyn(A,B,x5)-C);
f_ss = @(A,B,C)(n_ss_syncleft(A,B,x4) - A .* x5/diff_n-C);

% % % 6. Flux calculation on network 
Adj = readmatrix([matdir filesep 'mouse_adj_matrix_19_01.csv']);
nroi = size(Adj,1);
network_flux = zeros(nroi);
for i = 1:nroi
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

% Functions
    function nprime = ode_ss_axon(x,A,n) %#ok<INUSD> 
        nprime = (1/diff_n * -A ./ ip.Results.frac) + ...
            1/diff_n * ((1 - ip.Results.frac) ./ ip.Results.frac) .* n...
            .* (v_a*(1 + ip.Results.delta .* n) .* (1 -...
            ((ip.Results.gamma1 * ip.Results.epsilon .* n .^2) ./ ...
            (ip.Results.beta - ip.Results.gamma2 .* n))) - v_r);
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

% n_ss_postsyn=n_ss_postsyn(A_1,xmesh_postsyn);
% m_ss_postsyn=m_ss_postsyn(A_1,xmesh_postsyn);
% n_ss_postsyn(end);
% figure
% subplot(2,1,1)
% plot([xmesh_presyn xmesh_ais xmesh_axon xmesh_syncleft xmesh_postsyn],[ n_ss_presyn n_ss_ais n_ss_axon n_ss_syncleft n_ss_postsyn]);
% 
% xlabel('x');
% ylabel('n(x)');
% 
% subplot(2,1,2);
% plot([xmesh_presyn xmesh_ais xmesh_axon xmesh_syncleft xmesh_postsyn], [m_ss_presyn m_ss_ais m_ss_axon m_ss_syncleft  m_ss_postsyn]);
% 
% xlabel('x');
% ylabel('m(x)');