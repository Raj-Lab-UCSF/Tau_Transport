function network_flux = NetworkFluxCalculator(tau_x0,tau_xL,matdir,varargin)

if nargin < 3
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
beta_new = ip.Results.beta*ip.Results.time_scale;
gamma_new = ip.Results.gamma1*ip.Results.time_scale;
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

% % % 4. Steady State Calculation
% Follows the derivation of Michiel Bertsh

% % % 4a. Presynaptic somatodendritic compartment
presyn_mask = spatial_mask('presyn');
xmesh_presyn = xmesh(presyn_mask);
x1 = xmesh_presyn(end);

% % % 4b. Axon initial segment
ais_mask = spatial_mask('ais');
xmesh_ais = xmesh(ais_mask);
x2 = xmesh_ais(end);

% % % 4c. Axon
axon_mask = spatial_mask('axon');
xmesh_axon = xmesh(axon_mask);
x3 = xmesh_axon(end);
nx_init = @(A_,B) max((-A_.*x1/diff_n + B - A_.*(x2 - x1)/(ip.Results.lambda1*diff_n)),0);
options = odeset('RelTol',ip.Results.reltol,'AbsTol',ip.Results.abstol,'NonNegative',1:length(nx_init));
n_ss_axon = @(A,B) ode45(@(x,n)ode_ss_axon(x,A,n),[L1_new+L_ais_new,...
   L1_new+L_int_new-L_syn_new],nx_init(A,B),options);
n_ss_axon = @(A,B,x) deval(n_ss_axon(A,B),x);

% % % 4d. Synaptic cleft
syncleft_mask = spatial_mask('syncleft');
xmesh_syncleft = xmesh(syncleft_mask);
x4 = xmesh_syncleft(end);
n_ss_syncleft = @(A,B,x) max((n_ss_axon(A,B,x3) - A.*(x-x3)/(diff_n*ip.Results.lambda2)),0);

% % % 4e. Postsynaptic somatodendritic compartment
postsyn_mask = spatial_mask('postsyn');
xmesh_postsyn = xmesh(postsyn_mask);
x5 = xmesh_postsyn(end);
f_ss = @(A,B,C)(n_ss_syncleft(A,B,x4) - A .* x5/diff_n-C);

% % % 5. Flux calculation on network 
Adj = readmatrix([matdir filesep 'mouse_adj_matrix_19_01.csv']);
switch ip.Results.connectome_subset
    case 'Hippocampus'
        Adj = Adj([27:37 (27+213):(37+213)], [27:37 (27+213):(37+213)]);
    case 'RH'
        Adj = Adj(1:213,1:213);
    case 'LH'
        Adj = Adj(214:end,214:end);
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

% % % 6. Functions
    function nprime=ode_ss_axon(x,A,n) %#ok<INUSD> 
        nprime = -A/(ip.Results.frac*diff_n)+(1/diff_n)*((1-ip.Results.frac)./...
        ip.Results.frac).*n.*((v_a*(1+ip.Results.delta.*n).*(1-((gamma_new...
        *ip.Results.epsilon.*n.^2)./beta_new))-v_r));
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
