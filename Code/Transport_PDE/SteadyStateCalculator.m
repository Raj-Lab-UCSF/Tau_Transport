function [n_ss, m_ss, xmesh, B_ss] = SteadyStateCalculator(varargin)
% % %
% Calculates the steady state for the differential equations 

% % % 1. Set default values and initialize inputParser

alpha_ = 0; % Not currently explored in the model
beta_ = 1e-06;
gamma_ = 2e-05;
delta_ = 1;
epsilon_ = 0.01;
lambda1_ = 0.01;
lambda2_ = 0.01;
frac_ = 0.92; % Average fraction of n diffusing (Konsack 2007)
L_int_ = 1000; % in micrometers
L1_ = 200;
L2_ = 200; 
L_ais_ = 40;
L_syn_ = 40;
n0_ = zeros(1, L_int_);
m0_ = zeros(1, L_int_);
N1_0_ = 0; 
N2_0_ = 0; 
M1_0_ = 0; 
M2_0_ = 0;
resmesh_ = 'fine';
total_mass_ = 184;

ip = inputParser;
validScalar = @(x) isnumeric(x) && isscalar(x) && (x>=0);
validNonnegative = @(x) isnumeric(x) && (sum(x>=0) == length(x));
addParameter(ip, 'alpha', alpha_, validScalar);
addParameter(ip, 'beta', beta_, validScalar);
addParameter(ip, 'gamma', gamma_, validScalar);
addParameter(ip, 'delta', delta_, validScalar);
addParameter(ip, 'epsilon', epsilon_, validScalar);
addParameter(ip, 'frac', frac_, validScalar);
addParameter(ip, 'lambda1', lambda1_, validScalar);
addParameter(ip, 'lambda2', lambda2_, validScalar);
addParameter(ip, 'L_int', L_int_, validScalar);
addParameter(ip, 'L1', L1_, validScalar);
addParameter(ip, 'L2', L2_, validScalar);
addParameter(ip, 'n0', n0_, validNonnegative);
addParameter(ip, 'm0', m0_, validNonnegative);
addParameter(ip, 'N1_0', N1_0_, validScalar);
addParameter(ip, 'N2_0', N2_0_, validScalar);
addParameter(ip, 'M1_0', M1_0_, validScalar);
addParameter(ip, 'M2_0', M2_0_, validScalar);
addParameter(ip, 'resmesh', resmesh_);
addParameter(ip, 'L_ais', L_ais_);
addParameter(ip, 'L_syn', L_syn_);
addParameter(ip, 'total_mass',total_mass_)
parse(ip, varargin{:});

% The input parser creates a struct, ip, that has fields initialized with
% the default values indicated above (with an trailing underscore), which
% can be overwritten by passing a (parameter name (string), parameter
% value) pair to the function, where the allowable parameter names are
% given by the second arguments to the addParameter function. Allowable
% parameter values are constrained by validScalar and validNonnegative and
% the function will throw an error if these constraints aren't met.

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

% % % 4. Definition of ICs
n_init0 = [ip.Results.N1_0*ones(1,ip.Results.L1),ip.Results.n0,ip.Results.N2_0*ones(1,ip.Results.L2)];
m_init0 = [ip.Results.M1_0*ones(1,ip.Results.L1),ip.Results.m0,ip.Results.M2_0*ones(1,ip.Results.L2)];
rescale_factor = ip.Results.total_mass/(trapz(n_init0) + trapz(m_init0));
n_init0 = interp1(1:L_total,n_init0,xmesh,'makima'); n_init0(n_init0 < 0) = 0;
m_init0 = interp1(1:L_total,m_init0,xmesh,'makima'); m_init0(m_init0 < 0) = 0;
% n_init = n_init0 * rescale_factor; 
m_init = m_init0 * rescale_factor;

% % % 5. Definition of Steady-State Function
A = 0; % Specified to be 0 for Neumann BC
B0 = 0.04; % Initialization of B
% % % 5a. Steady State of Presynaptic Somatodendritic Compartment
% Follows the derivation of Michiel Bertsh
presyn_mask = spatial_mask('presyn');
% n0_presyn = n_init(presyn_mask); 
% m0_presyn = m_init(presyn_mask);
xmesh_presyn = xmesh(presyn_mask);
n_ss_presyn = @(B,x) (B - A*x/diff_n); % B = n_ss(0) and is an integration constant
m_ss_presyn = @(B,x) (ip.Results.gamma*n_ss_presyn(B,x)^2) / ...
                (ip.Results.beta - ip.Results.gamma*n_ss_presyn(B,x));

% % % 5b. Steady state of axon initial segment
ais_mask = spatial_mask('ais');
% n0_ais = n_init(ais_mask); 
% m0_ais = m_init(ais_mask);
xmesh_ais= xmesh(ais_mask);
x1 = xmesh_presyn(end);
n_ss_ais = @(B,x) (B - A*x1/diff_n - A*(x-x1)/(diff_n*ip.Results.lambda1));
m_ss_ais = @(B,x) (ip.Results.gamma*n_ss_ais(B,x)^2) / ...
                (ip.Results.beta - ip.Results.gamma*n_ss_ais(B,x));

% % % 5c. Steady state of axon
axon_mask = spatial_mask('axon');
% n0_axon = n_init(axon_mask); 
% m0_axon = m_init(axon_mask);
xmesh_axon= xmesh(axon_mask);
x2 = xmesh_ais(end);
nx_init = @(B_) (-A*x1 + diff_n*B_ - A*(x2 - x1)/ip.Results.lambda1)/diff_n;
% nx_init = nx2(B0);
% n_ss_axon_fun = @(nx,B,x) (x - x2)/diff_n - integral(@(nx_)(-A/ip.Results.frac +...
%     ((1 - ip.Results.frac).*nx_./ip.Results.frac) .* (v_a.*(1 + ip.Results.delta.*nx_).*...
%     (1 - (ip.Results.gamma.*ip.Results.epsilon.*nx_.^2)./(ip.Results.beta - ...
%     ip.Results.gamma.*nx_)) - v_r)).^(-1),nx2(B),nx);
% n_ss_axon = @(B,x) fsolve(@(nx)n_ss_axon_fun(nx,B,x),nx_init);
options = odeset('RelTol',1e-4,'AbsTol',1e-4);
odefun = @(x,nx_,B) (-A/ip.Results.frac +...
    ((1 - ip.Results.frac).*nx_./ip.Results.frac) .* (v_a.*(1 + ip.Results.delta.*nx_).*...
    (1 - (ip.Results.gamma.*ip.Results.epsilon.*nx_.^2)./(ip.Results.beta - ...
    ip.Results.gamma.*nx_)) - v_r))./diff_n;
n_ss_axon = @(B) ode45(@(x,nx)odefun(x,nx,B),[ip.Results.L1+ip.Results.L_ais,...
    ip.Results.L1+ip.Results.L_int-ip.Results.L_syn],nx_init(B),options);
n_ss_axon = @(B,x) deval(n_ss_axon(B),x);
m_ss_axon = @(B,x) (ip.Results.gamma*n_ss_axon(B,x)^2) / ...
                (ip.Results.beta - ip.Results.gamma*n_ss_axon(B,x));

% % % 5d. Steady state of synaptic cleft
syncleft_mask = spatial_mask('syncleft');
% n0_syncleft = n_init(syncleft_mask); 
% m0_syncleft = m_init(syncleft_mask);
xmesh_syncleft = xmesh(syncleft_mask);
x3 = xmesh_axon(end);
n_ss_syncleft = @(B,x) (n_ss_axon(B,x3) - A*(x-x3)/(diff_n*ip.Results.lambda2));
m_ss_syncleft = @(B,x) m_init(xmesh == x);
% m_ss_syncleft = @(B,x) 0;

% % % 5e. Steady state of postsynaptic somatodendritic compartment
postsyn_mask = spatial_mask('postsyn');
% n0_postsyn = n_init(postsyn_mask); 
% m0_postsyn = m_init(postsyn_mask);
xmesh_postsyn = xmesh(postsyn_mask);
x4 = xmesh_syncleft(end);
n_ss_postsyn = @(B,x) (n_ss_syncleft(B,x4) - A*x/diff_n); 
m_ss_postsyn = @(B,x) (ip.Results.gamma*n_ss_postsyn(B,x)^2) / ...
                (ip.Results.beta - ip.Results.gamma*n_ss_postsyn(B,x));

% % % 5f. Combining into a B-dependent vector over xmesh
n_ss_funs = cell(1,length(xmesh));
m_ss_funs = n_ss_funs;
for i = 1:length(xmesh)
    xi = xmesh(i);
    if ismember(xi,xmesh_presyn)
        n_ss_funs{i} = @(B) n_ss_presyn(B,xi);
        m_ss_funs{i} = @(B) m_ss_presyn(B,xi);
    elseif ismember(xi,xmesh_ais)
        n_ss_funs{i} = @(B) n_ss_ais(B,xi);
        m_ss_funs{i} = @(B) m_ss_ais(B,xi);
    elseif ismember(xi,xmesh_axon)
        n_ss_funs{i} = @(B) n_ss_axon(B,xi);
        m_ss_funs{i} = @(B) m_ss_axon(B,xi);
    elseif ismember(xi,xmesh_syncleft)
        n_ss_funs{i} = @(B) n_ss_syncleft(B,xi);
        m_ss_funs{i} = @(B) m_ss_syncleft(B,xi); 
    elseif ismember(xi,xmesh_postsyn)
        n_ss_funs{i} = @(B) n_ss_postsyn(B,xi);
        m_ss_funs{i} = @(B) m_ss_postsyn(B,xi);
    end
end

% % % % 6. Solving in terms of B and returning steady-state values
total_ss_funs = cell(1,length(n_ss_funs));
for i = 1:length(total_ss_funs)
    total_ss_funs{i} = @(B) n_ss_funs{i}(B) + m_ss_funs{i}(B);
end

tic;
% Testing with fmincon
numvals = 5;
int_total_ss = @(B) cell_integral(B,total_ss_funs,xmesh);
int_fun = @(B) 100*abs(ip.Results.total_mass - int_total_ss(B))/ip.Results.total_mass;
B_ss_test = linspace(0.5*ip.Results.beta/ip.Results.gamma,0.95*ip.Results.beta/ip.Results.gamma,numvals);
int_ss_test = zeros(1,length(B_ss_test));
for i = 1:length(B_ss_test)
    fprintf('Initialization iter %d/%d \n',i,length(B_ss_test))
    int_ss_test(i) = int_fun(B_ss_test(i));
end
[~,intind] = min(int_ss_test);
B_ss0 = B_ss_test(intind);
options = optimoptions('fmincon','Display','iter','OptimalityTolerance',1e-4,...
    'StepTolerance',1e-4,'ConstraintTolerance',1e-4);
B_ss = fmincon(int_fun, B_ss0,[],[],[],[],0,ip.Results.beta/ip.Results.gamma,[],options);

% Grid search
% numvals = 50;
% B_ss_test = linspace(0.5*ip.Results.beta/ip.Results.gamma,0.95*ip.Results.beta/ip.Results.gamma,numvals);
% int_ss_test = zeros(1,length(B_ss_test));
% for i = 1:length(B_ss_test)
%     fprintf('Iter %d/%d \n',i,length(B_ss_test))
%     int_ss_test(i) = int_fun(B_ss_test(i));
% end
% [intmin,intind] = min(int_ss_test);
% B_ss = B_ss_test(intind);
toc;
n_ss = zeros(1,length(n_ss_funs)); m_ss = n_ss;
for i = 1:length(n_ss)
    n_ss(i) = n_ss_funs{i}(B_ss);
    m_ss(i) = m_ss_funs{i}(B_ss);
end
fprintf('Percent Error = %0.2f \n',int_fun(B_ss))

% % % 7. Other functions
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

    function [int_trap] = cell_integral(B_,fun_cell,xmesh_)
        dxs = xmesh_(2:end) - xmesh_(1:(end-1));
        funcalls = zeros(1,length(xmesh_));
        for j = 1:length(funcalls)
            fun_cellj = @(b)fun_cell{j}(b);
            funcalls(j) = fun_cellj(B_);
        end
        funcalls_avg = (funcalls(1:(end-1)) + funcalls(2:end))/2;
        int_trap = sum(funcalls_avg .* dxs);
    end
end