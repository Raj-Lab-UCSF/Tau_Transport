function [n, m, xmesh, trange, jn, jm] = TauTransportPDE(varargin)
% % % 0. Introduction
% This is a function that creates a phenomenological model of pathological
% tau axonal transport and solves the following system of coupled 1D PDES:
% 
% Dn/Dt = D[frac*diff_n*Dn/Dx+v*(1-frac)*n]/Dx + alpha + beta*m - gamma*n*(n+m)
% Dm/Dt = diff_m*D^2m/Dx^2 - beta*m + gamma*n*(n+m)
% v(n,m;delta,epsilon) = v_a*(1 + delta*n)*(1 - epsilon*m) - v_r
% 
% The key aspect of this model will be to encode a SPATIAL MASK that takes  
% into account the different properties of WM and GM regions. Namely, while
% the equations above hold in white matter (WM) regions, there is no active
% transport expected in gray matter (GM) regions (v = 0). Additionally, the
% rate of diffusion between compartments is expected to be much slower than
% that within compartments due to the necessity of either trans-synaptic
% or axonal initial segment crossing. Therefore, the diffusion constant of
% soluble and insoluble tau (diff_n, diff_m) will be piecewise linear
% functions of x, as will v_a and v_r. The boundary conditions (bordering
% the two GM regions) can then be simply encoded as simple Neumann perfectly
% reflecting boundaries.
%
% Primary parameters (inputs)
% alpha: generative (linear) term of phosphorylated monomers (uM/s)
% beta: fragmentation rate of aggregates to form mobile 'monomers' (1/s)
% gamma: aggregation rate (1/(uM*s)), ~ 2 (Northup 1992)
% delta (contained in v): rate of anterograde transport enhancement by monomers
% (1/uM)
% epsilon (contained in v): rate of anterograde transport inhibition by aggs
% (1/uM)
% lambda: scaling constant of the rate of diffusion across GM-WM boundaries
%
% Secondary parameters (inputs)
% N1: Initial concentration of soluble tau in GM compartment 1
% N2: Initial concentration of soluble tau in GM compartment 2
% M1: Initial concentration of insoluble tau in GM compartment 1
% M2: Initial concentration of insoluble tau in GM compartment 2
% n0: Initial concentration profile of soluble tau in WM compartment
% (should be a 1D array of length L_int)
% m0: Initial concentration profile of insoluble tau in WM compartment
% (should be a 1D array of length (L_int)
% L1: Length of GM compartment 1 in um 
% L2: Length of GM compartment 2 in um 
% L_int: Length of the WM compartment in um
% T: Time of simulation in seconds
% tsteps: Number of time points to plot the PDEs over 
% (num_com: Number of compartments to evaluate the PDES over - should be
% internal to the function)
% 
% Other inputs:
% ndm: Boolean that determines the default parameter set (turns transport,
% 2-species, aggregation, etc. off if set to "1")
%
% Other constants (not inputs)
% diff_n: Diffusivity of soluble tau 
% diff_m: Diffusivity of insoluble tau
% v_a: Basal rate of anterograde tau transport
% v_r: Basal rate of retrograde tau transport
% frac: Net fraction of tau undergoing diffusion (i.e. unbound tau)
%
% This function relies upon the MATLAB builtin pdepe, which solves
% elliptical and parabolic partial differential equations, and inputParser,
% which allows for more flexible and readable inputs to the function.
%
% % % 1. Set default values and initialize inputParser

alpha_ = 0; % Not currently explored in the model
beta_ = 1e-06;
gamma_ = 2e-05;
delta_ = 1;
epsilon_ = 0.01;
lambda_ = 0.01;
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
T_ = 5e7;
tsteps_ = 1000;
resmesh_ = 'fine';

ip = inputParser;
validScalar = @(x) isnumeric(x) && isscalar(x) && (x>=0);
validNonnegative = @(x) isnumeric(x) && (sum(x>=0) == length(x));
addParameter(ip, 'alpha', alpha_, validScalar);
addParameter(ip, 'beta', beta_, validScalar);
addParameter(ip, 'gamma', gamma_, validScalar);
addParameter(ip, 'delta', delta_, validScalar);
addParameter(ip, 'epsilon', epsilon_, validScalar);
addParameter(ip, 'frac', frac_, validScalar);
addParameter(ip, 'lambda', lambda_, validScalar);
addParameter(ip, 'L_int', L_int_, validScalar);
addParameter(ip, 'L1', L1_, validScalar);
addParameter(ip, 'L2', L2_, validScalar);
addParameter(ip, 'n0', n0_, validNonnegative);
addParameter(ip, 'm0', m0_, validNonnegative);
addParameter(ip, 'N1_0', N1_0_, validScalar);
addParameter(ip, 'N2_0', N2_0_, validScalar);
addParameter(ip, 'M1_0', M1_0_, validScalar);
addParameter(ip, 'M2_0', M2_0_, validScalar);
addParameter(ip, 'T', T_, validScalar);
addParameter(ip, 'tsteps', tsteps_, validScalar);
addParameter(ip, 'resmesh', resmesh_);
addParameter(ip, 'L_ais', L_ais_);
addParameter(ip, 'L_syn', L_syn_);
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
diff_ratio = 0; % Diffusivity of m relative to n (0.01 is a heuristic)

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

% % % 4. Other pdepe inputs
m_pde = 0; % This tells pdepe to evaluate the PDE using 1D symmetry
trange = [linspace(0,990,50),exp(linspace(log(1000),log(ip.Results.T),ip.Results.tsteps-50))];
sol = pdepe(m_pde,@pdex4pde,@pdex4ic,@pdex4bc,xmesh,trange);
n = sol(:,:,1);
m = sol(:,:,2);
jn = zeros(length(trange),length(xmesh)); jm = jn;
for t = 1:length(trange)
    [jn_t,jm_t] = fluxcalc(n(t,:),m(t,:),xmesh);
    jn(t,:) = jn_t;
    jm(t,:) = jm_t;
end

% % --------------------------------------------------------------
    function maskval = axon_mask(x)
    % Defines the spatial mask as a function of x position, returning 1 if 
    % the compartment is axonal and 0 otherwise. Used for turning active
    % transport on and off.
    
    if x < ip.Results.L1 || x > ip.Results.L1+ip.Results.L_int
        maskval = 0;
    else
        maskval = 1;
    end
    end

    function maskval = aissyn_mask(x)
    % Defines the spatial mask as a function of x position, returning 1 if
    % the compartment is in the AIS/synapse and 0 elsewhere. Used for
    % controlling the diffusivity of soluble tau through lambda.
    
    if axon_mask(x)
        if  x < ip.Results.L1+ip.Results.L_ais || x > ip.Results.L1+ip.Results.L_int-ip.Results.L_syn
            maskval = 1;
        else
            maskval = 0;
        end
    else
        maskval = 0;
    end
    end

    function maskval = syn_mask(x)
    % Defines the spatial mask as a function of x position, returning 1 if
    % the compartment is in the AIS/synapse and 0 elsewhere. Used to
    % specify if interconversion is allowed or not.
    
    if axon_mask(x)
        if  x > ip.Results.L1+ip.Results.L_int-ip.Results.L_syn
            maskval = 1;
        else
            maskval = 0;
        end
    else
        maskval = 0;
    end
    end
    function [c,f,s] = pdex4pde(x,t,u,DuDx)
    % This creates a function handle to pass to the PDE solver for the coupled 
    % set of equations stated above.
    
    diff_op = @(z) (1-aissyn_mask(z))*diff_n + aissyn_mask(z)*ip.Results.lambda*diff_n; 
    c = [1; 1];
    v = (1-aissyn_mask(x))*axon_mask(x)*(v_a*(1+ip.Results.delta*u(1))*(1-ip.Results.epsilon*u(2)) - v_r);

    f_n = ip.Results.frac*diff_op(x)*DuDx(1) - v*(1-ip.Results.frac)*u(1);
    f_m = diff_ratio*diff_op(x)*DuDx(2);
    f = [f_n; f_m];

    s_n = (1-syn_mask(x))*(ip.Results.alpha + ip.Results.beta*u(2) - ip.Results.gamma*u(1)*(u(1)+u(2)));
    s_m = (1-syn_mask(x))*(ip.Results.gamma*u(1)*(u(1)+u(2)) - ip.Results.beta*u(2));
    s = [s_n; s_m];
    end
    % --------------------------------------------------------------
    function u0 = pdex4ic(x)
    % This creates a function handle to pass to the PDE solver for the
    % initial conditions, whose values are specified above as N1, n0 and N2
    % for soluble tau and M1, m0, and M2 for insoluble tau. In order to put
    % those values, which are defined per unit of length, into the xmesh
    % space used by pdepe, spline interpolation is used.
    n_init0 = [ip.Results.N1_0*ones(1,ip.Results.L1),ip.Results.n0,ip.Results.N2_0*ones(1,ip.Results.L2)];
    m_init0 = [ip.Results.M1_0*ones(1,ip.Results.L1),ip.Results.m0,ip.Results.M2_0*ones(1,ip.Results.L2)];
    n_init = interp1(1:L_total,n_init0,xmesh,'makima');
    m_init = interp1(1:L_total,m_init0,xmesh,'makima');
    icind = find(xmesh == x);
    u0 = [n_init(icind); m_init(icind)];
    u0(u0 < 0) = 0;
    end
    % --------------------------------------------------------------
    function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t)
    % This creates a function handle to pass to the PDE solver for the
    % boundary conditions. These should be perfectly insulating/reflective
    % for the purposes of this problem.
    
    pl = [0; 0];
    pr = [0; 0];
    ql = [1; 1]; 
    qr = [1; 1];
    end

    function [Jn,Jm] = fluxcalc(n,m,xrange)
    % Calculates the approximate flux for a given 1D concentration profile
    % (prof) over its domain (xrange), used for plotting.
    Jn = zeros(1,length(xrange)); Jm = Jn;
    diff_op_n = @(z)(1-aissyn_mask(z))*diff_n + aissyn_mask(z)*ip.Results.lambda*diff_n;
    diff_op_m = @(z)(1-aissyn_mask(z))*diff_n*diff_ratio + aissyn_mask(z)*ip.Results.lambda*diff_n*diff_ratio;
    for i = 1:length(xrange)
        if i == 1 || i == length(xrange)
            Jn(i) = 0; Jm(i) = 0;
        else
            v = (1-aissyn_mask(xrange(i)))*axon_mask(xrange(i))*(v_a*(1+ip.Results.delta*n(i))*(1-ip.Results.epsilon*m(i)) - v_r);
            DnDx = (n(i+1) - n(i-1))/(xrange(i+1) - xrange(i-1));
            DmDx = (m(i+1) - m(i-1))/(xrange(i+1) - xrange(i-1));
            Jn(i) = -ip.Results.frac*diff_op_n(xrange(i))*DnDx + v*(1-ip.Results.frac)*n(i);
            Jm(i) = -diff_op_m(xrange(i))*DmDx;
        end
    end
    end
end
