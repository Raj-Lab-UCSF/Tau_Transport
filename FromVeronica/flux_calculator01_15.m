function A_1=flux_calculator01_15(b,c,varargin)

alpha_ = 0; % Not currently explored in the model
beta_ = 1e-06;
gamma1_ = 2e-05;
gamma2_=2e-06;
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
n0_ = zeros(1, L_int_);
m0_ = zeros(1, L_int_);
N1_0_ = 0; 
N2_0_ = 0; 
M1_0_ = 0; 
M2_0_ = 0;
resmesh_ = 'coarse';
total_mass_=184;
ip = inputParser;
validScalar = @(x) isnumeric(x) && isscalar(x) && (x>=0);
validNonnegative = @(x) isnumeric(x) && (sum(x>=0) == length(x));
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
addParameter(ip, 'n0', n0_, validNonnegative);
addParameter(ip, 'm0', m0_, validNonnegative);
addParameter(ip, 'N1_0', N1_0_, validScalar);
addParameter(ip, 'N2_0', N2_0_, validScalar);
addParameter(ip, 'M1_0', M1_0_, validScalar);
addParameter(ip, 'M2_0', M2_0_, validScalar);
addParameter(ip, 'resmesh', resmesh_);
addParameter(ip, 'L_ais', L_ais_);
addParameter(ip, 'L_syn', L_syn_);
addParameter(ip, 'total_mass',total_mass_);
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

%Network stuff
Adj=readmatrix('mouse_adj_matrix_19_01.csv');
%Conn=readmatrix('mouse_connectome.csv');
nroi=size(Adj,1);


% % % 5a. Steady State of Presynaptic Somatodendritic Compartment
% Follows the derivation of Michiel Bertsh
presyn_mask = spatial_mask('presyn');

xmesh_presyn = xmesh(presyn_mask);
%n_ss_presyn = @(A,B,x) max((B - A.*x/diff_n),0); 
%m_ss_presyn = @(A,B,x) ((ip.Results.gamma*n_ss_presyn(A,B,x).^2)./(ip.Results.beta - n_ss_presyn(A,B,x)*ip.Results.gamma));
%m_ss_presyn = @(A,x) ((ip.Results.gamma*n_ss_presyn(A,x).^2)./(ip.Results.beta));


% % % 5b. Steady state of axon initial segment
ais_mask = spatial_mask('ais');

xmesh_ais= xmesh(ais_mask);
x1 = xmesh_presyn(end);
%n_ss_ais = @(A,B,x) max((B - A.*x1/diff_n - A.*(x-x1)/(diff_n*ip.Results.lambda1)),0);
%n_ss_ais=@(A,B,x)max((n_ss_presyn(A,B,x1)-A*(x-x1)/(diff_n*ip.Results.lambda1)),0);
%m_ss_ais = @(A,B,x) (ip.Results.gamma*n_ss_ais(A,B,x).^2)./(ip.Results.beta - n_ss_ais(A,B,x)*ip.Results.gamma);
      %m_ss_ais = @(A,x) (ip.Results.gamma*n_ss_ais(A,x).^2)./(ip.Results.beta);

 % % 5c. Steady state of axon
      axon_mask = spatial_mask('axon');

      xmesh_axon= xmesh(axon_mask);
      x2 = xmesh_ais(end);
      nx_init = @(A_,B) max((-A_.*x1/diff_n + B - A_.*(x2 - x1)/(ip.Results.lambda1*diff_n)),0);
      %nx_init=@(A_,B)  n_ss_ais(A_,B,x2);
      options = odeset('RelTol',1e-4,'AbsTol',1e-4,'NonNegative',1:length(nx_init));


     n_ss_axon = @(A,B) ode45(@(x,n)ode_ss_axon(x,A,n),[ip.Results.L1+ip.Results.L_ais,...
           ip.Results.L1+ip.Results.L_int-ip.Results.L_syn],nx_init(A,B),options);

      n_ss_axon= @(A,B,x) deval(n_ss_axon(A,B),x);
%       m_ss_axon = @(A,B,x) (ip.Results.gamma*n_ss_axon(A,B,x).^2)./ ...
%                   (ip.Results.beta - ip.Results.gamma*n_ss_axon(A,B,x));
      %m_ss_axon = @(A,x) (ip.Results.gamma*n_ss_axon(A,x).^2)./ ...
               % (ip.Results.beta);
      

      syncleft_mask = spatial_mask('syncleft');

      xmesh_syncleft = xmesh(syncleft_mask);
      x3 = xmesh_axon(end);
      n_ss_syncleft = @(A,B,x) max((n_ss_axon(A,B,x3) - A.*(x-x3)/(diff_n*ip.Results.lambda2)),0);
      %m_ss_syncleft = M1_0_*ones(1,length(xmesh_syncleft));
    

% % % 5e. Steady state of postsynaptic somatodendritic compartment
      postsyn_mask = spatial_mask('postsyn');

      xmesh_postsyn = xmesh(postsyn_mask);
      x4 = xmesh_syncleft(end);
      n_ss_postsyn = @(A,B,x) max((n_ss_syncleft(A,B,x4) - A.*x/diff_n),0); 
     % m_ss_postsyn = @(A,B,x) (ip.Results.gamma*n_ss_postsyn(A,B,x).^2)./ ...
               % (ip.Results.beta - ip.Results.gamma*n_ss_postsyn(A,B,x));
     % m_ss_postsyn = @(A,x) (ip.Results.gamma*n_ss_postsyn(A,x).^2)./ ...
      %          (ip.Results.beta);
      x5=xmesh_postsyn(end);

%f_ss=@(A)(n_ss_axon(A,x3)- (A/(diff_n*ip.Results.lambda2))*((1-ip.Results.lambda2)*x4+ip.Results.lambda2*x5-x3)-C_1);

 
      % f_ss=@(A,B,C)(n_ss_postsyn(A,B,x5)-C);
       f_ss=@(A,B,C)(n_ss_syncleft(A,B,x4) - A.*x5/diff_n-C);

     
 A_1=zeros(nroi);

for i=1:nroi
  
      Ad_in=logical(Adj(:,i));
      C_1=c(i);
    C_1=repmat(C_1,1,length(Ad_in));
   C_1=C_1.';
   i_app=logical(((b(:,i)>0)+(C_1>0)).*(Ad_in));
   B_1=b(i_app,i);
   C_1=C_1(i_app);
   if length(B_1)>0
    x0=zeros(length(B_1),1);
%     A_sample=(1e-8+1e-8).*rand(length(B_1),10)-1e-8; % N random numbers in the interval (a,b) with the formula r = a + (b-a).*rand(N,1)
%     A_sample=[A_0 A_sample];
%      f_test=zeros(11,1);
%       for m=1:11
%           A_test=A_sample(:,m);
%           %X=f_ss(A_test,B_2,C_2);
%           %size(X)
%           f_test(m)=norm(f_ss(A_test,B_1,C_1));
%       end
%       [m_test,i_test]=min(f_test,[],'all');
%       x0=A_sample(:,i_test);


%     if isempty(B_1>0)==1 && isempty(C_1>0)==1
%         A_1(i,:)=0;  
%     end

        
          %options=optimset('tolX',1e-6,'TolFun',1e-6);
          options=optimset('TolFun',1e-6);
         fun_ss=@(A)f_ss(A,B_1,C_1);
           A_1(i_app,i)=fsolve(fun_ss,x0,options);
   end
end
          

% n_ss_presyn=n_ss_presyn(A_1,xmesh_presyn);
% m_ss_presyn=m_ss_presyn(A_1, xmesh_presyn);
% n_ss_ais=n_ss_ais(A_1,xmesh_ais);
% m_ss_ais=m_ss_ais(A_1,xmesh_ais);
% n_ss_axon=n_ss_axon(A_1,xmesh_axon);
% m_ss_axon=m_ss_axon(A_1,xmesh_axon);
% n_ss_syncleft=n_ss_syncleft(A_1,xmesh_syncleft);
%m_ss_syncleft=m_ss_syncleft(xmesh_syncleft.');


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

 function nprime=ode_ss_axon(x,A,n) 
 nprime=1/diff_n*-A./ip.Results.frac+1/diff_n*((1-ip.Results.frac)./ip.Results.frac).*n.*(v_a*(1+ip.Results.delta.*n).*(1-((ip.Results.gamma1*ip.Results.epsilon.*n.^2)./(ip.Results.beta-ip.Results.gamma2.*n)))-v_r);
 %nprime=1/diff_n*-A./ip.Results.frac+1/diff_n*((1-ip.Results.frac)./ip.Results.frac).*n.*(v_a*(1+ip.Results.delta.*n).*(1-((ip.Results.gamma*ip.Results.epsilon.*n.^2)./(ip.Results.beta)))-v_r);
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