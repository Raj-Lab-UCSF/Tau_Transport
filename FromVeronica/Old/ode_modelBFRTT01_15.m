function ode_modelBFRTT01_15()

%load("220919_MouseDataforVeronica.mat");
load 220919_MouseDataforVeronica.mat  baseline_pathology; 
Conn=readmatrix('mouse_connectome_19_01.csv');
%Conn=csvread('mean80_fibercount_Ashish.csv', 1, 0);
%Adj=readmatrix('Adj_matrix_Ashish.csv');
nroi=size(Conn,1);
Adj=readmatrix('mouse_adj_matrix_19_01.csv');
beta = 1e-06;
gamma= 2e-05;
delta_ = 1; %the parameters delta,epsilon,lambda1, lambda2 are rewritten here just for plotting
epsilon_ = 0.01;
lambda1_ = 0.01;
lambda2_ = 0.01;%0.01;
initial_tau=baseline_pathology.Hurtado;
%init_tau=ones(86,1);
%init_tau([2,23:25, 45:48, 60, 79, 81] )=1;
%init_tau=1e-4*initial_tau;
init_tau=1e-4*initial_tau;
i_nonzero_init_tau=init_tau>0;
i_zero=logical(1-initial_tau);

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


 k=0.005;
% %
  t=0:k:0.05;
  nt=length(t);
  n=zeros(nroi,nt);
  m=zeros(nroi,nt);
  n(:,1)=init_tau(:);
 A_in= zeros([nroi,size(n)]);
 
 





for h=1:(nt-1)

 
 n_adj_in=n(:,h).*Adj; %incoming edges %for columns


  A_in(:,:,h)=flux_calculator01_15(n_adj_in,n(:,h));
  %A_0=A_in(:,:,h);
  %A_0=A_0(Adj>0);
  
 F_in=A_in(:,:,h);
F_out=F_in.';

 n(:,h+1)=n(:,h)+(diag((Conn.'*F_in)) - diag((Conn*F_out)))*k + beta*m(:,h)*k-gamma*(n(:,h).*n(:,h)+n(:,h).*m(:,h))*k;

%n(:,h+1)=max(n(:,h+1),1e-8*ones(nroi,1));

 m(:,h+1)=m(:,h)-beta*m(:,h)*k+gamma*(n(:,h).*(n(:,h)+m(:,h)))*k;
end

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
%end
figure
subplot(2,1,1)
plot(t,n);


xlabel('t');
ylabel('n(t)');
 title("n,m distributions on the network",'Fontsize',12);
 txt = ['$\mathbf{\lambda_1} = $' num2str(lambda1_) ',' '$\mathbf{\lambda_2} = $' num2str(lambda2_) ',' , '$\mathbf{\epsilon}=$' num2str(epsilon_) ',', '$\mathbf{\delta}=$' num2str(delta_) ];
 subtitle(txt,'Interpreter','latex');
subplot(2,1,2)
plot(t,m);
xlabel('t');
ylabel('m(t)');


%[M,i]=min(abs((n(:,2)-3.744e-05)),[],'all')
[M,i]=min(n,[],'all')
N=length(n(n>0))

figure
subplot(2,1,1)
plot(t,n(i_nonzero_init_tau,:));
subplot(2,1,2)
plot(t,m(i_nonzero_init_tau,:));


end
