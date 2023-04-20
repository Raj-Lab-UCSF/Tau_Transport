function model_BFRTT()

%load("220919_MouseDataforVeronica.mat");
load 220919_MouseDataforVeronica.mat  baseline_pathology; 
Conn=readmatrix('mouse_connectome.csv');
nroi=size(Conn,1);
Adj=readmatrix('mouse_adj_matrix.csv');
beta = 1e-06;
gamma= 2e-05;

init_tau=baseline_pathology.Hurtado;
init_tau=1e-4*init_tau;


Conn_in=(zeros(nroi));
for i=1:nroi

C_adj_in=Conn(:,i)>0;
CC_in=Conn(C_adj_in,i);
  for j=1:length(CC_in)
      Conn_in(i,j)=CC_in(j); %the incoming adiancency are set for rows
  end
end
Conn_out=(zeros(nroi));
for i=1:nroi

C_adj_out=Conn(i,:)>0;
CC_out=Conn(i,C_adj_out);
  for j=1:length(CC_out)
      Conn_out(i,j)=CC_out(j); %the outcoming adiancency are set for rows
  end
end
 k=0.001;
% %
  t=0:k:0.01;
  nt=length(t);
  n=zeros(nroi,nt);
  m=zeros(nroi,nt);
  n(:,1)=init_tau(:);
 A_in= zeros([nroi,size(n)]);
 A_out=zeros([nroi,size(n)]);

%h=1;
% 
for h=1:(nt-1)

 n_adj_out=n(:,h).*Adj.'; %outcoming edges   %for columns

 n_adj_in=n(:,h).*Adj; %incoming edges %for columns

 n_adj_out=n_adj_out.';

  A_in(:,:,h)=flux_calculator11_11(n_adj_in,n(:,h),1);

  A_out(:,:,h)=flux_calculator11_11(n(:,h),n_adj_out,0); 
 
 n(:,h+1)=n(:,h)+(diag((A_in(:,:,h)*Conn_in.')) - diag(A_out(:,:,h)*Conn_out.'))*k+ beta*m(:,h)*k-gamma*(n(:,h).*n(:,h)+n(:,h).*m(:,h))*k;

 
 m(:,h+1)=m(:,h)-beta*m(:,h)*k+gamma*(n(:,h).*(n(:,h)+m(:,h)))*k;
 end
 
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

subplot(2,1,2)
plot(t,m);
xlabel('t');
ylabel('m(t)');
end
