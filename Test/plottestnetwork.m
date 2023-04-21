N = sim_struct(4).Model_Outputs.Predicted.N;
M = sim_struct(4).Model_Outputs.Predicted.M;
t = 0:sim_struct(4).Model_Outputs.Sim.dt:sim_struct(1).Model_Outputs.Sim.T;
lambda1 = sim_struct(4).Model_Outputs.Parameters.lambda1;
lambda2 = sim_struct(4).Model_Outputs.Parameters.lambda2;
beta = sim_struct(4).Model_Outputs.Parameters.beta;
gamma1 = sim_struct(4).Model_Outputs.Parameters.gamma1;
delta = sim_struct(4).Model_Outputs.Parameters.delta;
epsilon = sim_struct(4).Model_Outputs.Parameters.epsilon;

figure
subplot(2,1,1)
plot(t,N);
xlabel('t');
ylabel('N(t)');
title("N,M distributions on the network",'Fontsize',12);
txt = ['$\mathbf{\lambda_1} = $' num2str(lambda1) ','...
    '$\mathbf{\lambda_2} = $' num2str(lambda2) ',' , ......
    '$\mathbf{\beta} = $' num2str(beta) ',' , ...
    '$\mathbf{\gamma_1}=$' num2str(gamma1) ',',...
    '$\mathbf{\epsilon}=$' num2str(epsilon) ',',...
    '$\mathbf{\delta}=$' num2str(delta)];
subtitle(txt,'Interpreter','latex');
subplot(2,1,2)
plot(t,M);
xlabel('t');
ylabel('M(t)');