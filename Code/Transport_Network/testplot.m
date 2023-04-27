cmap = lines(size(sim_struct(1).Model_Outputs.Predicted.N,1));
T = sim_struct(1).Model_Outputs.Sim.T; dt = sim_struct(1).Model_Outputs.Sim.dt;
trange = 0:dt:T;

for i = 1:9
    figure('Position',[0,0,500,1000]); tiledlayout(2,1);
    nexttile; hold on;
    N = sim_struct(i).Model_Outputs.Predicted.N;
    M = sim_struct(i).Model_Outputs.Predicted.M;
    for j = 1:size(N,1)
        plot(trange,N)
    end
    xlabel('t'); ylabel('N');
    txt = ['$\mathbf{\lambda_1} = $' num2str(sim_struct(i).Model_Outputs.Parameters.lambda1) ','...
        '$\mathbf{\lambda_2} = $' num2str(sim_struct(i).Model_Outputs.Parameters.lambda2) ',' , ......
        '$\mathbf{\beta} = $' num2str(sim_struct(i).Model_Outputs.Parameters.beta) ',' , ...
        '$\mathbf{\gamma1}=$' num2str(sim_struct(i).Model_Outputs.Parameters.gamma1) ',',...
        '$\mathbf{\epsilon}=$' num2str(sim_struct(i).Model_Outputs.Parameters.epsilon) ',',...
        '$\mathbf{\delta}=$' num2str(sim_struct(i).Model_Outputs.Parameters.delta)];
    title(txt,'Interpreter','latex');
    nexttile; hold on;
    for j = 1:size(M,1)
        plot(trange,M)
    end
    xlabel('t'); ylabel('M');
end
