cmap = lines(size(output_struct.Parameter_Grid,1));
T = output_struct.Simulations(1).Model_Outputs.Sim.T; 
dt = output_struct.Simulations(1).Model_Outputs.Sim.dt; 
trange = 0:dt:T;

for i = 1:size(output_struct.Parameter_Grid,1)
    figure('Position',[0,0,500,1000]); tiledlayout(2,1);
    nexttile; hold on;
    N = output_struct.Simulations(i).Model_Outputs.Predicted.N;
    M = output_struct.Simulations(i).Model_Outputs.Predicted.M;
    for j = 1:size(N,1)
        plot(trange,N)
    end
    xlabel('t'); ylabel('N');
    txt = ['$\mathbf{\lambda_1} = $' num2str(output_struct.Simulations(i).Model_Outputs.Parameters.lambda1) ','...
        '$\mathbf{\lambda_2} = $' num2str(output_struct.Simulations(i).Model_Outputs.Parameters.lambda2) ',' , ......
        '$\mathbf{\beta} = $' num2str(output_struct.Simulations(i).Model_Outputs.Parameters.beta) ',' , ...
        '$\mathbf{\gamma1}=$' num2str(output_struct.Simulations(i).Model_Outputs.Parameters.gamma1) ',',...
        '$\mathbf{\epsilon}=$' num2str(output_struct.Simulations(i).Model_Outputs.Parameters.epsilon) ',',...
        '$\mathbf{\delta}=$' num2str(output_struct.Simulations(i).Model_Outputs.Parameters.delta)];
    title(txt,'Interpreter','latex');
    nexttile; hold on;
    for j = 1:size(M,1)
        plot(trange,M)
    end
    xlabel('t'); ylabel('M');
end
