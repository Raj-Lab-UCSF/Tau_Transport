% Simulations for paper:
% (Rerun both hippocampome ones with same init path and t range)
% 
% 1. lambda_1 = lambda_2, but different values
% 2. anterograde vs. retrograde bias spread
% 3. different gamma1 (aggregation)
% 
% Plots for paper:
% 1a. Line plots (maybe remove legend because it looks bad)
% 1b. Heatmap plots 
% 1c. Brainframe images at select times
% 2. Correlation with seed plot
% 3. Pairwise correlations between simulations
% 4. Show comparison between this model and network diffusion model with
% different s parameters (maybe include alpha as well?)

figpath = '~/Documents/MATLAB/Tau_Transport/OutputFigures';
simpath = '~/Documents/MATLAB/Tau_Transport/SampleFiles';
loadpath = '~/Documents/MATLAB/Tau_Transport/MatFiles';
bfpath = '~/Documents/MATLAB/Brainframe-Dev/Brainframe';
simstr = 'hippocampome_final_lambda';
load([simpath filesep simstr '_1.mat'],'output_struct');

% [Co,Ci] = ConnectomePlot(simstr,1,loadpath,simpath,0,figpath);
TimeCoursePlot([simstr '_1'],21,'Heatmap',1,loadpath,simpath,0,figpath);
TimeCoursePlot([simstr '_1'],7,'Heatmap',1,loadpath,simpath,0,figpath);
% BrainframePlot(simstr,2,[1,10,20,37],1,0,loadpath,simpath,bfpath,figpath);
% BrainframePlot(simstr,3,[1,10,20,37],1,0,loadpath,simpath,bfpath,figpath);

TimeCoursePlot([simstr '_1'],21,'Line',0,loadpath,simpath,0,figpath);
TimeCoursePlot([simstr '_1'],7,'Line',0,loadpath,simpath,0,figpath);
% TimeCoursePlot([simstr '_2'],2,'Line',0,loadpath,simpath,0,figpath);
% TimeCoursePlot([simstr '_2'],3,'Line',0,loadpath,simpath,0,figpath);


inds = [1:24 26:30];
ts = output_struct.Simulations(2).Model_Outputs.Sim.trange * 180;
tot_ret = output_struct.Simulations(2).Model_Outputs.Predicted.N(inds,:) + ...
    output_struct.Simulations(2).Model_Outputs.Predicted.M(inds,:);
tot_ant = output_struct.Simulations(3).Model_Outputs.Predicted.N(inds,:) + ...
    output_struct.Simulations(3).Model_Outputs.Predicted.M(inds,:);
corrs_ret_Co = corr(Co,tot_ret);
corrs_ret_Ci = corr(Ci,tot_ret);
corrs_ant_Co = corr(Co,tot_ant);
corrs_ant_Ci = corr(Ci,tot_ant);


figure('Position',[0,0,1000,600]); 
subplot(1,2,1); hold on;
plot(ts,corrs_ant_Ci,'LineWidth',3,'Color','red');
plot(ts,corrs_ant_Co,'LineWidth',3,'Color','blue');
xlabel('Time Days'); ylabel("Pearson's R with C");
legend({'C_i_n, seed','C_o_u_t, seed'},'Location','northeast');
title('Ant, small gamma Simulation');
set(gca,'FontSize',16,'FontName','Times')
subplot(1,2,2); hold on;
plot(ts,corrs_ret_Ci,'LineWidth',3,'Color','red');
plot(ts,corrs_ret_Co,'LineWidth',3,'Color','blue');
xlabel('Time Days'); ylabel("Pearson's R with C");
legend({'C_i_n, seed','C_o_u_t, seed'},'Location','northeast');
title('Ret small gamma Simulation');
set(gca,'FontSize',16,'FontName','Times')

params = output_struct.Parameter_Grid;
params_ret = params((params(:,7) < params(:,8)),:);
params_ant = params((params(:,7) > params(:,8)),:);
retinds = find((params(:,7) < params(:,8)));
antinds = find((params(:,7) > params(:,8)));
cmap = lines(size(params_ret,1));
leg = cell(1,size(params_ret,1));
figure('Position',[0,0,600,600]); hold on;
for i = 1:size(params_ret,1)
    tot_ret = output_struct.Simulations(retinds(i)).Model_Outputs.Predicted.N(inds,:) + ...
        output_struct.Simulations(retinds(i)).Model_Outputs.Predicted.M(inds,:);
    tot_ant = output_struct.Simulations(antinds(i)).Model_Outputs.Predicted.N(inds,:) + ...
        output_struct.Simulations(antinds(i)).Model_Outputs.Predicted.M(inds,:);
    corrs = corr(tot_ret,tot_ant); corrs = diag(corrs);
    plot(ts,corrs,'LineWidth',3,'Color',cmap(i,:))
    leg{i} = ['$\gamma$ = ',num2str(params_ret(i,2),'%.2d')];
end
xlabel('Time (Days)'); ylabel("Pearson's R, Ant vs. Ret");
legend(leg,'Location','northeast','Interpreter','latex');
set(gca,'FontSize',16,'FontName','Times')

 