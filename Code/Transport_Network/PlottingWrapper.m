% Simulations for paper:
% (Rerun both hippocampome ones with same init path and t range)
% 
% 1. lambda_1 = lambda_2, but different values
% 2. anterograde vs. retrograde bias spread
% 3. different gamma1 (aggregation)
% 
% Plots for paper (old plan):
% 1a. Line plots (maybe remove legend because it looks bad)
% 1b. Heatmap plots 
% 1c. Brainframe images at select times
% 2. Correlation with seed plot
% 3. Pairwise correlations between simulations
% 4. Show comparison between this model and network diffusion model with
% different s parameters (maybe include alpha as well?)
%
% Plots for paper (new plan):
% 1. Show that steady state calculation lines up with numerical
% approximations from the paper (same parameter values as 2021)
% 2. Show high/low lambda for intermediate values of other parameters
% 3. Show high/low delta for intermediate values of other parameters
% 4. Show high/low epsilon for intermediate values of other parameters
% 5. Show high/low gamma1 for intermediate values of other parameters
% 6. Selected Brainframe images

%% 0. Load files, set paths, etc.
% figpath = '~/Documents/MATLAB/Tau_Transport_OtherFiles/Figures';
% simpath = '~/Documents/MATLAB/Tau_Transport_OtherFiles/FinalSimFiles';
figpath = '~/Documents/MATLAB/SampleFiles'; simpath = figpath;
loadpath = '~/Documents/MATLAB/Tau_Transport/MatFiles';
bfpath = '~/Documents/MATLAB/Brainframe-Dev/Brainframe';
simstr = 'hippocampome_final_round2';
load([simpath filesep simstr '.mat'],'output_struct');

load([loadpath filesep 'DefaultAtlas.mat'],'DefaultAtlas');
V_inv = 1./DefaultAtlas.volumes; V_inv = diag(V_inv);
load([loadpath filesep 'CCF_labels.mat'],'CCF_labels');
switch output_struct.Simulations(1).Model_Outputs.Sim.connectome_subset
    case 'Hippocampus'
        inds = ismember(CCF_labels(:,3),'Hippocampus');
    case 'Hippocampus+PC+RSP'
        inds_hipp = ismember(CCF_labels(:,3),'Hippocampus');
        inds_pc = ismember(CCF_labels(:,1),'Piriform area');
        inds_rsp = ismember(CCF_labels(:,3),'Retrosplenial Area');
        inds = logical(inds_hipp + inds_pc + inds_rsp);
    case 'RH'
        inds = ismember(CCF_labels(:,4),'Right Hemisphere');
    case 'LH'
        inds = ismember(CCF_labels(:,4),'Left Hemisphere');
end
C = output_struct.Simulations(1).Model_Outputs.Sim.C;
V_inv = V_inv(inds,inds);

%% 1. Check mass conservation
MassConservationPlot(simstr,loadpath,simpath,1,figpath);

%% 2. Time course plots
colind_lambda = find(ismember(output_struct.Parameter_Names,'lambda1'));
colind_gamma1 = find(ismember(output_struct.Parameter_Names,'gamma1'));
colind_delta = find(ismember(output_struct.Parameter_Names,'delta'));
colind_epsilon = find(ismember(output_struct.Parameter_Names,'epsilon'));

lambda_vals = unique(output_struct.Parameter_Grid(:,colind_lambda));
gamma1_vals = unique(output_struct.Parameter_Grid(:,colind_gamma1));
delta_vals = unique(output_struct.Parameter_Grid(:,colind_delta));
epsilon_vals = unique(output_struct.Parameter_Grid(:,colind_epsilon));

%% 2.1 lambda1/2
gamma1_midvals = gamma1_vals; gamma1_midvals([1,length(gamma1_vals)]) = [];
delta_midvals = delta_vals; delta_midvals([1,length(delta_vals)]) = [];
epsilon_midvals = epsilon_vals; epsilon_midvals([1,length(epsilon_vals)]) = [];
bool_midgamma1 = ismember(output_struct.Parameter_Grid(:,colind_gamma1),gamma1_midvals);
bool_middelta = ismember(output_struct.Parameter_Grid(:,colind_delta),delta_midvals);
bool_midepsilon = ismember(output_struct.Parameter_Grid(:,colind_epsilon),epsilon_midvals);
bool_allmids = bool_midgamma1 + bool_middelta + bool_midepsilon;
bool_allmids = bool_allmids == max(bool_allmids); 

%% 2.1.1 Maximum value
maxlambda = max(lambda_vals);
bool_maxlambda = ismember(output_struct.Parameter_Grid(:,colind_lambda),maxlambda);
bool_plotmax = bool_allmids + bool_maxlambda;
bool_plotmax = bool_plotmax == max(bool_plotmax); 
idxs_plotmax = find(bool_plotmax);

%% 2.1.2 Minimum value
minlambda = min(output_struct.Parameter_Grid(:,colind_lambda));
bool_minlambda = ismember(output_struct.Parameter_Grid(:,colind_lambda),minlambda);
bool_plotmin = bool_allmids + bool_minlambda;
bool_plotmin = bool_plotmin == max(bool_plotmin); 
idxs_plotmin = find(bool_plotmin);

%% 2.1.3 Plotting
for i = 1:length(idxs_plotmin)
%     TimeCoursePlot(simstr,idxs_plot(i),'Heatmap',1,loadpath,simpath,0,figpath);
    TimeCoursePlot(simstr,idxs_plotmax(i),'Line',0,loadpath,simpath,1,figpath);
    TimeCoursePlot(simstr,idxs_plotmin(i),'Line',0,loadpath,simpath,1,figpath);
end

%% 2.2 gamma1
lambda_midvals = lambda_vals; lambda_midvals([1,length(lambda_vals)]) = [];
delta_midvals = delta_vals; delta_midvals([1,length(delta_vals)]) = [];
epsilon_midvals = epsilon_vals; epsilon_midvals([1,length(epsilon_vals)]) = [];
bool_midlambda = ismember(output_struct.Parameter_Grid(:,colind_lambda),lambda_midvals);
bool_middelta = ismember(output_struct.Parameter_Grid(:,colind_delta),delta_midvals);
bool_midepsilon = ismember(output_struct.Parameter_Grid(:,colind_epsilon),epsilon_midvals);
bool_allmids = bool_midlambda + bool_middelta + bool_midepsilon;
bool_allmids = bool_allmids == max(bool_allmids); 

%% 2.2.1 Maximum value
maxgamma1 = max(gamma1_vals);
bool_maxgamma1 = ismember(output_struct.Parameter_Grid(:,colind_gamma1),maxgamma1);
bool_plotmax = bool_allmids + bool_maxgamma1;
bool_plotmax = bool_plotmax == max(bool_plotmax); 
idxs_plotmax = find(bool_plotmax);

%% 2.2.2 Minimum value
mingamma1 = min(output_struct.Parameter_Grid(:,colind_gamma1));
bool_mingamma1  = ismember(output_struct.Parameter_Grid(:,colind_gamma1),mingamma1);
bool_plotmin = bool_allmids + bool_mingamma1;
bool_plotmin = bool_plotmin == max(bool_plotmin); 
idxs_plotmin = find(bool_plotmin);

%% 2.2.3 Plotting
for i = 1
%     TimeCoursePlot(simstr,idxs_plot(i),'Heatmap',1,loadpath,simpath,0,figpath);
    TimeCoursePlot(simstr,idxs_plotmax(i),'Line',0,loadpath,simpath,1,figpath);
    TimeCoursePlot(simstr,idxs_plotmin(i),'Line',0,loadpath,simpath,1,figpath);
end

%% 3. Plot connectivity to/from seed region
SeedConnectivityPlot(simstr,1,loadpath,simpath,0,figpath);

%% 4. Plot correlations w.r.t. seed connectivity
lambda = 0.025;
beta = 1e-6;
frac = 0.92;
gamma1 = 2e-3;
gamma2 = 0;
delta = 10;
epsilon = 10;

%% 5. Brainframe plots
BrainframePlot(simstr,1,150,1,0,loadpath,simpath,bfpath,figpath)

%%
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

 