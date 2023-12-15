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
figpath = '~/Documents/MATLAB/OutputFigures/Tau_Transport/Network_Transport';
simpath = '~/Documents/MATLAB/Tau_Transport_OtherFiles/FinalSimFiles';
% figpath = '~/Documents/MATLAB/SampleFiles'; simpath = figpath;
loadpath = '~/Documents/MATLAB/Tau_Transport/MatFiles';
bfpath = '~/Documents/MATLAB/Brainframe-Dev/Brainframe';
simstr = 'hippocampome_final_round2_v2';
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

%% 1. Check that steady state matches end-timepoint tau distributions
% See SingleEdgeSims.m for simulation script
simstr_edge = 'constant_n0_L1000_fine_ant_neumann_pde_ss';
simstr_ss = '~/Documents/Tau_Transport/SampleFiles';
savenclose = 1;
SteadyStatePlotter(simstr_edge,simpath_ss,savenclose,figpath);

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
htmapplt = 1;
lnplt = 1;
savenclose = 1;
figpathlam = [figpath filesep 'LambdaCompare'];
for i = 1:length(idxs_plotmin)
    if lnplt
        TimeCoursePlot(simstr,idxs_plotmax(i),'Line',0,loadpath,simpath,savenclose,figpathlam); %#ok<UNRCH> 
        TimeCoursePlot(simstr,idxs_plotmax(i),'Line',1,loadpath,simpath,savenclose,figpathlam);
        TimeCoursePlot(simstr,idxs_plotmin(i),'Line',0,loadpath,simpath,savenclose,figpathlam);
        TimeCoursePlot(simstr,idxs_plotmin(i),'Line',1,loadpath,simpath,savenclose,figpathlam);
    end
    if htmapplt
        TimeCoursePlot(simstr,idxs_plotmax(i),'Heatmap',1,loadpath,simpath,savenclose,figpathlam); %#ok<UNRCH> 
        TimeCoursePlot(simstr,idxs_plotmin(i),'Heatmap',1,loadpath,simpath,savenclose,figpathlam);
    end
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
htmapplt = 1;
lnplt = 1;
savenclose = 1;
figpathgam = [figpath filesep 'Gamma1Compare'];
for i = 1:length(idxs_plotmin)
    if lnplt
        TimeCoursePlot(simstr,idxs_plotmax(i),'Line',0,loadpath,simpath,savenclose,figpathgam); %#ok<UNRCH> 
        TimeCoursePlot(simstr,idxs_plotmax(i),'Line',1,loadpath,simpath,savenclose,figpathgam);
        TimeCoursePlot(simstr,idxs_plotmin(i),'Line',0,loadpath,simpath,savenclose,figpathgam);
        TimeCoursePlot(simstr,idxs_plotmin(i),'Line',1,loadpath,simpath,savenclose,figpathgam);
    end
    if htmapplt
        TimeCoursePlot(simstr,idxs_plotmax(i),'Heatmap',1,loadpath,simpath,savenclose,figpathgam); %#ok<UNRCH> 
        TimeCoursePlot(simstr,idxs_plotmin(i),'Heatmap',1,loadpath,simpath,savenclose,figpathgam);
    end
end

%% 2.3 Anterograde/Retrograde
lambda_midvals = lambda_vals; lambda_midvals([1,length(lambda_vals)]) = [];
gamma1_midvals = gamma1_vals; gamma1_midvals([1,length(gamma1_vals)]) = [];
bool_midlambda = ismember(output_struct.Parameter_Grid(:,colind_lambda),lambda_midvals);
bool_midgamma1 = ismember(output_struct.Parameter_Grid(:,colind_gamma1),gamma1_midvals);
bool_allmids = bool_midlambda + bool_midgamma1;
bool_allmids = bool_allmids == max(bool_allmids); 

%% 2.3.1 Anterograde
maxdelta = max(delta_vals); minepsilon = min(epsilon_vals);
bool_maxdelta = ismember(output_struct.Parameter_Grid(:,colind_delta),maxdelta);
bool_minepsilon = ismember(output_struct.Parameter_Grid(:,colind_epsilon),minepsilon);
bool_plotant = bool_allmids + bool_maxdelta + bool_minepsilon;
bool_plotant = bool_plotant == max(bool_plotant); 
idxs_plotant = find(bool_plotant);

%% 2.3.2 Retrograde
mindelta = min(delta_vals); maxepsilon = max(epsilon_vals);
bool_mindelta = ismember(output_struct.Parameter_Grid(:,colind_delta),mindelta);
bool_maxepsilon = ismember(output_struct.Parameter_Grid(:,colind_epsilon),maxepsilon);
bool_plotret = bool_allmids + bool_mindelta + bool_maxepsilon;
bool_plotret = bool_plotret == max(bool_plotret); 
idxs_plotret = find(bool_plotret);

%% 2.3.3 Plotting
htmapplt = 1;
lnplt = 1;
savenclose = 1;
figpathantret = [figpath filesep 'AntRetCompare'];
for i = 1:length(idxs_plotant)
    if lnplt
        TimeCoursePlot(simstr,idxs_plotant(i),'Line',0,loadpath,simpath,savenclose,figpathantret); %#ok<UNRCH> 
        TimeCoursePlot(simstr,idxs_plotant(i),'Line',1,loadpath,simpath,savenclose,figpathantret);
        TimeCoursePlot(simstr,idxs_plotret(i),'Line',0,loadpath,simpath,savenclose,figpathantret);
        TimeCoursePlot(simstr,idxs_plotret(i),'Line',1,loadpath,simpath,savenclose,figpathantret);
    end
    if htmapplt
        TimeCoursePlot(simstr,idxs_plotant(i),'Heatmap',1,loadpath,simpath,savenclose,figpathantret); %#ok<UNRCH> 
        TimeCoursePlot(simstr,idxs_plotret(i),'Heatmap',1,loadpath,simpath,savenclose,figpathantret);
    end
end

%% 2.3.4 Brainframe plots for anterograde
ts = output_struct.Simulations(1).Model_Outputs.Sim.trange;
t1_ind = 1; 
[~,t2_ind] = min(abs(ts - (ts(end)/80)));
[~,t3_ind] = min(abs(ts - (ts(end)/40)));
[~,t4_ind] = min(abs(ts - (ts(end)/20)));
[~,t5_ind] = min(abs(ts - (ts(end)/10)));
[~,t6_ind] = min(abs(ts - (ts(end)/4)));
t7_ind = length(ts);
ts_inds_bf = [t1_ind,t2_ind,t3_ind,t4_ind,t5_ind,t6_ind,t7_ind];

gamma1_antret = 0.004;
lambda_antret = 0.075;
delta_ant = 100; delta_ret = 10;
epsilon_ant = 10; epsilon_ret = 100;

lambdabool = ismember(output_struct.Parameter_Grid(:,colind_lambda),lambda_antret); 
gamma1bool = ismember(output_struct.Parameter_Grid(:,colind_gamma1),gamma1_antret);
deltaantbool = ismember(output_struct.Parameter_Grid(:,colind_delta),delta_ant);
epsilonantbool = ismember(output_struct.Parameter_Grid(:,colind_epsilon),epsilon_ant);
antind = lambdabool + gamma1bool + deltaantbool + epsilonantbool; 
antind = find(antind == 4);

wflow = 1; savenclose = 0;
BrainframePlot(simstr,antind,ts_inds_bf,wflow,savenclose,loadpath,simpath,bfpath,figpath);

%% 2.3.5 Brainframe plots for retrograde
deltaretbool = ismember(output_struct.Parameter_Grid(:,colind_delta),delta_ret);
epsilonretbool = ismember(output_struct.Parameter_Grid(:,colind_epsilon),epsilon_ret);
retind = lambdabool + gamma1bool + deltaretbool + epsilonretbool; 
retind = find(retind == 4);

wflow = 1; savenclose = 0;
BrainframePlot(simstr,retind,ts_inds_bf,wflow,savenclose,loadpath,simpath,bfpath,figpath);

%% 3.1 Plot connectivity to/from seed region
savenclose = 1;
SeedConnectivityPlot(simstr,1,loadpath,simpath,savenclose,figpath);

%% 3.2 Correlation with seed connectivity, anterograde/retrograde
gamma1_antret = 0.004;
lambda_antret = 0.025;
delta_ant = 100; delta_ret = 10;
epsilon_ant = 10; epsilon_ret = 100;

lambdabool = ismember(output_struct.Parameter_Grid(:,colind_lambda),lambda_antret); 
gamma1bool = ismember(output_struct.Parameter_Grid(:,colind_gamma1),gamma1_antret);
deltaantbool = ismember(output_struct.Parameter_Grid(:,colind_delta),delta_ant);
epsilonantbool = ismember(output_struct.Parameter_Grid(:,colind_epsilon),epsilon_ant);
antind = lambdabool + gamma1bool + deltaantbool + epsilonantbool; 
antind = find(antind == 4);

deltaretbool = ismember(output_struct.Parameter_Grid(:,colind_delta),delta_ret);
epsilonretbool = ismember(output_struct.Parameter_Grid(:,colind_epsilon),epsilon_ret);
retind = lambdabool + gamma1bool + deltaretbool + epsilonretbool; 
retind = find(retind == 4);

savenclose = 1;
SimulationSeedCorrelationPlot(simstr,antind,retind,loadpath,simpath,savenclose,figpath)

%% S1. Check mass conservation
MassConservationPlot(simstr,loadpath,simpath,1,figpath);

%% Misc.
% params = output_struct.Parameter_Grid;
% params_ret = params((params(:,7) < params(:,8)),:);
% params_ant = params((params(:,7) > params(:,8)),:);
% retinds = find((params(:,7) < params(:,8)));
% antinds = find((params(:,7) > params(:,8)));
% cmap = lines(size(params_ret,1));
% leg = cell(1,size(params_ret,1));
% figure('Position',[0,0,600,600]); hold on;
% for i = 1:size(params_ret,1)
%     tot_ret = output_struct.Simulations(retinds(i)).Model_Outputs.Predicted.N(inds,:) + ...
%         output_struct.Simulations(retinds(i)).Model_Outputs.Predicted.M(inds,:);
%     tot_ant = output_struct.Simulations(antinds(i)).Model_Outputs.Predicted.N(inds,:) + ...
%         output_struct.Simulations(antinds(i)).Model_Outputs.Predicted.M(inds,:);
%     corrs = corr(tot_ret,tot_ant); corrs = diag(corrs);
%     plot(ts,corrs,'LineWidth',3,'Color',cmap(i,:))
%     leg{i} = ['$\gamma$ = ',num2str(params_ret(i,2),'%.2d')];
% end
% xlabel('Time (Days)'); ylabel("Pearson's R, Ant vs. Ret");
% legend(leg,'Location','northeast','Interpreter','latex');
% set(gca,'FontSize',16,'FontName','Times');


 