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
%% Load things
figpath = '~/Documents/MATLAB/Tau_Transport/OutputFigures';
simpath = '~/Documents/MATLAB/Tau_Transport/SampleFiles';
loadpath = '~/Documents/MATLAB/Tau_Transport/MatFiles';
bfpath = '~/Documents/MATLAB/Brainframe-Dev/Brainframe';
simstr = 'hippocampome_final_lambda';
simstr_2 = 3;
simstr = [simstr '_' num2str(simstr_2)];
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

%% Check mass conservation
masstots = NaN(27,91);
for i = 1:size(masstots,1)
    sim_ = output_struct.Simulations(i).Model_Outputs.Predicted;
    for j = 1:size(masstots,2)
        Ns = V_inv * sim_.N(:,j);
        Ms = V_inv * sim_.M(:,j);
        edges = C .* sim_.EdgeMass(:,:,j);
        masstots(i,j) = sum(Ns) + sum(Ms) + sum(edges(:));
    end
end
reldiffs = masstots - masstots(:,1);
reldiffs = reldiffs ./ masstots(:,1);
figure; hold on;
for i = 1:size(masstots,1)
    plot(output_struct.Simulations(i).Model_Outputs.Sim.trange*180,reldiffs(i,:)); 
end
xlabel('t (days)'); ylabel('Relative Difference w.r.t. t0'); set(gca,'FontSize',16)

%%
% [Co,Ci] = ConnectomePlot(simstr,1,loadpath,simpath,0,figpath);
% idxs = find(output_struct.Parameter_Grid(:,7)==100);
idxs = [3 12 21]+4;
for i = 1:length(idxs)
%     TimeCoursePlot(simstr,idxs(i),'Heatmap',1,loadpath,simpath,0,figpath);
    TimeCoursePlot(simstr,idxs(i),'Line',0,loadpath,simpath,0,figpath);
end
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

 