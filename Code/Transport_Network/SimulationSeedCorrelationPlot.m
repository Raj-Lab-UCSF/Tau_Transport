function SimulationSeedCorrelationPlot(simstr,idxant,idxret,loadpath,simpath,savenclose,figpath)

if nargin < 7
    figpath = '~/Documents/MATLAB/Tau_Transport/OutputFigures';
    if nargin < 6
        savenclose = 0;
        if nargin < 5
            simpath = '~/Documents/MATLAB/Tau_Transport/SampleFiles';
            if nargin < 4
                loadpath = '~/Documents/MATLAB/Tau_Transport/MatFiles';
            end
        end
    end
end

load([simpath filesep simstr '.mat'],'output_struct');
load([loadpath filesep 'DefaultAtlas.mat'],'DefaultAtlas');
V_inv = 1./DefaultAtlas.volumes; V_inv = diag(V_inv);
load([loadpath filesep 'CCF_labels.mat'],'CCF_labels');
switch output_struct.Simulations(idxant).Model_Outputs.Sim.connectome_subset
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
    otherwise
        inds = logical(ones(size(Conn,1),1)); %#ok<LOGL> 
end
if ~isfield(output_struct.Simulations(idxant).Model_Outputs.Sim,'C')   
    C = readmatrix([loadpath filesep 'mouse_connectome_19_01.csv']);
    C = C(inds,inds);
else
    C = output_struct.Simulations(idxant).Model_Outputs.Sim.C;
end
V_inv = V_inv(inds,inds);
C = V_inv * C;
% 
% if ~isfield(output_struct.Simulations(idxant).Model_Outputs.Sim,'region_names')
%     load([loadpath filesep 'CCF_labels.mat'],'CCF_labels');
%     switch output_struct.Simulations(idxant).Model_Outputs.Sim.connectome_subset
%         case 'Hippocampus'
%             inds = ismember(CCF_labels(:,3),'Hippocampus');
%         case 'Hippocampus+PC+RSP'
%             inds_hipp = ismember(CCF_labels(:,3),'Hippocampus');
%             inds_pc = ismember(CCF_labels(:,1),'Piriform area');
%             inds_rsp = ismember(CCF_labels(:,3),'Retrosplenial Area');
%             inds = logical(inds_hipp + inds_pc + inds_rsp);
%         case 'RH'
%             inds = ismember(CCF_labels(:,4),'Right Hemisphere');
%         case 'LH'
%             inds = ismember(CCF_labels(:,4),'Left Hemisphere');
%         otherwise
%             inds = logical(ones(size(Conn,1),1)); %#ok<LOGL> 
%     end
%     regnamecell = CCF_labels(inds,:);
%     regnames = cell(size(regnamecell,1),1);
%     for i = 1:length(regnames)
%         regname = regnamecell{i,1};
%         reghem = regnamecell{i,4};
%         if strcmp(reghem,'Right Hemisphere')
%             regnames{i} = [regname ' RH'];
%         else
%             regnames{i} = [regname ' LH'];
%         end
%     end
% else
%     regnames = output_struct.Simulations(idxant).Model_Outputs.Sim.region_names;
% end
initpath = output_struct.Simulations(idxant).Model_Outputs.Sim.init_path;

C(logical(eye(size(C)))) = 0; C = C/max(C(:));
C_out = C(initpath > 0,:).';
C_out(initpath > 0) = [];
outmax = max(C_out);
C_in = C(:,initpath > 0);
C_in(initpath > 0) = [];
inmax = max(C_in);
C_out = C_out/max([outmax,inmax]);
C_in = C_in/max([outmax,inmax]);

ts = output_struct.Simulations(idxant).Model_Outputs.Sim.trange * 180;
N_ant = output_struct.Simulations(idxant).Model_Outputs.Predicted.N;
M_ant = output_struct.Simulations(idxant).Model_Outputs.Predicted.M;
X_ant = N_ant + M_ant; X_ant(initpath > 0,:) = [];
N_ret = output_struct.Simulations(idxret).Model_Outputs.Predicted.N;
M_ret = output_struct.Simulations(idxret).Model_Outputs.Predicted.M;
X_ret = N_ret + M_ret; X_ret(initpath > 0,:) = [];

corrs_ant_Co = corr(C_out,X_ant);
corrs_ant_Ci = corr(C_in,X_ant);
corrs_ret_Co = corr(C_out,X_ret);
corrs_ret_Ci = corr(C_in,X_ret);

figure('Position',[0,0,900,450]); 
tiledlayout(1,2,"TileSpacing","compact");
nexttile; hold on;
plot(ts,corrs_ant_Ci,'LineWidth',3,'Color','red');
plot(ts,corrs_ant_Co,'LineWidth',3,'Color','blue');
xlabel('t (Days)'); ylabel("Pearson's R");
xmin = 0; xmax = ts(end); xlim([xmin, xmax]);
xticks([xmin, (xmin+xmax)/2, xmax]);
yticklabels({num2str(xmin), num2str((xmin+xmax)/2), num2str(xmax)})
ymin = 0.9*min([min(corrs_ant_Ci) min(corrs_ret_Ci) min(corrs_ant_Co) min(corrs_ret_Co)]);
ymax = 1;
ylim([ymin, ymax]); yticks([ymin, (ymin+ymax)/2, ymax]);
yticklabels({num2str(ymin,'%.2f'), num2str((ymin+ymax)/2,'%.2f'), num2str(ymax)})
legend({'C_i_n, LH ECl','C_o_u_t, LH ECl'},'Location','northeast');
title('Anterograde Condition');
set(gca,'FontSize',24,'FontName','Times')

nexttile; hold on;
plot(ts,corrs_ret_Ci,'LineWidth',3,'Color','red');
plot(ts,corrs_ret_Co,'LineWidth',3,'Color','blue');
xlabel('t (Days)'); 
xticks([xmin, (xmin+xmax)/2, xmax]);
yticklabels({num2str(xmin), num2str((xmin+xmax)/2), num2str(xmax)})
ylim([ymin, ymax]); yticks([ymin, (ymin+ymax)/2, ymax]);
yticklabels([])
title('Retrograde Condition');
set(gca,'FontSize',24,'FontName','Times')

if savenclose
    figstr = [simstr '_' 'sims' num2str(idxant) '_' num2str(idxret) '_' 'SeedConnCorrs'];
    print([figpath filesep figstr],'-dtiffn','-r600');
    close;
end
end