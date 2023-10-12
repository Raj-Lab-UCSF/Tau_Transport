function [C_out,C_in] = SeedConnectivityPlot(simstr,idx,loadpath,simpath,savenclose,figpath)

if nargin < 6
    figpath = '~/Documents/MATLAB/Tau_Transport/OutputFigures';
    if nargin < 5
        savenclose = 0;
        if nargin < 4
            simpath = '~/Documents/MATLAB/Tau_Transport/SampleFiles';
            if nargin < 3
                loadpath = '~/Documents/MATLAB/Tau_Transport/MatFiles';
            end
        end
    end
end

load([simpath filesep simstr '.mat'],'output_struct');
load([loadpath filesep 'DefaultAtlas.mat'],'DefaultAtlas');
V_inv = 1./DefaultAtlas.volumes; V_inv = diag(V_inv);
load([loadpath filesep 'CCF_labels.mat'],'CCF_labels');
switch output_struct.Simulations(idx).Model_Outputs.Sim.connectome_subset
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
if ~isfield(output_struct.Simulations(idx).Model_Outputs.Sim,'C')   
    C = readmatrix([loadpath filesep 'mouse_connectome_19_01.csv']);
    C = C(inds,inds);
else
    C = output_struct.Simulations(idx).Model_Outputs.Sim.C;
end
V_inv = V_inv(inds,inds);
C = V_inv * C;

if ~isfield(output_struct.Simulations(idx).Model_Outputs.Sim,'region_names')
    load([loadpath filesep 'CCF_labels.mat'],'CCF_labels');
    switch output_struct.Simulations(idx).Model_Outputs.Sim.connectome_subset
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
    regnamecell = CCF_labels(inds,:);
    regnames = cell(size(regnamecell,1),1);
    for i = 1:length(regnames)
        regname = regnamecell{i,1};
        reghem = regnamecell{i,4};
        if strcmp(reghem,'Right Hemisphere')
            regnames{i} = [regname ' RH'];
        else
            regnames{i} = [regname ' LH'];
        end
    end
else
    regnames = output_struct.Simulations(idx).Model_Outputs.Sim.region_names;
end
initpath = output_struct.Simulations(idx).Model_Outputs.Sim.init_path;
% if ~isfield(output_struct.Simulations(idx).Model_Outputs.Sim,'seed')
%     seedreg = regnames(initpath > 0);
% else
%     seedreg = output_struct.Simulations(idx).Model_Outputs.Sim.seed;
% end
% if length(seedreg) > 1
%     s2 = cell(1,length(seedreg)); s2(1:end-1) = {', '}; s2(end) = {''};
%     seedreg = [seedreg; s2]; seedreg = [seedreg{:}];
% else
%     seedreg = char(seedreg);
% end

C(logical(eye(size(C)))) = 0; C = C/max(C(:));
C_out = C(initpath > 0,:).';
C_out(initpath > 0) = [];
outmax = max(C_out);
C_in = C(:,initpath > 0);
C_in(initpath > 0) = [];
inmax = max(C_in);
C_out = C_out/max([outmax,inmax]);
C_in = C_in/max([outmax,inmax]);
regnames(initpath > 0) = [];

figure('Position',[0,0,600,800]); 
tiledlayout(1,2); nexttile;
imagesc(C_out, [0 1]); colormap('hot');
set(gca,'TickLength',[0 0],'XTickLabel',{},...
    'YTick',1:length(regnames),'YTickLabel',regnames,...
    'TickLabelInterpreter','tex','FontName','Times','FontSize',20)
title('Out');

nexttile;
imagesc(C_in,[0 1]); colormap('hot'); colorbar;
set(gca,'TickLength',[0 0],'XTickLabel',{},'YTickLabel',{},...
    'TickLabelInterpreter','tex','FontName','Times','FontSize',20)
title('In');

if savenclose
    figstr = [simstr '_' 'sim' num2str(idx) '_' 'Connectivity'];
    print([figpath filesep figstr],'-dtiffn','-r600');
    close;
end
end