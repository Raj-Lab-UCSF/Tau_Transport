function AMBCAHeatmap_transport(simstr,idx,loadpath,simpath,savenclose_,figdir_)

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
thresh = 80;
threshval = prctile(nonzeros(C(:)),thresh);
C_norm = C; C_norm(C_norm > threshval) = threshval;
C_norm(logical(eye(size(C_norm)))) = 0;
C_norm = (C_norm - min(C_norm(:)))/(max(C_norm(:)) - min(C_norm(:)));

cmap_ = [[ones(650,1), linspace(1,0.5,650).', linspace(1,0,650).'];...
        [ones(350,1), linspace(0.5,0,350).', 0*ones(350,1)]];
figure('Units','inches','Position',[0 0 18 18]);
imagesc(C_norm); colormap(cmap_); clrbr = colorbar; axis square;
set(clrbr,'YTick',0:0.25:1);
set(gca,'TickLength',[0 0],...
    'YTick',1:length(regnames),'YTickLabel',regnames,...
    'XTick',1:length(regnames),'XTickLabel',regnames,...
    'TickLabelInterpreter','tex','FontName','Times','FontSize',18)

if savenclose_
    print([figdir_ filesep 'AMBCA_Heatmap'],'-dtiffn','-r300'); close;
end
end