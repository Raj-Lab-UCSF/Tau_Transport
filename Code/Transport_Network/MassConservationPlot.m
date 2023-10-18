function MassConservationPlot(simstr,loadpath,simpath,savenclose,figpath)
if nargin < 5 
    figpath = cd;
    if nargin < 4
        savenclose = 0;
        if nargin < 3
            simpath = cd;
            if nargin < 2
                loadpath = cd;
            end
        end
    end
end

load([loadpath filesep 'DefaultAtlas.mat'],'DefaultAtlas');
load([simpath filesep simstr '.mat'],'output_struct');
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

ts = output_struct.Simulations(1).Model_Outputs.Sim.trange;
masstots = NaN(size(output_struct.Parameter_Grid,1),length(ts));
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

plottype = 'MassCons';
cmap = hsv(size(output_struct.Parameter_Grid,1));
figure; hold on;
for i = 1:size(masstots,1)
    plot(output_struct.Simulations(i).Model_Outputs.Sim.trange*180,reldiffs(i,:),'Color',cmap(i,:)); 
end
xlabel('t (Days)'); ylabel('Relative Difference w.r.t. t0'); xlim([0,max(ts)*180])
set(gca,'FontSize',20,'FontName','Times')

if savenclose
    figstr = [simstr '_' plottype];
    print([figpath filesep figstr],'-dtiffn','-r600');
    close;
end