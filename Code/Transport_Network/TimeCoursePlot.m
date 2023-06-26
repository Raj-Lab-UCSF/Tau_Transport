function TimeCoursePlot(simstr,idx,plottype,exclseed,loadpath,simpath,savenclose,figpath)

if nargin < 8
    figpath = '~/Documents/MATLAB/Tau_Transport/OutputFigures';
    if nargin < 7
        savenclose = 0;
        if nargin < 6
            simpath = '~/Documents/MATLAB/Tau_Transport/SampleFiles';
            if nargin < 5
                loadpath = '~/Documents/MATLAB/Tau_Transport/MatFiles';
            end
        end
    end
end

load([simpath filesep simstr '.mat'],'output_struct');
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
if ~isfield(output_struct.Simulations(idx).Model_Outputs.Sim,'seed')
    seedreg = regnames(initpath > 0);
else
    seedreg = output_struct.Simulations(idx).Model_Outputs.Sim.seed;
end
if length(seedreg) > 1
    s2 = cell(1,length(seedreg)); s2(1:end-1) = {', '}; s2(end) = {''};
    seedreg = [seedreg; s2]; seedreg = [seedreg{:}];
else
    seedreg = char(seedreg);
end

trange = 6*30*output_struct.Simulations(idx).Model_Outputs.Sim.trange;
N = output_struct.Simulations(idx).Model_Outputs.Predicted.N;
M = output_struct.Simulations(idx).Model_Outputs.Predicted.M;
if exclseed
    N(initpath > 0,:) = []; 
    M(initpath > 0,:) = [];
    regnames(initpath > 0) = [];
end

if strcmp(plottype,'Line')
    figure('Position',[0,0,500,1000]); tiledlayout(2,1,'TileSpacing','compact');
    nexttile; hold on;
    leghands = [];
    for j = 1:size(N,1)
        plot(trange,N(j,:));
    end
    xlim([0,trange(end)]); ylim([0,max(N(:))]);
    xticks([0,trange(end)/2,trange(end)]); yticks([0,max(N(:))/2,max(N(:))]);
    yticklabels({'0', num2str(max(N(:))/2,'%.1d'),num2str(max(N(:)),'%.1d')});
    ylabel('N');
    txt = {['$\mathbf{\lambda_1}~=~$' num2str(output_struct.Simulations(idx).Model_Outputs.Parameters.lambda1)],...
        ['$\mathbf{\lambda_2}~=~$' num2str(output_struct.Simulations(idx).Model_Outputs.Parameters.lambda2)],...
        ['$\mathbf{\beta}~=~$' num2str(output_struct.Simulations(idx).Model_Outputs.Parameters.beta)],...
        ['$\mathbf{\gamma_1}~=~$' num2str(output_struct.Simulations(idx).Model_Outputs.Parameters.gamma1)],...
        ['$\mathbf{\epsilon}~=~$' num2str(output_struct.Simulations(idx).Model_Outputs.Parameters.epsilon)],...
        ['$\mathbf{\delta}~=~$' num2str(output_struct.Simulations(idx).Model_Outputs.Parameters.delta)]};
    title(['Seed: ' seedreg]);
    text(0.6,0.775,txt,'FontSize',18,'FontName','Times','Interpreter','latex','Units','normalized');
    set(gca,'FontSize',20,'FontName','Times')
    nexttile; hold on;
    for j = 1:size(M,1)
        h = plot(trange,M(j,:));
        leghands = [leghands h];
    end
    tot_T = N(:,end) + M(:,end); 
    [~,sortinds] = sort(tot_T,'descend');
    sortinds = sortinds(1:5); leghands = leghands(sortinds);
    legend(leghands,regnames(sortinds),'Location','northeast','box','off','FontSize',18)
    xlim([0,trange(end)]); ylim([0,max(M(:))]);
    xticks([0,trange(end)/2,trange(end)]); yticks([0,max(M(:))/2,max(M(:))]);
    yticklabels({'0', num2str(max(M(:))/2,'%.1d'),num2str(max(M(:)),'%.1d')});
    xlabel('t (Days)'); ylabel('M');
    set(gca,'FontSize',20,'FontName','Times')
elseif strcmp(plottype,'Heatmap')
    figure('Position',[0,0,3000,800]); 
    t = tiledlayout(1,2); nexttile;
    imagesc(N); colormap('hot'); colorbar;
%     xlabs = cell(1,length(subclasses_));
%     for i = 1:length(subclasses_)
%         col = cmap_x(i,:);
%         xlabs{i} = sprintf('\\color[rgb]{%f,%f,%f}%s',col(1),col(2),col(3),subclasses_{i});
%     end
    trangecell = cell(1,length(trange));
    for i = 1:length(trange)
        trangecell{i} = num2str(trange(i),'%.1f');
    end
    set(gca,'TickLength',[0 0],'XTick',1:length(trange),...
        'XTickLabel',trangecell,...
        'YTick',1:length(regnames),'YTickLabel',regnames,...
        'TickLabelInterpreter','tex','FontName','Times','FontSize',20)
    title('N');

    nexttile;
    imagesc(M); colormap('hot'); colorbar;
    set(gca,'TickLength',[0 0],'XTick',1:length(trange),...
        'XTickLabel',trangecell,'YTickLabel',{},...
        'TickLabelInterpreter','tex','FontName','Times','FontSize',20)
    title('M');
    xlabel(t,'Time (Days)','FontName','Times','FontSize',20,'FontWeight','bold')
end

if savenclose
    figstr = [simstr '_' 'sim' num2str(idx) '_' plottype];
    print([figpath filesep figstr],'-dtiffn','-r600');
    close;
end
end