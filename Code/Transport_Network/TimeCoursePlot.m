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

trange = 6*30*output_struct.Simulations(idx).Model_Outputs.Sim.trange;
N = output_struct.Simulations(idx).Model_Outputs.Predicted.N;
M = output_struct.Simulations(idx).Model_Outputs.Predicted.M;
cmap_line = prism(size(N,1));
if exclseed
    N(initpath > 0,:) = []; 
    M(initpath > 0,:) = [];
    regnames(initpath > 0) = [];
    cmap_line(initpath > 0,:) = [];
end
X = N + M;

if strcmp(plottype,'Line')
    if exclseed
    end
    figure('Position',[0,0,600,500]); hold on;
%     tiledlayout(2,1,'TileSpacing','compact');
%     nexttile; hold on;
%     leghands = [];
%     for j = 1:size(N,1)
%         plot(trange,N(j,:));
%     end
    for j = 1:size(X,1)
        plot(trange,X(j,:),'LineWidth',2,'Color',cmap_line(j,:));
    end
    xlim([0,trange(end)]); ylim([0,max(X(:))]);
    xticks([0,trange(end)/2,trange(end)]); yticks([0,max(X(:))/2,max(X(:))]);
    yticklabels({'0', num2str(max(X(:))/2,'%.1d'),num2str(max(X(:)),'%.1d')});
    ylabel('[Tau] ($\mu$M)','interpreter','latex'); xlabel('t (Days)');
    set(gca,'FontName','Times','FontSize',28);
%     title(['Seed: ' seedreg]);
%     nexttile; hold on;
%     for j = 1:size(M,1)
%         h = plot(trange,M(j,:));
%         leghands = [leghands h];
%     end
%     tot_T = N(:,end) + M(:,end); 
    if ~exclseed
        txt1 = {['$\mathbf{\lambda_1}~=~$' num2str(output_struct.Simulations(idx).Model_Outputs.Parameters.lambda1)],...
            ['$\mathbf{\beta}~=~$' num2str(output_struct.Simulations(idx).Model_Outputs.Parameters.beta)],...
            ['$\mathbf{\delta}~=~$' num2str(output_struct.Simulations(idx).Model_Outputs.Parameters.delta)]};
        text(0.3,0.85,txt1,'FontSize',24,'FontName','Times','Interpreter','latex','Units','normalized');
        txt2 = {['$\mathbf{\lambda_2}~=~$' num2str(output_struct.Simulations(idx).Model_Outputs.Parameters.lambda2)],...
            ['$\mathbf{\gamma_1}~=~$' num2str(output_struct.Simulations(idx).Model_Outputs.Parameters.gamma1)],...
            ['$\mathbf{\epsilon}~=~$' num2str(output_struct.Simulations(idx).Model_Outputs.Parameters.epsilon)]};
        text(0.65,0.85,txt2,'FontSize',24,'FontName','Times','Interpreter','latex','Units','normalized');
    end
%     xlim([0,trange(end)]); ylim([0,max(M(:))]);
%     xticks([0,trange(end)/2,trange(end)]); yticks([0,max(M(:))/2,max(M(:))]);
%     yticklabels({'0', num2str(max(M(:))/2,'%.1d'),num2str(max(M(:)),'%.1d')});
%     xlabel('t (Days)'); ylabel('M');
elseif strcmp(plottype,'Heatmap')
    trange_int = linspace(0,180,181);
    N_int = zeros(size(N,1),length(trange_int));
    M_int = N_int;
    for i = 1:size(N,1)
        N_int(i,:) = interpn(trange,N(i,:),trange_int);
        M_int(i,:) = interpn(trange,M(i,:),trange_int);
    end
    X_int = N_int + M_int;

    figure('Position',[0,0,1200,1200]); 
%     t = tiledlayout(1,2); nexttile;
    maxX = 0.0025;
    cmap_ = [[ones(650,1), linspace(1,0.5,650).', linspace(1,0,650).'];...
            [ones(350,1), linspace(0.5,0,350).', 0*ones(350,1)]];
%     cmap_ = twocolor(cmap_(1,:),cmap_(2,:),1000);
    imagesc(X_int,[0,maxX]); colormap(cmap_); colorbar;
%     xlabs = cell(1,length(subclasses_));
%     for i = 1:length(subclasses_)
%         col = cmap_x(i,:);
%         xlabs{i} = sprintf('\\color[rgb]{%f,%f,%f}%s',col(1),col(2),col(3),subclasses_{i});
%     end
    trangecell = cell(1,length(trange_int));
    trangeinds = [1, (max(trange_int)/2)+1, max(trange_int)+1];
    for i = 1:length(trange_int)
        if ismember(i,trangeinds)
            trangecell{i} = num2str(trange_int(i),'%.0f');
        end
    end
    set(gca,'TickLength',[0 0],'XTick',1:length(trange_int),...
        'XTickLabel',trangecell,'XTickLabelRotation',0,...
        'YTick',1:length(regnames),'YTickLabel',regnames,...
        'TickLabelInterpreter','tex','FontName','Times','FontSize',26)
    title('[Tau] ($\mu$M)','interpreter','latex','FontSize',28); xlabel('t (Days)','FontSize',28);

%     nexttile;
%     imagesc(M_int); colormap('hot'); colorbar;
%     set(gca,'TickLength',[0 0],'XTick',1:length(trange_int),...
%         'XTickLabel',trangecell,'XTickLabelRotation',0,'YTickLabel',{},...
%         'TickLabelInterpreter','tex','FontName','Times','FontSize',20)
%     title('M');
%     xlabel(t,'Time (Days)','FontName','Times','FontSize',20,'FontWeight','bold')
end

if savenclose
    figstr = [simstr '_' 'sim' num2str(idx) '_' plottype];
    if exclseed
        figstr = [figstr '_' 'noseed'];
    else
        figstr = [figstr '_' 'withseed'];
    end
    print([figpath filesep figstr],'-dtiffn','-r600');
%     exportgraphics(gcf,[figpath filesep figstr '.tif'],'Resolution',600);
%     saveas(gcf,[figpath filesep figstr],'pdf');
    close;
end
end