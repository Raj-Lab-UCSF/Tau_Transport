clear; clc;
simstr = 'test_timescale_1';
matdir = '/Users/justintorok/Documents/MATLAB/Tau_Transport';
load([matdir filesep 'MatFiles' filesep 'CCF_labels.mat'],'CCF_labels');
load([matdir filesep 'SampleFiles' filesep simstr '.mat'],'output_struct')
cmap = lines(size(output_struct.Parameter_Grid,1));
trange = output_struct.Simulations(1).Model_Outputs.Sim.trange; 
% dt = output_struct.Simulations(1).Model_Outputs.Sim.dt; 
trange = 6*30 * trange;
rois = CCF_labels([27:37, 78:80, 147, 240:250, 291:293, 360],1);
for i = 1:length(rois)
    if i < 16
        rois{i} = [rois{i} ' RH'];
    else
        rois{i} = [rois{i} ' LH'];
    end
end

for i = 1:length(output_struct.Simulations)
    figure('Position',[0,0,500,1000]); tiledlayout(2,1,'TileSpacing','compact');
    nexttile; hold on;
    N = output_struct.Simulations(i).Model_Outputs.Predicted.N;
    M = output_struct.Simulations(i).Model_Outputs.Predicted.M;
    leghands = [];
    for j = 1:size(N,1)
        plot(trange,N(j,:));
    end
    xlim([0,trange(end)]); ylim([0,max(N(:))]);
    xticks([0,trange(end)/2,trange(end)]); yticks([0,max(N(:))/2,max(N(:))]);
    yticklabels({'0', num2str(max(N(:))/2,'%.1d'),num2str(max(N(:)),'%.1d')});
    ylabel('N');
    txt = {['$\mathbf{\lambda_1}~=~$' num2str(output_struct.Simulations(i).Model_Outputs.Parameters.lambda1)],...
        ['$\mathbf{\lambda_2}~=~$' num2str(output_struct.Simulations(i).Model_Outputs.Parameters.lambda2)],...
        ['$\mathbf{\beta}~=~$' num2str(output_struct.Simulations(i).Model_Outputs.Parameters.beta)],...
        ['$\mathbf{\gamma_1}~=~$' num2str(output_struct.Simulations(i).Model_Outputs.Parameters.gamma1)],...
        ['$\mathbf{\epsilon}~=~$' num2str(output_struct.Simulations(i).Model_Outputs.Parameters.epsilon)],...
        ['$\mathbf{\delta}~=~$' num2str(output_struct.Simulations(i).Model_Outputs.Parameters.delta)]};
    title(['Seed: ' rois{21}]);
    text(0.15,0.775,txt,'FontSize',18,'FontName','Times','Interpreter','latex','Units','normalized');
    set(gca,'FontSize',20,'FontName','Times')
    nexttile; hold on;
    for j = 1:size(M,1)
        h = plot(trange,M(j,:));
        leghands = [leghands h];
    end
%     tot_T = max(N + M,[],2); 
    tot_T = N(:,end) + M(:,end); 
    [~,sortinds] = sort(tot_T,'descend');
    sortinds = sortinds(1:5); leghands = leghands(sortinds);
    legend(leghands,rois(sortinds),'Location','northeast','box','off','FontSize',18)
    xlim([0,trange(end)]); ylim([0,max(M(:))]);
    xticks([0,trange(end)/2,trange(end)]); yticks([0,max(M(:))/2,max(M(:))]);
    yticklabels({'0', num2str(max(M(:))/2,'%.1d'),num2str(max(M(:)),'%.1d')});
    xlabel('t (Days)'); ylabel('M');
    set(gca,'FontSize',20,'FontName','Times')
    print([matdir filesep 'OutputFigures' filesep simstr '_' num2str(i)],'-dpng'); 
%     close;
end
