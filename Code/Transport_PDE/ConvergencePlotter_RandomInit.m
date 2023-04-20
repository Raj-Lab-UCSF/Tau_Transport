function ConvergencePlotter_RandomInit(savenclose,basepath)
% Takes a few minutes to run given the overhead with calculating pairwise
% relative error

if nargin < 2
    basepath = cd;
    if nargin < 1
        savenclose = 0;
    end
end

fpath = [basepath filesep 'SampleFiles'];
outpath = [basepath filesep 'OutputFigures'];

simstr = 'random_init';
biasstrs = {'_ant','_ret','_mixed'};
figure('Units', 'inch', 'Position', [0 0 30 10]);
% sgtitle('Random Initial Conditions Comparison','FontSize',32,'FontWeight','bold');

for k = 1:length(biasstrs)
    matstr = [simstr biasstrs{k}];
    load([fpath filesep matstr '.mat'],'n_mat','m_mat','numiters','trange');
    pairwise_err_cell = cell(1,numiters-1);
    while_n = 1;
    while_mat = cat(3,n_mat,m_mat);
    while while_n < numiters
        num_rep = size(while_mat,1) - 1;
        nm_i = while_mat(1,:,:);
        nm_i_rep = repmat(nm_i,num_rep,1,1);
        nm_rest = while_mat(2:end,:,:);
        pairwise_err = sum(abs(nm_i_rep - nm_rest),3);
        pairwise_err_cell{while_n} = pairwise_err./sum(nm_i,3);
        while_n = while_n+1;
        while_mat(1,:,:) = [];
    end
    pairwise_diff = [];
    for j = 1:length(pairwise_err_cell)
        pairwise_diff = cat(1,pairwise_diff,pairwise_err_cell{j});
    end

    rng(0);
    randominds = randperm(size(pairwise_diff,1));
    randominds = randominds(1:50);
    trange = trange/86400;
    subplot(1,3,k)
    for j = 1:length(randominds)
        semilogx(trange,pairwise_diff(randominds(j),:),'Color',[0.6,0.6,0.6],'LineStyle','-','LineWidth',0.5);
        hold on;
    end
    p1 = semilogx(trange,mean(pairwise_diff,1),'Color','r','LineWidth',3); hold on
    p2 = semilogx(trange,mean(pairwise_diff,1)+2*std(pairwise_diff),'Color','b','LineWidth',1.5,'LineStyle','--');
    p3 = semilogx(trange,mean(pairwise_diff,1)-2*std(pairwise_diff),'Color','b','LineWidth',1.5,'LineStyle','--');
    xlabel('Time (Days)'); 
    xlim([trange(2),max(trange)]); ylim([0,1.05]);
    xticks([0.01,1,100]); yticks([0,0.5,1]);
    if k == 1
        ylabel('Relative Pairwise Difference'); 
        legend([p1 p2],{'Mean', 'Mean \pm 2 Std'});
    else
        yticklabels({});
    end
    set(gca,'FontSize',24);
    if strcmp(biasstrs{k},'_ant')
        title('Anterograde-Biased Regime','FontSize',28);
    elseif strcmp(biasstrs{k},'_ret')
        title('Retrograde-Biased Regime','FontSize',28);
    elseif strcmp(biasstrs{k},'_mixed')
        title('Net-Unbiased Regime','FontSize',28);
    end
    
    hold off;
    clearvars -except savenclose fpath simstr biasstrs outpath
end
set(findall(gcf,'-property','FontName'),'FontName','Times')

if savenclose
    if ~isfolder(outpath)
        mkdir(outpath)
    end
    print([outpath filesep simstr],'-dtiffn','-r300');
    close;
end
end
