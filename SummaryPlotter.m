function SummaryPlotter(savenclose,basepath)

if nargin < 2
    basepath = cd;
    if nargin < 1
        savenclose = 0;
    end
end

fpath = [basepath filesep 'SampleFiles'];
outpath = [basepath filesep 'OutputFigures'];

filestr = 'constant_n0_L1000_coarse_'; % coarse xmesh gives smoother bias manifolds at early t
matstrs = {'ant','ret','nob'};
titles = {'Anterograde-Biased Regime','Retrograde-Biased Regime','Net Unbiased Regime'};
figure('Units', 'inch', 'Position', [0 0 28 8]);
for i = 1:length(matstrs)
    load([fpath filesep filestr matstrs{i} '.mat'],'N1','N2','M1','M2','trange');
    bias = (N2+M2-N1-M1)./(N1+N2+M1+M2+eps);
    trange = trange/86400; % convert seconds --> days

    subplot(1,3,i)
    yyaxis left
    hold off;
    semilogx(trange,N1,'Color',[1 0 0],'LineStyle','-','LineWidth',3);
    hold on;
    semilogx(trange,N2,'Color',[1 0 0],'LineStyle','-.','LineWidth',3); 
    semilogx(trange,M1,'Color',[0 0 1],'LineStyle','-','LineWidth',3);
    semilogx(trange,M2,'Color',[0 0 1],'LineStyle','-.','LineWidth',3);
    if i == 1
        ylabel('Mean Somatodendritic Tau (\muM)');
    end
    ylim([0,0.25])
    yticks([0,0.05,0.1,0.15,0.2]);
%     if i == 1
%         yticklabels({'0','0.05','0.1','0.15','0.2'});
%     else
%         yticklabels({''})
%     end
    set(gca,'ycolor','k');
    yyaxis right
    semilogx(trange,bias,'Color',[1 0 1],'LineStyle',':','LineWidth',3);
    xlabel('t (Days)');
    xticks([0.01,1,100])
    xlim([0,max(trange)])
    ylim([-0.85,0.85]);
    yticks([-0.75,0,0.75]);
    if i == 3
        ylabel('Bias Between S.D. Compartments')
    end
    title(titles{i});
    set(gca, 'FontSize', 24, 'ycolor', 'k')
    if i == 1
        legend({'Presynaptic Soluble','Postsynaptic Soluble',...
            'Presynaptic Insoluble','Postsynaptic Insoluble','Net Bias'},...
            'Location','northwest','FontSize',18)
    end
    clear N1 N2 M1 M2 trange bias
end
tightfig;
if savenclose
    if ~isfolder(outpath)
        mkdir(outpath)
    end
    print([outpath filesep 'summary_plots'],'-dpng','-r300');
    close;
end
end