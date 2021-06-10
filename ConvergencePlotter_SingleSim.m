function ConvergencePlotter_SingleSim(savenclose,basepath)
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

simstr = 'constant_n0_L1000';
biasstrs = {'_ant','_ret','_nob'};
colors = {'r','b','m'};
figure('Units', 'inch', 'Position', [0 0 20 10]);
sgtitle('Convergence to Steady State','FontSize',28,'FontWeight','bold');
latepts = 70;
for k = 1:length(biasstrs)
    matstr = [simstr biasstrs{k}];
    load([fpath filesep matstr '.mat'],'n','m','trange');
    diffs_n = abs(n(2:end,:) - n(1:(end-1),:));
    diffs_m = abs(m(2:end,:) - m(1:(end-1),:));
    trange = trange/86400;
    deltat = trange(2:end)-trange(1:(end-1));
    trange_diffs = trange(2:end);
    relerrs = (sum(diffs_n,2) + sum(diffs_m,2))./sum((n(2:end,:)+m(2:end,:)),2);
    relerrs_dt = relerrs.' ./ deltat;
    
    subplot(3,4,[1:3,5:7,9:11]); hold on;
    plot(trange_diffs,relerrs_dt,'Color',colors{k},'LineWidth',3);
end
xlabel('Time (Days)'); 
xlim([trange_diffs(1),max(trange)]); ylim([0,1.05*max(relerrs_dt)]);
xticks([0.01,1,100]); yticks([0,0.5*max(relerrs_dt),max(relerrs_dt)]);
ytickformat('%.1d');
ylabel('\Delta Tau Distribution / \Delta Time'); 
legend({'Anterograde-Biased', 'Retrograde-Biased', 'Net-Unbiased'},'Location','northeast');
set(gca,'FontSize',24,'xscale','log');
hold off;

for k = 1:length(biasstrs)
    matstr = [simstr biasstrs{k}];
    load([fpath filesep matstr '.mat'],'n','m','trange');
    diffs_n = abs(n(2:end,:) - n(1:(end-1),:));
    diffs_m = abs(m(2:end,:) - m(1:(end-1),:));
    trange = trange/86400;
    deltat = trange(2:end)-trange(1:(end-1));
    trange_diffs = trange(2:end);
    relerrs = (sum(diffs_n,2) + sum(diffs_m,2))./sum((n(2:end,:)+m(2:end,:)),2);
    relerrs_dt = relerrs.' ./ deltat;
    
    coord = 4*k;
    subplot(3,4,coord)
    logrelerrs_dt = log(relerrs_dt);
    plot(trange_diffs((end-latepts):end),logrelerrs_dt((end-latepts):end),[colors{k} 'o'],'MarkerSize',10); hold on;
    p = polyfit(trange_diffs((end-latepts):end),logrelerrs_dt((end-latepts):end),1);
    logrelerrs_fit = polyval(p,trange_diffs((end-latepts):end));
    plot(trange_diffs((end-latepts):end),logrelerrs_fit,'k:','LineWidth',3);
    mdl = fitlm(logrelerrs_fit,logrelerrs_dt((end-latepts):end));
    xlabel('Time (Days)'); 
    xlim([trange_diffs(end-latepts),max(trange)]); 
    ylim([1.05*min(logrelerrs_dt((end-latepts):end)),0.95*max(logrelerrs_dt((end-latepts):end))]);
    xticks([trange_diffs(end-latepts),...
            (trange_diffs(end-latepts)+trange_diffs(end))/2,...
            trange_diffs(end)]); 
    yticks([min(logrelerrs_dt((end-latepts):end)),...
            (max(logrelerrs_dt((end-latepts):end))+min(logrelerrs_dt((end-latepts):end)))/2,...
            max(logrelerrs_dt((end-latepts):end))]);
    ytickformat('%.1f'); xtickformat('%.1d');
    ylabel('ln(\Delta Tau / \Delta Time)');     
    if strcmp(biasstrs{k},'_ant')
%         legend({'Ant.-Biased',sprintf('R^2 = %.2f',mdl.Rsquared.Adjusted)},'Location','northeast');
        legend({'Ant.-Biased',sprintf('y = %.3f x + %.1f',p(1),p(2))},'Location','northeast');
    elseif strcmp(biasstrs{k},'_ret')
        legend({'Ret.-Biased',sprintf('y = %.3f x + %.1f',p(1),p(2))},'Location','northeast');
    elseif strcmp(biasstrs{k},'_nob')
        legend({'Net-Unb.',sprintf('y = %.3f x + %.1f',p(1),p(2))},'Location','northeast');
    end    
    text(trange_diffs(end-latepts)*10,...
        min(logrelerrs_dt((end-latepts):end))*0.85,...
        sprintf('R^2 = %.2f',mdl.Rsquared.Adjusted),'FontSize',16);
    set(gca,'FontSize',16);
end

if savenclose
    if ~isfolder(outpath)
        mkdir(outpath)
    end
    print([outpath filesep 'convergence_single'],'-dpng','-r300');
    close;
end



end
