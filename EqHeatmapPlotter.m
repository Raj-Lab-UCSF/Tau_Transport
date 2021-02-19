function EqHeatmapPlotter(paramname,savenclose,basepath)

if nargin < 3
    basepath = cd;
    if nargin < 2
        savenclose = 0;
    end
end

fpath = [basepath filesep 'SampleFiles'];
outpath = [basepath filesep 'OutputFigures'];

switch paramname
    case 'gamma'
        filenames = {'heatmap_plot_default_gamma','heatmap_plot_larger_gamma','heatmap_plot_smaller_gamma'};
        plottitles = {'Original Aggregation Rate','2x Higher Aggregation Rate','2x Lower Aggregation Rate'};
        textpos = [[7.5,15.5];[4.5,15.5];[12.5,15.5]];
    case 'beta'
        filenames = {'heatmap_plot_default_gamma','heatmap_plot_larger_beta','heatmap_plot_smaller_beta'};
        plottitles = {'Original Fragmentation Rate','2x Higher Fragmentation Rate','2x Lower Fragmentation Rate'};
        textpos = [[7.5,15.5];[12.5,15.5];[4.5,15.5]];
    case 'frac'
        filenames = {'heatmap_plot_frac_05','heatmap_plot_frac_075','heatmap_plot_frac_096'};
        plottitles = {'f = 0.50','f = 0.75','f = 0.96'};
        textpos = [[7.5,15.5];[7.5,15.5];[7.5,15.5]];
end

figure('Units', 'inch', 'Position', [0 0 30 11]);
for k = 1:length(filenames)
    subplot(1,3,k); hold off;
    load([fpath filesep filenames{k} '.mat'],'biasmat','deltalist','epsilonlist');
    imagesc(squeeze(biasmat(:,:,end)),[-1,1]); 
    colormap redblue
    set(gca,'DataAspectRatio', [1 1 1]);
    x0s = zeros(1,length(deltalist));
    for i = 1:length(deltalist)
        biasvals = squeeze(biasmat(i,:,end));
        p = polyfit(epsilonlist,biasvals,3);
        rs = roots(p);
        logtest = (rs == real(rs));
        rangetest = (rs <= max(epsilonlist)) + (rs >= min(epsilonlist));
        if sum(logtest) == 1
            x0s(i) = rs(logtest);
        elseif any(rangetest == 2)
            x0s(i) = rs(rangetest == 2);
        else
            mindists = abs(rs - min(epsilonlist));
            maxdists = abs(rs - max(epsilonlist));
            [~, minind] = min([mindists;maxdists]);
            rs_ = [rs;rs];
            x0s(i) = rs_(minind);
        end
    end
    p0 = polyfit(x0s,deltalist,1);
    xvals = 2*(0.5 + x0s*(length(epsilonlist)-1)/(max(epsilonlist)-min(epsilonlist)));
    yvals = 2*((1:length(deltalist)) - 0.5);
    hold on;
    plot(xvals, yvals, 'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 3);
    ylabel('Positive Feedback (\delta)','FontWeight','bold'); 
    xlabel('Negative Feedback (\epsilon)','FontWeight','bold');
    xticks([1 (1+length(epsilonlist))/2 length(epsilonlist)]); 
    yticks([1 (1+length(deltalist))/2 length(deltalist)]); 
    xticklabels([min(deltalist) mean(deltalist) max(deltalist)]);
    yticklabels([min(epsilonlist) mean(epsilonlist) max(epsilonlist)]);
    set(gca,'FontSize',24,'ydir','normal','TickLength',[0 0]);
    text(textpos(k,1),textpos(k,2),['\delta = ' sprintf('%.3f + %.3f',p0(2),p0(1)) '*\epsilon'],'FontSize',22)
    title(plottitles{k}); 
    clear biasmat;
end
tightfig;
if savenclose
    if ~isfolder(outpath)
        mkdir(outpath)
    end
    print([outpath filesep 'heatmap_plot_' paramname],'-dpng','-r300');
    close;
end

figure('Units', 'inch', 'Position', [0 0 1 7]);
imagesc(linspace(-1,1,1001).'); colormap redblue
set(gca,'FontSize',24,'ydir','normal','TickLength',[0 0]);
yticks([1 501 1001]); xticks([]); 
yticklabels([-1 0 1]);

if savenclose
    print([outpath filesep 'heatmap_colorbar_' paramname],'-dpng','-r300');
    close;
end
end