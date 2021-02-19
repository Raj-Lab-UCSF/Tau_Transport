function SnapshotPlotter(matstr,customt,savenclose,basepath)

% customt = [1,100,145,188,250]; retrograde default
% customt = [1,100,146,191,234]; mixed default
% customt = [1,100,146,179,234]; anterograde default

if nargin < 4
    basepath = cd;
    if nargin < 3
        savenclose = 0;
    end
end

fpath = [basepath filesep 'SampleFiles'];
outpath = [basepath filesep 'OutputFigures'];
load([fpath filesep matstr '.mat']);

rmbds = @(x)[x-1,x,x+1];
jn(:,[rmbds(GM1_bound),rmbds(GM2_bound),rmbds(ais_bound),rmbds(syn_bound)]) = []; 
n(:,[rmbds(GM1_bound),rmbds(GM2_bound),rmbds(ais_bound),rmbds(syn_bound)]) = []; 
m(:,[rmbds(GM1_bound),rmbds(GM2_bound),rmbds(ais_bound),rmbds(syn_bound)]) = [];
xmesh_ = xmesh;
xmesh_([rmbds(GM1_bound),rmbds(GM2_bound),rmbds(ais_bound),rmbds(syn_bound)]) = [];
thresh = 3*median(abs(jn(:)));
jn(jn < -thresh) = -thresh;
jn(jn > thresh) = thresh;
    
figure('Units', 'inch', 'Position', [0 0 5*length(customt) 5]);
subplotinds = 1:7;
figsize1 = length(subplotinds)*length(customt);
figsize2 = length(subplotinds)+2;
indcell = cell(1,2*length(customt));
indmat = reshape(1:(figsize1*figsize2),figsize1,figsize2).';
Lt = L1 + L_int + L2;
for i = 1:length(indcell)
    if ismember(i,1:length(customt))
        ind_ = indmat(1:(length(subplotinds)-1),(1:(length(subplotinds)-1))+(i-1)*length(subplotinds)).';
        indcell{i} = (ind_(:)).';
    else
        ind_ = indmat(length(subplotinds)+2,(1:(length(subplotinds)-1))+(i-(1+length(customt)))*length(subplotinds)).';
        indcell{i} = (ind_(:)).';
    end
end

indcell = reshape(indcell,length(customt),2).';

for i = 1:length(customt)
    subplot(figsize2,figsize1,indcell{1,i})
    hold off;
    plotmax = 1.2*max([n(:);m(:)]); plotmin = min([n(:);m(:)]); 
    h1 = fill([0,0,L1,L1],[plotmin,plotmax,plotmax,plotmin],...
        [211/255,211/255,211/255],'EdgeColor','none');  
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on;    
    h2 = fill([L1+L_int,L1+L_int,L2+L1+L_int,L2+L1+L_int],[plotmin,plotmax,plotmax,plotmin],...
        [211/255,211/255,211/255],'EdgeColor','none');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    h3 = fill([L1,L1,L1+L_ais,L1+L_ais],[plotmin,plotmax,plotmax,plotmin],...
        [222/255,222/255,222/255],'EdgeColor','none');  
    set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');   
    h4 = fill([L1+L_int-L_syn,L1+L_int-L_syn,L1+L_int,L1+L_int],[plotmin,plotmax,plotmax,plotmin],...
        [222/255,222/255,222/255],'EdgeColor','none');
    set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    plot(xmesh_,n(customt(i),:),'Marker','none','Color',[1 0 0],'LineWidth',3);
    plot(xmesh_,m(customt(i),:),'Marker','none','Color',[0 0 1],'LineWidth',3);
    ylim([0,plotmax]);
    xlim([0,max(xmesh_)]);
    if i == 1
        ylabel('Tau Concentration (\muM)','Color',[0.15 0.15 0.15],'FontWeight','bold')
        set(gca,'YColor',[0.15 0.15 0.15])
        legend({'Soluble','Insoluble'});
    else
        set(gca,'YTickLabels',{});
    end
    if trange(customt(i)) > 0
        title(sprintf('t = %.1f days', trange(customt(i))/86400));
    else
        title(sprintf('t = %.1d days', trange(customt(i))/86400));
    end
    set(gca, 'FontSize', 18);
    
    subplot(figsize2,figsize1,indcell{2,i})
    hold off
    jn_ = interp1(xmesh_, jn(customt(i),:), 1:Lt);
    jn_pos = zeros(1,length(jn_)); jn_neg = jn_pos;
    jn_pos(jn_ > 0) = jn_(jn_ > 0);
    jn_neg(jn_ < 0) = jn_(jn_ < 0);
    ratiomax = max([jn_pos, abs(jn_neg)])/(thresh);
    plotinds = 1:120:Lt;
    q1 = quiver(40*plotinds/mean(plotinds),ones(1,length(plotinds)),jn_pos(plotinds),zeros(1,length(plotinds)),'Color',[1 0 0],'LineWidth',1.5);
    set(q1,'AutoScale','on', 'AutoScaleFactor',0.2*ratiomax)
    hold on;
    q2 = quiver(40*plotinds/mean(plotinds),ones(1,length(plotinds)),jn_neg(plotinds),zeros(1,length(plotinds)),'Color',[1 0 0],'LineWidth',1.5);
    set(q2,'AutoScale','on', 'AutoScaleFactor',0.2*ratiomax)
    yticks([]); xticks([]);
    set(gca,'TickLength',[0 0]);
    ylim([0,2]); xlim([0,40*Lt/mean(plotinds)]);
    if i == 1
        ylabel('Tau Flux','FontWeight','bold');
    end
    set(gca, 'FontSize', 18);
end

if savenclose
    if ~isfolder(outpath)
        mkdir(outpath)
    end
    print([outpath filesep matstr],'-dpng','-r300');
    close;
end
end