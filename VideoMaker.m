function VideoMaker(matstr,basepath)

if nargin < 2
    basepath = cd;
end

fpath = [basepath filesep 'SampleResults'];
outpath = [basepath filesep 'OutputFigures'];
load([fpath filesep matstr '.mat']);

Lt = L1+L2+L_int;
trange = trange/86400;
rmbds = @(x)[x-1,x,x+1];
jn(:,[rmbds(GM1_bound),rmbds(GM2_bound),rmbds(ais_bound),rmbds(syn_bound)]) = []; 
n(:,[rmbds(GM1_bound),rmbds(GM2_bound),rmbds(ais_bound),rmbds(syn_bound)]) = []; 
m(:,[rmbds(GM1_bound),rmbds(GM2_bound),rmbds(ais_bound),rmbds(syn_bound)]) = [];
xmesh_ = xmesh;
xmesh_([rmbds(GM1_bound),rmbds(GM2_bound),rmbds(ais_bound),rmbds(syn_bound)]) = [];
thresh = 3*median(abs(jn(:)));
jn(jn < -thresh) = -thresh;
jn(jn > thresh) = thresh;

v = VideoWriter([outpath filesep matstr '.avi']);
open(v);
figure('Units', 'inch', 'Position', [0 0 18 15]);
subplotinds = 1:8;
figsize = length(subplotinds)*2 + 1;
indcell = cell(1,6);
indmat = reshape(1:(figsize^2),figsize,figsize).';
ind1 = indmat(1:length(subplotinds)-1,1:length(subplotinds)-1).'; indcell{1} = (ind1(:)).';
ind2 = indmat(length(subplotinds)+1,1:length(subplotinds)-1).'; indcell{2} = (ind2(:)).';
ind3 = indmat((2:length(subplotinds))+length(subplotinds)+1,1:length(subplotinds)-1).'; indcell{3} = (ind3(:)).';
ind4 = indmat(1:length(subplotinds),(2:length(subplotinds))+length(subplotinds)).'; indcell{4} = (ind4(:)).';
ind5 = indmat((1:length(subplotinds))+length(subplotinds)+1,(2:length(subplotinds))+length(subplotinds)).'; indcell{5} = (ind5(:)).';

for i = 1:length(trange)
    subplot(figsize,figsize,indcell{1})
    hold off;
    h1 = fill([0,0,L1,L1],[min(n(:)),max(n(:)),max(n(:)),min(n(:))],...
        [211/255,211/255,211/255],'EdgeColor','none');  
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on;    
    h2 = fill([L1+L_int,L1+L_int,L2+L1+L_int,L2+L1+L_int],[min(n(:)),max(n(:)),max(n(:)),min(n(:))],...
        [211/255,211/255,211/255],'EdgeColor','none');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    h3 = fill([L1,L1,L1+L_ais,L1+L_ais],[min(n(:)),max(n(:)),max(n(:)),min(n(:))],...
        [222/255,222/255,222/255],'EdgeColor','none');  
    set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');   
    h4 = fill([L1+L_int-L_syn,L1+L_int-L_syn,L1+L_int,L1+L_int],[min(n(:)),max(n(:)),max(n(:)),min(n(:))],...
        [222/255,222/255,222/255],'EdgeColor','none');
    set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    plot(xmesh_,n(i,:),'Marker','none','Color',[1 0 0],'LineWidth',3)
    ylim([0,max(n(:))]);
    xlim([0,max(xmesh_)]);
    ylabel('Concentration (\muM)','Color',[0.15 0.15 0.15],'FontWeight','bold')
    set(gca,'YColor',[0.15 0.15 0.15])  
    title(sprintf('Soluble Tau @ t = %.2f Days', trange(i)));
    set(gca, 'FontSize', 17);
    
    subplot(figsize,figsize,indcell{2})
    hold off
    jn_ = interp1(xmesh_, jn(i,:), 1:Lt);
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
    ylabel('Tau Flux','FontWeight','bold','FontSize',17);
    
    subplot(figsize,figsize,indcell{3})
    hold off;
    h1 = fill([0,0,L1,L1],[min(m(:)),max(m(:)),max(m(:)),min(m(:))],...
        [211/255,211/255,211/255],'EdgeColor','none');  
    set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    hold on;    
    h2 = fill([L1+L_int,L1+L_int,L2+L1+L_int,L2+L1+L_int],[min(m(:)),max(m(:)),max(m(:)),min(m(:))],...
        [211/255,211/255,211/255],'EdgeColor','none');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    h3 = fill([L1,L1,L1+40,L1+40],[min(m(:)),max(m(:)),max(m(:)),min(m(:))],...
        [222/255,222/255,222/255],'EdgeColor','none');  
    set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    h4 = fill([L1+L_int-40,L1+L_int-40,L1+L_int,L1+L_int],[min(m(:)),max(m(:)),max(m(:)),min(m(:))],...
        [222/255,222/255,222/255],'EdgeColor','none');
    set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    plot(xmesh_,m(i,:),'Marker','none','Color',[0 0 1],'LineWidth',3)
    ylim([0,max(m(:))]);
    xlim([0,max(xmesh_)]);
    xlabel('Axonal Position x (\mum)','FontWeight','bold')
    ylabel('Concentration (\muM)','Color',[0.15 0.15 0.15],'FontWeight','bold')
    set(gca,'YColor',[0.15 0.15 0.15])
    title(sprintf('Insoluble Tau @ t = %.2f Days', trange(i)));
    set(gca, 'FontSize', 17);
    
    subplot(figsize,figsize,indcell{4})
    semilogx(trange,[N1(1:i), NaN(1,(length(trange)-i))],'r-','LineWidth',3); 
    hold on;
    semilogx(trange,[N2(1:i), NaN(1,(length(trange)-i))],'r--','LineWidth',3); 
    ylabel('Mean Concentration (\muM)','FontWeight','bold');
    ylim([0,max([N1,N2])])
    xlim([0,max(trange)])
    title(sprintf('Soluble Somatodendritic Tau'));
    set(gca, 'FontSize', 17)
    legend({'Presynaptic S.D. Soluble', 'Postsynaptic S.D.Soluble'},...
        'FontSize',14,'Location','northwest'); 

    subplot(figsize,figsize,indcell{5})
    semilogx(trange,[M1(1:i), NaN(1,(length(trange)-i))],'b-','LineWidth',3); 
    hold on;
    semilogx(trange,[M2(1:i), NaN(1,(length(trange)-i))],'b--','LineWidth',3); 
    xlabel('t (Days)','FontWeight','bold'); 
    ylabel('Mean Concentration (\muM)','FontWeight','bold');
    ylim([0,max([M1,M2])])
    xlim([0,max(trange)])
    title(sprintf('Insoluble Somatodendritic Tau'));
    set(gca, 'FontSize', 17)
    legend({'Presynaptic S.D. Insoluble', 'Postsynaptic S.D. Insoluble'},...
        'FontSize',14,'Location','northwest'); 

    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v)

end