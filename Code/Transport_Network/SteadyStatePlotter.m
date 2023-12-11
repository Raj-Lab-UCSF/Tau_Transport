function SteadyStatePlotter(matstr,basepath,savenclose,figpath)

if nargin < 4
    figpath = cd;
    if nargin < 3
        savenclose = 0;
        if nargin < 2
            basepath = cd;
        end
    end
end
fpath = [basepath filesep 'SampleFiles'];
outpath = figpath;
load([fpath filesep matstr '.mat']); %#ok<LOAD>

rmbds = @(x)[x-1,x,x+1];
n(:,[rmbds(GM1_bound),rmbds(GM2_bound),rmbds(ais_bound),rmbds(syn_bound)]) = []; 
m(:,[rmbds(GM1_bound),rmbds(GM2_bound),rmbds(ais_bound),rmbds(syn_bound)]) = [];
n_ss(:,[rmbds(GM1_bound),rmbds(GM2_bound),rmbds(ais_bound),rmbds(syn_bound)]) = []; 
m_ss(:,[rmbds(GM1_bound),rmbds(GM2_bound),rmbds(ais_bound),rmbds(syn_bound)]) = [];
xmesh_ = xmesh;
xmesh_([rmbds(GM1_bound),rmbds(GM2_bound),rmbds(ais_bound),rmbds(syn_bound)]) = [];
% Lt = L1 + L_int + L2;
    
figure('Units', 'inch', 'Position', [0 0 14 6]);
tiledlayout(1,2,'Padding','tight'); nexttile; hold on;
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
plot(xmesh_,n(end,:),'Marker','none','Color',[1 0 0],'LineWidth',3);
plot(xmesh_,m(end,:),'Marker','none','Color',[0 0 1],'LineWidth',3);
ylim([0,plotmax]);
xlim([0,max(xmesh_)]);

ylabel('Tau Concentration (\muM)','Color',[0.15 0.15 0.15],'FontWeight','bold')
set(gca,'YColor',[0.15 0.15 0.15])
legend({'Soluble','Insoluble'});
title(sprintf('Single-edge simulation'));
set(gca, 'FontSize', 24,'box','off');

nexttile; hold on;
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
plot(xmesh_,n_ss,'Marker','none','Color',[1 0 0],'LineWidth',3);
plot(xmesh_,m_ss,'Marker','none','Color',[0 0 1],'LineWidth',3);
ylim([0,plotmax]);
xlim([0,max(xmesh_)]);
yticks([]);
title(sprintf('Single-edge steady-state'));
set(gca, 'FontSize', 24,'box','off');

set(findall(gcf,'-property','FontName'),'FontName','Times')
if savenclose
    print([outpath filesep matstr],'-dtiffn','-r600');
    close;
end
end