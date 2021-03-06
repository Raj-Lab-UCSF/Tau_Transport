function BiasComparisonPlotter(savenclose,basepath)

if nargin < 2
    basepath = cd;
    if nargin < 1
        savenclose = 0;
    end
end

fpath = [basepath filesep 'SampleFiles'];
outpath = [basepath filesep 'OutputFigures'];
load([fpath filesep 'Svals_From_DirNTResults.mat'],'svals_adl','tpts_adl',...
                'adl_names','svals_nadl','tpts_nadl','nadl_names');
figure('Units', 'inch', 'Position', [0 0 30 10]);
% sgtitle('Bias Comparisons','FontSize',32,'FontWeight','bold');

% AD-like studies
matstr = 'bias_comparison_simulation_adl.mat';
load([fpath filesep matstr],'N1','N2','M1','M2','trange','delta','epsilon');
bias = (N2+M2-N1-M1)./(N1+N2+M1+M2+eps);
trange = trange/86400; % convert seconds --> days

lgnd = {['Model Bias - \delta = ' num2str(delta) ', \epsilon = ' num2str(epsilon)]};
svals_adl_vec = []; tpts_adl_vec = [];
for i = 1:length(svals_adl)
    svals_adl{i} = 1 - 2*svals_adl{i};
    tpts_adl{i} = tpts_adl{i}*30;
    svals_adl_vec = [svals_adl_vec, svals_adl{i}];
    tpts_adl_vec = [tpts_adl_vec, tpts_adl{i}];
    lgnd{i+1} = adl_names{i};
end
% lgnd{end+1} = 'Study Mean';

tpts_adl_vec = tpts_adl_vec(~isnan(svals_adl_vec));
svals_adl_vec = svals_adl_vec(~isnan(svals_adl_vec));
% tpts_adl_vec_uni = unique(tpts_adl_vec);
% svals_adl_vec_uni = zeros(1,length(tpts_adl_vec_uni));
% 
% for i = 1:length(tpts_adl_vec_uni)
%     svals_adl_vec_uni(i) = mean(svals_adl_vec(tpts_adl_vec == tpts_adl_vec_uni(i)));
% end
% svals_adl_vec = svals_adl_vec_uni;
% tpts_adl_vec = tpts_adl_vec_uni;

svals_adl_bias = zeros(1,length(svals_adl_vec));
tpts_adl_bias = svals_adl_bias;
for i = 1:length(svals_adl_vec)
    tpt = tpts_adl_vec(i);
    timediff = abs(trange - repmat(tpt,1,length(trange)));
    [~,ind] = min(timediff);
    svals_adl_bias(i) = bias(ind);
    tpts_adl_bias(i) = trange(ind);
end
mdl_adl = fitlm(svals_adl_bias.',svals_adl_vec.');
gof_adl = mdl_adl.Rsquared.Adjusted;

subplot(1,2,1)
markers = {'o','s','d'};
hold off;
plot(trange,bias,'Color',[1 0 1],'LineStyle',':','LineWidth',4);
hold on;
for i = 1:length(svals_adl)
    scatter(tpts_adl{i},svals_adl{i},250,[0.5 0.5 0.5],markers{i},'filled'); 
end
ylabel('Bias');
xlabel('t (Days)');
ylim([-1,1])
xlim([0,max(trange)])
yticks(-1:0.5:1);
title('AD-Like Mouse Studies','FontSize',28);
set(gca, 'FontSize', 24)
[~,objh] = legend(lgnd,'Location','northeast','FontSize',24,'FontName','Times');
objhl = findobj(objh, 'Type', 'patch');
set(objhl, 'MarkerSize', 12);
text(300,-0.7,['R^{2} = ' sprintf('%.2f',gof_adl)],'FontSize',24)

% Non-AD-like studies
matstr = 'bias_comparison_simulation_nadl.mat';
load([fpath filesep matstr],'N1','N2','M1','M2','trange','delta','epsilon');
bias = (N2+M2-N1-M1)./(N1+N2+M1+M2+eps);
trange = trange/86400; % convert seconds --> days
includeinds = [];
excludenames = {'IbaHippInj','IbaStrInj'}; % P301L studies
lgnd = {['Model Bias - \delta = ' num2str(delta) ', \epsilon = ' num2str(epsilon)]};
for i = 1:length(svals_nadl)
    svals = svals_nadl{i};
    svals = svals(~isnan(svals));
    all0 = all(svals == 0);
    all1 = all(svals == 1);
    if ~ismember(nadl_names{i},excludenames) && ~all0 && ~all1 % Remove studies with constant bias
        includeinds = [includeinds,i];
    end
end
svals_nadl = svals_nadl(includeinds);
tpts_nadl = tpts_nadl(includeinds);
nadl_names = nadl_names(includeinds);
svals_nadl_vec = []; tpts_nadl_vec = [];
for i = 1:length(svals_nadl)
    svals_nadl{i} = 1 - 2*svals_nadl{i};
    tpts_nadl{i} = tpts_nadl{i}*30;
    svals_nadl_vec = [svals_nadl_vec, svals_nadl{i}];
    tpts_nadl_vec = [tpts_nadl_vec, tpts_nadl{i}];
    lgnd{i+1} = nadl_names{i};
end
% lgnd{end+1} = 'Study Mean';

tpts_nadl_vec = tpts_nadl_vec(~isnan(svals_nadl_vec));
svals_nadl_vec = svals_nadl_vec(~isnan(svals_nadl_vec));
% tpts_nadl_vec_uni = unique(tpts_nadl_vec);
% svals_nadl_vec_uni = zeros(1,length(tpts_nadl_vec_uni));
% 
% for i = 1:length(tpts_nadl_vec_uni)
%     svals_nadl_vec_uni(i) = mean(svals_nadl_vec(tpts_nadl_vec == tpts_nadl_vec_uni(i)));
% end
% svals_nadl_vec = svals_nadl_vec_uni;
% tpts_nadl_vec = tpts_nadl_vec_uni;

svals_nadl_bias = zeros(1,length(svals_nadl_vec));
tpts_nadl_bias = svals_nadl_bias;
for i = 1:length(svals_nadl_vec)
    tpt = tpts_nadl_vec(i);
    timediff = abs(trange - repmat(tpt,1,length(trange)));
    [~,ind] = min(timediff);
    svals_nadl_bias(i) = bias(ind);
    tpts_nadl_bias(i) = trange(ind);
end
mdl_nadl = fitlm(svals_nadl_bias.',svals_nadl_vec.');
gof_nadl = mdl_nadl.Rsquared.Adjusted;

subplot(1,2,2);
markers = {'o','p','s','d','^','v','>','<','h'};
hold off;
plot(trange,bias,'Color',[1 0 1],'LineStyle',':','LineWidth',4);
hold on;
for i = 1:length(svals_nadl)
    scatter(tpts_nadl{i},svals_nadl{i},250,[0.5 0.5 0.5],markers{i},'filled')
end
ylabel('Bias');
xlabel('t (Days)');
ylim([-1,1])
yticks(-1:0.5:1);
title('Non AD-Like Mouse Studies','FontSize',28);
set(gca, 'FontSize', 24)
[~,objh] = legend(lgnd,'Location','northeast','FontSize',24,'FontName','Times');
objhl = findobj(objh, 'Type', 'patch');
set(objhl, 'MarkerSize', 12);
text(300,-0.7,['R^{2} = ' sprintf('%.2f',gof_nadl)],'FontSize',24)
set(findall(gcf,'-property','FontName'),'FontName','Times')
if savenclose
    if ~isfolder(outpath)
        mkdir(outpath)
    end
    print([outpath filesep 'bias_comparison_vs_t'],'-dtiffn','-r300');
    close;
end

figure('Units','inch','Position',[0 0 12 10]); hold on;
plot(svals_adl_bias,svals_adl_vec,'ro','MarkerFaceColor','r','MarkerSize',10);
plot(svals_nadl_bias,svals_nadl_vec,'bs','MarkerFaceColor','b','MarkerSize',10);
h1 = lsline;
h1(1).LineWidth = 3; h1(1).LineStyle = '--';
h1(2).LineWidth = 3; h1(2).LineStyle = '--';
xlabel('Model Bias'); ylabel('Study Bias');
title('Bias Comparison');
yticks(-1:0.5:1); xticks(-1:0.25:1);
set(gca,'FontSize',28);
legend({['AD-Like Studies - R^{2} = ' sprintf('%.2f',gof_adl)],...
    ['Non-AD-Like Studies - R^{2} = ' sprintf('%.2f',gof_nadl)]},...
    'FontSize',22,'Location','northwest');
set(findall(gcf,'-property','FontName'),'FontName','Times')

if savenclose
    print([outpath filesep 'bias_comparison_gof'],'-dtiffn','-r300');
    close;
end
end