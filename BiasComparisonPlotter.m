function BiasComparisonPlotter(savenclose,basepath)

if nargin < 2
    basepath = cd;
    if nargin < 1
        savenclose = 0;
    end
end

fpath = [basepath filesep 'SampleResults'];
outpath = [basepath filesep 'OutputFigures'];
load([basepath filesep 'Svals_From_DirNTResults.mat'],'svals_adl','tpts_adl',...
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

svals_adl_vec = svals_adl_vec(~isnan(svals_adl_vec));
tpts_adl_vec = tpts_adl_vec(~isnan(svals_adl_vec));
tpts_adl_vec_uni = unique(tpts_adl_vec);
svals_adl_vec_uni = zeros(1,length(tpts_adl_vec_uni));

for i = 1:length(tpts_adl_vec_uni)
    svals_adl_vec_uni(i) = mean(svals_adl_vec(tpts_adl_vec == tpts_adl_vec_uni(i)));
end
svals_adl_vec = svals_adl_vec_uni;
tpts_adl_vec = tpts_adl_vec_uni;

svals_adl_bias = zeros(1,length(svals_adl_vec));
tpts_adl_bias = svals_adl_bias;
for i = 1:length(svals_adl_vec)
    tpt = tpts_adl_vec(i);
    timediff = abs(trange - repmat(tpt,1,length(trange)));
    [~,ind] = min(timediff);
    svals_adl_bias(i) = bias(ind);
    tpts_adl_bias(i) = trange(ind);
end
gof_adl = corr(svals_adl_bias.',svals_adl_vec.');

subplot(1,2,1)
markers = {'o','s','d'};
hold off;
plot(trange,bias,'Color',[1 0 1],'LineStyle',':','LineWidth',4);
hold on;
for i = 1:length(svals_adl)
    scatter(tpts_adl{i},svals_adl{i},250,['k' markers{i}],'filled'); 
end
ylabel('Bias');
xlabel('t (Days)');
ylim([-1,1])
xlim([0,max(trange)])
title('AD-Like Mouse Studies','FontSize',28);
set(gca, 'FontSize', 24)
[~,objh] = legend(lgnd,'Location','northeast','FontSize',24);
objhl = findobj(objh, 'Type', 'patch');
set(objhl, 'MarkerSize', 12);
text(300,-0.8,['Pearson R = ' sprintf('%.2f',gof_adl)],'FontSize',22)

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

svals_nadl_vec = svals_nadl_vec(~isnan(svals_nadl_vec));
tpts_nadl_vec = tpts_nadl_vec(~isnan(svals_nadl_vec));
tpts_nadl_vec_uni = unique(tpts_nadl_vec);
svals_nadl_vec_uni = zeros(1,length(tpts_nadl_vec_uni));

for i = 1:length(tpts_nadl_vec_uni)
    svals_nadl_vec_uni(i) = mean(svals_nadl_vec(tpts_nadl_vec == tpts_nadl_vec_uni(i)));
end
svals_nadl_vec = svals_nadl_vec_uni;
tpts_nadl_vec = tpts_nadl_vec_uni;

svals_nadl_bias = zeros(1,length(svals_nadl_vec));
tpts_nadl_bias = svals_nadl_bias;
for i = 1:length(svals_nadl_vec)
    tpt = tpts_nadl_vec(i);
    timediff = abs(trange - repmat(tpt,1,length(trange)));
    [~,ind] = min(timediff);
    svals_nadl_bias(i) = bias(ind);
    tpts_nadl_bias(i) = trange(ind);
end
gof_nadl = corr(svals_nadl_bias.',svals_nadl_vec.');

subplot(1,2,2);
markers = {'o','p','s','d','^','v','>','<','h'};
hold off;
plot(trange,bias,'Color',[1 0 1],'LineStyle',':','LineWidth',4);
hold on;
for i = 1:length(svals_nadl)
    scatter(tpts_nadl{i},svals_nadl{i},250,['k' markers{i}],'filled'); 
end
ylabel('Bias');
xlabel('t (Days)');
ylim([-1,1])
title('Non AD-Like Mouse Studies','FontSize',28);
set(gca, 'FontSize', 24)
[~,objh] = legend(lgnd,'Location','northeast','FontSize',24);
objhl = findobj(objh, 'Type', 'patch');
set(objhl, 'MarkerSize', 12);
text(300,-0.8,['Pearson R = ' sprintf('%.2f',gof_nadl)],'FontSize',22)

if savenclose
    if ~isfolder(outpath)
        mkdir(outpath)
    end
    print([outpath filesep 'bias_comparison'],'-dpng','-r300');
    close;
end
end