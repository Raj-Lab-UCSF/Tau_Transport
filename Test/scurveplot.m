load([cd filesep 'SampleFiles' filesep 'beta_gamma_curve_finerange.mat']);
figure('Units','inches','Position',[0 0 10 9]); hold on;
cmap = hsv(size(biasmat,2));
smat = 0.5 - 0.5*biasmat;
legcell = cell(1,length(fraclist));
for i = 1:length(fraclist)
    legcell{i} = sprintf('f = %.1f',fraclist(i));
end
for i = 1:size(biasmat,2)
    gbratio = gammalist(8:24)/beta;
    plot(gbratio,smat(8:24,i,end),'x-','Color',cmap(i,:),'LineWidth',3);
end
xlabel('\gamma/\beta'); ylabel('s'); title('Steady State Bias vs. Aggregation');
legend(legcell,'Location','southeast');
set(gca,'FontSize',20,'xscale','log');