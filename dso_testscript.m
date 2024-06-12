clear; clc;
figure('Units','inches','Position',[0 0 14 6.8]);
t = tiledlayout(1,2,'TileSpacing','loose');

os_nosr = load('/Users/justintorok/Documents/MATLAB/Tau_Transport/SampleFiles/hippocampome_ret_debug_no_sr.mat','output_struct');
os_withsr = load('/Users/justintorok/Documents/MATLAB/Tau_Transport/SampleFiles/hippocampome_ret_debug_with_sr.mat','output_struct');
% F_nosr = os_nosr.output_struct.Simulations.Model_Outputs.Predicted.F;
N_nosr = os_nosr.output_struct.Simulations.Model_Outputs.Predicted.N;
% TE_nosr = num2str(os_nosr.output_struct.Time_Elapsed/60,'%.0f');
% F_withsr = os_withsr.output_struct.Simulations.Model_Outputs.Predicted.F;
N_withsr = os_withsr.output_struct.Simulations.Model_Outputs.Predicted.N;
% TE_withsr = num2str(os_withsr.output_struct.Time_Elapsed/60,'%.0f');

inds = 1:30; inds(inds == 25) = [];
N_nosr_noseed = N_nosr(inds,:);
N_withsr_noseed = N_withsr(inds,:);

nexttile; hold on;
scatter(N_nosr_noseed,N_withsr_noseed);
plotmax = 1.05*max(N_nosr_noseed(:));
xlim([0,plotmax]); ylim([0,plotmax]);
plot([0 plotmax],[0 plotmax],'LineStyle','--','LineWidth',1.5)
ticklabs = (0:3)*10^-5;
Rval = corr(N_nosr_noseed(:),N_withsr_noseed(:));
text(0.7,0.1,sprintf('R^2 = %.2f',Rval^2),'Units','normalized',...
    'FontSize',22);
xticks(ticklabs); yticks(ticklabs); 
% title('Out-of-Sample (Time)')
set(gca,'FontSize',22,'box','on');

os_nosr = load('/Users/justintorok/Documents/MATLAB/Tau_Transport/SampleFiles/LH_CA1_ret_debug_no_sr.mat','output_struct');
os_withsr = load('/Users/justintorok/Documents/MATLAB/Tau_Transport/SampleFiles/LH_CA1_ret_debug_with_sr.mat','output_struct');
% F_nosr = os_nosr.output_struct.Simulations.Model_Outputs.Predicted.F;
N_nosr = os_nosr.output_struct.Simulations.Model_Outputs.Predicted.N;
% TE_nosr = num2str(os_nosr.output_struct.Time_Elapsed/60,'%.0f');
% F_withsr = os_withsr.output_struct.Simulations.Model_Outputs.Predicted.F;
N_withsr = os_withsr.output_struct.Simulations.Model_Outputs.Predicted.N;
% TE_withsr = num2str(os_withsr.output_struct.Time_Elapsed/60,'%.0f');

inds = 1:213; inds(inds == 30) = [];
N_nosr_noseed = N_nosr(inds,:);
N_withsr_noseed = N_withsr(inds,:);

nexttile; hold on;
scatter(N_nosr_noseed,N_withsr_noseed);
plotmax = 1.05*max(N_nosr_noseed(:));
xlim([0,plotmax]); ylim([0,plotmax]);
plot([0 plotmax],[0 plotmax],'LineStyle','--','LineWidth',1.5)
ticklabs = (0:0.5:2.5)*10^-4;
xticks(ticklabs); yticks(ticklabs);
Rval = corr(N_nosr_noseed(:),N_withsr_noseed(:));
text(0.7,0.1,sprintf('R^2 = %.2f',Rval^2),'Units','normalized',...
    'FontSize',22);
% title('Out-of-Sample (Space)');
set(gca,'FontSize',22,'box','on');
xlabel(t,'N_t_r_u_e','FontSize',28); 
ylabel(t,'N_D_S_O','FontSize',28);
print('DSOfig','-dtiff','-r300'); close;
% nonzeroinds = (F_nosr ~= 0);
% F_nosr_nonzeros = F_nosr(nonzeroinds);
% F_withsr_nonzeros = F_withsr(nonzeroinds);
% figure;
% scatter(F_nosr_nonzeros(:),F_withsr_nonzeros(:));
% % legend({'t0','t1','t2','t3','t4','t5','t6','t7','t8','t9','t10'},'Location','northwest')
% xlabel(['F_N_T_M, Time Elapsed = ' TE_nosr]); 
% ylabel(['F_D_S_O, Time Elapsed = ' TE_withsr]);
% set(gca,'FontSize',24,'FontName','Times');
