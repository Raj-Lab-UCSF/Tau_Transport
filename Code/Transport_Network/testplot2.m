figpath = '~/Documents/MATLAB/Tau_Transport/OutputFigures';
simpath = '~/Documents/MATLAB/Tau_Transport/SampleFiles';
loadpath = '~/Documents/MATLAB/Tau_Transport/MatFiles';
simstr = 'test_timescale_3c';
load([simpath filesep simstr '.mat'],'output_struct');

[Co,Ci] = ConnectomePlot('test_timescale_3c',6,loadpath,simpath,0,figpath);
TimeCoursePlot(simstr,6,'Heatmap',1,loadpath,simpath,0,figpath);
TimeCoursePlot(simstr,7,'Heatmap',1,loadpath,simpath,0,figpath);

inds = [1:24 26:30];
ts = output_struct.Simulations(2).Model_Outputs.Sim.trange * 180;
tot_ret = output_struct.Simulations(2).Model_Outputs.Predicted.N(inds,:) + ...
    output_struct.Simulations(2).Model_Outputs.Predicted.M(inds,:);
tot_ant = output_struct.Simulations(3).Model_Outputs.Predicted.N(inds,:) + ...
    output_struct.Simulations(3).Model_Outputs.Predicted.M(inds,:);
corrs_ret_Co = corr(Co,tot_ret);
corrs_ret_Ci = corr(Ci,tot_ret);
corrs_ant_Co = corr(Co,tot_ant);
corrs_ant_Ci = corr(Ci,tot_ant);


figure('Position',[0,0,1000,600]); 
subplot(1,2,1); hold on;
plot(ts,corrs_ant_Ci,'LineWidth',3,'Color','red');
plot(ts,corrs_ant_Co,'LineWidth',3,'Color','blue');
xlabel('Time Days'); ylabel("Pearson's R with C");
legend({'C_i_n, seed','C_o_u_t, seed'},'Location','northeast');
title('Ant, small gamma Simulation');
set(gca,'FontSize',16,'FontName','Times')
subplot(1,2,2); hold on;
plot(ts,corrs_ret_Ci,'LineWidth',3,'Color','red');
plot(ts,corrs_ret_Co,'LineWidth',3,'Color','blue');
xlabel('Time Days'); ylabel("Pearson's R with C");
legend({'C_i_n, seed','C_o_u_t, seed'},'Location','northeast');
title('Ret small gamma Simulation');
set(gca,'FontSize',16,'FontName','Times')