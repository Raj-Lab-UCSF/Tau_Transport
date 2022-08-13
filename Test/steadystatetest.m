matdir = '/Users/justintorok/Documents/MATLAB/Tau_Transport/SampleFiles';
biascalc = @(n,m) (n(end) + m(end) - n(1) - m(1))/...
    (n(end) + n(1) + m(end) + m(1));
% % % 
simstr = 'constant_n0_L1000_fine_ant';
load([matdir filesep simstr],'n','m','xmesh');
[nss,mss,xmeshss,Bss] = SteadyStateCalculator('resmesh','coarse','n0',...
    [zeros(1,40), 20*ones(1,920), zeros(1,40)]);

figure('Units','inches','Position',[0 0 15 7]); 
subplot(1,2,1); hold on;
plot(xmesh,n(end,:),'r','LineWidth',2); 
plot(xmesh,m(end,:),'b','LineWidth',2); 
ylim([0,max([mss,m(end,:)])]);
legend({'Steady State n', 'Steady State m'},'Location','northwest'); 
xlabel('x'); 
ylabel('Concentration'); 
title('t_f for PDE'); 
set(gca,'FontSize',16)
subplot(1,2,2); hold on;
plot(xmeshss,nss,'r','LineWidth',2); 
plot(xmeshss,mss,'b','LineWidth',2); 
ylim([0,max([mss,m(end,:)])]);
xlabel('x'); 
title(sprintf('Steady-state for B = %0.3f',Bss)); 
set(gca,'FontSize',16)
sgtitle('\delta = 1, \epsilon = 0.01','FontSize',24)
% % %
simstr = 'constant_n0_L1000_fine_ret';
load([matdir filesep simstr],'n','m','xmesh');
[nss,mss,xmeshss,Bss] = SteadyStateCalculator('resmesh','coarse','n0',...
    [zeros(1,40), 20*ones(1,920), zeros(1,40)],'delta',0.01,'epsilon',1);

figure('Units','inches','Position',[0 0 15 7]); 
subplot(1,2,1); hold on;
plot(xmesh,n(end,:),'r','LineWidth',2); 
plot(xmesh,m(end,:),'b','LineWidth',2); 
ylim([0,max([mss,m(end,:)])]);
legend({'Steady State n', 'Steady State m'},'Location','northeast'); 
xlabel('x'); 
ylabel('Concentration'); 
title('t_f for PDE'); 
set(gca,'FontSize',16)
subplot(1,2,2); hold on;
plot(xmeshss,nss,'r','LineWidth',2); 
plot(xmeshss,mss,'b','LineWidth',2); 
ylim([0,max([mss,m(end,:)])]);
xlabel('x'); 
title(sprintf('Steady-state for B = %0.3f',Bss)); 
set(gca,'FontSize',16)
sgtitle('\delta = 0.01, \epsilon = 1','FontSize',24)

% % %
simstr = 'constant_n0_L1000_fine_nob';
load([matdir filesep simstr],'n','m','xmesh');
[nss,mss,xmeshss,Bss] = SteadyStateCalculator('resmesh','coarse','n0',...
    [zeros(1,40), 20*ones(1,920), zeros(1,40)],'delta',1,'epsilon',0.35);
figure('Units','inches','Position',[0 0 15 7]); 
subplot(1,2,1); hold on;
plot(xmesh,n(end,:),'r','LineWidth',2); 
plot(xmesh,m(end,:),'b','LineWidth',2); 
ylim([0,max([mss,m(end,:)])]);
legend({'Steady State n', 'Steady State m'},'Location','northwest'); 
xlabel('x'); 
ylabel('Concentration'); 
title('t_f for PDE'); 
set(gca,'FontSize',16)
subplot(1,2,2); hold on;
plot(xmeshss,nss,'r','LineWidth',2); 
plot(xmeshss,mss,'b','LineWidth',2); 
ylim([0,max([mss,m(end,:)])]);
xlabel('x'); 
title(sprintf('Steady-state for B = %0.3f',Bss)); 
set(gca,'FontSize',16)
sgtitle('\delta = 1, \epsilon = 0.35','FontSize',24)

%% 
[nss1,mss1,xmeshss1,Bss1] = SteadyStateCalculator('resmesh','coarse','n0',...
    [zeros(1,40), 20*ones(1,920), zeros(1,40)]);
[nss2,mss2,xmeshss2,Bss2] = SteadyStateCalculator('resmesh','coarse','n0',...
    [zeros(1,40), 20*ones(1,920), zeros(1,40)],'total_mass',46);
[nss3,mss3,xmeshss3,Bss3] = SteadyStateCalculator('resmesh','coarse','n0',...
    [zeros(1,40), 20*ones(1,920), zeros(1,40)],'total_mass',736);


figure('Units','inches','Position',[0 0 20 7]); 
subplot(1,3,1); hold on;
plot(xmeshss1,nss1,'r','LineWidth',2); 
plot(xmeshss1,mss1,'b','LineWidth',2); 
ylim([0,max([nss1,mss1])]); xlim([0,1400])
legend({'Steady State n', 'Steady State m'},'Location','northwest'); 
xlabel('x'); 
ylabel('Concentration'); 
title(sprintf('Mass = 184, Bias = %0.3f', biascalc(nss1,mss1))); 
set(gca,'FontSize',16)
subplot(1,3,2); hold on;
plot(xmeshss2,nss2,'r','LineWidth',2); hold on;
plot(xmeshss2,mss2,'b','LineWidth',2); 
ylim([0,max([nss2,mss2])]); xlim([0,1400])
legend({'Steady State n', 'Steady State m'},'Location','northwest'); 
xlabel('x'); 
ylabel('Concentration'); 
title(sprintf('Mass = 46, Bias = %0.3f', biascalc(nss2,mss2))); 
set(gca,'FontSize',16)
subplot(1,3,3); hold on;
plot(xmeshss3,nss3,'r','LineWidth',2); 
plot(xmeshss3,mss3,'b','LineWidth',2); 
ylim([0,max(mss3)]);  xlim([0,1400])
legend({'Steady State n', 'Steady State m'},'Location','northwest'); 
xlabel('x'); 
ylabel('Concentration'); 
title(sprintf('Mass = 736, Bias = %0.3f', biascalc(nss3,mss3))); 
set(gca,'FontSize',16)
sgtitle('\delta = 1, \epsilon = 0.01','FontSize',24)