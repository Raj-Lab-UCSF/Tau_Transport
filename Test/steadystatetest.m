matdir = '/Users/justintorok/Documents/MATLAB/Tau_Transport/SampleFiles';

% % % 
simstr = 'constant_n0_L1000_fine_ant';
load([matdir filesep simstr],'n','m','xmesh');
[nss,mss,xmeshss,Bss] = SteadyStateCalculator('resmesh','coarse','n0',[zeros(1,40), 20*ones(1,920), zeros(1,40)]);

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