clear; clc;
%% Test for SingleEdgeCalc
matdir = '~/Documents/MATLAB/Tau_Transport_OtherFiles/FinalSimFiles';
loadpath = '~/Documents/MATLAB/Tau_Transport/MatFiles';
rng(0);
load([matdir filesep 'hippocampus_testmat_for_single_edge_debug.mat'],'output_struct');
ntests = 250;
Fs_true = NaN(ntests,1); Fs_test = Fs_true;
W1s_true = NaN(ntests,1); W1s_test = W1s_true;
nsims = size(output_struct.Parameter_Grid,1);
nroi = size(output_struct.Simulations(1).Model_Outputs.Predicted.N,1);
nt = size(output_struct.Simulations(1).Model_Outputs.Predicted.N,2);
for i = 1:ntests
    fprintf('Random iteration %d/%d\n',i,ntests)
    randsimno = datasample(1:nsims,1);
    % randt = datasample(1:nt,1);
    randt = 1; % Have to test at the first time point prior to weight correction on N
    F_i = squeeze(output_struct.Simulations(randsimno).Model_Outputs.Predicted.F(:,:,randt));
    nonzero_inds = find(F_i);
    randroi_ind = datasample(nonzero_inds,1);
    [randroi_1, randroi_2] = ind2sub([nroi,nroi],randroi_ind);
    Fs_true(i) = squeeze(output_struct.Simulations(randsimno).Model_Outputs.Predicted.F(randroi_1,randroi_2,randt));
    W1s_true(i) = squeeze(output_struct.Simulations(randsimno).Model_Outputs.Predicted.W1(randroi_1,randroi_2,randt));
    tau_x0_i = output_struct.Simulations(randsimno).Model_Outputs.Predicted.N(randroi_1,randt);
    tau_xL_i = output_struct.Simulations(randsimno).Model_Outputs.Predicted.N(randroi_2,randt);
    paramstruct_i = output_struct.Simulations(randsimno).Model_Outputs.Parameters;
    mdloutput = SingleEdgeFluxCalculator_WrapperFun(loadpath,'gamma1',paramstruct_i.gamma1,...
                                    'lambda1',paramstruct_i.lambda1,...
                                    'delta',paramstruct_i.delta,...
                                    'epsilon',paramstruct_i.epsilon,...
                                    'tau_x0',tau_x0_i,...
                                    'tau_xL',tau_xL_i,... % use paramlist(i) instead of paramnames 
                                    'beta',paramstruct_i.beta,...
                                    'frac',paramstruct_i.frac,...
                                    'lambda2',paramstruct_i.lambda2,... % same as lambda1
                                    'gamma2',paramstruct_i.gamma2);
    Fs_test(i) = mdloutput.Predicted.F;
    W1_i = SingleEdgeWCalculator(loadpath,Fs_true(i),...
                                    'gamma1',paramstruct_i.gamma1,...
                                    'lambda1',paramstruct_i.lambda1,...
                                    'delta',paramstruct_i.delta,...
                                    'epsilon',paramstruct_i.epsilon,...
                                    'tau_x0',tau_x0_i,...
                                    'tau_xL',tau_xL_i,... % use paramlist(i) instead of paramnames 
                                    'beta',paramstruct_i.beta,...
                                    'frac',paramstruct_i.frac,...
                                    'lambda2',paramstruct_i.lambda2,... % same as lambda1
                                    'gamma2',paramstruct_i.gamma2);
    W1s_test(i) = W1_i;
end

figure('Units','inches','Position',[0 0 15 16]); 
subplot(2,2,1); hold on;       
scatter(Fs_test,Fs_true,'bo');
xlabel('F_t_e_s_t'); ylabel('F_t_r_u_e');
xlim([1.1*min(Fs_test),1.1*max(Fs_test)]);
ylim([1.1*min(Fs_true),1.1*max(Fs_true)]);
h = lsline;
legend(h,{sprintf('R^2 = %.2f',corr(Fs_test,Fs_true)^2)},'Location','northwest');
set(gca,'FontSize',20,'FontName','Times','box','on');

subplot(2,2,2); hold on;       
scatter(W1s_test,W1s_true,'ro');
xlabel('W1_t_e_s_t'); ylabel('W1_t_r_u_e');
xlim([0.9*min(W1s_test),1.1*max(W1s_test)]);
ylim([0.9*min(W1s_true),1.1*max(W1s_true)]);
h = lsline;
legend(h,{sprintf('R^2 = %.2f',corr(W1s_test,W1s_true)^2)},'Location','northwest');
set(gca,'FontSize',20,'FontName','Times','box','on');

subplot(2,2,3); hold on;       
scatter(Fs_true,(Fs_test - Fs_true),'bo');
plot([1.1*min(Fs_true),1.1*max(Fs_true)],[0 0],'LineStyle','--','Color','k')
xlabel('F_t_r_u_e'); ylabel('F_t_e_s_t - F_t_r_u_e');
xlim([1.1*min(Fs_true),1.1*max(Fs_true)]);
ylim([1.1*min(Fs_test - Fs_true),1.1*max(Fs_test - Fs_true)]);
set(gca,'FontSize',20,'FontName','Times','box','on');

subplot(2,2,4); hold on;       
scatter(W1s_true,(W1s_test - W1s_true),'ro');
plot([1.1*min(W1s_true),1.1*max(W1s_true)],[0 0],'LineStyle','--','Color','k')
xlabel('W1_t_r_u_e'); ylabel('W1_t_e_s_t - W1_t_r_u_e');
xlim([1.1*min(W1s_true),1.1*max(W1s_true)]);
ylim([1.1*min(W1s_test - W1s_true),1.1*max(W1s_test - W1s_true)]);
set(gca,'FontSize',20,'FontName','Times','box','on');

%% Test for SR
matdir = '~/Documents/MATLAB/Tau_Transport/SampleFiles';
outputs_int = load([matdir filesep 'test_for_flux.mat'],'output_struct');
outputs_sr = load([matdir filesep 'test_for_flux_sr.mat'],'output_struct');
Fs_int = NaN(100,1);
Fs_sr = Fs_int;

for i = 1:100
    Fs_int(i) = outputs_int.output_struct.Simulations(i).Model_Outputs.Predicted.F;
    Fs_sr(i) = outputs_sr.output_struct.Simulations(i).Model_Outputs.Predicted.F;
end

figure; hold on;
scatter(Fs_sr,Fs_int,'ko');
h = lsline;
xlabel('F (SR)'); ylabel('F (Int)');
legend(h,{sprintf('R = %.2f',corr(Fs_sr,Fs_int))});
set(gca,'FontSize',16,'FontName','Times');