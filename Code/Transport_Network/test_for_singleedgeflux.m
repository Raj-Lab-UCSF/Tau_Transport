clear; clc;
%% Test for SingleEdgeCalc
matdir = '~/Documents/MATLAB/Tau_Transport_OtherFiles/FinalSimFiles';
loadpath = '~/Documents/MATLAB/Tau_Transport/MatFiles';
rng(0);
load([matdir filesep 'hippocampome_final_round2.mat'],'output_struct');
ntests = 100;
Fs_true = NaN(ntests,1); Fs_test = Fs_true;
nsims = size(output_struct.Parameter_Grid,1);
nroi = size(output_struct.Simulations(1).Model_Outputs.Predicted.N,1);
nt = size(output_struct.Simulations(1).Model_Outputs.Predicted.N,2);
for i = 1:ntests
    randsimno = datasample(1:nsims,1);
    % randt = datasample(1:nt,1);
    randt = 1; % Have to test at the first time point prior to weight correction on N
    F_i = squeeze(output_struct.Simulations(randsimno).Model_Outputs.Predicted.F(:,:,randt));
    nonzero_inds = find(F_i);
    randroi_ind = datasample(nonzero_inds,1);
    [randroi_1, randroi_2] = ind2sub([nroi,nroi],randroi_ind);
    Fs_true(i) = squeeze(output_struct.Simulations(randsimno).Model_Outputs.Predicted.F(randroi_1,randroi_2,randt));
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
end

figure; hold on;
scatter(Fs_test,Fs_true,'ko');
h = lsline;
xlabel('F (Test)'); ylabel('F (True)');
xlim([1.1*min(Fs_test),1.1*max(Fs_test)]);
ylim([1.1*min(Fs_true),1.1*max(Fs_true)]);
legend(h,{sprintf('R = %.2f',corr(Fs_test,Fs_true))},'Location','northwest');
set(gca,'FontSize',20,'FontName','Times');

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