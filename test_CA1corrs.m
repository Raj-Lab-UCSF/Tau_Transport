clear; clc;
loadpath = '/Users/justintorok/Documents/MATLAB/Tau_Transport_OtherFiles/FinalSimFiles';
loadpath2 = '/Users/justintorok/Documents/MATLAB/Tau_Transport/MatFiles';
load([loadpath2 filesep 'Mouse_Tauopathy_Data_HigherQ_CCF.mat'])
simno = 13;
% studynames = {'DS4','DS6','DS7','DS9','DS6_110','DS9_110'};
studynames = {'DS4','DS9'};
load([loadpath filesep 'hippocampome_CA1seed_dir_ant.mat'])
Nant = output_struct.Simulations(simno).Model_Outputs.Predicted.N;
Mant = output_struct.Simulations(simno).Model_Outputs.Predicted.M;
Xant = Nant+Mant;

load([loadpath filesep 'hippocampome_CA1seed_dir_ret.mat'])
Nret = output_struct.Simulations(simno).Model_Outputs.Predicted.N;
Mret = output_struct.Simulations(simno).Model_Outputs.Predicted.M;
Xret = Nret+Mret;
trange = output_struct.Simulations(simno).Model_Outputs.Sim.trange;

C = output_struct.Simulations(simno).Model_Outputs.Sim.C;
C_in = C(:,19);
C_out = C(19,:).';

Xant(19,:) = []; Xret(19,:) = []; 
C_in(19) = []; C_out(19) = [];

figure; hold on;
scatter(trange,diag(corr(Xant,Xret)),'gx');
xlabel('t (Model)'); ylabel('R (no seed)'); title(['Simulation: ' num2str(simno) ', Ant/Ret Corrs'])
set(gca,'FontSize',20','FontName','Times')

hippinds_ccf = [27:37,78:80,147];
hippinds_ccf = [hippinds_ccf,hippinds_ccf+213]; hippinds_ccf(hippinds_ccf==243) = [];

corrCin_ret = corr([C_in, Xret]); corrCin_ret = corrCin_ret(1,2:end);
corrCout_ret = corr([C_out, Xret]); corrCout_ret = corrCout_ret(1,2:end);
corrCin_ant = corr([C_in, Xant]); corrCin_ant = corrCin_ant(1,2:end);
corrCout_ant = corr([C_out, Xant]); corrCout_ant = corrCout_ant(1,2:end);

figure; hold on;
scatter(trange,corrCin_ret,'bo'); 
scatter(trange,corrCout_ret,'ro'); 
legend({'C_i_n','C_o_u_t'});
xlabel('t Model'); ylabel('R'); title(['Sim. Number: ' num2str(simno) ', Retrograde Condition']);
set(gca,'FontSize',20','FontName','Times')

figure; hold on;
scatter(trange,corrCin_ant,'bo'); 
scatter(trange,corrCout_ant,'ro'); 
legend({'C_i_n','C_o_u_t'});
xlabel('t Model'); ylabel('R'); title(['Sim. Number: ' num2str(simno) ', Anterograde Condition']);
set(gca,'FontSize',20','FontName','Times')

for i = 1:3
    tpt = i;
    for j = 1:length(studynames)
        DSdata = mousedata_struct_ccf.(studynames{j}).data(:,tpt);
        corrret = corr([DSdata(hippinds_ccf),Xret],'rows','complete'); corrret = corrret(1,2:end);
        corrant = corr([DSdata(hippinds_ccf),Xant],'rows','complete'); corrant = corrant(1,2:end);
        figure; hold on; scatter(trange, corrret,'bo'); scatter(trange, corrant,'ro'); legend({'ret','ant'})
        xlabel('t'); ylabel('R (no seed)'); 
        title(['Sim. Number: ' num2str(simno) ', Study: ' studynames{j} ', t = ' num2str(i)]);
        set(gca,'FontSize',20,'FontName','Times');
    end
end

for j = 1:length(studynames)
    figure; hold on;
    for i = 1:3
        tpt = i;
        DSdata = mousedata_struct_ccf.(studynames{j}).data(:,tpt);
        corr_Cin = corr(C_in,DSdata(hippinds_ccf),'rows','complete');
        corr_Cout = corr(C_out,DSdata(hippinds_ccf),'rows','complete');
        scatter(tpt,corr_Cin,'bo','filled'); scatter(tpt,corr_Cout,'ro','filled');
        legend({'C_i_n','C_o_u_t'});
    end
    xlabel('t'); ylabel('R (no seed)'); 
    title(['R with Conn to CA1, Study: ' studynames{j}]);
    set(gca,'FontSize',20,'FontName','Times');
end