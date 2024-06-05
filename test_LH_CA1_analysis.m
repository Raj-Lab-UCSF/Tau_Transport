%% Loading
clear; clc;
loadpath = '/Users/justintorok/Documents/MATLAB/Tau_Transport_OtherFiles/FinalSimFiles';
loadpath2 = '/Users/justintorok/Documents/MATLAB/Nexis_Project/Nexis/raw_data_mouse';
addpath('/Users/justintorok/Documents/MATLAB/Nexis_Project/Nexis/MATLAB/lib_NexIS');
load([loadpath2 filesep 'Mouse_Tauopathy_Data_HigherQ.mat'])
load([loadpath filesep 'LH_CA1seed_dir.mat'])

%% Defining X & calculating Ant/Ret corrs
simno_ret = 14;
Nret = output_struct.Simulations(simno_ret).Model_Outputs.Predicted.N;
Mret = output_struct.Simulations(simno_ret).Model_Outputs.Predicted.M;
Xret = Nret + Mret;
trange_ret = output_struct.Simulations(simno_ret).Model_Outputs.Sim.trange;

simno_ant = 3;
Nant = output_struct.Simulations(simno_ant).Model_Outputs.Predicted.N;
Mant = output_struct.Simulations(simno_ant).Model_Outputs.Predicted.M;
Xant = Nant + Mant;
trange_ant = output_struct.Simulations(simno_ant).Model_Outputs.Sim.trange;

seedbool = (Xret(:,1) > 0);
Xret_noseed = Xret; Xret_noseed(seedbool,:) = [];
Xant_noseed = Xant; Xant_noseed(seedbool,:) = [];

R_antret = diag(corr(Xret_noseed,Xant_noseed));

%% Plot Ant/Ret corrs over time
figure; hold on;
scatter(trange_ret,R_antret,'gx');
xlabel('t (Model)'); ylabel('R (no seed)'); title(['Simulations: ' num2str(simno_ret) ' & ' ...
    num2str(simno_ant) ', Ant/Ret Corrs'])
set(gca,'FontSize',20','FontName','Times')

%% Calculate corrs to incoming/outgoing seed connectivity
C = output_struct.Simulations(simno_ret).Model_Outputs.Sim.C;
C_in_seed = C(:,seedbool);
C_out_seed = C(seedbool,:).';
C_in_seed(seedbool) = []; C_out_seed(seedbool) = [];

corrs_Cin_ret = corr([C_in_seed, Xret_noseed]); corrs_Cin_ret = corrs_Cin_ret(1,2:end);
corrs_Cout_ret = corr([C_out_seed, Xret_noseed]); corrs_Cout_ret = corrs_Cout_ret(1,2:end);

corrs_Cin_ant = corr([C_in_seed, Xant_noseed]); corrs_Cin_ant = corrs_Cin_ant(1,2:end);
corrs_Cout_ant = corr([C_out_seed, Xant_noseed]); corrs_Cout_ant = corrs_Cout_ant(1,2:end);

%% Plot
figure; hold on;
scatter(trange_ret,corrs_Cin_ret,'bo'); 
scatter(trange_ret,corrs_Cout_ret,'ro'); 
legend({'C_i_n','C_o_u_t'});
xlabel('t Model'); ylabel('R'); title(['Sim. Number: ' num2str(simno_ret) ', Retrograde Condition']);
set(gca,'FontSize',20','FontName','Times')

figure; hold on;
scatter(trange_ret,corrs_Cin_ant,'bo'); 
scatter(trange_ret,corrs_Cout_ant,'ro'); 
legend({'C_i_n','C_o_u_t'});
xlabel('t Model'); ylabel('R'); title(['Sim. Number: ' num2str(simno_ant) ', Anterograde Condition']);
set(gca,'FontSize',20','FontName','Times')

%% Calculate correlations to data
studynames = {'DS4','DS6','DS6_110','DS7','DS9','DS9_110'};
corrs_data = struct;
Xret_ccf = [NaN(size(Xret)); Xret]; Xant_ccf = [NaN(size(Xant)); Xant];
for j = 1:length(studynames)
    timestamps_j =  mousedata_struct.(studynames{j}).time_stamps;
    corrs_data.(studynames{j}).data_times = timestamps_j;
    for i = 1:length(timestamps_j)
        data_ij = mousedata_struct.(studynames{j}).data(:,i);
        seed_ij = logical(mousedata_struct.(studynames{j}).seed);
        Xret_ij = CCFToData(Xret_ccf,studynames{j},loadpath2);
        Xant_ij = CCFToData(Xant_ccf,studynames{j},loadpath2);
        % if ~isnan(seed_ij)
        %     data_ij(seed_ij) = [];
        %     Xret_ij(seed_ij,:) = [];
        %     Xant_ij(seed_ij,:) = [];
        % end
        corrret = corr([data_ij,Xret_ij],'rows','complete');
        corrant = corr([data_ij,Xant_ij],'rows','complete');

        corrs_data.(studynames{j}).ret.model_times = trange_ret;
        corrs_data.(studynames{j}).ant.model_times = trange_ant;
        corrs_data.(studynames{j}).ret.(sprintf('R_%d',i)) = corrret(1,2:end);
        corrs_data.(studynames{j}).ant.(sprintf('R_%d',i)) = corrant(1,2:end);
    end
end

%% Plot
shapes_plot = {'-','--',':'};
legend_str = {'R_r_e_t, t_1', 'R_r_e_t, t_2', 'R_r_e_t, t_3',...
              'R_a_n_t, t_1', 'R_a_n_t, t_2', 'R_a_n_t, t_3'};
for j = 1:length(studynames)
    figure('Units','inches','Position',[0 0 10 10]); hold on;
    corrs_data_j = corrs_data.(studynames{j});
    for i = 1:length(shapes_plot)
        t_ret = corrs_data_j.ret.model_times;
        corrret_i = corrs_data_j.ret.(sprintf('R_%d',i));
        plot(t_ret,corrret_i,'b','LineStyle',shapes_plot{i},'LineWidth',2);
    end
    for i = 1:length(shapes_plot)
        t_ant = corrs_data_j.ant.model_times;
        corrant_i = corrs_data_j.ant.(sprintf('R_%d',i));
        plot(t_ant,corrant_i,'r','LineStyle',shapes_plot{i},'LineWidth',2); 
    end
    studylabel = strrep(studynames{j},'_',' ');
    legend(legend_str,'Location','southwest','FontName','Times',...
        'FontSize',18,'NumColumns',2);
    xlabel('Model Time'); ylabel('R (w/seed)');
    title(sprintf('Simulations %d & %d: Corrs with Study %s',simno_ret,...
        simno_ant,studylabel));
    set(gca,'FontName','Times','FontSize',24);
end

%% 
for i = 1:3
    tpt = i;
    for j = 1:length(studynames)
        figure; hold on; scatter(trange_ret, corrret,'bo'); scatter(trange_ret, corrant,'ro'); legend({'ret','ant'})
        xlabel('t'); ylabel('R (no seed)'); 
        title(['Sim. Number: ' num2str(simno) ', Study: ' studynames{j} ', t = ' num2str(i)]);
        set(gca,'FontSize',20,'FontName','Times');
    end
end

for j = 1:length(studynames)
    figure; hold on;
    for i = 1:3
        tpt = i;
        data_ij = mousedata_struct_ccf.(studynames{j}).data(:,tpt);
        corr_Cin = corr(C_in_seed,data_ij(hippinds_ccf),'rows','complete');
        corr_Cout = corr(C_out_seed,data_ij(hippinds_ccf),'rows','complete');
        scatter(tpt,corr_Cin,'bo','filled'); scatter(tpt,corr_Cout,'ro','filled');
        legend({'C_i_n','C_o_u_t'});
    end
    xlabel('t'); ylabel('R (no seed)'); 
    title(['R with Conn to CA1, Study: ' studylabel]);
    set(gca,'FontSize',20,'FontName','Times');
end