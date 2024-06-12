%% Loading
clear; clc;
loadpath = '/Users/justintorok/Documents/MATLAB/Tau_Transport_OtherFiles/FinalSimFiles';
loadpath2 = '/Users/justintorok/Documents/MATLAB/Nexis_Project/Nexis/raw_data_mouse';
loadpath3 = '~/Documents/MATLAB/Tau_Transport/MatFiles';
addpath('/Users/justintorok/Documents/MATLAB/Nexis_Project/Nexis/MATLAB/lib_NexIS');
addpath('/Users/justintorok/Documents/MATLAB/Nexis_Project/Nexis/TauDirectionality_Related');
load([loadpath2 filesep 'Mouse_Tauopathy_Data_HigherQ.mat'])
simstr = 'LH_CA1seed_dir';
load([loadpath filesep simstr '.mat'])
figpath = cd;

%% Defining X & calculating Ant/Ret corrs
simno_ret = 2;
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
figure('Units','inches','Position',[0 0 6.5 6]); hold on;
scatter(trange_ret,corrs_Cout_ret,'bo'); 
scatter(trange_ret,corrs_Cin_ret,'ro'); 
legend({'C_o_u_t','C_i_n'});
xlabel('t (Model)'); 
xticks([0 0.5 1]); yticks([0,0.3,0.6])
% title(['Sim. Number: ' num2str(simno_ret) ', Retrograde Condition']);
set(gca,'FontSize',24','FontName','Times','box','on')
print('Cinout_LH_ret','-dtiff','-r300'); close;

figure('Units','inches','Position',[0 0 6.5 6]); hold on;
scatter(trange_ret,corrs_Cout_ant,'bo'); 
scatter(trange_ret,corrs_Cin_ant,'ro'); 
legend({'C_o_u_t','C_i_n'});
xlabel('t (Model)'); ylabel("Pearson's R"); 
xticks([0 0.5 1]); yticks([-0.2,0.4,1])
% title(['Sim. Number: ' num2str(simno_ant) ', Anterograde Condition']);
set(gca,'FontSize',24','FontName','Times','box','on')
print('Cinout_LH_ant','-dtiff','-r300'); close;

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

%% Plot DS6/DS7, t = end
shapes_plot = {'-','--'};
studynames_plot = {'DS6','DS7'};
legend_str = {'R_a_n_t, t = 3 months','R_r_e_t, t = 3 months'};
for j = 1:length(studynames_plot)
    figure('Units','inches','Position',[0 0 6.5 6]); hold on;
    corrs_data_j = corrs_data.(studynames_plot{j});

    t_ant = corrs_data_j.ant.model_times;
    corrant_i = corrs_data_j.ant.(sprintf('R_%d',i));
    plot(t_ant,corrant_i,'b','LineStyle',shapes_plot{j},'LineWidth',2); 
    t_ret = corrs_data_j.ret.model_times;
    corrret_i = corrs_data_j.ret.(sprintf('R_%d',3));
    plot(t_ret,corrret_i,'r','LineStyle',shapes_plot{j},'LineWidth',2);

    legend(legend_str,'Location','southwest','FontName','Times');
    xlabel('t (Model)'); 
    if j == 1
        xticks([0 0.5 1]); yticks([-0.3,0.2,0.7])
        ylim([-0.45,0.85]);
        yticklabels({'0.3','0.2','0.7'})
    else
        xticks([0 0.5 1]); yticks([0,0.4,0.8])
        yticklabels({'0','-0.4','0.8'})
        ylim([-0.1,0.9]);
        ylabel("Pearson's R");
    end
    title(['Study ' studynames_plot{j}]);
    set(gca,'FontName','Times','FontSize',24,'box','on');
    print([studynames_plot{j} '_LHcorrs'],'-dtiff','-r300'); close;
end


%% Plot all
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

%% Brainframe
bfpath = '/Users/justintorok/Documents/MATLAB/Brainframe-Dev/Brainframe';
addpath(bfpath);

reggroups_ = zeros(213,1); %Chunk of code to define region_groups
amy = 1:11; cer = 12:23; sub = 24:26; hip = 27:37; hyp = 38:57;
ncx = 58:95; med = 96:120; mid = 121:141; olf = 142:149; pal = 150:157;
pon = 158:170; str = 171:178; tha = 179:213;
reggroups_(amy) = 1; reggroups_(cer) = 2; reggroups_(sub) = 3; 
reggroups_(hip) = 4; reggroups_(hyp) = 5; reggroups_(ncx) = 6;
reggroups_(med) = 7; reggroups_(mid) = 8; reggroups_(olf) = 9;
reggroups_(pal) = 10; reggroups_(pon) = 11; reggroups_(str) = 12;
reggroups_(tha) = 13;
reggroups_ = [reggroups_;reggroups_];

ts = output_struct.Simulations(simno_ant).Model_Outputs.Sim.trange;
studynames_bf = {'DS6','DS7'};
[~,tind1] = max(corrs_data.DS6.ret.R_3);
[~,tind2] = max(corrs_data.DS7.ant.R_3);
wflow = 0; savenclose = 1;
cmap_ = hsv(length(unique(reggroups_)));

% Transform data & pred to CCF space, obtain queried time point
datinput_data_DS7 = DataToCCF([],studynames_bf{2},loadpath2);
datinput_data_DS6 = DataToCCF([],studynames_bf{1},loadpath2);
isnans_data = isnan(datinput_data_DS7(:,1));
datinput_data_DS7(isnans_data,:) = 0;
datinput_data_DS7(1:213,:) = 0; % LH only
datinput_data_DS6(isnans_data,:) = 0;
datinput_data_DS6(1:213,:) = 0; % LH only

Xant_plot = [zeros(size(Xant,1),1); Xant(:,tind2)];
Xant_plot(isnans_data) = 0;

Xret_plot = Xret_ccf(:,tind1);
Xret_plot(isnans_data) = 0; Xret_plot(isnan(Xret_plot)) = 0;

reggroups_data_ = reggroups_;
cmap_data_ = cmap_;
tptsplotind_ = 3;
datinput_data_DS7 = datinput_data_DS7(:,tptsplotind_);
datinput_data_DS6 = datinput_data_DS6(:,tptsplotind_);
tpt = tptsplotind_;

% Thresholding
threshval = 50;
datathresh_DS6_val = prctile(nonzeros(datinput_data_DS6),threshval);
datathresh_DS7_val = prctile(nonzeros(datinput_data_DS7),threshval);
datathresh_ret_val = prctile(nonzeros(Xret_plot),threshval);
datathresh_ant_val = prctile(nonzeros(Xant_plot),threshval);

thresh_inds_DS6 = (datinput_data_DS6 >= datathresh_DS6_val);
datinput_data_DS6(~thresh_inds_DS6) = 0;
thresh_inds_DS7 = (datinput_data_DS7 >= datathresh_DS7_val);
datinput_data_DS7(~thresh_inds_DS7) = 0;
thresh_inds_ret = (Xret_plot >= datathresh_ret_val);
Xret_plot(~thresh_inds_ret) = 0;
thresh_inds_ant = (Xant_plot >= datathresh_ant_val);
Xant_plot(~thresh_inds_ant) = 0;

% Generate glass brains
imglabel = 'DS7_t3';
imgview = [-90,-18];
savenclose = 1;
input_struct_data = brainframe_inputs_mouse(bfpath,...
                                             'region_groups',reggroups_data_,...
                                             'cmap',cmap_data_,...
                                             'xfac',4.5,...
                                             'sphere',1,...
                                             'sphere_npts',50,...
                                             'pointsize',5,...
                                             'voxUreg',1,...
                                             'data',datinput_data_DS7,...
                                             'norm_method','max',...
                                             'bgcolor','w',...
                                             'img_labels',imglabel,...
                                             'img_format','tiffn',...
                                             'img_views',imgview,...
                                             'img_directory',figpath,...
                                             'savenclose',savenclose);
brainframe(input_struct_data);

imglabel = 'DS6_t3';
imgview = [-90,-18];
input_struct_data = brainframe_inputs_mouse(bfpath,...
                                             'region_groups',reggroups_data_,...
                                             'cmap',cmap_data_,...
                                             'xfac',3,...
                                             'sphere',1,...
                                             'sphere_npts',50,...
                                             'pointsize',5,...
                                             'voxUreg',1,...
                                             'data',datinput_data_DS6,...
                                             'norm_method','max',...
                                             'bgcolor','w',...
                                             'img_labels',imglabel,...
                                             'img_format','tiffn',...
                                             'img_views',imgview,...
                                             'img_directory',figpath,...
                                             'savenclose',savenclose);
brainframe(input_struct_data);

imglabel = 'Antmodel_Rmax';
imgview = [-90,-18];
input_struct_data = brainframe_inputs_mouse(bfpath,...
                                             'region_groups',reggroups_data_,...
                                             'cmap',cmap_data_,...
                                             'xfac',4.5,...
                                             'sphere',1,...
                                             'sphere_npts',50,...
                                             'pointsize',5,...
                                             'voxUreg',1,...
                                             'data',Xant_plot,...
                                             'norm_method','max',...
                                             'bgcolor','w',...
                                             'img_labels',imglabel,...
                                             'img_format','tiffn',...
                                             'img_views',imgview,...
                                             'img_directory',figpath,...
                                             'savenclose',savenclose);
brainframe(input_struct_data);

imglabel = 'Retmodel_Rmax';
imgview = [-90,-18];
input_struct_data = brainframe_inputs_mouse(bfpath,...
                                             'region_groups',reggroups_data_,...
                                             'cmap',cmap_data_,...
                                             'xfac',4.5,...
                                             'sphere',1,...
                                             'sphere_npts',50,...
                                             'pointsize',5,...
                                             'voxUreg',1,...
                                             'data',Xret_plot,...
                                             'norm_method','max',...
                                             'bgcolor','w',...
                                             'img_labels',imglabel,...
                                             'img_format','tiffn',...
                                             'img_views',imgview,...
                                             'img_directory',figpath,...
                                             'savenclose',savenclose);
brainframe(input_struct_data);

%%
% for i = 1:3
%     tpt = i;
%     for j = 1:length(studynames)
%         figure; hold on; 
%         corrret_i = corrs_data_j.ret.(sprintf('R_%d',i));
%         corrant_i = corrs_data_j.ant.(sprintf('R_%d',i));
%         scatter(trange_ret, corrret_i,'bo'); 
%         scatter(trange_ret, corrant_i,'ro'); 
%         legend({'ret','ant'})
%         xlabel('t'); ylabel('R (no seed)'); 
%         title(['Sim. Number: ' num2str(simno) ', Study: ' studynames{j} ', t = ' num2str(i)]);
%         set(gca,'FontSize',20,'FontName','Times');
%     end
% end
% 
% for j = 1:length(studynames)
%     figure; hold on;
%     for i = 1:3
%         tpt = i;
%         data_ij = mousedata_struct_ccf.(studynames{j}).data(:,tpt);
%         corr_Cin = corr(C_in_seed,data_ij(hippinds_ccf),'rows','complete');
%         corr_Cout = corr(C_out_seed,data_ij(hippinds_ccf),'rows','complete');
%         scatter(tpt,corr_Cin,'bo','filled'); scatter(tpt,corr_Cout,'ro','filled');
%         legend({'C_i_n','C_o_u_t'});
%     end
%     xlabel('t'); ylabel('R (no seed)'); 
%     title(['R with Conn to CA1, Study: ' studylabel]);
%     set(gca,'FontSize',20,'FontName','Times');
% end