function F_plot_net = BrainframePlot(simstr,simno,tpts,wflow,savenclose,loadpath_,simpath_,bfpath_,figpath_)

addpath(bfpath_)

% Unpack struct
load([simpath_ filesep simstr '.mat'],'output_struct');
% load([loadpath_ filesep 'CCF_labels.mat'],'CCF_labels');
regs = output_struct.Simulations(1).Model_Outputs.Sim.region_names;
% hems = cell(length(regs),1); regs_ = regs;
% for i = 1:length(regs)
%     hems{i} = regs{i}(end-1:end);
%     regs_{i} = regs{i}(1:end-3);
% end
% lh_strfun = @(x) strrep(x,'LH','Left Hemisphere');
% rh_strfun = @(x) strrep(x,'RH','Right Hemisphere');
% hems = cellfun(lh_strfun,hems,'UniformOutput',false);
% hems = cellfun(rh_strfun,hems,'UniformOutput',false);
% regs_426 = CCF_labels(:,1);
% hems_426 = CCF_labels(:,4);
% 
% reg_inds = ismember(regs_426,regs);

ts = output_struct.Simulations(simno).Model_Outputs.Sim.trange(tpts) * 180;
N = output_struct.Simulations(simno).Model_Outputs.Predicted.N(:,tpts);
M = output_struct.Simulations(simno).Model_Outputs.Predicted.M(:,tpts);
X = N + M;
X = DataToCCF_Transport(X,regs,1,loadpath_);
X_plot = X;
X_plot(isnan(X_plot)) = 0;
X_plot = X_plot / max(X(:));
X_plot = X_plot .^ 0.5; % normalization is off, check later 10/12/2023
if wflow
    C = output_struct.Simulations(simno).Model_Outputs.Sim.C;
    F = output_struct.Simulations(simno).Model_Outputs.Predicted.F(:,:,tpts);
    F = C .* F;
    F_plot = DataToCCF_Transport(F,regs,0,loadpath_);
%     F_plot = NaN(426);
%     for j = 1:length(tpts)
%         F_plot(reg_inds,reg_inds,j) = F(:,:,j); 
%     end
    F_plot(isnan(F_plot)) = 0;
    F_plot = F_plot / max(abs(F_plot(:)));
end
seedreg = DataToCCF_Transport(output_struct.Simulations(simno).Model_Outputs.Sim.init_path,regs,1,loadpath_);
seedreg(isnan(seedreg)) = 0; seedreg = logical(seedreg);

% Chunk of code to define region_groups
reggroups = zeros(213,1); 
amy = 1:11; cer = 12:23; sub = 24:26; hip = 27:37; hyp = 38:57;
ncx = 58:95; med = 96:120; mid = 121:141; olf = 142:149; pal = 150:157;
pon = 158:170; str = 171:178; tha = 179:213;
reggroups(amy) = 1; reggroups(cer) = 2; reggroups(sub) = 3; 
reggroups(hip) = 4; reggroups(hyp) = 5; reggroups(ncx) = 6;
reggroups(med) = 7; reggroups(mid) = 8; reggroups(olf) = 9;
reggroups(pal) = 10; reggroups(pon) = 11; reggroups(str) = 12;
reggroups(tha) = 13; 
reggroups = [reggroups;reggroups];
reggroups_conn = reggroups;
reggroups(seedreg) = 14;
cmap = hsv(length(unique(reggroups))-1); %Creating colormap
cmap = [cmap; [1 0 0]];
imgview = [-1.068505859375000e+02,-32.923828125000000];
for i = 1:length(ts)
    imglabel = sprintf([simstr '_t%d'],tpts(i));
    X_plot_i = X_plot(:,i);
    if wflow
        F_plot_i = squeeze(F_plot(:,:,i));
        F_plot_net = zeros(size(F_plot_i));
        for j = 1:size(F_plot_net,1)
            for k = j:size(F_plot_net,1)
                F_jk = F_plot_i(j,k);
                F_kj = F_plot_i(k,j);
                F_net = F_jk - F_kj;
                if F_net > 0 
                    F_plot_net(j,k) = F_net;
                elseif F_net < 0
                    F_plot_net(k,j) = abs(F_net);
                end
            end
        end
        F_plot_i = F_plot_net;
        F_plot_i(F_plot_i < prctile(nonzeros(F_plot_i(:)),95)) = 0;
        input_struct = brainframe_inputs_mouse(bfpath_,'conmat',F_plot_i,...
                                                     'region_groups',reggroups,...
                                                     'con_regiongroups',reggroups_conn,...
                                                     'cmap',cmap,...
                                                     'con_cmap',cmap,...
                                                     'xfac',5,...
                                                     'sphere',1,...
                                                     'sphere_npts',20,...
                                                     'pointsize',5,...
                                                     'voxUreg',1,...
                                                     'iscon',1,...
                                                     'conarrow_WL',[1.5 1],...
                                                     'data',X_plot_i,...
                                                     'norm_method','none',...
                                                     'bgcolor','w',...
                                                     'con_rescale',20,...
                                                     'img_labels',imglabel,...
                                                     'img_format','tiffn',...
                                                     'img_views',imgview,...
                                                     'img_directory',figpath_,...
                                                     'savenclose',savenclose,...
                                                     'con_arch',0.3);
    else
        input_struct = brainframe_inputs_mouse(bfpath_,'region_groups',reggroups,...
                                                     'cmap',cmap,...
                                                     'con_cmap',cmap,...
                                                     'xfac',3,...
                                                     'sphere',1,...
                                                     'sphere_npts',15,...
                                                     'pointsize',5,...
                                                     'voxUreg',1,...
                                                     'data',X_plot_i,...
                                                     'norm_method','max',...
                                                     'bgcolor','w',...
                                                     'con_rescale',20,...
                                                     'img_labels',imglabel,...
                                                     'img_format','tiffn',...
                                                     'img_views',imgview,...
                                                     'img_directory',figpath_,...
                                                     'savenclose',savenclose,...
                                                     'con_arch',0.3);
    end
    brainframe(input_struct);
end


% if ~isfield(output_struct.Simulations(idx).Model_Outputs.Sim,'region_names')
%     load([loadpath_ filesep 'CCF_labels.mat'],'CCF_labels');
%     switch output_struct.Simulations(idx).Model_Outputs.Sim.connectome_subset
%         case 'Hippocampus'
%             inds = ismember(CCF_labels(:,3),'Hippocampus');
%         case 'Hippocampus+PC+RSP'
%             inds_hipp = ismember(CCF_labels(:,3),'Hippocampus');
%             inds_pc = ismember(CCF_labels(:,1),'Piriform area');
%             inds_rsp = ismember(CCF_labels(:,3),'Retrosplenial Area');
%             inds = logical(inds_hipp + inds_pc + inds_rsp);
%         case 'RH'
%             inds = ismember(CCF_labels(:,4),'Right Hemisphere');
%         case 'LH'
%             inds = ismember(CCF_labels(:,4),'Left Hemisphere');
%         otherwise
%             inds = logical(ones(size(Conn,1),1)); %#ok<LOGL> 
%     end
%     regnamecell = CCF_labels(inds,:);
%     regnames = cell(size(regnamecell,1),1);
%     for i = 1:length(regnames)
%         regname = regnamecell{i,1};
%         reghem = regnamecell{i,4};
%         if strcmp(reghem,'Right Hemisphere')
%             regnames{i} = [regname ' RH'];
%         else
%             regnames{i} = [regname ' LH'];
%         end
%     end
% else
%     regnames = output_struct.Simulations(idx).Model_Outputs.Sim.region_names;
% end
% initpath = output_struct.Simulations(idx).Model_Outputs.Sim.init_path;
% if ~isfield(output_struct.Simulations(idx).Model_Outputs.Sim,'seed')
%     seedreg = regnames(initpath > 0);
% else
%     seedreg = output_struct.Simulations(idx).Model_Outputs.Sim.seed;
% end
% if length(seedreg) > 1
%     s2 = cell(1,length(seedreg)); s2(1:end-1) = {', '}; s2(end) = {''};
%     seedreg = [seedreg; s2]; seedreg = [seedreg{:}];
% else
%     seedreg = char(seedreg);
% end





end