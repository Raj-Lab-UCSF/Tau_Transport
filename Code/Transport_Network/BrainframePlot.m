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
X_plot(X_plot < 0.25*median(nonzeros(X_plot(:)))) = 0;
X_plot = X_plot .^ (1/3); 
if wflow
    C = output_struct.Simulations(simno).Model_Outputs.Sim.C;
    F = output_struct.Simulations(simno).Model_Outputs.Predicted.F(:,:,tpts);
    F = C .* F;
    F_plot = DataToCCF_Transport(F,regs,0,loadpath_);
    F_plot(isnan(F_plot)) = 0;
    F_plot = F_plot / max(abs(F_plot(:)));
%     F_plot(F_plot < 0.05*max(nonzeros(F_plot(:)))) = 0;
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
reggroups(seedreg) = 14;
reggroups_conn = reggroups;
cmap = lines(length(unique(reggroups))-1); %Creating colormap
cmap = [cmap; [0 1 0.5]];
imgview = [-1 0 0];
for i = 1:length(ts)
    imglabel = sprintf([simstr '_' num2str(simno) '_t%d'],tpts(i));
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
                                                     'xfac',4,...
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
                                                     'img_directory',figpath_,...
                                                     'img_renderer','painters',...
                                                     'img_view',imgview,...
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
                                                     'norm_method','none',...
                                                     'bgcolor','w',...
                                                     'con_rescale',20,...
                                                     'img_labels',imglabel,...
                                                     'img_format','tiffn',...
                                                     'img_directory',figpath_,...
                                                     'img_renderer','painters',...
                                                     'img_view',imgview,...
                                                     'savenclose',savenclose,...
                                                     'con_arch',0.3);
    end
    brainframe(input_struct);
end

end