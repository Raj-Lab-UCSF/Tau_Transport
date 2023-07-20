%% Nexis website
clear; clc;
matdir_nexis = '/Users/justintorok/Documents/MATLAB/Nexis_Project/Nexis/raw_data_mouse';
figdirectory = '/Users/justintorok/Documents/MATLAB/CellTypeVulnerability_Project/Figures';
load([matdir_nexis filesep 'Connectomes.mat'],'Connectomes');
load([matdir_nexis filesep 'CCF_labels.mat'],'CCF_labels');
C = Connectomes.default; C = C/max(C);
seedind = 30+213;
seed = zeros(426,1); seed(seedind) = 0.1;
tpts = linspace(0,3,90);
datinput = eNDM_general_dir(seed,tpts,C,zeros(426,1),1,0.5,1,0,0,0,'analytic',1);

brainframedir = '/Users/justintorok/Documents/MATLAB/Brainframe-Dev/Brainframe';
addpath(brainframedir)
savenclose = 0;
reggroups = zeros(213,1); %Chunk of code to define region_groups
amy = 1:11; cer = 12:23; sub = 24:26; hip = 27:37; hyp = 38:57;
ncx = 58:95; med = 96:120; mid = 121:141; olf = 142:149; pal = 150:157;
pon = 158:170; str = 171:178; tha = 179:213;
reggroups(amy) = 1; reggroups(cer) = 2; reggroups(sub) = 3; 
reggroups(hip) = 4; reggroups(hyp) = 5; reggroups(ncx) = 6;
reggroups(med) = 7; reggroups(mid) = 8; reggroups(olf) = 9;
reggroups(pal) = 10; reggroups(pon) = 11; reggroups(str) = 12;
reggroups(tha) = 13;
reggroups = [reggroups;reggroups];
cmap = hsv(length(unique(reggroups))); %Creating colormap

nframes = length(tpts);
totangle = 360;
inc = totangle/nframes;
inc1 = 0; inc2 = 0;
incs1 = zeros(1,length(nframes)); incs2 = incs1;

for i = 1:nframes
    if i < (nframes/4 + 1)
        inc1 = inc1 + inc; inc2 = inc2 + inc/6;
    elseif i < (nframes/2 + 1)
        inc1 = inc1 + inc; inc2 = inc2 - inc/6;
    elseif i < (3*nframes/4 + 1)
        inc1 = inc1 + inc; inc2 = inc2 + inc/6;
    else
        inc1 = inc1 + inc; inc2 = inc2 - inc/6;
    end
    incs1(i) = inc1; incs2(i) = inc2;
end

v = VideoWriter([figdirectory filesep 'testpath.avi']);
open(v);
for i = 1:length(tpts)
    i
    datinput_i = datinput(:,i);
    nany = isnan(datinput_i);
    datinput_i(nany) = 0;
    xfac = 0.25*sum(datinput_i)/sum(datinput(~nany,1));
    imglab = sprintf('website_nexis_spacefilling');
    input_struct_data = brainframe_inputs_mouse(brainframedir,'data',datinput_i,...
                                                'voxUreg',1,...
                                                'xfac',xfac,...
                                                'pointsize',1,...
                                                'norm_method','mean',...
                                                'bgcolor','w',...
                                                'img_format','tiffn',...
                                                'cmap',cmap,...
                                                'region_groups',reggroups,...
                                                'centered',[0 1],...
                                                'regsUbins',0,...
                                                'img_directory',figdirectory,...
                                                'img_labels',imglab,...
                                                'img_renderer','painters',...
                                                'savenclose',savenclose);
    brainframe(input_struct_data);
    view([90 90])
    camorbit(incs1(i),incs2(i),'data',[1 0 0])
    frame = getframe(gcf);
    writeVideo(v,frame);
    close;
end
close(v);

%% CT Website
clear; clc;
brainframedir = '/Users/justintorok/Documents/MATLAB/Brainframe-Dev/Brainframe';
figdirectory = '/Users/justintorok/Documents/MATLAB/CellTypeVulnerability_Project/Figures';
datapath = '/Users/justintorok/Documents/MATLAB/CellTypeVulnerability_Project/Large_MatFiles';
addpath(brainframedir)
load([datapath filesep 'Yao_Dependencies.mat']);
subclasses = classkey;
type = 'DG';
typeinds = find(ismember(subclasses,type));
savenclose = 0;
% Define colormap based on cell class
nonneuronal_inds = ismember(subclasses,{'Endo','Astro','Oligo',...
    'Micro-PVM','SMC-Peri','VLMC'});
gaba_inds = ismember(subclasses,{'Lamp5','Sncg','Meis2','CR',...
    'Pvalb','Sst','Sst Chodl','Vip'});
ctxtest = @(x) strcmp(x(end),'X');
glutctx_inds = logical(cell2mat(cellfun(ctxtest,subclasses,'UniformOutput',false)));
gluthipp_inds = ~logical(nonneuronal_inds + gaba_inds + glutctx_inds);
indcell = {glutctx_inds,gluthipp_inds,gaba_inds,nonneuronal_inds};
indtest = glutctx_inds + 2*gluthipp_inds + 3*gaba_inds + 4*nonneuronal_inds;
cmap_col = [[0 0.75 1]; [1 0 0.5]];
cmap = twocolor(cmap_col(1,:),cmap_col(2,:),length(indcell));
voxthreshes = 0.6; % Brainframe parameter

for i = 1:length(typeinds)
    ctdata = outstruct.corrB(:,typeinds(i));
    newVoxMap = zeros(size(GENGDmod));
    newVoxMap(nonzerovox) = ctdata;
    datinput = imresize3(newVoxMap,[133 81 115]);
    datinput(datinput < 0) = 0;
    col_base = cmap(indtest(typeinds(i)),:);
    col_min = (col_base + 1)/2;
    nbin = 10;
    cmap_i = twocolor(col_min,col_base,nbin);
    input_struct = brainframe_inputs_mouse(brainframedir,'data',datinput,...
                                                'voxthresh',voxthreshes(i),...
                                                'nbin',nbin,...
                                                'voxUreg',0,...
                                                'xfac',0.02,...
                                                'pointsize',0.1,...
                                                'bgcolor','w',...
                                                'img_format','tiffn',...
                                                'cmap',cmap_i,...
                                                'regsUbins',0,...
                                                'img_directory',figdirectory,...
                                                'img_labels',['Yao_' subclasses{typeinds(i)}],...
                                                'img_renderer','painters',...
                                                'savenclose',savenclose);
    brainframe(input_struct);
end

nframes = 180;
totangle = 360;
inc = totangle/nframes;
inc1 = 0; inc2 = 0;
incs1 = zeros(1,length(nframes)); incs2 = incs1;
for i = 1:nframes
    if i < (nframes/4 + 1)
        inc1 = inc1 + inc; inc2 = inc2 + inc/3;
    elseif i < (nframes/2 + 1)
        inc1 = inc1 + inc; inc2 = inc2 - inc/3;
    elseif i < (3*nframes/4 + 1)
        inc1 = inc1 + inc; inc2 = inc2 + inc/3;
    else
        inc1 = inc1 + inc; inc2 = inc2 - inc/3;
    end
    incs1(i) = inc1; incs2(i) = inc2;
end

v = VideoWriter([figdirectory filesep 'testct.avi']);
open(v);
for i = 1:nframes
    view([90 90])
    camorbit(incs1(i),incs2(i),'data',[1 0 0]);
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);

% view([90 90])
% nframes = 180;
% totangle = 360;
% inc = totangle/nframes;
% v = VideoWriter([figdirectory filesep 'testct.avi']);
% open(v);
% for i = 1:nframes
%     if i < (nframes/4 + 1)
%         camorbit(inc,inc/2.5,'data',[1 0 0])
%     elseif i < (nframes/2 + 1)
%         camorbit(inc,-inc/2.5,'data',[1 0 0])
%     elseif i < (3*nframes/4 + 1)
%         camorbit(inc,inc/2.5,'data',[1 0 0])
%     else
%         camorbit(inc,-inc/2.5,'data',[1 0 0])
%     end
%     frame = getframe(gcf);
%     writeVideo(v,frame);
% end
% close(v)