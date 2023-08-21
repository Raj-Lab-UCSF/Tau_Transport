% Wrapper script for generating output files for NetworkTransportModel for
% generating tau distributions on the network

%% 1. Define directories for saving outputs
clear; clc;
curpath = '/wynton/protected/home/rajlab/jtorok/MATLAB/Tau_Transport';
p = genpath(curpath);
addpath(p);
simpath = [curpath filesep 'SampleFiles'];
loadpath = [curpath filesep 'MatFiles'];
if ~isfolder(simpath)
    mkdir(simpath)
end
simstr = 'hippocampome_final_round2_3'; % for saving the outputs

%% 2. Parameter definitions
% 2a. Define actively tuned parameters as scalars or arrays to be explored
% on a grid search
inputparams = cell(2,8);
paramnames = {'beta','gamma1','gamma2','frac','lambda1','lambda2',...
    'delta','epsilon'};
inputparams(1,:) = paramnames;
inputparams{2,1} = 1e-6; % beta
inputparams{2,2} = [4e-3,8e-3]; % gamma1
inputparams{2,3} = 0; % gamma2
inputparams{2,4} = 0.92; % frac
inputparams{2,5} = 0.05; % lambda1
inputparams{2,6} = 0.05; % lambda2
inputparams{2,7} = [10,25,50,100]; % delta
inputparams{2,8} = [10,25,50,100]; % epsilon

% 2b. Create parameter array to grid search using allcomb()
paramgrid = allcomb(inputparams{2,1},...
                    inputparams{2,2},...
                    inputparams{2,3},...
                    inputparams{2,4},...
                    inputparams{2,5},...
                    inputparams{2,6},...
                    inputparams{2,7},...
                    inputparams{2,8});
paramnamescell = repmat(paramnames,size(paramgrid,1),1);

% 2c. Define other parameters
L_int = 1000; % default = 1000; in microns
L1 = 200; % default = 200
L2 = 200; % default = 200
L_ais = 40; % default = 40
L_syn = 40; % default = 40
T = []; % default = 0.05
dt = []; % default = 0.005
trange = [0:0.0025:0.1, 0.105:0.005:0.3, 0.31:0.01:1];
resmesh = 'coarse'; % 'fine' or 'coarse' - use 'coarse' for faster, less precise simulations
plotting = 0;
reltol = 1e-4;
abstol = 1e-4;
fsolvetol = 1e-6;
init_rescale = 0.02;
init_path = {'Entorhinal area, lateral part_L'};
study = 'Hurtado';
connectome_subset = 'Hippocampus+PC+RSP';
ncores = 32;

%% 3. Run NetworkTransportModel
output_struct = struct;
output_struct.Parameter_Grid = paramgrid;   
output_struct.Parameter_Names = inputparams(1,:);
sim_struct = struct;
parpool(ncores)
tic
% for i = 1:size(paramgrid,1)
parfor i = 1:size(paramgrid,1)
    fprintf('Simulation %d/%d \n',i,size(paramgrid,1))
    paramlist = paramgrid(i,:);
    paramnames_i = paramnamescell(i,:); % prevents broadcast warning message
    mdloutput = NetworkTransportModel(loadpath,paramnames_i{1},paramlist(1),...
                                paramnames_i{2},paramlist(2),...
                                paramnames_i{3},paramlist(3),...
                                paramnames_i{4},paramlist(4),...
                                paramnames_i{5},paramlist(5),...
                                paramnames_i{6},paramlist(6),...
                                paramnames_i{7},paramlist(7),...
                                paramnames_i{8},paramlist(8),...
                                'L_int',L_int,...
                                'L1',L1,...
                                'L2',L2,...
                                'L_ais',L_ais,...
                                'L_syn',L_syn,...
                                'T',T,...
                                'dt',dt,...
                                'trange',trange,...
                                'resmesh', resmesh,...
                                'plotting',plotting,...
                                'reltol',reltol,...
                                'abstol',abstol,...
                                'fsolvetol',fsolvetol,...
                                'init_rescale',init_rescale,...
                                'init_path',init_path,...
                                'study',study,...
                                'connectome_subset',connectome_subset,...
                                'sim_no',i);
    sim_struct(i).Model_Outputs = mdloutput;
end
output_struct.Simulations = sim_struct;
delete(gcp('nocreate'));
toc

% testplot

%% 4. Save output file
save([simpath filesep simstr '.mat'],'output_struct') 
clear