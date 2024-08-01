% Wrapper script for generating output files for NetworkTransportModel for
% generating tau distributions on the network

%% 1. Define directories for saving outputs
clear; clc;
curpath = '/wynton/protected/home/rajlab/jtorok/MATLAB/Tau_Transport'; % CHANGE THIS LINE TO WHERE Tau_Transport DIRECTORY IS
p = genpath(curpath);
addpath(p);
simpath = [curpath filesep 'SampleFiles']; % THIS IS WHERE THE OUTPUTS WILL BE SAVED
loadpath = [curpath filesep 'MatFiles'];
if ~isfolder(simpath)
    mkdir(simpath)
end
simstr = 'All_CA1seed_dir_Brandon'; % CHANGE THIS AS NEEDED, THIS IS THE NAME OF THE OUTPUT FILE

%% 2. Parameter definitions
% 2a. Define actively tuned parameters as scalars or arrays to be explored
% on a grid search
inputparams = cell(2,8);
paramnames = {'beta','gamma1','gamma2','frac','lambda1','lambda2',...
    'delta','epsilon'};
inputparams(1,:) = paramnames;

% CAN CHANGE THESE VALUES (LINES 25-32), BUT KEEP VALUES CLOSE. ALL CAN BE
% ARRAYS
inputparams{2,1} = 1e-6; % beta, KEEP SAME
inputparams{2,2} = [1e-3,2e-3,4e-3,8e-3]; % gamma1
inputparams{2,3} = 0; % gamma2, KEEP SAME
inputparams{2,4} = 0.92; % frac, KEEP SAME (though this one will be interesting at some point, maybe)
inputparams{2,5} = 0.05; % lambda1
inputparams{2,6} = 0.05; % lambda2, KEEP SAME AS lambda1
inputparams{2,7} = [1,100]; % delta
inputparams{2,8} = [1,100]; % epsilon

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

% DO NOT CHANGE THESE EXCEPT WHERE INDICATED!!!
L_int = 1000; % default = 1000; in microns
L1 = 200; % default = 200
L2 = 200; % default = 200
L_ais = 40; % default = 40
L_syn = 40; % default = 40
T = []; % default = 0.05
dt = []; % default = 0.005
trange = [0:0.001:0.05, 0.055:0.005:0.3, 0.31:0.01:1];
resmesh = 'coarse';
plotting = 0;
reltol = 1e-4;
abstol = 1e-4;
fsolvetol = 1e-6;
init_rescale = 0.2;
init_path = {'Field CA1_L'};
study = 'DS9';
connectome_subset = 'LH'; % YOU SHOULDN'T HAVE TO CHANGE THIS, BUT FOR WHOLE BRAIN IT WOULD BE 'WB'
conn_thresh = 0.8;
ncores = 16; % NUMBER OF CORES TO USE FOR PARFOR, IF RUNNING IN PARALLEL

%% 3. Run NetworkTransportModel
output_struct = struct;
output_struct.Parameter_Grid = paramgrid;   
output_struct.Parameter_Names = inputparams(1,:);
sim_struct = struct;
parpool(ncores) % IF DOING SERIALLY, COMMENT OUT
tic
% for i = 1:size(paramgrid,1) % IF DOING SERIALLY, UNCOMMENT
parfor i = 1:size(paramgrid,1) % IF DOING SERIALLY, COMMENT OUT
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
                                'sim_no',i,...
                                'conn_thresh',conn_thresh);
    sim_struct(i).Model_Outputs = mdloutput;
end
output_struct.Simulations = sim_struct;
delete(gcp('nocreate')); % IF DOING SERIALLY, COMMENT OUT
toc

%% 4. Save output file
save([simpath filesep simstr '.mat'],'output_struct') 
clear