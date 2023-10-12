output_struct_final = struct;
simstr_base = 'hippocampome_final_round2_';
simpath = '~/Documents/MATLAB/Tau_Transport/SampleFiles';
output_struct_final.Parameter_Names = [];
output_struct_final.Parameter_Grid = [];
output_struct_final.Simulations = [];
for i = 1:7
    load([simpath filesep simstr_base num2str(i) '.mat'],'output_struct');
    output_struct_final.Parameter_Names = output_struct.Parameter_Names;
    output_struct_final.Parameter_Grid = [output_struct_final.Parameter_Grid; output_struct.Parameter_Grid];
    if i ~= 1
        output_struct_final.Simulations = [output_struct_final.Simulations, output_struct.Simulations];
    else
        output_struct_final.Simulations = output_struct.Simulations;
    end
    clear output_struct
end