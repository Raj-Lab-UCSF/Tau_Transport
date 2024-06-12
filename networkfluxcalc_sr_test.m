%%% Test code for SR NetworkFluxCalculator
%% Dummy parameter values
clear; clc;
rng(0)
y1 = rand(1,1);
y2 = rand(1,1);
y3 = rand(1,1);
y4 = rand(1,1);
y5 = rand(1,1);
y6 = rand(1,1);
y7 = rand(1,1);
y8 = rand(1,1);
% y5 = rand(30,1);
% y6 = rand(30,1);
% y7 = rand(30,1);
% y8 = rand(30,1);

loaddir = '/Users/justintorok/Documents/MATLAB/Tau_Transport/MatFiles';

%% True expressions for given complexity
complexity_val_1 = 7;
f_flux_1 = @(y1_,y2_,y3_,y4_,y5_,y6_,y7_,y8_)...
    y2_ * y5_^2 * y7_;
flux_true_1 = f_flux_1(y1,y2,y3,y4,y5,y6,y7,y8);

complexity_val_2 = 6;
f_flux_2 = @(y1_,y2_,y3_,y4_,y5_,y6_,y7_,y8_)...
    -1.9244e-6 * y2_^2 * y7_ / (-y2_ + y7_);
flux_true_2 = f_flux_2(y1,y2,y3,y4,y5,y6,y7,y8);

%% Load in and parse csv of DSO outputs
dso_filename = [loaddir filesep 'dso_test_expr_output.csv'];
dso_table = readmatrix(dso_filename,'OutputType','char');
complexity_rows = cellfun(@(x)str2double(x),dso_table(:,1),'UniformOutput',false);
complexity_rows = cell2mat(complexity_rows);

row_ind_1 = find(complexity_rows == complexity_val_1);
expression_expr_1 = dso_table{row_ind_1,5};
row_ind_2 = find(complexity_rows == complexity_val_2);
expression_expr_2 = dso_table{row_ind_2,5};

%% 
syms f(x,y)
f(x,y) = x^2 * y;
f(1,2)

%% Parse traversal expression 1 and create anonymous function
exprcell = expression_expr_1(:);
exprfun = @(a2,a5,a7) 1;
i = 1;
while ~isempty(exprcell)
    if strcmp(exprcell(1),'x')
        switch str2double(exprcell(2))
            case 1
            case 2
                exprfun = @(a2,a5,a7) exprfun(a2,a5,a7) * a2;
                exprcell(1:2) = [];
            case 3
            case 4
            case 5
                exprfun = @(a2,a5,a7) exprfun(a2,a5,a7) * a5;
                exprcell(1:2) = [];
            case 6
            case 7
                exprfun = @(a2,a5,a7) exprfun(a2,a5,a7) * a7;
                exprcell(1:2) = [];
            case 8
        end
    elseif strcmp(exprcell(1),'*') && strcmp(exprcell(2),'*')
        exprcell(1) = [];
        exprcell(1) = '^';
    elseif strcmp(exprcell(1),'*') && ~strcmp(exprcell(2),'*')
        exprcell(1) = [];
        if strcmp(exprcell(1),'x')
            switch str2double(exprcell(2))
                case 1
                case 2
                    exprfun = @(a2,a5,a7) exprfun(a2,a5,a7) * a2;
                    exprcell(1:2) = [];
                case 3
                case 4
                case 5
                    exprfun = @(a2,a5,a7) exprfun(a2,a5,a7) * a5;
                    exprcell(1:2) = [];
                case 6
                case 7
                    exprfun = @(a2,a5,a7) exprfun(a2,a5,a7) * a7;
                    exprcell(1:2) = [];
                case 8
            end
        elseif ~isnan(str2double(exprcell(1)))
            exprfun = @(a2,a5,a7) exprfun(a2,a5,a7) * str2double(exprcell(1));
            exprcell(1) = [];
        end
    elseif strcmp(exprcell(1),'^')
        exprcell(1) = [];
        if strcmp(exprcell(1),'x')
            switch str2double(exprcell(2))
                case 1
                case 2
                    exprfun = @(a2,a5,a7) exprfun(a2,a5,a7) ^ a2;
                    exprcell(1:2) = [];
                case 3
                case 4
                case 5
                    exprfun = @(a2,a5,a7) exprfun(a2,a5,a7) ^ a5;
                    exprcell(1:2) = [];
                case 6
                case 7
                    exprfun = @(a2,a5,a7) exprfun(a2,a5,a7) ^ a7;
                    exprcell(1:2) = [];
                case 8
            end
        elseif ~isnan(str2double(exprcell(1)))
            exprfun = @(a2,a5,a7) exprfun(a2,a5,a7) ^ str2double(exprcell(1));
            str2double(exprcell(1))
            exprcell(1) = [];
        end
    end
    i = i + 1;
end
flux_calc_1 = exprfun(y2,y5,y7);
% expr_fun_1 = expression_parser(expression_expr_1);
% flux_expr_1 = expr_fun_1(y1,y2,y3,y4,y5,y6,y7,y8);

%% Parse traversal expression 2 and create anonymous function


%% Functions
function expr_fun = expression_parser(expression_expr_)
    expr_cell = strsplit(expression_expr_,',');
    expr_fun = @(z1_,z2_,z3_,z4_,z5_,z6_,z7_,z8_) 1;
    i = 1;
    while ~isempty(expr_cell)
        i = i + 1;
    end
end


