clearvars
close all
clc

delete(gcp('nocreate'));
parpool('threads');
disp('==============================')
disp('-- Optimization analysis ')

%path = "H:\Il mio Drive\Tesi magistrale\Python\CTM-s\data\opti_data_i-j-delta-beta-priority.csv";
%path = "C:\A_Tesi\Python\CTM-s\data\opti_data_i-j-delta-beta-priority.csv";
path = "C:\Users\adria\Documents\Uni\LM II anno\Tesi\python\CTM-s\data\opti_data_i-j-delta-beta-priority.csv";

read = 1;
generate = 1-read;

shallow = 1;
deep = 1-shallow;

output = single(1); % 1 = integral delta, 2 = pi greco

path_output=strcat(pwd,'\opti-ide-val-single-precision.xlsx');

%% extract data
T = readtable(path);

varname=T.Properties.VariableNames;
A = table2array(T);
cell_in = A(:, 1);
cell_out = A(:, 2);
delta = A(:, 3);
beta = A(:, 4);
priority = A(:,5);
integral=A(:,6);
pi_greco = A(:, 8);

%% Generate dataset
if (generate)
    aa = randperm(length(cell_in),length(cell_in)*0.15)';
    not_aa=setdiff(1:length(cell_in),aa)';
    A_val = A(aa, :);
    A_ide = A(not_aa, :);
    save_file = [pwd, '/A_val','.mat'];
    save(save_file,'A_val')
    save_file = [pwd, '/A_ide','.mat'];
    save(save_file,'A_ide')
end

%% Read dataset
if(read)
    path=strcat(pwd,'/A_val.mat');
    aaa = load(path, '*');
    A_val = aaa.A_val;
    path=strcat(pwd,'/A_ide.mat');
    aaa = load(path, '*');
    A_ide = aaa.A_ide;
    clear pi_greco integral priority beta delta aaa T A cell_out cell_in varname path_ctm

    I_val=[A_val(:,1) A_val(:,2) A_val(:,3) A_val(:,4)];
    I_ide=[A_ide(:,1) A_ide(:,2) A_ide(:,3) A_ide(:,4)];
    O_val = [A_val(:,6) A_val(:,8)];
    O_ide = [A_ide(:,6) A_ide(:,8)];
    clear A_ide A_val
end

%% Find best degree polynomial model
if(shallow)
    first_j=single(1);
    last_j=single(12);

    grado=ones(last_j-first_j+1, 1);
    R2_ide=ones(last_j-first_j+1, 1);
    R2a_ide=ones(last_j-first_j+1, 1);
    RMSE_ide=ones(last_j-first_j+1, 1);
    SSR_ide=ones(last_j-first_j+1, 1);
    Val_FPE=ones(last_j-first_j+1, 1);
    Val_AIC=ones(last_j-first_j+1, 1);
    Val_MDL=ones(last_j-first_j+1, 1);
    R2_val=ones(last_j-first_j+1, 1);
    R2a_val=ones(last_j-first_j+1, 1);
    RMSE_val=ones(last_j-first_j+1, 1);
    SSR_val=ones(last_j-first_j+1, 1);

    header = ["grado", "SSR_ide", "SSR_val", "Val_FPE", "Val_AIC", "Val_MDL", "R2_ide", "R2a_ide", "RMSE_ide", "R2_val", "R2a_val", "RMSE_val"];
    if (output == 1)
        sheet = 'shallow - integral delta';
    elseif (output == 2)
        sheet = 'shallow - pi';
    end
    writematrix(header,path_output,'Sheet',sheet);

    for j=first_j:last_j
        i=j-first_j+1;
        fprintf('... degree %d \n',j)
        grado(i) = j;
        n = length(O_ide(:,output));
        mdl_ide = polyfitn(I_ide, O_ide(:,output), i);
        ypred_ide = polyvaln(mdl_ide,I_ide);
        R2_ide(i) = mdl_ide.R2;
        R2a_ide(i) = mdl_ide.AdjustedR2;
        RMSE_ide(i) = mdl_ide.RMSE;
        q=length(mdl_ide.Coefficients); %number of parameters
        epsilonIde=O_ide(:,output)-ypred_ide; %residui identificazione
        SSR_ide(i)= epsilonIde'*epsilonIde;
        Val_FPE(i)=((n+q)/(n-q))*SSR_ide(i);
        Val_AIC(i)=2*q/n+log(SSR_ide(i));
        Val_MDL(i)=log(n)*q/n + log(SSR_ide(i));

        %% cross validation

        mdl_val = polyfitn(I_val, O_val(:,output), i);
        ypred_val = polyvaln(mdl_val,I_val);
        R2_val(i) = mdl_val.R2;
        R2a_val(i) = mdl_val.AdjustedR2;
        RMSE_val(i) = mdl_val.RMSE;

        epsilonVal=O_val(:,output)-ypred_val; %residui validazione
        SSR_val(i)= epsilonVal'*epsilonVal;

    end
    fprintf('Saving information in %s ...\n',path)
    tabella = [grado,SSR_ide,SSR_val,Val_FPE,Val_AIC,Val_MDL,R2_ide, R2a_ide, RMSE_ide, R2_val, R2a_val, RMSE_val];
    writematrix(tabella, path_output,'Sheet',sheet,'WriteMode','append');
    disp('==============================')
end

if(deep)
    header = ["grado", "tol", "SSR_ide", "SSR_val", "Val_FPE", "Val_AIC", "Val_MDL", "R2_ide", "R2a_ide", "RMSE_ide", "R2_val", "R2a_val", "RMSE_val"];
    if (output == 1)
        sheet = 'deep - integral delta';
    elseif (output == 2)
        sheet = 'deep - pi';
    end
    writematrix(header, path_output, 'Sheet', sheet);
    model = strcat('/mdl_ide_', '9', '.mat');
    load_path=strcat(pwd,model);
    aaa = load(load_path, '*');
    mdl_ide = aaa.mdl_ide;
    clear aaa

    for ijk=1:10
        tol=10^(-ijk);
        bin_index = [];

        ModelTerms_ide = mdl_ide.ModelTerms;
        Coefficients_ide = mdl_ide.Coefficients;
        for k=1:length(Coefficients_ide)
            if(abs(Coefficients_ide(k))<tol)
                bin_index = [bin_index k];
            end
        end
        ModelTerms_ide(bin_index, :) = [];
        Coefficients_ide(bin_index) = [];

        fprintf('... tolerance %d \n',tol)
        grado = 9;
        n = length(O_ide(:,output));
        mdl_ide = polyfitn(I_ide, O_ide(:,output), ModelTerms_ide);
        ypred_ide = polyvaln(mdl_ide,I_ide);
        R2_ide = mdl_ide.R2;
        R2a_ide = mdl_ide.AdjustedR2;
        RMSE_ide = mdl_ide.RMSE;
        q=length(mdl_ide.Coefficients); %number of parameters
        epsilonIde=O_ide(:,output)-ypred_ide; %residui identificazione
        SSR_ide= epsilonIde'*epsilonIde;
        Val_FPE=((n+q)/(n-q))*SSR_ide;
        Val_AIC=2*q/n+log(SSR_ide);
        Val_MDL=log(n)*q/n + log(SSR_ide);

        %% cross validation

        mdl_val = polyfitn(I_val, O_val(:,output), ModelTerms_ide);
        ypred_val = polyvaln(mdl_val,I_val);
        R2_val = mdl_val.R2;
        R2a_val = mdl_val.AdjustedR2;
        RMSE_val = mdl_val.RMSE;

        epsilonVal=O_val(:,output)-ypred_val; %residui validazione
        SSR_val= epsilonVal'*epsilonVal;
        
        fprintf('Saving information in %s ...\n',path_output)
        tabella = [grado,tol,SSR_ide,SSR_val,Val_FPE,Val_AIC,Val_MDL,R2_ide, R2a_ide, RMSE_ide, R2_val, R2a_val, RMSE_val];
        writematrix(tabella, path_output,'Sheet',sheet,'WriteMode','append');
    end   
end

