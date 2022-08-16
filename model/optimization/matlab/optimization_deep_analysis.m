clearvars
close all
clc

delete(gcp('nocreate'));
parpool('threads');
disp('==============================')
disp('-- Optimization analysis ')

path = "H:\Il mio Drive\Tesi magistrale\Python\CTM-s\data\opti_data_i-j-delta-beta-priority.csv";
path_ctm = "H:\Il mio Drive\Tesi magistrale\CTMs-identification\fnc\extracted_data\CTM_param_out_nice.xls";
%path = "C:\A_Tesi\Python\CTM-s\data\opti_data_i-j-delta-beta-priority.csv";
%path_ctm = "C:\A_Tesi\CTMs-identification\fnc\extracted_data\CTM_param_out_nice.xls";
%path = "C:\Users\adria\Documents\Uni\LM II anno\Tesi\python\CTM-s\data\opti_data_i-j-delta-beta-priority.csv";
%path_ctm = "C:\Users\adria\Documents\Uni\LM II anno\Tesi\CTMs-identification\fnc\extracted_data\CTM_param_out_nice.xls";
path=strcat(pwd,'/mdl_ide_9.mat');
aaa = load(path, '*');
mdl_ide_9 = aaa.mdl_ide;

path=strcat(pwd,'/mdl_ide_8.mat');
aaa = load(path, '*');
mdl_ide = aaa.mdl_ide;
clear aaa

bin_index = [];
ModelTerms_ide = mdl_ide.ModelTerms;
Coefficients_ide = mdl_ide.Coefficients;
for k=1:length(Coefficients_ide)
    if(abs(Coefficients_ide(k))<1e-5)
        bin_index = [bin_index k];
    end
end
ModelTerms_ide(bin_index, :) = [];
Coefficients_ide(bin_index) = [];

