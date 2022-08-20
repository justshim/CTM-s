clearvars
close all
clc

delete(gcp('nocreate'));
parpool('threads');
disp('==============================')
disp('-- Model extraction ')

load_path=strcat(pwd,'/data models/4 vars\mdl_ide_4_integral.mat');
aaa = load(load_path, '*');
mdl_ide_integral = aaa.mdl_ide;

load_path=strcat(pwd,'/data models/4 vars\mdl_ide_4_pi.mat');
aaa = load(load_path, '*');
mdl_ide_pi = aaa.mdl_ide;
clear aaa

ModelTerms_ide_integral = mdl_ide_integral.ModelTerms;
Coefficients_ide_integral = mdl_ide_integral.Coefficients;

ModelTerms_ide_pi = mdl_ide_pi.ModelTerms;
Coefficients_ide_pi = mdl_ide_pi.Coefficients;

tol = 10^(-2);
bin_index=[];
for k=1:length(Coefficients_ide_integral)
    if(abs(Coefficients_ide_integral(k))<tol)
        bin_index = [bin_index k];
    end
end
ModelTerms_ide_integral(bin_index, :) = [];
Coefficients_ide_integral(bin_index) = [];

tol = 10^(-2);
bin_index=[];
for k=1:length(Coefficients_ide_pi)
    if(abs(Coefficients_ide_pi(k))<tol)
        bin_index = [bin_index k];
    end
end
ModelTerms_ide_pi(bin_index, :) = [];
Coefficients_ide_pi(bin_index) = [];

save_file = [pwd, '/data models/4 vars/model_integral_coeff','.mat'];
save(save_file,'Coefficients_ide_integral')

save_file = [pwd, '/data models/4 vars/model_integral_terms','.mat'];
save(save_file,'ModelTerms_ide_integral')

save_file = [pwd, '/data models/4 vars/model_pi_coeff','.mat'];
save(save_file,'Coefficients_ide_pi')

save_file = [pwd, '/data models/4 vars/model_pi_terms','.mat'];
save(save_file,'ModelTerms_ide_pi')

disp('==============================')


