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

addpath(path)
addpath(path_ctm)
warning('off')

poly_fit = 0;
plots = 0;

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
% aa = randperm(length(cell_in),length(cell_in)*0.15)';
% not_aa=setdiff(1:length(cell_in),aa)';
% A_val = A(aa, :);
% A_ide = A(not_aa, :);
% save_file = [pwd, '/A_val','.mat'];
% save(save_file,'A_val')
% save_file = [pwd, '/A_ide','.mat'];
% save(save_file,'A_ide')

%% Read dataset
path=strcat(pwd,'/A_val.mat');
aaa = load(path, '*');
A_val = aaa.A_val;
path=strcat(pwd,'/A_ide.mat');
aaa = load(path, '*');
A_ide = aaa.A_ide;
clear pi_greco integral priority beta delta aaa T A cell_out cell_in varname path_ctm
path=strcat(pwd,'\opti-ide-val-single-precision.xlsx');
header = ["grado", "tol", "SSR_ide", "SSR_val", "Val_FPE", "Val_AIC", "Val_MDL", "R2_ide", "R2a_ide", "RMSE_ide", "R2_val", "R2a_val", "RMSE_val"];
writematrix(header,path, 'Sheet','deep - integral delta');

I_val=single([A_val(:,1) A_val(:,2) A_val(:,3) A_val(:,4)]);
I_ide=single([A_ide(:,1) A_ide(:,2) A_ide(:,3) A_ide(:,4)]);
O_val = single([A_val(:,6) A_val(:,8)]);
O_ide = single([A_ide(:,6) A_ide(:,8)]);
clear A_ide A_val

first_j=int8(8);
last_j=int8(8);
output = int8(1); % 1 = integral delta, 2 = pi greco

grado=ones(last_j-first_j+1, 1, 'single');
R2_ide=ones(last_j-first_j+1, 1, 'single');
R2a_ide=ones(last_j-first_j+1, 1, 'single');
RMSE_ide=ones(last_j-first_j+1, 1, 'single');
SSR_ide=ones(last_j-first_j+1, 1, 'single');
Val_FPE=ones(last_j-first_j+1, 1, 'single');
Val_AIC=ones(last_j-first_j+1, 1, 'single');
Val_MDL=ones(last_j-first_j+1, 1, 'single');
R2_val=ones(last_j-first_j+1, 1, 'single');
R2a_val=ones(last_j-first_j+1, 1, 'single');
RMSE_val=ones(last_j-first_j+1, 1, 'single');
SSR_val=ones(last_j-first_j+1, 1, 'single');

%% zozzerie
load_path=strcat(pwd,'/mdl_ide_8.mat');
aaa = load(load_path, '*');
mdl_ide_8 = aaa.mdl_ide;
clear aaa

for ijk=1:10
    tol=10^ijk;
    bin_index = [];

    ModelTerms_ide = mdl_ide_8.ModelTerms;
    Coefficients_ide = mdl_ide_8.Coefficients;
    for k=1:length(Coefficients_ide)
        if(abs(Coefficients_ide(k))<tol)
            bin_index = [bin_index k];
        end
    end
    ModelTerms_ide(bin_index, :) = [];
    Coefficients_ide(bin_index) = [];

    for j=first_j:last_j
        i=j-first_j+1;
        fprintf('... degree %d \n',j)
        grado(i) = j;
        n = length(O_ide(:,output));
        mdl_ide = polyfitn(I_ide, O_ide(:,output), ModelTerms_ide);
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

        mdl_val = polyfitn(I_val, O_val(:,output), ModelTerms_ide);
        ypred_val = polyvaln(mdl_val,I_val);
        R2_val(i) = mdl_val.R2;
        R2a_val(i) = mdl_val.AdjustedR2;
        RMSE_val(i) = mdl_val.RMSE;

        epsilonVal=O_val(:,output)-ypred_val; %residui validazione
        SSR_val(i)= epsilonVal'*epsilonVal;

    end
    fprintf('Saving information in %s ...\n',path)
    tabella = [grado,tol,SSR_ide,SSR_val,Val_FPE,Val_AIC,Val_MDL,R2_ide, R2a_ide, RMSE_ide, R2_val, R2a_val, RMSE_val];
    writematrix(tabella, path,'Sheet','deep - integral delta','WriteMode','append');
end
%% save file
fprintf('Saving information in %s ...\n',path)
tabella = [grado,tol,SSR_ide,SSR_val,Val_FPE,Val_AIC,Val_MDL,R2_ide, R2a_ide, RMSE_ide, R2_val, R2a_val, RMSE_val];
writematrix(tabella, path,'Sheet','deep - integral delta','WriteMode','append');
disp('==============================')



if(poly_fit)
    T = readtable(path_ctm);
    A = table2array(T);
    q_max = A(:, 5);
    rho_max = A(:, 6);
    fit_name=["poly11","poly12","poly21","poly22","poly13","poly31","poly23","poly32","poly33", ...
        "poly14","poly41","poly24","poly42","poly34","poly43","poly44",...
        "poly15","poly51","poly25","poly52","poly35","poly53","poly45","poly54","poly55"];

    %% 1 - Integral with different i and j, delta fixed
    disp("1 - Integral with different i and j, delta fixed")
    fixed_delta = 480;
    x = [];
    y = [];
    z = [];
    parfor i=1:length(cell_in)
        if(delta(i)==fixed_delta)
            x=[x cell_in(i)]; %i
            y=[y cell_out(i)]; %j
            z=[z integral(i)]; %integral
        end
    end
    x=x';
    y=y';
    z=z';

    [f, gof] = fit([x, y],z,fit_name(1),'Normalize', 'on', 'Robust','Bisquare');
    f_best=f;
    gof_best=gof;
    adjrsq_best = gof_best.adjrsquare;
    sse_best = gof_best.sse;
    for i=2:length(fit_name)
        [f, gof] = fit([x, y],z,fit_name(i), 'Normalize', 'on', 'Robust','Bisquare');
        if(gof.adjrsquare > adjrsq_best && gof.sse < sse_best)
            f_best=f;
            gof_best=gof;
            adjrsq_best = gof_best.adjrsquare;
            sse_best = gof_best.sse;
        end
    end

    if(plots>0)
        figure()
        plot(f_best,[x y],z, 'Exclude', z < 100)
        zlim([0 max(z)]);
        xlabel('x - Cell IN')
        ylabel('y - Cell OUT')
        zlabel("Integral")
        title("Integral VS i, j")
        legend(type(f_best))
        grid on
    end

    %% 1.1 - Integral with different rho out, delta and beta fixed
    disp("1.1 - Integral with different rho in and rho out, delta fixed")
    fixed_delta = 480;
    fixed_priority = 0.05;
    fixed_beta = 0.07;
    y = [];
    z = [];
    parfor i=1:length(cell_in)
        if(delta(i)==fixed_delta && beta(i)==fixed_beta)
            y=[y cell_out(i)]; %j
            z=[z integral(i)]; %integral
        end
    end
    y=y';
    z=z';

    for i=1:length(y)
        for k=1:length(rho_max)
            if(y(i)==k)
                y(i)=rho_max(k); %j -> rho max
            end
        end
    end
    y = sort(y);

    if(plots>0)
        figure()
        scatter(y,z)
        hold on
        plot(fit(y,z,'poly3'))
        xlabel('x - rho max cell OUT')
        ylabel('Integral')
        title("Integral VS rho max j")
        grid on
    end

    %% 1.2 - Integral with different q_max out, delta and beta fixed
    disp("1.2 - Integral with different q_max out, delta and beta fixed")
    fixed_delta = 480;
    fixed_priority = 0.05;
    fixed_beta = 0.07;
    y = [];
    z = [];
    parfor i=1:length(cell_in)
        if(delta(i)==fixed_delta && beta(i)==fixed_beta)
            y=[y cell_out(i)]; %j
            z=[z integral(i)]; %integral
        end
    end
    y=y';
    z=z';

    for i=1:length(y)
        for k=1:length(q_max)
            if(y(i)==k)
                y(i)=q_max(k); %j -> rho max
            end
        end
    end
    y = sort(y);

    if(plots>0)
        figure()
        scatter(y,z)
        hold on
        plot(fit(y,z,'poly3'))
        xlabel('x - q max cell OUT')
        ylabel('Integral')
        title("Integral VS q max j")
        grid on
    end

    %% 2 - Integral with different i and delta(10-120 min), j fixed
    disp('2 - Integral with different i and delta(10-120 min), j fixed')
    j = 2;
    x = [];
    y = [];
    z = [];
    parfor i=1:length(cell_in)
        if(cell_out(i)==cell_in(i)+j)
            x=[x cell_in(i)]; %i
            y=[y delta(i)]; %delta
            z=[z integral(i)]; %integral
        end
    end
    x=x';
    y=y';
    z=z';
    [f, gof] = fit([x, y],z,fit_name(1),'Normalize', 'on', 'Robust','Bisquare');
    f_best=f;
    gof_best=gof;
    adjrsq_best = gof_best.adjrsquare;
    sse_best = gof_best.sse;
    for i=2:length(fit_name)
        [f, gof] = fit([x, y],z,fit_name(i), 'Normalize', 'on', 'Robust','Bisquare');
        if(gof.adjrsquare > adjrsq_best && gof.sse < sse_best)
            f_best=f;
            gof_best=gof;
            adjrsq_best = gof_best.adjrsquare;
            sse_best = gof_best.sse;
        end
    end

    if(plots>0)
        figure()
        plot(f_best,[x y],z,'Exclude', z < 100)
        zlim([0 max(z)]);
        xlabel('x - Cell IN')
        ylabel('y - Delta (timesteps)')
        yticks(0:min(y):max(y))
        zlabel("Integral")
        title("Integral VS i, delta")
        legend(type(f_best))
        grid on
    end


    %% 3 - Pi with different i and delta(10-120 min), j fixed
    disp('3 - Pi with different i and delta(10-120 min), j fixed')
    j = 2;
    x = [];
    y = [];
    z = [];
    parfor i=1:length(cell_in)
        if(cell_out(i)==cell_in(i)+j)
            x=[x cell_in(i)]; %i
            y=[y delta(i)]; %delta
            z=[z pi_greco(i)]; %pi
        end
    end
    x=x';
    y=y';
    z=z';
    [f, gof] = fit([x, y],z,fit_name(1),'Exclude', z < 0);
    f_best=f;
    gof_best=gof;
    adjrsq_best = gof_best.adjrsquare;
    sse_best = gof_best.sse;
    for i=2:length(fit_name)
        [f, gof] = fit([x, y],z,fit_name(i), 'Exclude', z < 0);
        if(gof.adjrsquare > adjrsq_best && gof.sse < sse_best)
            f_best=f;
            gof_best=gof;
            adjrsq_best = gof_best.adjrsquare;
            sse_best = gof_best.sse;
        end
    end

    if(plots>0)
        figure()
        plot(f_best,[x y],z, 'Exclude', z < 0)
        zlim([-1 1]);
        xlabel('x - Cell IN')
        ylabel('y - Delta (timesteps)')
        yticks(0:min(y):max(y))
        zlabel("PI")
        title("PI VS i, delta")
        legend(type(f_best))
        grid on
    end

    %% 3.1 - Pi with different rho out, delta and beta fixed
    disp("3.1 - Pi with different rho in and rho out, delta fixed")
    fixed_delta = 480;
    fixed_priority = 0.05;
    fixed_beta = 0.07;
    y = [];
    z = [];
    parfor i=1:length(cell_in)
        if(delta(i)==fixed_delta && beta(i)==fixed_beta)
            y=[y cell_out(i)]; %j
            z=[z pi_greco(i)]; %integral
        end
    end
    y=y';
    z=z';

    for i=1:length(y)
        for k=1:length(rho_max)
            if(y(i)==k)
                y(i)=rho_max(k); %j -> rho max
            end
        end
    end
    y = sort(y);

    if(plots>0)
        figure()
        scatter(y,z)
        hold on
        plot(fit(y,z,'poly3'))
        xlabel('x - rho max cell OUT')
        ylabel('Pi')
        title("Pi VS rho max j")
        grid on
    end

    %% 3.2 - Pi with different q_max out, delta and beta fixed
    disp("3.2 - Pi with different q_max out, delta and beta fixed")
    fixed_delta = 480;
    fixed_priority = 0.05;
    fixed_beta = 0.07;
    y = [];
    z = [];
    parfor i=1:length(cell_in)
        if(delta(i)==fixed_delta && beta(i)==fixed_beta)
            y=[y cell_out(i)]; %j
            z=[z pi_greco(i)]; %integral
        end
    end
    y=y';
    z=z';

    for i=1:length(y)
        for k=1:length(q_max)
            if(y(i)==k)
                y(i)=q_max(k); %j -> rho max
            end
        end
    end
    y = sort(y);

    if(plots>0)
        figure()
        scatter(y,z)
        hold on
        plot(fit(y,z,'poly3'))
        xlabel('x - q max cell OUT')
        ylabel('Pi')
        title("Pi VS q max j")
        grid on
    end

    %% 4 - Integral with different i and beta, j fixed
    disp('4 - Integral with different i and beta, j fixed')
    j = 2;
    x = [];
    y = [];
    z = [];

    parfor i=1:length(cell_in)
        if(cell_out(i)==cell_in(i)+j)
            x=[x cell_in(i)]; %i
            y=[y beta(i)]; %beta
            z=[z integral(i)]; %integral
        end
    end
    x=x';
    y=y';
    z=z';
    [f, gof] = fit([x, y],z,fit_name(1));
    f_best=f;
    gof_best=gof;
    adjrsq_best = gof_best.adjrsquare;
    sse_best = gof_best.sse;
    for i=2:length(fit_name)
        [f, gof] = fit([x, y],z,fit_name(i));
        if(gof.adjrsquare > adjrsq_best && gof.sse < sse_best)
            f_best=f;
            gof_best=gof;
            adjrsq_best = gof_best.adjrsquare;
            sse_best = gof_best.sse;
        end
    end

    if(plots>0)
        figure()
        plot(f_best,[x y],z, 'Exclude', z < 100)
        zlim([0 max(z)]);
        xlabel('x - Cell IN')
        ylabel('y - Beta')
        zlabel("Integral")
        title("Integral VS i, beta")
        legend(type(f_best))
        grid on
    end
    %% 5 - Pi with different i and beta, j fixed
    disp('5 - Pi with different i and beta, j fixed')
    j = 3;
    x = [];
    y = [];
    z = [];
    parfor i=1:length(cell_in)
        if(cell_out(i)==cell_in(i)+j)
            x=[x cell_in(i)]; %i
            y=[y beta(i)]; %beta
            z=[z pi_greco(i)]; %pi
        end
    end
    x=x';
    y=y';
    z=z';
    [f, gof] = fit([x, y],z,fit_name(1));
    f_best=f;
    gof_best=gof;
    adjrsq_best = gof_best.adjrsquare;
    sse_best = gof_best.sse;
    for i=2:length(fit_name)
        [f, gof] = fit([x, y],z,fit_name(i));
        if(gof.adjrsquare > adjrsq_best && gof.sse < sse_best)
            f_best=f;
            gof_best=gof;
            adjrsq_best = gof_best.adjrsquare;
            sse_best = gof_best.sse;
        end
    end


    if(plots>0)
        figure()
        plot(f_best,[x y],z, 'Exclude', z < 0)
        zlim([-1 1]);
        xlabel('x - Cell IN')
        ylabel('y - Beta')
        zlabel("PI")
        title("PI VS i, beta")
        legend(type(f_best))
        grid on
    end

    % %% 6 - Pi with different delta and beta
    % disp('6 - Pi with different delta and beta')
    %
    % x = [];
    % y = [];
    % z = [];
    % parfor i=1:length(cell_in)
    %     x=[x delta(i)]; %delta
    %     y=[y beta(i)]; %beta
    %     z=[z pi_greco(i)]; %pi
    % end
    % x=x';
    % y=y';
    % z=z';
    % [f, gof] = fit([x, y],z,fit_name(1),'Normalize', 'on', 'Robust','Bisquare');
    % f_best=f;
    % gof_best=gof;
    % adjrsq_best = gof_best.adjrsquare;
    % sse_best = gof_best.sse;
    % for i=2:length(fit_name)
    %     [f, gof] = fit([x, y],z,fit_name(i), 'Normalize', 'on', 'Robust','Bisquare');
    %     if(gof.adjrsquare > adjrsq_best && gof.sse < sse_best)
    %         f_best=f;
    %         gof_best=gof;
    %         adjrsq_best = gof_best.adjrsquare;
    %         sse_best = gof_best.sse;
    %     end
    % end
    %
    % if(plots>0)
    %     figure()
    %     plot(f_best,[x y],z, 'Exclude', z < 0)
    %     zlim([-1 1]);
    %     xlabel('x - Delta')
    %     xticks(0:min(x):max(x))
    %     ylabel('y - Beta')
    %     zlabel("PI")
    %     title("PI VS Delta, beta")
    %     legend(type(f_best))
    %     grid on
    % end
    % %% 7 - Integral with different delta and beta
    % disp('7 - Integral with different delta and beta ')
    %
    % x = [];
    % y = [];
    % z = [];
    % parfor i=1:length(cell_in)
    %     x=[x delta(i)]; %delta
    %     y=[y beta(i)]; %beta
    %     z=[z integral(i)]; %integral
    % end
    % x=x';
    % y=y';
    % z=z';
    %
    % [f, gof] = fit([x, y],z,fit_name(1),'Normalize', 'on', 'Robust','Bisquare');
    % f_best=f;
    % gof_best=gof;
    % adjrsq_best = gof_best.adjrsquare;
    % sse_best = gof_best.sse;
    % for i=2:length(fit_name)
    %     [f, gof] = fit([x, y],z,fit_name(i), 'Normalize', 'on', 'Robust','Bisquare');
    %     if(gof.adjrsquare > adjrsq_best && gof.sse < sse_best)
    %         f_best=f;
    %         gof_best=gof;
    %         adjrsq_best = gof_best.adjrsquare;
    %         sse_best = gof_best.sse;
    %     end
    % end
    %
    % if(plots>0)
    %     figure()
    %     plot(f_best,[x y],z, 'Exclude', z < 100)
    %      zlim([0 max(z)]);
    %     xlabel('x - Delta')
    %     xticks(0:min(x):max(x))
    %     ylabel('y - Beta')
    %     zlabel("Integral")
    %     title("Integral VS Delta, beta")
    %     legend(type(f_best))
    %     grid on
    % end
    % %% 8 - Pi with different beta, priority
    % disp('8 - Pi with different beta, priority ')
    %
    % x = [];
    % y = [];
    % z = [];
    % parfor i=1:length(cell_in)
    %     x=[x beta(i)]; %beta
    %     y=[y priority(i)]; %priority
    %     z=[z pi_greco(i)]; %pi
    % end
    % x=x';
    % y=y';
    % z=z';
    % [f, gof] = fit([x, y],z,fit_name(1),'Normalize', 'on', 'Robust','Bisquare','Exclude', z < 0);
    % f_best=f;
    % gof_best=gof;
    % adjrsq_best = gof_best.adjrsquare;
    % sse_best = gof_best.sse;
    % for i=2:length(fit_name)
    %     [f, gof] = fit([x, y],z,fit_name(i), 'Normalize', 'on', 'Robust','Bisquare','Exclude', z < 0);
    %     if(gof.adjrsquare > adjrsq_best && gof.sse < sse_best)
    %         f_best=f;
    %         gof_best=gof;
    %         adjrsq_best = gof_best.adjrsquare;
    %         sse_best = gof_best.sse;
    %     end
    % end
    %
    % if(plots>0)
    %     figure()
    %     plot(f_best,[x y],z, 'Exclude', z < 0)
    %     zlim([-1 1]);
    %     xlabel('x - beta')
    %     ylabel('y - priority')
    %     zlabel("PI")
    %     title("PI VS Beta, priority")
    %     legend(type(f_best))
    %     grid on
    % end
    % %% 9 - Integral with different beta and priority
    % disp('9 - Integral with different beta and priority ')
    %
    % x = [];
    % y = [];
    % z = [];
    % parfor i=1:length(cell_in)
    %     x=[x beta(i)]; %beta
    %     y=[y priority(i)]; %priority
    %     z=[z integral(i)]; %integral
    % end
    % x=x';
    % y=y';
    % z=z';
    %
    % [f, gof] = fit([x, y],z,fit_name(1),'Normalize', 'on', 'Robust','Bisquare');
    % f_best=f;
    % gof_best=gof;
    % adjrsq_best = gof_best.adjrsquare;
    % sse_best = gof_best.sse;
    % for i=2:length(fit_name)
    %     [f, gof] = fit([x, y],z,fit_name(i), 'Normalize', 'on', 'Robust','Bisquare');
    %     if(gof.adjrsquare > adjrsq_best && gof.sse < sse_best)
    %         f_best=f;
    %         gof_best=gof;
    %         adjrsq_best = gof_best.adjrsquare;
    %         sse_best = gof_best.sse;
    %     end
    % end
    %
    % if(plots>0)
    %     figure()
    %     plot(f_best,[x y],z, 'Exclude', z < 100)
    %      zlim([0 max(z)]);
    %     xlabel('x - beta')
    %     ylabel('y - priority')
    %     zlabel("Integral")
    %     title("Integral VS Beta, priority")
    %     legend(type(f_best))
    %     grid on
    % end
    %%
    disp('==============================')
end