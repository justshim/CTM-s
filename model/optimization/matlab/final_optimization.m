clearvars
close all
clc

delete(gcp('nocreate'));
parpool('threads');

%% Load files
load_path=strcat(pwd,'\model_integral_coeff.mat');
aaa = load(load_path, '*');
model_integral_coeff = aaa.Coefficients_ide_integral;

load_path=strcat(pwd,'\model_integral_terms.mat');
aaa = load(load_path, '*');
model_integral_terms = aaa.ModelTerms_ide_integral;

load_path=strcat(pwd,'\model_pi_coeff.mat');
aaa = load(load_path, '*');
model_pi_coeff = aaa.Coefficients_ide_pi;

load_path=strcat(pwd,'\model_pi_terms.mat');
aaa = load(load_path, '*');
model_pi_terms = aaa.ModelTerms_ide_pi;
clear aaa

%% Build model functions
model_integral = struct();
model_integral.Coefficients = model_integral_coeff;
model_integral.ModelTerms = model_integral_terms;

fun_integral = polyn2sym_mod(model_integral);

model_pi = struct();
model_pi.Coefficients = model_pi_coeff;
model_pi.ModelTerms = model_pi_terms;

fun_pi = polyn2sym_mod(model_pi);
fun_pi = @(x) -fun_pi(x);

% w = 2460;
w = (8991+14)/2;
fun_both = @(x) (fun_integral(x) + w*fun_pi(x));

%% Find minima
% x(1) = i; x(2) = j; x(3) = delta; x(4) = beta

% linear constraints
A = [1 1 0 0];  %A*x < b
b = 3;
Aeq = [];       %Aeq*x = beq
beq = [];

% nonlinear constraints
lb = [1 2 60 0.01];
ub = [12 13 720 0.2];
x0 = (lb + ub)/2;

fprintf("\n");
fprintf("=== Optimizing integral Delta ===");
sol_integral = fmincon(fun_integral,x0,A,b,Aeq,beq,lb,ub);
fprintf("\n");
fprintf("=== Optimizing Pi ===");
sol_pi = fmincon(fun_pi,x0,A,b,Aeq,beq,lb,ub);
fprintf("\n");
fprintf("=== Optimizing both ===");
sol_both = fmincon(fun_both,x0,A,b,Aeq,beq,lb,ub);

