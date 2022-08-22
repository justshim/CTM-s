clearvars
close all
clc

delete(gcp('nocreate'));
parpool('threads');

%% Load files
load_path=strcat(pwd,'/data models/4 vars\model_integral_coeff.mat');
min = load(load_path, '*');
model_integral_coeff = min.Coefficients_ide_integral;

load_path=strcat(pwd,'/data models/4 vars\model_integral_terms.mat');
min = load(load_path, '*');
model_integral_terms = min.ModelTerms_ide_integral;

load_path=strcat(pwd,'/data models/4 vars\model_pi_coeff.mat');
min = load(load_path, '*');
model_pi_coeff = min.Coefficients_ide_pi;

load_path=strcat(pwd,'/data models/4 vars\model_pi_terms.mat');
min = load(load_path, '*');
model_pi_terms = min.ModelTerms_ide_pi;
clear min

%% Build model functions
model_integral = struct();
model_integral.Coefficients = model_integral_coeff;
model_integral.ModelTerms = model_integral_terms;

fun_integral = polyn2sym_mod(model_integral);

model_pi = struct();
model_pi.Coefficients = model_pi_coeff;
model_pi.ModelTerms = model_pi_terms;

fun_pi = polyn2sym_mod(model_pi);
f_int = string(char(fun_pi));
f_int = replace(f_int,["x(3)","x(4)"],["240","0.07"]);
f_int = replace(f_int,["x(1)","x(2)"],["x","y"]);
f_int = replace(f_int,"@(x)","@(x,y)");
f_int = str2func(f_int);
figure()
fmesh(f_int, [0 13])



%fun_pi = @(x) -fun_pi(x);
f = string(char(fun_pi));
f = replace(f,["x(3)","x(4)"],["240","0.07"]);
f = replace(f,["x(1)","x(2)"],["x","y"]);
f = replace(f,"@(x)","@(x,y)");
f = str2func(f);

f = @(x,y) -f(x,y);
figure()
fmesh(f, [1 13])


% w = 2460;
%w = (8991+14)/2;
%fun_both = @(x) (fun_integral(x) + w*fun_pi(x));

%% Find minima
% x(1) = i; x(2) = j; x(3) = delta; x(4) = beta

% linear constraints
A = [-1, 1, 0, 0;
     1, -1, 0, 0;];  %A*x <= b
b = [3; -1]; %station lunga almeno 1 e massimo 3: -x1 +x2<=3 && -(-x1 + x2)<=-1
Aeq = [];       %Aeq*x = beq
beq = [];
%nonlcon = [];
% nonlinear constraints
lb = [1 2 60 0.01];
ub = [12 13 720 0.2];
x0 = (lb + ub)/2;
% lb = [];
% ub = [];
%x0=[0 0 0 0];
%options = optimoptions('fmincon','Algorithm','interior-point');


fprintf("\n");
fprintf("=== Optimizing integral Delta ===");
[sol_integral,y_sol_integral,exitflag,output]  = fmincon(fun_integral,x0,A,b,Aeq,beq,lb,ub);

fprintf("\n");
fprintf("=== Optimizing Pi ===");
[sol_pi,y_sol_pi,exitflag,output] = fmincon(fun_pi,x0,A,b,Aeq,beq,lb,ub);

% fprintf("\n");
% fprintf("=== Optimizing both ===");
% [sol_both,y_sol_both,exitflag,output] = fmincon(fun_both,x0,A,b,Aeq,beq,lb,ub);

