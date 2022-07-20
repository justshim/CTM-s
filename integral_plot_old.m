close all
clearvars
clc

path = readtable('C:/A_Tesi/Python/CTM-s/opti_data.xls');
%% Different i, delta (5-60 min); j = i+2

x = table2array(path(:, "i"));
y = table2array(path(:, "delta"));
z = table2array(path(:, "integral"));
n=length(z);
x_grid=linspace(0, 11, 50);% fare dinamico
y_grid=linspace(0, 400, 50);

[xq, yq]=meshgrid(x_grid,y_grid); %crea una griglia cartesiana 2D o 3D
xq_vec=xq(:);
yq_vec=yq(:); %serve per poter visualizzare la nostra stima

Phi2=[ones(n,1), x, y, x.^2, y.^2, x.*y, x.^3, y.^3, (x.^2).*y, x.*(y.^2)...
    x.^4, y.^4, (x.^3).*y, x.*(y.^3), (x.^2).*(y.^2)];


[theta2, std_theta2]=lscov(Phi2, z); %calcola il vettore dei parametri e la sua std dev

q2=length(theta2);
y_hat2=Phi2*theta2; %stima
epsilon2=z-y_hat2; %residui
SSR2=epsilon2'*epsilon2

Phi2_grid=[ones(length(xq_vec), 1), xq_vec, yq_vec, xq_vec.^2, yq_vec.^2, xq_vec.*yq_vec, xq_vec.^3, yq_vec.^3,  (xq_vec.^2).*yq_vec, xq_vec.*(yq_vec.^2)...
    xq_vec.^4, yq_vec.^4, (xq_vec.^3).*yq_vec, xq_vec.*(yq_vec.^3), (xq_vec.^2).*(yq_vec.^2)];

curva2=Phi2_grid*theta2;
curva2_matrix=reshape(curva2, size(xq));

figure(2)
mesh(xq, yq, curva2_matrix)
hold on
plot3(x,y,z, '*')
xlabel('Cell IN')
ylabel('delta (timesteps)')
zlabel("Integral Delta")
title("Integral Delta VS i, delta")
grid on
legend('polinomio di grado 4', 'dati')

FPE2=((n+q2)/(n-q2))*SSR2
AIC2=2*q2/n +log(SSR2)
MDL2=log(n)*q2/n + log(SSR2)