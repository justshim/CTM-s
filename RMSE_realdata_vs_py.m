clearvars
close all
clc

%% matlab data
path=strcat('H:\Il mio Drive\Tesi magistrale\CTMs-identification\fnc\extracted_data\data.mat');
aa = load(path, '*');
data = aa.sensor_sum;
clear aa;

%% python data
D = readtable("C:\A_Tesi\Python\CTM-s\model\density.csv");
dens_py = table2array(D);
F = readtable("C:\A_Tesi\Python\CTM-s\model\flow.csv");
flow_py = table2array(F);

%% plots and RMSE
RMSE = [];
last_fig_num = get(gcf,'Number');
figure(last_fig_num)

for CELL = 1:13
    y = data(CELL).flow';
    yhat = flow_py(:,CELL);
    yhat = [yhat; 0];
    RMSE = [RMSE; sqrt(mean((y - yhat).^2))];  % Root Mean Squared Error
    n_row = 3;

    x = linspace(1, 24, 8640);

    subplot(n_row,ceil((13)/n_row),CELL)
    plot(x, y)
    hold on
    plot(x, yhat)
%     ax = gca;
%     ax.YLim = [0,2500];
    title_str1 = ['Cell - ',CELL];
    title(title_str1)
    grid on
%     f1 = figure;
%     scatter(sensor_sum(10).density,sensor_sum(10).vehicle_speed,[],'filled')
%     hold on
%     grid on
%     f1.WindowState = 'maximized';
%     ax = gca();
%     font_sz = 25;
%     ax.XAxis.FontSize = font_sz; ax.XAxis.TickLabelInterpreter = 'latex';
%     ax.YAxis.FontSize = font_sz; ax.YAxis.TickLabelInterpreter = 'latex';
%     ax.XAxis.Label.String = '$\rho[veh/km]$'; ax.XAxis.Label.FontSize = font_sz;
%     ax.XAxis.Label.Interpreter = 'latex';
%     ax.YAxis.Label.String = '$speed[km/h]$'; ax.YAxis.Label.FontSize = font_sz;
%     ax.YAxis.Label.Interpreter = 'latex';
%     exportgraphics(f1,['density_speed.pdf'],...
%                    'BackgroundColor','none');
%     exportgraphics(f1,['density_speed.eps'],...
%                    'BackgroundColor','none');


end

