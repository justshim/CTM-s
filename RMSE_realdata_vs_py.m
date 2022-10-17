clearvars
close all
clc

%% matlab data
path=strcat('H:\Il mio Drive\Tesi magistrale\CTMs-identification\fnc\extracted_data\data.mat');
aa = load(path, '*');
data = aa.sensor_sum;
clear aa;

%% python data
%D = readtable("H:\Il mio Drive\Tesi magistrale\Python\CTM-s\model\density_station.csv");
D = readtable("H:\Il mio Drive\Tesi magistrale\Python\CTM-s\model\density_no_station.csv");
dens_py = table2array(D);
%F = readtable("H:\Il mio Drive\Tesi magistrale\Python\CTM-s\model\flow_station.csv");
F = readtable("H:\Il mio Drive\Tesi magistrale\Python\CTM-s\model\flow_no_station.csv");
flow_py = table2array(F);

%% plots and RMSE
RMSE = [];
last_fig_num = get(gcf,'Number');
figure(last_fig_num)

for CELL = 1:13
    y = data(CELL).density'; % realdata
    yhat = dens_py(:,CELL); % python
    
    yhat = [yhat; 0];
    
    if (sum(isnan(y))~=0)
        y(isnan(y))=0;
    end
    RMSE = [RMSE; sqrt(mean((y - yhat).^2))];  % Root Mean Squared Error
    n_row = 3;

    x = linspace(1, 24, 8640);
y_mm=movmean(y, 90*8640/1440);


 f1 = figure;
    h = plot(x,y);
    h.LineWidth = 3;
    h.Color=[0.4, 0.9, 1];
    hold on
    grid on
    h = plot(x,y_mm);
    h.LineWidth = 3;
    h.Color="blue";
    h = plot(x,yhat);
    h.LineWidth = 3;
    h.Color="red";
    f1.WindowState = 'maximized';
    ax = gca();
    %ax.YLim = [0; 2500];
    font_sz = 25;
    ax.XAxis.FontSize = font_sz; ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.FontSize = font_sz; ax.YAxis.TickLabelInterpreter = 'latex';
    ax.XAxis.Label.String = '$Time \mathrm{[h]}$'; ax.XAxis.Label.FontSize = font_sz;
    ax.XAxis.Label.Interpreter = 'latex';
    ax.YAxis.Label.String = '$\rho \mathrm{[veh/km]}$'; ax.YAxis.Label.FontSize = font_sz;
    ax.YAxis.Label.Interpreter = 'latex';
    str_pdf=strcat('density_',num2str(CELL), '.pdf');
    str_eps=strcat('density_',num2str(CELL), '.eps');
    exportgraphics(f1,str_pdf,'BackgroundColor','none');
    exportgraphics(f1,str_eps,'BackgroundColor','none');


end

