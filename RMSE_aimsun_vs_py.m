clearvars
close all
clc

% sec_phi_1_24h_realsmooth_incr
% sec_phi_1_24h_realsmooth
% sec_phi_1_24h_real
% sec_phi_1_24h_doublepeak
% sec_phi_1_24h_singlepeak

% _no_station

file_name = "sec_phi_1_24h_realsmooth";
n_lanes = 1;
if(n_lanes == 1)
    path=strcat('C:\A_Tesi\Aimsun\CTM-s\A2 (1 lane)\Model\Resources\Scripts\', file_name,'.csv');
end

%% extract data
import_raw = importdata(path, ';');
cell_raw = import_raw.textdata;
index=1;
for i = 1:size(cell_raw,2)
    if(isempty(cell_raw{2,i}))
        cell_raw(2:end,i) = num2cell(import_raw.data(:,index));
        index=index+1;
    end
end

% create an empty structure
data = struct();

%% Reading data: Header
%The first row is the header, thus it is used to create the structure
%fields
header = cell_raw(1,:);

% create all the fields
for j = 1 : length(header)
    field_name = char(header(1,j));
    data.(field_name) = [];
end

%% Reading data: Data rows
% fill all the fields row by row
data_field_name = fieldnames(data);
for i = 1:numel(data_field_name)
    % get the data field
    field_i_name = char(data_field_name(i));
    % select the data in the row
    data.(field_i_name) = [cell_raw(2:end, i)];
end

%% Find sensor names
id_section = unique(string(data.name));

%% Extract useful data
% extract from the whole data only the values that interest us
% find the indeces associated with the different sensors

section(length(id_section)) = struct(); %preallocate space for speed-up
for j = 1:length(id_section)
    section_index = strcmp(string(data.name),id_section(j));
    section(j).travel_time_avg = cell2mat(data.travel_time_avg(section_index));
    section(j).id = id_section(j);
    section(j).flow = cell2mat(data.flow(section_index));
    section(j).count = cell2mat(data.count(section_index));
    section(j).delay_time_avg = cell2mat(data.delay_time_avg(section_index));
    section(j).veh_avg_speed = cell2mat(data.speed(section_index));
    section(j).density = cell2mat(data.density(section_index));
    section(j).ending_s_time = data.finish_time(section_index);
    section(j).starting_s_time = data.init_time(section_index);
    section(j).stop_time_avg = cell2mat(data.stop_time_avg(section_index));
    section(j).num_stops = cell2mat(data.num_stops(section_index));
    section(j).queue_avg = cell2mat(data.queue_avg(section_index));
    section(j).queue_max = cell2mat(data.queue_max(section_index));
    section(j).virtual_queue_avg = cell2mat(data.virtual_queue_avg(section_index));
end
%% cleaning of data
for j = 1:length(id_section)
    bin = isnan(section(j).num_stops);
    section(j).num_stops(bin) = 0;
    for k = 1:length(section(j).delay_time_avg)
        if (section(j).delay_time_avg(k) == -1)
            section(j).delay_time_avg(k) = 0;
        end
        if (section(j).stop_time_avg(k) == -1)
            section(j).stop_time_avg(k) = 0;
        end
    end
end
%% Divide data between main sections and service station sections
jj = 1;
kk = 1;
for j = 1:length(id_section)
    if (contains(section(j).id, "Cell"))
        section(j).id = str2double(erase(section(j).id, "Cell "));
        section_main_temp(jj)=section(j);
        jj = jj+1;
    elseif (contains(section(j).id, "Station"))
        section_station(kk)=section(j);
        kk = kk+1;
    end

end
jj = 1;
for j = 1:length(section_main_temp)
    for k = 1:length(section_main_temp)
        if(section_main_temp(k).id == j)
            section_main(jj)=section_main_temp(k);
            jj = jj+1;
            break;
        end
    end
end
clear section_main_temp
for j = 1:length(section_main)
    section_main(j).id = num2str(section_main(j).id);
end
%% python data
D = readtable("C:\A_Tesi\Python\CTM-s\model\density.csv");
dens = table2array(D);
F = readtable("C:\A_Tesi\Python\CTM-s\model\flow.csv");
flow = table2array(F);


%% plots and RMSE
RMSE = [];

last_fig_num = get(gcf,'Number');
figure(last_fig_num)
 f1 = figure;
for CELL = 1:13
    y = section_main(CELL).flow; %aimsun
    y_mm=movmean(y, 121);
    yhat_temp = flow(:,CELL); %py
    yhat_temp = [yhat_temp; 0];

    z=1;
    yhat=zeros(1440,1);
    for p=1:6:length(yhat_temp)
        yhat(z) = (yhat_temp(p) + yhat_temp(p+1) + yhat_temp(p+2) + yhat_temp(p+3) + yhat_temp(p+4) + yhat_temp(p+5))/6;
        z=z+1;
    end

    RMSE = [RMSE; sqrt(mean((y - yhat).^2))];  % Root Mean Squared Error


    n_row = 3;

    x = linspace(1, 24, 1440);

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
    font_sz = 25;
    ax.XAxis.FontSize = font_sz; ax.XAxis.TickLabelInterpreter = 'latex';
    ax.YAxis.FontSize = font_sz; ax.YAxis.TickLabelInterpreter = 'latex';
    ax.XAxis.Label.String = '$Time [h]$'; ax.XAxis.Label.FontSize = font_sz;
    ax.XAxis.Label.Interpreter = 'latex';
    ax.YAxis.Label.String = '$Flow [veh/h]$'; ax.YAxis.Label.FontSize = font_sz;
    ax.YAxis.Label.Interpreter = 'latex';
    str_pdf=strcat('flow_',num2str(CELL), '.pdf');
    str_eps=strcat('flow_',num2str(CELL), '.eps');
%     exportgraphics(f1,str_pdf,'BackgroundColor','none');
%     exportgraphics(f1,str_eps,'BackgroundColor','none');

end

