function [Analysis_Output] = APDAnalysis(file_folder, file_name)
%APDAnalysis Analysis of the optical mapping APD signals
%   APDAnalysis(x) performs the analysis of the file named x and gives
%   several outputs:
%   global maximum - length_max, t_length_max, 
%   baseline maximum - length_max1, t_length_max1, 
%   global minimum - length_min, t_length_min, 
%   fractional shortening (FS) - FS,
%   maximal shortening velocity - dlength_dt_min, length_dmin, t_length_dmin, 
%   maximal relaxation velocity - dlength_dt_max, length_dmax, t_length_dmax
%
% Analysis_Output = [length_max, t_length_max, length_max1, t_length_max1, length_min, t_length_min, FS, dlength_dt_min, length_dmin, t_length_dmin, dlength_dt_max, length_dmax, t_length_dmax]
%
% [length_max, t_length_max, length_max1, t_length_max1, length_min, t_length_min, FS, dlength_dt_min, length_dmin, t_length_dmin, dlength_dt_max, length_dmax, t_length_dmax]
% [length_max,t_length_max,length_max1,t_length_max1,length_min,t_length_min,FS,dlength_dt_min,length_dmin,t_length_dmin,dlength_dt_max,length_dmax,t_length_dmax]
%
%Written by Klemen Ziberna
%
%Update: includes fractional shortening output

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load data

file_name2 = file_name(1:end-4);

cell_data = importdata([file_folder,file_name]);

t = transpose(cell_data.data(:,1));
length = transpose(cell_data.data(:,2));


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical Analysis

% Numerical first derivative
dlength_dt = diff(length)./diff(t);

% Associated time values (central difference)
tt = t(1:end-1)+diff(t)./2;

% Maximum and minimum values

t_max = max(t);

[length_max, length_max_i] = max(length);
[length_min, length_min_i] = min(length);
[length_max1, length_max1_i] = max(length(1:length_min_i)); %max length before the min (i.e. baseline)

t_length_max = t(length_max_i);
t_length_min = t(length_min_i);
t_length_max1 = t(length_max1_i);

[dlength_dt_max, dlength_dt_max_i] = max(dlength_dt); %max relaxation velocity
[dlength_dt_min, dlength_dt_min_i] = min(dlength_dt); %max shortening velocity

%Corresponding length and t values at max/min dlength_dt (max relaxation and shortening velocity)

length_dmax = (length(dlength_dt_max_i)+length(dlength_dt_max_i+1))/2; %length at max dlength_dt (max relaxation velocity)
length_dmin = (length(dlength_dt_min_i)+length(dlength_dt_min_i+1))/2; %length at min dlength_dt (max shortening velocity)

t_length_dmax = tt(dlength_dt_max_i); %time at max dlength_dt (max relaxation velocity)
t_length_dmin = tt(dlength_dt_min_i); %time at min dlength_dt (max shortening velocity)

% Fractional shortening (FS = (length_max1 - length_min)./length_max1)

FS = ((length_max1 - length_min)./length_max1);


% Create output of all the calculated values

Analysis_Output = [length_max, t_length_max, length_max1, t_length_max1, length_min, t_length_min, FS, dlength_dt_min, length_dmin, t_length_dmin, dlength_dt_max, length_dmax, t_length_dmax];


%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot graphs

%Figure1: Data with tangents for max, min, minslope, maxslope
%Figure2: First derivative


% Calculate figure's y_min and y_max to set y-axis dimensions

%t_max = max(t);

%length_max = max(length);
%length_min = min(length);
%dlength_dt_max = max(dlength_dt);
%dlength_dt_min = min(dlength_dt);

% Draw plots


cell_fig = figure(1);

subplot(2,1,1);
set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(t,length, 'bo-'), grid; %plot length(t)
line([t(length_max1_i)-0.01, t(length_max1_i)+0.01], [length_max1, length_max1], 'LineWidth', 2, 'Color', [0 0.6 0]); %max baseline length line
line([t(length_max_i)-0.01, t(length_max_i)+0.01], [length_max, length_max], 'LineWidth', 2, 'Color', [0 0.6 0]); %max length line
line([t(length_min_i)-0.01, t(length_min_i)+0.01], [length_min, length_min], 'LineWidth', 2, 'Color', [0 0.6 0]); %min length line
line([t_length_dmin-0.005, t_length_dmin+0.005], [length_dmin-dlength_dt_min*(0.005), length_dmin-dlength_dt_min*(-0.005)], 'LineWidth', 2, 'Color', [0 0.6 0]); %max shortening velocity (dlength_dt_min)
line([t_length_dmax-0.005, t_length_dmax+0.005], [length_dmax-dlength_dt_max*(0.005), length_dmax-dlength_dt_max*(-0.005)], 'LineWidth', 2, 'Color', [0 0.6 0]); %max relaxation velocity (dlength_dt_max)

title([file_name2, ' - Cell Shortening'],'Interpreter', 'none');
xlim([0,t_max]);
ylim([length_min-0.5, length_max+0.5]);
xlabel 'Time (s)', ylabel 'Length (AU)';

subplot(2,1,2);
set(gcf,'Visible','off'); %prevent figures to pop up on screen
plot(tt,dlength_dt, 'rd-'), grid;
line([tt(dlength_dt_max_i)-0.01, tt(dlength_dt_max_i)+0.01], [dlength_dt_max, dlength_dt_max], 'LineWidth', 2, 'Color', [0 0.6 0]); %max length line
line([tt(dlength_dt_min_i)-0.01, tt(dlength_dt_min_i)+0.01], [dlength_dt_min, dlength_dt_min], 'LineWidth', 2, 'Color', [0 0.6 0]); %min length line
%title([file_name2, ' - Derivative']);
xlim([0,t_max]);
ylim([dlength_dt_min-100, dlength_dt_max+100]);
xlabel 'Time (s)', ylabel 'Velocity (AU/s)';


% Save the figure in the Matlab format

saveas(cell_fig, [file_name2], 'fig'); 

% Export plot to jpg

print('-djpeg', '-r300', [file_name2,'.jpg']);





end

