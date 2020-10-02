function plot_PracId(virus_ind)
% % 
% characterizing practical identifiability (PI) using pairwise correlation 
% between the corresponding parameter values estimated in the final
% iteration of iABC. See the SI file of the manuscript for more details.

warning off;
indexing; % info about parameter, variable and virus indexing

%% figure formatting
frame_x0 = 0.1; frame_y0 = 0.1; frame_x1 = 0.575; frame_y1 = 0.8;
fs = 32;

%% files and information

img_file_name = [virus_name{1, virus_ind}, '_PI']; % storing the PI image

% extracting the parameter combinations selected in the final iteration of
% iABC. PRM_nest: stores the selections for each iteration
load(['session_main_', virus_name{1, virus_ind}, '.mat'], 'PRM_nest');
PRM_last = PRM_nest(:, :, size(PRM_nest, 3));

%% calculation
rho = round(corr(PRM_last), 2); % estimating correlation

% "rho1" is the lower triangular part of "rho". note "rho" is symmetric
for x_ind = 1:size(rho, 1)-1
    for y_ind = x_ind:size(rho, 2)-1
        rho1(x_ind, y_ind) = NaN;
    end
    for y_ind = 1:x_ind
        rho1(x_ind, y_ind) = rho(x_ind+1, y_ind);
    end    
end

%% making the figure

% rows and cols name: parameter
xvalues = prm_name(1:size(rho, 1)-1);
yvalues = prm_name(2:size(rho, 1));

% initializing the figure
f = figure('Units', 'normalized', 'Position',[frame_x0 frame_y0 frame_x1 frame_y1]);

% displaying the corr matrix as a heatmap
h = heatmap(f, xvalues, yvalues, rho1, 'MissingDataLabel', {},...
    'GridVisible','off',...
    'MissingDataColor',[0.95, 0.95, 0.95]);
% heatmap formatting
h.ColorLimits = [-1 1]; % corr values range
h.Title = {['Pairwise correlation (', virus_name{1, virus_ind}, ')']}; % title
h.FontSize = fs; % font size

% user defined colormap for the heatmap
col_red = [linspace(0, 0, 8), linspace(0, 1, 24), linspace(1, 0.5, 32)]';
col_green = [linspace(0, 0, 8), linspace(0, 1, 24), linspace(1, 0, 24), linspace(0, 0, 8)]';
col_blue = [linspace(0.5, 1, 32), linspace(1, 0, 24), linspace(0, 0, 8)]';
h.Colormap = [col_red, col_green, col_blue];

print(f, img_file_name,'-djpeg','-r960') 
% close all;
end