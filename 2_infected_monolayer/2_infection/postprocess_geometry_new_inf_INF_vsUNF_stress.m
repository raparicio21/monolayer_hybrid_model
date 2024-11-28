%% MONOLAYER POSTPROCESS: GEOMETRY

clear all;
clc;
close all;

% Initialize structures

% Load first step to retrieve some parameters
load ./results/raw/data_00001.mat

case_name        = 'infect';
make_video       = true;
make_figure      = true;
total_time_steps = 200;
delta_time       = 10; % minutes
sim_time_min     = 0:delta_time:total_time_steps*delta_time-delta_time; % minutes
sim_time_hou     = sim_time_min/60; % hours
number_of_cells  = numel(id_out);


% Load parameters
for i = 1:total_time_steps

    file_name = sprintf('./results/raw/data_%05i.mat', i); 
    load(file_name)

end


clear area_out pol_sizes_out perimeter_out location_out disp_vec_out   ...
      velo_vec_out tot_dist_out tot_disp_out t_dis_ve_out id_out       ...
      init_pos_out radius_out curr_pos_out tot_force_out int_force_out ...
      con_force_out rand_walk_out pro_force_out inf_crawl_out          ...
      orientation_min_s_I_vect orientation_min_s_II_vect               ...
      surround_out surround_2nd_out surround_3rd_out area_change       ...
      cellshape_out infected_out weighted_force_I weighted_force_II;

%% Data from mesh

% Column 1 : Integration_point_coord_x 
% Column 2 : Integration_point_coord_y
% Column 3 : S_x 
% Column 4 : S_y 
% Column 5 : S_xy 
% Column 6 : S_maxprin or sigma_I  (range [0  s_I_max]) positive values
% Column 7 : S_minprin or sigma_II (range [s_II_max 0]) negative values
% Column 8 : S_maxprin (all range)
% Column 9 : S_minprin (all range)


mean_sI         = zeros(1,total_time_steps);
mean_sII        = zeros(1,total_time_steps);
mean_LANS       = zeros(1,total_time_steps);
mean_MLSS       = zeros(1,total_time_steps);
weighted_sI     = zeros(1,total_time_steps);
weighted_sII    = zeros(1,total_time_steps);
weighted_LANS   = zeros(1,total_time_steps);
weighted_MLSS   = zeros(1,total_time_steps);
% Local average normal stress (LANS) = (SigmaI + SigmaII) /2
% Maximum local shear stress (MLSS) = (SigmaI - SigmaII) /2
stress_time     =  cell(1,total_time_steps);

infection_state =  cell(1,total_time_steps);


% Load parameters
for i = 1:total_time_steps

    file_name = sprintf('./results/raw_fem/fem_model_%05i.mat', i); 
    load(file_name)

    sI_i             = stress_solution(:,8)      ;
    sII_i            = stress_solution(:,9)      ;
    LANS_i           = (sI_i + sII_i)/2          ;
    MLSS_i           = (sI_i - sII_i)/2          ;


    %infection_triang  = triang_inf_id             ;

    stress_time{i}   = [sI_i sII_i LANS_i MLSS_i];
    % Column 1 : sI
    % Column 2 : sII
    % Column 3 : LANS 
    % Column 4 : MLSS 

    mean_sI(i)       = mean(sI_i)                ;
    mean_sII(i)      = mean(sII_i)               ;
    mean_LANS(i)     = mean(LANS_i)              ;
    mean_MLSS(i)     = mean(MLSS_i)              ;

    weighted_sI(i)   = sum(sI_i   .*triang_areas) / sum(triang_areas);
    weighted_sII(i)  = sum(sII_i  .*triang_areas) / sum(triang_areas);
    weighted_LANS(i) = sum(LANS_i .*triang_areas) / sum(triang_areas);
    weighted_MLSS(i) = sum(MLSS_i .*triang_areas) / sum(triang_areas);


    infection_state{i} = triang_inf_id;

end

clear area_triangles displacements_solution infection_triang LANS_i MLSS_i ...
    sI_i sII_i triang_inf_id triang_areas stress_solution


%% Postprocess data

meanforceAll = zeros(total_time_steps,1);
meanforceInf = zeros(total_time_steps,1);
meanforceSur = zeros(total_time_steps,1);

for i = 91:total_time_steps

    inf_state_iter = infection_state{i};
    MLSS_iter      = stress_time{i}(:,4);
    
    % All
    N_elements_all  = numel(inf_state_iter);
    MLSS_all        = MLSS_iter;
    meanforceAll(i) = sum(MLSS_all)/N_elements_all;

    % Infected
    N_elements_inf  = sum(inf_state_iter);
    MLSS_inf        = inf_state_iter.*MLSS_iter;
    meanforceInf(i) = sum(MLSS_inf)/N_elements_inf;

    % Uninfected
    Sur_Elem        = imcomplement(inf_state_iter);
    N_elements_Sur  = sum(Sur_Elem);
    MLSS_Sur        = Sur_Elem.*MLSS_iter;
    meanforceSur(i) = sum(MLSS_Sur)/N_elements_Sur;




end

%% Plots

%............................  MLSS  ...................................
% Plot MLSS mean over time
% between this two frame points
time_1 = 91;
time_2 = 200;

%plot(sim_time_min(time_1:time_2) , meanforceAll(time_1:time_2),'LineWidth', 3)
hold on
plot(sim_time_min(time_1:time_2) , meanforceInf(time_1:time_2),'LineWidth', 3)
plot(sim_time_min(time_1:time_2) , meanforceSur(time_1:time_2),'LineWidth', 3)
title ('mean \sigma_I [Pa]'    );
ylabel('mean \sigma_I [Pa]'    );
xlabel('Time [hours]'          );
% ylim([20,50])
% xlim([0,4000])
%xlim([800,1550])
picname    = "mean_sI_time.svg";

ax = gca;
ax.LineWidth = 2; % Set line width to 2 points
ax.XAxis.FontSize = 20; % Set font size of x axis label to 12 points
ax.YAxis.FontSize = 20; % Set font size of y axis label to 12 points
ax.XAxis.FontName = 'Arial';
ax.YAxis.FontName = 'Arial';
ax.Box            = 'off'  ;
hold off


saveas(gcf, picname);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT DYNAMICS FEM %%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Plot mean weighted values overtime %%%%%%%%%%%%%%%%%%%%%%%%
% .................. s_I ...................
plot(sim_time_min , weighted_sI)
title ('mean s_I [Pa]'         );
ylabel('s_I [Pa]'              );
xlabel('Simulation time [min]' );
% ylim([20,50])
% xlim([0,4000])
hold on
plot(sim_time_min , mean_sI)
legend('mean s_I weighted','mean s_I');
hold off

if make_figure == true
    picname    = "sI_" + case_name + ".png";
    saveas(gcf, picname);
end



%% Plot paper
%............................  sigma_I  ...................................
% Plot s_I over time
% between this two frame points
time_1 = 88;
time_2 = 164;

plot(sim_time_min(time_1:time_2) , mean_sI(time_1:time_2),'LineWidth', 3)
title ('mean \sigma_I [Pa]'    );
ylabel('mean \sigma_I [Pa]'    );
xlabel('Time [hours]'          );
% ylim([20,50])
% xlim([0,4000])
%xlim([800,1550])
picname    = "mean_sI_time.svg";

ax = gca;
ax.LineWidth = 2; % Set line width to 2 points
ax.XAxis.FontSize = 20; % Set font size of x axis label to 12 points
ax.YAxis.FontSize = 20; % Set font size of y axis label to 12 points
ax.XAxis.FontName = 'Arial';
ax.YAxis.FontName = 'Arial';
ax.Box            = 'off'  ;


saveas(gcf, picname);




%%
% .................. s_II ...................
plot(sim_time_min , weighted_sII)
legend('mean s_{II} weighted'   );
title ('mean s_{II} [Pa]'       );
ylabel('s_{II} [Pa]'            );
xlabel('Simulation time [min]'  );
hold on
plot(sim_time_min , mean_sII)
legend('mean s_II weighted','mean s_II');
hold off

if make_figure == true
    picname    = "sII_" + case_name + ".png";
    saveas(gcf, picname);
end

% .................. LANS ...................
plot(sim_time_min , weighted_LANS)
legend('mean LANS weighted'     );
title ('mean LANS [Pa]'         );
ylabel('s_{LANS} [Pa]'          );
xlabel('Simulation time [min]'  );
hold on
plot(sim_time_min , mean_LANS)
legend('mean LANS weighted','mean LANS');
hold off

if make_figure == true
    picname    = "sLANS_" + case_name + ".png";
    saveas(gcf, picname);
end

% .................. MLSS ...................
plot(sim_time_min , weighted_MLSS)
legend('mean MLSS weighted'     );
title ('mean MLSS [Pa]'         );
ylabel('s_{MLSS} [Pa]'          );
xlabel('Simulation time [min]'  );
hold on
plot(sim_time_min , mean_MLSS)
legend('mean MLSS weighted','mean MLSS');
hold off

if make_figure == true
    picname    = "sMLSS_" + case_name + ".png";
    saveas(gcf, picname);
end


%% Plot paper
%............................  MLSS  ...................................
% Plot MLSS over time
% between this two frame points
time_1 = 88;
time_2 = 164;

plot(sim_time_min(time_1:time_2) , mean_MLSS(time_1:time_2),'LineWidth', 3)
title ('mean MLSS [Pa]'    );
ylabel('mean MLSS [Pa]'    );
xlabel('Time [hours]'          );
% ylim([20,50])
% xlim([0,4000])
%xlim([800,1550])
picname    = "mean_MLSS_time.svg";

ax = gca;
ax.LineWidth = 2; % Set line width to 2 points
ax.XAxis.FontSize = 20; % Set font size of x axis label to 12 points
ax.YAxis.FontSize = 20; % Set font size of y axis label to 12 points
ax.XAxis.FontName = 'Arial';
ax.YAxis.FontName = 'Arial';
ax.Box            = 'off'  ;


saveas(gcf, picname);

%% MLSS paper
% Compare distribution of two frames (for example initial vs final)
edges             = 0:2:100; % histogram edges (range)
frames_to_compare = [88 164];

% Get maps
sI_map   = cell(1,length(frames_to_compare));
sII_map  = cell(1,length(frames_to_compare));
LANS_map = cell(1,length(frames_to_compare));
MLSS_map = cell(1,length(frames_to_compare));

% Load stress maps
for i = 1:length(frames_to_compare)

    file_name = sprintf('./results/raw_fem/fem_model_%05i.mat', frames_to_compare(i)); 
    load(file_name)

    sI_i             = stress_solution(:,8)      ;
    sII_i            = stress_solution(:,9)      ;
    LANS_i           = (sI_i + sII_i)/2          ;
    MLSS_i           = (sI_i - sII_i)/2          ;

    sI_map{i}   = sI_i  ;
    sII_map{i}  = sII_i ;
    LANS_map{i} = LANS_i;
    MLSS_map{i} = MLSS_i;



end


for i = 1:length(frames_to_compare)

    hold on

    hist_data          = sI_map{i}(:);
    histogram(hist_data, 'BinEdges', edges);

    hold off

end


xlim  ([0, 100]            ) ;     ylim  ([0, 400]                 );
xlabel('MLSS [Pa]'   ) ;     ylabel('frequency of appearence' );
name1= 'primero';
name2= 'segundo';
%name3= ['t: ' num2str(time_min(frames_to_compare(3))) ' min'];
legend(name1,name2)
title('MLSS distribution')

ax = gca;
ax.LineWidth = 2; % Set line width to 2 points
ax.XAxis.FontSize = 16; % Set font size of x axis label to 12 points
ax.YAxis.FontSize = 16; % Set font size of y axis label to 12 points
ax.XAxis.FontName = 'Arial';
ax.YAxis.FontName = 'Arial';

picname    = "compare_dist_s_MLSS.svg";
saveas(gcf, picname,'svg');


%%
%%%%%%%%%%%%%%% Plot distribution of stresses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sI distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "sI_" + case_name + ".mp4";
    vidObj           = VideoWriter(video_name,'MPEG-4');
    vidObj.FrameRate = 10;    % Set the frame rate to 10 frames per second
    open(vidObj);
end
edges = -70:5:120; % histogram edges (range)

for i = 1:total_time_steps

    hist_data          = stress_time{i}(:,1);

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([-70, 120]         ) ;     ylim  ([0, 1000]                  );
    xlabel('s_{I} [Pa]'       ) ;     ylabel('frequency of appearence' );

    txt = ['t: ' num2str(sim_time_min(i)) ' min'];
    title({['s_I distribution - ' txt], ''})

    % Add frame to video object
    frame       = getframe(gcf);
    if make_video == true
        writeVideo(vidObj,frame);
    end
    % Clear figure for next frame
    clf;
end

% Close video object
if make_video == true
    close(vidObj);
end









% sII distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "sII_" + case_name + ".mp4";
    vidObj           = VideoWriter(video_name,'MPEG-4');
    vidObj.FrameRate = 10;    % Set the frame rate to 10 frames per second
    open(vidObj);
end
edges = -140:5:70; % histogram edges (range)

for i = 1:total_time_steps

    hist_data          = stress_time{i}(:,2);

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([-140, 70]          ) ;     ylim  ([0, 1000]                  );
    xlabel('s_{II} [Pa]'       ) ;     ylabel('frequency of appearence' );

    txt = ['t: ' num2str(sim_time_min(i)) ' min'];
    title({['s_{II} distribution - ' txt], ''})

    % Add frame to video object
    frame       = getframe(gcf);
    if make_video == true
        writeVideo(vidObj,frame);
    end
    % Clear figure for next frame
    clf;
end

% Close video object
if make_video == true
    close(vidObj);
end







% sLANS distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "sLANS_" + case_name + ".mp4"   ;
    vidObj           = VideoWriter(video_name,'MPEG-4');
    vidObj.FrameRate = 10;    % Set the frame rate to 10 frames per second
    open(vidObj);
end
edges = -140:5:140; % histogram edges (range)

for i = 1:total_time_steps

    hist_data          = stress_time{i}(:,3);

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([-140, 140]        ) ;     ylim  ([0, 1000]                  );
    xlabel('s_{LANS} [Pa]'    ) ;     ylabel('frequency of appearence' );

    txt = ['t: ' num2str(sim_time_min(i)) ' min'];
    title({['s_{LANS} distribution - ' txt], ''})

    % Add frame to video object
    frame       = getframe(gcf);
    if make_video == true
        writeVideo(vidObj,frame);
    end
    % Clear figure for next frame
    clf;
end

% Close video object
if make_video == true
    close(vidObj);
end







%%

% sMLSS distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "sMLSS_" + case_name + ".mp4"   ;
    vidObj           = VideoWriter(video_name,'MPEG-4');
    vidObj.FrameRate = 10;    % Set the frame rate to 10 frames per second
    open(vidObj);
end
edges = -10:5:140; % histogram edges (range)

for i = 1:total_time_steps

    hist_data          = stress_time{i}(:,4);

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([-10, 140]          ) ;     ylim  ([0, 2500]                );
    xlabel('s_{MLSS} [Pa]'    ) ;     ylabel('frequency of appearence' );

    txt = ['t: ' num2str(sim_time_min(i)) ' min'];
    title({['s_{MLSS} distribution - ' txt], ''})

    % Add frame to video object
    frame       = getframe(gcf);
    if make_video == true
        writeVideo(vidObj,frame);
    end
    % Clear figure for next frame
    clf;
end

% Close video object
if make_video == true
    close(vidObj);
end





%%

% Save postrpocess parameters ............................................
filename = "var_"+ case_name;
save(filename, 'area_time', 'polygon_size_time', 'perimeter_time', 'cell_shape_time'           ...
             , 'location_time', 'infection_time', 'prctile_area_step','median_pol_size_step'   ...
             , 'prctile_perim_step', 'prctile_cell_shape_step', 'total_time_steps'             ...
             , 'number_of_cells', 'displacement_time_magn','prctile_disp_vec_magn_step'        ...
             , 'prctile_velocity_magn_step', 'prctile_tot_dist_step', 'prctile_tot_disp_step'  ...
             , 'prctile_dir_ratio_step');
