%% MONOLAYER POSTPROCESS: GEOMETRY

clear all;
clc;
close all;

% Initialize structures

% Load first step to retrieve some parameters
load ./results/raw/data_00001.mat

case_name        = 'prueba';
make_video       = true;
make_figure      = true;
total_time_steps = 200;
delta_time       = 20; % minutes
sim_time_min     = 0:delta_time:total_time_steps*delta_time-delta_time; % minutes
sim_time_hou     = sim_time_min/60; % hours
number_of_cells  = numel(id_out);

% Cell geometry
area_time                    = zeros(number_of_cells,total_time_steps);
area_change_time             = zeros(number_of_cells,total_time_steps);
polygon_size_time            = zeros(number_of_cells,total_time_steps);
perimeter_time               = zeros(number_of_cells,total_time_steps);
cell_shape_time              = zeros(number_of_cells,total_time_steps);

% Cell state
location_time                = zeros(number_of_cells,total_time_steps);
infection_time               = zeros(number_of_cells,total_time_steps);

% Cell kinematics
displacement_time_x          = zeros(number_of_cells,total_time_steps);
displacement_time_y          = zeros(number_of_cells,total_time_steps);
displacement_time_magn       = zeros(number_of_cells,total_time_steps);

% Cell total kinematics
total_distance               = zeros(number_of_cells,total_time_steps);
total_displacement           = zeros(number_of_cells,total_time_steps);
directionality_ratio         = zeros(number_of_cells,total_time_steps);

% .........................................................................
% Cell forces (at the center of the cell/nuclei)
% .........................................................................
% Cell-cell interaction forces
interaction_time_x           = zeros(number_of_cells,total_time_steps);
interaction_time_y           = zeros(number_of_cells,total_time_steps);
% Cell contraction 
contraction_time_x           = zeros(number_of_cells,total_time_steps);
contraction_time_y           = zeros(number_of_cells,total_time_steps);
% Cell protrusion
protrusion_time_x            = zeros(number_of_cells,total_time_steps);
protrusion_time_y            = zeros(number_of_cells,total_time_steps);
% Cell crawling: infection
% inf_craw_time_x              = zeros(number_of_cells,total_time_steps);
% inf_craw_time_y              = zeros(number_of_cells,total_time_steps);
% Random walk: not included anymore
% Cell resultant force
tot_forc_time_x              = zeros(number_of_cells,total_time_steps);
tot_forc_time_y              = zeros(number_of_cells,total_time_steps);


% Load parameters
for i = 1:total_time_steps

    file_name = sprintf('./results/raw/data_%05i.mat', i); 
    load(file_name)

    % Cell geometry
    area_time                    (:,i) = area_out                            ;
    area_change_time             (:,i) = area_change                         ;
    polygon_size_time            (:,i) = pol_sizes_out                       ;
    perimeter_time               (:,i) = perimeter_out                       ;
    cell_shape_time              (:,i) = perimeter_out ./ (area_out).^(0.5)  ;

    % Cell state
    location_time                (:,i) = location_out                        ;
    infection_time               (:,i) = infected_out                        ;
    
    % Cell kinematics
    displacement_time_x          (:,i) = disp_vec_out(:,1)                   ;
    displacement_time_y          (:,i) = disp_vec_out(:,2)                   ;
    displacement_time_magn       (:,i) = (displacement_time_x(:,i).^2  + ...
                                          displacement_time_y(:,i).^2).^(0.5);
    
    % Cell total kinematics
    total_distance               (:,i) = tot_dist_out                        ;
    total_displacement           (:,i) = tot_disp_out                        ;
    directionality_ratio         (:,i) = tot_disp_out./tot_dist_out          ;

    % Cell forces
    interaction_time_x           (:,i) = int_force_out(:,1)                  ;
    interaction_time_y           (:,i) = int_force_out(:,2)                  ;
    contraction_time_x           (:,i) = con_force_out(:,1)                  ;
    contraction_time_y           (:,i) = con_force_out(:,2)                  ;
    protrusion_time_x            (:,i) = pro_force_out(:,1)                  ;
    protrusion_time_y            (:,i) = pro_force_out(:,2)                  ;
    % inf_craw_time_x              (:,i) = inf_crawl_out(:,1)                  ;
    % inf_craw_time_y              (:,i) = inf_crawl_out(:,2)                  ;
    tot_forc_time_x              (:,i) = tot_force_out(:,1)                  ;
    tot_forc_time_y              (:,i) = tot_force_out(:,2)                  ;

    

end

velocity_magn  =  displacement_time_magn/delta_time;

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

end

clear area_triangles displacements_solution infection_triang LANS_i MLSS_i ...
    sI_i sII_i triang_inf_id triang_areas stress_solution

%%
% PLOT AREA ...............................................................

% Area .........
% Exclude cells on the border for the analysis
prctile_area_step        = zeros(total_time_steps,3);

for i = 1 : total_time_steps

    id_cells_in            = ~location_time(:,i)              ;
    area_all_cells         = area_time(:,i)                   ;
    area_cells_in          = area_all_cells(id_cells_in)      ;
    prctile_area_step(i,1) = prctile(area_cells_in,25)        ;
    prctile_area_step(i,2) = prctile(area_cells_in,50)        ;
    prctile_area_step(i,3) = prctile(area_cells_in,75)        ;
    
end
hold on
plot(sim_time_min , prctile_area_step(1:end,1))
plot(sim_time_min , prctile_area_step(1:end,2))
plot(sim_time_min , prctile_area_step(1:end,3))
legend('25th perc','50th perc','75th perc');
title ('Area [\mum^2]' );
ylabel('Area [\mum^2]' );
xlabel('Simulation time [min]'         );
hold off
% nbins = 10;
% histogram(area_time(:,1),nbins);
if make_figure == true
    picname    = "prctile_area_" + case_name + ".png";
    saveas(gcf, picname);
end


% Area distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "area_" + case_name + ".mp4";
    vidObj           = VideoWriter(video_name,'MPEG-4');
    vidObj.FrameRate = 10;    % Set the frame rate to 10 frames per second
    open(vidObj);
end
edges = 200:20:800; % histogram edges (range)

for i = 1:total_time_steps

    data_to_analyse    = area_time(:,i)                   ;
    id_cells_in        = ~location_time(:,i)              ;
    hist_data          = data_to_analyse(id_cells_in)     ;

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([200, 800]         ) ;     ylim  ([0, 100]                );
    xlabel('area [\mum^2]'   ) ;     ylabel('frequency of appearence');

    txt = ['t: ' num2str(sim_time_min(i)) ' min'];
    title({['area distribution - ' txt], ''})

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
% PLOT POLYGON SIZE .......................................................

% Polygon size .........
% Exclude cells on the border for the analysis
median_pol_size_step        = zeros(total_time_steps,1);

for i = 1 : total_time_steps

    id_cells_in             = ~location_time(:,i)                  ;
    pol_size_all_cells      = polygon_size_time(:,i)               ;
    pol_size_cells_in       = pol_size_all_cells(id_cells_in)      ;
    
    median_pol_size_step(i) = median(pol_size_cells_in)            ;
    
end
% plot(median_pol_size_step(1:end))


% Polygon size distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "pol_size_" + case_name + ".mp4";
    vidObj           = VideoWriter(video_name,'MPEG-4');
    vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
    open(vidObj);
end


for i = 1:total_time_steps
    
    data_to_analyse = polygon_size_time(:,i)           ;
    id_cells_in     = ~location_time(:,i)              ;
    hist_data       = data_to_analyse(id_cells_in)     ;

    histogram(hist_data,'BinMethod','integers');
    xlim  ([2, 10]                    );    ylim  ([0, 400]                 );
    xlabel('polygon size'             );    ylabel('frequency of appearence');

    txt = ['t: ' num2str(sim_time_min(i)) ' min'];
    title({['polygon size distribution - ' txt], ''})

    % Add frame to video object
    frame       = getframe(gcf);
    if make_video == true
        writeVideo(vidObj,frame);
    end
    % Clear figure for next frame
    clf;
end

if make_video == true
    % Close video object
    close(vidObj);
end









%%
% PLOT PERIMETER ..........................................................

% Perimeter .........
% Exclude cells on the border for the analysis
prctile_perim_step          = zeros(total_time_steps,3);
for i = 1 : total_time_steps

    id_cells_in             = ~location_time(:,i)              ;
    perim_all_cells         = perimeter_time(:,i)              ;
    perim_cells_in          = perim_all_cells(id_cells_in)     ;
    prctile_perim_step(i,1) = prctile(perim_cells_in,25)       ;
    prctile_perim_step(i,2) = prctile(perim_cells_in,50)       ;
    prctile_perim_step(i,3) = prctile(perim_cells_in,75)       ;

end

hold on
plot(sim_time_min , prctile_perim_step(1:end,1))
plot(sim_time_min , prctile_perim_step(1:end,2))
plot(sim_time_min , prctile_perim_step(1:end,3))
legend('25th perc','50th perc','75th perc');
title ('Perimeter [\mum]' );
ylabel('Perimeter [\mum]' );
xlabel('Simulation time [min]'         );
hold off
% nbins = 10;
% histogram(perimeter_time(:,1),nbins);
if make_figure == true
    picname    = "prctile_perim_" + case_name + ".png";
    saveas(gcf, picname);
end




% Perimeter distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "perim_" + case_name + ".mp4";
    vidObj           = VideoWriter(video_name,'MPEG-4');
    vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
    open(vidObj);
end

edges = 30:2:80;  % histogram edges (range)
% Loop over frames
for i = 1:total_time_steps

    data_to_analyse   = perimeter_time(:,i)             ;
    id_cells_in       = ~location_time(:,i)             ;
    hist_data         = data_to_analyse(id_cells_in)    ;

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([40, 80]               );   ylim  ([0, 350]                 );
    xlabel('perimeter [\mum]'      );   ylabel('frequency of appearence');

    txt = ['t: ' num2str(sim_time_min(i)) ' min'];
    title({['Perimeter distribution - ' txt], ''})

    % Add frame to video object
    frame       = getframe(gcf);
    if make_video == true
        writeVideo(vidObj,frame);
    end
    
    % Clear figure for next frame
    clf;
end

if make_video == true
    % Close video object
    close(vidObj);
end






%%
% PLOT CELL SHAPE .........................................................

% Cell shape .........
% Exclude cells on the border for the analysis
prctile_cell_shape_step         = zeros(total_time_steps,3);
for i = 1 : total_time_steps

    id_cells_in                  = ~location_time(:,i)              ;
    cshape_all_cells             = cell_shape_time(:,i)             ;
    cshape_cells_in              = cshape_all_cells(id_cells_in)    ;
    prctile_cell_shape_step(i,1) = prctile(cshape_cells_in,25)      ;
    prctile_cell_shape_step(i,2) = prctile(cshape_cells_in,50)      ;
    prctile_cell_shape_step(i,3) = prctile(cshape_cells_in,75)      ;

end

hold on
plot(sim_time_min , prctile_cell_shape_step(1:end,1))
plot(sim_time_min , prctile_cell_shape_step(1:end,2))
plot(sim_time_min , prctile_cell_shape_step(1:end,3))
legend('25th perc','50th perc','75th perc');
title ('Cell shape index [-]' );
ylabel('Cell shape index [-]' );
xlabel('Simulation time [min]'         );
hold off
% nbins = 10;
% histogram(cell_shape_index(:,1),nbins);
if make_figure == true
    picname    = "prctile_cshape_" + case_name + ".png";
    saveas(gcf, picname);
end



% Cell shape distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "cshape_" + case_name + ".mp4";
    vidObj           = VideoWriter(video_name,'MPEG-4');
    vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
    open(vidObj);
end

edges = 3.51:0.05:5; % histogram edges (range)

for i = 1:total_time_steps

    data_to_analyse = cell_shape_time(:,i)           ;
    id_cells_in     = ~location_time(:,i)            ;
    hist_data       = data_to_analyse(id_cells_in)   ;

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([3.4, 5]               );    ylim  ([0, 400]                 );
    xlabel('Cell shape [-]'       );    ylabel('frequency of appearence');

    txt = ['t: ' num2str(sim_time_min(i)) ' min'];
    title({['Cell shape distribution - ' txt], ''})

    % Print shape index limit
    hold on;
    xline = 3.81;
    line([xline, xline], ylim, 'Color', 'r');
    hold off

    % Add frame to video object
    frame       = getframe(gcf);
    if make_video == true
        writeVideo(vidObj,frame);
    end
    
    % Clear figure for next frame
    clf;
end

if make_video == true
    % Close video object
    close(vidObj);
end









%%
% PLOT DISPLACEMENT VECTOR ................................................

% Displacement vector .........
% Exclude cells on the border for the analysis
prctile_disp_vec_magn_step     = zeros(total_time_steps,3);
for i = 1 : total_time_steps

    id_cells_in                     = ~location_time(:,i)                ;
    disp_mag_all_cells              = displacement_time_magn(:,i)        ;
    disp_mag_cells_in               = disp_mag_all_cells(id_cells_in)    ;
    prctile_disp_vec_magn_step(i,1) = prctile(disp_mag_cells_in,25)      ;
    prctile_disp_vec_magn_step(i,2) = prctile(disp_mag_cells_in,50)      ;
    prctile_disp_vec_magn_step(i,3) = prctile(disp_mag_cells_in,75)      ;

end

hold on
plot(sim_time_min(2:end) , prctile_disp_vec_magn_step(2:end,1))
plot(sim_time_min(2:end) , prctile_disp_vec_magn_step(2:end,2))
plot(sim_time_min(2:end) , prctile_disp_vec_magn_step(2:end,3))
legend('25th perc','50th perc','75th perc');
title ('Cell displacement [\mum]' );
ylabel('Cell displacement [\mum]' );
xlabel('Simulation time [min]'         );
hold off
% nbins = 10;
% histogram(displacement_time_magn(:,2),nbins);
if make_figure == true
    picname    = "prctile_disp_vec_" + case_name + ".png";
    saveas(gcf, picname);
end



% Displacement vector distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "disp_magn_" + case_name + ".mp4";
    vidObj           = VideoWriter(video_name,'MPEG-4');
    vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
    open(vidObj);
end

edges                = 0:0.02:2; % histogram edges (range)
% Loop over frames
for i = 2:total_time_steps
    
    data_to_analyse = displacement_time_magn(:,i)   ;
    id_cells_in     = ~location_time(:,i)           ;
    hist_data       = data_to_analyse(id_cells_in)  ;

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([0, 2]                              );   ylim  ([0, 150]                  ); 
    xlabel('Cell displacement magnitude [\mum]');   ylabel('frequency of appearence');

    txt = ['t: ' num2str(sim_time_min(i)) ' min'];
    title({['Cell displacement magn distribution - ' txt], ''})
    

    % Add frame to video object
    frame       = getframe(gcf);
    if make_video == true
        writeVideo(vidObj,frame);
    end
    
    % Clear figure for next frame
    clf;
end

if make_video == true
    % Close video object
    close(vidObj);
end






%%
% PLOT VELOCITY VECTOR ................................................

% Velocity vector .........
% Exclude cells on the border for the analysis
prctile_velocity_magn_step     = zeros(total_time_steps,3);
for i = 1 : total_time_steps

    id_cells_in                     = ~location_time(:,i)                ;
    vel_mag_all_cells               = velocity_magn(:,i)                 ;
    vel_mag_cells_in                = vel_mag_all_cells(id_cells_in)     ;
    prctile_velocity_magn_step(i,1) = prctile(vel_mag_cells_in,25)       ;
    prctile_velocity_magn_step(i,2) = prctile(vel_mag_cells_in,50)       ;
    prctile_velocity_magn_step(i,3) = prctile(vel_mag_cells_in,75)       ;

end

hold on
plot(sim_time_min(2:end) , prctile_velocity_magn_step(2:end,1))
plot(sim_time_min(2:end) , prctile_velocity_magn_step(2:end,2))
plot(sim_time_min(2:end) , prctile_velocity_magn_step(2:end,3))
legend('25th perc','50th perc','75th perc');
title ('Cell nuclei velocity [\mum/min]' );
ylabel('Cell velocity [\mum/min]' );
xlabel('Simulation time [min]'         );
hold off
% nbins = 10;
% histogram(displacement_time_magn(:,2),nbins);
if make_figure == true
    picname    = "prctile_veloc_" + case_name + ".png";
    saveas(gcf, picname);
end



% Cell velocity vector distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "vel_magn_" + case_name + ".mp4";
    vidObj           = VideoWriter(video_name,'MPEG-4');
    vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
    open(vidObj);
end

edges                = 0:0.02:2; % histogram edges (range)
% Loop over frames
for i = 2:total_time_steps
    
    data_to_analyse = velocity_magn(:,i)            ;
    id_cells_in     = ~location_time(:,i)           ;
    hist_data       = data_to_analyse(id_cells_in)  ;

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([0, 2]                              );   ylim  ([0, 150]                  ); 
    xlabel('Cell velocity magnitude [\mum/min]');   ylabel('frequency of appearence');
    title ('Cell velocity magnitude'           );

    txt = ['t: ' num2str(sim_time_min(i)) ' min'];
    title({['Cell velocity distribution - ' txt], ''})
    

    % Add frame to video object
    frame       = getframe(gcf);
    if make_video == true
        writeVideo(vidObj,frame);
    end
    
    % Clear figure for next frame
    clf;
end

if make_video == true
    % Close video object
    close(vidObj);
end






%%
% PLOT TOTAL DISTANCE ................................................

% Total distance walked by the cell .........
% Exclude cells on the border for the analysis
prctile_tot_dist_step     = zeros(total_time_steps,3);
for i = 1 : total_time_steps

    id_cells_in                     = ~location_time(:,i)                ;
    tot_dis_all_cells               = total_distance(:,i)                ;
    tot_dis_cells_in                = tot_dis_all_cells(id_cells_in)     ;
    prctile_tot_dist_step(i,1) = prctile(tot_dis_cells_in,25)       ;
    prctile_tot_dist_step(i,2) = prctile(tot_dis_cells_in,50)       ;
    prctile_tot_dist_step(i,3) = prctile(tot_dis_cells_in,75)       ;

end

hold on
plot(sim_time_min(2:end) , prctile_tot_dist_step(2:end,1))
plot(sim_time_min(2:end) , prctile_tot_dist_step(2:end,2))
plot(sim_time_min(2:end) , prctile_tot_dist_step(2:end,3))
legend('25th perc','50th perc','75th perc');
title ('Cell total distance [\mum]' );
ylabel('Cell total distance [\mum]' );
xlabel('Simulation time [min]'         );
hold off
% nbins = 10;
% histogram(displacement_time_magn(:,2),nbins);
if make_figure == true
    picname    = "prctile_tot_dist_" + case_name + ".png";
    saveas(gcf, picname);
end



% Cell total distance distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "tot_dis_" + case_name + ".mp4";
    vidObj           = VideoWriter(video_name,'MPEG-4');
    vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
    open(vidObj);
end

edges                = 0:0.5:100; % histogram edges (range)
% Loop over frames
for i = 2:total_time_steps
    
    data_to_analyse = total_distance(:,i)           ;
    id_cells_in     = ~location_time(:,i)           ;
    hist_data       = data_to_analyse(id_cells_in)  ;

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([0, 100]                              );   ylim  ([0, 100]                  ); 
    xlabel('Cell total distance [\mum]'        );   ylabel('frequency of appearence');

    txt = ['t: ' num2str(sim_time_min(i)) ' min'];
    title({['Cell total distance - ' txt], ''})
    

    % Add frame to video object
    frame       = getframe(gcf);
    if make_video == true
        writeVideo(vidObj,frame);
    end
    
    % Clear figure for next frame
    clf;
end

if make_video == true
    % Close video object
    close(vidObj);
end




%%
% PLOT TOTAL DISPLACEMENT ................................................

% Total displacement walked by the cell .........
% Exclude cells on the border for the analysis
prctile_tot_disp_step     = zeros(total_time_steps,3);
for i = 1 : total_time_steps

    id_cells_in                     = ~location_time(:,i)                 ;
    tot_disp_all_cells              = total_displacement(:,i)             ;
    tot_disp_cells_in               = tot_disp_all_cells(id_cells_in)     ;
    prctile_tot_disp_step(i,1) = prctile(tot_disp_cells_in,25)       ;
    prctile_tot_disp_step(i,2) = prctile(tot_disp_cells_in,50)       ;
    prctile_tot_disp_step(i,3) = prctile(tot_disp_cells_in,75)       ;

end

hold on
plot(sim_time_min(2:end) , prctile_tot_disp_step(2:end,1))
plot(sim_time_min(2:end) , prctile_tot_disp_step(2:end,2))
plot(sim_time_min(2:end) , prctile_tot_disp_step(2:end,3))
legend('25th perc','50th perc','75th perc');
title ('Cell total displacement [\mum]' );
ylabel('Cell total displacement [\mum]' );
xlabel('Simulation time [min]'         );
hold off
% nbins = 10;
% histogram(displacement_time_magn(:,2),nbins);
if make_figure == true
    picname    = "prctile_tot_disp_" + case_name + ".png";
    saveas(gcf, picname);
end



% Cell total displacement distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "tot_disp_" + case_name + ".mp4";
    vidObj           = VideoWriter(video_name,'MPEG-4');
    vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
    open(vidObj);
end

edges                = 0:0.5:50; % histogram edges (range)
% Loop over frames
for i = 2:total_time_steps
    
    data_to_analyse = total_displacement(:,i)       ;
    id_cells_in     = ~location_time(:,i)           ;
    hist_data       = data_to_analyse(id_cells_in)  ;

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([0, 50]                              );    ylim  ([0, 150]                 ); 
    xlabel('Cell total displacement [\mum]'      );   ylabel('frequency of appearence');

    txt = ['t: ' num2str(sim_time_min(i)) ' min'];
    title({['Cell total displacement - ' txt], ''})
    

    % Add frame to video object
    frame       = getframe(gcf);
    if make_video == true
        writeVideo(vidObj,frame);
    end
    
    % Clear figure for next frame
    clf;
end

if make_video == true
    % Close video object
    close(vidObj);
end






%%
% PLOT DIRECTIONALITY RATIO................................................

% Directionality ratio .........
% Exclude cells on the border for the analysis
prctile_dir_ratio_step     = zeros(total_time_steps,3);
for i = 1 : total_time_steps

    id_cells_in                     = ~location_time(:,i)                 ;
    dir_ratio_all_cells             = directionality_ratio(:,i)           ;
    dir_ratio_cells_in              = dir_ratio_all_cells(id_cells_in)    ;
    prctile_dir_ratio_step(i,1)     = prctile(dir_ratio_cells_in,25)      ;
    prctile_dir_ratio_step(i,2)     = prctile(dir_ratio_cells_in,50)      ;
    prctile_dir_ratio_step(i,3)     = prctile(dir_ratio_cells_in,75)      ;

end

hold on
plot(sim_time_min(3:end) , prctile_dir_ratio_step(3:end,1))
plot(sim_time_min(3:end) , prctile_dir_ratio_step(3:end,2))
plot(sim_time_min(3:end) , prctile_dir_ratio_step(3:end,3))
legend('25th perc','50th perc','75th perc');
title ('Directionality ratio [-]' );
ylabel('Directionality ratio [-]' );
xlabel('Simulation time [min]'    );
hold off
% nbins = 10;
% histogram(displacement_time_magn(:,2),nbins);
if make_figure == true
    picname    = "prctile_dir_rat_" + case_name + ".png";
    saveas(gcf, picname);
end



% Directionality ratio distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "dir_ratio_" + case_name + ".mp4";
    vidObj           = VideoWriter(video_name,'MPEG-4');
    vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
    open(vidObj);
end

edges                = 0:0.02:1.1; % histogram edges (range)
% Loop over frames
for i = 3:total_time_steps
    
    data_to_analyse = directionality_ratio(:,i)     ;
    id_cells_in     = ~location_time(:,i)           ;
    hist_data       = data_to_analyse(id_cells_in)  ;

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([0, 1.1]                              );    ylim  ([0, 150]                 ); 
    xlabel('Directionality ratio [-]'            );    ylabel('frequency of appearence');

    txt = ['t: ' num2str(sim_time_min(i)) ' min'];
    title({['Cell directionality ratio - ' txt], ''})
    

    % Add frame to video object
    frame       = getframe(gcf);
    if make_video == true
        writeVideo(vidObj,frame);
    end
    
    % Clear figure for next frame
    clf;
end

if make_video == true
    % Close video object
    close(vidObj);
end



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
    xlim  ([-70, 120]         ) ;     ylim  ([0, 500]                  );
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
    xlim  ([-140, 70]          ) ;     ylim  ([0, 500]                  );
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
    xlim  ([-140, 140]        ) ;     ylim  ([0, 500]                  );
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
    xlim  ([-10, 140]          ) ;     ylim  ([0, 1000]                );
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
% Boxplots

% Area
data_vector          = reshape(area_time, [], 1);
mean_area_cell       = mean(area_time,2); % average each cell over time

boxplot(area_time);
boxplot(data_vector);
ylabel('Area [\mum^2]');
if make_figure == true
    picname    = "boxplot_area_" + case_name + ".png";
    saveas(gcf, picname);
end
boxplot(mean_area_cell);

% Polygon size
% data_vector          = reshape(polygon_size_time, [], 1);
% median_pol_size_cell = median(polygon_size_time,2); % median of each cell over time
% 
% boxplot(polygon_size_time);
% boxplot(data_vector);
% boxplot(median_pol_size_cell);

% Perimeter
data_vector          = reshape(perimeter_time, [], 1);
mean_perim_cell      = mean(perimeter_time,2); % average each cell over time

boxplot(perimeter_time);
boxplot(data_vector);
ylabel('perimeter [\mum]');
if make_figure == true
    picname    = "boxplot_perim_" + case_name + ".png";
    saveas(gcf, picname);
end
boxplot(mean_perim_cell);

% Cell shape
data_vector          = reshape(cell_shape_time, [], 1);
mean_cshape_cell     = mean(cell_shape_time,2); % average each cell over time

boxplot(cell_shape_time);
boxplot(data_vector);
ylabel('Cell shape [-]');
if make_figure == true
    picname    = "boxplot_cshape_" + case_name + ".png";
    saveas(gcf, picname);
end
boxplot(mean_cshape_cell);

% Cell displacement vector
data_vector          = reshape(displacement_time_magn, [], 1);
mean_cell_displ_cell = mean(displacement_time_magn,2); % average each cell over time

boxplot(displacement_time_magn(2:end,:));      % skip first step
boxplot(data_vector(number_of_cells+1:end));   % skip cells from the first step
ylabel('Cell displacements [\mum]');
if make_figure == true
    picname    = "boxplot_celldisp_" + case_name + ".png";
    saveas(gcf, picname);
end
boxplot(mean_cell_displ_cell(2:end,:));

% Cell velocity
data_vector          = reshape(velocity_magn, [], 1);
mean_cell_veloc_cell = mean(velocity_magn,2); % average each cell over time

boxplot(velocity_magn(2:end,:));      % skip first step
boxplot(data_vector(number_of_cells+1:end));   % skip cells from the first step
ylabel('Cell velocity [\mum/min]');
if make_figure == true
    picname    = "boxplot_cellveloc_" + case_name + ".png";
    saveas(gcf, picname);
end
boxplot(mean_cell_veloc_cell(2:end,:));

% Cell total distance
data_vector          = reshape(total_distance, [], 1);
mean_tot_dist_cell   = mean(total_distance,2); % average each cell over time

boxplot(total_distance(2:end,:));              % skip first step
boxplot(data_vector(number_of_cells+1:end));   % skip cells from the first step
ylabel('Total distance [\mum]');
if make_figure == true
    picname    = "boxplot_totdist_" + case_name + ".png";
    saveas(gcf, picname);
end
boxplot(mean_tot_dist_cell(2:end,:));

% Cell total displacement
data_vector          = reshape(total_displacement, [], 1);
mean_tot_disp_cell   = mean(total_displacement,2); % average each cell over time

boxplot(total_displacement(2:end,:));              % skip first step
boxplot(data_vector(number_of_cells+1:end));       % skip cells from the first step
ylabel('Total displacement [\mum^2]');
if make_figure == true
    picname    = "boxplot_totdisp_" + case_name + ".png";
    saveas(gcf, picname);
end
boxplot(mean_tot_disp_cell(2:end,:));

% Cell directionality ratio


% Save postrpocess parameters ............................................
filename = "var_"+ case_name;
save(filename, 'area_time', 'polygon_size_time', 'perimeter_time', 'cell_shape_time'           ...
             , 'location_time', 'infection_time', 'prctile_area_step','median_pol_size_step'   ...
             , 'prctile_perim_step', 'prctile_cell_shape_step', 'total_time_steps'             ...
             , 'number_of_cells', 'displacement_time_magn','prctile_disp_vec_magn_step'        ...
             , 'prctile_velocity_magn_step', 'prctile_tot_dist_step', 'prctile_tot_disp_step'  ...
             , 'prctile_dir_ratio_step');
