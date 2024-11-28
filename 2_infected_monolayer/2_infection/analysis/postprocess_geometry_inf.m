%% MONOLAYER POSTPROCESS: GEOMETRY

clear all;
clc;
close all;

% Initialize structures

% Load first step to retrieve some parameters
load ./results/raw/data_00001.mat

case_name        = 'case_a2';
make_video       = false;
total_time_steps = 250;
number_of_cells  = numel(id_out);

% Cell geometry
area_time                    = zeros(number_of_cells,total_time_steps);
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

% Load parameters
for i = 1:total_time_steps

    file_name = sprintf('./results/raw/data_%05i.mat', i); 
    load(file_name)

    % Cell geometry
    area_time                    (:,i) = area_out                          ;
    polygon_size_time            (:,i) = pol_sizes_out                     ;
    perimeter_time               (:,i) = perimeter_out                     ;
    cell_shape_time              (:,i) = perimeter_out ./ (area_out).^(0.5);

    % Cell state
    location_time                (:,i) = location_out                      ;
    infection_time               (:,i) = infected_out                      ;
    
    % Cell kinematics
    displacement_time_x          (:,i) = disp_vec_out(:,1);
    displacement_time_y          (:,i) = disp_vec_out(:,2);
    displacement_time_magn       (:,i) = (displacement_time_x(:,i).^2  + ...
                                          displacement_time_y(:,i).^2).^(0.5);

end

clear area_out pol_sizes_out perimeter_out location_out disp_vec_out   ...
      velo_vec_out tot_dist_out tot_disp_out t_dis_ve_out id_out       ...
      init_pos_out radius_out curr_pos_out tot_force_out int_force_out ...
      con_force_out rand_walk_out;










% PLOT AREA ...............................................................

% Area .........
% Exclude cells on the border for the analysis
prctile_area_step         = zeros(total_time_steps,3);
prctile_area_inf_step     = zeros(total_time_steps,3);
prctile_area_uninf_step   = zeros(total_time_steps,3);

for i = 1 : total_time_steps

    id_cells_in                  = ~location_time(:,i)                   ;
    id_cells_out                 =  logical(location_time(:,i))          ;
    area_all_cells               = area_time(:,i)                        ;
    area_cells_in                = area_all_cells(id_cells_in)           ;
    prctile_area_step(i,1)       = prctile(area_cells_in,25)             ;
    prctile_area_step(i,2)       = prctile(area_cells_in,50)             ;
    prctile_area_step(i,3)       = prctile(area_cells_in,75)             ;

    % Infected versus non-infected
    % Infected:
    id_cells_infected            = logical(infection_time(:,i))          ;
    area_cells_infected          = area_all_cells(id_cells_infected)     ;
    prctile_area_inf_step(i,1)   = prctile(area_cells_infected,25)       ;
    prctile_area_inf_step(i,2)   = prctile(area_cells_infected,50)       ;
    prctile_area_inf_step(i,3)   = prctile(area_cells_infected,75)       ;
    % Uninfected: remove cells at the boundary
    id_cells_uninfected          = ~infection_time(:,i)                  ;
    id_cells_uninfected          = id_cells_uninfected-id_cells_out      ; % remove cells at the boundary
    id_cells_uninfected          = logical(id_cells_uninfected)          ; 
    area_cells_uninfected        = area_all_cells(id_cells_uninfected)   ;
    prctile_area_uninf_step(i,1) = prctile(area_cells_uninfected,25)     ;
    prctile_area_uninf_step(i,2) = prctile(area_cells_uninfected,50)     ;
    prctile_area_uninf_step(i,3) = prctile(area_cells_uninfected,75)     ;
    
end

hold on
plot(prctile_area_step(1:end,1))
plot(prctile_area_step(1:end,2))
plot(prctile_area_step(1:end,3))
legend('25th perc','50th perc','75th perc');
title ('Area [\mum^2]' );
ylabel('Area [\mum^2]' );
xlabel('Steps'         );
hold off
% nbins = 10;
% histogram(area_time(:,1),nbins);

% Plot area infected vs uninfected
figure()
hold on
name1 = "25th perc uninf "+ case_name; name4 = "25th perc inf "+ case_name       ;
name2 = "50th perc uninf "+ case_name; name5 = "50th perc inf "+ case_name       ;
name3 = "75th perc uninf "+ case_name; name6 = "75th perc inf "+ case_name       ;
p1    = plot(prctile_area_uninf_step(1:end,1), 'LineStyle', ':', 'Color', 'blue');
p2    = plot(prctile_area_uninf_step(1:end,2),                   'Color', 'blue');
p3    = plot(prctile_area_uninf_step(1:end,3), 'LineStyle', ':', 'Color', 'blue');
p4    = plot(prctile_area_inf_step  (1:end,1), 'LineStyle', ':', 'Color', 'red' );
p5    = plot(prctile_area_inf_step  (1:end,2),                   'Color', 'red' );
p6    = plot(prctile_area_inf_step  (1:end,3), 'LineStyle', ':', 'Color', 'red' );

lgd = legend([p2, p5], {name2, name5});
title ('Area [\mum^2]' );
ylabel('Area [\mum^2]' );
xlabel('Steps'         );
hold off
picname    = "area_inf_vs_uninf_comparison_" + case_name + ".png";
saveas(gcf, picname);






% Area distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "area_" + case_name + ".avi";
    vidObj           = VideoWriter(video_name);
    vidObj.FrameRate = 10;    % Set the frame rate to 10 frames per second
    open(vidObj);
end
edges = 120:20:800; % histogram edges (range)

for i = 1:total_time_steps

    data_to_analyse    = area_time(:,i)                   ;
    id_cells_in        = ~location_time(:,i)              ;
    hist_data          = data_to_analyse(id_cells_in)     ;

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([0, 800]          ) ;     ylim  ([0, 100]                 );
    xlabel('area [\mum^2]'   ) ;     ylabel('frequency of appearence');
    title ('area distribution');

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




% Area infected vs uninfected
video_name       = "area_infvsuninf" + case_name + ".avi";
vidObj           = VideoWriter(video_name);
vidObj.FrameRate = 10;    % Set the frame rate to 10 frames per second
open(vidObj);

edges = 120:20:800; % histogram edges (range)

for i = 1:total_time_steps

    hold on

    id_cells_out                 =  logical(location_time(:,i))          ;
    area_all_cells               = area_time(:,i)                        ;
    % Infected versus non-infected
    % Infected:
    id_cells_infected            = logical(infection_time(:,i))          ;
    area_cells_infected          = area_all_cells(id_cells_infected)     ;
    hist_data                    = area_cells_infected                   ;
    histogram(hist_data, 'BinEdges', edges);
    
    % Uninfected: remove cells at the boundary
    id_cells_uninfected          = ~infection_time(:,i)                  ;
    id_cells_uninfected          = id_cells_uninfected-id_cells_out      ; % remove cells at the boundary
    id_cells_uninfected          = logical(id_cells_uninfected)          ; 
    area_cells_uninfected        = area_all_cells(id_cells_uninfected)   ;
    hist_data                    = area_cells_uninfected;
    histogram(hist_data, 'BinEdges', edges);

    xlim  ([0, 800]          ) ;     ylim  ([0, 100]                 );
    xlabel('area [\mum^2]'   ) ;     ylabel('frequency of appearence');
    title ('area distribution');     legend('infected','uninfected')  ;

    hold off

    % Add frame to video object
    frame       = getframe(gcf);    
    writeVideo(vidObj,frame);
    % Clear figure for next frame
    clf;
end

% Close video object
close(vidObj);









% PLOT POLYGON SIZE .......................................................

% Polygon size .........
% Exclude cells on the border for the analysis
median_pol_size_step        = zeros(total_time_steps,1);
median_pol_size_inf_step    = zeros(total_time_steps,1);
median_pol_size_uninf_step  = zeros(total_time_steps,1);


for i = 1 : total_time_steps

    id_cells_in                     = ~location_time(:,i)                     ;
    id_cells_out                    =  logical(location_time(:,i))            ;
    pol_size_all_cells              = polygon_size_time(:,i)                  ;
    pol_size_cells_in               = pol_size_all_cells(id_cells_in)         ; 
    median_pol_size_step(i)         = median(pol_size_cells_in)               ;
    
    % Infected versus non-infected
    % Infected:
    id_cells_infected               = logical(infection_time(:,i))            ;
    pol_size_infected               = pol_size_all_cells(id_cells_infected)   ;
    median_pol_size_inf_step(i)     = median(pol_size_infected)               ;
   
    % Uninfected: remove cells at the boundary
    id_cells_uninfected             = ~infection_time(:,i)                    ;
    id_cells_uninfected             = id_cells_uninfected-id_cells_out        ; % remove cells at the boundary
    id_cells_uninfected             = logical(id_cells_uninfected)            ; 
    pol_size_uninfected             = pol_size_all_cells(id_cells_uninfected) ;
    median_pol_size_uninf_step(i,1) = median(pol_size_uninfected)             ;
   
end

% plot(median_pol_size_step(1:end))
% plot(median_pol_size_inf_step(1:end))
% plot(median_pol_size_uninf_step(1:end))




% Polygon size distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "pol_size" + case_name + ".avi";
    vidObj           = VideoWriter(video_name);
    vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
    open(vidObj);
end


for i = 1:total_time_steps
    
    data_to_analyse = polygon_size_time(:,i)           ;
    id_cells_in     = ~location_time(:,i)              ;
    hist_data       = data_to_analyse(id_cells_in)     ;

    histogram(hist_data,'BinMethod','integers');
    xlim  ([2, 10]                    );    ylim  ([0, 600]                 );
    xlabel('polygon size'             );    ylabel('frequency of appearence');
    title ('polygon size distribution');

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



% Pol size infected vs uninfected
video_name       = "pol_size_infvsuninf" + case_name + ".avi";
vidObj           = VideoWriter(video_name);
vidObj.FrameRate = 10;    % Set the frame rate to 10 frames per second
open(vidObj);

for i = 1:total_time_steps

    hold on

    id_cells_out                 =  logical(location_time(:,i))            ;
    pol_size_all_cells           = polygon_size_time(:,i)                  ;
    % Infected versus non-infected
    % Infected:
    id_cells_infected            = logical(infection_time(:,i))            ;
    pol_size_cells_infected      = pol_size_all_cells(id_cells_infected)   ;
    hist_data                    = pol_size_cells_infected                 ;
    histogram(hist_data,'BinMethod','integers');
    
    % Uninfected: remove cells at the boundary
    id_cells_uninfected          = ~infection_time(:,i)                    ;
    id_cells_uninfected          = id_cells_uninfected-id_cells_out        ; % remove cells at the boundary
    id_cells_uninfected          = logical(id_cells_uninfected)            ; 
    pol_size_cells_uninfected    = pol_size_all_cells(id_cells_uninfected) ;
    hist_data                    = pol_size_cells_uninfected;
    histogram(hist_data,'BinMethod','integers');

    xlim  ([2, 10]                    );    ylim  ([0, 400]                 );
    xlabel('polygon size'             );    ylabel('frequency of appearence');
    title ('polygon size distribution');    legend('infected','uninfected'  );

    hold off

    % Add frame to video object
    frame       = getframe(gcf);    
    writeVideo(vidObj,frame);
    % Clear figure for next frame
    clf;


end

% Close video object
close(vidObj);










% PLOT PERIMETER ..........................................................

% Perimeter .........
% Exclude cells on the border for the analysis
prctile_perim_step         = zeros(total_time_steps,3);
prctile_perim_inf_step     = zeros(total_time_steps,3);
prctile_perim_uninf_step   = zeros(total_time_steps,3);

for i = 1 : total_time_steps

    id_cells_in                   = ~location_time(:,i)              ;
    perim_all_cells               = perimeter_time(:,i)              ;
    perim_cells_in                = perim_all_cells(id_cells_in)     ;
    prctile_perim_step(i,1)       = prctile(perim_cells_in,25)       ;
    prctile_perim_step(i,2)       = prctile(perim_cells_in,50)       ;
    prctile_perim_step(i,3)       = prctile(perim_cells_in,75)       ;

    % Infected versus non-infected
    % Infected:
    id_cells_infected             = logical(infection_time(:,i))          ;
    peri_cells_infected           = perim_all_cells(id_cells_infected)    ;
    prctile_perim_inf_step(i,1)   = prctile(peri_cells_infected,25)       ;
    prctile_perim_inf_step(i,2)   = prctile(peri_cells_infected,50)       ;
    prctile_perim_inf_step(i,3)   = prctile(peri_cells_infected,75)       ;
    % Uninfected: remove cells at the boundary
    id_cells_uninfected           = ~infection_time(:,i)                  ;
    id_cells_uninfected           = id_cells_uninfected-id_cells_out      ; % remove cells at the boundary
    id_cells_uninfected           = logical(id_cells_uninfected)          ; 
    peri_cells_uninfected         = perim_all_cells(id_cells_uninfected)  ;
    prctile_perim_uninf_step(i,1) = prctile(peri_cells_uninfected,25)     ;
    prctile_perim_uninf_step(i,2) = prctile(peri_cells_uninfected,50)     ;
    prctile_perim_uninf_step(i,3) = prctile(peri_cells_uninfected,75)     ;

end

hold on
plot(prctile_perim_step(1:end,1))
plot(prctile_perim_step(1:end,2))
plot(prctile_perim_step(1:end,3))
legend('25th perc','50th perc','75th perc');
title ('Perimeter [\mum]' );
ylabel('Perimeter [\mum]' );
xlabel('Steps'         );
hold off
% nbins = 10;
% histogram(perimeter_time(:,1),nbins);

% Plot perimeter infected vs uninfected
figure()
hold on
name1 = "25th perc uninf "+ case_name; name4 = "25th perc inf "+ case_name       ;
name2 = "50th perc uninf "+ case_name; name5 = "50th perc inf "+ case_name       ;
name3 = "75th perc uninf "+ case_name; name6 = "75th perc inf "+ case_name       ;
p1    = plot(prctile_perim_uninf_step(1:end,1), 'LineStyle', ':', 'Color', 'blue');
p2    = plot(prctile_perim_uninf_step(1:end,2),                   'Color', 'blue');
p3    = plot(prctile_perim_uninf_step(1:end,3), 'LineStyle', ':', 'Color', 'blue');
p4    = plot(prctile_perim_inf_step  (1:end,1), 'LineStyle', ':', 'Color', 'red' );
p5    = plot(prctile_perim_inf_step  (1:end,2),                   'Color', 'red' );
p6    = plot(prctile_perim_inf_step  (1:end,3), 'LineStyle', ':', 'Color', 'red' );

lgd = legend([p2, p5], {name2, name5});
title ('Perimeter [\mum]' );
ylabel('Perimeter [\mum]' );
xlabel('Steps'         );
hold off
picname    = "perim_inf_vs_uninf_comparison_" + case_name + ".png";
saveas(gcf, picname);







% Perimeter distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "perim_" + case_name + ".avi";
    vidObj           = VideoWriter(video_name);
    vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
    open(vidObj);
end

edges = 50:2:120;  % histogram edges (range)
% Loop over frames
for i = 1:total_time_steps

    data_to_analyse   = perimeter_time(:,i)             ;
    id_cells_in       = ~location_time(:,i)             ;
    hist_data         = data_to_analyse(id_cells_in)    ;

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([40, 120]               );   ylim  ([0, 300]                 );
    xlabel('perimeter [\mum]'      );   ylabel('frequency of appearence');
    title ('Perimeter distribution');

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


% Perimeter infected vs uninfected
video_name       = "perim_infvsuninf" + case_name + ".avi";
vidObj           = VideoWriter(video_name);
vidObj.FrameRate = 10;    % Set the frame rate to 10 frames per second
open(vidObj);

edges = 50:2:120; % histogram edges (range)

for i = 1:total_time_steps

    hold on

    id_cells_out                 =  logical(location_time(:,i))          ;
    perim_all_cells              = perimeter_time(:,i)                   ;
    % Infected versus non-infected
    % Infected:
    id_cells_infected            = logical(infection_time(:,i))          ;
    peri_cells_infected          = perim_all_cells(id_cells_infected)    ;
    hist_data                    = peri_cells_infected                   ;
    histogram(hist_data, 'BinEdges', edges);
    
    % Uninfected: remove cells at the boundary
    id_cells_uninfected          = ~infection_time(:,i)                  ;
    id_cells_uninfected          = id_cells_uninfected-id_cells_out      ; % remove cells at the boundary
    id_cells_uninfected          = logical(id_cells_uninfected)          ; 
    peri_cells_uninfected        = perim_all_cells(id_cells_uninfected)  ;
    hist_data                    = peri_cells_uninfected;
    histogram(hist_data, 'BinEdges', edges);

    xlim  ([40, 120]          );     ylim  ([0, 100]                 );
    xlabel('perimeter [\mum]'   ) ;     ylabel('frequency of appearence');
    title ('perimeter distribution');     legend('infected','uninfected')  ;

    hold off

    % Add frame to video object
    frame       = getframe(gcf);    
    writeVideo(vidObj,frame);
    % Clear figure for next frame
    clf;
end

% Close video object
close(vidObj);









% PLOT CELL SHAPE .........................................................

% Cell shape .........
% Exclude cells on the border for the analysis
prctile_cell_shape_step         = zeros(total_time_steps,3);
prctile_cell_shape_inf_step     = zeros(total_time_steps,3);
prctile_cell_shape_uninf_step   = zeros(total_time_steps,3);

for i = 1 : total_time_steps

    id_cells_in                        = ~location_time(:,i)                   ;
    cshape_all_cells                   = cell_shape_time(:,i)                  ;
    cshape_cells_in                    = cshape_all_cells(id_cells_in)         ;
    prctile_cell_shape_step(i,1)       = prctile(cshape_cells_in,25)           ;
    prctile_cell_shape_step(i,2)       = prctile(cshape_cells_in,50)           ;
    prctile_cell_shape_step(i,3)       = prctile(cshape_cells_in,75)           ;

    % Infected versus non-infected
    % Infected:
    id_cells_infected                  = logical(infection_time(:,i))          ;
    cshape_cells_infected              = cshape_all_cells(id_cells_infected)   ;
    prctile_cell_shape_inf_step(i,1)   = prctile(cshape_cells_infected,25)     ;
    prctile_cell_shape_inf_step(i,2)   = prctile(cshape_cells_infected,50)     ;
    prctile_cell_shape_inf_step(i,3)   = prctile(cshape_cells_infected,75)     ;
    % Uninfected: remove cells at the boundary
    id_cells_uninfected                = ~infection_time(:,i)                  ;
    id_cells_uninfected                = id_cells_uninfected-id_cells_out      ; % remove cells at the boundary
    id_cells_uninfected                = logical(id_cells_uninfected)          ; 
    cshape_cells_uninfected            = cshape_all_cells(id_cells_uninfected) ;
    prctile_cell_shape_uninf_step(i,1) = prctile(cshape_cells_uninfected,25)   ;
    prctile_cell_shape_uninf_step(i,2) = prctile(cshape_cells_uninfected,50)   ;
    prctile_cell_shape_uninf_step(i,3) = prctile(cshape_cells_uninfected,75)   ;

end

hold on
plot(prctile_cell_shape_step(1:end,1))
plot(prctile_cell_shape_step(1:end,2))
plot(prctile_cell_shape_step(1:end,3))
legend('25th perc','50th perc','75th perc');
title ('Cell shape index [-]' );
ylabel('Cell shape index [-]' );
xlabel('Steps'         );
hold off
% nbins = 10;
% histogram(cell_shape_index(:,1),nbins);


% Plot cshape infected vs uninfected
figure()
hold on
name1 = "25th perc uninf "+ case_name; name4 = "25th perc inf "+ case_name       ;
name2 = "50th perc uninf "+ case_name; name5 = "50th perc inf "+ case_name       ;
name3 = "75th perc uninf "+ case_name; name6 = "75th perc inf "+ case_name       ;
p1    = plot(prctile_cell_shape_uninf_step(1:end,1), 'LineStyle', ':', 'Color', 'blue');
p2    = plot(prctile_cell_shape_uninf_step(1:end,2),                   'Color', 'blue');
p3    = plot(prctile_cell_shape_uninf_step(1:end,3), 'LineStyle', ':', 'Color', 'blue');
p4    = plot(prctile_cell_shape_inf_step  (1:end,1), 'LineStyle', ':', 'Color', 'red' );
p5    = plot(prctile_cell_shape_inf_step  (1:end,2),                   'Color', 'red' );
p6    = plot(prctile_cell_shape_inf_step  (1:end,3), 'LineStyle', ':', 'Color', 'red' );

lgd = legend([p2, p5], {name2, name5});
title ('Cell shape [-]' );
ylabel('Cell shape [-]' );
xlabel('Steps'         );
hold off
picname    = "cshape_inf_vs_uninf_comparison_" + case_name + ".png";
saveas(gcf, picname);





% Cell shape distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "cshape_" + case_name + ".avi";
    vidObj           = VideoWriter(video_name);
    vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
    open(vidObj);
end

edges = 3.51:0.05:5.1; % histogram edges (range)

for i = 1:total_time_steps

    data_to_analyse = cell_shape_time(:,i)           ;
    id_cells_in     = ~location_time(:,i)            ;
    hist_data       = data_to_analyse(id_cells_in)   ;

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([3.4, 5.2]               );    ylim  ([0, 300]                 );
    xlabel('Cell shape [\mu]'       );    ylabel('frequency of appearence');
    title ('Cell shape distribution');

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


% Cshape infected vs uninfected
video_name       = "cshape_infvsuninf" + case_name + ".avi";
vidObj           = VideoWriter(video_name);
vidObj.FrameRate = 10;    % Set the frame rate to 10 frames per second
open(vidObj);

edges = 3.51:0.05:5.1; % histogram edges (range)

for i = 1:total_time_steps

    hold on

    id_cells_out                 =  logical(location_time(:,i))          ;
    cshape_all_cells             = cell_shape_time(:,i)                  ;
    % Infected versus non-infected
    % Infected:
    id_cells_infected            = logical(infection_time(:,i))          ;
    cshape_cells_infected        = cshape_all_cells(id_cells_infected)   ;
    hist_data                    = cshape_cells_infected                 ;
    histogram(hist_data, 'BinEdges', edges);
    
    % Uninfected: remove cells at the boundary
    id_cells_uninfected          = ~infection_time(:,i)                  ;
    id_cells_uninfected          = id_cells_uninfected-id_cells_out      ; % remove cells at the boundary
    id_cells_uninfected          = logical(id_cells_uninfected)          ; 
    cshape_cells_uninfected      = cshape_all_cells(id_cells_uninfected) ;
    hist_data                    = cshape_cells_uninfected;
    histogram(hist_data, 'BinEdges', edges);

    xlim  ([3.4, 5.2]               );     ylim  ([0, 100]                 );
    xlabel('Cell shape [-]'         );     ylabel('frequency of appearence');
    title ('Cell shape distribution');     legend('infected','uninfected')  ;

    hold off

    % Print shape index limit
    hold on;
    xline = 3.81;
    line([xline, xline], ylim, 'Color', 'r');
    hold off

    % Add frame to video object
    frame       = getframe(gcf);    
    writeVideo(vidObj,frame);
    % Clear figure for next frame
    clf;
end

% Close video object
close(vidObj);







% PLOT DISPLACEMENT VECTOR ................................................

% Displacement vector .........
% Exclude cells on the border for the analysis
prctile_disp_vec_magn_step     = zeros(total_time_steps,3);
prctile_disp_vec_inf_step      = zeros(total_time_steps,3);
prctile_disp_vec_uninf_step    = zeros(total_time_steps,3);

for i = 1 : total_time_steps

    id_cells_in                      = ~location_time(:,i)                       ;
    disp_mag_all_cells               = displacement_time_magn(:,i)               ;
    disp_mag_cells_in                = disp_mag_all_cells(id_cells_in)           ;
    prctile_disp_vec_magn_step(i,1)  = prctile(disp_mag_cells_in,25)             ;
    prctile_disp_vec_magn_step(i,2)  = prctile(disp_mag_cells_in,50)             ;
    prctile_disp_vec_magn_step(i,3)  = prctile(disp_mag_cells_in,75)             ;

    % Infected versus non-infected
    % Infected:
    id_cells_infected                = logical(infection_time(:,i))              ;
    disp_magn_cells_infected         = disp_mag_all_cells(id_cells_infected)     ;
    prctile_disp_vec_inf_step(i,1)   = prctile(disp_magn_cells_infected,25)      ;
    prctile_disp_vec_inf_step(i,2)   = prctile(disp_magn_cells_infected,50)      ;
    prctile_disp_vec_inf_step(i,3)   = prctile(disp_magn_cells_infected,75)      ;
    % Uninfected: remove cells at the boundary
    id_cells_uninfected              = ~infection_time(:,i)                      ;
    id_cells_uninfected              = id_cells_uninfected-id_cells_out          ; % remove cells at the boundary
    id_cells_uninfected              = logical(id_cells_uninfected)              ; 
    disp_magn_cells_uninfected       = disp_mag_all_cells(id_cells_uninfected)   ;
    prctile_disp_vec_uninf_step(i,1) = prctile(disp_magn_cells_uninfected,25)    ;
    prctile_disp_vec_uninf_step(i,2) = prctile(disp_magn_cells_uninfected,50)    ;
    prctile_disp_vec_uninf_step(i,3) = prctile(disp_magn_cells_uninfected,75)    ;

end

hold on
plot(prctile_disp_vec_magn_step(2:end,1))
plot(prctile_disp_vec_magn_step(2:end,2))
plot(prctile_disp_vec_magn_step(2:end,3))
legend('25th perc','50th perc','75th perc');
title ('Cell displacement [\mum]' );
ylabel('Cell displacement [\mum]' );
xlabel('Steps'         );
hold off
% nbins = 10;
% histogram(displacement_time_magn(:,2),nbins);


% Plot displacement magnitude infected vs uninfected
figure()
hold on
name1 = "25th perc uninf "+ case_name; name4 = "25th perc inf "+ case_name       ;
name2 = "50th perc uninf "+ case_name; name5 = "50th perc inf "+ case_name       ;
name3 = "75th perc uninf "+ case_name; name6 = "75th perc inf "+ case_name       ;
p1    = plot(prctile_disp_vec_uninf_step(2:end,1), 'LineStyle', ':', 'Color', 'blue');
p2    = plot(prctile_disp_vec_uninf_step(2:end,2),                   'Color', 'blue');
p3    = plot(prctile_disp_vec_uninf_step(2:end,3), 'LineStyle', ':', 'Color', 'blue');
p4    = plot(prctile_disp_vec_inf_step  (2:end,1), 'LineStyle', ':', 'Color', 'red' );
p5    = plot(prctile_disp_vec_inf_step  (2:end,2),                   'Color', 'red' );
p6    = plot(prctile_disp_vec_inf_step  (2:end,3), 'LineStyle', ':', 'Color', 'red' );

lgd = legend([p2, p5], {name2, name5});
title ('|Cell displacement| [\mum]' );
ylabel('|Cell displacement| [\mum]' );
xlabel('Steps'         );
hold off
picname    = "disp_magn_inf_vs_uninf_comparison_" + case_name + ".png";
saveas(gcf, picname);


% Displacement vector distribution .........
% Set up parameters
if make_video == true
    % Create video object
    video_name       = "disp_magn_" + case_name + ".avi";
    vidObj           = VideoWriter(video_name);
    vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
    open(vidObj);
end

edges                = 0:0.02:2; % histogram edges (range)
% Loop over frames
for i = 1:total_time_steps
    
    data_to_analyse = displacement_time_magn(:,i)   ;
    id_cells_in     = ~location_time(:,i)           ;
    hist_data       = data_to_analyse(id_cells_in)  ;

    histogram(hist_data, 'BinEdges', edges);
    xlim  ([0, 2]                              );   ylim  ([0, 80]                  ); 
    xlabel('Cell displacement magnitude [\mum]');   ylabel('frequency of appearence');
    title ('Cell displacement magnitude'       );

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





% Boxplots

% Area
data_vector          = reshape(area_time, [], 1);
mean_area_cell       = mean(area_time,2); % average each cell over time

boxplot(area_time);
boxplot(data_vector);
boxplot(mean_area_cell);

% Polygon size
data_vector          = reshape(polygon_size_time, [], 1);
median_pol_size_cell = median(polygon_size_time,2); % median of each cell over time

boxplot(polygon_size_time);
boxplot(data_vector);
boxplot(median_pol_size_cell);

% Perimeter
data_vector          = reshape(perimeter_time, [], 1);
mean_perim_cell      = mean(perimeter_time,2); % average each cell over time

boxplot(perimeter_time);
boxplot(data_vector);
boxplot(mean_perim_cell);

% Cell shape
data_vector          = reshape(cell_shape_time, [], 1);
mean_cshape_cell     = mean(cell_shape_time,2); % average each cell over time

boxplot(cell_shape_time);
boxplot(data_vector);
boxplot(mean_cshape_cell);

% Cell displacement vector
data_vector          = reshape(displacement_time_magn, [], 1);
mean_cell_displ_cell = mean(displacement_time_magn,2); % average each cell over time

boxplot(displacement_time_magn(2:end,:));      % skip first step
boxplot(data_vector(number_of_cells+1:end));   % skip cells from the first step
boxplot(mean_cell_displ_cell(2:end,:));


% Save postrpocess parameters ............................................
filename = "var_"+ case_name;
save(filename, 'area_time', 'polygon_size_time', 'perimeter_time', 'cell_shape_time'                 ...
             , 'location_time', 'infection_time', 'prctile_area_step','median_pol_size_step'         ...
             , 'prctile_perim_step', 'prctile_cell_shape_step', 'total_time_steps'                   ...
             , 'number_of_cells', 'displacement_time_magn','prctile_disp_vec_magn_step'              ...
             , 'prctile_area_inf_step','prctile_area_uninf_step','median_pol_size_inf_step'          ...
             , 'median_pol_size_uninf_step', 'prctile_perim_uninf_step','prctile_perim_inf_step'     ...
             , 'prctile_cell_shape_uninf_step', 'prctile_cell_shape_inf_step'                        ...
             , 'prctile_disp_vec_uninf_step', 'prctile_disp_vec_inf_step');
