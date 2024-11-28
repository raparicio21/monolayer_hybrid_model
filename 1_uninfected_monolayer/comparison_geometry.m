%% MONOLAYER POSTPROCESS: GEOMETRY COMPARISON

clear all;
clc;
close all;
name_case_1           = 'case_a2';
name_case_2           = 'case_b2';
name_case_3           = 'case_c2';
composed_name         = 'abc2';
%-------------------------------
string_name           = "var_"+name_case_1+".mat"  ;
load (string_name)
numb_steps            = total_time_steps           ;
numb_cells            = number_of_cells            ;

loc_time_1            = location_time              ;
area_1                = area_time                  ;
pol_size_1            = polygon_size_time          ;
perim_1               = perimeter_time             ;
cell_shape_1          = cell_shape_time            ;
disp_magn_1           = displacement_time_magn     ;

prctile_area_1        = prctile_area_step          ;
med_pol_size_1        = median_pol_size_step       ;
prctile_per_1         = prctile_perim_step         ;
prctile_cshape_1      = prctile_cell_shape_step    ;
prctile_cdispl_magn_1 = prctile_disp_vec_magn_step ;
%-------------------------------
string_name           = "var_"+name_case_2+".mat"  ;
load (string_name)

loc_time_2            = location_time              ;
area_2                = area_time                  ;
pol_size_2            = polygon_size_time          ;
perim_2               = perimeter_time             ;
cell_shape_2          = cell_shape_time            ;
disp_magn_2           = displacement_time_magn     ;

prctile_area_2        = prctile_area_step          ;
med_pol_size_2        = median_pol_size_step       ;
prctile_per_2         = prctile_perim_step         ;
prctile_cshape_2      = prctile_cell_shape_step    ;
prctile_cdispl_magn_2 = prctile_disp_vec_magn_step ;
%-------------------------------
string_name           = "var_"+name_case_3+".mat"  ;
load (string_name)

loc_time_3            = location_time              ;
area_3                = area_time                  ;
pol_size_3            = polygon_size_time          ;
perim_3               = perimeter_time             ;
cell_shape_3          = cell_shape_time            ;
disp_magn_3           = displacement_time_magn     ;

prctile_area_3        = prctile_area_step          ;
med_pol_size_3        = median_pol_size_step       ;
prctile_per_3         = prctile_perim_step         ;
prctile_cshape_3      = prctile_cell_shape_step    ;
prctile_cdispl_magn_3 = prctile_disp_vec_magn_step ;


clear area_time polygon_size_time perimeter_time cell_shape_time prctile_area_step ...
    median_pol_size_step prctile_perim_step prctile_cell_shape_step string_name    ...
    displacement_time_magn prctile_disp_vec_magn_step location_time;


% Area comparison
hold on
plot(prctile_area_1(1:end,1), 'LineStyle', ':', 'Color', 'blue')
plot(prctile_area_1(1:end,2),                   'Color', 'blue')
plot(prctile_area_1(1:end,3), 'LineStyle', ':', 'Color', 'blue')
plot(prctile_area_3(1:end,1), 'LineStyle', ':', 'Color', 'red' )
plot(prctile_area_3(1:end,2),                   'Color', 'red' )
plot(prctile_area_3(1:end,3), 'LineStyle', ':', 'Color', 'red' )
name1 = "25th perc "+ name_case_1; name4 = "25th perc "+ name_case_3;
name2 = "50th perc "+ name_case_1; name5 = "50th perc "+ name_case_3;
name3 = "75th perc "+ name_case_1; name6 = "75th perc "+ name_case_3;
legend(name1,name2,name3,name4,name5,name6);
title ('Area [\mum^2]' );
ylabel('Area [\mum^2]' );
xlabel('Steps'         );
hold off
picname    = "prctile_area_comparison_" + composed_name + ".png";
saveas(gcf, picname);


% Area distribution .........
% Create video object
video_name       = "area_comparison_" + composed_name + ".avi";
vidObj           = VideoWriter(video_name);
vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
open(vidObj);
edges = 120:20:800; % histogram edges (range)

for i = 1:total_time_steps

    hold on

    data_to_analyse    = area_1(:,i)                   ;
    id_cells_in        = ~loc_time_1(:,i)              ;
    hist_data          = data_to_analyse(id_cells_in)  ;
    histogram(hist_data, 'BinEdges', edges);

%     data_to_analyse    = area_2(:,i)                   ;
%     id_cells_in        = ~loc_time_2(:,i)              ;
%     hist_data          = data_to_analyse(id_cells_in)  ;
%     histogram(hist_data, 'BinEdges', edges);

    data_to_analyse    = area_3(:,i)                   ;
    id_cells_in        = ~loc_time_3(:,i)              ;
    hist_data          = data_to_analyse(id_cells_in)  ;
    histogram(hist_data, 'BinEdges', edges);

    xlim  ([0, 800]          ) ;     ylim  ([0, 100]                 );
    xlabel('area [\mum^2]'   ) ;     ylabel('frequency of appearence');
    title ('area distribution');     legend(name_case_1,name_case_3)  ;

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
% the median does not give you relevant information in this case


% Polygon size distribution .........
video_name       = "pol_size_comparison" + composed_name + ".avi";
vidObj           = VideoWriter(video_name);
vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
open(vidObj);


for i = 1:total_time_steps

    hold on

    data_to_analyse = pol_size_1(:,i)               ;
    id_cells_in     = ~loc_time_1(:,i)              ;
    hist_data       = data_to_analyse(id_cells_in)  ;
    histogram(hist_data,'BinMethod','integers')     ;

%     data_to_analyse = pol_size_2(:,i)               ;
%     id_cells_in     = ~loc_time_2(:,i)              ;
%     hist_data       = data_to_analyse(id_cells_in)  ;
%     histogram(hist_data,'BinMethod','integers')     ;

    data_to_analyse = pol_size_3(:,i)               ;
    id_cells_in     = ~loc_time_3(:,i)              ;
    hist_data       = data_to_analyse(id_cells_in)  ;
    histogram(hist_data,'BinMethod','integers')     ;

    xlim  ([2, 10]                    );    ylim  ([0, 600]                 );
    xlabel('Polygon size'             );    ylabel('frequency of appearence');
    title ('Polygon size distribution');    legend(name_case_1,name_case_3)  ;

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
hold on
plot(prctile_per_1(1:end,1), 'LineStyle', ':', 'Color', 'blue')
plot(prctile_per_1(1:end,2),                   'Color', 'blue')
plot(prctile_per_1(1:end,3), 'LineStyle', ':', 'Color', 'blue')
plot(prctile_per_3(1:end,1), 'LineStyle', ':', 'Color', 'red' )
plot(prctile_per_3(1:end,2),                   'Color', 'red' )
plot(prctile_per_3(1:end,3), 'LineStyle', ':', 'Color', 'red' )
name1 = "25th perc "+ name_case_1; name4 = "25th perc "+ name_case_3;
name2 = "50th perc "+ name_case_1; name5 = "50th perc "+ name_case_3;
name3 = "75th perc "+ name_case_1; name6 = "75th perc "+ name_case_3;
legend(name1,name2,name3,name4,name5,name6);
title ('Perimeter [\mum]' );
ylabel('Perimeter [\mum]' );
xlabel('Steps'            );
hold off
picname    = "prctile_perim_comparison_" + composed_name + ".png";
saveas(gcf, picname);




% Perimeter distribution .........
video_name       = "perimeter_comparison" + composed_name + ".avi";
vidObj           = VideoWriter(video_name);
vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
open(vidObj);

edges            = 50:2:140;  % histogram edges (range)

for i = 1:total_time_steps

    hold on

    data_to_analyse    = perim_1(:,i)                   ;
    id_cells_in        = ~loc_time_1(:,i)              ;
    hist_data          = data_to_analyse(id_cells_in)  ;
    histogram(hist_data, 'BinEdges', edges);

%     data_to_analyse    = perim_2(:,i)                   ;
%     id_cells_in        = ~loc_time_2(:,i)              ;
%     hist_data          = data_to_analyse(id_cells_in)  ;
%     histogram(hist_data, 'BinEdges', edges);

    data_to_analyse    = perim_3(:,i)                   ;
    id_cells_in        = ~loc_time_3(:,i)              ;
    hist_data          = data_to_analyse(id_cells_in)  ;
    histogram(hist_data, 'BinEdges', edges);

    xlim  ([40, 120]               );   ylim  ([0, 300]                 );
    xlabel('perimeter [\mum]'      );   ylabel('frequency of appearence');
    title ('Perimeter distribution');   legend(name_case_1,name_case_3)  ;

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
hold on
plot(prctile_cshape_1(1:end,1), 'LineStyle', ':', 'Color', 'blue')
plot(prctile_cshape_1(1:end,2),                   'Color', 'blue')
plot(prctile_cshape_1(1:end,3), 'LineStyle', ':', 'Color', 'blue')
plot(prctile_cshape_3(1:end,1), 'LineStyle', ':', 'Color', 'red' )
plot(prctile_cshape_3(1:end,2),                   'Color', 'red' )
plot(prctile_cshape_3(1:end,3), 'LineStyle', ':', 'Color', 'red' )
name1 = "25th perc "+ name_case_1; name4 = "25th perc "+ name_case_3;
name2 = "50th perc "+ name_case_1; name5 = "50th perc "+ name_case_3;
name3 = "75th perc "+ name_case_1; name6 = "75th perc "+ name_case_3;
legend(name1,name2,name3,name4,name5,name6);
title ('Cell shape index [-]' );
ylabel('Cell shape index [-]' );
xlabel('Steps'                );
hold off
picname    = "prctile_cshape_comparison_" + composed_name + ".png";
saveas(gcf, picname);




% Cell shape distribution .........
video_name       = "cell_shape_comparison" + composed_name + ".avi";
vidObj           = VideoWriter(video_name);
vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
open(vidObj);

edges = 3.51:0.05:5.1; % histogram edges (range)

for i = 1:total_time_steps

    hold on

    data_to_analyse    = cell_shape_1(:,i)             ;
    id_cells_in        = ~loc_time_1(:,i)              ;
    hist_data          = data_to_analyse(id_cells_in)  ;
    histogram(hist_data, 'BinEdges', edges)            ;

%     data_to_analyse    = cell_shape_2(:,i)             ;
%     id_cells_in        = ~loc_time_2(:,i)              ;
%     hist_data          = data_to_analyse(id_cells_in)  ;
%     histogram(hist_data, 'BinEdges', edges)            ;

    data_to_analyse    = cell_shape_3(:,i)             ;
    id_cells_in        = ~loc_time_3(:,i)              ;
    hist_data          = data_to_analyse(id_cells_in)  ;
    histogram(hist_data, 'BinEdges', edges)            ;

    xlim  ([3.4, 5.2]               );    ylim  ([0, 300]                 );
    xlabel('Cell shape [-]'         );    ylabel('frequency of appearence');
    title ('Cell shape distribution');    legend(name_case_1,name_case_3)  ;

    % Print shape index limit
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









% PLOT CELL DISPLACEMENTS .................................................

% Cell displacement .........
hold on
plot(prctile_cdispl_magn_1(2:end,1), 'LineStyle', ':', 'Color', 'blue')
plot(prctile_cdispl_magn_1(2:end,2),                   'Color', 'blue')
plot(prctile_cdispl_magn_1(2:end,3), 'LineStyle', ':', 'Color', 'blue')
plot(prctile_cdispl_magn_3(2:end,1), 'LineStyle', ':', 'Color', 'red' )
plot(prctile_cdispl_magn_3(2:end,2),                   'Color', 'red' )
plot(prctile_cdispl_magn_3(2:end,3), 'LineStyle', ':', 'Color', 'red' )
name1 = "25th perc "+ name_case_1; name4 = "25th perc "+ name_case_3;
name2 = "50th perc "+ name_case_1; name5 = "50th perc "+ name_case_3;
name3 = "75th perc "+ name_case_1; name6 = "75th perc "+ name_case_3;
legend(name1,name2,name3,name4,name5,name6);
title ('Cell displacement [\mum]' );
ylabel('Cell displacement [\mum]' );
xlabel('Steps'                    );
hold off
picname    = "prctile_cdispl_magn_comparison_" + composed_name + ".png";
saveas(gcf, picname);




% Cell displacement distribution .........
video_name       = "cell_disp_magn_comparison" + composed_name + ".avi";
vidObj           = VideoWriter(video_name);
vidObj.FrameRate = 10; % Set the frame rate to 10 frames per second
open(vidObj);

edges            = 0:0.02:2; % histogram edges (range)

for i = 2:total_time_steps

    hold on

    data_to_analyse    = disp_magn_1(:,i)              ;
    id_cells_in        = ~loc_time_1(:,i)              ;
    hist_data          = data_to_analyse(id_cells_in)  ;
    histogram(hist_data, 'BinEdges', edges)            ;

%     data_to_analyse    = disp_magn_2(:,i)              ;
%     id_cells_in        = ~loc_time_2(:,i)              ;
%     hist_data          = data_to_analyse(id_cells_in)  ;
%     histogram(hist_data, 'BinEdges', edges)            ;

    data_to_analyse    = disp_magn_3(:,i)              ;
    id_cells_in        = ~loc_time_3(:,i)              ;
    hist_data          = data_to_analyse(id_cells_in)  ;
    histogram(hist_data, 'BinEdges', edges)            ;

    xlim  ([0, 2]                              );   ylim  ([0, 80]                  ); 
    xlabel('Cell displacement magnitude [\mum]');   ylabel('frequency of appearence');
    title ('Cell displacement magnitude'       );   legend(name_case_1,name_case_3)  ;


    % Add frame to video object
    frame       = getframe(gcf);
    writeVideo(vidObj,frame);
    % Clear figure for next frame
    clf;
    
end

% Close video object
close(vidObj);

%
% Box plot

% Area
% data_vector          = reshape(area_1, [], 1);
median_area_cell_1      = median(area_1,2); % median for each cell over time
median_area_cell_2      = median(area_2,2); % median for each cell over time
median_area_cell_3      = median(area_3,2); % median for each cell over time

median_area_cell        = [ median_area_cell_1; median_area_cell_2; median_area_cell_3];

g1 = repmat({name_case_1},numb_cells,1);
g2 = repmat({name_case_2},numb_cells,1);
g3 = repmat({name_case_3},numb_cells,1);
g  = [g1; g2; g3];

boxplot(median_area_cell,g)
ylabel('Area (median for each cell) [\mum^2]');
picname    = "boxplot_area_foreachcell_" + composed_name + ".png";
saveas(gcf, picname);



% Polygon size
% data_vector          = reshape(pol_size_1, [], 1);
median_pol_siz_cell_1      = median(pol_size_1,2); % median for each cell over time
median_pol_siz_cell_2      = median(pol_size_2,2); % median for each cell over time
median_pol_siz_cell_3      = median(pol_size_3,2); % median for each cell over time

median_pol_siz_cell        = [ median_pol_siz_cell_1; median_pol_siz_cell_2; median_pol_siz_cell_3];

g1 = repmat({name_case_1},numb_cells,1);
g2 = repmat({name_case_2},numb_cells,1);
g3 = repmat({name_case_3},numb_cells,1);
g  = [g1; g2; g3];

boxplot(median_pol_siz_cell,g)
ylabel('Polygon size (median for each cell)');
picname    = "boxplot_pol_size_foreachcell_" + composed_name + ".png";
saveas(gcf, picname);



% Perimeter
% data_vector          = reshape(area_1, [], 1);
median_perim_cell_1      = median(perim_1,2); % median for each cell over time
median_perim_cell_2      = median(perim_2,2); % median for each cell over time
median_perim_cell_3      = median(perim_3,2); % median for each cell over time

median_perim_cell        = [ median_perim_cell_1; median_perim_cell_2; median_perim_cell_3];

g1 = repmat({name_case_1},numb_cells,1);
g2 = repmat({name_case_2},numb_cells,1);
g3 = repmat({name_case_3},numb_cells,1);
g = [g1; g2; g3];

boxplot(median_perim_cell,g)
ylabel('Perimeter (median for each cell) [\mum]');
picname    = "boxplot_perim_foreachcell_" + composed_name + ".png";
saveas(gcf, picname);





% Cell shape
% data_vector          = reshape(area_1, [], 1);
median_cshape_cell_1     = median(cell_shape_1,2); % median for each cell over time
median_cshape_cell_2     = median(cell_shape_2,2); % median for each cell over time
median_cshape_cell_3     = median(cell_shape_3,2); % median for each cell over time

median_cshape_cell       = [ median_cshape_cell_1; median_cshape_cell_2; median_cshape_cell_3];

g1 = repmat({name_case_1},numb_cells,1);
g2 = repmat({name_case_2},numb_cells,1);
g3 = repmat({name_case_3},numb_cells,1);
g = [g1; g2; g3];

boxplot(median_cshape_cell,g)
ylabel('Cell shape (median for each cell) [-]');
picname    = "boxplot_cshape_foreachcell_" + composed_name + ".png";
saveas(gcf, picname);




% Cell displacement magnitude
% data_vector           = reshape(area_1, [], 1);
median_cdisp_magn_cell_1  = median(disp_magn_1,2); % median for each cell over time
median_cdisp_magn_cell_2  = median(disp_magn_2,2); % median for each cell over time
median_cdisp_magn_cell_3  = median(disp_magn_3,2); % median for each cell over time

median_cdisp_magn_cell    = [ median_cdisp_magn_cell_1; median_cdisp_magn_cell_2; median_cdisp_magn_cell_3];

g1 = repmat({name_case_1},numb_cells,1);
g2 = repmat({name_case_2},numb_cells,1);
g3 = repmat({name_case_3},numb_cells,1);
g  = [g1; g2; g3];

boxplot(median_cdisp_magn_cell,g)
ylabel('|Cell displacement| (median for each cell) [\mum]');
picname    = "boxplot_cdisp_foreachcell_" + composed_name + ".png";
saveas(gcf, picname);
