classdef vtk_mat_Parser

    methods

        %................... Object initialization ........................
        function obj = vtk_mat_Parser()
            addpath(fullfile('inputs')); 
            addpath(fullfile('src'));
        end




        %........... Write monolayer data to VTK file format ..............
        %..................................................................
        function obj = writeVTK_mat_file(obj,monolayer,stress_results)

            obj.writeData(monolayer,stress_results);

            
        end






        function writeData(~,monolayer,stress_results)

            % Load coordinates
            x = monolayer.monolayer_cell_points(:,1);
            y = monolayer.monolayer_cell_points(:,2);
            z = zeros(size(monolayer.monolayer_cell_points(:,1)));

            % Choose variables to print
%             % Cell geometry
%             variables.id                        = inputs.output_files_parameters_id;
%             variables.radius                    = inputs.output_files_parameters_radius;
%             variables.area                      = inputs.output_files_parameters_area;
%             variables.polygon_size              = inputs.output_files_parameters_polygon_size;
%             variables.perimeter                 = inputs.output_files_parameters_perimeter;
%             variables.location                  = inputs.output_files_parameters_location;
%     
%             % Cell kinematics
%             variables.initial_position          = inputs.output_files_parameters_initial_position;
%             variables.current_position          = inputs.output_files_parameters_current_position;
%             variables.displacement_vector       = inputs.output_files_parameters_displacement_vector;
%             variables.velocity_vector           = inputs.output_files_parameters_velocity_vector;
%     
%             % Cell total kinematics
%             variables.total_distance            = inputs.output_files_parameters_total_distance;
%             variables.total_displacement        = inputs.output_files_parameters_total_displacement; 
%             variables.total_displacement_vector = inputs.output_files_parameters_total_displacement_vector;
%     
%             % Cell forces       
%             variables.forces                    = inputs.output_files_parameters_forces;
%             variables.cell_cell_int_forces      = inputs.output_files_parameters_cell_cell_int_forces;
%             variables.cell_cont_cent_forces     = inputs.output_files_parameters_cell_cont_cent_forces;

            % Initialize output variables
            % Cell geometry
            %id_out        = zeros(numel(monolayer.cell_id),1);
            %radius_out    = zeros(numel(monolayer.cell_id),1);
            area_out      = zeros(numel(monolayer.cell_id),1);
            pol_sizes_out = zeros(numel(monolayer.cell_id),1);
            perimeter_out = zeros(numel(monolayer.cell_id),1);
            %cellshape_out = zeros(numel(monolayer.cell_id),1);

            % Cell state
            %location_out  = zeros(numel(monolayer.cell_id),1);
            %infected_out  = zeros(numel(monolayer.cell_id),1);
            %surround_out  = zeros(numel(monolayer.cell_id),1);

            % Cell kinematics
            %init_pos_out  = zeros(numel(monolayer.cell_id),3);
            %curr_pos_out  = zeros(numel(monolayer.cell_id),3);
            %disp_vec_out  = zeros(numel(monolayer.cell_id),3);
            %velo_vec_out  = zeros(numel(monolayer.cell_id),3);

            % Cell total kinematics
            %tot_dist_out  = zeros(numel(monolayer.cell_id),1); 
            %tot_disp_out  = zeros(numel(monolayer.cell_id),1); 
            %t_dis_ve_out  = zeros(numel(monolayer.cell_id),3);

            % Cell forces
            %tot_force_out = zeros(numel(monolayer.cell_id),3);
            %int_force_out = zeros(numel(monolayer.cell_id),3);
            %con_force_out = zeros(numel(monolayer.cell_id),3); % centroid
            %rand_walk_out = zeros(numel(monolayer.cell_id),3); 
            
            
            % Update output variables .....................................

            % Cell geometry
            id_out           =        monolayer.cell_id                                              ;
            radius_out       =        monolayer.cell_radius                                          ;

            % Cell state
            location_out     =        monolayer.cellLocationList                                     ;
            infected_out     =        monolayer.cell_infection                                       ;
            surround_out     =        monolayer.surrounder_cell                                      ;
            surround_2nd_out =        monolayer.surrounder_cell_2nd                                  ;
            surround_3rd_out =        monolayer.surrounder_cell_3rd                                  ;

            % Cell kinematics
            init_pos_out     =        monolayer.cell_init_position                                   ;
            curr_pos_out     =        monolayer.cell_position                                        ;

            boundary_condition = inputs.getBoundaryCondition(monolayer.step)                         ;

            if strcmp(boundary_condition,'PBC') || strcmp(boundary_condition,'material_around')

                disp_vec_out     =        monolayer.cell_PBC_position_PBC - monolayer.cell_previous_position;
            else
                disp_vec_out     =        curr_pos_out - monolayer.cell_previous_position            ;
            end
            
            velo_vec_out     =        disp_vec_out ./ inputs.general_time_step_duration              ;

            % Cell total kinematics
            tot_dist_out     =        monolayer.cell_total_distance  ; % +  sqrt(sum(disp_vec_out.^2, 2)); It is already computed when you move the cells at the end of the loop
            tot_disp_out     =        sqrt(   sum(  (curr_pos_out - init_pos_out).^2   , 2)         );
            t_dis_ve_out     =        curr_pos_out - init_pos_out                                    ;

            % Cell forces
            tot_force_out    =        monolayer.force_list                                           ; 
            int_force_out    =        monolayer.force_cell_cell_int_list                             ;
            con_force_out    =        monolayer.force_cell_cont_cent_list                            ;  
            rand_walk_out    =        monolayer.force_random_walk                                    ; 
            inf_crawl_out    =        monolayer.force_inf_crawl_cent_list                            ;
            pro_force_out    =        monolayer.force_cell_prot_cent_list                            ; 

            extruding_cell   =        monolayer.cell_extrusion                                       ;
            for i = 1 : numel(monolayer.cell_id)

                if extruding_cell(i) == 0 % Non-extruding cell
                
                    % Cell geometry
                    area_out      (i)   =        monolayer.cellPolyList{i}.area                        ;
                    pol_sizes_out (i)   = length(monolayer.cellPolyList{i}.Vertices(:,1))              ;
                    perimeter_out (i)   =        monolayer.cellPolyList{i}.perimeter                   ;
                end

            end  

            % Cell geometry
            cellshape_out     = perimeter_out ./ sqrt(area_out)   ;

            % Cell average weighted stresses from previous step
            step              = monolayer.step                             ;
            if step == 1
                % weighted_sigma_I                         = zeros(numel(monolayer.cell_id),1)      ;
                % weighted_sigma_II                        = zeros(numel(monolayer.cell_id),1)      ;
                weighted_force_I                         = zeros(numel(monolayer.cell_id),1)      ;
                weighted_force_II                        = zeros(numel(monolayer.cell_id),1)      ;

                orientation_min_s_I                      = zeros(numel(monolayer.cell_id),1)      ;
                orientation_min_s_II                     = zeros(numel(monolayer.cell_id),1)      ;

                %accumulated_strain                       = zeros(numel(monolayer.cell_id),1)      ;
            else
                % Variables saved on fem_results.average_cell_stress
                % Column 1 : Cell averaged weighted stress sigma_I  (all range         ) 
                % Column 2 : Cell averaged weighted stress sigma_II (all range         )
                % Column 3 : Cell averaged weighted stress sigma_I  (range [0  s_I_max]) positive values
                % Column 4 : Cell averaged weighted stress sigma_II (range [s_II_max 0]) negative values 
                % Column 5 : Column 1 * cell area = averaged cell force_I  (all range    )
                % Column 6 : Column 2 * cell area = averaged cell force_II (all range    )
                % Column 7 : Column 3 * cell area = averaged cell force_I  (limited range)
                % Column 8 : Column 4 * cell area = averaged cell force_II (limited range)
                % Column 9 : averaged cell force_I  from the resultant stress tensor (principal stress I )
                % Column 10: averaged cell force_II from the resultant stress tensor (principal stress II)


                % Finally not computed
                % Column 11: weighted e_I
                % Column 12: weighted e_II
                % Column 13: weighted e_III
                % Column 14: weighted e invariant 1 (volumetric strain)


                % weighted_sigma_I         = stress_results.average_cell_stress(:,3)      ;
                % weighted_sigma_II        = stress_results.average_cell_stress(:,4)      ;
                weighted_force_I         = stress_results.average_cell_stress(:,9)       ;
                weighted_force_II        = stress_results.average_cell_stress(:,10)      ;

                orientation_min_s_I      = monolayer.orientation_min_sigma_I_prev_i      ;
                orientation_min_s_II     = monolayer.orientation_min_sigma_II_prev_i     ;

                %accumulated_strain       = monolayer.accumulated_strain_invariant_1*1000 ; 




            end

            if step <= 19
                area_change = zeros(numel(monolayer.cell_id),1)   ;
            elseif step > 19
                area_change = monolayer.area_change               ;
            end

            % Create vectors for the minimal stress direction in order to
            % represent it
            orientation_min_s_I_vect       = zeros(numel(monolayer.cell_id),3)    ;
            orientation_min_s_II_vect      = zeros(numel(monolayer.cell_id),3)    ;

            magnitude                      = 15                                   ;
            x_component                    = magnitude * cos(orientation_min_s_I) ;
            y_component                    = magnitude * sin(orientation_min_s_I) ;
            orientation_min_s_I_vect(:,1)  = x_component                          ;
            orientation_min_s_I_vect(:,2)  = y_component                          ;

            magnitude                      = 15                                   ;
            x_component                    = magnitude * cos(orientation_min_s_II);
            y_component                    = magnitude * sin(orientation_min_s_II);
            orientation_min_s_II_vect(:,1) = x_component                          ;
            orientation_min_s_II_vect(:,2) = y_component                          ;



            %.................... Write mat files .........................
            % In case there was extrusion, we save all variables for the
            % analysis. We will take care of it later on when we do the
            % postprocess
            if inputs.output_files_mat_file==true
                filename = sprintf('./results/raw/data_%05i.mat', step);
                save(filename, 'id_out'             , 'radius_out'       , 'area_out'          , 'pol_sizes_out'           , 'perimeter_out'   ,           ...
                               'cellshape_out'      ,                                                                                                      ...
                               'location_out'       , 'infected_out'     , 'surround_out'      , 'surround_2nd_out'        , 'surround_3rd_out',           ...
                               'init_pos_out'       , 'curr_pos_out'     , 'disp_vec_out'      , 'velo_vec_out'            ,                               ...
                               'tot_dist_out'       , 'tot_disp_out'     , 't_dis_ve_out'      ,                                                           ...
                               'tot_force_out'      , 'int_force_out'    , 'con_force_out'     , 'rand_walk_out'           , 'inf_crawl_out'             , ...
                               'pro_force_out'      , 'weighted_force_I' , 'weighted_force_II' , 'orientation_min_s_I_vect', 'orientation_min_s_II_vect' , ...
                               'area_change'                                                                                                             );
            end
           






            %.................... Write VTK files .........................
            if inputs.output_files_vtk_file==true

                filename = sprintf('./results/cells/cell_%05i.vtk', step); 

                % Avoid the representation of extruding cells
                id_extruding_cells     = monolayer.index_extruding_cells;
                if any(id_extruding_cells)

                    x(id_extruding_cells) = []; y(id_extruding_cells) = []; z(id_extruding_cells) = [];
    
                                      id_out(id_extruding_cells  ) = [];                       radius_out(id_extruding_cells  ) = []; 
                                    area_out(id_extruding_cells  ) = [];                    pol_sizes_out(id_extruding_cells  ) = []; 
                               perimeter_out(id_extruding_cells  ) = [];                     location_out(id_extruding_cells  ) = []; 
                                init_pos_out(id_extruding_cells,:) = [];                     curr_pos_out(id_extruding_cells,:) = []; 
                                disp_vec_out(id_extruding_cells,:) = [];                     velo_vec_out(id_extruding_cells,:) = []; 
                                tot_dist_out(id_extruding_cells  ) = [];                     tot_disp_out(id_extruding_cells  ) = []; 
                                t_dis_ve_out(id_extruding_cells,:) = [];                    tot_force_out(id_extruding_cells,:) = []; 
                               int_force_out(id_extruding_cells,:) = [];                    con_force_out(id_extruding_cells,:) = []; 
                               rand_walk_out(id_extruding_cells,:) = [];                     infected_out(id_extruding_cells  ) = []; 
                                surround_out(id_extruding_cells  ) = [];                    inf_crawl_out(id_extruding_cells,:) = []; 
                               cellshape_out(id_extruding_cells  ) = [];                    pro_force_out(id_extruding_cells,:) = []; 
                            surround_2nd_out(id_extruding_cells  ) = [];                 surround_3rd_out(id_extruding_cells  ) = []; 
                            weighted_force_I(id_extruding_cells  ) = [];                weighted_force_II(id_extruding_cells  ) = []; 
                    orientation_min_s_I_vect(id_extruding_cells,:) = [];        orientation_min_s_II_vect(id_extruding_cells,:) = []; 
                                 area_change(id_extruding_cells,:) = []; 
    
                end
    
                % Write cells agent-based model
                vtkwrite_cells(   filename                                   , 'POLYDATA'                      , x, y, z                                                           , ...
                                 'int'     , 'id'                            , id_out                   ,    'scalars' , 'radius'                      , radius_out                , ...
                                 'scalars' , 'area'                          , area_out                 ,    'int'     , 'polygon_size'                , pol_sizes_out             , ... 
                                 'scalars' , 'perimeter'                     , perimeter_out            ,    'int'     , 'location'                    , location_out              , ...
                                 'vectors' , 'initial_position'              , init_pos_out             ,    'vectors' , 'current_position'            , curr_pos_out              , ...
                                 'vectors' , 'displacement_vector'           , disp_vec_out             ,    'vectors' , 'velocity_vector'             , velo_vec_out              , ...
                                 'scalars' , 'total_distance'                , tot_dist_out             ,    'scalars' , 'total_displacement'          , tot_disp_out              , ...
                                 'vectors' , 'total_displacement_vector'     , t_dis_ve_out             ,    'vectors' , 'forces'                      , tot_force_out             , ...
                                 'vectors' , 'cell_cell_interaction_forces'  , int_force_out            ,    'vectors' , 'cell_cont_cent_forces'       , con_force_out             , ...
                                 'vectors' , 'rand_walk_forces'              , rand_walk_out            ,    'int'     , 'infected'                    , infected_out              , ...
                                 'int'     , '1st_line'                      , surround_out             ,    'vectors' , 'cell_inf_crawl_forces'       , inf_crawl_out             , ...  
                                 'scalars' , 'cell_shape'                    , cellshape_out            ,    'vectors' , 'cell_prot_cent_forces'       , pro_force_out             , ...
                                 'int'     , '2nd_line'                      , surround_2nd_out         ,    'int'     , '3rd_line'                    , surround_3rd_out          , ...   
                                 'scalars' , 'weighted_force_I_prev_step'    , weighted_force_I         ,    'scalars' , 'weighted_force_II_prev_step' , weighted_force_II         , ... 
                                 'vectors' , 'minim_s_I_dir_prev_step'       , orientation_min_s_I_vect ,    'vectors' , 'minim_s_II_dir_prev_step'    , orientation_min_s_II_vect , ...
                                 'scalars' , 'area_change'                   , area_change              ,                                                                            ...
                                 'Precision', 2); 
    
                % Write Voronoi agent-based model
                filename = sprintf('./results/voronoi/voronoi_%05i.vtk', step);
                v_x      = monolayer.voronoi_vertices(:,1);
                v_y      = monolayer.voronoi_vertices(:,2);
                v_z      = zeros(length(v_x),1);
                vtkwrite_polygons( filename , 'POLYGONS' , v_x, v_y, v_z, pol_sizes_out, monolayer.voronoi_region                                                                   , ...
                                 'int'     , 'id'                            , id_out                    ,    'scalars' , 'radius'                      , radius_out                , ...
                                 'scalars' , 'area'                          , area_out                  ,    'int'     , 'polygon_size'                , pol_sizes_out             , ... 
                                 'scalars' , 'perimeter'                     , perimeter_out             ,    'int'     , 'location'                    , location_out              , ...
                                 'vectors' , 'initial_position'              , init_pos_out              ,    'vectors' , 'current_position'            , curr_pos_out              , ...
                                 'vectors' , 'displacement_vector'           , disp_vec_out              ,    'vectors' , 'velocity_vector'             , velo_vec_out              , ...
                                 'scalars' , 'total_distance'                , tot_dist_out              ,    'scalars' , 'total_displacement'          , tot_disp_out              , ...
                                 'vectors' , 'total_displacement_vector'     , t_dis_ve_out              ,    'vectors' , 'forces'                      , tot_force_out             , ...
                                 'vectors' , 'cell_cell_interaction_forces'  , int_force_out             ,    'vectors' , 'cell_cont_cent_forces'       , con_force_out             , ...
                                 'vectors' , 'rand_walk_forces'              , rand_walk_out             ,    'int'     , 'infected'                    , infected_out              , ...
                                 'int'     , '1st_line'                      , surround_out              ,    'vectors' , 'cell_inf_crawl_forces'       , inf_crawl_out             , ...
                                 'scalars' , 'cell_shape'                    , cellshape_out             ,    'vectors' , 'cell_prot_cent_forces'       , pro_force_out             , ...
                                 'int'     , '2nd_line'                      , surround_2nd_out          ,    'int'     , '3rd_line'                    , surround_3rd_out          , ...  
                                 'scalars' , 'weighted_force_I_prev_step'    , weighted_force_I          ,    'scalars' , 'weighted_force_II_prev_step' , weighted_force_II         , ...  
                                 'vectors' , 'minim_s_I_dir_prev_step'       , orientation_min_s_I_vect  ,    'vectors' , 'minim_s_II_dir_prev_step'    , orientation_min_s_II_vect , ...
                                 'scalars' , 'area_change'                   , area_change               ,                                                                            ...
                                 'Precision', 2);
    
                % Write alpha_shape
%                 filename  = sprintf('./results/alpha_shape/alpha_shape_%05i.vtk', step);
%                 v_alpha_x = monolayer.alpha_shape(:,1);
%                 v_alpha_y = monolayer.alpha_shape(:,2);
%                 v_alpha_z = zeros(length(v_alpha_x),1);
    
                % Modify code when there is more than one alpha shape
                % (infection for example) CAREFUL
%                 vtkwrite_alphashape( filename , 'ALPHASHAPE' , v_alpha_x, v_alpha_y, v_alpha_z);

            end







            %.......... Write elipsoids (major and minor axis) ............

            % Load centroids
            list_centroids_x       = monolayer.centroid_of_elipsoid_step_i(:,1);
            list_centroids_y       = monolayer.centroid_of_elipsoid_step_i(:,2);
            list_centroids_z       = monolayer.centroid_of_elipsoid_step_i(:,3);

            % To show major or minor axis, plot two arrows from the
            % centroid, so we need to duplicate the centroids to associate
            % the first and second arrow
            list_centroids_x       = [list_centroids_x ; list_centroids_x];
            list_centroids_y       = [list_centroids_y ; list_centroids_y];
            list_centroids_z       = [list_centroids_z ; list_centroids_z];

            % Project axis to get x and y components

            % Major axis
            % Define the angle and magnitude of the vector
            angle                  = monolayer.orientation_axis_elipsoid_i;   % Angle in rad
            magnitude              = monolayer.major_axis_length_i/2;         % Major axis/2
                       
            % Compute the x and y coordinates
            major_axis_x_1st_arrow = magnitude .* cos(angle);
            major_axis_y_1st_arrow = magnitude .* sin(angle);
            major_axis_z_1st_arrow = zeros(size(major_axis_x_1st_arrow,1),1);

            major_axis_x_2nd_arrow = - major_axis_x_1st_arrow;
            major_axis_y_2nd_arrow = - major_axis_y_1st_arrow;
            major_axis_z_2nd_arrow = zeros(size(major_axis_x_2nd_arrow,1),1);

            major_axis_x_together  = [major_axis_x_1st_arrow ; major_axis_x_2nd_arrow];
            major_axis_y_together  = [major_axis_y_1st_arrow ; major_axis_y_2nd_arrow];
            major_axis_z_together  = [major_axis_z_1st_arrow ; major_axis_z_2nd_arrow];

            major_axis_alltogether = [major_axis_x_together major_axis_y_together major_axis_z_together];



            % Minor axis
            % Define the angle and magnitude of the vector
            angle                  = monolayer.orientation_axis_elipsoid_i + 90*pi/180;   % Angle in rad
            magnitude              = monolayer.minor_axis_length_i/2;                     % Minor axis/2
                       
            % Compute the x and y coordinates
            minor_axis_x_1st_arrow = magnitude .* cos(angle);
            minor_axis_y_1st_arrow = magnitude .* sin(angle);
            minor_axis_z_1st_arrow = zeros(size(minor_axis_x_1st_arrow,1),1);

            minor_axis_x_2nd_arrow = - minor_axis_x_1st_arrow;
            minor_axis_y_2nd_arrow = - minor_axis_y_1st_arrow;
            minor_axis_z_2nd_arrow = zeros(size(minor_axis_x_2nd_arrow,1),1);

            minor_axis_x_together  = [minor_axis_x_1st_arrow ; minor_axis_x_2nd_arrow];
            minor_axis_y_together  = [minor_axis_y_1st_arrow ; minor_axis_y_2nd_arrow];
            minor_axis_z_together  = [minor_axis_z_1st_arrow ; minor_axis_z_2nd_arrow];

            minor_axis_alltogether = [minor_axis_x_together minor_axis_y_together minor_axis_z_together];


            filename               = sprintf('./results/elipsoids/elipsoids_%05i.vtk', step); 
    
            % Write elipsoid major and minor axis (info from agent based model)
            vtkwrite_cells(   filename                                   , 'POLYDATA'                , list_centroids_x, list_centroids_y, list_centroids_z , ...
                             'vectors' , 'major_axis'                    , major_axis_alltogether    , ...
                             'vectors' , 'minor_axis'                    , minor_axis_alltogether    , ...
                             'Precision', 2);





            if inputs.forces_contraction_active == true

                %............. Write vertex contraction forces ................
                list_vert_x  = monolayer.contraction_vertices_coord(:,1);
                list_vert_y  = monolayer.contraction_vertices_coord(:,2);
                list_vert_z  = monolayer.contraction_vertices_coord(:,3);
    
                vert_force   = monolayer.contraction_vertices_force;
    
                filename = sprintf('./results/force_vert/cont_f_vert_%05i.vtk', step); 
        
                % Write vertex contraction forces for every vertex and cell
                vtkwrite_cells(   filename                                   , 'POLYDATA'                , list_vert_x, list_vert_y, list_vert_z , ...
                                 'vectors' , 'vertex_force'                  , vert_force                , ...
                                 'Precision', 2);

            end





%             %%%%%%% Write vertex contraction forces (resultant) %%%%%%%%%%%
%             list_vert_x          = monolayer.voronoi_vertices(:,1);
%             list_vert_y          = monolayer.voronoi_vertices(:,2);
%             list_vert_z          = zeros(size(list_vert_x,1),1);
% 
%             vert_force_resultant = monolayer.vertices_forces;
% 
%             filename = sprintf('./results/force_vert/cont_f_vert_result%05i.vtk', step); 
%     
%             % Write vertex contraction force (resultant)
%             vtkwrite_cells(   filename                                   , 'POLYDATA'                , list_vert_x, list_vert_y, list_vert_z , ...
%                              'vectors' , 'vertex_force'                  , vert_force_resultant      , ...
%                              'Precision', 2);





            %............. Write vertex protrusive forces .................

            if inputs.forces_protrusion_active == true

                if monolayer.step > 2
    
                    list_vert_x  = monolayer.protrusion_vertices_coord(:,1);
                    list_vert_y  = monolayer.protrusion_vertices_coord(:,2);
                    list_vert_z  = monolayer.protrusion_vertices_coord(:,3);
        
                    vert_force   = monolayer.protrusion_vertices_force;
        
                    filename = sprintf('./results/force_vert/protr_f_vert_%05i.vtk', step); 
            
                    % Write vertex protrusive forces 
                    vtkwrite_cells(   filename                                   , 'POLYDATA'                , list_vert_x, list_vert_y, list_vert_z , ...
                                     'vectors' , 'vertex_force'                  , vert_force                , ...
                                     'Precision', 2);
    
                else
                    % The first two steps there are no protrusive forces. To
                    % avoid delayed representations in paraview, we copy
                    % contraction forces for the first two steps
    
                    filename = sprintf('./results/force_vert/protr_f_vert_%05i.vtk', step); 
        
                    % Write vertex contraction forces for every vertex and cell
                    vtkwrite_cells(   filename                                   , 'POLYDATA'                , list_vert_x, list_vert_y, list_vert_z , ...
                                     'vectors' , 'vertex_force'                  , vert_force                , ...
                                     'Precision', 2);
    
                end

            end



            








         
        end

    end




end