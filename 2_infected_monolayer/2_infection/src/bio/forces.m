classdef forces


    methods

        %................... Object initialization ........................
        function obj = forces()
            addpath(fullfile('inputs')); 
            addpath(fullfile('src')); 
        end





        %......... Calculate cell-cell interaction forces .................
        function cell_int_force = calculateInteraction(~,iterator,neighbours_id,mon_cell_radius,mon_cell_pos,mon_cell_infection,next_extruding_cell)

            tolerance      = inputs.general_tolerance;
            cell_int_force = [0 0 0];


            for i = 1 : length(neighbours_id)

                % Avoid forces related with the next extruding cell
                if neighbours_id(i) ~= next_extruding_cell

                    % Define equilibrium distance used in interaction forces
                    eq_distance     = mon_cell_radius(iterator)    + ...
                                      mon_cell_radius(neighbours_id(i));
    
                    % Alter equilibrium distance considering cell polarity (not
                    % implemented here)
    
                    % Define interaction mode: right now only Lennard Jones
                    % interaction_mode = ...;
                    
                    % Relative distance between cell A and cell neighbour B
                    position_A      = mon_cell_pos(iterator,:);
                    position_B      = mon_cell_pos(neighbours_id(i),:);
                    distance_vector = position_B - position_A;
                    distance_module = norm(distance_vector);
    
                    % Lennard Jones interaction forces
                    if distance_module > tolerance
    
                        
    
    
                        if inputs.monolayer_infection_active == true 
                            % Only if infection == true, we ask this question:
            
                            % Cell-cell interactions magnitude depending on infection state
                            if inputs.forces_interaction_infection == true
                
                                cell_status_inf       = mon_cell_infection(iterator);
                                if cell_status_inf == 0
                                    Eps_constant  = inputs.forces_interaction_eps_uninf;
                                elseif cell_status_inf == 1
                                    Eps_constant  = inputs.forces_interaction_eps_inf;
                                end
                            else
                                Eps_constant = inputs.forces_interaction_lennard_constant;
                            end
            
            
                        else
            
                            Eps_constant = inputs.forces_interaction_lennard_constant;
            
                        end
    
    
    
    
                        force_module = -12 * Eps_constant *              ...
                            ((eq_distance.^12 ./ distance_module.^13) -  ...
                             (eq_distance.^6 ./ distance_module.^7));
    
                        % Normalize direction using relative position module
                        force = force_module / distance_module * distance_vector;
    
                        cell_int_force = cell_int_force + force;
    
                    end
                end

            end

            % Limit cell-cell interaction forces
            force_module = norm(cell_int_force)           ;
            force_limit  = inputs.forces_interaction_limit;
            if force_module > force_limit

                scale_ratio    = force_limit / force_module  ;
                cell_int_force = scale_ratio * cell_int_force;                    

            end

        end









        %.............. Calculate cell contraction forces .................
        % Contraction forces calculation. Total contraction force is 
        % distributed proportionally on the cell edges using edge length. 
        % Edge force is then distributed on its vertices          
        function [centroid_cont_force, vert_forces, cont_vert_coor, cont_vert_forc] = calculateContraction(~,iterator, ...
                   inf_state, inf_cells_border, step, cell_pos, mon_polygons, mon_regions, mon_vertices, vert_forces, cont_vert_coor, cont_vert_forc ...
                   , first_extrusion, random_contraction, average_cell_stress)


            %%%%%% Contraction mode :  1.- Default contraction
            if inputs.forces_contraction_mode == 1

                if inputs.monolayer_infection_active == true 
                    % Only if infection == true, we ask this question:
    
                    % Contraction magnitude depending on infection state
                    if inputs.forces_contraction_infection == true
        
                        cell_status       = inf_state(iterator);
                        if cell_status == 0
                            contraction_magn  = inputs.forces_contraction_magn_uninf;
                        elseif cell_status == 1
                            contraction_magn  = inputs.forces_contraction_magn_inf;
                        end
                    else
                        contraction_magn  = inputs.forces_contraction_magnitude;
                    end
    
    
                else
    
                    contraction_magn  = inputs.forces_contraction_magnitude;
    
                end





            %%%%%% Contraction mode :  2.- Random contraction 
            elseif inputs.forces_contraction_mode == 2

                contraction_magn = random_contraction(iterator);





            %%%%%% Contraction mode :  3.- Contraction taking into account principal stress
            elseif inputs.forces_contraction_mode == 3

                % At the beginning forces might be high, wait for 15 steps
                if step < 15
                    contraction_magn  = inputs.forces_contraction_magnitude;
                else
                    contraction_magn  = average_cell_stress(iterator)                    ;
                    contraction_magn  = abs(contraction_magn)                            ;
                    contraction_magn  = contraction_magn * inputs.force_contraction_scale;
                    
                end

                % If the cell is infected, increase contraction magnitude
                cell_status       = inf_state(iterator);
                % Remove infected cells at the border
                cell_status       = cell_status - inf_cells_border(iterator);
                factor            = 3;
                
                if cell_status == 1 %&& first_extrusion == false
                    %contraction_magn = contraction_magn * factor;
                    % Constant contraction for infected cells that are not
                    % at the border
                    contraction_magn = inputs.force_contraction_infection;                    
                end

                if contraction_magn > inputs.forces_contraction_limit
                    contraction_magn = inputs.forces_contraction_limit;
                end


            end




                

            
            cell_centroid       = cell_pos(iterator,:)   ;
            cell_polygon        = mon_polygons{iterator} ;
            cell_perimeter      = perimeter(cell_polygon);
            cell_vertices_id    = mon_regions{iterator}  ;

            centroid_cont_force = [0 0]                  ;

            for i = 1 : length(cell_vertices_id)
                % Find edges related to the vertex and sum their lengths
                point_i_id        = cell_vertices_id(i);

                if i == 1
                    point_A_id = cell_vertices_id(end);
                    point_B_id = cell_vertices_id(i+1);
                    
                elseif i == length(cell_vertices_id)
                    point_A_id = cell_vertices_id(i-1);
                    point_B_id = cell_vertices_id(1);
                else
                    point_A_id = cell_vertices_id(i-1);
                    point_B_id = cell_vertices_id(i+1);        
                end

                coord_point_A = mon_vertices(point_A_id,:);
                coord_point_B = mon_vertices(point_B_id,:);
                coord_point_i = mon_vertices(point_i_id,:);

                length_1          = norm(coord_point_A-coord_point_i);
                length_2          = norm(coord_point_i-coord_point_B);
                two_edges_length  = length_1 + length_2;

                % Force is distributed on each vertex considering total 
                % perimeter and the half of the sum of the edge lenghts
                force_module      = (contraction_magn * two_edges_length) / ...
                                    ( 2 * cell_perimeter);

                % Limit forces
                if force_module > inputs.forces_contraction_limit
                    force_module = inputs.forces_contraction_limit;
                end

                force_direction   = cell_centroid(1:2) - coord_point_i;

                % Normalize direction. There are centroids at the edge that
                % are considered vertices. Avoid the vector direction [0 0]
                if isequal(force_direction, [0 0])
                    contraction_force = [0 0];
                else
                    force_direction   = force_direction/norm(force_direction);
                    contraction_force = force_module*force_direction; 
                end

                
                % Update global vertices forces
                vert_forces(point_i_id,1:2) = ...
                    vert_forces(point_i_id,1:2) + contraction_force;

                % Update cell centroid forces
                centroid_cont_force = centroid_cont_force - contraction_force;


                % Save variables to plot later contraction forces (we need 
                % to repeat points to plot every single force)
                coord_point      = [coord_point_i 0]                   ;
                cont_force       = [contraction_force 0]               ;
                cont_vert_coor   = [cont_vert_coor   ; coord_point    ];
                cont_vert_forc   = [cont_vert_forc   ; cont_force     ];

               
            end

        end






        %.............. (1) Calculate cell protrusive forces ..............
        % - Protruding along the major axis of the cell elipsoid (both 
        %   directions +- 180 degrees) 
        % - Only for vertices aligned to the ellipsoid major axis of the 
        %   cell. 
        % - The protrusive force is distributed proportionally on the cell 
        %   edges using edge length.

        function [centroid_prot_force,vert_forces_2, prot_vert_coor, prot_vert_forc] = calculateProtrusion(~, iterator, mon_step, aver_cell_stress, mon_cell_pos...
                                , mon_polygons, mon_regions, orient_axis_elips, mon_vertices, vert_forces_2, prot_vert_coor, prot_vert_forc )

            % Protrusion: magnitude considering principal stress sigma 
            % At the beginning forces mighbt be high, wait for 15 steps
            if mon_step < 15
                protrusion_magn  = inputs.forces_protrusion_magnitude;
            else

                % Considering sigma I ([0 s_Imax]) tensile stresses
                protrusion_magn  = abs(aver_cell_stress(iterator))                 ;
                protrusion_magn  = protrusion_magn * inputs.forces_protrusion_scale;

                if protrusion_magn > inputs.forces_protrusion_limit
                    protrusion_magn = inputs.forces_protrusion_limit;
                end
            end

        
            
            
            cell_centroid       = mon_cell_pos(iterator,:);
            cell_polygon        = mon_polygons{iterator}  ;
            cell_perimeter      = perimeter(cell_polygon) ;
            cell_vertices_id    = mon_regions{iterator}   ;
            centroid_prot_force = [0 0]                   ;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Cell protrusion only for vertices aligned with the major
            % axis of the ellipsoid formed by the cell!!
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Major axis direction (load value)
            axis_direction          = orient_axis_elips(iterator);
            axis_degrees            = axis_direction*180/pi      ;


            % Axis direction goes from -180 to 180 degrees, take it into
            % account so it goes  from    0 to 360 degrees

            if axis_direction < 0
                % Case [-180, 0] degrees
                axis_direction          = axis_direction + 2*pi;
                axis_degrees            = axis_direction*180/pi;

                % We need to check the other orientation of the major
                % axis
                % Difference with the opposite direction (-180 degrees)
                opp_orientation_rad     = axis_direction - pi;
                opp_orientation_deg     = opp_orientation_rad*180/pi;
                
            else
                % Case [0, 180] degrees
                % Difference with the opposite direction (+180 degrees)
                opp_orientation_rad     = axis_direction + pi;
                opp_orientation_deg     = opp_orientation_rad*180/pi;

            end

            % Loop over all vertices of the cell
            for i = 1 : length(cell_vertices_id)


                % Compute vector orientation: centroid-->vertex
                % Find vertex
                point_i_id          = cell_vertices_id(i)       ;
                coord_point_i       = mon_vertices(point_i_id,:);

                % Compute vector: centroid -> vertex
                vector_cent_vert    = coord_point_i - cell_centroid(1:2);

                
                if     vector_cent_vert(1) >= 0 && vector_cent_vert(2) >= 0     % First quadrant

                    direction_rad       = atan(vector_cent_vert(2)/vector_cent_vert(1));
                    direction_deg       = direction_rad*180/pi;

                elseif vector_cent_vert(1) <  0 && vector_cent_vert(2) >  0     % Second quadrant

                    direction_rad       = atan(vector_cent_vert(2)/vector_cent_vert(1)) + pi;
                    direction_deg       = direction_rad*180/pi;

                elseif vector_cent_vert(1) <  0 && vector_cent_vert(2) <  0     % Third quadrant

                    direction_rad       = atan(vector_cent_vert(2)/vector_cent_vert(1)) + pi;
                    direction_deg       = direction_rad*180/pi;

                elseif vector_cent_vert(1) >  0 && vector_cent_vert(2) <  0     % Fourth quadrant

                    direction_rad       = atan(vector_cent_vert(2)/vector_cent_vert(1))  + 2*pi; % It gives [-90,0] degrees, sum 360
                    direction_deg       = direction_rad*180/pi;
                end
                

                % Compute alignment between 
                %               (1) direction centroid -> vertex
                %               (2) direction major axis (two orientations)
                % Range of alignment +- 30 degrees

                % Angle difference
                angle_diff              = abs(direction_rad - axis_direction)     ;
                angle_diff_degrees      = angle_diff*180/pi                       ;
                opposite_angle_diff     = abs(direction_rad - opp_orientation_rad);
                opposite_angle_diff_deg = opposite_angle_diff*180/pi              ;
       
                % Desired range in degrees
                range_deg               = 30                                                             ;  
                range_rad               = deg2rad(range_deg)                                             ;
                is_aligned              = (angle_diff <= range_rad) || (opposite_angle_diff <= range_rad);
    
                % If there is alignment, the protrusive force is calculated
                if is_aligned == true
                
                    % Find edges related to the vertex and sum their lengths
                    point_i_id        = cell_vertices_id(i);
    
                    if i == 1
                        point_A_id = cell_vertices_id(end);
                        point_B_id = cell_vertices_id(i+1);
                        
                    elseif i == length(cell_vertices_id)
                        point_A_id = cell_vertices_id(i-1);
                        point_B_id = cell_vertices_id(1)  ;
                    else
                        point_A_id = cell_vertices_id(i-1);
                        point_B_id = cell_vertices_id(i+1);        
                    end
    
                    coord_point_A     = mon_vertices(point_A_id,:)       ;
                    coord_point_B     = mon_vertices(point_B_id,:)       ;
                    coord_point_i     = mon_vertices(point_i_id,:)       ;
    
                    length_1          = norm(coord_point_A-coord_point_i);
                    length_2          = norm(coord_point_i-coord_point_B);
                    two_edges_length  = length_1 + length_2              ;
    
                    % Force is distributed on each vertex considering total 
                    % perimeter and the half of the sum of the edge lenghts
                    force_module      = (protrusion_magn * two_edges_length) / ...
                                        ( 2 * cell_perimeter);

                    % Limit forces
                    if force_module > inputs.forces_protrusion_limit
                        force_module = inputs.forces_protrusion_limit;
                    end
    
                    force_direction   = cell_centroid(1:2) - coord_point_i;
    
                    % Normalize direction. There are centroids at the edge that
                    % are considered vertices. Avoid the vector direction [0 0]
                    if isequal(force_direction, [0 0])
                        protrusive_force = [0 0];
                    else
                        force_direction   = force_direction/norm(force_direction);
                        protrusive_force  = force_module*force_direction         ;
                    end

                    
                     % Update global vertices forces
                    vert_forces_2(point_i_id,1:2) = ...
                        vert_forces_2(point_i_id,1:2) + protrusive_force;
    
                    % Update cell centroid forces
                    centroid_prot_force = centroid_prot_force - protrusive_force;
    
    
                    % Save variables to plot later protrusive forces (we need 
                    % to repeat points to plot every single force)
                    coord_point       = [coord_point_i 0]                   ;
                    prot_force        = [protrusive_force 0]                ;
                    prot_vert_coor    = [prot_vert_coor   ; coord_point    ];
                    prot_vert_forc    = [prot_vert_forc   ; prot_force     ];

                end

                

            end


        end






        %.............. (2) Calculate cell protrusive forces ..............
        % - Protrusion along the direction of minimal stress sigma I.
        % - Only for vertices that are aligned with that orientation.
        % - The protrusive force is distributed proportionally on the cell 
        %   edges using edge length. 
        % - Nothing related with the cell ellipsoid major axis.

        function [centroid_prot_force,vert_forces_2, prot_vert_coor, prot_vert_forc] = calculateProtrusion_dir_minimal_sigma_I(~, iterator, mon_step, aver_cell_stress, mon_cell_pos...
                                , mon_polygons, mon_regions, orient_sigma_I, mon_vertices, vert_forces_2, prot_vert_coor, prot_vert_forc )

            % Protrusion: magnitude considering principal stress sigma I
            % At the beginning forces mighbt be high, wait for 15 steps
            if mon_step < 15
                protrusion_magn  = inputs.forces_protrusion_magnitude;
            else

                % Considering sigma I ([0 s_Imax]) tensile stresses
                protrusion_magn  = abs(aver_cell_stress(iterator))                 ;
                protrusion_magn  = protrusion_magn * inputs.forces_protrusion_scale;

                if protrusion_magn > inputs.forces_protrusion_limit
                    protrusion_magn = inputs.forces_protrusion_limit;
                end
            end

        
            
            
            cell_centroid       = mon_cell_pos(iterator,:);
            cell_polygon        = mon_polygons{iterator}  ;
            cell_perimeter      = perimeter(cell_polygon) ;
            cell_vertices_id    = mon_regions{iterator}   ;
            centroid_prot_force = [0 0]                   ;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Cell protrusion only for vertices aligned with the direction
            % of minimal stress sigma I!!
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % HERE the difference with the function calculateProtrusion
            % There is only one orientation where protrusion will be
            % applied

            % Direction of minimal stress sigma I(load value)
            axis_direction          = orient_sigma_I(iterator);
            axis_degrees            = axis_direction*180/pi   ;

            % We have already computed this direction within the range [0 2pi]


            % Loop over all vertices of the cell
            for i = 1 : length(cell_vertices_id)


                % Compute vector orientation: centroid-->vertex
                % Find vertex
                point_i_id          = cell_vertices_id(i)        ;
                coord_point_i       = mon_vertices(point_i_id,:) ;

                % Compute vector: centroid -> vertex
                vector_cent_vert    = coord_point_i - cell_centroid(1:2);

                
                if     vector_cent_vert(1) >= 0 && vector_cent_vert(2) >= 0     % First quadrant

                    direction_rad       = atan(vector_cent_vert(2)/vector_cent_vert(1));
                    direction_deg       = direction_rad*180/pi;

                elseif vector_cent_vert(1) <  0 && vector_cent_vert(2) >  0     % Second quadrant

                    direction_rad       = atan(vector_cent_vert(2)/vector_cent_vert(1)) + pi;
                    direction_deg       = direction_rad*180/pi;

                elseif vector_cent_vert(1) <  0 && vector_cent_vert(2) <  0     % Third quadrant

                    direction_rad       = atan(vector_cent_vert(2)/vector_cent_vert(1)) + pi;
                    direction_deg       = direction_rad*180/pi;

                elseif vector_cent_vert(1) >  0 && vector_cent_vert(2) <  0     % Fourth quadrant

                    direction_rad       = atan(vector_cent_vert(2)/vector_cent_vert(1))  + 2*pi; % It gives [-90,0] degrees, sum 360
                    direction_deg       = direction_rad*180/pi;
                end
                

                % Compute alignment between 
                %               (1) direction centroid -> vertex
                %               (2) direction major specific axis (one orientation)
                % Range of alignment +- 40 degrees

                % Angle difference
                angle_diff              = abs(direction_rad - axis_direction);
                angle_diff_degrees      = angle_diff*180/pi                  ;
                    
                % Desired range in degrees
                range_deg               = 40                                 ;  
                range_rad               = deg2rad(range_deg)                 ;
                is_aligned              = (angle_diff <= range_rad)          ;
    
                % If there is alignment, the protrusive force is calculated
                if is_aligned == true
                
                    % Find edges related to the vertex and sum their lengths
                    point_i_id        = cell_vertices_id(i);
    
                    if i == 1
                        point_A_id = cell_vertices_id(end);
                        point_B_id = cell_vertices_id(i+1);
                        
                    elseif i == length(cell_vertices_id)
                        point_A_id = cell_vertices_id(i-1);
                        point_B_id = cell_vertices_id(1);
                    else
                        point_A_id = cell_vertices_id(i-1);
                        point_B_id = cell_vertices_id(i+1);        
                    end
    
                    coord_point_A = mon_vertices(point_A_id,:);
                    coord_point_B = mon_vertices(point_B_id,:);
                    coord_point_i = mon_vertices(point_i_id,:);
    
                    length_1          = norm(coord_point_A-coord_point_i);
                    length_2          = norm(coord_point_i-coord_point_B);
                    two_edges_length  = length_1 + length_2;
    
                    % Force is distributed on each vertex considering total 
                    % perimeter and the half of the sum of the edge lenghts
                    force_module      = (protrusion_magn * two_edges_length) / ...
                                        ( 2 * cell_perimeter);

                    % Limit forces
                    if force_module > inputs.forces_protrusion_limit
                        force_module = inputs.forces_protrusion_limit;
                    end
    
                    force_direction   = cell_centroid(1:2) - coord_point_i;
    
                    % Normalize direction. There are centroids at the edge that
                    % are considered vertices. Avoid the vector direction [0 0]
                    if isequal(force_direction, [0 0])
                        protrusive_force = [0 0];
                    else
                        force_direction   = force_direction/norm(force_direction);
                        protrusive_force  = force_module*force_direction         ;
                    end

                    
                    % Update global vertices forces
                    vert_forces_2(point_i_id,1:2) = ...
                        vert_forces_2(point_i_id,1:2) + protrusive_force;
    
                    % Update cell centroid forces
                    centroid_prot_force = centroid_prot_force - protrusive_force;
    
    
                    % Save variables to plot later contraction forces (we need 
                    % to repeat points to plot every single force)
                    coord_point      = [coord_point_i 0]                   ;
                    prot_force       = [protrusive_force 0]               ;
                    prot_vert_coor   = [prot_vert_coor   ; coord_point    ];
                    prot_vert_forc   = [prot_vert_forc   ; prot_force     ];

    
                end
               
            end

        end


        
        %.............. Calculate cell random walk ........................      
        function monolayer = cell_random_walk(~,monolayer)

            % Customize random generator

            if strcmp(inputs.forces_random_walk_generator,'shuffle')

                rng('shuffle'); % Initializes generator based on the current 
                                % time, resulting in a different sequence of 
                                % random numbers after each call to rng. 
                                % choose default for a fixed random
                                % distribution. shuffle default
            end

            
            if monolayer.step == 1                
                % Random walk
                numb_of_cells                     = numel(monolayer.cell_id);
                random_force                      = randn(numb_of_cells, 2);
                rand_for_normalized               = random_force ./ sqrt(sum(random_force.^2, 2));
                cte_random_walk_force             = rand_for_normalized * inputs.forces_random_walk_magnitude;

                % Apply a factor so not all the cells have the same force
                % magnitude
                minValue                          = 0.7;  % Minimum value
                maxValue                          = 1.3;  % Maximum value
                factor                            = (minValue + (maxValue - minValue) * rand(1, numb_of_cells));
                random_walk_force                 = factor'.*cte_random_walk_force;
                z_component                       = zeros(numb_of_cells,1);
                random_walk_force                 = [random_walk_force z_component];
            
                max_value                         = inputs.forces_random_walk_max_persistance;
                random_persistance                = randi([1 max_value], [numb_of_cells 1]); % randi([min max], size) 

                % Exclude random walk for infected cells
                if inputs.forces_random_walk_infected_cells == false
                    list_infected_cells   = monolayer.cell_infection;   % 1 infected 0 uninfected
                    list_uninfected_cells = ~list_infected_cells    ;   % 0 infected 1 uninfected

                    % Remove values from infected cells
                    random_walk_force     = random_walk_force .*list_uninfected_cells;
                    random_persistance    = random_persistance.*list_uninfected_cells;

                    % Update persistance of infected cells, so it is always
                    % ==1 and it will never go to 0
                    random_persistance    = random_persistance + list_infected_cells*2;
                end
                
                monolayer.random_walk_persistance = random_persistance-1;
                monolayer.force_random_walk       = random_walk_force;
            
                monolayer.force_list              = monolayer.force_list + random_walk_force;
            else
            
                persistance_vec      = monolayer.random_walk_persistance;
                rand_force_vec       = monolayer.force_random_walk;
            
                if any(persistance_vec == 0)
            
                    idx_cero                         = find(persistance_vec == 0); % get index where persistance is zero
            
                    numb_cells_to_change             = numel(idx_cero);
                    random_force                     = randn(numb_cells_to_change, 2);
                    rand_for_normalized              = random_force ./ sqrt(sum(random_force.^2, 2));
                    cte_random_walk_force            = rand_for_normalized * inputs.forces_random_walk_magnitude;

                    % Apply a factor so not all the cells have the same force
                    % magnitude
                    minValue                          = 0.7;  % Minimum value
                    maxValue                          = 1.3;  % Maximum value
                    factor                            = (minValue + (maxValue - minValue) * rand(1, numb_cells_to_change));
                    random_walk_force                 = factor'.*cte_random_walk_force;
                    z_component                       = zeros(numb_cells_to_change,1);
                    random_walk_force                 = [random_walk_force z_component];

                    
                    max_value                        = inputs.forces_random_walk_max_persistance;
                    random_persistance               = randi([1 max_value], [numb_cells_to_change 1]); % randi([min max], size) 
            
                    persistance_vec(idx_cero)        = random_persistance;
                    rand_force_vec(idx_cero,:)       = random_walk_force ;
                    monolayer.force_random_walk      = rand_force_vec    ;
                end

                % Exclude random walk for infected cells
                if inputs.forces_random_walk_infected_cells == false
                    list_infected_cells   = monolayer.cell_infection;   % 1 infected 0 uninfected

                    % Update persistance of infected cells, so it is always
                    % ==1 and it will never go to 0
                    persistance_vec       = persistance_vec + list_infected_cells;
                end
            
                monolayer.random_walk_persistance    = persistance_vec-1;
                monolayer.force_list                 = monolayer.force_list + monolayer.force_random_walk;
                
            
            end
            


        end






        %.............. Calculate crawling infection forces ...............       
        function monolayer = calculateInfectionCrawling(~,iterator,monolayer,mon_regions,border_vertices,factor)

            crawling_magn     = inputs.forces_infection_magnitude * factor   ;
            cell_centroid     = monolayer.cell_position(iterator,:)          ;
            cell_polygon      = monolayer.cellPolyList{iterator}             ;
            cell_perimeter    = perimeter(cell_polygon)                      ;
            cell_vertices_id  = mon_regions{iterator}                        ;


            % Check what are the vertices that are on the border between 
            % infected and uninfected cells

            % (1) Cell vertices
            coord_cell_vertices  = monolayer.voronoi_vertices(cell_vertices_id,:);

            % (2) Retrieve border domain (infection or first lines of
            % cells)
            inf_focus_vertices   = border_vertices                    ;
            
            % (3) Check whether there are vertices on the edge of
            % the infection (1)intersection with (2) 
            [~,on]      = inpolygon(coord_cell_vertices(:,1),coord_cell_vertices(:,2), ... 
                                     inf_focus_vertices(:,1), inf_focus_vertices(:,2));

            % (4) Plot vertices and polygons

%             plot(cell_polygon)
%             hold on
%             plot(infected_polygon)
%             for i = 1 : length(on)
%                 if on(i) == 1
%                     plot(coord_cell_vertices(i,1), coord_cell_vertices(i,2), 'r*', 'MarkerSize', 10);
%                 else
%                     plot(coord_cell_vertices(i,1), coord_cell_vertices(i,2), 'bo', 'MarkerSize', 10);
%                 end
%             end
%             hold off

            % (5) Indices vertices that are on the edge of the infection
            indices_surround_vertices = find(on==1);


            % (6) Finally, the computation of infective crawling forces
            % Compute forces only on the vertices on the edge of the
            % infection

            for i = 1 : length(cell_vertices_id)

                if ismember(i, indices_surround_vertices)
                
                
                    % Find edges related to the vertex and sum their lengths
                    point_i_id        = cell_vertices_id(i);
    
                    if i == 1
                        point_A_id = cell_vertices_id(end);
                        point_B_id = cell_vertices_id(i+1);
                        
                    elseif i == length(cell_vertices_id)
                        point_A_id = cell_vertices_id(i-1);
                        point_B_id = cell_vertices_id(1);
                    else
                        point_A_id = cell_vertices_id(i-1);
                        point_B_id = cell_vertices_id(i+1);        
                    end
    
                    coord_point_A = monolayer.voronoi_vertices(point_A_id,:);
                    coord_point_B = monolayer.voronoi_vertices(point_B_id,:);
                    coord_point_i = monolayer.voronoi_vertices(point_i_id,:);
    
                    length_1          = norm(coord_point_A-coord_point_i);
                    length_2          = norm(coord_point_i-coord_point_B);
                    two_edges_length  = length_1 + length_2;
    
                    % Force is distributed on each vertex considering total 
                    % perimeter and the half of the sum of the edge lenghts
                    force_module      = (crawling_magn * two_edges_length) / ...
                                        ( 2 * cell_perimeter);
    
                    force_direction   = cell_centroid(1:2) - coord_point_i;
    
                    % Normalize direction. There are centroids at the edge that
                    % are considered vertices. Avoid the vector direction [0 0]
                    if isequal(force_direction, [0 0])
                        infective_crawling_force = [0 0];
                    else
                        force_direction   = force_direction/norm(force_direction);
                        infective_crawling_force = force_module*force_direction; 
                    end
    
                    % Update global vertices forces
                    monolayer.vertices_forces(point_i_id,1:2) = ...
                        monolayer.vertices_forces(point_i_id,1:2) + infective_crawling_force;
    
                    % Update cell centroid forces
                    centroid_infective_crawling_force = -infective_crawling_force;
    
                    monolayer.force_list(iterator,1:2)                = ...
                        monolayer.force_list(iterator,1:2) + centroid_infective_crawling_force;
    
                    monolayer.force_inf_crawl_cent_list(iterator,1:2) = ...
                        monolayer.force_inf_crawl_cent_list(iterator,1:2) + centroid_infective_crawling_force;

                end


            end


        end








        %.............. Calculate extrusion forces ........................       
        function monolayer = calculateExtrusionForces(~,monolayer,inf_cell_pos,iterator,extruding_cell_pos)

            extrusion_magn    = inputs.forces_extrusion_magnitude   ;
            tolerance         = inputs.general_tolerance            ;

            % Relative distance between cell A and extruding cell B
            position_A      = inf_cell_pos            ;
            position_B      = extruding_cell_pos      ;
            distance_vector = position_B - position_A ;
            distance_module = norm(distance_vector)   ;

            % Compute force towards the extrunding cell
            if distance_module > tolerance

                force_module    = extrusion_magn                                  ;
                % Normalize direction using relative position module
                extrusion_force = force_module / distance_module * distance_vector;

            end

            % Update cell extruding forces
            monolayer.force_list(iterator,1:2)           = monolayer.force_list(iterator,1:2)           + extrusion_force;

            monolayer.force_extruding_cell(iterator,1:2) = monolayer.force_extruding_cell(iterator,1:2) + extrusion_force;



        end




        
    end




end