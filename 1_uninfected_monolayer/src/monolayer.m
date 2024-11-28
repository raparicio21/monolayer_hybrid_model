classdef monolayer

    properties
        step                  = 0         % Simulation step
        monolayer_size        = [0 0 0]   % Monolayer size may differ from scaffold size

        monolayer_cell_points             % List of cell centroid points
        alpha_shape                       % Alpha shape vertices
        voronoi_vertices                  % Voronoi vertices
        voronoi_region                    % Voronoi region
        vertices_forces                   % Forces in global vertices

        cellPolyList                      % List of cell Voronoi polygons 
        cellLocationList                  % List of cell locations 1 = exterior 0 = interior

        previous_force_list               % List of previous cell centroid forces
        force_list                        % List of cell centroid forces
        force_cell_cell_int_list          % List of cell_cell interaction forces              (AUXILIAR VARIABLE)
        force_cell_cont_cent_list         % List of cell contraction forces (centroid)        (AUXILIAR VARIABLE)
        force_cell_prot_cent_list         % List of cell protrusive forces  (centroid)        (AUXILIAR VARIABLE)
        force_random_walk                 % List of cell random walk forces (centroid)
        random_walk_persistance           % List of cell random walk persistance
        force_inf_crawl_cent_list         % List of cell infective crawling forces (centroid) (AUXILIAR VARIABLE)
        %contraction_persistance           % List of cell persistance for cell contraction magnitude
        contraction_per_cell              % List of cell random contraction

                                          % The next variables are coupled to plot them:
        contraction_vertices_coord        % List of vertex coordinates for every cell vertex  (AUXILIAR VARIABLE)
        contraction_vertices_force        % List of contraction forces for every vertex       (AUXILIAR VARIABLE)
        protrusion_vertices_coord         % List of vertex coordinates for protrusion         (AUXILIAR VARIABLE)
        protrusion_vertices_force         % List of protrusive forces (vertex)                (AUXILIAR VARIABLE)

        cell_radius                       % List of cell virtual radius (agent-based model)
        cell_id                           % List of cell identification numbers
        cell_position                     % List of cell positions, vector          [x y z]
        cell_init_position                % List of cell initial positions, vector  [x y z]
        cell_previous_position            % List of cell previous positions, vector [x y z]
        cell_PBC_position_PBC             % List of cell positions considering PBC  [x y z] taking into account periodic boundary conditions (PBC)
                                          % (cells can move from one side to the other one)

        cell_total_distance               % List of cell total distances walked by every cell (all segments)
        cell_total_displacement           % Position - Initial_position
        cell_velocity                     % List of cell velocities

        cell_extrusion                    % List of cell extrusion:     0 = not extruding, 1 = extruding cell
        index_extruding_cells             % Id cells that have been extruded
        cell_infection                    % List of infected cells:     0 = not extruding, 1 = extruding cell
        infected_polygon                  % Polygon object of infected cells
        infected_polygon_1st              % Polygon object of infected cells + 1st line of neigbouring uninfected cells
        infected_polygon_2nd              % Polygon object of infected cells + 1st line + 2nd line of neigbouring uninfected cells
        surrounder_cell                   % List of neighbouring uninfected cells (1st line)
        surrounder_cell_2nd               % List of neighbouring uninfected cells (2nd line)
        surrounder_cell_3rd               % List of neighbouring uninfected cells (3rd line)
        infection_age                     % Infection age (steps)

        centroid_of_elipsoid_step_i       % List of centroids for ellipsoid around the cell in the current step i
        orientation_axis_elipsoid_i       % List of orientation axis for ellipsoides in the current step i 
        major_axis_length_i               % List of major axis length for the ellipsoid in the current step i
        minor_axis_length_i               % List of minor axis length for the ellipsoid in the current step i

        orientation_axis_elipsoid_prev_i  % List of orientation axis for ellipsoides in the previous step i 
        major_axis_length_prev_i          % List of major axis length for the ellipsoid in the previous step i
        minor_axis_length_prev_i          % List of minor axis length for the ellipsoid in the previous step i

        orientation_min_sigma_I_prev_i    % List of orientation vectors of minimal stress (s_I)  from previous step i
        orientation_min_sigma_II_prev_i   % List of orientation vectors of minimal stress (s_II) from previous step i
        %accumulated_strain_invariant_1    % List of weigthed accumulated strains in every cell   from previous step i
        initial_area                      % List of initial cell area at step 20 (or another one chosen in the code)
        area_change                       % List of area change normalized with respect to the initial cell area

    end




    methods

        %................... Object initialization ........................
        function obj = monolayer()
            addpath(fullfile('src','geom','distribution')); 
            addpath(fullfile('src','geom','shape')       );
            addpath(fullfile('src','geom','cgal')        );
            addpath(fullfile('src','bio')                );
        end


        %................... Set monolayer size ...........................
        function obj = setMonolayerSize(obj,x,y,z)
            arguments
                obj 
                x   {mustBeNumeric}
                y   {mustBeNumeric}
                z   {mustBeNumeric}
            end
            % Check object format
            if ~isa(obj, 'monolayer')
                error('Error: the input is not a monolayer object');
            end

            obj.monolayer_size = [x y z];
        end



        %................... Initialize system............................. 
        function obj = initialize_monolayer(obj)

           my_cell_distribution = cell_distribution();
           obj                  = my_cell_distribution.distribute(obj);
           
        end




        %.......... Wrapper function to generate monolayer geometry ....... 
        function obj = generateGeometry(obj)

            %tic
            %disp('Updating list monolayer points');
            % Update list of monolayer_points
            obj.monolayer_cell_points      = zeros(length(obj.cell_id),2);
            obj.monolayer_cell_points(:,1) = obj.cell_position(:,1)      ;
            obj.monolayer_cell_points(:,2) = obj.cell_position(:,2)      ;
            %toc

            % Update cell virtual shapes and compute alpha shapes (not used
            % anymore, now the border of the monolayer is a rectangle 
            % instead of thevirtual cell shapes)
            %disp('\tCGAL: Generating cell virtual shapes and alpha shape');
            %my_cell_shapes = cell_shapes();
            %obj            = my_cell_shapes.update_shape(obj);

            %tic
            % Update cell extrusion list
            if inputs.monolayer_extrusion_active == true

                if obj.step >= inputs.monolayer_extrusion_initial_step
                    
                    disp('Info: Updating cell extrusion')      ;
                    list_inf_cells  = obj.cell_infection       ;

                    % Only if there are any infected cells at the moment
                    if any(list_inf_cells)
                        infection_index = find(list_inf_cells == 1);

                        if obj.step >= 45
                            obj.cell_extrusion(925) = 1;
                        end

                    end


                end

            end
            %toc

            %tic
            % Generate Voronoi diagram
            disp('Info: Generating Voronoi diagram');
            my_voronoi     = voronoi_object();
            obj            = my_voronoi.generateVoronoi(obj);
            %toc

            %tic
            % Propagate infection
            if inputs.monolayer_infection_propagation == true && obj.step >= inputs.monolayer_infection_step_propagation
                
                disp('Info: Propagate infection');
                % When does infection start
                mon_step             = obj.step                                              ;
                mon_step_start_infec = inputs.monolayer_infection_step_propagation           ;
                step_of_infection    = mon_step - mon_step_start_infec                       ;
                MPI                  = step_of_infection * inputs.general_time_step_duration ; % Minutes Post Infection (MPI)
                % Variables to be used in the function 
                mon_cell_points      = obj.monolayer_cell_points                             ;
                
                % Propagate infection
                my_polygons          = polygons()                                            ;
                [list_inf_cells]     = my_polygons.propagateInfection(mon_cell_points, MPI)  ;


                % Keep old ones in the vector (infected cells), even they 
                % move away from the infection focus area

                % Initialize vector to retain ones
                retain_ones          = false(numel(obj.cell_id), 1)                          ;
                % Update the retain_ones logical vector with the old
                % infection vector
                retain_ones          = retain_ones | (obj.cell_infection == 1)               ;
                
                % Update the column_vector with new ones and retained ones
                list_inf_cells       = list_inf_cells | retain_ones                          ;
                % Update list of infected cells
                obj.cell_infection   = list_inf_cells                                        ;

            end
            %toc


            %tic
            % Convert geometry to polygons
            disp('Info: Generating polygons');
            my_polygons    = polygons();
            obj            = my_polygons.generatePolygons(obj);
            %toc

           
        end



        %.......... Function to generate nuclei plots for PIV ............. 
        function obj = generateNucleiPlots(obj)
            
            % Plot cell centroids
            % Set the visibility of figures to 'off'
            set(groot, 'DefaultFigureVisible', 'off');
        
            FigH           = figure('Position', get(0, 'Screensize'));
            x_nuclei_coord = obj.monolayer_cell_points(:,1);
            y_nuclei_coord = obj.monolayer_cell_points(:,2);
            markerSize     = 75; % Increase or decrease this value to change the marker size
            scatter(x_nuclei_coord, y_nuclei_coord, markerSize, 'filled','k');
            axis equal; % Set equal aspect ratio for the axes
            axis off
            set(gca, 'xticklabels', []);
            set(gca, 'yticklabels', []);
            monolayer_x_size = inputs.monolayer_size_x;
            monolayer_y_size = inputs.monolayer_size_y;
            xlim([0 monolayer_x_size]); ylim([0 monolayer_y_size]);
            
            folderPath     = './results/nuclei/'; % Replace with your desired folder path
            fileName       = sprintf('./nuclei_%04i.tif', obj.step);
        
            %set(gca, 'units', 'normalized'); %Just making sure it's normalized
            Tight = get(gca, 'TightInset');  %Gives you the bording spacing between plot box and any axis labels
                                             %[Left Bottom Right Top] spacing
            NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
            set(gca, 'Position', NewPos);
            exportgraphics(gca,fullfile(folderPath, fileName),'Resolution',100) 
            %saveas(gca, fullfile(folderPath, fileName));
        
            % Restore the visibility of figures to 'on' (optional)
            set(groot, 'DefaultFigureVisible', 'on');
            
           
        end




        %................... Calculate forces .............................
        function obj = calculateForces(obj,FEM_results)

            %tic
            % Compute neighbours (Voronoi adjacency matrix)
            list_neighbours                = obj.getNeighbours();
            %toc

            %tic
            my_forces                      = forces()           ; 

            % Reset forces (force lists)
            obj.previous_force_list        = zeros(size(obj.monolayer_cell_points,1),3);
            obj.force_list                 = zeros(size(obj.monolayer_cell_points,1),3);
            obj.force_cell_cell_int_list   = zeros(size(obj.monolayer_cell_points,1),3);
            obj.force_cell_cont_cent_list  = zeros(size(obj.monolayer_cell_points,1),3);
            obj.force_cell_prot_cent_list  = zeros(size(obj.monolayer_cell_points,1),3);
            obj.force_inf_crawl_cent_list  = zeros(size(obj.monolayer_cell_points,1),3);
            % Vertices forces already reset when computing voronoi_object
            % force_random_walk has persistance in time
            % Reset vertex information (info used for plots in paraview)
            obj.contraction_vertices_coord = [];
            obj.contraction_vertices_force = [];
            obj.protrusion_vertices_coord  = [];
            obj.protrusion_vertices_force  = [];



            % .............................................................
            % Compute cell-cell interaction forces
            % .............................................................          
            number_cells         = length(obj.cell_id);
            % Get variables from object monolayer before starting the loop
            mon_cell_radius      = obj.cell_radius    ;
            mon_cell_pos         = obj.cell_position  ;
            mon_cell_infection   = obj.cell_infection ;
                        
            % Output variables
            vector_f_interaction = zeros(number_cells,3)    ;
            
            % Compute interaction forces
            extruding_cells      = obj.cell_extrusion       ;
            for i = 1:number_cells

                if extruding_cells(i) == 0    % Non-extruding cell            
                    neighbours_id              = list_neighbours{i}'                                                                            ;
                    vector_f_interaction(i,:)  = my_forces.calculateInteraction(i,neighbours_id,mon_cell_radius,mon_cell_pos,mon_cell_infection); 
                end
            end

            % Update cell-cell interaction forces            
            obj.force_list               = vector_f_interaction;
            obj.force_cell_cell_int_list = vector_f_interaction;





            % .............................................................
            % Compute contraction forces
            % .............................................................
            % Modes of contraction:
            %           - (1) Default contraction:
            %                 Contraction magnitude constant for all the
            %                 cells
            %
            %           - (2) Random contraction:
            %                 Contraction magnitude random within a range
            %
            %           - (3) Principal stress proportional to contraction:
            %                 Contraction magnitude depending on the
            %                 cell average principal stress

            % Get variables from object monolayer before starting the loop
            mon_step            = obj.step            ;
            mon_polygons        = obj.cellPolyList    ;
            mon_regions         = obj.voronoi_region  ;
            mon_vertices        = obj.voronoi_vertices;

            %%%% UPDATE MON_REGIONS IN CASE OF CELL EXTRUSION %%%%
            id_extruding_cells     = obj.index_extruding_cells;
            if any(id_extruding_cells)

                num_extruding_cells = numel(id_extruding_cells);
                for k = 1:num_extruding_cells        
                    pos  = id_extruding_cells(k);
                    % Move list
                    mon_regions(pos+1:end+1) = mon_regions(pos:end)    ;
                    mon_regions{pos} = 0;
                    
                end

            end

            % Initialize variables
            vert_forces         = obj.vertices_forces ; % Now it is empty but it has the correct size
            cont_vert_coor      = [];
            cont_vert_forc      = [];

            % Output variables
            vector_f_contraction = zeros(number_cells,2);


            if inputs.forces_contraction_active == true

                %%%%%% Contraction mode :  1.- Default contraction %%%%%%%%
                if inputs.forces_contraction_mode == 1

                    % Compute contraction forces
                    for i = 1:number_cells  
                        if extruding_cells(i) == 0    % Non-extruding cell 
                            [vector_f_contraction(i,:),vert_forces, cont_vert_coor, cont_vert_forc] = my_forces.calculateContraction(i, mon_cell_infection, mon_step, mon_cell_pos, ...
                                    mon_polygons, mon_regions, mon_vertices, vert_forces, cont_vert_coor, cont_vert_forc);
                        end
                    end

                



                %%%%%% Contraction mode :  2.- Random contraction %%%%%%%%%%   
                elseif inputs.forces_contraction_mode == 2

                    if obj.step == 1 
                        % Generate list of contraction forces
                        contraction_magnitude = inputs.forces_contraction_magnitude       ;
                        number_cells          = numel(obj.cell_id)                        ;
                        cont_magn_list        = ones(number_cells,1)*contraction_magnitude;
    
                        % Apply a factor so not all the cells have the same force
                        % magnitude
                        minValue                          = 0.5;  % Minimum value
                        maxValue                          = 1.5;  % Maximum value
                        factor                            = (minValue + (maxValue - minValue) * rand(1, number_cells));
                        contr_force_magnitude             = factor'.*cont_magn_list;
                            
                        % Update value
                        obj.contraction_per_cell          = contr_force_magnitude  ;
                    end

                    random_contraction_cell = obj.contraction_per_cell;

                    % Compute contraction forces 
                    for i = 1:number_cells  
                        if extruding_cells(i) == 0    % Non-extruding cell 
                            [vector_f_contraction(i,:),vert_forces, cont_vert_coor, cont_vert_forc] = my_forces.calculateContraction(i, mon_cell_infection, mon_step, mon_cell_pos, ...
                                    mon_polygons, mon_regions, mon_vertices, vert_forces, cont_vert_coor, cont_vert_forc ...
                                    , random_contraction_cell);
                        end
                    end  



                %%%%%% Contraction mode : 3.- Contraction taking into account principal stress
                elseif inputs.forces_contraction_mode == 3

                    random_contraction_cell = []; % It is not going to be used here but needed
                    % At the beginning forces might be high, wait for 15 steps
                    if obj.step < 15
                        average_cell_stress  = []; % There are no results in the first step
                    else
%                         % Considering sigma I  (within all range)
%                         average_cell_stress     = FEM_results.average_cell_stress(:,1);
%                         % Considering sigma II (within all range)
%                         average_cell_stress     = FEM_results.average_cell_stress(:,2);
%                         % Considering sigma I  ([0 s_Imax])      tensile
%                         average_cell_stress     = FEM_results.average_cell_stress(:,3);
%                         % Considering sigma II ([s_II_max 0])      compressive
%                         average_cell_stress     = FEM_results.average_cell_stress(:,4);

%                         % The following four variables are called stress but they are forces = stress * area
%                         % Considering force_I  = sigma I  * cell_area  (within all range)
%                         average_cell_stress     = FEM_results.average_cell_stress(:,5);
%                         % Considering force_II = sigma II * cell_area  (within all range)
%                         average_cell_stress     = FEM_results.average_cell_stress(:,6);
%                         % Considering force_I  = sigma I  * cell_area   ([0 s_Imax])      tensile
%                         average_cell_stress     = FEM_results.average_cell_stress(:,7);
%                         % Considering force_II = sigma II * cell_area  ([s_II_max 0])      compressive
%                         average_cell_stress     = FEM_results.average_cell_stress(:,8);
% 
                        % Considering averaged cell force_I  from the resultant stress tensor (principal stress  I)
                        average_cell_stress_I      = FEM_results.average_cell_stress(:,9);
                        % Considering averaged cell force_II from the resultant stress tensor (principal stress II)
                        average_cell_stress_II     = FEM_results.average_cell_stress(:,10);
                        % Combine the two last vectors so all cells have
                        % the posibility to contract
                        abs_vector_I        = abs(average_cell_stress_I)                   ;
                        abs_vector_II       = abs(average_cell_stress_II)                  ;
                        average_cell_stress = max(abs_vector_I,abs_vector_II)              ;
                        
                    end

                    % Compute contraction forces 
                    for i = 1:number_cells  
                        if extruding_cells(i) == 0    % Non-extruding cell 
                            [vector_f_contraction(i,:),vert_forces, cont_vert_coor, cont_vert_forc] = my_forces.calculateContraction(i, mon_cell_infection, mon_step, mon_cell_pos, ...
                                    mon_polygons, mon_regions, mon_vertices, vert_forces, cont_vert_coor, cont_vert_forc ...
                                    , random_contraction_cell, average_cell_stress);
                        end
                    end  

                
                end


            end



            % Update global vertices forces
            obj.vertices_forces                    = vert_forces   ;

            % Update cell centroid forces
            obj.force_list(:,1:2)                  = obj.force_list(:,1:2)                + vector_f_contraction;
            obj.force_cell_cont_cent_list(:,1:2)   = obj.force_cell_cont_cent_list(:,1:2) + vector_f_contraction;

            % Update list of vertices and their forces to plot them later
            obj.contraction_vertices_coord         = cont_vert_coor;
            obj.contraction_vertices_force         = cont_vert_forc;








            % .............................................................
            % Compute protrusive forces
            % .............................................................               
            if inputs.forces_protrusion_active == true

               
                % 1st step there are no stress from the previous step
                % 2nd step there might be high forces
                if obj.step > 2

                    % Get variables from object monolayer before starting the loop
%                     % Considering sigma I                        ([0 s_Imax]) tensile stresses
%                     aver_cell_stress  = FEM_results.average_cell_stress(:,3) ;
%                     % Considering force I = sigma I * cell_area  ([0 s_Imax]) tensile stresses
%                     aver_cell_stress  = FEM_results.average_cell_stress(:,7) ; % the variable is called stress but it is force
                     % Considering averaged cell force_I  from the resultant stress tensor (principal stress  I)
                     aver_cell_stress  = FEM_results.average_cell_stress(:,9) ; % the variable is called stress but it is force
                     orient_sigma_I    = obj.orientation_min_sigma_I_prev_i   ;
%                    % Considering averaged cell force_II  from the resultant stress tensor (principal stress  II)
                     %aver_cell_stress  = FEM_results.average_cell_stress(:,10) ; % the variable is called stress but it is force
                     %orient_sigma_I    = obj.orientation_min_sigma_II_prev_i  ; % It is sigma_II but in this way we do not change the name in the next function
                     

                    orient_axis_elips = obj.orientation_axis_elipsoid_prev_i ;
                    
                     
                    

                    % Initialize variables
                    vert_forces_2       = obj.vertices_forces ; % Now it is empty but it has the correct size
                    prot_vert_coor      = [];
                    prot_vert_forc      = [];
        
                    % Output variables
                    vector_f_protrusion = zeros(number_cells,2);
                    
                    % Protrusion law only for uninfected cells
                    cell_status_inf       = obj.cell_infection;

                    for i = 1 : number_cells
                        if cell_status_inf(i) == 0 && extruding_cells(i) == 0 % Uninfected an non-extruding cell
                            % Compute protrusive forces, different ways (do
                            % not use them at the same time, comment!!!)

%                             % (1) Original way, protruding along the major 
%                             % axis of the elipsoid (both directions +- 180
%                             % degrees)
%                             [vector_f_protrusion(i,:),vert_forces_2, prot_vert_coor, prot_vert_forc] = my_forces.calculateProtrusion(i, mon_step, ...
%                                 aver_cell_stress, mon_cell_pos, mon_polygons, mon_regions, orient_axis_elips, ...
%                                 mon_vertices, vert_forces_2, prot_vert_coor, prot_vert_forc);

                            % (2) Protruding along the direction of minimal
                            % stress sigma I, nothing related with the
                            % axis of the ellipsoid
                            [vector_f_protrusion(i,:),vert_forces_2, prot_vert_coor, prot_vert_forc] = my_forces.calculateProtrusion_dir_minimal_sigma_I(i, mon_step, ...
                                aver_cell_stress, mon_cell_pos, mon_polygons, mon_regions, orient_sigma_I, ...
                                mon_vertices, vert_forces_2, prot_vert_coor, prot_vert_forc);
                        end
                    end

                    % Update global vertices forces
                    obj.vertices_forces                    = vert_forces_2   ;
        
                    % Update cell centroid forces
                    obj.force_list(:,1:2)                  = obj.force_list(:,1:2)                + vector_f_protrusion;
                    obj.force_cell_prot_cent_list(:,1:2)   = obj.force_cell_prot_cent_list(:,1:2) + vector_f_protrusion;
        
                    % Update list of vertices and their forces to plot them later
                    obj.protrusion_vertices_coord          = prot_vert_coor;
                    obj.protrusion_vertices_force          = prot_vert_forc;

                end
            end







            % .............................................................
            % Compute cell random walk
            % ............................................................. 
            %tic
            if inputs.forces_random_walk_active == true

                obj           = my_forces.cell_random_walk(obj);

            else

                numb_cells            = numel(obj.cell_id) ;
                list_force_zeros      = zeros(numb_cells,3);
                obj.force_random_walk = list_force_zeros   ;

            end
            %toc





            % .............................................................
            % Compute infection crawling forces
            % ............................................................. 

            %tic
            if inputs.forces_crawling_infection_active == true

                list_surr_cells = obj.surrounder_cell       ;
                surr_indices    = find(list_surr_cells == 1);
                surr_length     = length(surr_indices)      ;

                for i = 1 : surr_length

                    cell_number = surr_indices(i);
                    % Compute infection forces
                    obj = my_forces.calculateInfectionCrawling(cell_number,obj,mon_regions);
                end

            end
            %toc



%             % Plot force magnitude for each polygon
%             hold on
%             for i = 1:length(obj.cell_container)
%                 % Print cell cell interaction forces
%                 X = obj.cell_container(i).polygon.Vertices(:,1);
%                 Y = obj.cell_container(i).polygon.Vertices(:,2);
%                 patch(X,Y,norm(obj.cell_container(i).force_cell_cell_int),'FaceColor','flat')
%                 colormap('parula');
%                 colorbar
%             end            
%             plot(obj.monolayer_cell_points(:,1),obj.monolayer_cell_points(:,2),'.')
%             hold off

        end







        %................... Get cell neighbours ..........................
        function neighbour_list = getNeighbours(monolayer)

            % Load Voronoi regions (rows = cells, columns = vertex that
            % belong to that cell)
            C          = monolayer.voronoi_region       ;

            % Find the cell with more neighbours in the system
            max_length = max(cellfun(@length, C))       ;
            
            % Convert the structure of Voronoi regions (cell structure) 
            % into a matrix. To do that, we need to write zeros in some 
            % positions because every cell has different number of
            % neighbours

            % Initialize structure
            vor_region = zeros(length(C), max_length)   ;

            % Convert cell structure into a matrix structure
            for i = 1:length(C)
                vector     = C{i}                       ; % i_th cell neighbours
                vor_region(i, 1:length(vector)) = vector; % Copy them into the new matrix
            end


            % Initialize list of neighbours
            neighbour_list         = cell(length(C),1)              ;
            id_extruding_cells     = monolayer.index_extruding_cells;
            
            for i = 1: length(C) % Over all the cells
                % Initialize list of neighbours id
                neighbours = [];
                for j = 1:max_length % Over all the vertices
                    % Get id of the vertex in that cell
                    id_vertex = vor_region(i,j);
                    % Remember that there are vertices with id 0 when we
                    % converted the structure Voronoi regions into a
                    % matrix. id = 0 means there is no vertex, we skip
                    % those positions
                    if id_vertex ~= 0 
                        rows_with_that_id = any(vor_region == id_vertex, 2);
                        positions         = find(rows_with_that_id)        ;
                        neighbours        = [neighbours; positions]        ;
                    end
                end
                % Clean neighours that are duplicated
                neighbours             = unique(neighbours)              ;
                % Remove the id corresponding to the same cell
                positions_to_remove    = neighbours == i                 ;
                new_vector             = neighbours(~positions_to_remove);




                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                %%%%%%% Update indices in case there was extrusion %%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                if any(id_extruding_cells)
    
                    num_extruding_cells = numel(id_extruding_cells);
    
                    for k = 1:num_extruding_cells
                        
                        pos                         = id_extruding_cells(k)          ;
                        find_ids_higher             = (new_vector>=pos)              ;
                        if any(find_ids_higher)
                            new_vector(find_ids_higher) = new_vector(find_ids_higher) + 1;
                        end
       
                    end
                end





                neighbour_list{i} = new_vector;

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            %%%%%%%% Update dimension in case there was extrusion %%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if any(id_extruding_cells)
                num_extruding_cells = numel(id_extruding_cells);

                for k = 1:num_extruding_cells
                    
                    pos                         = id_extruding_cells(k)      ;
                    % Move list
                    neighbour_list(pos+1:end+1) = neighbour_list(pos:end)    ;
                    neighbour_list{pos} = 0;
                    
   
                end
            end
           
        
        end




        %................... Set displacements ............................
        function monolayer = setDisplacements(monolayer,fem_results,mesh_object)

            cell_number          = size(monolayer.monolayer_cell_points,1);
            if strcmp(inputs.fem_software,'pde_matlab')
                node_displacements_x = fem_results.results_pde.Displacement.x;
                node_displacements_y = fem_results.results_pde.Displacement.y;
                
            elseif strcmp(inputs.fem_software,'abaqus')
                node_displacements_x = fem_results.displ_results(:,2);
                node_displacements_y = fem_results.displ_results(:,3);
            end

            node_displacement    = [node_displacements_x node_displacements_y];
            node_list            = mesh_object.triangulation_obj.Points       ;

            % Load extruding cell list
            list_extruding_cells = monolayer.cell_extrusion;
            
            for i = 1:cell_number


                if list_extruding_cells(i) == 0 % Non-extruding cells

                    % Cell current position
                    position           = monolayer.monolayer_cell_points(i,:) ;
        
                    % Store position so one can use it in the following step
                    monolayer.cell_previous_position(i,:) = [position 0]      ;
                    
    
                    % Get node displacements
                    [list, index]      = ismember(position, node_list, 'rows');
                    if ~any(list)
                        error('Node coordinates not found while passing AMB-FEM data');
                    end
    
                    % Calculate new position
                    local_displacement = node_displacement(index,:)                       ;
                    new_position       = position + local_displacement                    ;
                    
                    tolerance          = inputs.general_tolerance                         ;
                    boundary           = [inputs.monolayer_size_x inputs.monolayer_size_y];
    
    
    
                    % In case of Periodic Bounday Condition, cells can move
                    % from one side to the opposite side, material_around condition too
                    boundary_condition = inputs.getBoundaryCondition(monolayer.step)      ;
    
                    if strcmp(boundary_condition,'PBC') || strcmp(boundary_condition,'material_around')
    
                        monolayer.cell_PBC_position_PBC(i,:) = [new_position 0] ;
                        out_of_border_position               = new_position     ;
    
                        % X POSITION
                        if new_position(1) < 0                           
                                new_position(1) = boundary(1) - abs(new_position(1) );
                        elseif new_position(1) > boundary(1)
                                new_position(1) = new_position(1)-boundary(1);
                        end
                        
                        % Y POSITION
                        if new_position(2) < 0
                                new_position(2) = boundary(2) - abs(new_position(2) );
                        elseif new_position(2) > boundary(2)
                                new_position(2) = new_position(2)-boundary(2);
                        end
                        
    
                    end
                    
    
                    for j = 1:numel(new_position)
                        if new_position(j) < tolerance
                            new_position(j) = 0;
                        elseif boundary(j) - new_position(j) < tolerance
                            new_position(j) = boundary(j);
                        end
                    end

                else % Extruding cell

                    % Cell current position
                    position           = monolayer.monolayer_cell_points(i,:) ;
        
                    % Store position so one can use it in the following step
                    monolayer.cell_previous_position(i,:) = [position 0]      ;

                    % We keep the same cell position
                    new_position       = position                             ;

                end
                

                % Update position
                monolayer.cell_position(i,:) = [new_position 0];
                %monolayer.monolayer_cell_points(i,:) = new_position;


                % Calculate new total displacement and store total distance
                if strcmp(boundary_condition,'PBC') || strcmp(boundary_condition,'material_around') 
                    % Some cells might go to one side to another side

                    initial_position   = monolayer.cell_init_position(i,1:2)    ;
                    total_displacement = norm(new_position - initial_position)  ;
    
                    old_total_distance = monolayer.cell_total_distance(i)       ;
                    real_displacement  = norm(out_of_border_position - position); 
                    new_total_distance = old_total_distance + real_displacement ;

                else

                    initial_position   = monolayer.cell_init_position(i,1:2)    ;
                    total_displacement = norm(new_position - initial_position)  ;
    
                    old_total_distance = monolayer.cell_total_distance(i)       ;
                    new_displacement   = norm(new_position - position)          ;
                    new_total_distance = old_total_distance + new_displacement  ;

                end

                % Update values
                monolayer.cell_total_displacement(i) = total_displacement;
                monolayer.cell_total_distance(i)     = new_total_distance;





                
            end

        end





        %................... Compute Average Stress .......................
        function [monolayer,fem_results] = computeAverageStress(monolayer,fem_results,mesh_object)

            numb_cells                 = numel(monolayer.cell_id);
            mean_sI_cell_all_range     = zeros(numb_cells,1)     ;
            mean_sII_cell_all_range    = zeros(numb_cells,1)     ;
            mean_sI_cell_positive_val  = zeros(numb_cells,1)     ;
            mean_sII_cell_negative_val = zeros(numb_cells,1)     ;
            mean_fI_cell_all_range     = zeros(numb_cells,1)     ;
            mean_fII_cell_all_range    = zeros(numb_cells,1)     ;
            mean_fI_cell_positive_val  = zeros(numb_cells,1)     ;
            mean_fII_cell_negative_val = zeros(numb_cells,1)     ;

            resultant_fI_cell          = zeros(numb_cells,1)     ;
            resultant_fII_cell         = zeros(numb_cells,1)     ;

            % weighted_e_I               = zeros(numb_cells,1)     ;
            % weighted_e_II              = zeros(numb_cells,1)     ;
            % weighted_e_III             = zeros(numb_cells,1)     ;
            % weighted_inv               = zeros(numb_cells,1)     ;
            cell_reference_area        = zeros(numb_cells,1)     ; 
            change_of_area             = zeros(numb_cells,1)     ;

            step_i           =  monolayer.step    ;
            if step_i > 18
                A_reference = monolayer.initial_area;
            end

            % Reset orientation of minimal stress from previous step
            monolayer.orientation_min_sigma_I_prev_i  = zeros(numb_cells,1);
            monolayer.orientation_min_sigma_II_prev_i = zeros(numb_cells,1);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            %%%%%%%%% Update indices in case there was extrusion %%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            % Cell ID might need to be updated because of extruding cells
            % Load FEM element label related to the ID of the cell
            list_element_labels  = mesh_object.elem_cell_id       ;
            % Load extruding cell list
            list_extruding_cells = monolayer.cell_extrusion       ;
            id_extruding_cells   = monolayer.index_extruding_cells;

            if(any(list_extruding_cells))

                num_extruding_cells = sum(list_extruding_cells);

                for k = 1:num_extruding_cells
                    
                    pos                         = id_extruding_cells(k)              ;
                    find_ids_higher             = (list_element_labels>=pos)         ;
                    if any(find_ids_higher)
                        list_element_labels(find_ids_higher) = list_element_labels(find_ids_higher) + 1;
                    end
   
                end

            end

                           
            infection_state             = monolayer.cell_infection           ;
            poisson_ratio               = inputs.fem_properties_poisson_ratio;

            for i = 1:numb_cells

                if list_extruding_cells(i) == 0 % Non-extruding cell

                    % Find elements that belong to that cell
                    positions      = find(list_element_labels==i)       ;
    
                    % Get area related to those elements
                    element_areas  = mesh_object.triang_areas(positions);
                    sum_areas_cell = sum(element_areas)                 ;
    
    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % (A) Magnitude of principal stresses
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % There are two ways to compute the magnitude:
                    %   - (A.1) Get the principal values of all the triangles
                    %   of the cell and average and weight them. 
                    %   - (A.2) Sum all the components (sigma_x sigma_y tau_xy)
                    %   of every triangle of the cell and then compute the
                    %   principal values of the resultant
    
    
                    % Get stress in the i_th cell
                    % Column 1 : Integration_point_coord_x 
                    % Column 2 : Integration_point_coord_y
                    % Column 3 : S_x 
                    % Column 4 : S_y 
                    % Column 5 : S_xy 
                    % Column 6 : S_maxprin or sigma_I  (range [0  s_I_max]) positive values
                    % Column 7 : S_minprin or sigma_II (range [s_II_max 0]) negative values
                    % Column 8 : S_maxprin (all range)
                    % Column 9 : S_minprin (all range)
    
                    %%%%%%%%%%%%%%%%%%%%%%%%%% (A.1) %%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Within all range of principal stresses
                    sigma_I_cell                  = fem_results.stress_results(positions,8)     ;
                    sigma_II_cell                 = fem_results.stress_results(positions,9)     ;
                    weighted_sigma_I_cell         = sigma_I_cell.*element_areas                 ;
                    weighted_sigma_II_cell        = sigma_II_cell.*element_areas                ;
                    sI_average_cell               = sum(weighted_sigma_I_cell )/sum_areas_cell  ;
                    sII_average_cell              = sum(weighted_sigma_II_cell)/sum_areas_cell  ;                               
                    % Old mean without weights (areas):
                    %sI_average_cell               = mean(sigma_I_cell) ; 
                    %sII_average_cell              = mean(sigma_II_cell);
                    % Save values
                    mean_sI_cell_all_range(i)     = sI_average_cell                             ;
                    mean_sII_cell_all_range(i)    = sII_average_cell                            ;
                    % Stress = Force/area
                    % Here stress is saved, but we might be interested on the
                    % cell averaged force, let us save  force = stress*area as
                    % well
                    mean_fI_cell_all_range(i)     = sum(weighted_sigma_I_cell )                 ;
                    mean_fII_cell_all_range(i)    = sum(weighted_sigma_II_cell)                 ;
    
    
    
    
                    % Only positive (sigma_I) and negative (sigma_II) values of principal stresses
                    sigma_I_cell                  = fem_results.stress_results(positions,6)     ;
                    sigma_II_cell                 = fem_results.stress_results(positions,7)     ;
                    weighted_sigma_I_cell         = sigma_I_cell.*element_areas                 ;
                    weighted_sigma_II_cell        = sigma_II_cell.*element_areas                ;
                    sI_average_cell               = sum(weighted_sigma_I_cell )/sum_areas_cell  ;
                    sII_average_cell              = sum(weighted_sigma_II_cell)/sum_areas_cell  ;
                     % Old mean without weights (areas):
                    %sI_average_cell               = mean(sigma_I_cell) ;
                    %sII_average_cell              = mean(sigma_II_cell);
                    % Save values
                    mean_sI_cell_positive_val(i)  = sI_average_cell                             ;
                    mean_sII_cell_negative_val(i) = sII_average_cell                            ;
                    % Stress = Force/area
                    % Here stress is saved, but we might be interested on the
                    % cell averaged force, let us save  force = stress*area as
                    % well
                    mean_fI_cell_positive_val(i)  = sum(weighted_sigma_I_cell )                 ;
                    mean_fII_cell_negative_val(i) = sum(weighted_sigma_II_cell)                 ;
    
    
    
    
                    %%%%%%%%%%%%%%%%%%%%%%%%%% (A.2) %%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Get local sigma values from all the triangles that belong
                    % to that i_th cell
                    sigma_x_cell                  = fem_results.stress_results(positions,3)     ;
                    sigma_y_cell                  = fem_results.stress_results(positions,4)     ;
                    sigma_xy_cell                 = fem_results.stress_results(positions,5)     ;
    
                    % Weight them
                    weighted_sigma_x_cell         = sigma_x_cell .*element_areas                ;
                    weighted_sigma_y_cell         = sigma_y_cell .*element_areas                ;
                    weighted_sigma_xy_cell        = sigma_xy_cell.*element_areas                ;
    
                    total_weighted_sigma_x_cell   = sum(weighted_sigma_x_cell )/sum_areas_cell  ;
                    total_weighted_sigma_y_cell   = sum(weighted_sigma_y_cell )/sum_areas_cell  ;
                    total_weighted_sigma_xy_cell  = sum(weighted_sigma_xy_cell)/sum_areas_cell  ;
    
                    % Compute eigen vectors 
                    stress_matrix                 = [total_weighted_sigma_x_cell    total_weighted_sigma_xy_cell  ;  total_weighted_sigma_xy_cell  total_weighted_sigma_y_cell];
                    [~, prin_values]              = eig(stress_matrix,'vector')                 ;
    
                    % Sort eigenvectors
                    if prin_values(1) >= prin_values(2)
                        sigma_I_resultant  = prin_values(1);
                        sigma_II_resultant = prin_values(2);
                    else
                        sigma_I_resultant  = prin_values(2);
                        sigma_II_resultant = prin_values(1);
                    end

                    % % Compute principal strains
                    % inf_label     =  infection_state(i)                 ;
                    % 
                    % if inf_label == 0
                    %     young_mod             = inputs.fem_properties_E_uninfected;
                    % elseif inf_label == 1
                    %     young_mod             = inputs.fem_properties_E_infected  ;
                    % end
                    % 
                    % e_I               = (1/young_mod) * (sigma_I_resultant  - poisson_ratio * sigma_II_resultant)                  ;
                    % e_II              = (1/young_mod) * (sigma_II_resultant - poisson_ratio * sigma_I_resultant )                  ;
                    % e_III             = (1/young_mod) * (-poisson_ratio) * (sigma_I_resultant + sigma_II_resultant)                ;
                    % invariant_1       = e_I + e_II                                                                                 ; % area change
                    % %invariant_1       = e_I + e_II + e_III                                                                         ; % volumetric change
                    % 
                    % weighted_e_I(i)   = e_I        ;
                    % weighted_e_II(i)  = e_II       ;
                    % weighted_e_III(i) = e_III      ;
                    % weighted_inv(i)   = invariant_1;
    
                    if sigma_I_resultant < 0
                        sigma_I_resultant = 0;
                    end
                    if sigma_II_resultant > 0
                        sigma_II_resultant = 0;
                    end
    
                    f_I_resultant                 = sigma_I_resultant  *sum_areas_cell          ;
                    f_II_resultant                = sigma_II_resultant *sum_areas_cell          ;
    
                    resultant_fI_cell(i)          = f_I_resultant                               ;
                    resultant_fII_cell(i)         = f_II_resultant                              ;


                    % Compute area change from step 18
                    step_i           =  monolayer.step    ;
                    if step_i == 18
                        cell_reference_area(i) = sum_areas_cell;

                    elseif step_i > 18
                        change_of_area(i)      = sum_areas_cell / A_reference(i);
                    end
    
          
    
    
    
    
    
    
    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % (B) Direction of minimal principal stress
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % .... Compute the direction in which the      .....
                    % .... cell experiences the minimum σI stress  .....
    
                    % Get the cell triangle(s) with minimal stress values
                    sigma_I_cell                  = fem_results.stress_results(positions,6)     ;
                    minValue                      = min(sigma_I_cell)                           ;
                    minIndices                    = sigma_I_cell == minValue                    ;
    
                    cell_centroid                 = monolayer.cell_position(i,:)                                             ;
                    global_index_min_values       = positions(minIndices)                                                    ;
                    triangle_id_points            = mesh_object.triangulation_obj.ConnectivityList(global_index_min_values,:);
    
                    % Get orientation vector between centroid and those
                    % triangles with minimal stress
                    
                    if size(triangle_id_points(:,1)) == 1
                        % If there is only one triangle with minimal stress:
    
                        % Get the vertices of the triangle (1st column is the
                        % centroid, 2nd and 3rd columns are the vertices)
                        id_vert_A      = triangle_id_points(2)                            ;
                        id_vert_B      = triangle_id_points(3)                            ;
                        coor_vert_A    = mesh_object.triangulation_obj.Points(id_vert_A,:);
                        coor_vert_B    = mesh_object.triangulation_obj.Points(id_vert_B,:);
    
                        % Compute middle point between both vertices
                        middle_point   = (coor_vert_A + coor_vert_B)/2                    ;
                        % Compute vector: centroid -> middle point
                        vect_centr_mid = middle_point - cell_centroid(1:2)                ;
                    
                    else
    
                         % If there are more than one triangle with minimal stress:
                         
                         % Get the vertices of the triangles (1st column is the
                         % centroid, 2nd and 3rd columns are the vertices)
                         id_vert_A      = triangle_id_points(:,2)                          ;
                         id_vert_B      = triangle_id_points(:,3)                          ;
                         coor_vert_A    = mesh_object.triangulation_obj.Points(id_vert_A,:);
                         coor_vert_B    = mesh_object.triangulation_obj.Points(id_vert_B,:);
                        
                         % Compute middle points between both vertices for
                         % every triangle
                         middle_point   = (coor_vert_A + coor_vert_B)/2                    ;
    
                         % Compute mean of middle points to get an average point of
                         % minimal stress
                         average_point  = mean(middle_point)                               ;
    
                         % Compute vector: centroid -> middle average point
                         vect_centr_mid = average_point - cell_centroid(1:2)               ;
    
                    end
    
                    % Get orientation (angle)
                    if     vect_centr_mid(1) >= 0 && vect_centr_mid(2) >= 0     % First quadrant
                    
                     direction_rad       = atan(vect_centr_mid(2)/vect_centr_mid(1));
                     direction_deg       = direction_rad*180/pi;
                    
                    elseif vect_centr_mid(1) <  0 && vect_centr_mid(2) >  0     % Second quadrant
                    
                     direction_rad       = atan(vect_centr_mid(2)/vect_centr_mid(1)) + pi;
                     direction_deg       = direction_rad*180/pi;
                    
                    elseif vect_centr_mid(1) <  0 && vect_centr_mid(2) <  0     % Third quadrant
                    
                     direction_rad       = atan(vect_centr_mid(2)/vect_centr_mid(1)) + pi;
                     direction_deg       = direction_rad*180/pi;
                    
                    elseif vect_centr_mid(1) >  0 && vect_centr_mid(2) <  0     % Fourth quadrant
                    
                     direction_rad       = atan(vect_centr_mid(2)/vect_centr_mid(1))  + 2*pi; % It gives [-90,0] degrees, sum 360
                     direction_deg       = direction_rad*180/pi;
                    end
                    
                    monolayer.orientation_min_sigma_I_prev_i(i) = direction_rad;
    
    
    
    
    
    
    
                    % .... Compute the direction in which the       .....
                    % .... cell experiences the minimum σII stress  .....
    
                    % Get the cell triangle(s) with minimal stress values
                    sigma_II_cell                 = abs(fem_results.stress_results(positions,7)) ;
                    minValue                      = min(sigma_II_cell)                           ;
                    minIndices                    = sigma_II_cell == minValue                    ;
    
                    cell_centroid                 = monolayer.cell_position(i,:)                                             ;
                    global_index_min_values       = positions(minIndices)                                                    ;
                    triangle_id_points            = mesh_object.triangulation_obj.ConnectivityList(global_index_min_values,:);
    
                    % Get orientation vector between centroid and those
                    % triangles with minimal stress
                    
                    if size(triangle_id_points(:,1)) == 1
                        % If there is only one triangle with minimal stress:
    
                        % Get the vertices of the triangle (1st column is the
                        % centroid, 2nd and 3rd columns are the vertices)
                        id_vert_A      = triangle_id_points(2)                            ;
                        id_vert_B      = triangle_id_points(3)                            ;
                        coor_vert_A    = mesh_object.triangulation_obj.Points(id_vert_A,:);
                        coor_vert_B    = mesh_object.triangulation_obj.Points(id_vert_B,:);
    
                        % Compute middle point between both vertices
                        middle_point   = (coor_vert_A + coor_vert_B)/2                    ;
                        % Compute vector: centroid -> middle point
                        vect_centr_mid = middle_point - cell_centroid(1:2)                ;
                    
                    else
    
                         % If there are more than one triangle with minimal stress:
                         
                         % Get the vertices of the triangles (1st column is the
                         % centroid, 2nd and 3rd columns are the vertices)
                         id_vert_A      = triangle_id_points(:,2)                          ;
                         id_vert_B      = triangle_id_points(:,3)                          ;
                         coor_vert_A    = mesh_object.triangulation_obj.Points(id_vert_A,:);
                         coor_vert_B    = mesh_object.triangulation_obj.Points(id_vert_B,:);
                        
                         % Compute middle points between both vertices for
                         % every triangle
                         middle_point   = (coor_vert_A + coor_vert_B)/2                    ;
    
                         % Compute mean of middle points to get an average point of
                         % minimal stress
                         average_point  = mean(middle_point)                               ;
    
                         % Compute vector: centroid -> middle average point
                         vect_centr_mid = average_point - cell_centroid(1:2)               ;
    
                    end
    
                    % Get orientation (angle)
                    if     vect_centr_mid(1) >= 0 && vect_centr_mid(2) >= 0     % First quadrant
                    
                     direction_rad       = atan(vect_centr_mid(2)/vect_centr_mid(1));
                     direction_deg       = direction_rad*180/pi;
                    
                    elseif vect_centr_mid(1) <  0 && vect_centr_mid(2) >  0     % Second quadrant
                    
                     direction_rad       = atan(vect_centr_mid(2)/vect_centr_mid(1)) + pi;
                     direction_deg       = direction_rad*180/pi;
                    
                    elseif vect_centr_mid(1) <  0 && vect_centr_mid(2) <  0     % Third quadrant
                    
                     direction_rad       = atan(vect_centr_mid(2)/vect_centr_mid(1)) + pi;
                     direction_deg       = direction_rad*180/pi;
                    
                    elseif vect_centr_mid(1) >  0 && vect_centr_mid(2) <  0     % Fourth quadrant
                    
                     direction_rad       = atan(vect_centr_mid(2)/vect_centr_mid(1))  + 2*pi; % It gives [-90,0] degrees, sum 360
                     direction_deg       = direction_rad*180/pi;
                    end
                    
                    monolayer.orientation_min_sigma_II_prev_i(i) = direction_rad;

                end



            end

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


            fem_results.average_cell_stress = [mean_sI_cell_all_range    mean_sII_cell_all_range    ...
                                               mean_sI_cell_positive_val mean_sII_cell_negative_val ...
                                               mean_fI_cell_all_range    mean_fII_cell_all_range    ...
                                               mean_fI_cell_positive_val mean_fII_cell_negative_val ...
                                               resultant_fI_cell         resultant_fII_cell         ];%...
                                               % weighted_e_I              weighted_e_II              ...
                                               % weighted_e_III            weighted_inv                 ];

            % step_i                   = monolayer.step                                ;
            % if step_i < 2
            %     old_accumulated_strain   = zeros(numb_cells,1)                           ;
            % else
            %     old_accumulated_strain   = monolayer.accumulated_strain_invariant_1      ; 
            % end
            % 
            % new_invariant            = abs(fem_results.average_cell_stress(:,14))    ;
            % accumulated_strain       = old_accumulated_strain + new_invariant        ;
            % monolayer.accumulated_strain_invariant_1 = accumulated_strain            ; 
            

            % Save major and minor axis of elipsoids for next step
            monolayer.orientation_axis_elipsoid_prev_i = monolayer.orientation_axis_elipsoid_i;
            monolayer.major_axis_length_prev_i         = monolayer.major_axis_length_i        ;    
            monolayer.minor_axis_length_prev_i         = monolayer.minor_axis_length_i        ;  

            % Save reference area and normalized area change
            step_i           =  monolayer.step    ;
            if step_i == 18
                monolayer.initial_area = cell_reference_area;
            elseif step_i > 18
                monolayer.area_change  = change_of_area     ;
            end




        end





    end





end