classdef cell_distribution

    methods (Static) % Static because cell_distribution contains a set of 
                     % functions but it does not change parameters of that
                     % class. In this way, the object cell_distribution is
                     % not exposed when calling its functions

                     

        %................... Object initialization.........................
        function obj = cell_distribution()
            addpath(fullfile('inputs')); 
            addpath(fullfile('src'));
        end



        %................... Select distribution mode......................
        % Select distribution mode and apply shake cell, hole and infection
        % to the cell monolayer
        function monolayer = distribute(monolayer)

            % Distribution mode
            if strcmp(inputs.monolayer_distribution_mode,'fill')
                monolayer = cell_distribution.fillDistribution(monolayer);
            else
                error('Error: Invalid cell distribution mode')
            end

            % Distort the current distribution using the magnitude given
            if inputs.monolayer_shake_active == true
                monolayer = cell_distribution.shakeAllCells(monolayer);
            end

            % Infect cells
            if inputs.monolayer_infection_active == true
                if inputs.monolayer_infection_propagation == false
                    % Fixed infection, applied during all the steps
                    monolayer = cell_distribution.infectCells(monolayer);
                else
                    % Infection will start after some steps, so no need to
                    % infect now the cells
                    monolayer.cell_infection = zeros(numel(monolayer.cell_id),1);
                end
            else
                monolayer.cell_infection = zeros(numel(monolayer.cell_id),1);
            end


            % Initialize list of extruding cells and other structures
            monolayer.cell_extrusion                 = zeros(numel(monolayer.cell_id),1);
            monolayer.infection_age                  = zeros(numel(monolayer.cell_id),1);
            monolayer.extrusion_counter              = -1                               ; % Initialize counter
            monolayer.first_extrusion                = false                            ; % there is no previous extrusion at the beginning
            monolayer.first_evaluation_extrusion     = false                            ; % you have not evaluated extrusion yet
            monolayer.total_infection_area           = 0                                ; % initialize infection area

        end




        %................... Distribute cells to fill all the domain.......
        
        function monolayer = fillDistribution(monolayer)

            % Get configuration data
            radius           = inputs.geometry_radius;
            % Number of cells for each dimension
            number_cells     = [0 0];
            number_cells(1)  = floor( monolayer.monolayer_size(1) / (2*radius) );
            number_cells(2)  = floor( monolayer.monolayer_size(2) / (2*radius) );
            % Set coordinates to position conversion step
            position_step    = [0 0];
            position_step(1) = monolayer.monolayer_size(1) / number_cells(1);
            position_step(2) = monolayer.monolayer_size(2) / number_cells(2);
            % Create cell objects knowing total cell number
            total_cells      = number_cells(1) * number_cells(2);
            monolayer        = cell_distribution.addCells(total_cells,monolayer);
            disp(['Info: Filling surface with ' num2str(total_cells) ' cells']);
            % Set position by ID in the grid
            for i = 1:total_cells

                cell_id     = monolayer.cell_id(i);
                % Compute coordinates in the grid
                column      = floor( (cell_id-1) / number_cells(1) );
                row         = cell_id - column * number_cells(1);
                % Set position from grid coordinates
                position    = [0 0];
                offset      = 0.5;
                if inputs.monolayer_distribution_hexagonal == true
                    offset    = 0.16666666;
                end
                position(1) = position_step(1) * ( row    - offset );
                position(2) = position_step(2) * ( column + 0.5    );
                % Move cells in odd rows
                if mod(column, 2) ~= 0 && inputs.monolayer_distribution_hexagonal == true
                    position(1) = position(1) - monolayer.cell_radius(i);
                end

                monolayer.cell_position(i,:)          = [position 0];
                monolayer.cell_init_position(i,:)     = [position 0];
                monolayer.cell_previous_position(i,:) = [position 0];

                if strcmp(inputs.fem_BC_type,'PBC') || strcmp(inputs.fem_BC_type,'material_around')
                    % Boundary conditions: PBC or material around
                    % Cells can move to one side to another side, we need a
                    % copy of the position
                    monolayer.cell_PBC_position_PBC(i,:) = [position 0];
                end

            end

            % Update list of monolayer_points
            for i = 1:total_cells
                monolayer.monolayer_cell_points(i,1) = monolayer.cell_position(i,1);
                monolayer.monolayer_cell_points(i,2) = monolayer.cell_position(i,2);               
            end
                     
        end





        %............. Add cells to the cell monolayer object .............
        function monolayer = addCells(quantity,monolayer)

            % Properties:
            %       - cell_radius             : cell virtual radii (ABM)
            %       - cell_id                 : cell label
            %       - cell_total_distance     : cell total distance (trajectory)
            %       - cell_total_displacement : current position - initial_position

            monolayer.cell_radius             = ones(quantity,1)*inputs.geometry_radius;
            monolayer.cell_id                 = (1:quantity)';
            monolayer.cell_total_distance     = zeros(quantity,1);
            monolayer.cell_total_displacement = zeros(quantity,1);
            
        end



        %................... Shake cells using the magnitude given.........
        % This magnitude is multiplied by radius and a random number [-1, 1]
        function monolayer = shakeAllCells(monolayer)

            rng('default'); % Initializes generator based on the current 
                            % time, resulting in a different sequence of 
                            % random numbers after each call to rng. 
                            % choose default for a fixed random
                            % distribution. Modes: shuffle default
            magnitude = inputs.monolayer_shake_magnitude;
            for i = 1:length(monolayer.cell_id)
                new_position    = monolayer.cell_position(i,:);
                radius          = monolayer.cell_radius(i);
                for j = 1:length(monolayer.monolayer_size)-1
                    % rand() generates a number between [0,1). To generate
                    % a uniformly distributed random number in the interval
                    % between (-1,1): 2*rand()-1.
                    new_position(j) = new_position(j) + radius * magnitude * (2*rand()-1);
                    if new_position(j)<0
                        new_position(j)=0;
                    elseif new_position(j)> monolayer.monolayer_size(j)
                        new_position(j)= monolayer.monolayer_size(j);

                    end

                monolayer.cell_position(i,:)          = new_position;
                monolayer.cell_init_position(i,:)     = new_position;
                monolayer.cell_previous_position(i,:) = new_position;

                if strcmp(inputs.fem_BC_type,'PBC') || strcmp(inputs.fem_BC_type,'material_around')
                    % Boundary conditions: PBC or material around
                    % Cells can move to one side to another side, we need a
                    % copy of the position
                    monolayer.cell_PBC_position_PBC(i,:) = new_position;
                end

                end
                 
            end

            % Update list of monolayer_points
            for i = 1:length(monolayer.cell_id)
                monolayer.monolayer_cell_points(i,1) = monolayer.cell_position(i,1);
                monolayer.monolayer_cell_points(i,2) = monolayer.cell_position(i,2);               
            end

        end




        %................... Infect cells .................................
        function monolayer = infectCells(monolayer)

            % Create an area (polygon) where there is infection
            x_center        = inputs.monolayer_size_x/2;
            y_center        = inputs.monolayer_size_y/2;
            cir_radius      = inputs.monolayer_inection_radius;
            num_points      = 100;

            theta           = (0:num_points-1)*(2*pi/num_points);
            x_poly          = x_center + cir_radius*cos(theta);
            y_poly          = y_center + cir_radius*sin(theta);
            infected_circle = polyshape(x_poly,y_poly);
         
            x_cell_points   = monolayer.monolayer_cell_points(:,1);
            y_cell_points   = monolayer.monolayer_cell_points(:,2);

            inf_cells       = isinterior(infected_circle,x_cell_points,y_cell_points);

            % Update infection list
            monolayer.cell_infection = inf_cells;

%             list_infected_cells_x   = [];
%             list_infected_cells_y   = [];
%             list_uninfected_cells_x = [];
%             list_uninfected_cells_y = [];
% 
%             for i = 1:length(monolayer.cell_id)
% 
%                 if inf_cells(i) == 1
%                     list_infected_cells_x   = [list_infected_cells_x x_cell_points(i)];
%                     list_infected_cells_y   = [list_infected_cells_y y_cell_points(i)];
%                 else
%                     list_uninfected_cells_x = [list_uninfected_cells_x x_cell_points(i)];
%                     list_uninfected_cells_y = [list_uninfected_cells_y y_cell_points(i)];
%                 end
%             end
%            
%             plot(infected_circle)
%             axis equal
%             hold on
%             plot(list_infected_cells_x,list_infected_cells_y,'r*')
%             plot(list_uninfected_cells_x,list_uninfected_cells_y,'g*')
%             hold off

        
        end

    end




end