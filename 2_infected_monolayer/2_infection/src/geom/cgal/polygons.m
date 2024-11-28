classdef polygons

    properties %(Access = public)
        polygon
    end


    methods


        %................... Object initialization.........................
        function obj = polygons()
            addpath(fullfile('inputs')); 
            addpath(fullfile('src'));
        end


        
        
        %................. Generate polygons ..............................
        function [monolayer] = generatePolygons(~,monolayer)

            % Load cell vertices and regions
            vertices            = monolayer.voronoi_vertices;
            regions             = monolayer.voronoi_region  ;

            % Load alpha shape and create a polygon or load boundaries
            % (rectangle)
            %alphaVertices       = monolayer.alpha_shape; % Not used anymore
            alphaVertices       = [0 0; inputs.monolayer_size_x 0;                ...
                                  inputs.monolayer_size_x inputs.monolayer_size_y;...
                                  0 inputs.monolayer_size_y];

%             alpha_shape_polygon = polyshape(alphaVertices,'KeepCollinearPoints',true);
%             figure()
%             plot(alphaVertices(:,1),alphaVertices(:,2),'-k');
%             figure()
%             plot(alpha_shape_polygon)

            % Initialize 
            numRegions           = size(regions, 1)           ;  % Number of polygons
            polygonList          =  cell(1, numRegions)       ;  % List of polygons
            locationList         = zeros(1, numRegions)       ;  % Label border cell = 1, interior cell = 0
            elipsoid_centroid_x  = zeros(1, numRegions)       ;  % Cell ellipsoid, centroid coordinates 
            elipsoid_centroid_y  = zeros(1, numRegions)       ;
            elipsoid_centroid_z  = zeros(1, numRegions)       ;
            elipsoid_orientation = zeros(1, numRegions)       ;  % Cell ellipsoid, axis orientation 
            elipsoid_major_axis  = zeros(1, numRegions)       ;  % Cell ellipsoid, major axis length
            elipsoid_minor_axis  = zeros(1, numRegions)       ;  % Cell ellipsoid, minor axis length

            boundary_x           = alphaVertices(:,1)         ;
            boundary_y           = alphaVertices(:,2)         ;

            tol                  = inputs.general_tolerance   ;
            monolayer_size_x     = monolayer.monolayer_size(1);
            monolayer_size_y     = monolayer.monolayer_size(2);

%             figure()
%             hold on

            % THIS LOOP IS OVER THE NUMBER OF REGIONS, NOT INCLUDING
            % POSSIBLE EXTRUDING CELLS (IT COMES FROM VORONOI). CORRECT ID
            % NUMBERS AFTER THE LOOP
            for i = 1:numRegions                
                
                % Create a polygon for each cell
                coord_x_vert         = vertices(regions{i},1)            ;
                coord_y_vert         = vertices(regions{i},2)            ;
                my_polygon           = polyshape(...
                    coord_x_vert,coord_y_vert,'KeepCollinearPoints',true);

                polygonList(i)       = {my_polygon}                      ;
%                 plot(polygonList{i})

                % Small numbers can lead to misclassifications (for
                % instance, 1.77e-15), make sure there is no such problem
                % at the border.
                
                coord_x_vert(abs(coord_x_vert) < tol)                    = 0;
                coord_y_vert(abs(coord_y_vert) < tol)                    = 0               ;
                coord_x_vert(abs(coord_x_vert - monolayer_size_x) < tol) = monolayer_size_x;
                coord_y_vert(abs(coord_y_vert - monolayer_size_y) < tol) = monolayer_size_y;

                % Get whether there are vertices on the edge of the alpa
                % shape
                [~,on]               = inpolygon(coord_x_vert,coord_y_vert, ... 
                                                   boundary_x,boundary_y        );

                % If there is any vertex on the edge, polygon = exterior 
                locationList(i)      = any(on)                                   ;

%                 if  locationList(i) == 1  % exterior
%                     color = 'r';
%                 else
%                     color = 'b';
%                 end
%                 plot(polygonList{i},'FaceColor',color)



                % Compute major and minor axis associated to the cell
                % polygon

                %plot(my_polygon)
                % Convert polygon to binarize image
                BW    = poly2mask(coord_x_vert,coord_y_vert,monolayer_size_x,monolayer_size_y);
                %imagesc(BW)

                % Compute the properties of the binary image
                props = regionprops(BW,'Centroid','Orientation','MajorAxisLength','MinorAxisLength')                ;

                % Get the centroid, orientation and axis lengths
                cx    = props.Centroid(1)          ;
                cy    = props.Centroid(2)          ;
                theta = deg2rad(-props.Orientation);
                L1    = props.MajorAxisLength      ;
                L2    = props.MinorAxisLength      ;

                % Plot them on Matlab
                % Calculate the endpoints of the major axis
                %x1    = cx + L1/2 * cos(theta);
                %y1    = cy + L1/2 * sin(theta);
                %x2    = cx - L1/2 * cos(theta);
                %y2    = cy - L1/2 * sin(theta);

                % Draw the major axis
                %%line([x1 x2],[y1 y2],'Color','r','LineWidth',2)
                
                % Calculate the endpoints of the minor axis
                %x3    = cx + L2/2 * cos(theta+90*pi/180);
                %y3    = cy + L2/2 * sin(theta+90*pi/180);
                %x4    = cx - L2/2 * cos(theta+90*pi/180);
                %y4    = cy - L2/2 * sin(theta+90*pi/180);

                % Draw the minor axis
                %line([x3 x4],[y3 y4],'Color','g','LineWidth',2)
                
                % Use the parametric equation of an ellipse to generate points
                
                %t     = linspace(0,2*pi,50);
                %a     = L1/2;
                %b     = L2/2;
                %x     = cx + a*cos(t)*cos(theta) - b*sin(t)*sin(theta);
                %y     = cy + a*cos(t)*sin(theta) + b*sin(t)*cos(theta);
                %plot(x,y,'r','Linewidth',2)

                % Save variables
                elipsoid_centroid_x(i)  = cx;
                elipsoid_centroid_y(i)  = cy;
                elipsoid_orientation(i) = theta;
                elipsoid_major_axis(i)  = L1;
                elipsoid_minor_axis(i)  = L2;
                
                


            end
%             hold off
            
            % Clean structures
            monolayer.cellPolyList                = [];
            monolayer.cellLocationList            = [];
            monolayer.centroid_of_elipsoid_step_i = [];
            monolayer.orientation_axis_elipsoid_i = [];
            monolayer.major_axis_length_i         = [];
            monolayer.minor_axis_length_i         = [];

            % Update properties
            monolayer.cellPolyList                = polygonList;
            monolayer.cellLocationList            = locationList;

            monolayer.centroid_of_elipsoid_step_i = [elipsoid_centroid_x ; elipsoid_centroid_y ; elipsoid_centroid_z]';
            monolayer.orientation_axis_elipsoid_i = elipsoid_orientation';
            monolayer.major_axis_length_i         = elipsoid_major_axis';
            monolayer.minor_axis_length_i         = elipsoid_minor_axis';

            % CORRECT THE FACT THAT THERE ARE EXTRUDING CELLS THAT HAVE NOT
            % BEEN TAKEN INTO ACCOUNT IN THE VORONOI CALCULATION
            id_extruding_cells = monolayer.index_extruding_cells;
            if any(id_extruding_cells)

                num_extruding_cells = numel(id_extruding_cells);

                for k = 1:num_extruding_cells

                    pos                                   = id_extruding_cells(k);

                    monolayer.cellPolyList                = [               monolayer.cellPolyList(1:pos-1  )   0                  monolayer.cellPolyList(pos:end  ) ];
                    monolayer.cellLocationList            = [           monolayer.cellLocationList(1:pos-1  )   1              monolayer.cellLocationList(pos:end  ) ];
                    centroid_1                            = [monolayer.centroid_of_elipsoid_step_i(1:pos-1,1) ; 0 ; monolayer.centroid_of_elipsoid_step_i(pos:end,1) ];
                    centroid_2                            = [monolayer.centroid_of_elipsoid_step_i(1:pos-1,2) ; 0 ; monolayer.centroid_of_elipsoid_step_i(pos:end,2) ];
                    centroid_3                            = [monolayer.centroid_of_elipsoid_step_i(1:pos-1,3) ; 0 ; monolayer.centroid_of_elipsoid_step_i(pos:end,3) ];
                    monolayer.centroid_of_elipsoid_step_i = [centroid_1 centroid_2 centroid_3];
                    monolayer.orientation_axis_elipsoid_i = [monolayer.orientation_axis_elipsoid_i(1:pos-1  ) ; 0 ;  monolayer.orientation_axis_elipsoid_i(pos:end  )];
                    monolayer.major_axis_length_i         = [        monolayer.major_axis_length_i(1:pos-1  ) ; 0 ;          monolayer.major_axis_length_i(pos:end  )];
                    monolayer.minor_axis_length_i         = [        monolayer.minor_axis_length_i(1:pos-1  ) ; 0 ;          monolayer.minor_axis_length_i(pos:end  )];

                end

            end



            % ........... Identify uninfected surrounding cells ...........
            % .............................................................

            % ...... (1) First line of neighbouring uninfected cells ......
            % Identify first line of surrounding cells close to the
            % infection border
            
            % Initialize surrounder cell list, even if there is no infection,
            % so one can plot it later on paraview anyways
            numb_cells                = length(monolayer.cell_id) ;
            cell_surrounder_list      = zeros(1, numb_cells      );
            monolayer.surrounder_cell = cell_surrounder_list      ;
            infection_list            = monolayer.cell_infection  ;
            extruding_cells           = monolayer.cell_extrusion  ;

            % get right polygonlist in case there was extrusion
            polygonList               = monolayer.cellPolyList    ;



            if inputs.monolayer_infection_active == true

                % This part is computed in two different scenarios
                %     (1) When the infection is fixed (infection over all
                %     steps)
                %     (2) When we simulate propagation, but first we let
                %     the program run a certain number of steps to avoid
                %     high forces during the first steps

                if inputs.monolayer_infection_propagation == false || ...
                        (inputs.monolayer_infection_propagation == true && monolayer.step >= inputs.monolayer_infection_step_propagation && ...
                        any(infection_list)~=false)
    
                    % Get infected polygon ......................               
                    infected_polygon = [];
                    counter          = 1 ;
        
                    % Go over every infected cell and get its polyshape union with
                    % the previous infected cells
                    for i = 1:numb_cells
        
                        if infection_list(i) == 1 && extruding_cells(i) == 0 % Infected and not extruding cell
                            if counter == 1
                                infected_polygon = polygonList{i};
                            else
                                infected_polygon = union(infected_polygon,polygonList{i});
                            end
                            counter                    = counter + 1;
                        end
                                      
                    end
%                     plot(infected_polygon)
    
                    % Update infected polygons
                    monolayer.infected_polygon = []              ;
                    monolayer.infected_polygon = infected_polygon;
    
    
    
    
                    % Retrieve and tag cells that intersect with the infected
                    % polygon
    
                    inf_focus_vertices   = infected_polygon.Vertices;
    
    %                 figure();
    %                 hold on
    
                    for i = 1:numb_cells
    
                        if infection_list(i) == 0 && extruding_cells(i) == 0 % Uninfected and not extruding cell
    
                            uninfected_polygon = polygonList{i}             ;
                            uninf_pol_vertices = uninfected_polygon.Vertices;
    
                            % Check whether there are vertices on the edge of
                            % the infection
                            [~,on]      = inpolygon(uninf_pol_vertices(:,1),uninf_pol_vertices(:,2), ... 
                                                    inf_focus_vertices(:,1),inf_focus_vertices(:,2));
                            % If there is any vertex on the edge, polygon = exterior 
                            cell_surrounder_list(i) = any(on);
                            
    %                         if  cell_surrounder_list(i) == 1  % exterior
    %                             color = 'r';
    %                         else
    %                             color = 'b';
    %                         end
    %                         plot(polygonList{i},'FaceColor',color)
    
                        end
    
                    end
    %                 hold off     
    
                    % Update list
                    monolayer.surrounder_cell = [];
                    monolayer.surrounder_cell = cell_surrounder_list;

                end

            end


            % ..... (2) Second line of neighbouring uninfected cells ......
            % Identify second line of surrounding cells close to the
            % infection border
            
            % Initialize surrounder cell list, even if there is no infection,
            % so one can plot it later on paraview anyways
            cell_surrounder_list_2nd      = zeros(1, numb_cells      )    ;
            monolayer.surrounder_cell_2nd = cell_surrounder_list_2nd      ;



            if inputs.monolayer_infection_active == true

                % This part is computed in two different scenarios
                %     (1) When the infection is fixed (infection over all
                %     steps)
                %     (2) When we simulate propagation, but first we let
                %     the program run a certain number of steps to avoid
                %     high forces during the first steps

                if inputs.monolayer_infection_propagation == false || ...
                        (inputs.monolayer_infection_propagation == true && monolayer.step >= inputs.monolayer_infection_step_propagation && ...
                        any(infection_list)~=false)
    
                    % Get infected polygon + first line of uninfected cells               
                    infected_polygon_1stline = [];
                    counter                  = 1 ;
        


                    % Go over every uninfected neigbouring cell (1st line) 
                    % and get its polyshape union with the previous 
                    % infected cell polygon
                    for i = 1:numb_cells
        
                        if cell_surrounder_list(i) == 1 && extruding_cells(i) == 0 % Infected and not extruding cell
                            if counter == 1
                                infected_polygon_1stline = union(infected_polygon        ,polygonList{i}); 
                            else
                                infected_polygon_1stline = union(infected_polygon_1stline,polygonList{i}); 
                            end
                            counter = counter + 1;
                            
                        end
                                      
                    end
%                     plot(infected_polygon_1stline)
    
                    % Update infected cell + 1st line of neighbouring
                    % uninfected cell polygon
                    monolayer.infected_polygon_1st = []                      ;
                    monolayer.infected_polygon_1st = infected_polygon_1stline;
    
    
    
    
                    % Retrieve and tag cells that intersect with the 
                    % polygon of infected cells + 1st line of neigh. uninf.
                    % cells
    
                    polygon_vertices   = infected_polygon_1stline.Vertices;
    
%                     figure();
%                     hold on
    
                    for i = 1:numb_cells
    
                        if infection_list(i) == 0 && extruding_cells(i) == 0 % Uninfected and not extruding cell

                            if cell_surrounder_list(i) == 0
    
                                uninfected_polygon = polygonList{i}             ;
                                uninf_pol_vertices = uninfected_polygon.Vertices;
        
                                % Check whether there are vertices on the edge of
                                % the infection
                                [~,on]      = inpolygon(uninf_pol_vertices(:,1),uninf_pol_vertices(:,2), ... 
                                                          polygon_vertices(:,1), polygon_vertices(:,2));
                                % If there is any vertex on the edge, polygon = exterior 
                                cell_surrounder_list_2nd(i) = any(on);
                                
%                                 if  cell_surrounder_list(i) == 1  % exterior
%                                     color = 'r';
%                                 else
%                                     color = 'b';
%                                 end
%                                 plot(polygonList{i},'FaceColor',color)

                            end
    
                        end
    
                    end
%                     hold off     
    
                    % Update list
                    monolayer.surrounder_cell_2nd = [];
                    monolayer.surrounder_cell_2nd = cell_surrounder_list_2nd;

                end

            end





            % ..... (3) Third line of neighbouring uninfected cells .......
            % Identify third line of surrounding cells close to the
            % infection border
            
            % Initialize surrounder cell list, even if there is no infection,
            % so one can plot it later on paraview anyways
            cell_surrounder_list_3rd      = zeros(1, numb_cells      )    ;
            monolayer.surrounder_cell_3rd = cell_surrounder_list_3rd      ;



            if inputs.monolayer_infection_active == true

                % This part is computed in two different scenarios
                %     (1) When the infection is fixed (infection over all
                %     steps)
                %     (2) When we simulate propagation, but first we let
                %     the program run a certain number of steps to avoid
                %     high forces during the first steps

                if inputs.monolayer_infection_propagation == false || ...
                        (inputs.monolayer_infection_propagation == true && monolayer.step >= inputs.monolayer_infection_step_propagation && ...
                        any(infection_list)~=false)
    
                    % Get infected polygon + first and second line of 
                    % uninfected cells               
                    infected_polygon_2ndline = [];
                    counter                  = 1 ;
        


                    % Go over every uninfected neigbouring cell (2nd line) 
                    % and get its polyshape union with the previous 
                    % infected cell polygon + 1st neighbouring line
                    for i = 1:numb_cells
        
                        if cell_surrounder_list_2nd(i) == 1 && extruding_cells(i) == 0 % surr. cell and not extruding cell
                            if counter == 1
                                infected_polygon_2ndline = union(infected_polygon_1stline, polygonList{i});
                            else
                                infected_polygon_2ndline = union(infected_polygon_2ndline, polygonList{i}); 
                            end
                            counter = counter + 1;
                            
                        end
                                      
                    end
%                     plot(infected_polygon_2ndline)
    
                    % Update infected cell + 1st and 2nd line of 
                    % neighbouring uninfected cell polygon
                    monolayer.infected_polygon_2nd = []                      ;
                    monolayer.infected_polygon_2nd = infected_polygon_2ndline;
    
    
    
    
                    % Retrieve and tag cells that intersect with the 
                    % polygon of infected cells + 1st and 2nd line of neigh
                    % uninfected cells
                    polygon_vertices   = infected_polygon_2ndline.Vertices;
    
%                     figure();
%                     hold on
    
                    for i = 1:numb_cells
    
                        if infection_list(i) == 0 && extruding_cells(i) == 0 % Infected and not extruding cell

                            if cell_surrounder_list_2nd(i) == 0 && cell_surrounder_list(i) == 0
    
                                uninfected_polygon = polygonList{i}             ;
                                uninf_pol_vertices = uninfected_polygon.Vertices;
        
                                % Check whether there are vertices on the edge of
                                % the infection
                                [~,on]      = inpolygon(uninf_pol_vertices(:,1),uninf_pol_vertices(:,2), ... 
                                                          polygon_vertices(:,1), polygon_vertices(:,2));
                                % If there is any vertex on the edge, polygon = exterior 
                                cell_surrounder_list_3rd(i) = any(on);
                                
%                                 if  cell_surrounder_list_3rd(i) == 1  % exterior
%                                     color = 'r';
%                                 else
%                                     color = 'b';
%                                 end
%                                 plot(polygonList{i},'FaceColor',color)

                            end
    
                        end
    
                    end
%                     hold off     
    
                    % Update list
                    monolayer.surrounder_cell_3rd = [];
                    monolayer.surrounder_cell_3rd = cell_surrounder_list_3rd;

                end






            % ........... Identify exterior infected cells ................
            % .............................................................

            % ...... (1) First line of exterior INnfected cells ...........
            % Identify first line of exterior infected cells in the
            % infection
            
            % Initialize surrounder cell list, even if there is no infection,
            % so one can plot it later on paraview anyways
            numb_cells                = length(monolayer.cell_id) ;
            inf_exterior_cell         = zeros(1, numb_cells      );
            monolayer.exterior_infected_cells = inf_exterior_cell ;
            infection_list            = monolayer.cell_infection  ;
            extruding_cells           = monolayer.cell_extrusion  ;

            % get right polygonlist in case there was extrusion
            polygonList               = monolayer.cellPolyList    ;



            if inputs.monolayer_infection_active == true

                % This part is computed in two different scenarios
                %     (1) When the infection is fixed (infection over all
                %     steps)
                %     (2) When we simulate propagation, but first we let
                %     the program run a certain number of steps to avoid
                %     high forces during the first steps

                if inputs.monolayer_infection_propagation == false || ...
                        (inputs.monolayer_infection_propagation == true && monolayer.step >= inputs.monolayer_infection_step_propagation && ...
                        any(infection_list)~=false)
    
                    % Get infected polygon ......................    
                    % Previously computed
                    infected_polygon = monolayer.infected_polygon;

                    % Retrieve and tag cells that intersect with the infected
                    % polygon
                    inf_focus_vertices   = infected_polygon.Vertices;
    
                    % figure();
                    % hold on
    
                    for i = 1:numb_cells
    
                        if infection_list(i) == 1 && extruding_cells(i) == 0 % Infected and not extruding cell
    
                            infected_polygon = polygonList{i}             ;
                            inf_pol_vertices = infected_polygon.Vertices  ;
    
                            % Check whether there are vertices on the edge of
                            % the infection
                            [~,on]      = inpolygon(inf_pol_vertices(:,1),inf_pol_vertices(:,2), ... 
                                                    inf_focus_vertices(:,1),inf_focus_vertices(:,2));
                            % If there is any vertex on the edge, polygon = exterior 
                            inf_exterior_cell(i) = any(on);
                            
                            % if  inf_exterior_cell(i) == 1  % exterior
                            %     color = 'r';
                            % else
                            %     color = 'b';
                            % end
                            % plot(polygonList{i},'FaceColor',color)
    
                        end
    
                    end
                    % hold off     
    
                    % Update list
                    monolayer.exterior_infected_cells = [];
                    monolayer.exterior_infected_cells = inf_exterior_cell;

                end

            end

















                % ........... Compute infection focus area ................
                % .........................................................
                infection_list            = monolayer.cell_infection  ;
                extruding_cells           = monolayer.cell_extrusion  ;
                list_of_polygons          = monolayer.cellPolyList    ;

                % Extruding cells are only infected cells. If they have
                % been already extruded, they are not considered in the
                % calculation of the area of infection. So, we need to
                % remove them from the infection list here:

                combined_list = infection_list - extruding_cells;
                inf_indices   = find(combined_list == 1)        ;
                inf_area      = 0                               ;

                if ~isempty(inf_indices)
                
                    for i = 1 : length(inf_indices) 
    
                        index        = inf_indices(i)         ;
                        inf_polygon  = list_of_polygons{index};
                        area_inf_pol = area(inf_polygon)      ;
                        inf_area     = inf_area + area_inf_pol;
    
    
                    end

                    monolayer.total_infection_area = inf_area;
                    disp("Infection area")
                    disp(inf_area)

                    step_here = monolayer.step;

                    fid = fopen('./results/area_inf.txt', 'a'); % write infection area and append
                    if fid == -1
                        error('cannot open the file area_inf.txt');
                    end

                    % Write values   
                    fprintf(fid, '%f %f\n', step_here, inf_area);
                    fclose(fid); 

                    
                else

                    monolayer.total_infection_area = 0;
                    disp("Infection area")
                    disp("0")

                    step_here = monolayer.step;
                    inf_area = 0;

                    fid = fopen('./results/area_inf.txt', 'a'); % write infection area
                    if fid == -1
                        error('cannot open the file area_inf.txt');
                    end

                    % Write values   
                    fprintf(fid, '%f %f\n', step_here, inf_area);
                    fclose(fid); 

                end



            end










       
        end







        %................... Propagate infection ..........................
        function [list_inf_cells] = propagateInfection(~,mon_cell_points, mins_post_inf)

            % Create an area (polygon) where there is infection
            x_center        = inputs.monolayer_size_x/2         ;
            y_center        = inputs.monolayer_size_y/2         ;
            %cir_radius      = inputs.monolayer_inection_radius  ;
            num_points      = 100                               ;
            

            % Assume the area of infection increases: area = a*time+ b
            a = inputs.monolayer_infection_propg_constant       ;
            b = -6000                                           ;
            % make sure there is at least one infected cell when the
            % infection starts (choose a right value for b)

            % Compute area of infection, assuming the focus is a circle
            area_infection  = a * (mins_post_inf+310) + b       ;
            cir_radius      = (area_infection/pi)^0.5           ;

            theta           = (0:num_points-1)*(2*pi/num_points);
            x_poly          = x_center + cir_radius*cos(theta)  ;
            y_poly          = y_center + cir_radius*sin(theta)  ;
            infected_circle = polyshape(x_poly,y_poly)          ;
         
            x_cell_points   = mon_cell_points(:,1)              ;
            y_cell_points   = mon_cell_points(:,2)              ;

            inf_cells       = isinterior(infected_circle,x_cell_points,y_cell_points);

            % Update infection list
            list_inf_cells  = inf_cells                         ;

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