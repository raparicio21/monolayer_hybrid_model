classdef create_triang

    methods


        %................... Object initialization.........................
        function obj = create_triang()
            addpath(fullfile('inputs')); 
            addpath(fullfile('src'));
        end


        
        
        %.......... Generate Triangulation ..............................
        function [TR, infection_label, id_cell_element, triangle_areas, save_id_removed_vertices] = createTriangulation(~,monolayer,ecm)

            % Register cells (cell centroids) and polygon vertices
            list_centroids = monolayer.monolayer_cell_points;
            list_vertices  = monolayer.voronoi_vertices;
            list_regions   = monolayer.voronoi_region;

            % If there are extruding cells, list_centroids includes them.
            % Update it to remove them here
            list_extruding_cells = monolayer.cell_extrusion;
            if(any(list_extruding_cells))
                indices                   = find(list_extruding_cells == 1);
                list_centroids(indices,:) = []                             ;
            end

            % Remove border centroids that are on the vertex list
%             indices = ismember(list_vertices,list_centroids,'rows');
%             list_vertices(indices,:) = [];
            

            % Limit position with tolerance
            var1       = ecm.ecm_size(1) - list_centroids(:,1);
            positions  = find(  var1 < inputs.general_tolerance);
            list_centroids(positions,1) = ecm.ecm_size(1);

            var2       = ecm.ecm_size(2) - list_centroids(:,2);
            positions  = find(  var2 < inputs.general_tolerance);
            list_centroids(positions,2) = ecm.ecm_size(2);

            var3       = list_centroids(:,1);
            positions  = find(  var3 < inputs.general_tolerance);
            list_centroids(positions,1) = 0;

            var4       = list_centroids(:,2);
            positions  = find(  var4 < inputs.general_tolerance);
            list_centroids(positions,2) = 0;

            % Create list        
            global_list              = [list_centroids; list_vertices];
            connect_list             = []                             ;
            infection_label          = []                             ;
            id_cell_element          = []                             ;
            save_id_removed_vertices = []                             ;
            numb_cells               = numel(monolayer.cell_id)       ;
            %mon_size_x      = inputs.monolayer_size_x;
            %mon_size_y      = inputs.monolayer_size_y;

            % In case of extruding cells, numb_cells is not updated. Update
            % it here. We adapt all the structures to work only with
            % non-extruding cells. Regarding that, we adapt the list of
            % infected cells that is used in the following for loop. 
            new_list_inf_cells       = monolayer.cell_infection       ;
            if(any(list_extruding_cells))

                numb_extruding_cells          = sum(list_extruding_cells)        ;
                numb_cells                    = numb_cells - numb_extruding_cells;

                
                indices                       = find(list_extruding_cells == 1)  ;
                new_list_inf_cells(indices,:) = []                               ;

            end

            for i = 1  : numb_cells

                  % Load parameters
                  cent                          = list_centroids(i,:)     ;
                  vert_id                       = list_regions{i}         ; 
                  vert_coord                    = list_vertices(vert_id,:);
                
                  % Remove (locally) vertices that are the centroid
                  isCloseVertices               = abs     (vert_coord(:,1) - cent(1)) < 1e-3 & abs(vert_coord(:,2) - cent(2)) < 1e-3;
                  vertices_to_remove            = vert_id (isCloseVertices,:)                                                       ;
                  save_id_removed_vertices      = [save_id_removed_vertices ; vertices_to_remove]                                   ;
                  vert_coord(isCloseVertices,:) = [];
                  vert_id   (isCloseVertices,:) = [];
                
                  % Create connectivity
                  connect                       = [];
                  
                  for j = 1 : length(vert_id)
                                
                          if j == length(vert_id)
                              local_triangle = [i   vert_id(j)+numb_cells    vert_id(1)+numb_cells  ];
                              % There are some cases where this triangle
                              % can be wrong (negative area, cells at the
                              % border for example, since we cut the cell,
                              % the centroid might be on the corner and the
                              % triangle should not exist, Get those cases
                              % so we do not consider the triangle
                              centroid           = cent                                             ;
                              vertex_A           = vert_coord(j,:)                                  ;
                              vertex_B           = vert_coord(1,:)                                  ;

                              vector1            = vertex_A - centroid                              ;
                              vector2            = vertex_B - centroid                              ;
                              cross_product      = vector1(1) * vector2(2) - vector1(2) * vector2(1);
                              isCounterClockwise = (cross_product > 0)                              ;
                          else  
                              local_triangle = [i   vert_id(j)+numb_cells    vert_id(j+1)+numb_cells];

                              % Same problem as the other case on the "if"
                              % loop
                              centroid           = cent                                             ;
                              vertex_A           = vert_coord(j,:)                                  ;
                              vertex_B           = vert_coord(j+1,:)                                ;

                              vector1            = vertex_A - centroid                              ;
                              vector2            = vertex_B - centroid                              ;
                              cross_product      = vector1(1) * vector2(2) - vector1(2) * vector2(1);
                              isCounterClockwise = (cross_product > 0)                              ;                         
                          end

                          % Only if the area is positive (counter clockwise
                          % triangle) we include the triangle
                          if isCounterClockwise == 1
                    
                              connect            = [connect; local_triangle]             ;                  
                        
                              % Label element depending on infection. 
                              % 0 = uninfected
                              % 1 = infected                  
                              infection_id_label = new_list_inf_cells(i)                 ;
                              infection_label    = [infection_label; infection_id_label] ;
    
                              % Label element depending on the cell it belongs
                              % to
                              id_cell_element    = [id_cell_element; i]                  ;
                          end
                  end         

                  connect_list = [connect_list; connect]                   ;      

            end
              

            % Create a background mesh (NOT IMPLEMENTED)


            P                   = global_list                    ;
            T                   = connect_list                   ;

        
            % Solve the problem there are centroids that are on the list of
            % vertices and update region matrix accordingly
            vertices_unused     = setdiff(1:size(P,1), unique(T));
            num_vertices_unused = numel(vertices_unused)         ;

            if num_vertices_unused ~= 0
                for i=1:size(T,1)
                    for j=1:3
                        id = T(i,j);
                        if id < min(vertices_unused)
                            continue
                        else
                            for k = num_vertices_unused:-1:1
                                if id > vertices_unused(k)
                                    T(i,j) = id - k;
                                    break
                                end
                            end
                        end
                    end
                end
                P(vertices_unused,:) = [];
            end

            % Check whether there are negative areas in the triangulation
            % (for ex. exception centroids on the border)
            triangle_id    = [];
            triangle_areas = zeros(size(T,1),1);
            
            for i = 1:size(T,1)
                v1 = P(T(i,1),:);
                v2 = P(T(i,2),:);
                v3 = P(T(i,3),:);
                area = polyarea([v1(1), v2(1), v3(1)], [v1(2), v2(2), v3(2)]);
                if area <= inputs.general_tolerance
                    triangle_id = [triangle_id i];
                end
                triangle_areas(i) = area;
            end

            % Remove triangle and id infection
            T(triangle_id,:)             = [];
            infection_label(triangle_id) = [];
            id_cell_element(triangle_id) = [];
            triangle_areas(triangle_id)  = [];
            

            %warning off
            TR = triangulation(T,P);
            %warning on

            % Plot triangulation
%             triplot(TR)
          
                         
        end
        
        



    end

end