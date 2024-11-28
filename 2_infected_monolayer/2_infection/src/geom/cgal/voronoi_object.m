classdef voronoi_object

    methods


        %................... Object initialization.........................
        function obj = voronoi_object()
            addpath(fullfile('inputs')); 
            addpath(fullfile('src'));
        end


        
        
        %.......... Generate Voronoi diagram ..............................
        function monolayer = generateVoronoi(~,monolayer)

            % ........... Load x, y coordinates (cell centroids)...........
            x       = monolayer.monolayer_cell_points(:,1);
            y       = monolayer.monolayer_cell_points(:,2);

            % Remove extruding cell points
            list_extruding_cells = monolayer.cell_extrusion;
            if any(list_extruding_cells)
                indices      = find(list_extruding_cells == 1);
                x(indices)   = [];
                y(indices)   = [];
                % Save indices
                monolayer.index_extruding_cells = indices;

            end

            % ........... Compute Voronoi vertices (plot) .................
            % ........... without boundary intersection
%             [VX,VY] = voronoi(x,y);
%             plot(VX,VY,'-b',x,y,'.r');
%             xlim([-20, inputs.substrate_size_x + 20])
%             ylim([-20, 220                         ])
%             % Assign labels to the points X
%             nump    = size(x);
%             plabels = arrayfun(@(n) {sprintf('X%d', n)}, (1:nump)');
%             hold on
%             text(x, y+0.2, plabels, 'color', 'r', ...
%             'FontWeight', 'bold', 'HorizontalAlignment',...
%             'center', 'BackgroundColor', 'none');
%             hold off

            % ........... Compute the Voronoi diagram .....................
            % ........... without boundary intersection
%             dt                   = delaunayTriangulation(x,y);
%             [V,R]                = voronoiDiagram(dt);
            % V : is a matrix representing the coordinates of the Voronoi
            %     vertices (the vertices are the end pointsof the
            %     Voronoiedges). By convection the first vertex in V is the
            %     infinite vertex.
            % R : is a vector cell array representing the Voronoi region
            %     associated with each point. Hence, the Voronoi region
            %     associated with the point [x(i),y(i)] is R{i}

            % Assign labels to the Voronoi vertices V (plot)
            % By convention the first vertex is at infinity
%             numv = size(V,1);
%             vlabels = arrayfun(@(n) {sprintf('V%d', n)}, (2:numv)');
%             hold on
%             Hpl = text(V(2:end,1), V(2:end,2)+.2, vlabels, ...
%                   'FontWeight', 'bold', 'HorizontalAlignment',...
%                   'center', 'BackgroundColor', 'none');
%             hold off



            
            % ........... Compute Voronoi diagram .........................
            %         (intersect with monolayer border) 
            
            % Define parameters to crop Voronoi           
            %bs_ext = monolayer.alpha_shape; % not used anymore
            % Exterior boundary: cell monolayer size
            bs_ext = [0 inputs.monolayer_size_x inputs.monolayer_size_x 0;...
                        inputs.monolayer_size_y inputs.monolayer_size_y 0 0]'; 
            X      = [x y];                 % Data points

            % VoronoiLimit function (used to crop voronoi diagram, some
            % points are in the infinite)
            [V_new,C_new,~] = VoronoiLimit(X(:,1),X(:,2),'bs_ext',bs_ext,'figure','off');
%             [V_new,C_new,XY_new] = VoronoiLimit(X(:,1),X(:,2),'bs_ext',bs_ext); % Plot
% 
% 
%             % Plot result
%             hold on
%             % Assign labels to the points.
%             nump    = size(XY_new,1);
%             plabels = arrayfun(@(n) {sprintf('X%d', n)}, (1:nump)');
%             Hpl     = text(XY_new(:,1), XY_new(:,2), plabels, 'FontWeight', ...
%                       'bold', 'HorizontalAlignment','center', ...
%                       'BackgroundColor', 'none');
%             
%             % Assign labels to the Voronoi vertices V.
%             % By convention the first vertex is at infinity.
%             numv    = size(V_new,1);
%             vlabels = arrayfun(@(n) {sprintf('V%d', n)}, (2:numv)');
%             hold on
%             Hpl     = text(V_new(2:end,1), V_new(2:end,2)+.2, vlabels, ...
%                      'FontWeight', 'bold', 'HorizontalAlignment',...
%                      'center', 'BackgroundColor', 'none');
%             hold off

            % Update voronoi intersected with alpha shapes
            monolayer.voronoi_vertices = [];
            monolayer.voronoi_region   = [];

            monolayer.voronoi_vertices = V_new;
            monolayer.voronoi_region   = C_new;

            % Clear vertices forces and update size
            monolayer.vertices_forces = zeros(length(V_new),3);      
            

%             [V_new,C_new,XY_new] = VoronoiLimit(X(:,1),X(:,2),'bs_ext',bs_ext); % Plot
        end
        



    end

end