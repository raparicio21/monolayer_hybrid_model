classdef cell_shapes


    methods

        %................... Object initialization.........................
        function obj = cell_shapes()
            addpath(fullfile('inputs')); 
            addpath(fullfile('src'));
        end

        
        %................... Cell  shape calculation.......................
        %................... Alpha shape calculation.......................
        % Generate a discretized circumference around each cell through
        % cell radius and compute alpha shape

        function monolayer = update_shape(~,monolayer)

            % Load parameters
            resolution        = inputs.geometry_shape_resolution;
            shape_mode        = inputs.geometry_shape_mode;
            %shape_orientation = inputs.geometry_shape_orientation;


            % Load cell center coordinates 
            x = monolayer.monolayer_cell_points(:,1);
            y = monolayer.monolayer_cell_points(:,2);

            if strcmp(shape_mode,'circle')

                angle_step = 2*pi/resolution;
                theta      = 0:angle_step:2*pi;
                r          = inputs.geometry_radius;
                ECM_x_size = inputs.substrate_size_x;
                ECM_y_size = inputs.substrate_size_y;

                % Case 1: consider all the discretized points = cell
                % centroids + discretized circumferences
                additional_points_x = [];
                additional_points_y = [];

                for i = 1:length(monolayer.cell_container)

                    % Load cell centroid coordinates
                    x0                  = x(i);
                    y0                  = y(i);

                    % Generate discretized circumference                    
                    x_new               = r * cos(theta) + x0;
                    y_new               = r * sin(theta) + y0;

                    % Remove points out of the cell monolayer border
                    % CONSIDER DOING THIS WITH FUNCTION INPOLYGON TO SPEED
                    % UP
                    index_x_down        = find( x_new < 0          );
                    index_x_up          = find( x_new > ECM_x_size );
                    index_x             = [index_x_down; index_x_up];
                    x_new(index_x)      = [];
                    y_new(index_x)      = [];
    
                    index_y_down        = find( y_new < 0          );
                    index_y_up          = find( y_new > ECM_y_size );
                    index_y             = [index_y_down; index_y_up];
                    x_new(index_y)      = [];
                    y_new(index_y)      = [];
    
                    additional_points_x = [additional_points_x  x_new];
                    additional_points_y = [additional_points_y  y_new];

                end

                
                % Case 2: case 1 but only considering border points (20%)
                % Set limits
%                 size_rect_x  = inputs.substrate_size_x;
%                 size_rect_y  = inputs.substrate_size_y;
%                 x_down_limit = size_rect_x*0.2;
%                 x_up_limit   = size_rect_x - x_down_limit;
%                 y_down_limit = size_rect_y*0.2;
%                 y_up_limit   = size_rect_y - y_down_limit;
%     
%                 additional_points_x = [];
%                 additional_points_y = [];
% 
%                 for i = 1:length(monolayer.cell_container)
% 
%                     % Load cell centroid coordinates
%                     x0                  = x(i);
%                     y0                  = y(i);
%                     
%                     % Create circles only at the points that meet the 
%                     % requirements
%                     if x0<x_down_limit || x0>x_up_limit || y0<y_down_limit || y0>y_up_limit
% 
%                         % Generate discretized circumference     
%                         x_new               = r * cos(theta) + x0;
%                         y_new               = r * sin(theta) + y0;
% 
%                         % Remove points out of the cell monolayer border    
%                         index_x_down        = find( x_new < 0  );
%                         index_x_up          = find( x_new > 200);
%                         index_x             = [index_x_down; index_x_up];
%                         x_new(index_x)      = [];
%                         y_new(index_x)      = [];
%     
%                         index_y_down        = find( y_new < 0  );
%                         index_y_up          = find( y_new > 200);
%                         index_y             = [index_y_down; index_y_up];
%                         x_new(index_y)      = [];
%                         y_new(index_y)      = [];
% 
%                         additional_points_x = [additional_points_x  x_new];
%                         additional_points_y = [additional_points_y  y_new];
%     
%                     end
%                 end

                % Add new points to cell centroids
                additional_points_all_x = [x ; additional_points_x'];
                additional_points_all_y = [y ; additional_points_y'];
    
                % Plot points
%                 figure()
%                 hold on
%                 plot(x(:),y(:),'.')
%                 plot(additional_points_x,additional_points_y,'.')
%                 hold off
    
                % Removed duplicated points
                additional_points_all = [additional_points_all_x additional_points_all_y];
                additional_points_all = unique(additional_points_all, 'rows');

                % Compute alpha shape
                shp       = alphaShape(additional_points_all);
                shp.Alpha = shp.Alpha + r;
%                 figure()
%                 plot(shp)

                % Extract the contour
                % p contains the points, bf the connectivity
                [~, p]            = shp.boundaryFacets();
%                 [bf, p]            = shp.boundaryFacets();
                
%                 figure()
%                 plot(p(:,1),p(:,2),'-k');

                % If you want to close the alpha shape for plotting:
%                 bf       = bf';
%                 xp       = p(:,1);
%                 yp       = p(:,2);
%                 xp       = xp(bf);
%                 yp       = yp(bf);
%                 xp       = xp(:);
%                 yp       = yp(:); 
%                 numedges = size(bf,2); % This is useful when there is more than one alpha shape
%                 xp       = [xp; NaN(1, numedges)];
%                 yp       = [yp; NaN(1, numedges)];
%                 figure()
%                 plot(xp,yp,'-k');


                % Update alpha shape
                monolayer.alpha_shape        = p;

                
            else
                error('Error: unsupported geometry shape mode!'); 
            end

        end





        
    end

end