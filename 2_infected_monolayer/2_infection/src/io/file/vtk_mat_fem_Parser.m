classdef vtk_mat_fem_Parser

    methods

        %................... Object initialization ........................
        function obj = vtk_mat_fem_Parser()
            addpath(fullfile('inputs')); 
            addpath(fullfile('src'));
        end




        %................... Write FEM data to VTK/mat file format ......
        function [obj,fem_results] = write_fem_VTK_mat_file(obj,fem_results,monolayer,mesh_obj)

            [fem_results] = obj.writeFEMData(fem_results,monolayer,mesh_obj);
            
        end






        function [fem_results] = writeFEMData(~,fem_results,monolayer,mesh_obj)

            step              = monolayer.step                                ;
                            
            nodes_x           = mesh_obj.triangulation_obj.Points(:,1)        ;             
            nodes_y           = mesh_obj.triangulation_obj.Points(:,2)        ;
            nodes_z           = zeros(size(nodes_x,1),1)                      ;

            mesh_connectivity = mesh_obj.triangulation_obj.ConnectivityList   ;

            if strcmp(inputs.fem_software,'pde_matlab')
    
                u_x     = fem_results.results_pde.Displacement.x          ;
                u_y     = fem_results.results_pde.Displacement.y          ;
                u_z     = zeros(size(u_x,1),1)                            ;
                u_disp  = [u_x, u_y, u_z]                                 ;
                %u_magn = fem_results.results_pde.Displacement.Magnitude;
    
                s_x     = fem_results.results_pde.Stress.xx               ;
                s_y     = fem_results.results_pde.Stress.yy               ;
                s_xy    = fem_results.results_pde.Stress.xy               ;
    
                pStress = evaluatePrincipalStress(fem_results.results_pde);
                s_I     = pStress.s1                                      ;
                s_II    = pStress.s2                                      ;
            
            elseif strcmp(inputs.fem_software,'abaqus')

                u_x         = fem_results.displ_results(:,2)            ;
                u_y         = fem_results.displ_results(:,3)            ;
                u_z         = zeros(size(u_x,1),1)                      ;
                u_disp      = [u_x, u_y, u_z]                           ;

                int_point_x = fem_results.stress_results(:,1)           ;
                int_point_y = fem_results.stress_results(:,2)           ;
                int_point_z = zeros(size(int_point_x,1),1)              ;
                int_point   = [int_point_x int_point_y int_point_z]     ;
                s_x         = fem_results.stress_results(:,3)           ;
                s_y         = fem_results.stress_results(:,4)           ;
                s_xy        = fem_results.stress_results(:,5)           ;
                s_I         = fem_results.stress_results(:,6)           ;
                s_II        = fem_results.stress_results(:,7)           ;
    

            end

            stress_I_projection_x  = zeros(length(s_I),1);
            stress_I_projection_y  = zeros(length(s_I),1);
            stress_II_projection_x = zeros(length(s_I),1);
            stress_II_projection_y = zeros(length(s_I),1);

            % Compute eigen vectors 
            for i = 1:length(s_I)

                stress_matrix           = [s_x(i)    s_xy(i)  ;  s_xy(i)  s_y(i)];
                [prin_dir, prin_values] = eig(stress_matrix,'vector')            ;

                % Sort eigenvectors
                if prin_values(1) > prin_values(2)
                    [prin_values, ind] = sort(prin_values);
                    prin_dir           = prin_dir(:, ind) ;
                end

%                 if prin_values(1)>0
%                     prin_values(1) = 0;
%                 end
% 
%                 if prin_values(2)<0
%                     prin_values(2) = 0;
%                 end
                            
                % Project principal components in order to plot principal 
                % stresses in X,Y coordinates. Make use of scalar product and 
                % vector projection in 2D   
                stress_I_projection_x(i)  = abs(prin_values(2)) * prin_dir(1,2);  % stress_I   vector projected in X axis
	            stress_I_projection_y(i)  = abs(prin_values(2)) * prin_dir(2,2);  % stress_I   vector projected in Y axis
	            stress_II_projection_x(i) = abs(prin_values(1)) * prin_dir(1,1);  % stress_II  vector projected in X axis
	            stress_II_projection_y(i) = abs(prin_values(1)) * prin_dir(2,1);  % stress_II  vector projected in Y axis

                % Abaqus gives you s_I within the range [0 s_Imax] and s_II
                % within the range [s_IImax 0]. To retrieve all the values,
                % I save the principal values on another variable here to
                % use it later on
                fem_results.stress_results(i,8) = prin_values(2);
                fem_results.stress_results(i,9) = prin_values(1);

            end

            stress_I_projection_z  = zeros(numel(stress_I_projection_x),1)                                 ;
            stress_II_projection_z = zeros(numel(stress_II_projection_x),1)                                ;

            stress_I_components    = [stress_I_projection_x  stress_I_projection_y  stress_I_projection_z ];
            stress_II_components   = [stress_II_projection_x stress_II_projection_y stress_II_projection_z];



            % ........ Write data ........
            if inputs.output_files_vtk_file==true

                if strcmp(inputs.fem_software,'pde_matlab')

                    filename = sprintf('./results/fem/solution_%05i.vtk', step); 
        
                    % Write fem variables (pde)
                    vtkwrite_fem(     filename , 'UNSTRUCTURED_GRID'   , nodes_x, nodes_y, nodes_z             , mesh_connectivity                        , ...
                                     'vectors' , 'displacement'        , u_disp                 ,    'scalars' , 'sigma_x_Pa'    , s_x                    , ...
                                     'scalars' , 'sigma_y_Pa'          , s_y                    ,    'scalars' , 'sigma_xy_Pa'   , s_xy                   , ... 
                                     'scalars' , 'sigma_I_Pa'          , s_I                    ,    'scalars' , 'sigma_II_Pa'   , s_II                   , ...
                                     'scalars' , 'sigma_Ix_Pa'         , stress_I_projection_x  ,    'scalars' , 'sigma_Ix_Pa'   , stress_I_projection_y  , ...
                                     'scalars' , 'sigma_IIx_Pa'        , stress_II_projection_x ,    'scalars' , 'sigma_IIx_Pa'  , stress_II_projection_y , ...
                                    'Precision', 2);

                elseif strcmp(inputs.fem_software,'abaqus')
    
                    filename = sprintf('./results/fem/solution_%05i.vtk', step); 
        
                    % Write fem: displacements at nodes
                    vtkwrite_fem(     filename , 'UNSTRUCTURED_GRID'   , nodes_x, nodes_y, nodes_z             , mesh_connectivity                        , ...
                                      'vectors' , 'displacement'       , u_disp, 'Precision', 2);


                    filename = sprintf('./results/fem/stress_%05i.vtk', step); 
        
                    % Write fem: stresses, value per element
                    vtkwrite_stress(   filename , 'POLYGONS'   , nodes_x, nodes_y, nodes_z             , mesh_connectivity                                 , ...
                                      'scalars' , 'sigma_x_Pa'         , s_x                     , 'scalars' , 'sigma_y_Pa'          , s_y                 , ...
                                      'scalars' , 'sigma_xy_Pa'        , s_xy                    , 'scalars' , 'sigma_I_Pa'          , s_I                 , ...
                                      'scalars' , 'sigma_II_Pa'        , s_II                    , 'vectors' , 'sigma_I_dir'         , stress_I_components , ...
                                      'vectors' , 'sigma_II_dir'       , stress_II_components    , ...
                                      'Precision', 2);

                end
    
            end

            % Save raw data
            if inputs.output_files_mat_file==true

                if strcmp(inputs.fem_software,'pde_matlab')
                
                    pde_model_results = fem_results.results_pde;
                    filename = sprintf('./results/raw_fem/pde_model_%05i.mat', step);
                    save(filename, 'pde_model_results')
                    
    
                 elseif strcmp(inputs.fem_software,'abaqus')

                    mesh_triangulation     = mesh_obj.triangulation_obj; 
                    triang_inf_id          = mesh_obj.elem_infect_id   ;  
                    triang_areas           = mesh_obj.triang_areas     ;
                    displacements_solution = fem_results.displ_results ;
                    stress_solution        = fem_results.stress_results;
                    filename = sprintf('./results/raw_fem/fem_model_%05i.mat', step);
                    save(filename, 'displacements_solution','stress_solution','mesh_triangulation', ...
                                   'triang_inf_id','triang_areas');

                end

            end






         
        end

    end




end