classdef inp_file

    methods

        %................... Object initialization ........................
        function obj = inp_file()
            addpath(fullfile('inputs')); 
            addpath(fullfile('src'));           
        end




        %................... Write inp ....................................
        function inp_obj = write_inp(inp_obj,pde_model,monolayer,mesh_obj)

            step     = monolayer.step                             ;

            filename = sprintf('./results/inp/inp_%05i.inp', step);

            %........ Open file ........
            fout     = fopen(filename, 'wt')                      ;


            % ........Heading ........
            fprintf(fout,'*Heading\n')                                            ;
            fprintf(fout,'** Job name: monolayer Model name: monolayer_model\n')  ;
            fprintf(fout,'** Generated by: monolayer_model matlab\n')             ;
            fprintf(fout,'*Preprint, echo=NO, model=NO, history=NO, contact=NO\n');
            fprintf(fout,'** PART: 2DMONOLAYER\n')                                ;
            fprintf(fout,'*FILE FORMAT, ASCII\n')                                 ;
            

            % ........ Choose linear o quadratic mesh ........
            geometric_order = inputs.fem_geometric_mesh_order                     ;
            % Generate :
            %           - same mesh if order = linear
            %           - new mesh if order = quadratic (adding middle edge
            %           points)
            %generateMesh(pde_model,"Hmin",100, "GeometricOrder",geometric_order) ;
            %generateMesh(pde_model,"Hmin",200, "GeometricOrder",geometric_order) ;




            % ........ Write nodes ........
            fprintf(fout,'*Node\n')                            ;           
            mesh_nodes   = mesh_obj.triangulation_obj.Points'  ; % nodes
            node_id      = 1:size(mesh_nodes,2)                ;
            var1         = mesh_nodes(1,:)                     ;
            var2         = mesh_nodes(2,:)                     ;
            output       = [node_id; var1; var2]               ;            
            fprintf(fout,'%d, %f, %f\n', output)               ;




            % ........ Write elements ........
            if strcmp(inputs.fem_geometric_mesh_order,'linear')
                fprintf(fout,'*Element, type=CPS3\n') ; 
            elseif strcmp(inputs.fem_geometric_mesh_order,'quadratic')    
                fprintf(fout,'*Element, type=CPS6M\n');
            end
            mesh_elements = mesh_obj.triangulation_obj.ConnectivityList';
            element_id    = 1:size(mesh_elements,2)                     ;


            if strcmp(inputs.fem_geometric_mesh_order,'linear')
                var1          = mesh_elements(1,:)                              ;
                var2          = mesh_elements(2,:)                              ;
                var3          = mesh_elements(3,:)                              ;
                output        = [element_id; var1; var2; var3]                  ; 
                fprintf(fout,'%d, %d, %d, %d\n', output)                        ;
            elseif strcmp(inputs.fem_geometric_mesh_order,'quadratic')
                var1          = mesh_elements(1,:)                              ;
                var2          = mesh_elements(2,:)                              ;
                var3          = mesh_elements(3,:)                              ;
                var4          = mesh_elements(4,:)                              ;
                var5          = mesh_elements(5,:)                              ;
                var6          = mesh_elements(6,:)                              ;
                output        = [element_id; var1; var2; var3; var4; var5; var6]; 
                fprintf(fout,'%d, %d, %d, %d, %d, %d, %d\n', output)            ;
            end




            % ........ Write Nsets and Elsets ........                
            fprintf(fout,'*Nset, nset=ALL_CELL_NODES, generate\n')     ;
            output     = [node_id(1); node_id(end)]                    ;
            fprintf(fout,'%d, %d, 1\n', output)                        ;  
            % generate means node_id(1):1:node_id(end), 


            fprintf(fout,'*Elset, elset=ALL_CELL_ELEMENTS, generate\n');
            output     = [element_id(1); element_id(end)]              ;
            fprintf(fout,'%d, %d, 1\n', output)                        ;  
            % generate means el_id(1):1:el_id(end)


            % This part is computed in two different scenarios
            %     (1) When the infection is fixed (infection over all
            %     steps)
            %     (2) When we simulate propagation, but first we let
            %     the program run a certain number of steps to avoid
            %     high forces during the first steps

            if inputs.monolayer_infection_active == true
                list_inf_cells       = monolayer.cell_infection;

                % Remove extruding cells from list_inf_cells
                list_extruding_cells = monolayer.cell_extrusion;
                if(any(list_extruding_cells))
                    indices                 = find(list_extruding_cells == 1);
                    list_inf_cells(indices) = []                             ;
                end

                if inputs.monolayer_infection_propagation == false || ...
                        (inputs.monolayer_infection_propagation == true && monolayer.step >= inputs.monolayer_infection_step_propagation && ...
                        any(list_inf_cells)~=false)
                    id_infected   = find(mesh_obj.elem_infect_id)                                      ;
                    id_uninfected = find(~mesh_obj.elem_infect_id)                                     ;
    
                    fprintf(fout,'*Elset, elset=UNINFECT_CELLS\n')                                     ;
                    fprintf(fout,'%d, %d, %d, %d, %d, %d\n', id_uninfected)                            ;  
                    fprintf(fout,'\n*Solid Section, elset=UNINFECT_CELLS, material=UNINFECTED_CELL\n') ;
                    fprintf(fout,',\n')                                                                ;
    
                    fprintf(fout,'*Elset, elset=INFECT_CELLS\n')                                       ;
                    fprintf(fout,'%d, %d, %d, %d, %d, %d\n', id_infected)                              ;
                    fprintf(fout,'\n*Solid Section, elset=INFECT_CELLS, material=INFECTED_CELL\n')     ;
                    fprintf(fout,',\n')                                                                ;

                else
                    fprintf(fout,'*Solid Section, elset=ALL_CELL_ELEMENTS, material=UNINFECTED_CELL\n');
                    fprintf(fout,',\n')                                                                ;
                end

            else
                    fprintf(fout,'*Solid Section, elset=ALL_CELL_ELEMENTS, material=UNINFECTED_CELL\n');
                    fprintf(fout,',\n')                                                                ;
            end

            


            fprintf(fout,'*Nset, nset=CELL_CENTROIDS, generate\n')   ;
            number_of_cells      = numel(monolayer.cell_id)          ;
            % Update number_of_cells in case of extruding cells
            list_extruding_cells = monolayer.cell_extrusion          ;
            if(any(list_extruding_cells))

                numb_extruding_cells   = sum(list_extruding_cells)             ;
                number_of_cells        = number_of_cells - numb_extruding_cells;

            end
            output          = [node_id(1); node_id(number_of_cells)] ;
            fprintf(fout,'%d, %d, 1\n', output)                      ;  
            % generate means node_id(1):1:node_id(end), 






            % ........ Write boundary set ........               
            boundary     = [0 0; inputs.monolayer_size_x 0;                ...
                           inputs.monolayer_size_x inputs.monolayer_size_y;...
                           0 inputs.monolayer_size_y];
            boundary_x   = boundary(:,1)             ;
            boundary_y   = boundary(:,2)             ;
            coord_x_node = mesh_nodes(1,:)'          ;
            coord_y_node = mesh_nodes(2,:)'          ;

            % Small numbers can lead to misclassifications (for
            % instance, 1.77e-15), make sure there is no such problem
            % at the border
            tol                                                                 = inputs.general_tolerance   ;
            coord_x_node(abs(coord_x_node) < tol)                               = 0                          ;
            coord_y_node(abs(coord_y_node) < tol)                               = 0                          ;
            coord_x_node(abs(coord_x_node - monolayer.monolayer_size(1)) < tol) = monolayer.monolayer_size(1);
            coord_y_node(abs(coord_y_node - monolayer.monolayer_size(2)) < tol) = monolayer.monolayer_size(2);
            

            % Tag the boundary nodes according to the chosen boundary 
            % conditions

            % ........ Find points that are on the border ........
            if strcmp(inputs.fem_BC_type,'fixed')

                [~,on]      = inpolygon(coord_x_node,coord_y_node, ... 
                                    boundary_x,boundary_y);
                bound_nodes = find(on == 1);
                
                fprintf(fout,'*Nset, nset=Boundary_nodes\n');
                fprintf(fout,'%d, %d, %d, %d\n', bound_nodes);

            elseif strcmp(inputs.fem_BC_type,'symmetric') || strcmp(inputs.fem_BC_type,'PBC') ||...
                   strcmp(inputs.fem_BC_type,'springs')   || strcmp(inputs.fem_BC_type,'material_around')
            
                % Find points that are on the left and right borders
                left_border       = abs(coord_x_node - boundary_x(1)) < tol ;
                right_border      = abs(coord_x_node - boundary_x(2)) < tol ;
    
                on_left_border    = find(left_border  == 1)                 ;
                on_right_border   = find(right_border == 1)                 ;                
                
    
                % Find points that are on the upper and bottom borders
                bottom_border     = abs(coord_y_node - boundary_y(1)) < tol ;
                top_border        = abs(coord_y_node - boundary_y(3)) < tol ;
    
                on_bottom_border  = find(bottom_border == 1)                ;
                on_top_border     = find(top_border    == 1)                ;

                % Extract the corner nodes separately
                left_bottom_corner  = intersect(on_left_border  , on_bottom_border);
                right_bottom_corner = intersect(on_right_border , on_bottom_border);
                left_upper_corner   = intersect(on_left_border  , on_top_border   );
                right_upper_corner  = intersect(on_right_border , on_top_border   );
               

                % ........ Print the nodes ........
                fprintf(fout, '*Nset, nset=Left_border_nodes\n')   ;
                fprintf(fout, '%d, %d, %d, %d\n', on_left_border)  ;
                fprintf(fout, '\n')                                ;
                
                fprintf(fout, '*Nset, nset=Right_border_nodes\n')  ;
                fprintf(fout, '%d, %d, %d, %d\n', on_right_border) ;
                fprintf(fout, '\n')                                ;
                
                fprintf(fout, '*Nset, nset=Bottom_border_nodes\n') ;
                fprintf(fout, '%d, %d, %d, %d\n', on_bottom_border);
                fprintf(fout, '\n')                                ;
                
                fprintf(fout, '*Nset, nset=Top_border_nodes\n')    ;
                fprintf(fout, '%d, %d, %d, %d\n', on_top_border)   ;
                fprintf(fout, '\n')                                ;

                % ........ Corners ........
                fprintf(fout, '*Nset, nset=left_bottom_corner\n')  ;
                fprintf(fout, '%d\n', left_bottom_corner)          ;

                fprintf(fout, '*Nset, nset=right_bottom_corner\n') ;
                fprintf(fout, '%d\n', right_bottom_corner)         ;

                fprintf(fout, '*Nset, nset=left_upper_corner\n')   ;
                fprintf(fout, '%d\n', left_upper_corner)           ;

                fprintf(fout, '*Nset, nset=right_upper_corner\n')  ;
                fprintf(fout, '%d\n', right_upper_corner)          ;

            end



            % Create sets depending on the boundary condition chosen
            boundary_condition = inputs.getBoundaryCondition(monolayer.step)       ;
            if strcmp(boundary_condition,'PBC')
                % Define surface to use it later on the contact (PBC)
                % Left edge
                fprintf(fout, '*Surface, type=NODE, name=master_surf_left\n')      ;
                fprintf(fout, 'Left_border_nodes\n')                               ;
    
                % Right edge
                fprintf(fout, '*Surface, type=NODE, name=master_surf_right\n')     ;
                fprintf(fout, 'Right_border_nodes\n')                              ;

                % Upper edge
                fprintf(fout, '*Surface, type=NODE, name=master_surf_upper\n')     ;
                fprintf(fout, 'Top_border_nodes\n')                                ;
    
                % Bottom edge
                fprintf(fout, '*Surface, type=NODE, name=master_surf_bottom\n')    ;
                fprintf(fout, 'Bottom_border_nodes\n')                             ;
        

                % Include PBC file input
                fprintf(fout, '** include\n')                                      ;
                fprintf(fout, '*Include, input=PBC_mesh.inp\n')                    ;

    
                % Contact between cell monolayer and PBC boundaries
                fprintf(fout, '** Constraint: contact_left_border_PBC\n')          ;
                fprintf(fout, '*Tie, name=contact_left_border_PBC, adjust=no\n')   ;
                fprintf(fout, 'slave_surf_left, master_surf_left\n')               ;
    
                fprintf(fout, '** Constraint: contact_right_border_PBC\n')         ;
                fprintf(fout, '*Tie, name=contact_right_border_PBC, adjust=no\n')  ;
                fprintf(fout, 'slave_surf_right, master_surf_right\n')             ;

                fprintf(fout, '** Constraint: contact_upper_border_PBC\n')         ;
                fprintf(fout, '*Tie, name=contact_upper_border_PBC, adjust=no\n')  ;
                fprintf(fout, 'slave_surf_upper, master_surf_upper\n')             ;
    
                fprintf(fout, '** Constraint: contact_bottom_border_PBC\n')        ;
                fprintf(fout, '*Tie, name=contact_bottom_border_PBC, adjust=no\n') ;
                fprintf(fout, 'slave_surf_bottom, master_surf_bottom\n')           ;

            end

            if strcmp(inputs.fem_BC_type,'springs')
                % Define surface to use it later on the BC definition
                % Left edge
                fprintf(fout, '*Surface, type=NODE, name=left_edge_surface\n')     ;
                fprintf(fout, 'Left_border_nodes\n')                               ;

                % Right edge
                fprintf(fout, '*Surface, type=NODE, name=right_edge_surface\n')    ;
                fprintf(fout, 'Right_border_nodes\n')                              ;
    
                % Upper edge
                fprintf(fout, '*Surface, type=NODE, name=upper_edge_surface\n')    ;
                fprintf(fout, 'Top_border_nodes\n')                                ;

                % Bottom edge
                fprintf(fout, '*Surface, type=NODE, name=bottom_edge_surface\n')   ;
                fprintf(fout, 'Bottom_border_nodes\n')                             ;


            end

            if strcmp(boundary_condition,'material_around')
                % Define surface to use it later on the contact
                % Left edge
                fprintf(fout, '*Surface, type=NODE, name=master_surf_left\n')      ;
                fprintf(fout, 'Left_border_nodes\n')                               ;
    
                % Right edge
                fprintf(fout, '*Surface, type=NODE, name=master_surf_right\n')     ;
                fprintf(fout, 'Right_border_nodes\n')                              ;

                % Upper edge
                fprintf(fout, '*Surface, type=NODE, name=master_surf_upper\n')     ;
                fprintf(fout, 'Top_border_nodes\n')                                ;
    
                % Bottom edge
                fprintf(fout, '*Surface, type=NODE, name=master_surf_bottom\n')    ;
                fprintf(fout, 'Bottom_border_nodes\n')                             ;
        

                % Include PBC file input
                fprintf(fout, '** include\n')                                      ;
                fprintf(fout, '*Include, input=grid500.inp\n')                     ;

    
                % Contact between cell monolayer and material around
                fprintf(fout, '** Constraint: contact_left_border\n')              ;
                fprintf(fout, '*Tie, name=contact_left_border_PBC, adjust=no\n')   ;
                fprintf(fout, 'left_surface_out, master_surf_left\n')              ;
    
                fprintf(fout, '** Constraint: contact_right_border\n')             ;
                fprintf(fout, '*Tie, name=contact_right_border_PBC, adjust=no\n')  ;
                fprintf(fout, 'right_surface_out, master_surf_right\n')            ;

                fprintf(fout, '** Constraint: contact_upper_border\n')             ;
                fprintf(fout, '*Tie, name=contact_upper_border_PBC, adjust=no\n')  ;
                fprintf(fout, 'upper_surface_out, master_surf_upper\n')            ;
    
                fprintf(fout, '** Constraint: contact_bottom_border\n')            ;
                fprintf(fout, '*Tie, name=contact_bottom_border_PBC, adjust=no\n') ;
                fprintf(fout, 'bottom_surface_out, master_surf_bottom\n')          ;

            end




            % ........ Write materials ........
            uninf_young = inputs.fem_properties_E_uninfected    ;
            inf_young   = inputs.fem_properties_E_infected      ;
            poisson_mod = inputs.fem_properties_poisson_ratio   ;

            fprintf(fout,'** MATERIALS'+"\n")                   ;
            fprintf(fout,'*Material, name=UNINFECTED_CELL'+"\n");
            fprintf(fout,'*Elastic'+"\n")                       ;           
            fprintf(fout,'%f,%f\n',uninf_young,poisson_mod)     ;

            fprintf(fout,'*Material, name=INFECTED_CELL'+"\n")  ;
            fprintf(fout,'*Elastic'+"\n")                       ;           
            fprintf(fout,'%f,%f\n',inf_young,poisson_mod)       ;





            % ........ Write boundary conditions ........
            if strcmp(inputs.fem_BC_type,'fixed')

                fprintf(fout,'** BOUNDARY CONDITIONS\n')                                   ;
                fprintf(fout,'** Name: BC_exterior Type: Symmetry/Antisymmetry/Encastre\n');
                fprintf(fout,'*Boundary\n')                                                ;
                fprintf(fout,'Boundary_nodes, ENCASTRE\n')                                 ;

            elseif strcmp(boundary_condition,'symmetric')

                fprintf(fout,'** BOUNDARY CONDITIONS\n')                                   ;
                %fprintf(fout,'** Name: BC_symmetry Type: Symmetry\n');
                %fprintf(fout,'*Boundary\n');
                %fprintf(fout,'Right_border_nodes, SYMMETRY\n' );
                %fprintf(fout,'Left_border_nodes, SYMMETRY\n'  );
                %fprintf(fout,'Top_border_nodes, SYMMETRY\n'   );
                %fprintf(fout,'Bottom_border_nodes, SYMMETRY\n');


                fprintf(fout,'** Name: BC_bottom Type: Symmetry/Antisymmetry/Encastre\n')  ;
                fprintf(fout,'*Boundary\n')                                                ;
                fprintf(fout,'Bottom_border_nodes, YSYMM\n')                               ;
                fprintf(fout,'** Name: BC_top Type: Symmetry/Antisymmetry/Encastre\n'   )  ;
                fprintf(fout,'*Boundary\n')                                                ;
                fprintf(fout,'Top_border_nodes, YSYMM\n'   )                               ;


                fprintf(fout,'** Name: BC_left Type: Symmetry/Antisymmetry/Encastre\n'  )  ;
                fprintf(fout,'*Boundary\n')                                                ;
                fprintf(fout,'Left_border_nodes, XSYMM\n' )                                ;
                fprintf(fout,'** Name: BC_right Type: Symmetry/Antisymmetry/Encastre\n' )  ;
                fprintf(fout,'*Boundary\n')                                                ;
                fprintf(fout,'Right_border_nodes, XSYMM\n' )                               ;



            elseif strcmp(boundary_condition,'PBC')

                fprintf(fout,'** BOUNDARY CONDITIONS\n')                                              ;

                % When one applies PBC, the BC at the corners must be specified
                fprintf(fout,'** Name: left_bottom_corner_x Type: Symmetry/Antisymmetry/Encastre\n')  ;
                fprintf(fout,'*Boundary\n')                                                           ;
                fprintf(fout,'left_bottom_corner, XSYMM\n')                                           ;

                fprintf(fout,'** Name: left_bottom_corner_y Type: Symmetry/Antisymmetry/Encastre\n')  ;
                fprintf(fout,'*Boundary\n')                                                           ;
                fprintf(fout,'left_bottom_corner, YSYMM\n')                                           ;

                fprintf(fout,'** Name: right_bottom_corner_x Type: Symmetry/Antisymmetry/Encastre\n') ;
                fprintf(fout,'*Boundary\n')                                                           ;
                fprintf(fout,'right_bottom_corner, XSYMM\n')                                          ;

                fprintf(fout,'** Name: right_bottom_corner_y Type: Symmetry/Antisymmetry/Encastre\n') ;
                fprintf(fout,'*Boundary\n')                                                           ;
                fprintf(fout,'right_bottom_corner, YSYMM\n')                                          ;

                fprintf(fout,'** Name: upper_left_corner_x Type: Symmetry/Antisymmetry/Encastre\n')   ;
                fprintf(fout,'*Boundary\n')                                                           ;
                fprintf(fout,'left_upper_corner, XSYMM\n')                                            ;

                fprintf(fout,'** Name: upper_left_corner_y Type: Symmetry/Antisymmetry/Encastre\n')   ;
                fprintf(fout,'*Boundary\n')                                                           ;
                fprintf(fout,'left_upper_corner, YSYMM\n')                                            ;

                fprintf(fout,'** Name:upper_right_corner_x Type: Symmetry/Antisymmetry/Encastre\n')   ;
                fprintf(fout,'*Boundary\n')                                                           ;
                fprintf(fout,'right_upper_corner, XSYMM\n')                                           ;

                fprintf(fout,'** Name: upper_right_corner_y Type: Symmetry/Antisymmetry/Encastre\n')  ;
                fprintf(fout,'*Boundary\n')                                                           ;
                fprintf(fout,'right_upper_corner, YSYMM\n')                                           ;

            elseif strcmp(inputs.fem_BC_type,'springs')

                fprintf(fout,'**\n')                                                               ;
                fprintf(fout,'** INTERACTIONS\n')                                                  ;
                fprintf(fout,'**\n')                                                               ;
                fprintf(fout,'** Interaction: springs\n')                                          ;
                fprintf(fout,'*Foundation\n')                                                      ;
                fprintf(fout,'left_edge_surface, F1, %d.\n', inputs.fem_properties_E_uninfected)   ;
                fprintf(fout,'upper_edge_surface, F2, %d.\n', inputs.fem_properties_E_uninfected)  ;
                fprintf(fout,'right_edge_surface, F1, %d.\n', inputs.fem_properties_E_uninfected)  ;
                fprintf(fout,'bottom_edge_surface, F2, %d.\n', inputs.fem_properties_E_uninfected) ;

            end



            % ........ Write step ........
            fprintf(fout,'**STEP: Step-1\n')               ;
            fprintf(fout,'*Step, name=Step-1, nlgeom=NO\n');
            fprintf(fout,'*Static\n')                      ;
            fprintf(fout,'%d,%d,%d,%d\n',1.,1.,1e-05,1.)   ;







            % ........ Write concentrated forces ........

            % Parse forces: initial model (vertices, centroids id) might
            % not be the same as the final model when generating the mesh
            % with pde (for example, when you choose quadratic mesh)

            % Centroid forces
            list_centroids       = monolayer.monolayer_cell_points  ;
            numb_cells           = numel(monolayer.cell_id)         ;
            centroid_nodes       = (1:numb_cells)'                  ;
            cent_force_map       = monolayer.force_list(:,1:2)      ;

            % Update lists in case of extruding cells
            list_extruding_cells = monolayer.cell_extrusion;
            if(any(list_extruding_cells))
                indices                   = find(list_extruding_cells == 1)                ;
                list_centroids(indices,:) = []                                             ;
                cent_force_map(indices,:) = []                                             ;

                numb_extruding_cells      = sum(list_extruding_cells)                      ;
                numb_cells                = numel(monolayer.cell_id) - numb_extruding_cells;
                centroid_nodes            = (1:numb_cells)'                                ;

            end

            % Vertices forces
            list_vertices  = monolayer.voronoi_vertices       ;
            % (1) Remove vertices that have been removed for the triangulation,
            % eg. vertices that are close to the centroid on the monolayer border
            % (they generate negative areas)
            removed_vector = mesh_obj.vertex_removed          ;

            if ~isempty(removed_vector)
                list_vertices(removed_vector,:) = [];
            end

            global_list    = [list_centroids; list_vertices]  ;

            numb_vertices  = numel(list_vertices(:,1))        ;    

            vertex_nodes   = (numb_cells + 1:     ...
                              numb_cells +        ...
                              numb_vertices)'                 ;
            vert_force_map = monolayer.vertices_forces(:,1:2) ;

            % Same as in (1)
            if ~isempty(removed_vector)
                vert_force_map(removed_vector,:) = [];
            end

            % Assembly force maps (according to initial mesh)
            id_nodes       = [centroid_nodes; vertex_nodes]   ;
            force_map      = [cent_force_map; vert_force_map] ;

%             % Parse cells to the model according to the new mesh 
%             % organization (it might be different from the initial one)
%             B = global_list;                           % initial mesh/configuration
%             A = pde_model.Mesh.Nodes';                 % new     mesh/configuration 
%             [~,Locb] = ismember(A,B,'rows');
% 
%             % Order new force_map according to the new mesh
%             new_force_map = force_map(Locb,:);


            fprintf(fout,'** LOADS\n')                                    ;
            fprintf(fout,'** Name: Load_model Type: Concentrated force\n');
            fprintf(fout,'*Cload\n')                                      ;

            output_x = [id_nodes'; force_map(:,1)']                       ;
            output_y = [id_nodes'; force_map(:,2)']                       ;
            fprintf(fout,'%d, 1, %f\n', output_x)                         ;  
            fprintf(fout,'%d, 2, %f\n', output_y)                         ;

            


            % ........ Output variables ........
            fprintf(fout,'** OUTPUT REQUESTS\n')                                         ;
            fprintf(fout,'*Restart, write, frequency=0\n')                               ;
            fprintf(fout,'** FIELD OUTPUT: F-Output-1n\n')                               ;
%             fprintf(fout,'*Output, field, variable=PRESELECT\n')                         ;
%             fprintf(fout,'*Output, history, variable=PRESELECT\n')                       ;
            


            fprintf(fout,'*Output, field\n')                                             ;
            fprintf(fout,'*Node Output\n')                                               ;
            fprintf(fout,'CF, U\n')                                                      ;
            %fprintf(fout,'*Element Output, POSITION=INTEGRATION POINTS, directions=NO\n');
            fprintf(fout,'*Element Output, POSITION=CENTROID, directions=NO\n')          ;
            % print integration points/ centroids, requested on element output:
            fprintf(fout,'COORD, S,\n')                                                  ;  
            %fprintf(fout,'COORD, EE, S,\n');  
            fprintf(fout,'**\n')                                                         ;
            fprintf(fout,'** HISTORY OUTPUT: H-Output-1\n')                              ;
            fprintf(fout,'**\n')                                                         ;
            fprintf(fout,'*Output, history, variable=PRESELECT\n')                       ;
%             fprintf(fout,'* Node print, nset=CELL_CENTROIDS, frequency=1\n'); % FOR THE .DAT
%             fprintf(fout,'COOR1, COOR2, U1, U2,\n')                                      ;
%             fprintf(fout,'* EL print, POSITION=INTEGRATION POINTS , frequency=1'+"\n")   ;
%             fprintf(fout,'COORD,'+"\n")                                                  ;
%             fprintf(fout,'EE11, EE22, EE12,'+"\n")                                       ;
            fprintf(fout,'*End Step\n')                                                  ;

            % *Element Output, POSITION = values are being written at:
            %                             - INTEGRATION POINTS (default)
            %                             - CENTROIDAL:  centroid of the element
            %                             - NODES     :  extrapolated to the nodes
            %                  DIRECTIONS = write the element material directions to the output database
            %                             - YES
            %                             - NO
            %                  VARIABLE  = 
            %                             - ALL      :  all element variables
            %                             - PRESELECT


            fclose(fout);

          
        end
        
    end




end