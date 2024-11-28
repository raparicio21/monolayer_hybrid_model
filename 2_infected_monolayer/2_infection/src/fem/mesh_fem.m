classdef mesh_fem

    properties
        triangulation_obj        % triangulation object
        elem_infect_id           % id to label infection in the elements
        elem_cell_id             % id to label cells from ABM in the elements
        triang_areas             % list of areas of mesh triangles
        vertex_removed           % list of vertex to be removed (eg vertex that are really close to the centroid on the border)
    end

    methods

        %................... Object initialization ........................
        function obj = mesh_fem()
            addpath(fullfile('inputs')); 
            addpath(fullfile('src'));            
        end




        %.................. Initialize the mesh ...........................
        function mesh_obj = initialize(mesh_obj,monolayer,ecm)

            % Create triangulation from Monolayer (ABM)
            my_triang        = create_triang();
            [triang, inf_id,cell_id, triangles_areas, vertex_to_be_removed] = my_triang.createTriangulation(monolayer,ecm);

            % Update triangulation
            mesh_obj.triangulation_obj = [];
            mesh_obj.triangulation_obj = triang;

            mesh_obj.elem_infect_id    = [];
            mesh_obj.elem_infect_id    = inf_id;

            mesh_obj.elem_cell_id      = [];
            mesh_obj.elem_cell_id      = cell_id;

            mesh_obj.triang_areas      = [];
            mesh_obj.triang_areas      = triangles_areas;

            mesh_obj.vertex_removed    = [];
            mesh_obj.vertex_removed    = vertex_to_be_removed;
       
        end


        %............... Create periodic boundary mesh ....................
        function mesh_obj = period_boundary_mesh(mesh_obj,monolayer)

            % (1) Get monolayer size
            x_size = monolayer.monolayer_size(1);
            y_size = monolayer.monolayer_size(2);

            % (2) Define the dimensions of the new rectangle/grid
            rectangle_height   = y_size;
            rectangle_width    = 0.5   ;
            element_size       = 10    ; % microns
            dragging_dimension = 0.01  ;

            % (3) Define node coordinates of the rectangle on the left side
            % = rectangle 1
            node_Ycoord_rect            = 0:element_size:rectangle_height                               ;

            node_Xcoord_rect_left_side  = 0 - dragging_dimension - rectangle_width                      ;
            node_Xcoord_rect_right_side = 0 - dragging_dimension                                        ;

            node_Xcoord_rect_left_side  = node_Xcoord_rect_left_side  * ones(1,length(node_Ycoord_rect));
            node_Xcoord_rect_right_side = node_Xcoord_rect_right_side * ones(1,length(node_Ycoord_rect));

            rect_1_left_side_coord      = [node_Xcoord_rect_left_side ; node_Ycoord_rect]'              ;
            rect_1_right_side_coord     = [node_Xcoord_rect_right_side; node_Ycoord_rect]'              ;

            % (4) Define node coordinates of the rectangle on the right
            % side = rectangle 2
            node_Xcoord_rect_right_side = x_size + dragging_dimension + rectangle_width                 ;
            node_Xcoord_rect_left_side  = x_size + dragging_dimension                                   ; 

            node_Xcoord_rect_left_side  = node_Xcoord_rect_left_side  * ones(1,length(node_Ycoord_rect));
            node_Xcoord_rect_right_side = node_Xcoord_rect_right_side * ones(1,length(node_Ycoord_rect));

            rect_2_left_side_coord      = [node_Xcoord_rect_left_side ; node_Ycoord_rect]'              ;
            rect_2_right_side_coord     = [node_Xcoord_rect_right_side; node_Ycoord_rect]'              ;

            

            % Now upper and bottom rectangles
            % (2') Define the dimensions of the new rectangle/grid
            rectangle_height   = 0.5   ;
            rectangle_width    = x_size;
            element_size       = 10    ; % microns
            dragging_dimension = 0.01  ;

            % (3') Define node coordinates of the rectangle on the upper side
            % = rectangle 3
            node_Xcoord_rect_3          = 0:element_size:rectangle_width                                   ;

            node_Ycoord_rect_upper_side = y_size + dragging_dimension + rectangle_height                   ;
            node_Ycoord_rect_bott_side  = y_size + dragging_dimension                                      ;

            node_Ycoord_rect_upper_side = node_Ycoord_rect_upper_side  * ones(1,length(node_Xcoord_rect_3));
            node_Ycoord_rect_bott_side  = node_Ycoord_rect_bott_side   * ones(1,length(node_Xcoord_rect_3));

            rect_3_upper_side_coord     = [node_Xcoord_rect_3 ; node_Ycoord_rect_upper_side ]'             ;
            rect_3_bott_side_coord      = [node_Xcoord_rect_3 ; node_Ycoord_rect_bott_side  ]'             ;

            % (4') Define node coordinates of the rectangle on the bottom
            % side = rectangle 4
            node_Ycoord_rect_bott_side  = -dragging_dimension - rectangle_height                           ; 
            node_Ycoord_rect_upper_side = -dragging_dimension                                              ; 

            node_Ycoord_rect_bott_side  = node_Ycoord_rect_bott_side  * ones(1,length(node_Xcoord_rect_3)) ;
            node_Ycoord_rect_upper_side = node_Ycoord_rect_upper_side * ones(1,length(node_Xcoord_rect_3)) ;

            rect_4_bott_side_coord      = [node_Xcoord_rect_3 ; node_Ycoord_rect_bott_side ]'              ; 
            rect_4_upper_side_coord     = [node_Xcoord_rect_3 ; node_Ycoord_rect_upper_side]'              ;

            % (5) Gather all nodes
            mesh_nodes                  = [rect_1_left_side_coord ; rect_1_right_side_coord ;        ...
                                           rect_2_left_side_coord ; rect_2_right_side_coord ;        ...
                                           rect_3_upper_side_coord; rect_3_bott_side_coord  ;        ...
                                           rect_4_upper_side_coord; rect_4_bott_side_coord   ]'           ;
            

            % (6) Generate the connectivity
            % Left rectangle
            first_element_node      =  1    :   length(node_Ycoord_rect)-1;
            second_element_node     = (length(node_Ycoord_rect)+1)  : (length(node_Ycoord_rect)*2-1);
            third_element_node      = (length(node_Ycoord_rect)+2)  : (length(node_Ycoord_rect)*2  );
            fourth_element_node     =  1+1  : length(node_Ycoord_rect)    ;

            mesh_connectivity_left  = [first_element_node;second_element_node;third_element_node;fourth_element_node];
            element_id_left         = 1:size(mesh_connectivity_left,2);

            % Right rectangle
            number_nodes_left       = size(rect_1_left_side_coord,1);
            mesh_connectivity_right = mesh_connectivity_left + number_nodes_left*2;
            element_id_right        = 1:size(mesh_connectivity_right,2); 
            element_id_right        = element_id_right + numel(element_id_left);

            % Upper rectangle
            number_nodes_so_far     = size(rect_1_left_side_coord,1)*4;  
            fourth_element_node     = number_nodes_so_far + 1 : number_nodes_so_far + size(rect_3_upper_side_coord,1)-1;
            third_element_node      = fourth_element_node + 1 ;
            first_element_node      = fourth_element_node + size(rect_3_upper_side_coord,1);
            second_element_node     = first_element_node  + 1 ;

            mesh_connectivity_up    = [first_element_node;second_element_node;third_element_node;fourth_element_node];
            element_id_up           = 1:size(mesh_connectivity_up,2); 
            element_id_up           = element_id_up + numel(element_id_left) + numel(element_id_right);

            
            % Bottom rectangle
            number_nodes_up         = size(rect_3_upper_side_coord,1);
            mesh_connectivity_down  = mesh_connectivity_up + number_nodes_up*2;        
            element_id_down         = 1:size(mesh_connectivity_down,2); 
            element_id_down         = element_id_down + numel(element_id_left) + numel(element_id_right) + numel(element_id_up);



            % All rectangles
            mesh_connectivity       = [ mesh_connectivity_left mesh_connectivity_right mesh_connectivity_up mesh_connectivity_down];
            element_id              = [ element_id_left element_id_right element_id_up element_id_down];

            



 


            % Write inp file for the new mesh (periodic boundary
            % conditions PBC)
            filename = sprintf('./results/inp/PBC_mesh.inp');

            % Open file
            fout = fopen(filename, 'wt');

            overlap_number = 1000000; % add 1 million to not overlap nodes
            % Write nodes
            
            fprintf(fout,'*Node\n');           
            node_id    = 1:size(mesh_nodes,2) ;  
            var1       = mesh_nodes(1,:);
            var2       = mesh_nodes(2,:);
            output     = [node_id + overlap_number; var1; var2];            
            fprintf(fout,'%d, %f, %f\n', output);

            % Write elements
            fprintf(fout,'*Element, type=CPS4\n');

            % Write elements
            var1          = mesh_connectivity(1,:) + overlap_number;
            var2          = mesh_connectivity(2,:) + overlap_number;
            var3          = mesh_connectivity(3,:) + overlap_number;
            var4          = mesh_connectivity(4,:) + overlap_number;
            output        = [element_id + overlap_number; var1; var2; var3; var4]; 
            fprintf(fout,'%d, %d, %d, %d, %d\n', output);

            % Write Nsets and Elsets 
            % All nodes
            fprintf(fout,'*Nset, nset=ALL_PBC_NODES, generate\n');
            output = [node_id(1); node_id(end)];
            output = output + overlap_number;
            fprintf(fout,'%d, %d, 1\n', output);  
            % generate means node_id(1):1:node_id(end), 

            % Bottom nodes rectangles 1 and 2 (left and right)
            fprintf(fout, '*Nset, nset=PBC_bottom_nodes\n');
            node_1_left_rect  = 1                                        + overlap_number;
            node_2_left_rect  = length(node_Ycoord_rect)+1               + overlap_number;
            node_1_right_rect = node_1_left_rect + number_nodes_left*2                   ;
            node_2_right_rect = node_2_left_rect + number_nodes_left*2                   ;
            bottom_nodes      = [node_1_left_rect node_2_left_rect node_1_right_rect node_2_right_rect];
            fprintf(fout, '%d, %d, %d, %d\n', bottom_nodes);

            % Top nodes rectangles 1 and 2 (left and right)
            fprintf(fout, '*Nset, nset=PBC_top_nodes\n');
            node_1_left_rect  = length(node_Ycoord_rect)                 + overlap_number;
            node_2_left_rect  = length(node_Ycoord_rect)*2               + overlap_number;
            node_1_right_rect = node_1_left_rect + number_nodes_left*2                   ;
            node_2_right_rect = node_2_left_rect + number_nodes_left*2                   ;
            top_nodes         = [node_1_left_rect node_2_left_rect node_1_right_rect node_2_right_rect];
            fprintf(fout, '%d, %d, %d, %d\n', top_nodes);

            % Left nodes rectangles 3 and 4 (up and down)
            fprintf(fout, '*Nset, nset=PBC_left_nodes\n');
            node_1_upp_rect  = number_nodes_so_far + 1                         + overlap_number;
            node_2_upp_rect  = node_1_upp_rect + length(node_Xcoord_rect_3)                    ;
            node_1_down_rect = node_1_upp_rect + length(node_Xcoord_rect_3)  *2                ;
            node_2_down_rect = node_2_upp_rect + length(node_Xcoord_rect_3)  *2                ;
            left_nodes       = [node_1_upp_rect node_2_upp_rect node_1_down_rect node_2_down_rect];
            fprintf(fout, '%d, %d, %d, %d\n', left_nodes);

            % Right nodes rectangles 3 and 4 (up and down)
            fprintf(fout, '*Nset, nset=PBC_right_nodes\n');
            node_1_upp_rect  = number_nodes_so_far + length(node_Xcoord_rect_3) + overlap_number;
            node_2_upp_rect  = node_1_upp_rect + length(node_Xcoord_rect_3)                     ;
            node_1_down_rect = node_1_upp_rect + length(node_Xcoord_rect_3)  *2                 ;
            node_2_down_rect = node_2_upp_rect + length(node_Xcoord_rect_3)  *2                 ;
            right_nodes      = [node_1_upp_rect node_2_upp_rect node_1_down_rect node_2_down_rect];
            fprintf(fout, '%d, %d, %d, %d\n', right_nodes);

            % Elset all nodes
            fprintf(fout,'*Elset, elset=ALL_PBC_ELEMENTS, generate\n');
            output = [element_id(1); element_id(end)];
            output = output + overlap_number;
            fprintf(fout,'%d, %d, 1\n', output);  
            % generate means el_id(1):1:el_id(end)

            fprintf(fout,'*Solid Section, elset=ALL_PBC_ELEMENTS, material=PBC_material\n');
            fprintf(fout,',\n');




            % Write nodes involved in the periodic boundary condition for
            % the equations (NSET needed)
            % Left border
            for i = 1:length(node_Ycoord_rect)
                fprintf(fout, '*Nset, nset=N%d\n', i + overlap_number);         
                fprintf(fout, '%d,\n', i + overlap_number);
            end

            % Right border
            % Bottom
            node_2_left_rect         = length(node_Ycoord_rect)+1              + overlap_number;
            node_2_right_rect_bottom = node_2_left_rect + number_nodes_left*2                  ;
            %Top
            node_2_left_rect         = length(node_Ycoord_rect)*2              + overlap_number;
            node_2_right_rect_top    = node_2_left_rect + number_nodes_left*2                  ;

            for i = node_2_right_rect_bottom:node_2_right_rect_top
                fprintf(fout, '*Nset, nset=N%d\n', i);         
                fprintf(fout, '%d,\n', i);
            end
            

            % Top border
            first_node        = number_nodes_so_far + 1                   + overlap_number;
            last_node         = first_node + length(node_Xcoord_rect_3) -1                ;

            for i = first_node:last_node
                fprintf(fout, '*Nset, nset=N%d\n', i);         
                fprintf(fout, '%d,\n', i);
            end

            % Bottom border
            first_node        = node_1_upp_rect + 1 + length(node_Xcoord_rect_3)  *2      ;
            last_node         = first_node + length(node_Xcoord_rect_3) -1                ;

            for i = first_node:last_node
                fprintf(fout, '*Nset, nset=N%d\n', i);         
                fprintf(fout, '%d,\n', i);
            end



            % Define surface to use it later on the contact
            % Left edge
            node1 = length(node_Ycoord_rect)+1       + overlap_number;
            node2 = length(node_Ycoord_rect)*2       + overlap_number;
            nodes = [node1 node2];
            fprintf(fout, '*Nset, nset=left_slave_surf, generate\n');
            fprintf(fout, '%d, %d, 1\n', nodes);
            fprintf(fout, '*Surface, type=NODE, name=slave_surf_left\n');
            fprintf(fout, 'left_slave_surf\n');

            % Right edge
            node1 = length(node_Ycoord_rect)*2+1     + overlap_number;
            node2 = length(node_Ycoord_rect)*3       + overlap_number;
            nodes = [node1 node2];
            fprintf(fout, '*Nset, nset=right_slave_surf, generate\n');
            fprintf(fout, '%d, %d, 1\n', nodes);
            fprintf(fout, '*Surface, type=NODE, name=slave_surf_right\n');
            fprintf(fout, 'right_slave_surf\n');

            % Upper edge
            node1 = number_nodes_so_far + 1 + length(node_Xcoord_rect_3)  + overlap_number;
            node2 = node1 + length(node_Xcoord_rect_3) - 1;
            nodes = [node1 node2];
            fprintf(fout, '*Nset, nset=upper_slave_surf, generate\n');
            fprintf(fout, '%d, %d, 1\n', nodes);
            fprintf(fout, '*Surface, type=NODE, name=slave_surf_upper\n');
            fprintf(fout, 'upper_slave_surf\n');

            % Bottom edge
            node1 = number_nodes_so_far + 1 + length(node_Xcoord_rect_3)*2  + overlap_number;
            node2 = node1 + length(node_Xcoord_rect_3) -1 ;
            nodes = [node1 node2];
            fprintf(fout, '*Nset, nset=bottom_slave_surf, generate\n');
            fprintf(fout, '%d, %d, 1\n', nodes);
            fprintf(fout, '*Surface, type=NODE, name=slave_surf_bottom\n');
            fprintf(fout, 'bottom_slave_surf\n');





            % Define constrains left and right rectangle
            node_ref_right = length(node_Ycoord_rect)*3 + 1  + overlap_number;
            node_ref_left  = 1                               + overlap_number;
            for i = 2:length(node_Ycoord_rect)
                fprintf(fout, '** Constraint: eq%dx\n', i);
                fprintf(fout, '*Equation\n');
                fprintf(fout, '4\n');
                node1 = length(node_Ycoord_rect)*3+i         + overlap_number;
                node2 = i                                    + overlap_number;

                fprintf(fout, 'N%d, 1, 1.\n', node1);
                fprintf(fout, 'N%d, 1, -1.\n', node2);
                fprintf(fout, 'N%d, 1, -1.\n', node_ref_right);
                fprintf(fout, 'N%d, 1, 1.\n', node_ref_left);

                fprintf(fout, '** Constraint: eq%dy\n', i);
                fprintf(fout, '*Equation\n');
                fprintf(fout, '4\n');
                fprintf(fout, 'N%d, 2, 1.\n', node1);
                fprintf(fout, 'N%d, 2, -1.\n', node2);
                fprintf(fout, 'N%d, 2, -1.\n', node_ref_right);
                fprintf(fout, 'N%d, 2, 1.\n', node_ref_left);


            end



            % Define constrains upper and bottom rectangle
            node_ref_up      = number_nodes_so_far + 1         + overlap_number     ;
            node_ref_down    = node_1_upp_rect + 1 + length(node_Xcoord_rect_3)  *2 ;

            for i = 1  :  length(node_Xcoord_rect_3)-1
                fprintf(fout, '** Constraint: eq_up_%dx\n', i);
                fprintf(fout, '*Equation\n');
                fprintf(fout, '4\n');

                node1 = node_ref_up   + i ;
                node2 = node_ref_down + i ;

                fprintf(fout, 'N%d, 1, 1.\n', node1);
                fprintf(fout, 'N%d, 1, -1.\n', node2);
                fprintf(fout, 'N%d, 1, -1.\n', node_ref_up);
                fprintf(fout, 'N%d, 1, 1.\n', node_ref_down);

                fprintf(fout, '** Constraint: eq_up_%dy\n', i);
                fprintf(fout, '*Equation\n');
                fprintf(fout, '4\n');
                fprintf(fout, 'N%d, 2, 1.\n', node1);
                fprintf(fout, 'N%d, 2, -1.\n', node2);
                fprintf(fout, 'N%d, 2, -1.\n', node_ref_up);
                fprintf(fout, 'N%d, 2, 1.\n', node_ref_down);


            end





            % Write materials
            uninf_young = inputs.fem_properties_E_uninfected;
            PBC_young   = uninf_young/100;
            poisson_mod = inputs.fem_properties_poisson_ratio;
            
            fprintf(fout,'** MATERIALS'+"\n");
            fprintf(fout,'*Material, name=PBC_material'+"\n");
            fprintf(fout,'*Elastic'+"\n");           
            fprintf(fout,'%f,%f\n',PBC_young,poisson_mod);

            % Write boundary conditions
            fprintf(fout,'** BOUNDARY CONDITIONS\n');
            fprintf(fout,'** Name: PBC_bottom Type: Symmetry/Antisymmetry/Encastre\n');
            fprintf(fout,'*Boundary\n');
            fprintf(fout,'PBC_bottom_nodes, YSYMM\n');

            fprintf(fout,'** Name: PBC_top Type: Symmetry/Antisymmetry/Encastre\n');
            fprintf(fout,'*Boundary\n');
            fprintf(fout,'PBC_top_nodes, YSYMM\n');

            fprintf(fout,'** Name: PBC_left Type: Symmetry/Antisymmetry/Encastre\n');
            fprintf(fout,'*Boundary\n');
            fprintf(fout,'PBC_left_nodes, XSYMM\n');

            fprintf(fout,'** Name: PBC_right Type: Symmetry/Antisymmetry/Encastre\n');
            fprintf(fout,'*Boundary\n');
            fprintf(fout,'PBC_right_nodes, XSYMM\n');




            fprintf(fout,'** Name: x_constrained_bottom Type: Symmetry/Antisymmetry/Encastre\n');
            fprintf(fout,'*Boundary\n');
            fprintf(fout,'PBC_bottom_nodes, XSYMM\n');

            fprintf(fout,'** Name: y_constrained_bottom Type: Symmetry/Antisymmetry/Encastre\n');
            fprintf(fout,'*Boundary\n');
            fprintf(fout,'PBC_left_nodes, YSYMM\n');


            fclose(fout);
            


       
        end






        
        
    end




end