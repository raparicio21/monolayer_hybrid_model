classdef fem_analysis

    properties
        pde_model             % partial diferential equation model/object
        results_pde           % fem pde toolbox results
        displ_results         % displacements solutions at nodes from abaqus
        stress_results        % stress solutions at integration points from abaqus
        average_cell_stress   % average principal stresses per cell 
                              % Column 1 mean s_I  per cell within all range
                              % Column 2 mean s_II per cell within all range
                              % Column 3 mean s_I  per cell within the range [0 s_Imax]
                              % Column 4 mean s_II per cell within the range [s_II_max 0]
    end

    methods

        %................... Object initialization ........................
        function obj = fem_analysis()
            addpath(fullfile('inputs')); 
            addpath(fullfile('src'));    
            addpath(fullfile('src','fem'));  
        end




        %................... FEM analysis .................................
        function FEM_obj = run(FEM_obj,mesh_obj,monolayer_obj)


            %..................... PREPROCESS .............................

            % Define geometry
            %disp('geometry')
            %tic
            %FEM_obj = FEM_obj.defineGeometry(mesh_obj);
            %toc

            if strcmp(inputs.fem_software,'pde_matlab')

                % Define material properties
                %disp('material')
                %tic
                FEM_obj = FEM_obj.defineMaterial();
                %toc
    
                % Define boundary conditions
                %disp('BC')
                %tic
                FEM_obj = FEM_obj.defineBC();
                %toc
    
                % Define quadratic mesh
                %disp('mesh')
                %tic
                FEM_obj = FEM_obj.defineMesh();
                %toc
    
                % Introduce concentrated force
                %disp('force')
                %tic
                FEM_obj = FEM_obj.defineConcentratedForce(monolayer_obj);
                %toc

            elseif strcmp(inputs.fem_software,'abaqus')

                %tic                
                % When running Abaqus == true, create file.inp
                %disp('write inp')
                my_inp = inp_file();
                my_inp.write_inp(FEM_obj.pde_model,monolayer_obj,mesh_obj);
                %toc
            end

            %..................... ANALYSIS ..............................
            if strcmp(inputs.fem_software,'pde_matlab')
                % Solve the FEM
                %disp('solve')
                %tic
                structuralresults = solve(FEM_obj.pde_model);
                %toc
                  
                % Update results
                %disp('update')
                %tic
                FEM_obj.results_pde   = structuralresults;
                %toc

            elseif strcmp(inputs.fem_software,'abaqus')

                %tic
                % Set directory tu run .inp model
                current_directory = pwd;
                cd('./results/inp/')   ;
                
                % Run abaqus
                step          = monolayer_obj.step                                         ; 
                command_1     = sprintf('abaqus job=inp_%05i.inp cpus=1 interactive', step);
                system(command_1)                                                          ;
                %toc

                %tic
                % Extract output variables on .dat file through python
                % library (odbAccess)                
                odb_name      = sprintf('abaqus viewer noGUI=extract_variables.py -- inp_%05i.odb %i', step,step);
                %odb_name      = str(odb_name);
                system(odb_name);
                %toc

                %tic
                % Read output files (displacement and stress)
                displ_file   = sprintf('displacements_%05i.dat', step);
                data_displ   = load(displ_file)                       ;

                stress_file  = sprintf('stress_%05i.dat', step)       ;
                data_stress  = load(stress_file)                      ;
                %toc

                %tic
                % Column 1 : Node_id 
                % Column 2 : Ux 
                % Column 3 : Uy 
                FEM_obj.displ_results    = data_displ;
                % Column 1 : Integration_point_coord_x 
                % Column 2 : Integration_point_coord_y
                % Column 3 : S_x 
                % Column 4 : S_y 
                % Column 5 : S_xy 
                % Column 6 : S_maxprin from abaqus (range [0 s_I_max] )
                % Column 7 : S_minprin from abaqus (range [s_II_max 0])
                % Column 8 : S_maxprin from abaqus (all range         ) (added later)
                % Column 9 : S_minprin from abaqus (all range         ) (added later)
                FEM_obj.stress_results   = data_stress;

                % Abaqus gives you s_I within the range [0 s_Imax] and s_II
                % within the range [s_IImax 0]. To retrieve all the principal values,
                % we add two more columns we will compute later on
                column_zeros             = zeros(size(data_stress,1),1)                      ;
                FEM_obj.stress_results   = [FEM_obj.stress_results column_zeros column_zeros];
                
                %toc


                %tic                
                % Remove some output files                
                filename     = cell(5, 1)                             ;
                filename{1}  = sprintf('inp_%05i.com'          , step);
                filename{2}  = sprintf('inp_%05i.msg'          , step);
                filename{3}  = sprintf('inp_%05i.prt'          , step);
                filename{4}  = sprintf('inp_%05i.sim'          , step);
                filename{5}  = sprintf('inp_%05i.sta'          , step);
                filename{6}  = sprintf('inp_%05i.odb'          , step);
                filename{7}  = sprintf('inp_%05i.dat'          , step);
                filename{8}  = sprintf('displacements_%05i.dat', step);
                filename{9}  = sprintf('stress_%05i.dat'       , step);
                filename{10} = sprintf('abaqus.rpy')                  ;
                delete(filename{:})                                   ;

                cd(current_directory)                                 ;
                %toc

                
            end


            %..................... POSTPROCESS ............................

            % Manual postprocess while implementing the code
            %FEM_obj           = FEM_obj.manualPostprocess();
       
        end





        %................... Model Geometry ...............................
        function FEM_obj = defineGeometry(FEM_obj,mesh_obj)

            % Initialize pde model            
            model     = createpde("structural","static-planestress");
            TR        = mesh_obj.triangulation_obj;
            tnodes    = TR.Points';
            telements = TR.ConnectivityList';
            elem_id   = (1:numel(telements(1,:)));
            geometryFromMesh(model,tnodes,telements,elem_id);

            % Update model
            FEM_obj.pde_model = model;

%             figure
%             pdeplot(model)
%             
%             figure
%             pdegplot(model,"FaceLabels","on") 
% 
%             figure
%             pdegplot(model,"EdgeLabels","on") 
%             
%             figure
%             pdemesh(model.Mesh,"NodeLabels","on")
%             hold on
%             plot(model.Mesh.Nodes(1,:), ...
%                  model.Mesh.Nodes(2,:), ...
%                  "ok","MarkerFaceColor","g")  
%             hold off
% 
%             figure
%             pdegplot(model,"VertexLabels","on")        
        end




        %................... Material properties ..........................
        function FEM_obj = defineMaterial(FEM_obj)

            youngs_modulus = inputs.fem_properties_E_uninfected; % MPa
            poissonsratio  = inputs.fem_properties_poisson_ratio; 

            % Material properties
            structuralProperties(FEM_obj.pde_model,              ...
                                 "YoungsModulus",youngs_modulus, ... 
                                 "PoissonsRatio",poissonsratio);
      
        end





        %................... Boundary Conditions ..........................
        function FEM_obj = defineBC(FEM_obj)

            % To compute freeBoundaries, the triangulation is needed

            % The triangulation of the model cannot be accesed directly, so
            % one needs to compute the  triangulation given the new
            % organization of the vertices. Careful->

            % While loading the model with GeometryFromMesh, the vertices
            % of the model are not organised as the initial configuration.
            % Vertex_ID =! Mesh_nodes vertices

            % Translate mesh vertices/nodes to the new vertex list:


            % Load model
            model            = FEM_obj.pde_model;

            % Load geometry and mesh data from model
            model_vertices   = model.Geometry.Vertices(:,1:2); % [X Y]
            model_nodes      = model.Mesh.Nodes';
            model_elements   = model.Mesh.Elements; % mesh connectivity

            % Parse vertex information:  
            % find vertices from mesh_node in the new list of vertices of 
            % the model
            [~,loc_id] = ismember(model_nodes,model_vertices,'rows');
            % Update connectivity list
            model_elements=loc_id(model_elements);

            % Create triangulation
            TR_model = triangulation(model_elements',model_vertices);
%             figure
%             hold on
%             triplot(TR_model)
%             hold off

            % Find freeBoundary vertices
            free_edges = freeBoundary(TR_model);
%             hold on
%             plot(model_vertices(free_edges,1),model_vertices(free_edges,2),'-r','LineWidth',2)
%             hold off
    
            % Get boundary edges list
            Boundary_Edges = [];
            for i = 1:numel(free_edges(:,1))
                edges_vertices = model_vertices(free_edges(i,:),:);
                midpoint       = (edges_vertices(1,:) + edges_vertices(2,:)) / 2;
                EdgeID         = nearestEdge(model.Geometry,midpoint);
                Boundary_Edges = [Boundary_Edges, EdgeID];
            
            end

%             figure
%             pdegplot(model,"EdgeLabels","on")

            % Apply BC to the model
            type_BC = inputs.fem_BC_type;
            structuralBC(FEM_obj.pde_model,"Edge",Boundary_Edges,...
                                             "Constraint", type_BC); 
      
        end


        %................... Mesh parameters ..............................
        function FEM_obj = defineMesh(FEM_obj)

            % The order of aproximation of the initial mesh is linear. In
            % order to generate the same mesh but quadratic, one sets the
            % Hmin of the element high so all the elements will remain as
            % the initial mesh

            geometric_order = inputs.fem_geometric_mesh_order;

            new_mesh = generateMesh(FEM_obj.pde_model,"Hmin",100, ...
                               "GeometricOrder",geometric_order);

%             pdemesh(new_mesh,"NodeLabels","on")
%             hold on
%             plot(new_mesh.Nodes(1,:), ...
%                  new_mesh.Nodes(2,:), ...
%                  "ok","MarkerFaceColor","g")  
%             hold off

            % Evaluate shape quality of mesh elements
%             Q       = meshQuality(new_mesh);
%             elemIDs = find(Q < 0.2);
% 
%             pdemesh(new_mesh,"FaceAlpha",0.5)
%             hold on
%             pdemesh(new_mesh.Nodes,new_mesh.Elements(:,elemIDs), ...
%                               "EdgeColor","blue")
%             hold off
% 
%             figure
%             hist(Q)
%             xlabel("Element Shape Quality","fontweight","b")
%             ylabel("Number of Elements","fontweight","b")
% 
%             Qworst  = min(Q)
%             elemIDs = find(Q==Qworst)
            

        end


        %................... Concentrated forces ..........................
        function FEM_obj = defineConcentratedForce(FEM_obj,monolayer)

            % Introduce concentrated forces on vertex. On nodes is not
            % possible, that is the reason one eneds to load all the nodes
            % as vertices in your model

            % All forces are defined as the initial triangulation / 
            % configuration, but the model has changed the order of
            % vertices. One needs to parse the forces from the initial
            % configuration to the current one.

            % Load forces from monolayer object

            % Centroid forces
            list_centroids = monolayer.monolayer_cell_points;
            centroid_nodes = (1:numel(monolayer.cell_id))';
            %cent_force_map = zeros(numel(monolayer.cell_id),3);

            cent_force_map = monolayer.force_list;

            % Vertices forces
            list_vertices  = monolayer.voronoi_vertices;
            global_list    = [list_centroids; list_vertices];

            vertex_nodes   = (numel(monolayer.cell_id) + 1: ...
                              numel(monolayer.cell_id) + ...
                              numel(list_vertices(:,1)))';
            vert_force_map = monolayer.vertices_forces;


            % Assembly force maps (according to initial mesh)
            id_nodes  = [centroid_nodes; vertex_nodes];
            force_map = [cent_force_map; vert_force_map];

            % Add new vertices 
            %loadedVertex = addVertex(model.Geometry,"Coordinates",global_list);

            % Parse cells to the model according to the new mesh 
            % organization
            B = global_list;                                 % initial mesh/configuration
            A = FEM_obj.pde_model.Geometry.Vertices(:,1:2);  % new     mesh/configuration
            [~,Locb] = ismember(A,B,'rows');

            % Order new force_map according to the new mesh
            new_force_map = force_map(Locb,:);

            % Add forces to the model (to the vertices)
            for i = 1:numel(id_nodes)
                structuralBoundaryLoad(FEM_obj.pde_model, ...
                       "Vertex", id_nodes(i)            , ...
                        "Force", new_force_map(i,1:2))  ;
            end

        end



        %................... Manual postprocess ...........................
        function FEM_obj = manualPostprocess(FEM_obj)

            model             = FEM_obj.pde_model;
            structuralresults = FEM_obj.results_pde;

            %nodes_prueba = model.Mesh.findNodes('region','Vertex',5);

            figure
            pdeplot(model, "Deformation",structuralresults.Displacement)

            figure
            u_magn = structuralresults.Stress.xx;
            pdeplot(model,"XYData",u_magn,"ColorMap","jet")

            figure
            u_magn = structuralresults.Stress.yy;
            pdeplot(model,"XYData",u_magn,"ColorMap","jet")

            figure
            u_magn = structuralresults.Stress.xy;
            pdeplot(model,"XYData",u_magn,"ColorMap","jet")
            
            figure
            u_magn = structuralresults.Displacement.Magnitude;
            pdeplot(model,"XYData",u_magn,"ColorMap","jet")
            
            figure
            u = structuralresults.Displacement.x;
            pdeplot(model,"XYData",u,"ColorMap","jet")
            
            figure
            u = structuralresults.Displacement.y;
            pdeplot(model,"XYData",u,"ColorMap","jet")
            
            
            
            figure
            u = structuralresults.NodalSolution;
            pdeplot(model,"XYData",u,"ColorMap","jet")
            
            figure
            pdeplot(model,"XYData",u,"ZData",u,"ColorMap","jet")
            
            figure
            ux = structuralresults.XGradients;
            uy = structuralresults.YGradients;
            pdeplot(model,"FlowData",[ux,uy])
            
            figure
            u = structuralresults.NodalSolution;
            ux = structuralresults.XGradients;
            uy = structuralresults.YGradients;
            h = pdeplot(model,"XYData",u,"ZData",u, ...
                              "FaceAlpha",0.5, ...
                              "FlowData",[ux,uy], ...
                              "ColorMap","jet", ...
                              "Mesh","on");
            
            % Plot the deformed shape using the default scale factor. By default, 
            % pdeplot internally determines the scale factor based on the dimensions 
            % of the geometry and the magnitude of deformation.
            
            figure
            pdeplot(model, ...
                    "XYData",structuralresults.VonMisesStress, ...
                    "Deformation",structuralresults.Displacement, ...
                    "ColorMap","jet")
            
            % Plot the deformed shape with the scale factor 500.
            figure
            pdeplot(model, ...
                    "XYData",structuralresults.VonMisesStress, ...
                    "Deformation",structuralresults.Displacement, ...
                    "DeformationScaleFactor",0.1,...
                    "ColorMap","jet")
            
            % Plot the deformed shape without scaling.
            figure
            pdeplot(model,"XYData",structuralresults.VonMisesStress, ...
                                    "ColorMap","jet")
            
            % To plot the kth component of a solution to a PDE system, extract the 
            % relevant part of the solution. For example, when using a PDEModel object,
            % specify:
            figure
            u = structuralresults.NodalSolution; % each column of u has one component of u
            pdeplot(model,"XYData",u(:,1)) % data for column k

    
      
        end




        
    end




end