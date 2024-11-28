%%%%%%%%%%%%%%%%%%%%%%%% MONOLAYER MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computational model to simulate cell monolayers in two-dimensions (2-D)
% Hybrid model combining:
%                           - Agent Based Models     (ABM)
%                           - Finite Element Methods (FEM)
% By Raul Aparicio Yuste
% PhD student University of Zaragoza & University of Tübingen
% E-mail : raparicio@unizar.es
% Based on the original code written in C++ from Ismael Gonzalez Valverde
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
clc;
close all;

% ..................... Print version......................................
addpath inputs                       ;   % Folder with input parameters      
addpath src                          ;   % Source code
addpath ext_func                     ;   % External functions from others
addpath(fullfile('src','io','file')) ;   % Code related to write files
addpath(fullfile('src','fem'))       ;   % Code related to finite element m.
disp(inputs.version)                 ;   % Version display

if inputs.general_run_cluster == false
       
    % ..................... Create directories.................................
    % Create subdirectories for results:
    %   - cells       : cell centroids with ABM variables
    %   - fem         : mesh generated with FEM variables
    %   - voronoi     : voronoi distribution (polygons) with ABM variables
    %   - raw         : raw data from ABM for analysis
    %   - fem         : raw data from FEM for analysis
    %   - inp         : input files used during the simulation in Abaqus
    %   - nuclei      : image of nuclei to compute trackmate and coordination
    %   - elipsoids   : maps to plot cell ellipsoids in ParaView (Glyph mode)
    %   - force_vert  : maps of forces at the vertex of the cell 
    %   - alpha_shape : map of alpha shape (monolayer border), not used anymore
    
    rmdir('results', 's');                   % Remove old folder with results
    mkdir('results')     ;                   % Make new directory results
    
    subdirectories = {'cells', 'fem', 'voronoi','raw','raw_fem','inp','nuclei','elipsoids', 'force_vert'}; % 'alpha_shape'
    
    for i = 1:length(subdirectories)
        mkdir(fullfile('results', subdirectories{i}));
    end
    
    clear subdirectories ;                   % Clear variable
    % Copy python file in order to extract variables from odb file (Abaqus):
    copyfile('extract_variables.py', fullfile('results','inp'));
    
    
    
    % ..................... Define parameters .................................
    inputs.check_variables
    disp('Info: Configuring ECM');
    my_ecm       = ecm();
    my_ecm       = my_ecm.setECMSize(inputs.substrate_size_x,inputs.substrate_size_y,0);
    % ECM is not used in this simulation so far
    
    disp('Info: Configuring monolayer');
    my_monolayer = monolayer();
    my_monolayer = my_monolayer.setMonolayerSize(inputs.monolayer_size_x,inputs.monolayer_size_y,0);
    my_monolayer = my_monolayer.initialize_monolayer();
    
    % In case of PBC, we need to create the external mesh
    if strcmp(inputs.fem_BC_type,'PBC')
        disp('Info: Define Periodic Boundary Condition');
        my_mesh      = mesh_fem();
        my_mesh      = my_mesh.period_boundary_mesh(my_monolayer);
    end
    
    my_FEM         = fem_analysis();      % Initialize structure
    text1          = 'Step '       ;
    
    
    
    % ..................... Calculation loop ..................................
    disp('Info: Starting calculation loop');
    for i = 1:inputs.general_time_steps
        %tic
        % Set step
        my_monolayer.step = i;
        disp([text1, num2str(my_monolayer.step), '/', num2str(inputs.general_time_steps)]);
    
        % Life cycle
        %disp('Info: Computing cell life cycle');
        %..............
        %tic
    
        % Generate monolayer geometry
        disp('Info: Generating monolayer geometry');
        my_monolayer   = my_monolayer.generateGeometry();
        %toc
    
        % Generate plots for PIV analysis and TrackMate (Fiji)
        my_monolayer   = my_monolayer.generateNucleiPlots();
    
        %tic
        % Calculate Forces
        disp('Info: Calculating interaction forces');
        my_monolayer   = my_monolayer.calculateForces(my_FEM);
        %toc
    
        %tic
        % Write data to VTK and mat file
        disp('Info: Writing agent-based model results');
        my_data        = vtk_mat_Parser();
        my_data        = my_data.writeVTK_mat_file(my_monolayer,my_FEM);
        %toc
    
        %tic
        % Create a Mesh from monolayer and ECM data
        disp('Info: Generating mesh');
        my_mesh        = mesh_fem();
        my_mesh        = my_mesh.initialize(my_monolayer,my_ecm);
        %toc
    
        %tic
        % Run FE analysis
        disp('Info: Solving FEM');
        my_FEM         = fem_analysis();
        my_FEM         = my_FEM.run(my_mesh,my_monolayer);
        %toc
    
        %tic
        % Write FEM data to VTK and mat file
        disp('Info: Writing FEM results');
        my_FEM_results          = vtk_mat_fem_Parser();
        [my_FEM_results,my_FEM] = my_FEM_results.write_fem_VTK_mat_file(my_FEM,my_monolayer,my_mesh);
        %toc
    
        %tic
        % Pass data to discrete cell model
        disp('Info: Moving cells');
        my_monolayer          = my_monolayer.setDisplacements(my_FEM,my_mesh)    ;
        [my_monolayer,my_FEM] = my_monolayer.computeAverageStress(my_FEM,my_mesh);
        %toc
        %toc
        
    end
    
    
    % % Plot cell centers
    % for i = 1:length(my_monolayer.cell_container)
    % 
    %     x(i)= my_monolayer.cell_container(i).position(1);
    %     y(i)= my_monolayer.cell_container(i).position(2);
    %     figure(1)
    %     hold on
    %     plot(x(i), y(i), '.', 'markersize', 8)
    %     hold off
    % 
    % end

elseif inputs.general_run_cluster == true
     
    try    
        
        % ..................... Create directories.................................
        % Create subdirectories for results:
        %   - cells       : cell centroids with ABM variables
        %   - fem         : mesh generated with FEM variables
        %   - voronoi     : voronoi distribution (polygons) with ABM variables
        %   - raw         : raw data from ABM for analysis
        %   - fem         : raw data from FEM for analysis
        %   - inp         : input files used during the simulation in Abaqus
        %   - nuclei      : image of nuclei to compute trackmate and coordination
        %   - elipsoids   : maps to plot cell ellipsoids in ParaView (Glyph mode)
        %   - force_vert  : maps of forces at the vertex of the cell 
        %   - alpha_shape : map of alpha shape (monolayer border), not used anymore
        
        rmdir('results', 's');                   % Remove old folder with results
        mkdir('results')     ;                   % Make new directory results
        
        subdirectories = {'cells', 'fem', 'voronoi','raw','raw_fem','inp','nuclei','elipsoids', 'force_vert'}; % 'alpha_shape'
        
        for i = 1:length(subdirectories)
            mkdir(fullfile('results', subdirectories{i}));
        end
        
        clear subdirectories ;                   % Clear variable
        % Copy python file in order to extract variables from odb file (Abaqus):
        copyfile('extract_variables.py', fullfile('results','inp'));
        
        
        
        % ..................... Define parameters .................................
        inputs.check_variables
        disp('Info: Configuring ECM');
        my_ecm       = ecm();
        my_ecm       = my_ecm.setECMSize(inputs.substrate_size_x,inputs.substrate_size_y,0);
        % ECM is not used in this simulation so far
        
        disp('Info: Configuring monolayer');
        my_monolayer = monolayer();
        my_monolayer = my_monolayer.setMonolayerSize(inputs.monolayer_size_x,inputs.monolayer_size_y,0);
        my_monolayer = my_monolayer.initialize_monolayer();
        
        % In case of PBC, we need to create the external mesh
        if strcmp(inputs.fem_BC_type,'PBC')
            disp('Info: Define Periodic Boundary Condition');
            my_mesh      = mesh_fem();
            my_mesh      = my_mesh.period_boundary_mesh(my_monolayer);
        end
        
        my_FEM         = fem_analysis();      % Initialize structure
        text1          = 'Step '       ;
        
        
        
        % ..................... Calculation loop ..................................
        disp('Info: Starting calculation loop');
        for i = 1:inputs.general_time_steps
            %tic
            % Set step
            my_monolayer.step = i;
            disp([text1, num2str(my_monolayer.step), '/', num2str(inputs.general_time_steps)]);
        
            
            % Life cycle
            %disp('Info: Computing cell life cycle');
            %..............
            %tic
        
            % Generate monolayer geometry
            disp('Info: Generating monolayer geometry');
            my_monolayer   = my_monolayer.generateGeometry();
            %toc
        
            % Generate plots for PIV analysis and TrackMate (Fiji)
            my_monolayer   = my_monolayer.generateNucleiPlots();
        
            %tic
            % Calculate Forces
            disp('Info: Calculating interaction forces');
            my_monolayer   = my_monolayer.calculateForces(my_FEM);
            %toc
        
            %tic
            % Write data to VTK and mat file
            disp('Info: Writing agent-based model results');
            my_data        = vtk_mat_Parser();
            my_data        = my_data.writeVTK_mat_file(my_monolayer,my_FEM);
            %toc
        
            %tic
            % Create a Mesh from monolayer and ECM data
            disp('Info: Generating mesh');
            my_mesh        = mesh_fem();
            my_mesh        = my_mesh.initialize(my_monolayer,my_ecm);
            %toc
        
            %tic
            % Run FE analysis
            disp('Info: Solving FEM');
            my_FEM         = fem_analysis();
            my_FEM         = my_FEM.run(my_mesh,my_monolayer);
            %toc
        
            %tic
            % Write FEM data to VTK and mat file
            disp('Info: Writing FEM results');
            my_FEM_results          = vtk_mat_fem_Parser();
            [my_FEM_results,my_FEM] = my_FEM_results.write_fem_VTK_mat_file(my_FEM,my_monolayer,my_mesh);
            %toc
        
            %tic
            % Pass data to discrete cell model
            disp('Info: Moving cells');
            my_monolayer          = my_monolayer.setDisplacements(my_FEM,my_mesh)    ;
            [my_monolayer,my_FEM] = my_monolayer.computeAverageStress(my_FEM,my_mesh);
            %toc
            %toc
            
        end
        
        
        % % Plot cell centers
        % for i = 1:length(my_monolayer.cell_container)
        % 
        %     x(i)= my_monolayer.cell_container(i).position(1);
        %     y(i)= my_monolayer.cell_container(i).position(2);
        %     figure(1)
        %     hold on
        %     plot(x(i), y(i), '.', 'markersize', 8)
        %     hold off
        % 
        % end
        exit(0);
    catch exception

        % Capture the error message and stack trace
        error_message = getReport(exception, 'extended', 'hyperlinks', 'on');

        % Display the error message
        disp(['Error: ' error_message]);

        % Save the error message to a file for later inspection
        error_log_file = fopen('error_log.txt', 'w');
        fprintf(error_log_file, '%s', error_message);
        fclose(error_log_file);

        exit(1);
    
    end
end




