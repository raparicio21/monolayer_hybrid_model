classdef inputs
   properties (Constant)

       version = 'v11_18_08_2023';


       % ..................... General parameters .........................
       general_time_steps                   = 200       ;
       general_vtk_reduction                = 1         ;     % Reduction of output VTK files
       general_tolerance                    = 1e-6      ;     % Tolerance for non-exact geometric operations
       general_time_step_duration           = 10        ;     % Time step [min]
       general_run_cluster                  = true      ;     % run code on cluster true/false


       % ..................... Monolayer initialization ...................
       monolayer_size_x                     = 450       ;     % [µm]
       monolayer_size_y                     = 450       ;     % [µm]

       monolayer_distribution_mode          = 'fill'    ;     % fill / lattice / random
       monolayer_distribution_cells         = 5000      ;     % random and lattice
       monolayer_distribution_hexagonal     = true      ;     % only fill mode
       monolayer_distribution_variability   = 0.005     ;     % only lattice mode

       monolayer_shake_active               = true      ;
       monolayer_shake_magnitude            = 0.8       ;

       monolayer_hole_active                = false     ;
       monolayer_hole_mode                  = 'interior';     % interior / exterior
       monolayer_hole_size_x                = [100  150];     % [µm]
       monolayer_hole_size_y                = [100  150];     % [µm]

       monolayer_infection_active           = false     ;
       monolayer_infection_mode             = 'interior';     % interior / exterior
       monolayer_inection_radius            = 30        ;     % Fixed infection radius [µm] (when there is no propagation)
       monolayer_infection_propagation      = false     ;     % Propagate infection
       monolayer_infection_step_propagation = 20        ;     % Apply propagation after this number of steps (reach an equilibrium first)
       monolayer_infection_propg_constant   = 8         ;     % [µm^2/min] 1.25 before

       monolayer_extrusion_active           = false     ;
       monolayer_extrusion_initial_step     = 21        ;


       % ..................... Substrate parameters .......................
       substrate_size_x                     = 450       ;     % [µm]
       substrate_size_y                     = 450       ;     % [µm]

       substrate_grid_active                = true      ;
       substrate_grid_size                  = 50        ;     % [µm]

       substrate_properties_E               = 3000      ;     % [Pa]
       substrate_properties_poisson         = 0.49      ;     % [-]
       % Stiffness gradient: f(x) = E
       substrate_properties_gradient_active = false     ;
       substrate_properties_gradient_m      = 0         ;
       substrate_properties_gradient_c      = 1         ;


       % ..................... Cell geometry ..............................
       geometry_radius                      = 6.5       ;     % [µm]
       % Alpha shapes parameter
       geometry_alpha                       = 120       ;

       geometry_shape_resolution            = 8         ;
       geometry_shape_orientation           = 'random'  ;     % random / fixed
       geometry_shape_mode                  = 'circle'  ;     % circle / sinusoid /random

       geometry_shape_sinusoid_frequency    = 4         ;
       geometry_shape_sinusoid_amplitude    = 1.0       ;

       geometry_shape_random_range          = [0.5 2]   ;


       % ..................... Cell life cycle ............................
       % Not implemented here


       % ..................... Cell forces (ABM) ..........................

       % _________________ (1) Cell-Cell interaction forces _______________
       forces_interaction_mode              = 'lennard' ;     % lennard / custom
       forces_interaction_limit             = 3000      ;     % [pN]    before:  1000 pN

       forces_interaction_polarity_active   = false     ;     % effect of cell polarity
       forces_interaction_polarity_factor   = -0.5      ;

       forces_interaction_lennard_constant  = 15        ;     % [pN]    
       forces_interaction_infection         = false     ;     % Differenciate infected vs uninfected
       forces_interaction_eps_uninf         = 30        ;     % Epsilon constant of uninfected cells
       forces_interaction_eps_inf           = 30        ;     % Epsilon constant of infected cells

       % _________________ (2) Cell contraction forces ____________________
       forces_contraction_active            = true      ;
       forces_contraction_mode              = 3         ;     % Contraction mode:
                                                              % 1.- Default contraction (constant)
                                                              % 2.- Random  contraction
                                                              % 3.- Contraction taking into account principal stress
       forces_contraction_limit             = 3000      ;     % [pN] Only for contraction mode = 3
       force_contraction_scale              = 0.2       ;     % Scale factor, only for contraction mode = 3

       % Set contraction magnitude (only cases 1.- default and 2.- random contraction)
       forces_contraction_magnitude         = 30        ;     % [pN]    
       forces_contraction_infection         = false     ;     % Differenciate infected vs uninfected
       forces_contraction_magn_uninf        = 30        ;     % Contraction magnitude of uninfected cells
       forces_contraction_magn_inf          = 30        ;     % Contraction magnitude of infected cells

       % __________________ (3) Cell protrusive forces ____________________
       forces_protrusion_active             = true      ;
       forces_protrusion_magnitude          = 30        ;     % [pN] Only first steps  
       forces_protrusion_limit              = 7000      ;     % [pN] 
       forces_protrusion_scale              = 1.5       ;     % Scale factor
       
       % ____________ (4) Cell infection forces (polarization) ____________
       %           (at the border of infected vs uninfected cells)
       forces_crawling_infection_active     = false     ;
       forces_infection_magnitude           = 3000      ;     % [pN]   

       % _____________________ (5) Cell random walk _______________________
       forces_random_walk_active            = false     ;
       forces_random_walk_generator         = 'default' ;     % default/shuffle  defaults means same rand distribution
       forces_random_walk_magnitude         = 50        ;     % [pN]   
       forces_random_walk_max_persistance   = 20        ;     % [steps] 
       forces_random_walk_infected_cells    = true      ;     % false = abolish random walk on infected cells


       % ..................... FEM paremeters .............................
       fem_properties_E_infected            = 250       ;     % [Pa]    
       fem_properties_E_uninfected          = 500       ;     % [Pa]    
       fem_properties_poisson_ratio         = 0.49      ;     % [-]

       % Boundary conditions, type:
       %        - fixed           : encastre all borders
       %        - symmetric       : XY symmetry
       %        - free            : free borders
       %        - PBC             : Periodic Boundary Conditions
       %        - springs         : springs on the border
       %        - material_around : another material around the monolayer
       % PDE matlab only fixed recommended, symmetric gives problems in the
       % monolayer vertices (abaqus is fine)
       fem_BC_type                          = "symmetric";
       fem_geometric_mesh_order             = "linear"   ;    % quadratic/linear
       fem_software                         = "abaqus"   ;    % abaqus/pde_matlab


       % ..................... Output printing parameters..................
       % Print files (agent-based model)
       output_files_vtk_file = true;
       output_files_mat_file = true;
       
%        % Cell geometry
%        output_files_parameters_id                        = true;
%        output_files_parameters_radius                    = true;
%        output_files_parameters_area                      = true;
%        output_files_parameters_polygon_size              = true;
%        output_files_parameters_perimeter                 = true;
%        output_files_parameters_location                  = true;
% 
%        % Cell kinematics
%        output_files_parameters_initial_position          = true;
%        output_files_parameters_current_position          = true;
%        output_files_parameters_displacement_vector       = true;
%        output_files_parameters_velocity_vector           = true;
% 
%        % Cell total kinematics
%        output_files_parameters_total_distance            = true;
%        output_files_parameters_total_displacement        = true; 
%        output_files_parameters_total_displacement_vector = true;
% 
%        % Cell forces       
%        output_files_parameters_forces                    = true;
%        output_files_parameters_cell_cell_int_forces      = true;
%        output_files_parameters_cell_cont_cent_forces     = true;
       
   end

   methods (Static)
       % Static methods are associated with a class, but not with specific 
       % instances of that class. These methods do not require an object of
       % the class as an input argument. Therefore, you can call static 
       % methods without creating an object of the class.


       function check_variables()
           
           % Cell circumference resolution
           while inputs.geometry_shape_resolution < 3
               error('Error: cell circumference resolution');
           end

           % Check dimensions
           while (inputs.monolayer_size_x > inputs.substrate_size_x || ... 
                  inputs.monolayer_size_y > inputs.substrate_size_y )
               error('Error: monolayer is bigger than substrate!');               
           end

           % Check infection size
           while (inputs.monolayer_size_x/2 < inputs.monolayer_inection_radius || ... 
                  inputs.monolayer_size_y/2 < inputs.monolayer_inection_radius )
               error('Error: the infection domain is bigger than the monolayer!');               
           end

           % Check distribution mode
           while (~strcmp(inputs.monolayer_distribution_mode,'random')  && ...
                  ~strcmp(inputs.monolayer_distribution_mode,'lattice') && ...
                  ~strcmp(inputs.monolayer_distribution_mode,'fill'))
               error('Error: unsupported distribution mode!');            
           end

           % Check background grid dimensions
           if inputs.substrate_grid_active == 1
               while (mod(inputs.substrate_size_x,inputs.substrate_grid_size)~=0) || ...
                     (mod(inputs.substrate_size_y,inputs.substrate_grid_size)~=0) 
                   error('Error: Scaffold size must be divisible by grid size!');
               end
           end

           % Check external stiffness  E(x) = m*x + c
           external_stiffness = inputs.substrate_properties_gradient_m * ...
               inputs.substrate_size_x + inputs.substrate_properties_gradient_c;
           while external_stiffness < 1e-6
               error('Error: Negative external stiffness!');
           end

           % Check cell shape modes
           while (~strcmp(inputs.geometry_shape_mode,'circle')   && ...
                  ~strcmp(inputs.geometry_shape_mode,'sinusoid') && ...
                  ~strcmp(inputs.geometry_shape_mode,'random'))
               error('Error: unsupported shape mode!');            
           end

           while (~strcmp(inputs.geometry_shape_orientation,'random')   && ...
                  ~strcmp(inputs.geometry_shape_orientation,'fixed'))
               error('Error: unsupported orientation mode!');            
           end

           % Missing: Check interaction mode


       end






       % Function useful to apply PBC when the equilibrium is achieved
       function boundary_condition = getBoundaryCondition(iterator)

           % If BC = PBC, apply it after 50 steps, first steps BC =
           % symmetric
           
           if strcmp(inputs.fem_BC_type,'PBC') || strcmp(inputs.fem_BC_type,'material_around')
               if iterator < 50
                    boundary_condition = 'symmetric';         % Set alternate value
                else
                    boundary_condition = inputs.fem_BC_type;  % Get the default value
               end
           end

           % For the the rest of conditions, do not do anything
           if strcmp(inputs.fem_BC_type,'fixed')
               boundary_condition = 'fixed';
           elseif strcmp(inputs.fem_BC_type,'symmetric')
               boundary_condition = 'symmetric';
           elseif strcmp(inputs.fem_BC_type,'free')
               boundary_condition = 'free';
           elseif strcmp(inputs.fem_BC_type,'springs')
               boundary_condition = 'springs';
%            elseif strcmp(inputs.fem_BC_type,'material_around')
%                boundary_condition = 'material_around';
           end

       end

   end
end