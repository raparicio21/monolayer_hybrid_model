%% Correlation length simulations
clear all;
clc;
close all;

addpath 'analysis'

% Load first step to retrieve some parameters
load ./results/raw/data_00001.mat

case_name                      = '1_high_conf_fI_1_5';
make_video                     = true;
make_figure                    = true;
fcal                           = 0.28;        % um/pixel  (SIMULATION VALUE!)
monolayer_length               = 450;      % um
% grid to interpolate displacements from nuclei of the simulation:
[xq,yq]                        = meshgrid(0:6.9:monolayer_length); 
total_time_steps               = 200;
delta_time                     = 20; % minutes
sim_time_min                   = 0:delta_time:total_time_steps*delta_time-delta_time; % minutes
sim_time_hou                   = sim_time_min/60; % hours
number_of_cells                = numel(id_out);
Imaging_minutes_post_infection = 0;

% Cell kinematics
position_x                     = zeros(number_of_cells,total_time_steps);
position_y                     = zeros(number_of_cells,total_time_steps);
displacement_time_x            = zeros(number_of_cells,total_time_steps);
displacement_time_y            = zeros(number_of_cells,total_time_steps);
% displacement_time_magn       = zeros(number_of_cells,total_time_steps);

% Load parameters
for i = 1:total_time_steps

    file_name = sprintf('./results/raw/data_%05i.mat', i); 
    load(file_name)

    % Cell kinematics
    position_x                   (:,i) = curr_pos_out(:,1)                   ;
    position_y                   (:,i) = curr_pos_out(:,2)                   ;
    displacement_time_x          (:,i) = disp_vec_out(:,1)                   ;
    displacement_time_y          (:,i) = disp_vec_out(:,2)                   ;
    % displacement_time_magn       (:,i) = (displacement_time_x(:,i).^2  + ...
    %                                       displacement_time_y(:,i).^2).^(0.5);

end


clear area_out pol_sizes_out perimeter_out location_out disp_vec_out   ...
      velo_vec_out tot_dist_out tot_disp_out t_dis_ve_out id_out       ...
      init_pos_out radius_out curr_pos_out tot_force_out int_force_out ...
      con_force_out rand_walk_out pro_force_out inf_crawl_out          ...
      orientation_min_s_I_vect orientation_min_s_II_vect               ...
      surround_out surround_2nd_out surround_3rd_out area_change       ...
      cellshape_out infected_out weighted_force_I weighted_force_II;


% Interpolate and create a matrix for the displacements (they are scattered
% points right now, a cloud of points, and we want to get a matrix of
% displacements to compute the correlation length)
matrix_u     =  cell(1,total_time_steps);
matrix_v     =  cell(1,total_time_steps);

for i = 1:total_time_steps

    x   = position_x(:,i)           ;
    y   = position_y(:,i)           ;
    u_x = displacement_time_x(:,i)  ;
    u_y = displacement_time_y(:,i)  ;

    % Create a scatteredInterpolant
    F   = scatteredInterpolant(x, y, u_x, 'nearest', 'nearest');
      
    u_interp = F(xq, yq);   
    % % Plot the original points
    % scatter(x, y, 30, u_x, 'filled');
    % hold on;
    % colorbar;
    % % Plot the interpolated mesh
    % mesh(xq, yq, u_interp, 'EdgeColor', 'interp', 'FaceColor', 'none');    
    % xlabel('X (microns)');
    % ylabel('Y (microns)');
    % zlabel('displ_x');
    % title('Interpolated displacement Mesh');  
    % hold off;

    % Create a scatteredInterpolant
    F   = scatteredInterpolant(x, y, u_y, 'nearest', 'nearest');   
    v_interp = F(xq, yq);
    % % Plot the original points
    % scatter(x, y, 30, u_y, 'filled');
    % hold on;
    % colorbar;
    % % Plot the interpolated mesh
    % mesh(xq, yq, v_interp, 'EdgeColor', 'interp', 'FaceColor', 'none');    
    % xlabel('X (microns)');
    % ylabel('Y (microns)');
    % zlabel('displ_y');
    % title('Interpolated displacement Mesh');  
    % hold off;

    u_interp_1 = u_interp(5:end-5,5:end-5);
    v_interp_1 = v_interp(5:end-5,5:end-5);

    % Save values
    matrix_u{i} = u_interp_1;
    matrix_v{i} = v_interp_1;

end

%%

set(0, 'DefaultFigureVisible', 'off') %'on' or 'off' to make figures pop up or not
%ScaleArrow = 5;
%overlap    = 24;
overlap    = 24; % Simulation, we assume overlap = 1
Fontsize   = 20;

% Script to plot vectors on nuclei AND calculate correlation length analysis (kymograph)
kin=2;       kfin = total_time_steps ;  % First and last frames
speed     = zeros(1,total_time_steps);  % Prepare array for speed calculation
R_max     = 25                       ;  % Maximum radius (pixels)
kymograph = zeros(R_max+1,kfin)      ;  % 26 because of R_max=25 and R_step=1



for k= kin:kfin  

     disp(['k = ', num2str(k)]);
    % Load displacements
    u = matrix_u{k};
    v = matrix_v{k};

    speed(1,k)= mean(mean(sqrt((u*fcal).^2+(v*fcal).^2)))/delta_time*60; %um/h 
    %%%%%%%%%%%%%%%%%%% Correlacion Angelini et al. (v1) %%%%%%%%%%%%%%%%%%
    % Subtraction of mean
    mean_u = mean(mean(u));
    mean_v = mean(mean(v));
    u      = u-mean_u     ;
    v      = v-mean_v     ;
    
    % Generate matrix with distances
    Generar_distancias_i = 1:1:length(u)                                              ;
    Generar_distancias_i = repmat(Generar_distancias_i,length(Generar_distancias_i),1);
    Generar_distancias_j = (1:1:length(u))'                                           ;
    Generar_distancias_j = repmat(Generar_distancias_j,1,length(Generar_distancias_j));
    
    %Initilization of matrices
    P_escalar_superior  = zeros(length(u));
    P_escalar_inferior  = zeros(length(u));
    P_escalar           = []              ;
    numero_R            = 1               ;
    R_ant               = 0               ;
    bordes              = 0               ; % Changing this value we can remove the correlation of the edges. (pixels)
    R_step              = 1               ; % Step of the radius discretization (pixels)
    %R_max               = 34              ; % Maximum radius (pixels). This corresponds to distance R_max*overal*fcal

    for R = 0 : R_step : R_max

        for i = 1+bordes:length(u) - bordes

            for j = 1+bordes:length(u) - bordes

                Distancias                  =  sqrt((i-Generar_distancias_i).^2+(j-Generar_distancias_j).^2);
                [Posiciones_y,Posiciones_x] =  find( Distancias <= R   & Distancias >= R_ant)               ;
                vector_0                    =  [u(i,j) v(i,j)]                                              ;
                
                for escalar = 1:length(Posiciones_x)

                    vector_Posicion        = [u(Posiciones_x(escalar),Posiciones_y(escalar)) v(Posiciones_x(escalar),Posiciones_y(escalar))]';
                    P_escalar(i,j,escalar) = vector_0*vector_Posicion                                       ;
               
                end

                P_escalar_superior(i,j) = mean( P_escalar(i,j,find(P_escalar(i,j,:)~=0 & ~isnan(P_escalar(i,j,:)))));
                P_escalar_inferior(i,j) = vector_0 * vector_0';

            end

        end

        % Correlation for each R iteration for each image(k):
        Coef(numero_R, k)       = sum(sum(P_escalar_superior))  /  sum(sum(P_escalar_inferior)); 
        % Radius for each iteration for each image(k)(the same for all
        % images):
        valor_Radio(numero_R,k) = R         ; 
        numero_R                = numero_R+1;

        R_ant                   = R         ;

        %%%%
        %Image of the correlation for a given radius (R)
        %Be careful (and comment) if R_step and/or k (number of images) is too large!!!
        % Coef_2D=P_escalar_superior./P_escalar_inferior; %Matrix of the individual correlation of each pixel for the radius R
        % % Plot the 2D map of the correlation coefficinet for each R considered
        % figure
        % imagesc([1:length(Coef_2D)]*24*fcal,[1:length(Coef_2D)]*24*fcal, Coef_2D); caxis([-1 1]);colorbar;axis image
        % set(gca,'FontSize',18,'DefaultAxesFontName', 'Arial');title(R);xlabel('um');ylabel('um');

    end
    
    valor_Radio(:,k) = valor_Radio(:,k)*overlap*fcal; % Change of units from pixels to um
    time(k)          = k*delta_time                 ; % time in 'delta' units
    

    % Length from equation
    %Fit equation
    eqexpon        = fittype('exp(-R/R0)', 'independent',{'R'}, 'coefficients',{'R0'})            ;
    myfit          = fit(valor_Radio(:,k),Coef(:,k),eqexpon,'startpoint', 50, 'maxiter',20000);
    R0_equation(k) = myfit.R0*2                                                                   ; %Value of R0 from equation
    
    % Plot Corr Coefficient over radius/distance per timestep
    fig0           = figure;
    plot(myfit,valor_Radio(:,k),Coef(:,k))
    xlabel('Radius (um)','FontSize',Fontsize)
    ylabel('C(Radius)','FontSize',Fontsize) % Correlation coefficient depending on radius
    ylim([-0.2 1])
    title([num2str(round(((k*delta_time+Imaging_minutes_post_infection)/60),2)) ' hpi'], 'FontWeight','normal', 'FontSize',Fontsize)
    %print(fig,'-dtiff',name);
    Save_figure = getframe(fig0); %or Save_figure = getframe(gcf); 
    imwrite(Save_figure.cdata, './CorrCoefTimesteps.tif', 'WriteMode', 'append', 'Compression','none');
    
    
    kymograph(:,k)=Coef(:,k);
    
    % Length from minimum value
    [min_cor, position] = min(Coef(:,k))                 ;
    R0_minimum(k)       = valor_Radio(position,k)        ;
    %close all
    save('./CorrCoefKymograph.mat', 'kymograph')
    Fig2_Correlation    = [time; R0_equation; R0_minimum];
    save('./Correlation.mat','Fig2_Correlation' )
    
    % Save whole workspace
    %save([Dir 'Results_summary/Save_Workspace/CorrLength_Pos' num2str(Pos) '.mat']) 
    
end 
    
%clear -Dir % Clear everything except Dir

Fig2Data    = load('Correlation.mat')       ;
Fig3Data    = load('CorrCoefKymograph.mat') ;
time        = Fig2Data.Fig2_Correlation(1,:);
R0_equation = Fig2Data.Fig2_Correlation(2,:);
%R0_minimum  = Fig2Data.Fig2_Correlation(3,:);
kymograph   = Fig3Data.kymograph;

% Plot of evolution of mean correlation length in each frame
fig2=figure;
hold on
plot((time(2:end)+Imaging_minutes_post_infection)/60,R0_equation(2:end))
%plot(time,R0_minimum)
xlabel('Time (hpi)','FontSize'             ,Fontsize)
ylabel('Correlation length (um)','FontSize',Fontsize)
title('Pos','FontSize',Fontsize, 'FontWeight','normal')
ylim([0 120])
%legend('equation','minimum value')
hold off
Save_figure = getframe(fig2); imwrite(Save_figure.cdata, 'Correlation.tif'); 
%print(fig2,'-dtiff',[Dir 'Urapiv/Pos' num2str(Pos) '/Pos' num2str(Pos) '_Correlation']);

% Kymograph where you can see how the correlation coefficient changes as a
% function of time and space
fig3=figure;
imagesc((time+Imaging_minutes_post_infection)/60,(1:size(kymograph,1))*overlap*fcal, kymograph); colormap redblue;Col=colorbar; Col.Label.String='Correlation coefficient';
clim([-1 1])
xlabel('Time post-infection (h)','FontSize',Fontsize)
ylabel('Radius (um)','FontSize'            ,Fontsize)
title('Pos', 'FontWeight','normal')
set(gca,'FontSize',Fontsize);
Save_figure = getframe(fig3); imwrite(Save_figure.cdata, 'CorrelationCoef.tif'); 
%print(fig3,'-dtiff',[Dir 'Results_summary/Pos' num2str(Pos) '_CorrelationCoef']);

% %         %%% Sort and delete files
% %         DeleteFolder = DirPIVresults; 
% %         status = rmdir(DeleteFolder, 's'); % remove the directory DeleteFolder, the option 's' indicates to also delete all included files, status=1=succesful deletion, status=0=error 

save('speed.mat', 'speed');

%%% Clear all unneccessary variables (Otherwise code gets slower after each Position)
%clearvars -except Pos
% close all

%end

set(0, 'DefaultFigureVisible', 'on')