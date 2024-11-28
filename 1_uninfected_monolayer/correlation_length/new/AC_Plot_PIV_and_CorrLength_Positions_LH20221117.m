% 10.11.2022 LH
% Adapted from github for infection focus

% readu_correlation_length_WT - Script to read txt files with cellular
% displacement vectors strored, plot them on top of the phase contrast image of cells (or image of their 
% Hoechst-stained nuclei) and calculate the correlation length of movement between neighbouring 
% cells same as in Angenlini et al. PRL 104, 168104 (2010)
%
% 1. Read the phase contrast images and plot on top of them the
% displacement vectors that cells undergo (frames taken every 10 min)
% 2. Construct kymographs of the average with respect to the wound
% displacement vectors, uy as a function of time, t (min) and vertical position, y (μm)  
% 3. Calculation of mean uy along a distance of 100 μm away from the wound 
% and calculated immediately after ablation
% 4. Storage of the results as a .mat file
% Last modified: F. Serrano-Alcalde and E. Bastounis: 2020-06-10

clear all


NumberPos=1;


%for Pos=1:NumberPos %Iterate through all positions
%for Pos = [2,3,5,7,13,18,19,21,24,26] %%% For IF8 analysis
Pos = 1;

clearvars -except Pos; 
close all; clc
addpath 'Functions'
set(0, 'DefaultFigureVisible', 'on') %'on' or 'off' to make figures pop up or not


%% ADAPT FOR EACH DATASET
% load('ExperimentVariables/Experiment_IF1_20220630.mat'); % INFECTION FOCUS 1 30.06.2022
% load('ExperimentVariables/Experiment_IF2_20220725.mat'); % INFECTION FOCUS 2 25.07.2022 
% load('ExperimentVariables/Experiment_IF3_20220811.mat'); % INFECTION FOCUS 3 11.08.2022
% load('ExperimentVariables/Experiment_IF4_20220826.mat'); % INFECTION FOCUS 4 26.08.2022
% load('ExperimentVariables/Experiment_IF5_20220922.mat'); % INFECTION FOCUS 5 22.09.2022
% load('ExperimentVariables/Experiment_Bleach1_20221017.mat'); % Bleach 1      17.10.2022
% load('ExperimentVariables/Experiment_Bleach2_20221018.mat'); % Bleach 2      18.10.2022
% load('ExperimentVariables/Experiment_Bleach3_20221019.mat'); % Bleach 3      19.10.2022
% load('ExperimentVariables/Experiment_IF6_20221125.mat'); % INFECTION FOCUS 6 25.11.2022
% load('ExperimentVariables/Experiment_TFM1_20221127.mat'); % TFM 1             27.11.2022
% load('ExperimentVariables/Experiment_TFM2_20221130.mat'); % TFM 2             30.11.2022
% load('ExperimentVariables/Experiment_IF8_20221213.mat');  % INFECTION FOCUS 8 13.12.2022
%load('ExperimentVariables/Experiment_TFM3_20221221.mat'); % TFM 3             21.12.2022
load('ExperimentVariables/Experiment_pos8_TFM6_2023_04_05.mat'); % TFM 6 pos_8 05_04_2023, analysed 06/11/2023



%%% Homeoffice
% load('\\storage.ag-bastounis.imit.uni-tuebingen.de\mn1111/lara_hundsdorfer/Code/Code_Lara/ExperimentVariables/Experiment_Bleach1_20221017.mat'); % Bleach 1      17.10.2022
% Dir='C:\Daten\PhD_2022\Homeoffice_22021108\DataHO/';
% Dir= 'D:\imaging_data_harddrive\20221213_FRET_MDCK_WT_EcatKO_no_and_with_EKAREV-NLS_glass_JAT607_PD_EGFRi_40x_op1_028umpx/';

ScaleArrow = 5;
overlap    = 24;
Fontsize   = 20;
%% Script to plot vectors on nuclei AND calculate correlation length analysis (kymograph)
kin=2; kfin=timesteps;    % First and last frames
speed = zeros(1,timesteps); % Prepare array for speed calculation

%%% TEST
% NumberPos = 2;
%kin = 2; kfin=4;


if ~exist([Dir 'Results_summary/Save_Workspace'],'dir') 
    mkdir([Dir 'Results_summary/Save_Workspace']) 
end



%for Pos=[2]
    tic
    %if ~exist([Dir 'Urapiv/Pos' num2str(Pos) '/Pos' num2str(Pos) '_speed.mat'],'file')
    if ~exist([Dir 'Urapiv_raul/Pos' num2str(Pos) '/Pos' num2str(Pos) '_PIVxxx.tif'],'file')
        disp(Pos)
%        DirPIVresults = [Dir 'Urapiv/Pos' num2str(Pos) '/PIV_Results_temporary']; mkdir(DirPIVresults)
        fil_txt = [Dir 'Urapiv_raul/Pos' num2str(Pos) '/Txt/' folder_nuclei]; %txt files from PIV
        all_velocities=[];
        kymograph=zeros(26,kfin);% 26 because of R_max=25 and R_step=1
    
        % Plotting vectors on nuclei
        
        for k= kin:kfin  
            Fontsize   = 20;
            fprintf('Pos%2.0f, timestep%2.0f \n', Pos, k)
%             fig1=figure; sgtitle(['Pos' num2str(Pos) ' ' num2str(round(((k*delta+Imaging_minutes_post_infection)/60),2)) ' hpi'])  % Title for whole figure
%             fig1.Position(3:4) = [750 400]; % Adapt size of the figure, here 5 columns and 2 rows thus a 5:2 ratio is useful
%             subplot(1,2,1); imagesc(img_NC);colormap gray;axis image;hold on;

            fig0=figure; 
            %img_NC  = double(imread([Dir 'Raw_Tif_separated/Pos' num2str(Pos) '/Pos' num2str(Pos) '_' folder_nuclei '.tif'],k)); % image of nuclei
            img_NC  = double(imread([Dir 'Pos' num2str(Pos) '_' 'TFM6_' folder_nuclei '.tif'],k)); % image of nuclei
            %%% fixed size white canvas
            [imgSize,~] = size(img_NC); % size of image
            canvas = imgSize*fcal;               % size of white area/ background = same size as image in um
            patch([0 0 canvas canvas],[0 canvas canvas 0],[1 1 1],'EdgeColor',[1 1 1])
            axis([0 canvas 0 canvas])
            hold on;
            imagesc([1:length(img_NC)]*fcal,[1:length(img_NC)]*fcal,img_NC);
            colormap gray; axis square;
            
            filename = ([fil_txt sprintf('%3.3d.txt', k ) ]);
            outl = 10;
            umax=100;
            s6nl= 0; %normally 0
            vec = load(filename); 
               n1 = min(find(diff(vec(:,1))<0));
               n2 = length(vec)/n1;
    
               vec=reshape(vec,n1,n2,5);
    
               x    = vec(:,:,1);
               y    = vec(:,:,2);
               u    = vec(:,:,3);
               v    = vec(:,:,4);
               s6n  = vec(:,:,5);
    
               if s6nl>0
    
                  % signal to noise check
    
                  s6n_low = find(s6n<s6nl);
    
                  w = 3; rad = 2;
                  kernel = zeros(2*w+1);
                  for i=1:2*w+1
                      kernel(i,:) =  exp(-((i-w-1)^2+[-w:w].^2)/rad^2);
                  end
                  kernel(w+1,w+1)=0;
                  kernel = kernel/sum(sum(kernel));
    
                  u(s6n_low) = 0;
                  v(s6n_low) = 0;
    
                  tmpv = (conv2(v,kernel,'same'));
                  tmpu = (conv2(u,kernel,'same'));
    
                  % Let's throw the outlayers out:
    
                  u(s6n_low) = tmpu(s6n_low); 
                  v(s6n_low) = tmpv(s6n_low); 
                  u(s6n<s6nl)=NaN;
                  v(s6n<s6nl)=NaN;
    
                 uu=sqrt(u.^2+v.^2);
                 u(uu>umax)=NaN;
                 v(uu>umax)=NaN;
              end 
    
              if outl>0
    
                  % Adaptive Local Median filtering
    
                  w = 2; rad = 1;
                  kernel = zeros(2*w+1);
                  for i=1:2*w+1
                      kernel(i,:) =  exp(-((i-w-1)^2+[-w:w].^2)/rad^2);
                  end
                  kernel(w+1,w+1)=0;
                  kernel = kernel/sum(sum(kernel));
    
                  tmpv = (conv2(v,kernel,'same'));
                  tmpu = (conv2(u,kernel,'same'));
    
                  lmtv_p = mean(mean(tmpv(2:end-1,2:end-1))) + ...
                         outl*std(reshape(tmpv(2:end-1,2:end-1),(n1-2)*(n2-2),1));
                  lmtv_m = mean(mean(tmpv(2:end-1,2:end-1))) - ...
                         outl*std(reshape(tmpv(2:end-1,2:end-1),(n1-2)*(n2-2),1));
                  lmtu_p = mean(mean(u(2:end-1,2:end-1))) + ...
                         outl*std(reshape(u(2:end-1,2:end-1),(n1-2)*(n2-2),1));
                  lmtu_m = mean(mean(u(2:end-1,2:end-1))) - ...
                         outl*std(reshape(u(2:end-1,2:end-1),(n1-2)*(n2-2),1));
    
                  u_out_p = find(u>lmtu_p);
                  u_out_m = find(u<lmtu_m);
                  v_out_p = find(v>lmtv_p);
                  v_out_m = find(v<lmtv_m);
    
                  % Let's throw the outlayers out:
    
                  u(u_out_m) = tmpu(u_out_m); 
                  u(v_out_m) = tmpu(v_out_m); 
                  v(u_out_m) = tmpv(u_out_m); 
                  v(v_out_m) = tmpv(v_out_m); 
    
                  u(u_out_p) = tmpu(u_out_p); 
                  u(v_out_p) = tmpu(v_out_p); 
                  v(u_out_p) = tmpv(u_out_p); 
                  v(v_out_p) = tmpv(v_out_p); 
    
              end  
    
            qq=1;
%             hold on
%             quiver(x(1:qq:end,1:qq:end),y(1:qq:end,1:qq:end),u(1:qq:end,1:qq:end)*10*fcal,v(1:qq:end,1:qq:end)*10*fcal,'AutoScale','off','color',[0 1 0]);%htt=text(110,860,'20 \mum','fontsize',24); 
            hold on
            quiver(x(1:qq:end,1:qq:end)*fcal,y(1:qq:end,1:qq:end)*fcal,u(1:qq:end,1:qq:end)*ScaleArrow*fcal,v(1:qq:end,1:qq:end)*ScaleArrow*fcal,'AutoScale','off','color','g');%htt=text(110,860,'20 \mum','fontsize',24); 
            axis off

            u=u';v=v';
            %MeanDisplacement = mean(mean(sqrt((u*fcal).^2+(v*fcal).^2)))
            speed(1,k)= mean(mean(sqrt((u*fcal).^2+(v*fcal).^2)))/delta*60; %um/h 
            
            u(isnan(u))=0;
            v(isnan(v))=0;
            %%% Save results of PIV vectors, plot them
% % %             print(fig0,'-dtiff',[DirPIVresults '/PIV_result_' int2str(k)]);  % Try getframe instead%%%%%%%% ???
% % %             img_PIV = imread([DirPIVresults '/PIV_result_' int2str(k) '.tif']); %%% Open the .tif per timestep and append one multi tif
% % %             imwrite(img_PIV,[Dir 'Urapiv/Pos' num2str(Pos) '/Pos' num2str(Pos) '_PIV.tif'], 'WriteMode', 'append',  'Compression','none') % saves image stack e.g. Pos1_PIV.tif            
% % %              
            Save_figure = getframe(fig0); %or Save_figure = getframe(gcf); 
            imwrite(Save_figure.cdata, [Dir 'Urapiv_raul/Pos' num2str(Pos) '/Pos' num2str(Pos) '_PIV.tif'], 'WriteMode', 'append',  'Compression','none');% saves image stack e.g. Pos1_PIV.tif  

    
            %% Correlacion Angelini et al. (v1)
            %Subtract of mean
            mean_u=mean(mean(u));
            mean_v=mean(mean(v));
            u=u-mean_u;
            v=v-mean_v;
    
            %Generate matrix with distances
            Generar_distancias_i=1:1:length(u);
            Generar_distancias_i=repmat(Generar_distancias_i,length(Generar_distancias_i),1);
            Generar_distancias_j=(1:1:length(u))';
            Generar_distancias_j=repmat(Generar_distancias_j,1,length(Generar_distancias_j));
    
            %Initilization of matrices
            P_escalar_superior=zeros(length(u));
            P_escalar_inferior=zeros(length(u));
            P_escalar=[];
            numero_R=1;
            R_ant=0;
            bordes=0; % Changing this value we can remove the correlation of the edges. (pixels)
            R_step=1; % Step of the radius discretization (pixels)
            R_max=25; % Maximum radius (pixels). This corresponds to distance R_max*overal*fcal
            for R=0:R_step:R_max
                for i=1+bordes:length(u)-bordes
                    for j=1+bordes:length(u)-bordes
                        Distancias=sqrt((i-Generar_distancias_i).^2+(j-Generar_distancias_j).^2);
                        [Posiciones_y,Posiciones_x]=find(Distancias<=R & Distancias>=R_ant);
                        vector_0=[u(i,j) v(i,j)];
                        for escalar=1:length(Posiciones_x)
                            vector_Posicion=[u(Posiciones_x(escalar),Posiciones_y(escalar)) v(Posiciones_x(escalar),Posiciones_y(escalar))]';
                            P_escalar(i,j,escalar)=vector_0*vector_Posicion;
                        end
                        P_escalar_superior(i,j)=mean(P_escalar(i,j,find(P_escalar(i,j,:)~=0 & ~isnan(P_escalar(i,j,:)))));
                        P_escalar_inferior(i,j)=vector_0*vector_0';
                    end
                end
                Coef(numero_R, k)=sum(sum(P_escalar_superior))/sum(sum(P_escalar_inferior)); %Correlation for each R iteration for each image(k)
                valor_Radio(numero_R,k)=R; % Radius for each iteration for each image(k)(the same for all images)
                numero_R=numero_R+1;
    
                R_ant=R;
    
                %%%%
                %Image of the correlation for a given radius (R)
                %Be careful (and comment) if R_step and/or k (number of images) is too large!!!
                Coef_2D=P_escalar_superior./P_escalar_inferior; %Matrix of the individual correlation of each pixel for the radius R
                % Plot the 2D map of the correlation coefficinet for each R considered
            %     figure
            %     imagesc([1:length(Coef_2D)]*24*fcal,[1:length(Coef_2D)]*24*fcal, Coef_2D); caxis([-1 1]);colorbar;axis image
            %     set(gca,'FontSize',18,'DefaultAxesFontName', 'Arial');title(R);xlabel('um');ylabel('um');
                %%%%
            end
    
            valor_Radio(:,k)=valor_Radio(:,k)*overlap*fcal; %Change of units from pixels to um
            time(k)=k*delta; %time in 'delta' units
    
            %% Length from equation
            %Fit equation
            eqexpon=fittype('exp(-R/R0)', 'independent',{'R'}, 'coefficients',{'R0'});
            myfit=fit(valor_Radio(:,k),Coef(:,k),eqexpon,'startpoint', [50], 'maxiter',[20000]);
            R0_equation(k)=myfit.R0*2; %Value of R0 from equation
            
            % Plot Corr Coefficient over radius/distance per timestep
            fig0=figure;
            plot(myfit,valor_Radio(:,k),Coef(:,k))
            xlabel('Radius (um)','FontSize',Fontsize)
            ylabel('C(Radius)','FontSize',Fontsize) % Correlation coefficient depending on radius
            ylim([-0.2 1])
            title(['Pos' num2str(Pos) ' ' num2str(round(((k*delta+Imaging_minutes_post_infection)/60),2)) ' hpi'], 'FontWeight','normal', 'FontSize',Fontsize)
            %print(fig,'-dtiff',name);
            Save_figure = getframe(fig0); %or Save_figure = getframe(gcf); 
            imwrite(Save_figure.cdata, [Dir 'Urapiv_raul/Pos' num2str(Pos) '/Pos' num2str(Pos) '_CorrCoefTimesteps.tif'], 'WriteMode', 'append', 'Compression','none');
            
            

            kymograph(:,k)=Coef(:,k);
    
            %% Length from minimum value
            [min_cor, position]=min(Coef(:,k));
            R0_minimum(k)=valor_Radio(position,k);
            %close all
            save([Dir 'Urapiv_raul/Pos' num2str(Pos) '/Pos' num2str(Pos) '_CorrCoefKymograph.mat'], 'kymograph')
            Fig2_Correlation=[time; R0_equation; R0_minimum];
            save([Dir 'Urapiv_raul/Pos' num2str(Pos) '/Pos' num2str(Pos) '_Correlation.mat'],'Fig2_Correlation' )
            
            % Save whole workspace
            %save([Dir 'Results_summary/Save_Workspace/CorrLength_Pos' num2str(Pos) '.mat']) 
    
         end 
    
        %clear -Dir % Clear everything except Dir
    
        Fig2Data = load([Dir 'Urapiv_raul/Pos' num2str(Pos) '/Pos' num2str(Pos) '_Correlation.mat']);
        Fig3Data = load([Dir 'Urapiv_raul/Pos' num2str(Pos) '/Pos' num2str(Pos) '_CorrCoefKymograph.mat']);
        time        = Fig2Data.Fig2_Correlation(1,:);
        R0_equation = Fig2Data.Fig2_Correlation(2,:);
        %R0_minimum  = Fig2Data.Fig2_Correlation(3,:);
        kymograph   = Fig3Data.kymograph;
    
        % Plot of evolution of mean correlation length in each frame
        fig2=figure;
        hold on
        plot((time(2:end)+Imaging_minutes_post_infection)/60,R0_equation(2:end))
        %plot(time,R0_minimum)
        xlabel('Time (hpi)','FontSize',Fontsize)
        ylabel('Correlation length (um)','FontSize',Fontsize)
        title(['Pos' num2str(Pos)],'FontSize',Fontsize, 'FontWeight','normal')
        ylim([0 100])
        %legend('equation','minimum value')
        hold off
        Save_figure = getframe(fig2); imwrite(Save_figure.cdata, [Dir 'Urapiv_raul/Pos' num2str(Pos) '/Pos' num2str(Pos) '_Correlation.tif']); 
        %print(fig2,'-dtiff',[Dir 'Urapiv/Pos' num2str(Pos) '/Pos' num2str(Pos) '_Correlation']);
    
        % Kymograph where you can see how the correlation coefficient changes as a
        % function of time and space
        fig3=figure;
        imagesc((time+Imaging_minutes_post_infection)/60,[1:size(kymograph,1)]*overlap*fcal, kymograph); colormap redblue;Col=colorbar; Col.Label.String='Correlation coefficient';
        caxis([-1 1])
        xlabel('Time post-infection (h)','FontSize',Fontsize)
        ylabel('Radius (um)','FontSize',Fontsize)
        title(['Pos' num2str(Pos)], 'FontWeight','normal')
        set(gca,'FontSize',Fontsize);
        Save_figure = getframe(fig3); imwrite(Save_figure.cdata, [Dir 'Results_summary/Pos' num2str(Pos) '_CorrelationCoef.tif']); 
        %print(fig3,'-dtiff',[Dir 'Results_summary/Pos' num2str(Pos) '_CorrelationCoef']);
    
% %         %%% Sort and delete files
% %         DeleteFolder = DirPIVresults; 
% %         status = rmdir(DeleteFolder, 's'); % remove the directory DeleteFolder, the option 's' indicates to also delete all included files, status=1=succesful deletion, status=0=error 
    else
        fprintf('For Pos%2.0f file Pos_PIV.tif already exists \n', Pos)
    end
    save([Dir 'Urapiv_raul/Pos' num2str(Pos) '/Pos' num2str(Pos) '_speed.mat'], 'speed');
    
    %%% Clear all unneccessary variables (Otherwise code gets slower after each Position)
    clearvars -except Pos
    % close all

    %end
    toc

set(0, 'DefaultFigureVisible', 'on')