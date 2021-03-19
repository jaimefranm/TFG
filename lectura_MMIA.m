%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %                   UNIVERSITAT DE VALÃˆNCIA                            %
% %              Image Processing Laboratory (IPL)                        %
% %               ASDC (Asim Science Data Center)                         %
% %                                                                       %
% %                Javier Navarro-GonzÃ¡lez               (nagonja@uv.es) %
% %                                                                       %
% %                      Dec.      2018                                   %
% %                                                                       %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Script is an example for reading MMIA level0 cdf New Version Files
% (After 2018DOY280)
% This scrip is preconfigured with the trigger #13 UV or #132 Bergen 
% TGF List DATE&TIME 2018 DOY 283 (l0 time 13:01:33.085950
%
% The Software makes an quick overview of the MMIA file
%
% 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % GO INTO SCRIPT AND SET: %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path to the files
% dirandfile='/home/navarro/TGF_20_NOV/mmiatriggeredobservationtm_observation_time_2018-10-10_13-01-32_observation_id_55591.cdf';%283
% 
% Path to the cdf_Matlab_patch
%  addpath '/matlab_cdf363_patch-64'
% 

% MODIFICACIÓN JESÚS LÓPEZ 

% CONCATENA TODOS LOS FICHEROS DE UN DIRECTORIO AL IGUAL QUE TODOS LOS
% FRAMES DE LOS MISMOS FICHEROS



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% COMENTARIO 2020-10-16 
% ESTA VERSION EMPLEA EL "t_corrected_l1" y completa el tiempo según el
% sample rate de 100 Khz. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Reset Values:
clearvars
clear all
close all

% SET cdf MATLAB patch
% if isunix
%     addpath '../matlab_cdf363_patch-64'
% else
     addpath '/Users/jaimemorandominguez/Desktop/TFG/Lectura_MMIA/cdf_patch/matlab_cdf370_patch'
% end

% Directorio datosa de MMIA a unir
% str='E:\Doctorado\Academico\ASIM\MMIA\ASIM_mmia_UV_overview_2_IAA_v0\MMIA_DATA_new\191030\06_05_46\';
str='/Users/jaimemorandominguez/Desktop/TFG/Lectura_MMIA/190722/09_03_17';

% ---------->LMA=load('E:\Doctorado\Academico\LMA\LMA viewer program\saved\LMA_191030_06_05_46.txt'); 
% f_LMA=find(LMA(:,3)==0);
% LMA(f_LMA,:)=[]; 
% % GLM info
% ---------->GLM=load('E:\UPC\ASIM\MMIA_Positive_CGs\Events\201111\GLM\output\GLM.txt'); 

% % LIS INFO -  Opción 3 programa Icar.
% ---------->load('E:\UPC\ASIM\MMIA_Positive_CGs\Events\201111\LIS\ISS_LIS_SC_V1.0_20201111_NQC_21914.mat');
% load('F:\Doctorado\Academico\LMA\LMA viewer program\LMA_LIS\LIS_data\output\ISS_LIS_SC_V1.0_20200919_NQC_21105.mat');


% %%%%%%%%%%%%%%%
% % LINET
% ---------->filename = 'E:\UPC\ASIM\MMIA_Positive_CGs\Events\201111\Linet\201111_00_43_17_Linet.txt';


% EN CASO QUE LOS DATOS PROCESADOS DESDE PROGRAMA DE ÍCAR NO ESTÉN ORDENADOS EN TIEMPO 
[~,idx] = sort(LIS_all(:,1));
sor_LIS= LIS_all(idx,:);
clear LIS_all
LIS_all=sor_LIS;


t_ajuste=0.15; % Tiempo de búsqueda = tiempo MMIA +/- t_ajuste

d_lat_lon=10/111.111; % Delta de búsqueda alrededor del rayos del LMA. 


% lis_sort para flashes
% Flash 1-> lis_sort=0
% Flash 2-> lis_sort=0
% Flash 3-> lis_sort=0
% Flash 4-> lis_sort=1

lis_sort=0; % Fuente de LIS para ajuste. 0=Valor de mayor radiance, =1, siguiente valor....
glm_sort=0;  % Fuente de GLM para ajuste. 0=Valor de mayor energía, 1 siguiente valor mayor....
mmia_sort=0; 

tresh_frame=100; % UMBRAL DE TRESHOLD PARA CADA FOTÓMETRO

folder_name=(str);
files=dir(fullfile(folder_name,'*.cdf'));
zz=1;
n_files=length(files); 

% Predefino las metrices de los fotones. 
PHOT1Data_all_tmp=[];
PHOT2Data_all_tmp=[];
PHOT3Data_all_tmp=[];
t_vector_tmp=[];

while zz<=n_files
    dirandfile=[str,files(zz).name];
    only_chu=0;     
    INFO=spdfcdfinfo(dirandfile);
    L=spdfcdfinfo(dirandfile);
    DATA=spdfcdfread(dirandfile);
    
    [a,numvars]=size(DATA);
    
    for nv=1:numvars
        eval([L.Variables{nv} '= double(DATA{nv});' ])
    end
    
    
    
    % Datos de la ISS
    fprintf('ISS data     lat: %s, lon: %s \n',num2str(latitude(1)),num2str(longitude(1)))
    fprintf('ISS attitude yaw: %s, pitch: %s, roll: %s, Meridian_angle: %s\n ',...
        num2str(iss_yaw(1)),num2str(iss_pitch(1)),num2str(iss_roll(1)),num2str(iss_angle_from_meridian(1)))
    
    
    %Frame en el que hay un trigger con MXGS
    frame=1;%find(MXGSTrigger);
    
    % Data time range:
    Init_T_File_=datestr(raw_datetime(1),'mmmm dd, yyyy HH:MM:SS.FFF');
    End_T_File_=str2double(datestr(raw_datetime(1),'SS.FFF'))+length(raw_datetime)*8320*10/1e6;
    fprintf('Data available from: %s till %f\n',Init_T_File_,End_T_File_)

    %Extraemos las curvas de los fotometros para 1 frame
    
    format long
    % CORRECCIÓN DEL TIEMPO frame_time_phot según el que se espcifica en
    % corrected_datetime_level1.
    t_corrected_l1=(str2num(datestr(corrected_datetime_level1(1,1),'SS.FFF'))-str2num(datestr(frame_time_phot(1,1),'SS.FFF'))) +  str2num(datestr(frame_time_phot(1,:),'SS.FFF'));
    
    % Difenrencia de tiempo entre "frames"
    t_btw_corrected_frame=t_corrected_l1(2:end,1) - t_corrected_l1(1:end-1,1);
    if (t_btw_corrected_frame(:)~=0) %
        
        PHOT1Data_all = reshape(PHOT1_photon_flux',[],1);
        PHOT1Data_all=PHOT1Data_all';
        PHOT2Data_all = reshape(PHOT2_photon_flux',[],1);
        PHOT2Data_all=PHOT2Data_all';
        PHOT3Data_all = reshape(PHOT3_photon_flux',[],1);
        PHOT3Data_all=PHOT3Data_all';
        
        tresh_phot1=find(PHOT1Data_all>=tresh_frame);
        tresh_phot2=find(PHOT2Data_all>=tresh_frame);
        tresh_phot3=find(PHOT3Data_all>=tresh_frame);
        
        frames_tresh=find(tresh_phot1(1,2:end)-tresh_phot1(1,1:end-1)>2);

        PHOT1Data_all(tresh_phot3)=[];
        PHOT2Data_all(tresh_phot3)=[];
        PHOT3Data_all(tresh_phot3)=[];
        
        
        PHOT1Data_all_tmp=[PHOT1Data_all_tmp,PHOT1Data_all];
        PHOT2Data_all_tmp=[PHOT2Data_all_tmp,PHOT2Data_all];
        PHOT3Data_all_tmp=[PHOT3Data_all_tmp,PHOT3Data_all];
        h_asim_corr=str2num(datestr(corrected_datetime_level1(1),'HH'))*3600;
        m_asim_corr=str2num(datestr(corrected_datetime_level1(1),'MM'))*60;
        s_asim_corr=str2num(datestr(corrected_datetime_level1(1),'SS.FFF'));
        t_ini_asim_corr=h_asim_corr+m_asim_corr+s_asim_corr;
        sample_r=1e-5;
        
        trigger_length=length(PHOT1Data_all(:)); %%DEBE SER ESTE FOTÓMETRO YA QUE SE ELIMINA LOS FAKES FRAMES
        t_asim_corr(1:trigger_length,1)=0;
        t_asim_corr(1)=t_ini_asim_corr;
        
        for i=2:trigger_length
            t_asim_corr(i,1)= t_asim_corr(i-1,1)+sample_r;
        end
        
        t_vectorL1=[t_vector_tmp;t_asim_corr];
        t_vector_tmp=t_vectorL1;

        % GRÁFICO DE LOS CHU
        CHU1Data=CHU1_photon_flux;
        CHU2Data=CHU2_photon_flux;
        CHU1_pixel_longitude;
        CHU1_pixel_latitude;
        CHU2_pixel_latitude;
        CHU2_pixel_longitude;
      end % Fin del ciclo con solo CHU!    
        for frame=1:length(CHU1Data_exists)
            
            if CHU1Data_exists(frame)==1 && CHU2Data_exists(frame)==1
                indiim=frame;
                
                %Reconstruct Images cut:
                dim_row=chu_maximum_row(indiim)-chu_minimum_row(indiim)+1;
                dim_column=chu_maximum_column(indiim)-chu_minimum_column(indiim)+1;
                       % GRÁFICO DE LOS CHU
                CHU1Data=CHU1_photon_flux;
                CHU2Data=CHU2_photon_flux;
                CHU1_pixel_longitude;
                CHU1_pixel_latitude;
                CHU2_pixel_latitude;
                CHU2_pixel_longitude;
                A_CHU1=reshape(CHU1Data(indiim,:),dim_column,dim_row);
                clc=reshape(CHU2Data(indiim,:),dim_column,dim_row);
                lat_chu1=reshape(CHU1_pixel_latitude(indiim,:),dim_column,dim_row);
                lon_chu1=reshape(CHU1_pixel_longitude(indiim,:),dim_column,dim_row);
                
                lat_chu2=reshape(CHU2_pixel_latitude(indiim,:),dim_column,dim_row);
                lon_chu2=reshape(CHU2_pixel_longitude(indiim,:),dim_column,dim_row);
                
                %         figure(44)
                %         contourf(lon_chu1,lat_chu1,A_CHU1,50,'LineStyle','none')
                %         title 'CHU 1 (337.0/ 5nm)'
                %         daspect([1 1 1])
                
               
                %         h = pcolor(A_CHU2);
                %         set(h, 'LineStyle','none');
                %         title 'CHU 2 (777.4)'
                %         caxis([300 2000])
                %         axis off
                
                flis_frame=find(LIS_all(:,1)>=t_asim_corr(1) & LIS_all(:,1)<=t_asim_corr(end));
                fglm_frame=find(GLM(:,1)>=t_asim_corr(1) & GLM(:,1)<=t_asim_corr(end));
                figure(46)
                axis off
                h = pcolor(lon_chu1 ,lat_chu1 ,A_CHU1);
                set(h, 'LineStyle','none');
                title 'PHOT 1 (337.0)'
                hold on
                plot(LIS_all(flis_frame,4),LIS_all(flis_frame,3),'or')
                plot(GLM(fglm_frame,3),GLM(fglm_frame,2),'+k')
                
                if exist('A_CHU2')
                    figure(45)
                    h = pcolor(lon_chu2 ,lat_chu2 ,A_CHU2);
                    hold on
                    plot(LIS_all(flis_frame,4),LIS_all(flis_frame,3),'or')
                    plot(GLM(fglm_frame,3),GLM(fglm_frame,2),'+k')
                    set(h, 'LineStyle','none');
                    title 'CHU 2 (777.4)'
                end

                

                imagenFOV1=zeros(1026,1056);
                imagenFOV2=zeros(1026,1056);
            end
        end
                
                
                if (t_btw_corrected_frame(:)~=0)

                    
                end
                


        zz=zz+1;  %% Finaliza concatenación de los fiecheros MMIA
    

end



MMIA_all(:,1)=t_vectorL1;
MMIA_all(:,2)=PHOT1Data_all_tmp';
MMIA_all(:,3)=PHOT2Data_all_tmp';
MMIA_all(:,4)=PHOT3Data_all_tmp';




% Caso Colombia LMA
% minlat=min(LMA(:,4))-d_lat_lon;
% maxlat=max(LMA(:,4))+d_lat_lon;
% minlon=min(LMA(:,5))-d_lat_lon;
% maxlon=max(LMA(:,5))+d_lat_lon; 


% minlat=latitude(1)-3;
% maxlat=latitude(1)+3;
% minlon=longitude(1)-3;
% maxlon=longitude(1)+3; 

minlat=7.58-3.5;
maxlat=7.58+3.5;
minlon=-72.6-3.5;
maxlon=-72.6+3.5; 

% % 
% 
% minlat=5;
% maxlat=9;
% minlon=-75.3;
% maxlon=-72.3; 
% 
% 
% 
%%

% % %NOTA: DESACTIVAR SI SE DESEA PLOTEAR TODOS LOS DATOS DE LIS 

ind_tmmia1=find( (LIS_all(:,1) <= (MMIA_all(1,1) - t_ajuste)) );
ind_tmmia2=find( (LIS_all(:,1) >= (MMIA_all(end,1) + t_ajuste)) );
ind_tmmia=[ind_tmmia1;ind_tmmia2];
LIS_original=LIS_all; 
LIS_all(ind_tmmia,:)=[]; % Elimino los datos de LIS fuera del tiempo de MMIA; 
% 




% %%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ENCUENTRA AJUSTE POR CADA FRAME
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [t_max_frame,LIS_tmp2]=t_frameMAX(t_matrix_L1, PHOT1_photon_flux, PHOT2_photon_flux, PHOT3_photon_flux,LIS_all,t_ajuste); 
% 
% LIS_all=LIS_tmp2; % LIS CON TIEMPO CORREGIDOS POR CADA !DT!

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%
% % GLM
% %%%%%%%%%%%%%%%%%
% 
% 
ind_tglm1=find( (GLM(:,1) <= (MMIA_all(1,1) - t_ajuste)) );
ind_tglm2=find( (GLM(:,1) >= (MMIA_all(end,1) + t_ajuste)) );
% ind_tglm=find( GLM(:,1)<= (MMIA_all(1,1)-t_ajuste) | GLM(:,1)>= (MMIA_all(end,1)+t_ajuste) ); 
ind_tglm=[ind_tglm1;ind_tglm2];
GLM(ind_tglm,:)=[]; % Elimino los datos de LIS fuera del tiempo de MMIA; 


% % Filtro para los datos dentro del área de LMA/MMIA 
ind_inLIS=find( LIS_all(:,3) <= maxlat & LIS_all(:,3) >= minlat & LIS_all(:,4) <= maxlon & LIS_all(:,4) >= minlon);
ind_inGLM=find( GLM(:,5) <= maxlat & GLM(:,5) >= minlat & GLM(:,6) <= maxlon & GLM(:,6) >= minlon);

% GLM Y LIS PARA EL TIMEPO DE MMIA+/-dt_ajust Y LA REGIÓN DE LMA + 10 km AL
% REDEDOR
GLM_in=GLM(ind_inGLM,:);
LIS_in=LIS_all(ind_inLIS,:);


% t_ajuste=0; 
t1=MMIA_all(1,1)-t_ajuste;
t_end=MMIA_all(end,1)+t_ajuste; 

% Habilitar si se quiere integrar o suavizar MMIA. 
[int_MMIA]=integral_MMIA2(MMIA_all);
figure
plot(MMIA_all(:,1),MMIA_all(:,4),'b')
hold on
plot(int_MMIA(:,1),(int_MMIA(:,3)),'r')



figure
plot(MMIA_all(:,1),MMIA_all(:,2),'b')
hold on
plot(int_MMIA(:,1),(int_MMIA(:,2)),'r')


% Integral de LIS cada 2 ms
[int_LIS,int_GLM]=integral_GLM_LIS2(GLM_in,LIS_in,t1,t_end);

% %%%%%%%%%%% % % Max MMIA 337 ((MMIA_all(:,2)) o 777 ((MMIA_all(:,4))
% % MAX MMIA 777
[vmax_mmia posmax_mmia]=max(MMIA_all(:,4));
t_mmia_max=MMIA_all(posmax_mmia,1);


% [~,idx_mmia] = sort(int_MMIA(:,2)); % Ordena según max. radiance
% t_mmia_max=int_MMIA(idx_mmia(end-mmia_sort),1);

% [~,idx_mmia] = sort(MMIA_all(:,4)); % Ordena según max. radiance
% t_mmia_max=MMIA_all(idx_mmia(end-mmia_sort),1);


clear posmax_glm
% [~,idx] = sort(int_LIS(:,2)); % Ordena según integral max. radiance
% t_lis_max=int_LIS(idx(end-lis_sort),1);
% dt_lis=t_lis_max-t_mmia_max; % Corrección de tiempo para MMIA según LIS

%Nota: Para el flash 4; lis_sort=1;

lis_sort=0; 
[~,idx] = sort(LIS_in(:,5)); % Ordena según max. radiance
t_lis_max=LIS_in(idx(end-lis_sort),1);
dt_lis=t_lis_max-t_mmia_max; % Corrección de tiempo para MMIA según LIS

% Max Energy GLM
% [~, posmax_glm]=sort(int_GLM(:,2));
% t_glm_max=int_GLM(posmax_glm(end-glm_sort),1);
% dt_glm=t_glm_max-t_mmia_max; % Corrección de tiempo para MMIA según LIS

glm_sort=0;
[~, posmax_glm]=sort(GLM_in(:,7));
t_glm_max=GLM_in(posmax_glm(end-glm_sort),1);
dt_glm=t_glm_max-t_mmia_max; % Corrección de tiempo para MMIA según LIS

dt_glm_lis=dt_glm-dt_lis;

linet = import_NewLinet2(filename); 



ind_tlinet1=find( (linet(:,11) <= (t1)) );
ind_tlinet2=find( (linet(:,11) >= (t_end)) );
ind_tlinet=[ind_tlinet1;ind_tlinet2];
linet(ind_tlinet,:)=[]; % Elimino los datos de LIS fuera del tiempo de MMIA; 

ind_inLINET=find( linet(:,5) <= maxlat & linet(:,5) >= minlat & linet(:,6) <= maxlon & linet(:,6) >= minlon);


max_phot_mmia337=max(MMIA_all(:,2)); % Max MMIA 
max_phot_mmia180=max(MMIA_all(:,3)); % Max MMIA 
max_phot_mmia777=max(MMIA_all(:,4)); % Max MMIA 
max_rad_LIS=max(LIS_in,5); % Max radiance LIS
max_ener_GLM=max(GLM_in,7); % Max radiance LIS
max_h_LMA=max(LMA(:,6)); % Max H LMA
max_pow_LMA=max(LMA(:,7)); % Max Pow LMA


close all
n_bins=60;
% n_bins=20;
% n_bins=30;%length(LMA(:,1));
% dt_mmia_linet=t_mmia_max-linet(2,11);
% dt_dt=0;
dt_dt=dt_glm;
%Flash 3 y 4. Promedio de las detecciones. 
% dt_dt=(dt_glm+dt_lis)/2;
%% AUSTES ADICIONAL DE TIEMPO EMPLEANDO DATOS DE LINET PARA CADA FLASH
% Flash 1
% t_ajuste_cg=0;
% dt_dt=dt_dt-t_ajuste_cg; 
%%%%%%%%%%
% Flash 2 
% t_ajuste_cg=0;
% dt_dt=dt_dt+t_ajuste_cg;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
% FLash 3
% t_ajuste_cg=0.005806825396576; % Reajuste empleando primer fuente LMA 
% dt_dt=dt_dt-t_ajuste_cg;
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%
% Flash4.
% t_ajuste_cg=0;% Reajuste primera fuenta LMA. 
% dt_dt=dt_dt-t_ajuste_cg;
% % 0.000538461699761;
% % 0.001461538300646;  
% % 0.001259615302843
%%%%%%%%%%%%%%%%%%%%%
%%

dt_plot=(dt_dt+10^-2);


%%
%%
% Potencia
figure
height_sub=0.2;
hAxis(1)=subplot(4,1,1); 
% p1=plot(LMA(:,1)-dt_dt,(LMA(:,6))/1000,'.');
% p1=bar(LMA(:,1),(LMA(:,7)),n_bins,'FaceColor',[0.255 0.255 0.204],'BarWidth',4); %,'BarWidth',1
p1=plot(MMIA_all(:,1)+dt_dt,MMIA_all(:,3),'k');
hold on
if ~isempty(ind_inLINET)
    dt_lma_linet=LMA(1,1)- linet(ind_inLINET(1),11);
    for li=1:length(linet(ind_inLINET))
        if linet(ind_inLINET(li),8)>1 % Rayo IC
        p2=plot(linet(ind_inLINET(li),11),0,'.k','MarkerSize',14);
        else
                p3=plot(linet(ind_inLINET(li),11),0,'xk','MarkerSize',12,'LineWidth',1);
        end
    end
end

% Plot CG GLD360 Flash 1
% p_ext=plot( (6*3600+5*60+11.994904999998113),0,'^k','MarkerSize',8);

% Plot CG GLD360 flash2
% p_ext=plot((6*3600+5*60+28.030931528),0,'^k','MarkerSize',8);


% % Plot CG GLD360 Flash4
% p_ext=plot((6*3600+5*60+47.006312914),0,'^k','MarkerSize',8);


% ylabel('H [m]') %'Power [dBW]'
ylabel('Power [dBW]') %
set(gca, 'YGrid', 'off', 'XGrid', 'on')
xlim([MMIA_all(1,1) MMIA_all(end,1)])
pos1 = get( hAxis(1), 'Position' );
pos1(4)=height_sub; % Increase height.
set( hAxis(1), 'Position', pos1) ;
set(gca,'xticklabel',[])
%Flash 1:
% legend([p1 p2 p3 p_ext],'LMA','IC(LINET)','CG(LINET)','CG(GLD360)')

% Flash 2:
% legend([p1 p2 p3 p_ext],'LMA','IC(LINET)','CG(LINET)','CG(GLD360)')

%Flash 3:
% legend(p1,'LMA')

%Flash 4:
% legend([p1 p3 p_ext],'LMA','IC(LINET)','CG(GLD360)')


hAxis(2)=subplot(4,1,2); 
plot(MMIA_all(:,1)+dt_dt,MMIA_all(:,2),'b');
hold on
if ~isempty(linet)
    for li=1:length(linet(ind_inLINET))
        if linet(ind_inLINET(li),8)>1 % Rayo IC
            plot(linet(ind_inLINET(li),11),0,'.k','MarkerSize',14);
        else
              plot(linet(ind_inLINET(li),11),0,'xk','MarkerSize',12,'LineWidth',1);
        end
    end
end

% Plot CG GLD360 Flash 1
% p_ext=plot( (6*3600+5*60+11.994904999998113),0,'^k','MarkerSize',8);

% Plot CG GLD360 flash2
% p_ext=plot((6*3600+5*60+28.030931528),0,'^k','MarkerSize',8);


% Plot CG GLD360 Flash4
% p_ext=plot((6*3600+5*60+47.006312914),0,'^k','MarkerSize',8);

%title 'PHOT 1 (337.0/ 5nm)'
ylabel('Energy [\muWm^-^2]') 
set(gca, 'YGrid', 'off', 'XGrid', 'on')
xlim([MMIA_all(1,1) MMIA_all(end,1)])
pos2=get(hAxis(2), 'Position' );
pos2(4)=height_sub;
set(hAxis(2), 'Position', pos2) ;
set(gca,'xticklabel',[])
legend('337\pm4 nm')

hAxis(3)=subplot(4,1,3); 
plot(MMIA_all(:,1)+dt_dt,MMIA_all(:,4),'r');
hold on
if ~isempty(linet)
    for li=1:length(linet(ind_inLINET))
        if linet(ind_inLINET(li),8)>1 % Rayo IC
            plot(linet(ind_inLINET(li),11),0,'.k','MarkerSize',14);
        else
            plot(linet(ind_inLINET(li),11),0,'xk','MarkerSize',12,'LineWidth',1);
        end
    end
end

% Plot CG GLD360 Flash 1
% p_ext=plot( (6*3600+5*60+11.994904999998113),0,'^k','MarkerSize',8);

% Plot CG GLD360 flash2
% p_ext=plot((6*3600+5*60+28.030931528),0,'^k','MarkerSize',8);

% Plot CG GLD360 Flash4
% p_ext=plot((6*3600+5*60+47.006312914),0,'^k','MarkerSize',8);

%title 'PHOT 1 (337.0/ 5nm)'
ylabel('Energy [\muWm^-^2]') 
set(gca, 'YGrid', 'off', 'XGrid', 'on')
xlim([MMIA_all(1,1) MMIA_all(end,1)])
pos3=get( hAxis(3), 'Position' );
pos3(4)=height_sub;
set( hAxis(3), 'Position', pos3) ;
set(gca,'xticklabel',[])
legend('777.4\pm5 nm')


hAxis(4)=subplot(4,1,4); 
yyaxis left
plot(int_LIS(:,1)+dt_lis,(int_LIS(:,2)))
ylabel('LIS [\muJ m^-^2 sr^-^1 \mum^-^1]')
yyaxis right
plot(int_GLM(:,1),(int_GLM(:,2))*10^15)
ylabel('GLM [pJ]')

hold on
if ~isempty(linet)
    for li=1:length(linet(ind_inLINET))
        if linet(ind_inLINET(li),8)>1 % Rayo IC
            plot(linet(ind_inLINET(li),11),0,'.k','MarkerSize',14);
        else
            plot(linet(ind_inLINET(li),11),0,'xk','MarkerSize',12,'LineWidth',1);
        end
    end
end
set(gca, 'YGrid', 'off', 'XGrid', 'on')
xlim([MMIA_all(1,1) MMIA_all(end,1)])
pos4=get( hAxis(4), 'Position' );
pos4(4)=height_sub;
set(hAxis(4), 'Position',pos4) ;
xlabel('Time [s]')
linkaxes([hAxis(1),hAxis(2),hAxis(3),hAxis(4)],'x');


% Plot CG GLD360 Flash 1
% p_ext=plot( (6*3600+5*60+11.994904999998113),0,'^k','MarkerSize',8);
% Plot CG GLD360 flash2
% p_ext=plot((6*3600+5*60+28.030931528),0,'^k','MarkerSize',8);


% Plot CG GLD360 Flash4
% p_ext=plot((6*3600+5*60+47.006312914),0,'^k','MarkerSize',8);




%%

% Código para grabar imágenes en PDF segun tamaño original
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(gcf,'06_05_11','-dpdf','-r0')

%%
% lma_sta=load('F:\Doctorado\Academico\LMA\LMA viewer program\map\stations_barranca_v2.txt');
% lma_buffer=importdata('F:\Doctorado\Academico\LMA\LMA viewer program\map\balma_buffer.txt',','); 
% figure
% plot(LMA(:,5),LMA(:,4),'.k','LineWidth',2,'MarkerSize',7)
% hold on
% plot(lma_sta(:,2),lma_sta(:,1),'^k')
% plot(GLM_in(:,3),GLM_in(:,2),'*r')
% plot(LIS_in(:,4),LIS_in(:,3),'ob') 
% plot(linet(ind_inLINET,6),linet(ind_inLINET,5),'Xm','MarkerSize',12,'LineWidth',2)
% xlabel('Longitude [º]')
% ylabel('Latitude [º]')
% legend('LMA','LMA stations','GLM','LIS','LINET')
% axis equal
% 
% dist_ISS_LMA=deg2km(distance(LMA(:,4),LMA(:,5),iss_lon_lat(2),iss_lon_lat(1))); 
% dist_ISS_GLM=deg2km(distance(GLM_in(:,2),GLM_in(:,3),iss_lon_lat(2),iss_lon_lat(1))); 
% dist_ISS_LIS=deg2km(distance(LIS_in(:,3),LIS_in(:,4),iss_lon_lat(2),iss_lon_lat(1))); 
% 
% figure
% plot(LMA(:,1)-dt_dt,dist_ISS_LMA,'.k','MarkerSize',7)
% hold on
% plot(GLM_in(:,1)-dt_dt,dist_ISS_GLM,'*r')
% plot(LIS_in(:,1)-dt_dt,dist_ISS_LIS,'ob')
% xlim([MMIA_all(1,1) MMIA_all(end,1)])
% xlabel('UT-time (s)')
% ylabel('Distance from ISS (km)')
% grid on 
% legend('LMA','GLM','LIS')
% 
% 
% % 
% % 
% % 
% % figure(2021)
% % height_sub=0.26;
% % hAxis(1)=subplot(3,1,1);
% % % bar(LMA(:,1)-dt_lis,(LMA(:,6)/1000),60,'BarWidth',1,'FaceColor',[0.255 0.255 0.204],'EdgeColor',[0.5 0.5 0.5],'LineWidth',0.1);
% % 
% % stairs(max_LIS(:,1)-dt_lis,max_LIS(:,3));
% % hold on
% % stairs(max_GLM(:,1)-dt_lis,max_GLM(:,3),'r');
% % plot(linet(ind_inLINET,11)-dt_lis,0,'Xk','MarkerSize',12)
% % ylabel('Energy [fJ]')
% % set(gca, 'YGrid', 'off', 'XGrid', 'on')
% % xlim([MMIA_all(1,1) MMIA_all(end,1)])
% % pos1 = get( hAxis(1), 'Position' );
% % % pos1(2)=0.23; % Shift down.
% % pos1(4)=height_sub; % Increase height.
% % set( hAxis(1), 'Position', pos1) ;
% % set(gca,'xticklabel',[])
% % 
% % hAxis(2)=subplot(3,1,2); 
% % plot(MMIA_all(:,1),MMIA_all(:,4),'r')
% % hold on
% % plot(linet(ind_inLINET,11)-dt_lis,0,'Xk','MarkerSize',12)
% % %title 'PHOT 3 (777.4/ 5nm)'
% % ylabel('Energy [uWm^-2]')
% % set(gca, 'YGrid', 'off', 'XGrid', 'on')
% % xlim([MMIA_all(1,1) MMIA_all(end,1)])
% % pos2=get( hAxis(2), 'Position' );
% % pos2(4)=height_sub;
% % set(hAxis(2), 'Position', pos2) ;
% % set(gca,'xticklabel',[])
% % 
% % hAxis(3)=subplot(3,1,3); 
% % plot(MMIA_all(:,1),MMIA_all(:,2))
% % hold on
% % plot(linet(ind_inLINET,11)-dt_lis,0,'Xk','MarkerSize',12)
% % %title 'PHOT 1 (337.0/ 5nm)'
% % ylabel('Energy [uWm^-2]')
% % set(gca, 'YGrid', 'off', 'XGrid', 'on')
% % xlim([MMIA_all(1,1) MMIA_all(end,1)])
% % pos3=get( hAxis(3), 'Position' );
% % pos3(4)=height_sub;
% % set( hAxis(3), 'Position', pos3) ;
% % xlabel('Time [s]')
% % 
% 
% 
% 
% % 
% % 
% figure(2021)
% hAxis(1)=subplot(3,1,1);
% bar(LMA(:,1)-dt_glm,LMA(:,7)./((LMA(:,6))/1000),n_bins,'BarWidth',1,'FaceColor',[0.255 0.255 0.204],'EdgeColor',[0.5 0.5 0.5],'LineWidth',0.1);
% ylabel('Height [km]')
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% xlim([MMIA_all(1,1) MMIA_all(end,1)])
% 
% hAxis(2)=subplot(3,1,2); 
% plot(MMIA_all(:,1),MMIA_all(:,4),'r')
% %title 'PHOT 3 (777.4/ 5nm)'
% ylabel('Energy [uWm^-2]')
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% xlim([MMIA_all(1,1) MMIA_all(end,1)])
% 
% hAxis(3)=subplot(3,1,3); 
% plot(MMIA_all(:,1),MMIA_all(:,2))
% %title 'PHOT 1 (337.0/ 5nm)'
% ylabel('Energy [uWm^-2]')
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% xlim([MMIA_all(1,1) MMIA_all(end,1)])
% 
% 
% 
% 
% % figure(1010)
% % hAxis(1)=subplot(3,1,1);
% % plot(MMIA_all(:,1),MMIA_all(:,2))
% % %title 'PHOT 1 (337.0/ 5nm)'
% % xlabel('Time [s]')
% % % 
% % % 
% % hAxis(2)=subplot(3,1,2);
% % plot(MMIA_all(:,1),MMIA_all(:,3))
% % %title 'PHOT 2 (180-230nm)'
% % xlabel('Time [s]')
% % % 
% % % 
% % hAxis(3)=subplot(3,1,3); 
% % plot(MMIA_all(:,1),MMIA_all(:,4))
% % %title 'PHOT 3 (777.4/ 5nm)'
% % xlabel('Time [s]')
% % % % PHOT 1
% % % 
% 
% 
% figure
% hold on
% plot(LIS_all(ind_inLIS,4),LIS_all(ind_inLIS,3),'or') % SE GRAFICA LAS FUENTES AL INTERIOR DE ZONA DE COBERTURA LMA 
% plot(LMA(:,5),LMA(:,4),'.');
% plot(GLM(ind_inGLM,3),GLM(ind_inGLM,2),'*g');
% plot(linet(ind_inLINET,6),linet(ind_inLINET,5),'X','MarkerSize',11,'LineWidth',2)
% grid on
% xlabel('Longitude [º]')
% ylabel('Latitude [º]')
% legend('LIS','LMA','GLM','LINET')
% 
% figure 
% plot(MMIA_all(:,1),MMIA_all(:,2)/max_phot_mmia337)
% hold on
% plot( GLM(ind_inGLM,1),GLM(ind_inGLM,7)/max_ener_GLM,'*g')
% plot(LIS_original(ind_inLIS,1),LIS_original(ind_inLIS,5)/max_rad_LIS,'or')
% plot(LMA(:,1),LMA(:,6)/max_h_LMA,'.m');
% grid on
% xlabel('Time [s]')
% ylabel('Normalized')
% title('Before time corrected V4')
% if ~isempty(linet)
%     for li=1:length(linet(ind_inLINET))
%         if linet(ind_inLINET(li),8)>1 % Rayo IC
%             plot(linet(ind_inLINET(li),11),0.2,'Xm','MarkerSize',14)
%         elseif linet(ind_inLINET(li),9)>0
%             plot(linet(ind_inLINET(li),11),0,'+r','MarkerSize',14,'LineWidth',2)
%         else
%             plot(linet(ind_inLINET(li),11),0,'xk','MarkerSize',14,'LineWidth',2)
%         end
%     end
% end
% legend('MMIA777','GLM','LIS','LMA','Linet')
% dt_LIS_GLM=abs(dt_lis - dt_glm); % Tiempo de diferencia LIS y GLM
% % hold on
% % plot(MMIA_all(tresh_phot3(frames_tresh)),0.1,'xr')
% 
% figure
% plot(MMIA_all(:,1),MMIA_all(:,4)/max_phot_mmia777)
% hold on
% grid on
% plot(GLM(ind_inGLM,1)-dt_glm,GLM(ind_inGLM,7)/max_ener_GLM,'*g');
% plot(LIS_all(ind_inLIS,1)-dt_lis,LIS_all(ind_inLIS,5)/max_rad_LIS,'or') % SE GRAFICA LAS FUENTES AL INTERIOR DE ZONA DE COBERTURA LMA 
% plot(LMA(:,1)-dt_lis,LMA(:,6)/max_h_LMA,'.');
% title 'PHOT 3 (777.4/ 5nm)'
% xlabel('Time [s]')
% if ~isempty(linet)
%     for li=1:length(linet(ind_inLINET))
%         if linet(ind_inLINET(li),8)>1 % Rayo IC
%             plot(linet(ind_inLINET(li),11)-dt_lis,0.2,'Xm','MarkerSize',14)
%         elseif linet(ind_inLINET(li),9)>0
%             plot(linet(ind_inLINET(li),11)-dt_lis,0,'+r','MarkerSize',14,'LineWidth',2)
%         else
%             plot(linet(ind_inLINET(li),11)-dt_lis,0,'xk','MarkerSize',14,'LineWidth',2)
%         end
%     end
% end
% legend('MMIA777','GLM','LIS','LMA','Linet')
% title('After time corrected V4')
% 
% 
% 
% 
% 
% figure
% % hAx(1)=gca;
% % hAx(2)=axes('Position',hAx(1).Position,'XAxisLocation','top','YAxisLocation','right','color','none');
% % hold(hAx(1),'on')
% % plot(hAx(1),MMIA_all(:,1),MMIA_all(:,4)/max_phot_mmia777)
% plot(MMIA_all(:,1),MMIA_all(:,2)/max_phot_mmia337)
% hold on
% grid on
% plot(GLM(ind_inGLM,1)-dt_glm,GLM(ind_inGLM,7)/max_ener_GLM,'*g');
% plot(LIS_all(ind_inLIS,1)-dt_lis,LIS_all(ind_inLIS,5)/max_rad_LIS,'or') % SE GRAFICA LAS FUENTES AL INTERIOR DE ZONA DE COBERTURA LMA 
% plot(LMA(:,1)-dt_lis,LMA(:,6)/max_h_LMA,'.');
% xlabel('Time [s]')
% if ~isempty(linet)
%     for li=1:length(linet(ind_inLINET))
%         if linet(ind_inLINET(li),8)>1 % Rayo IC
%             plot(linet(ind_inLINET(li),11)-dt_lis,0.2,'Xm','MarkerSize',14)
%         elseif linet(ind_inLINET(li),9)>0
%             plot(linet(ind_inLINET(li),11)-dt_lis,0,'+r','MarkerSize',14,'LineWidth',2)
%         else
%             plot(linet(ind_inLINET(li),11)-dt_lis,0,'xk','MarkerSize',14,'LineWidth',2)
%         end
%     end
% end
% legend('MMIA337','GLM','LIS','LMA','Linet')
% title('After time corrected V4')



% Encuentra los datos de GLM más cercanos a LINET. 
%  [~, vidx_glm] = histc(linet(ind_inLINET,11), GLM(:,1));
%    if vidx_glm ~= 0
%      result_GLM = GLM(vidx_glm);
%    end
% 
% 
% 
% [~,idx] = sort(LIS_all(:,1)); % sort just the first column
% LIS_all = LIS_all(idx,:); 
% 
%  [~, vidx_lis] = histc(linet(ind_inLINET,11), LIS_all(:,1) );
%    if vidx_lis ~= 0
%      result_LIS = LIS_all(vidx_lis);
%    end

% if ~isempty(linet)
%     
%     dt_LMmia=linet(ind_inLINET,11)-t_mmia_max; 
%     dt_LGlm=linet(ind_inLINET,11)-result_GLM;
%     dt_LLis=linet(ind_inLINET,11)-result_LIS; 
%     
%     dt_glm=result_GLM(1,1)-t_mmia_max; 
%     
% end


% figure
% plot(MMIA_all(:,1),MMIA_all(:,2)/max_phot_mmia337)
% hold on
% plot(LIS_all(ind_inLIS,1)-dt_lis,LIS_all(ind_inLIS,5)/max_rad_LIS,'or') % SE GRAFICA LAS FUENTES AL INTERIOR DE ZONA DE COBERTURA LMA 
% plot(LMA(:,1)-dt_lis,LMA(:,6)/max_h_LMA,'.');
% plot(GLM(ind_inGLM,1)-dt_glm,GLM(ind_inGLM,7)/max_ener_GLM,'*g');
% title 'PHOT 1 (337/ 5nm)'
% xlabel('level 1 time (s)')
% legend('MMIA337','LIS','LMA','GLM')



% 
% figure
% plot(MMIA_all(:,1),MMIA_all(:,4)/max_phot_mmia777)
% hold on
% plot(LIS_all(ind_inLIS,1)-dt_lis,LIS_all(ind_inLIS,5)/max_rad_LIS,'or') % SE GRAFICA LAS FUENTES AL INTERIOR DE ZONA DE COBERTURA LMA 
% plot(LMA(:,1)-dt_lis,LMA(:,6)/max_h_LMA,'.');
% plot(GLM(ind_inGLM,1)-dt_glm,GLM(ind_inGLM,7)/max_ener_GLM,'*g');
% title 'PHOT 1 (337/ 5nm)'
% xlabel('level 1 time (s)')

% figure
% plot(MMIA_all(:,1),MMIA_all(:,4)/max_phot_mmia777)
% hold on
% plot(LIS_all(ind_inLIS,1)-dt_lis,LIS_all(ind_inLIS,5)/max_rad_LIS,'or') % SE GRAFICA LAS FUENTES AL INTERIOR DE ZONA DE COBERTURA LMA 
% plot(LMA(:,1),LMA(:,6)/max_h_LMA,'.');
% plot(GLM(ind_inGLM,1)-dt_lis,GLM(ind_inGLM,7)/max_ener_GLM,'*g');
% title 'PHOT 3 (777.4/ 5nm)'
% xlabel('level 1 time (s)')


% figure
% plot(MMIA_all(:,1),MMIA_all(:,3)/max_phot_mmia180)
% hold on
% plot(LIS_all(ind_inLIS,1)-dt_lis,LIS_all(ind_inLIS,5)/max_rad_LIS,'or') % SE GRAFICA LAS FUENTES AL INTERIOR DE ZONA DE COBERTURA LMA 
% plot(LMA(:,1)-dt_lis,LMA(:,6)/max_h_LMA,'.');
% plot(GLM(ind_inGLM,1)-dt_lis,GLM(ind_inGLM,7)/max_ener_GLM,'*g');
% title 'PHOT 3 (180-230nm)'
% xlabel('level 1 time (s)')





% 
% 
% 
% 
% %  
% % id_flashLIS=unique(LIS_all_inLMA(:,2)); 
% figure
% 
% % stairs(max_MMIA(:,1)+dt,max_MMIA(:,2)/max_phot_mmia777,'r')
% 
% plot(MMIA_all(:,1)+dt,MMIA_all(:,4)/max_phot_mmia777,'b')
% 
% hold on
% stairs(max_LIS(:,1),max_LIS(:,2)/max_rad_LIS,'r')
% stairs(max_GLM(:,1)-dt_LIS_GLM,max_GLM(:,2)/max_ener_GLM,'g')
% stairs(max_LMA(:,1),max_LMA(:,2)/max_h_LMA,'g')
% xlabel('Time [s]'); 
% ylabel('P.U') 
% figure
% subplot(4,1,1)
% plot(MMIA_all(:,1)+dt,MMIA_all(:,4)/max_phot_mmia777,'r')
% title 'PHOT 3 (777.4/ 5nm)'
% xlim([MMIA_all(1,1) MMIA_all(end,1)])
% 
% subplot(4,1,2)
% stairs(max_LIS(:,1),max_LIS(:,2)/max_rad_LIS,'g')
% title 'LIS (2-ms bins)'
% xlim([MMIA_all(1,1) MMIA_all(end,1)])
% 
% 
% subplot(4,1,3)
% stairs(max_GLM(:,1)-dt_LIS_GLM,max_GLM(:,2)/max_ener_GLM,'k')
% title 'GLM (2-ms bins)'
% xlim([MMIA_all(1,1) MMIA_all(end,1)])
% 
% subplot(4,1,4)
% plot(LMA(:,1),LMA(:,6)/max_h_LMA,'.b')
% title 'LMA'
% xlim([MMIA_all(1,1) MMIA_all(end,1)])
% 
% 
% 




% 

% xlabel('Longitude [degrees]')
% ylabel('Latitude [degrees]')
% 
% 
% 
% [goes,R] = geotiffread('goes_tiff2.tif');
% figure
% mapshow(goes,R);


% for i = 1:length(id_flashLIS)
%     f_id_LIS=find(LIS_all_inLMA(:,2) == id_flashLIS(i) );
%     LIS_id{i} =  LIS_all_inLMA(f_id_LIS,:);
% end
% 
% id_flashGLM=unique(GLM(:,4)); 
% for i = 1:length(id_flashGLM)
%     f_id_GLM=find(GLM(:,4) == id_flashGLM(i) );
%     GLM_id{i} =  GLM(f_id_GLM,:);
% end
% 
% 
% 

% 
% 
% figure(5)
% plot(LMA(:,5),LMA(:,4),'.r')
% hold on
% %Para separar por ID de agrupación LIS
% col=['g' 'b' 'r'];
% for i=1:length(id_flashLIS)
% plot(LIS_id{i}(:,4),LIS_id{i}(:,3),'o')
% hold on
% end
% 
% figure(5)
% hold on
% for i=1:length(id_flashGLM)
% plot(GLM_id{i}(:,3),GLM_id{i}(:,2),'*')
% hold on
% end
% xlabel('Longitude [º]')
% ylabel('Latitude [º]')
% grid on
% 
% 





% Para separar por ID de agrupación LIS
% col=['g' 'b' 'r'];
% for i=1:length(id_flashLIS)
% plot(LIS_id{i}(:,4),LIS_id{i}(:,3),'*')
% hold on
% end






% 
% figure(4)
% plot(MMIA_all(:,1)+dt2,MMIA_all(:,2))
% hold on
% plot(LIS_all(ind_in,1),LIS_all(ind_in,5),'or') % SE GRAFICA LAS FUENTES AL INTERIOR DE ZONA DE COBERTURA LMA 
% plot(LMA(:,1),LMA(:,6)/1000,'.');
% title 'PHOT 1 (337.0/ 5nm)'
% xlabel('level 0 time (s)')
% 
% 
% figure(5)
% plot(MMIA_all(:,1)+dt2,MMIA_all(:,3))
% hold on
% plot(LIS_all(ind_in,1),LIS_all(ind_in,5),'or') % SE GRAFICA LAS FUENTES AL INTERIOR DE ZONA DE COBERTURA LMA 
% plot(LMA(:,1),LMA(:,6)/1000,'.');
% title 'PHOT 2 (180-230nm)'
% xlabel('level 0 time (s)')
% 
% 
% figure(6)
% plot(MMIA_all(:,1)+dt2,MMIA_all(:,4))
% hold on
% plot(LIS_all(ind_in,1),LIS_all(ind_in,5),'or') % SE GRAFICA LAS FUENTES AL INTERIOR DE ZONA DE COBERTURA LMA 
% plot(LMA(:,1),LMA(:,6)/1000,'.');
% title 'PHOT 3 (777.4/ 5nm)'
% xlabel('level 0 time (s)')
% 



% 
% subplot(3,1,1)
% plot(h_mmia+dt,PHOT1Data_all)
% hold on
% plot(LIS_all(ind_in,1),LIS_all(ind_in,5),'or')%referencia tiempo
% title 'PHOT 1 (337.0/ 5nm)'
% xlabel('level 0 time (s)')
% 
% %xlabel('level 0 time (s)')
% subplot(3,1,2)
% plot(h_mmia+dt,PHOT2Data_all)
% hold on
% plot(LIS_all(ind_in,1),LIS_all(ind_in,5),'or')%referencia tiempo
% title 'PHOT 2 (180-230nm)'
% %xlabel('level 0 time (s)')
% subplot(3,1,3)
% plot(h_mmia+dt,PHOT3Data_all)
% hold on
% plot(LIS_all(ind_in,1),LIS_all(ind_in,5),'or')%referencia tiempo
% title 'PHOT 3 (777.4/ 5nm)'
% xlabel('level 0 time (s)')
% 
% figure(10)
% plot(h_mmia+dt2,PHOT3Data_all)
% hold on
% plot(LIS_all(ind_in,1),LIS_all(ind_in,5),'or')%referencia tiempo
% title 'PHOT 3 (777.4/ 5nm)'
% xlabel('level 0 time (s)')
% 
% figure(11)
% plot(h_mmia+dt,PHOT3Data_all)
% hold on
% % plot(LIS_all(ind_in,1),LIS_all(ind_in,5),'or')%referencia tiempo
% title 'PHOT 3 (777.4/ 5nm)'
% xlabel('level 0 time (s)')
% 
% 
% figure(12)
% plot(h_mmia+dt,PHOT2Data_all)
% hold on
% plot(LIS_all(ind_in,1),LIS_all(ind_in,5),'or')%referencia tiempo
% title 'PHOT 2 (180-230nm)'
% xlabel('level 0 time (s)')
% 
% 
% figure(13)
% plot(h_mmia+dt,PHOT1Data_all)
% hold on
% plot(LIS_all(ind_in,1),LIS_all(ind_in,5),'or')%referencia tiempo
% title 'PHOT 1 (337.0/ 5nm)'
% xlabel('level 0 time (s)')
% 
% 
% %Photometers Plot
% frame=1; 
% figure(3)
% subplot(3,1,1)
% plot(tvector_L0,PHOT1Data(frame,1:length(tvector_L0)))
% hold on
% plot(T_Trigg,min(PHOT1Data(frame,1:length(tvector_L0))),'ro')%referencia tiempo
% title 'PHOT 1 (337.0/ 5nm)'
% %xlabel('level 0 time (s)')
% subplot(3,1,2)
% plot(tvector_L0,PHOT2Data(frame,1:length(tvector_L0)))
% hold on
% plot(T_Trigg,min(PHOT2Data(frame,1:length(tvector_L0))),'ro')%referencia tiempo
% title 'PHOT 2 (180-230nm)'
% %xlabel('level 0 time (s)')
% subplot(3,1,3)
% plot(tvector_L0,PHOT3Data(frame,1:length(tvector_L0)))
% hold on
% plot(T_Trigg,min(PHOT3Data(frame,1:length(tvector_L0))),'ro')%referencia tiempo
% title 'PHOT 1 (777.4/ 5nm)'
% xlabel('level 0 time (s)')

% Extraemos CHU1_CHU2 si lo hay
%CHU1/2_META are the sums(1026+1056=2082 )
%CHU1/2Data Have IMAGE information

%CHU 1 centroid coordinates: (Column 523.92 , Row 512.99)
%CHU 2 centroid coordinates: (Column 530.56 , Row 515.69)


% LMA sources

% xticklabels(time_label)


% [v_LMA_hmax pos_LMA_hmax]=max(LMA(:,5)); 
% dt3=abs(t_mmia_max-LMA(pos_LMA_hmax,1)); 
% dt3=abs(t_lis_max-LMA(3,1)); 

% figure(13)
% plot(h_mmia-dt3,PHOT3Data_all)
% hold on
% plot(LIS_all(ind_in,1),LIS_all(ind_in,5),'or')%referencia tiempo
% title 'PHOT 1 (777.4/ 5nm)'
% xlabel('level 0 time (s)')
% plot(LMA(:,1),LMA(:,6)/1000,'.');


%%
% figure(2022)
% height_sub=0.2;
% hAxis(1)=subplot(4,1,1); 
% bar(LMA(:,1)-dt_dt,(LMA(:,6))/1000,n_bins,'FaceColor',[0.255 0.255 0.204]); %,'BarWidth',1
% hold on
% if ~isempty(linet)
%     dt_lma_linet=LMA(1,1)- linet(ind_inLINET(1),11)
%     for li=1:length(linet(ind_inLINET))
%         if linet(ind_inLINET(li),8)>1 % Rayo IC
%             plot(linet(ind_inLINET(li),11)-dt_dt,linet(ind_inLINET(li),7),'*m','MarkerSize',12)
%         else
%             if linet(ind_inLINET(li),9)>0 % +CG
%                 plot(linet(ind_inLINET(li),11)-dt_dt,0,'+m','MarkerSize',12,'LineWidth',1)
%             else
%                 plot(linet(ind_inLINET(li),11)-dt_dt,0,'xk','MarkerSize',12,'LineWidth',1)
%             end
%             
%         end
%     end
% end
% ylabel('Height [km]')
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% xlim([MMIA_all(1,1) MMIA_all(end,1)])
% pos1 = get( hAxis(1), 'Position' );
% pos1(4)=height_sub; % Increase height.
% set( hAxis(1), 'Position', pos1) ;
% set(gca,'xticklabel',[])
% 
% hAxis(2)=subplot(4,1,2); 
% plot(MMIA_all(:,1),MMIA_all(:,4),'r')
% hold on
% if ~isempty(linet)
%     for li=1:length(linet(ind_inLINET))
%         if linet(ind_inLINET(li),8)>1 % Rayo IC
%             plot(linet(ind_inLINET(li),11)-dt_dt,linet(ind_inLINET(li),10),'*m','MarkerSize',12)
%         else
%             if linet(ind_inLINET(li),9)>0 % +CG
%                 plot(linet(ind_inLINET(li),11)-dt_dt,0,'+m','MarkerSize',12,'LineWidth',1)
%             else
%                 plot(linet(ind_inLINET(li),11)-dt_dt,0,'xk','MarkerSize',12,'LineWidth',1)
%             end
%         end
%     end
% end
% %title 'PHOT 1 (337.0/ 5nm)'
% ylabel('Energy [uWm^-^2]') 
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% xlim([MMIA_all(1,1) MMIA_all(end,1)])
% pos2=get(hAxis(2), 'Position' );
% pos2(4)=height_sub;
% set(hAxis(2), 'Position', pos2) ;
% set(gca,'xticklabel',[])
% legend('Phot.777','Linet')
% 
% hAxis(3)=subplot(4,1,3); 
% plot(MMIA_all(:,1),MMIA_all(:,2))
% hold on
% if ~isempty(linet)
%     for li=1:length(linet(ind_inLINET))
%         if linet(ind_inLINET(li),8)>1 % Rayo IC
%             plot(linet(ind_inLINET(li),11)-dt_dt,linet(ind_inLINET(li),10),'*m','MarkerSize',12)
%         else
%             if linet(ind_inLINET(li),9)>0 % +CG
%                 plot(linet(ind_inLINET(li),11)-dt_dt,0,'+m','MarkerSize',12,'LineWidth',1)
%             else
%                 plot(linet(ind_inLINET(li),11)-dt_dt,0,'xk','MarkerSize',12,'LineWidth',1)
%             end
%         end
%     end
% end
% %title 'PHOT 1 (337.0/ 5nm)'
% ylabel('Energy [uWm^-^2]') 
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% xlim([MMIA_all(1,1) MMIA_all(end,1)])
% pos3=get( hAxis(3), 'Position' );
% pos3(4)=height_sub;
% set( hAxis(3), 'Position', pos3) ;
% set(gca,'xticklabel',[])
% legend('Phot.337','Linet')
% 
% 
% hAxis(4)=subplot(4,1,4); 
% plot(int_GLM(:,1)-dt_glm,(int_GLM(:,2))/(max(int_GLM(:,2))),'r','LineWidth',2);
% hold on
% plot(int_LIS(:,1)-dt_dt,(int_LIS(:,2))/(max(int_LIS(:,2))));
% if ~isempty(linet)
%     for li=1:length(linet(ind_inLINET))
%         if linet(ind_inLINET(li),8)>1 % Rayo IC
%             plot(linet(ind_inLINET(li),11)-dt_dt,0,'*m','MarkerSize',12)
%         else
%             if linet(ind_inLINET(li),9)>0 % +CG
%                 plot(linet(ind_inLINET(li),11)-dt_dt,0,'+m','MarkerSize',12,'LineWidth',1)
%             else
%                 plot(linet(ind_inLINET(li),11)-dt_dt,0,'xk','MarkerSize',12,'LineWidth',1)
%             end
%         end
%     end
% end
% %title 'PHOT 3 (777.4/ 5nm)'
% ylabel('Normalized') 
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% xlim([MMIA_all(1,1) MMIA_all(end,1)])
% pos4=get( hAxis(4), 'Position' );
% pos4(4)=height_sub;
% set(hAxis(4), 'Position',pos4) ;
% legend('GLM','LIS','Linet')
% xlabel('Time [s]')
% linkaxes([hAxis(1),hAxis(2),hAxis(3),hAxis(4)],'x');

% ACTIVAR PARA PLOTS DEFINITIVOS
% Freedman-Diaconis rule
% bin_=2*(iqr( LMA(:,7) ) )* (length(LMA(:,1)))^(-1/3);
% n_bins=(max(LMA(:,7))-min(LMA(:,7)) )/bin_; 
% figure(2020)
% height_sub=0.28;
% hAxis(1)=subplot(3,1,1); 
% bar(LMA(:,1)-dt_glm,(LMA(:,7)),n_bins,'FaceColor',[0.255 0.255 0.204]); %,'BarWidth',1,
% % bar(LMA(:,1)-dt_lis,(LMA(:,6))/1000,n_bins,'BarWidth',1,'FaceColor',[0.255 0.255 0.204],'EdgeColor',[0.5 0.5 0.5],'LineWidth',0.1);
% % area(LMA(:,1)-dt_lis,abs(LMA(:,7)));
% hold on
% if ~isempty(linet)
%     for li=1:length(linet(ind_inLINET))
%         if linet(ind_inLINET(li),8)>1 % Rayo IC
%             plot(linet(ind_inLINET(li),11)-dt_glm,linet(ind_inLINET(li),7),'*m','MarkerSize',12)
%         else
%             plot(linet(ind_inLINET(li),11)-dt_glm,0,'xk','MarkerSize',14,'LineWidth',2)
%         end
%     end
% end
% ylabel('Power [dBW]')
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% xlim([MMIA_all(1,1) MMIA_all(end,1)])
% pos1 = get( hAxis(1), 'Position' );
% % pos1(2)=0.23; % Shift down.
% pos1(4)=height_sub; % Increase height.
% set( hAxis(1), 'Position', pos1) ;
% set(gca,'xticklabel',[])
% 
% hAxis(2)=subplot(3,1,2); 
% plot(MMIA_all(:,1),MMIA_all(:,4),'r')
% hold on
% if ~isempty(linet)
%     for li=1:length(linet(ind_inLINET))
%         if linet(ind_inLINET(li),8)>1 % Rayo IC
%             plot(linet(ind_inLINET(li),11)-dt_glm,linet(ind_inLINET(li),7),'*m','MarkerSize',12)
%         else
%             plot(linet(ind_inLINET(li),11)-dt_glm,0,'xk','MarkerSize',14,'LineWidth',2)
%         end
%     end
% end
% %title 'PHOT 3 (777.4/ 5nm)'
% ylabel('Energy [uWm^-^2]') 
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% xlim([MMIA_all(1,1) MMIA_all(end,1)])
% pos2=get( hAxis(2), 'Position' );
% pos2(4)=height_sub;
% set(hAxis(2), 'Position', pos2) ;
% set(gca,'xticklabel',[])
% 
% hAxis(3)=subplot(3,1,3); 
% plot(MMIA_all(:,1),MMIA_all(:,2))
% hold on
% if ~isempty(linet)
%     for li=1:length(linet(ind_inLINET))
%         if linet(ind_inLINET(li),8)>1 % Rayo IC
%             plot(linet(ind_inLINET(li),11)-dt_glm,linet(ind_inLINET(li),10),'*m','MarkerSize',12)
%         else
%             plot(linet(ind_inLINET(li),11)-dt_glm,0,'xk','MarkerSize',14,'LineWidth',2)
%         end
%     end
% end
% %title 'PHOT 1 (337.0/ 5nm)'
% ylabel('[uWm^-^2]') 
% set(gca, 'YGrid', 'off', 'XGrid', 'on')
% xlim([MMIA_all(1,1) MMIA_all(end,1)])
% pos3=get( hAxis(3), 'Position' );
% pos3(4)=height_sub;
% set( hAxis(3), 'Position', pos3) ;
% xlabel('Time [s]')
% linkaxes([hAxis(1),hAxis(2),hAxis(3)],'x');
% legend('Phot.337','Linet')
% 
% 

