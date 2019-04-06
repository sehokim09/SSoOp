% Monte Carlo Campaign

close all, clear all
% ... Initialization and Setup
BasePath = './Baseline';
Nrun = 1;
MaxCore = 1;

Nbatch = ceil(Nrun/MaxCore);

% ... Simulation Setup
Duration = 6048890.4; % 67-Nodal day cycle
DTSIM = 0.1; % Step size
DTOUT = 0.1; % File Output Interval

% lonRef = -180; % Longitude of 1km Grid
% latMax = 80; % Maximum Latitude of 1km Grid
% latMin = -80; % Minimum Latitude of 1km Grid

Reflectivity = 0.1;

% ... Campaign Loop
for Ibatch = 1:Nbatch
    Ncore = min(MaxCore,Nrun-(Ibatch-1)*MaxCore);
    BatchPath = '';
    for Icore = 1:Ncore
        Irun = (Ibatch-1)*MaxCore + Icore;

        % Copy Base folder to Run folder
        RunPath = sprintf('Run%03d',Irun);
        system(['rm -rf ',RunPath]);
        system(['cp -r ',BasePath,' ',RunPath]);
        BatchPath = [BatchPath,' ',RunPath];

        % Modify Input Files
        s = sprintf('%f   %1.1f                    !  Sim Duration, Step Size [sec]\n', Duration, DTSIM);
        OverwriteLineInFile([RunPath,'/Inp_Sim.txt'],4,s);
               
        s = sprintf('%1.1f                             !  File Output Interval [sec]\n', DTOUT);
        OverwriteLineInFile([RunPath,'/Inp_Sim.txt'],5,s);

%         s = sprintf('%d                             !  Longitude Reference [deg]\n', lonRef);
%         OverwriteLineInFile([RunPath,'/Inp_SoOp.txt'],3,s);
%         lonRef = lonRef + 60;
%         
%         s = sprintf('%d                           !  Latitude Maximum [deg]\n', latMax);
%         OverwriteLineInFile([RunPath,'/Inp_SoOp.txt'],4,s);
%         
%         s = sprintf('%d                            !  Latitude Minimum [deg]\n', latMin);
%         OverwriteLineInFile([RunPath,'/Inp_SoOp.txt'],5,s);
%         
      %  Reflectivity = Reflectivity + 0.1;
      %  Inp_Reflectivity(Irun) = Reflectivity;
        s = sprintf('%1.1f                           !  Reflectivity\n', Reflectivity);
        OverwriteLineInFile([RunPath,'/Inp_SoOp.txt'],7,s); 
    end

    % Perform Runs
    fprintf('Starting Runs %d - %d\n',(Ibatch-1)*MaxCore + 1, (Ibatch-1)*MaxCore+Ncore);
    system(['./BatchRun.sh ',BatchPath]);

    % Process Batch Results
    for Icore = 1:Ncore
        Irun = (Ibatch-1)*MaxCore + Icore;
        RunPath = sprintf('Run%03d',Irun);

        % Process Run in Matlab
        LoadString = [RunPath,'/MUOS.42'];
        load(LoadString);
        LoadString = [RunPath,'/Orbcomm.42'];
        load(LoadString);
        LoadString = [RunPath,'/GPS.42'];
        load(LoadString);
%         AbsTime(:,Irun) = MUOS(:,1);
%         Lat_MUOS(:,Irun) = MUOS(:,2);
%         Lon_MUOS(:,Irun) = MUOS(:,3);
%         SNR_D_MUOS(:,Irun) = MUOS(:,4);
%         SNR_R_MUOS(:,Irun) = MUOS(:,5);
%         az_MUOS(:,Irun) = MUOS(:,6);
%         
%         
%         Lat_Orbcomm(:,Irun) = Orbcomm(:,2);
%         Lon_Orbcomm(:,Irun) = Orbcomm(:,3);
%         SNR_D_Orbcomm(:,Irun) = Orbcomm(:,4);
%         SNR_R_Orbcomm(:,Irun) = Orbcomm(:,5);
%         IncAng_Orbcomm(:,Irun) = Orbcomm(:,6);

        % Save only selected *.42 files
        SaveList = {'MUOS','Orbcomm','GPS'};
        Nsave = length(SaveList);
        for Isave = 1:Nsave
          system(['mv ',RunPath,'/',SaveList{Isave},'.42 ',RunPath,'/',SaveList{Isave},'.tmp']);
        end
        system(['rm ',RunPath,'/*.42']);
        for Isave = 1:Nsave
          system(['mv ',RunPath,'/',SaveList{Isave},'.tmp ',RunPath,'/',SaveList{Isave},'.42']);
        end

        % Delete unaltered *.txt files
        system(['rm ',RunPath,'/Inp_Cmd.txt']);
        system(['rm ',RunPath,'/Inp_FOV.txt']);
        system(['rm ',RunPath,'/Inp_Graphics.txt']);
        system(['rm ',RunPath,'/Inp_IPC.txt']);
        system(['rm ',RunPath,'/Inp_Region.txt']);
        system(['rm ',RunPath,'/Inp_Sim.txt']);
        system(['rm ',RunPath,'/Inp_TDRS.txt']);
        system(['rm ',RunPath,'/Inp_NOS3.txt']);
        system(['rm ',RunPath,'/Inp_TLE.txt']);
    end

end
% 
% % ... Post-processing
% land_MUOS = landmask(Lat_MUOS(:,1), Lon_MUOS(:,1));
% 
% Lat_MUOS_land = Lat_MUOS(land_MUOS);
% Lon_MUOS_land = Lon_MUOS(land_MUOS);
% SNR_D_MUOS_land = SNR_D_MUOS(land_MUOS,1);
% SNR_R_MUOS_land1 = SNR_R_MUOS(land_MUOS,1);
% SNR_R_MUOS_land2 = SNR_R_MUOS(land_MUOS,2);
% SNR_R_MUOS_land3 = SNR_R_MUOS(land_MUOS,3);
% SNR_R_MUOS_land4 = SNR_R_MUOS(land_MUOS,4);
% SNR_R_MUOS_land5 = SNR_R_MUOS(land_MUOS,5);
% 
% land_Orbcomm = landmask(Lat_Orbcomm(:,1), Lon_Orbcomm(:,1));
% 
% Lat_Orbcomm_land = Lat_Orbcomm(land_Orbcomm);
% Lon_Orbcomm_land = Lon_Orbcomm(land_Orbcomm);
% SNR_D_Orbcomm_land = SNR_D_Orbcomm(land_Orbcomm,1);
% SNR_R_Orbcomm_land1 = SNR_R_Orbcomm(land_Orbcomm,1);
% SNR_R_Orbcomm_land2 = SNR_R_Orbcomm(land_Orbcomm,2);
% SNR_R_Orbcomm_land3 = SNR_R_Orbcomm(land_Orbcomm,3);
% SNR_R_Orbcomm_land4 = SNR_R_Orbcomm(land_Orbcomm,4);
% SNR_R_Orbcomm_land5 = SNR_R_Orbcomm(land_Orbcomm,5);
% 
% % SNR on Ground Track
% WorldMap = imread('BigBlueMarble.ppm');
% 
% figure(1)
% image([-180 180], [-90 90], flipud(WorldMap))
% set(gca,'YDir','normal')
% hold on; 
% scatter(Lon_MUOS_land, Lat_MUOS_land, [], SNR_D_MUOS_land, '.');
% cb = colorbar; % create colorbar
% title('MUOS SNR_D [dB]')
% grid on
% xlabel('Longitude [deg]')
% ylabel('Latitude [deg]')
% set(gcf, 'Position',  [200, 200, 720, 360])
% 
% figure(2)
% image([-180 180], [-90 90], flipud(WorldMap))
% set(gca,'YDir','normal')
% hold on; 
% scatter(Lon_MUOS_land, Lat_MUOS_land, [], SNR_R_MUOS_land1, '.');
% cb = colorbar;
% title('MUOS SNR_R [dB]')
% grid on
% xlabel('Longitude [deg]')
% ylabel('Latitude [deg]')
% set(gcf, 'Position',  [200, 200, 720, 360])
% 
% figure(3)
% image([-180 180], [-90 90], flipud(WorldMap))
% set(gca,'YDir','normal')
% hold on; 
% scatter(Lon_Orbcomm_land, Lat_Orbcomm_land, [], SNR_D_Orbcomm_land, '.');
% cb = colorbar;
% title('Orbcomm SNR_D [dB]')
% grid on
% xlabel('Longitude [deg]')
% ylabel('Latitude [deg]')
% set(gcf, 'Position',  [200, 200, 720, 360])
% 
% figure(4)
% image([-180 180], [-90 90], flipud(WorldMap))
% set(gca,'YDir','normal')
% hold on; 
% scatter(Lon_Orbcomm_land, Lat_Orbcomm_land, [], SNR_R_Orbcomm_land1, '.');
% cb = colorbar;
% title('Orbcomm SNR_R [dB]')
% grid on
% xlabel('Longitude [deg]')
% ylabel('Latitude [deg]')
% set(gcf, 'Position',  [200, 200, 720, 360])
% 
% % Incidence Angle PDF
% figure(5)
% histogram(Az_MUOS,'Normalization','probability')
% title('PDF of MUOS Specular Points Incidence Angle')
% xlabel('Incidence Angle [deg]')
% ylabel('Probability Density')
% 
% figure(6)
% histogram(IncAng_Orbcomm,'Normalization','probability')
% title('PDF of Orbcomm Specular Points Incidence Angle')
% xlabel('Incidence Angle [deg]')
% ylabel('Probability Density')
% 
% % SNR v.s. Reflectivity
% SNR_R_MUOS_m(1) = mean(SNR_R_MUOS_land1);
% SNR_R_MUOS_m(2) = mean(SNR_R_MUOS_land2);
% SNR_R_MUOS_m(3) = mean(SNR_R_MUOS_land3);
% SNR_R_MUOS_m(4) = mean(SNR_R_MUOS_land4);
% SNR_R_MUOS_m(5) = mean(SNR_R_MUOS_land5);
% 
% SNR_R_Orbcomm_m(1) = mean(SNR_R_Orbcomm_land1);
% SNR_R_Orbcomm_m(2) = mean(SNR_R_Orbcomm_land2);
% SNR_R_Orbcomm_m(3) = mean(SNR_R_Orbcomm_land3);
% SNR_R_Orbcomm_m(4) = mean(SNR_R_Orbcomm_land4);
% SNR_R_Orbcomm_m(5) = mean(SNR_R_Orbcomm_land5);
% 
% figure(7)
% plot(Inp_Reflectivity, SNR_R_MUOS_m, '-o')
% title('MUOS SNR_R vs Reflectivity')
% grid on
% xlabel('Reflectivity')
% ylabel('Mean SNR_R [dB]')
% 
% figure(8)
% plot(Inp_Reflectivity, SNR_R_Orbcomm_m, '-o')
% title('Orbcomm SNR_R vs Reflectivity')
% grid on
% xlabel('Reflectivity')
% ylabel('Mean SNR_R [dB]')

% ... Cleanup

