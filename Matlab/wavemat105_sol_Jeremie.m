function [soilCommonStateParams,soilInnerStateParams,soilOuterStateParams] =...
         wavemat105_sol_Jeremie(cauliConsParams,cauliStateIni,climateConsParams,fileSettings,leekConsParams,...
         leekStateIni,managementSettings,simulationSettings,soilCommonStateParams,soilConsParams,soilInnerStateParams,...
         soilOuterStateParams,parallelFlag)
close all force
if ~parallelFlag
    WAITBAR = waitbar(0,"Initializing parameters and initial conditions",'Name',"ECOFERT Simulation",...
        'CreateCancelBtn','setappdata(gcbf,"canceling",1)');
    set(WAITBAR,'units','normalized','Position',[0.4 0.8 0.2 0.15])
    setappdata(WAITBAR,'canceling',0);
end
%% CHECK RESULTS FOLDER 

ResultsPath = fileSettings.ResultsPath;
% check if directory already exists and choose what to do 
isDir = (7==exist(fileSettings.ResultsPath,'dir')); % check if folder exist
isEmpty = size(dir(fileSettings.ResultsPath),1)<=2; % check if folder is empty
if isDir
    if isEmpty
        disp(['Specified folder is still empty, thus the inputfiles are saved in ',...
            ResultsPath])
    else
        overrideFolderFlag = questdlg('The specified folder already exists. Would you like to override the current folder and its results?', ...
            'Save Folder', ...
            'Yes','No','No');
        switch overrideFolderFlag
            case 'No'
                % find the largest number with the same folder name
                folder_idx = strfind(ResultsPath,'\');
                folder_idx = folder_idx(end-1);
                folderPath = ResultsPath(1:folder_idx); % path where all different folders with inputfiles are stored
                saveFolder = ResultsPath(folder_idx+1:end-1);
                d = dir(folderPath);
                folderNames = {d(3:end).name}; % all the names of the folders containing inputfiles
                folderNumber = 1;
                for name = folderNames
                    matchingFolder = strfind(lower(name{1}),lower(saveFolder));
                    if matchingFolder
                        if strcmpi(name{1},saveFolder)
                            tmp_folderNumber = str2double(name{1}(end));
                        else
                            tmp_folderNumber = str2double(name{1}(length(saveFolder)+1:end));
                        end
                        if ~isempty(tmp_folderNumber)
                            if tmp_folderNumber>=folderNumber 
                                folderNumber = tmp_folderNumber+1;
                            end
                        end
                    end
                end
                ResultsPath = [ResultsPath(1:end-2),num2str(folderNumber),'\'];
                mkdir(ResultsPath)
                disp(['Results will be saved in ', ResultsPath])
            case 'Yes'
                
                % delete all files in de existing folder
                d = dir(ResultsPath);
                d = d(3:end);
                
                for i = 1:length(d)
                    baseFileName = d(i).name;
                    fullFileName = fullfile(ResultsPath, baseFileName);
                    if d(i).isdir
                        rmdir(fullFileName,'s')
                    else
                        delete(fullFileName)
                    end
                end
                disp(['Results in ', ResultsPath,' will be overwritten'])
        end
    end    
else 
    mkdir(ResultsPath)
    disp(['Folder is made and results will be saved in ',ResultsPath])
end 

%% MODEL INPUT AND INITIALIZATION
% In this module the inputfiles created with initialParameters.m are loaded

%Load climate data, convert to useful format and add necessary information
climateState = struct();
climateState.climate = dlmread(fileSettings.ClimatePath);
climateState = convert_climate_to_matlab(climateState);
climateState =calc_ET0(climateConsParams,climateState,leekConsParams); % add potential evapotranspiration to the climate data
climateState = convert_climate_to_daily(climateState,soilConsParams); % convert hourly data to daily data

% convert days after startdate to julian date for harvest and plant dates
managementSettings.pDateLeek = managementSettings.pDateLeek + simulationSettings.tstart;
managementSettings.pDateCauli = managementSettings.pDateCauli + simulationSettings.tstart;
managementSettings.hDateLeek = managementSettings.hDateLeek + simulationSettings.tstart;
managementSettings.hDateCauli = managementSettings.hDateCauli + simulationSettings.tstart;

%Also possible to define boundary conditions via in_climatic_data
[soilConsParams,soilInnerStateParams] = In_Boundary_conditions(2,climateState,...
    soilConsParams,soilInnerStateParams);
[~,soilOuterStateParams] = In_Boundary_conditions(2,climateState,soilConsParams,...
    soilOuterStateParams);

% soil organic matter
soilOuterStateParams.om_appl = om_applic(managementSettings,soilOuterStateParams);

soilInnerStateParams.soil_om = soilConsParams.soil_om;
soilOuterStateParams.soil_om = soilConsParams.soil_om;

% Set initial pressure head profile in inner and outer cilinder of quasi 2D
% WAVE
[soilInnerStateParams,soilOuterStateParams] = In_initial_pressure_head(simulationSettings,...
                                                                                soilConsParams,...
                                                                                soilInnerStateParams,...
                                                                                soilOuterStateParams);
if simulationSettings.simNitroFlag
    % total organic carbon concentration in the innder cilinder
    soilInnerStateParams.tcorgs = sum((soilInnerStateParams.soil_om(:,1)) +...
                                       soilInnerStateParams.soil_om(:,3) +...
                                       soilInnerStateParams.soil_om(:,5));
    % total organic nitrogen concentration in the inner cilinder
    soilInnerStateParams.tnorgs = sum((soilInnerStateParams.soil_om(:,2)) +...
                                       soilInnerStateParams.soil_om(:,4) +...
                                       soilInnerStateParams.soil_om(:,6));
    % total organic carbon concentration in the outer cilinder
    soilOuterStateParams.tcorgs = sum((soilInnerStateParams.soil_om(:,1)) +...
                                       soilInnerStateParams.soil_om(:,3) +...
                                       soilInnerStateParams.soil_om(:,5));
    % total organic nitrogen concetration in the outer cilinder
    soilOuterStateParams.tnorgs = sum((soilInnerStateParams.soil_om(:,2)) +...
                                       soilInnerStateParams.soil_om(:,4) +...
                                       soilInnerStateParams.soil_om(:,6));
end

% Set initial water content stock in the soil
soilInnerStateParams.dt = simulationSettings.dt_start;
soilOuterStateParams.dt = simulationSettings.dt_start;
soilInnerStateParams = calc_stock(...
    simulationSettings,soilConsParams,soilInnerStateParams,0);
soilOuterStateParams = calc_stock(...
    simulationSettings,soilConsParams,soilOuterStateParams,0);
soilInnerStateParams.stock_max=sum(moist_ret(zeros(1,simulationSettings.nComp),...
                        soilInnerStateParams,soilConsParams,0)*simulationSettings.dx);
soilOuterStateParams.stock_max=sum(moist_ret(zeros(1,simulationSettings.nComp),...
                        soilInnerStateParams,soilConsParams,0)*simulationSettings.dx);                    
soilCommonStateParams.water_storage = soilInnerStateParams.stock_initial;

        
plant_date = sort([managementSettings.pDateCauli;managementSettings.pDateLeek]);
harvest_date = sort([managementSettings.hDateCauli;managementSettings.hDateLeek]);

    
%% MODEL SIMULATION
calculationtime_day = [];
while simulationSettings.t < simulationSettings.tmax
    tic
    if simulationSettings.t-simulationSettings.tstart == 91
        disp('')
    end
    % calculate hydraulic properties based on the Weynats pedo transfer function
    if simulationSettings.PTFWeynants
        soilConsParams = PTFWeynants(simulationSettings,soilCommonStateParams,soilConsParams);
    end
    
    progress = (simulationSettings.t-simulationSettings.tstart)/...
        (simulationSettings.tmax-simulationSettings.tstart)*100;  % calculate the progress of the simulation
    day = simulationSettings.t-simulationSettings.tstart +1;
    % determine the number of the rotation and whether soil AND plant
    % should be simulated or only soil
    simulationSettings = find_cropnumber(managementSettings,simulationSettings);

    %% PLANT MODULE
    
    if simulationSettings.flag_double_sim
        %% UPDATE WAITBAR
        
        if ~parallelFlag
            waitbar(progress/100,WAITBAR,[string(['simulating day ',num2str(day)]);strcat(num2str(progress),'% of the simulation is completed');...
                "simulating plant processes"],'Name',"ECOFERT simulation",'height',500)
            if getappdata(WAITBAR,'canceling')
                delete(WAITBAR)
                error(['The simulation has been stopped manually the results or',...
                    ' progress are stored in ',char(ResultsPath)])
            end
        end
        
        if simulationSettings.t == plant_date(simulationSettings.ncrop)      
            if strcmp(managementSettings.croptype(simulationSettings.ncrop),"Cauli")
                % Initialize cauliflower growth
                
                cropStateParams = cauliStateIni;
                cropConsParams = cauliConsParams;

            elseif strcmp(managementSettings.croptype(simulationSettings.ncrop),"Leek")
                % Initialize leek growth
                
                cropStateParams = leekStateIni;    
                cropConsParams = leekConsParams;
            end
        end
        
        %% SELECT CROP FOR SIMULATION
        if  strcmp(managementSettings.croptype(simulationSettings.ncrop),"Cauli")
            
            %Run cauli simulation for one day
            [cropStateParams,soilInnerStateParams,soilOuterStateParams] =...
                                            run_model_cauli(cauliConsParams,...
                                            climateConsParams,climateState,...
                                            cropStateParams,fileSettings,...
                                            managementSettings,...
                                            simulationSettings,soilCommonStateParams,...
                                            soilInnerStateParams,soilOuterStateParams);
            
            if any(managementSettings.hDateCauli == simulationSettings.t)
                date = datevec(simulationSettings.t);
                if date(1) == 1971
                    disp('')
                end
                cropStateParams = dataLog(cropStateParams, "Cauli", fileSettings, "Force");
            else
                cropStateParams = dataLog(cropStateParams, "Cauli", fileSettings);
            end
            
        elseif strcmp(managementSettings.croptype(simulationSettings.ncrop),"Leek")
            
            % run leek simulation model for 1 day
            [cropStateParams,soilInnerStateParams,soilOuterStateParams] =...
                                            run_model_leek(climateConsParams,...
                                                           climateState,...
                                                           cropStateParams,...
                                                           leekConsParams,...
                                                           fileSettings,...
                                                           managementSettings,...
                                                           simulationSettings,...
                                                           soilCommonStateParams,...
                                                           soilInnerStateParams,...
                                                           soilOuterStateParams);
            
            if any(managementSettings.hDateLeek == simulationSettings.t)
                date = datevec(simulationSettings.t);
                if date(1) == 1971
                    disp('')
                end
                
                cropStateParams = dataLog(cropStateParams, "Leek", fileSettings, "Force");
            else
                cropStateParams = dataLog(cropStateParams, "Leek", fileSettings);
            end
        end
        
        % apply irrigation to boundary conditions
        i=max(find(soilOuterStateParams.BOUNDARY_CONDITIONS_MATRIX(:,1)<=simulationSettings.t));
        soilOuterStateParams.BOUNDARY_CONDITIONS_MATRIX(i,3) = soilOuterStateParams.BOUNDARY_CONDITIONS_MATRIX(i,3) + ...
            soilOuterStateParams.epa + soilOuterStateParams.esa - managementSettings.irri; %no plant simulation
        soilInnerStateParams.BOUNDARY_CONDITIONS_MATRIX(i,3) = soilInnerStateParams.BOUNDARY_CONDITIONS_MATRIX(i,3) + ...
            soilInnerStateParams.epa + soilInnerStateParams.esa - managementSettings.irri; %plant simulation
        
    elseif ~simulationSettings.flag_double_sim
        
        % soil simulation when no plant is on the field
        cropConsParams.crop_type = NaN;
        cropConsParams.rorad = cauliConsParams.rorad;
        cropConsParams.rd0 = cauliConsParams.rd0;
        cropConsParams.g = cauliConsParams.g;
        cropConsParams.rdens0 = cauliConsParams.rdens0;
        cropStateParams.TotRootDepth = 0;
        index = climateState.ET0_time== simulationSettings.t;
        soilInnerStateParams.esa = 0;   %mm to cm
        soilInnerStateParams.epa = 0; 
        soilOuterStateParams.esa = climateState.ET0_cm_per_day(index);
        soilOuterStateParams.epa = 0;
        cropStateParams.ShootNDemand = 0;
        soilOuterStateParams.tot_upt = 0;
        soilInnerStateParams.tot_upt = 0;
        cropStateParams.RLengthLA = 0;
        cropStateParams.fraction_plant = 0;
        cropStateParams.fraction_soil = 1;
        cropStateParams.fraction_plantp = 0;
        cropStateParams.fraction_soilp = 1;
        cropStateParams.Rdens = NaN;
        cropStateParams.RLDRCil = zeros(85,1);
    end
    
    %% SOIL MODULE
    % set variables to stop simulation if no convergence can be reached
    soilInnerStateParams.noConvCounter = 0;
    soilOuterStateParams.noConvCounter = 0;
    soilOuterStateParams.STOP = 0;
    soilInnerStateParams.STOP = 0;
    
    % simulate soil processes in the unrooted zone
    soilOuterStateParams.JulianDay(end+1) = simulationSettings.t;
    soilCommonStateParams.JulianDay(end+1) = simulationSettings.t;
    
    while soilOuterStateParams.numboolean == 0 && simulationSettings.t<simulationSettings.tmax
        if round(simulationSettings.t)== simulationSettings.t
            
            
            
            %Store results before timestep for use in second cilinder
            %Water
            php = soilOuterStateParams.ph;
            wat_flxsap = soilOuterStateParams.wat_flxsa;
            pvelap = soilOuterStateParams.pvela;
            pvelahp = soilOuterStateParams.pvelah;
            pvelohp = soilOuterStateParams.pveloh;
            pvelop = soilOuterStateParams.pvelo;
            wciop = soilOuterStateParams.wcio; 
            wcmap = soilOuterStateParams.wcma;
            wcmop = soilOuterStateParams.wcmo;
            wcop = soilOuterStateParams.wco;
            diffusp = soilOuterStateParams.diffus; 
            
            %Temperature
            tempp = soilOuterStateParams.temp;
            temptoppp = soilOuterStateParams.temptopp;
            
            %Solute
            csolp = soilOuterStateParams.csol; 
            cimp = soilOuterStateParams.cim; 
            cmp = soilOuterStateParams.cm;
            csp = soilOuterStateParams.cs;
            tcsolop = soilOuterStateParams.tcsolo;
            dtp = soilOuterStateParams.dt;
            soil_omp = soilOuterStateParams.soil_om;
        end
        
        % simulation of soil processes in unrooted cilinder
        simulationSettings.simPlantFlag = 0; % in the unrooted cilinder no plant processes should be simulated
        if ~parallelFlag
            waitbar(progress/100,WAITBAR,[string(['simulating day ',num2str(day)]);strcat(num2str(progress),'% of the simulation is completed');...
                "SIMULATING UNROOTED SOIL PROCESSES"],'Name',"ECOFERT simulation",'height',500)
            if getappdata(WAITBAR,'canceling')
                delete(WAITBAR)
                error(['The simulation has been stopped manually the results or',...
                    ' progress are stored in ',char(ResultsPath)])
            end
        end
        
        [climateState,managementSettings,simulationSettings,soilConsParams,...
            soilOuterStateParams] = runwave(climateState,cropStateParams,cropConsParams,...
            managementSettings,simulationSettings,soilConsParams,soilOuterStateParams);
        

        
        if round(simulationSettings.t) == simulationSettings.t
            
            soilOuterStateParams.numboolean = 1;          
            
            
            % store time series of the soil parameters in the unrooted zone
            soilOuterStateParams = saveWAVEresults(climateState,...
                                                   fileSettings,...
                                                   soilOuterStateParams,...
                                                   simulationSettings,...
                                                   'Soil_unrooted');
        end
        if soilOuterStateParams.STOP
            break
        end
        
    end
    
    %reset settings for soil simulation of the next day
    soilOuterStateParams.numboolean = 0;
    soilOuterStateParams.initsol = 1;
    simulationSettings.simPlantFlag = 1;
    
    % check stop statement
    if soilOuterStateParams.STOP
        break
    end
    
    % simulate soil processes in the rooted soil zone
    if simulationSettings.flag_double_sim == 1 &...
            (simulationSettings.t-1) <= harvest_date(simulationSettings.ncrop)
        
        simulationSettings.t = simulationSettings.t-1; % reset simulation time to the previous day to simulate inner and outer soil processes on the same day
        soilInnerStateParams.JulianDay(end+1) = simulationSettings.t;
        
        if simulationSettings.t == plant_date(simulationSettings.ncrop)      
            
            soilInnerStateParams.ph = php;
            soilInnerStateParams.temp = tempp;
            soilInnerStateParams.temptopp = temptoppp;
            soilInnerStateParams.cim = cimp; 
            soilInnerStateParams.cm = cmp;
            soilInnerStateParams.csol = csolp;
            soilInnerStateParams.cs = csp; 
            soilInnerStateParams.tcsolo_ini = tcsolop;
            soilInnerStateParams.tcsolo = soilInnerStateParams.tcsolo_ini; 
            soilInnerStateParams.soil_om = soil_omp;
                        
            soilInnerStateParams.wat_flxsa = wat_flxsap;
            soilInnerStateParams.pvela = wat_flxsap;
            soilInnerStateParams.pvelah = wat_flxsap;
            soilInnerStateParams.pveloh = pvelohp;
            soilInnerStateParams.pvelo = pvelop;
            soilInnerStateParams.pvela = pvelap;
            soilInnerStateParams.pvelah = pvelahp;
            soilInnerStateParams.wcio = wciop;
            soilInnerStateParams.wcma = wcmap;
            soilInnerStateParams.wcmo = wcmop;
            soilInnerStateParams.wco = wcop;
            soilInnerStateParams.diffus = diffusp;
            
            if simulationSettings.simNitroFlag
                soilInnerStateParams.tcorgs = sum(soilInnerStateParams.soil_om(:,1) +...
                                                  soilInnerStateParams.soil_om(:,3) +...
                                                  soilInnerStateParams.soil_om(:,5));
                soilInnerStateParams.tnorgs = sum(soilInnerStateParams.soil_om(:,2) +...
                                                  soilInnerStateParams.soil_om(:,4) +...
                                                  soilInnerStateParams.soil_om(:,6));
            
            end
            soilInnerStateParams = calc_stock(simulationSettings,soilConsParams,soilInnerStateParams,0);
            soilInnerStateParams.stock_max = sum(moist_ret(zeros(1,simulationSettings.nComp),...
                soilInnerStateParams,soilConsParams,0)*simulationSettings.dx);
            soilInnerStateParams.initsol = 2;
        end
        
       
        while  soilInnerStateParams.numboolean ==0
            if ~parallelFlag
                waitbar(progress/100,WAITBAR,[string(['simulating day ',num2str(day)]);strcat(num2str(progress),'% of the simulation is completed');...
                    "SIMULATING ROOTED SOIL PROCESSES"],'Name',"ECOFERT simulation",'height',500)
                if getappdata(WAITBAR,'canceling')
                    delete(WAITBAR)
                    error(['The simulation has been stopped manually the results or',...
                        ' progress are stored in ',char(ResultsPath)])
                end
            end
            
            [climateState,managementSettings,simulationSettings,soilConsParams,...
                soilInnerStateParams] = runwave(climateState,cropStateParams,...
                cropConsParams,managementSettings,simulationSettings,soilConsParams,soilInnerStateParams);
                             
            soilCommonStateParams.und_time(end+1) =  soilInnerStateParams.uptake_matrix(3);
            soilCommonStateParams.unc_time(end+1) =  soilInnerStateParams.uptake_matrix(2);
                        
            if round(simulationSettings.t) == simulationSettings.t
                
                soilInnerStateParams.numboolean = 2;
                soilInnerStateParams = saveWAVEresults(climateState,...
                                                       fileSettings,...
                                                       soilInnerStateParams,...
                                                       simulationSettings,...
                                                       'Soil_rooted');
            end
            
            soilInnerStateParams.initsol = 1;
            if soilInnerStateParams.STOP
                break
            end
        end
        
        soilInnerStateParams.numboolean = 0 ;
        
        if soilInnerStateParams.STOP
            break
        end
        
        if simulationSettings.t - 1 == harvest_date(simulationSettings.ncrop)
            % the crop is harvested at the start of harvest date so the
            fractSoil_new = 1;
            fractPlant_new = 0;
            fractSoil = cropStateParams.fraction_soil(end);
            fractPlant = cropStateParams.fraction_plant(end);
                                         
        else
            
            fractPlant = cropStateParams.fraction_plant(end-1);
            fractSoil = cropStateParams.fraction_soil(end-1);
            fractPlant_new = cropStateParams.fraction_plant(end);
            fractSoil_new = cropStateParams.fraction_soil(end);
            
        end
        
    elseif simulationSettings.flag_double_sim == 0
        
        fractSoil = 1; 
        fractPlant = 0;
        fractSoil_new = 1; 
        fractPlant_new = 0;      
    end
    
    % Output of soil model to crop model
    if ~parallelFlag
        waitbar(progress/100,WAITBAR,[string(['simulating day ',num2str(day)]);strcat(num2str(progress),'% of the simulation is completed');...
            "CALCULATING SOIL OUTPUT TOWARDS CROP MODEL"],'Name',"ECOFERT simulation",'height',500)
        
        if getappdata(WAITBAR,'canceling')
            delete(WAITBAR)
            error(['The simulation has been stopped manually the results or',...
                ' progress are stored in ',char(ResultsPath)])
        end
    end

    crop_type = cropConsParams.crop_type; % select the right plant density based on the loaded crop parameters
    
    if crop_type == 1
        PLM2 = managementSettings.PLM2Cauli;
        soilCommonStateParams.NSupply = soilInnerStateParams.tot_upt*fractPlant_new/PLM2*10^4;

    elseif crop_type == 2
        PLM2 = managementSettings.PLM2Leek;
        soilCommonStateParams.NSupply = soilInnerStateParams.tot_upt*fractPlant_new/PLM2*10^4;

    end
    
	%reset tot_upt for the next day
	soilInnerStateParams.tot_upt = 0;
    soilOuterStateParams.tot_upt = 0;
    
    %% Calculate common soil properties based on (un)rooted soil    
    [soilCommonStateParams, soilInnerStateParams] = mergeInnerOuter(...
                                                 fractSoil, fractPlant,...
                                                 simulationSettings,...
                                                 soilConsParams,...
                                                 soilCommonStateParams,...
                                                 soilInnerStateParams,...
                                                 soilOuterStateParams);
    
    %Save data of commonstate parameters to external storage
    fieldNames = {'und_time','unc_time', 'water_storage',...
                  'WCSoil_log', 'UreumSoil_log','NH4Soil_log',...
                  'NO3Soil_log', 'UreumBalance_log',...
                  'NH4Balance_log', 'NO3Balance_log', 'N_reaction_balance',...
                  'daily_mineral','JulianDay','TSoil', 'water_balance'};
    
    if simulationSettings.t == simulationSettings.tmax | (soilInnerStateParams.STOP | soilOuterStateParams.STOP)
        soilCommonStateParams = dataLog(soilCommonStateParams, "SoilCommonStateParams",...
            fileSettings, 'Force', fieldNames);
    else
       soilCommonStateParams = dataLog(soilCommonStateParams, "SoilCommonStateParams",...
            fileSettings, false, fieldNames);
    end
        
    %% Exchange information between cilinders
    
    %Water
    [wcmo1a, wcmo2a] = mix_variables(fractSoil, fractPlant, fractSoil_new, fractPlant_new,...
        soilOuterStateParams.wcmo, soilInnerStateParams.wcmo);
    
    [soilOuterStateParams.wco, soilInnerStateParams.wco] = mix_variables(fractSoil,fractPlant,fractSoil_new,...
        fractPlant_new,soilOuterStateParams.wco,soilInnerStateParams.wco);
    
    %Water storage in the profile
    soilInnerStateParams.extra_water(end+1) = sum(wcmo2a)-sum(soilInnerStateParams.wcmo); %Extra water added during this time step
    soilOuterStateParams.extra_water(end+1) = sum(wcmo1a)-sum(soilOuterStateParams.wcmo);%Extra water added during this time step
    
    soilInnerStateParams.wcmo = wcmo2a;
    soilOuterStateParams.wcmo = wcmo1a;
    
    [soilInnerStateParams, soilOuterStateParams] = change_state_var(fractSoil_new,...
                                                                    fractPlant_new,...
                                                                    soilConsParams,...
                                                                    soilInnerStateParams,...
                                                                    soilOuterStateParams);    
    
    for xx=1:7
        var1 = soilOuterStateParams.soil_om(:,xx);
        var2 = soilInnerStateParams.soil_om(:,xx);
        [var1a, var2a] = mix_variables(fractSoil,fractPlant, fractSoil_new, fractPlant_new, var1,var2);
        soilOuterStateParams.soil_om(1:simulationSettings.nComp,xx)= var1a;
        soilInnerStateParams.soil_om(1:simulationSettings.nComp,xx) = var2a;
    end
    
    %temperature
    [soilOuterStateParams.temp, soilInnerStateParams.temp] = mix_variables(fractSoil,fractPlant,...
        fractSoil_new, fractPlant_new, soilOuterStateParams.temp,soilInnerStateParams.temp);
    
    % solute
    [soilOuterStateParams.cim, soilInnerStateParams.cim] = mix_variables(fractSoil,...
        fractPlant, fractSoil_new, fractPlant_new,soilOuterStateParams.cim,soilInnerStateParams.cim);
    
    [soilOuterStateParams.acsolmo, soilInnerStateParams.acsolmo] = mix_variables(fractSoil,...
        fractPlant, fractSoil_new, fractPlant_new, soilOuterStateParams.acsolmo,soilInnerStateParams.acsolmo);
    
    [soilOuterStateParams.acsolio, soilInnerStateParams.acsolio] = mix_variables(fractSoil,...
        fractPlant, fractSoil_new, fractPlant_new, soilOuterStateParams.acsolio, soilInnerStateParams.acsolio);    
    
    [soilOuterStateParams.tcsoloa, soilInnerStateParams.tcsoloa] = mix_variables(fractSoil,...
        fractPlant, fractSoil_new, fractPlant_new, soilOuterStateParams.tcsolo, soilInnerStateParams.tcsolo);  
    
    % calculate added solute due to increasing plant size
    if isempty(soilInnerStateParams.tcsolo)
        soilInnerStateParams.extra_sol(:,end+1) = [0 0 0]';
    else
        soilInnerStateParams.extra_sol(:,end+1) = [sum(soilInnerStateParams.tcsoloa) - ...
            sum(soilInnerStateParams.tcsolo)]'; %Extra solute added during this time step
    end
     
    soilOuterStateParams.tcsolo = soilOuterStateParams.tcsoloa;
    soilInnerStateParams.tcsolo = soilInnerStateParams.tcsoloa;
    soilOuterStateParams = rmfield(soilOuterStateParams,'tcsoloa');
    soilInnerStateParams = rmfield(soilInnerStateParams,'tcsoloa');
    
    if simulationSettings.flag_double_sim == 1
        soilInnerStateParams.cm = (soilInnerStateParams.tcsolo - soilInnerStateParams.acsolmo.*...
                               repmat(soilConsParams.bulkDensity, 1, managementSettings.nsol).*...
                               repmat(simulationSettings.dx,simulationSettings.nComp, managementSettings.nsol))./...
                               repmat(soilInnerStateParams.wcmo', 1, managementSettings.nsol);
    end
                       
    soilOuterStateParams.cm = (soilOuterStateParams.tcsolo - soilOuterStateParams.acsolmo.*...
                               repmat(soilConsParams.bulkDensity, 1, managementSettings.nsol).*...
                               repmat(simulationSettings.dx,simulationSettings.nComp, managementSettings.nsol))./...
                               repmat(soilOuterStateParams.wcmo', 1, managementSettings.nsol);

    % ADD PART THAT RESET ORGANIC N on January 1st of each year - JVL - 2018/08/21
    JANTEST = datevec(simulationSettings.t);
    
    if sum(JANTEST(2:6)) == 2
        soilOuterStateParams.soil_om = soilConsParams.soil_om;
        soilInnerStateParams.soil_om = soilConsParams.soil_om;
        
        soilOuterStateParams.tcorgs = sum(soilOuterStateParams.soil_om(:,1) +...
                                          soilOuterStateParams.soil_om(:,3) +...
                                          soilOuterStateParams.soil_om(:,5));
        soilOuterStateParams.tnorgs = sum(soilOuterStateParams.soil_om(:,2) +...
                                          soilOuterStateParams.soil_om(:,4) +...
                                          soilOuterStateParams.soil_om(:,6));
        soilInnerStateParams.tcorgs = sum(soilInnerStateParams.soil_om(:,1) +...
                                          soilInnerStateParams.soil_om(:,3) +...
                                          soilInnerStateParams.soil_om(:,5));
        soilInnerStateParams.tnorgs = sum(soilInnerStateParams.soil_om(:,2) +...
                                          soilInnerStateParams.soil_om(:,4) +...
                                          soilInnerStateParams.soil_om(:,6));
        disp('Organic Matter was reset..');
    end
    
    %% MANAGEMENT MODULE
    if simulationSettings.flag_double_sim % only calculate management practices when plant are simulated
        
        %AUTOIRRIGATIE-----------------------------------------------------
        if ~parallelFlag
            waitbar(progress/100,WAITBAR,[string(['simulating day ',num2str(day)]);strcat(num2str(progress),'% of the simulation is completed');...
                "APPLYING MANAGEMENT PRACTICES"],'Name',"ECOFERT simulation",'height',500)
            
            if getappdata(WAITBAR,'canceling')
                
                delete(WAITBAR)
                error(['The simulation has been stopped manually the results or',...
                    ' progress are stored in ',char(ResultsPath)])                
            end
        end
        
        soilCommonStateParams = calc_FC(cropStateParams,soilCommonStateParams,soilConsParams); %optimal water
        soilCommonStateParams.WCFC = soilCommonStateParams.WCFC * 0.5;
        soilCommonStateParams.waterinrootzone =...
            sum(soilCommonStateParams.WCSoil_log(1:cropStateParams.TotRootDepth(end),end)); %available water
        if soilCommonStateParams.WCFC > soilCommonStateParams.waterinrootzone
            managementSettings.irri = soilCommonStateParams.WCFC -...
                soilCommonStateParams.waterinrootzone;
        else
            managementSettings.irri = 0;
        end
        %------------------------------------------------------------------
        
        % Calculation of fertilizer advice
        % still to be implemented
        
    end
    if (any(simulationSettings.t == sort(reshape(managementSettings.fertMomentCauli,...
            numel(managementSettings.fertMomentCauli),1))+simulationSettings.tstart-1)...
            | any(simulationSettings.t == sort(reshape(managementSettings.fertMomentLeek,...
            numel(managementSettings.fertMomentLeek),1))+simulationSettings.tstart-1))...
            & simulationSettings.calcAddFertFlag % check wheter fertilizer should be applied the next day AND if the user enabled calcAddfert 
        
        
        % calculate remaining amonium and nitrogen in the soil and apply the
        % required amount the next to fulfill the target values depending on
        % the method of calculating the target values
        [soilInnerStateParams, soilOuterStateParams] = calc_addFert(cropStateParams,...
                                                                  managementSettings,...
                                                                  simulationSettings,...
                                                                  soilCommonStateParams,...
                                                                  soilInnerStateParams,...
                                                                  soilOuterStateParams);
    end
    
    calculationtime_day(end+1) = toc;
    if length(calculationtime_day) == fileSettings.DataLog_size
               
        fileName = join([fileSettings.ResultsPath,'timing_datalog_size_',...
                         num2str(fileSettings.DataLog_size),'.txt'],'');
        data = calculationtime_day(1:end-2);
        try            
            dlmwrite(fileName, data','-append','newline','pc')
        catch
            mkdir(fileSettings.ResultsPath)
            dlmwrite(fileName, data', '-append','newline','pc')
        end
        
        calculationtime_day = calculationtime_day(end-1:end);
    end
    
end


fileName = join([fileSettings.ResultsPath,'timing_datalog_size_',...
    num2str(fileSettings.DataLog_size),'.txt'],'');
data = calculationtime_day;
try
    dlmwrite(fileName, data','-append','newline','pc')
catch
    mkdir(fileSettings.ResultsPath)
    dlmwrite(fileName, data', '-append','newline','pc')
end


if soilInnerStateParams.STOP | soilOuterStateParams.STOP
    warning('The code is ended because no convergence can be reached')
end

if ~parallelFlag
    waitbar(progress/100,WAITBAR,["ECOFERT simulation completed";...
        "SAVING DATA"],'Name',"ECOFERT simulation",'height',500)
    if getappdata(WAITBAR,'canceling')
        
        delete(WAITBAR)
        error(['The simulation has been stopped manually the results or',...
            ' progress are stored in ',char(ResultsPath)])
        
    end
    delete(WAITBAR)
end
% clear variables except state and constant parameters
clearvars -except cauliConsParams climatConsParams fileSettings leekConsParams...
    managementSettings simulationSettings soilCommonStateParams...
    soilConsParams soilInnerStateParams soilOuterStateParams cropOutput

