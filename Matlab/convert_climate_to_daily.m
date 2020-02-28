function climateState = convert_climate_to_daily(climateState, fileSettings, soilConsParams)
%% FUNCTION INPUT
%climatedata
crop_climate = climateState.crop_climate; % climate data [matlabTime,radiation,temperature,relative humidity,rain]
ET0_cm_per_day = climateState.ET0_cm_per_day; % potential evapotranspiration following Penman-Monteith equation

%soil constant parameters
irri = soilConsParams.irri;

if isfield(fileSettings, 'irri_file')
	irri = load(fileSettings.irri_file);
else
    irri = zeros(1, 3);
end


%% FUNCTION MAIN BODY 
%Summarize data files to daily files
%JD H Rad T RH
uniquetimes = unique(crop_climate(:,1));
for i=min(uniquetimes):max(uniquetimes)
    index = find(crop_climate(:,1) ==i);
    dailyrain(i-min(uniquetimes)+1,1) = sum(crop_climate(index,6));
    temp_max(i-min(uniquetimes)+1,1) = max(crop_climate(index,4));
    temp_min(i-min(uniquetimes)+1,1) = min(crop_climate(index,4));
end

%Export to usable tables
theta_table = (min(uniquetimes):max(uniquetimes))';
theta_table(:,2) = ET0_cm_per_day';
theta_table(:,3) = dailyrain;
theta_table(:,4) = 0;

if length(irri) == 1
    theta_table(:,4) = irri(3);
else
    theta_table(:, 4) = irri(:, 3);
end



temp_table =(min(uniquetimes):max(uniquetimes))';
temp_table(:,3) = temp_max;
temp_table(:,2) = temp_min;

%% FUNCTION OUTPUT

climateState.climateDaily = theta_table; %daily climate
climateState.MinMaxTemp = temp_table;   %table with daily min and max temperatures
