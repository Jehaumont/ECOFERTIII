function climateState = convert_climate_to_matlab(climateState)
%% Function input
crop_climate = climateState.climate;

%% Function main body
%Convert the julian day in the inputfile to matlabtime
matlabtime = datenum(crop_climate(:,1),1,crop_climate(:,2),0,0,0);
crop_climate = [matlabtime crop_climate(:,3:7)];
%% Function output
climateState.crop_climate = crop_climate;