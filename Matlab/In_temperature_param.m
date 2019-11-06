function [part, alamo, geomq, geomos,geomom,bd_water,heatcapsolids_corr, heatcapwater_corr]= In_temperature_param(units,ncomp)

%Input parameters for the temperature module as by Wierenga (and Wave
%Fortran)

%mainly fixed values
%Parameters are a fixed value,based on literature.
%Components: 
%quartz| organicmatter| other solids| air| water
part = [0.54 0.045 0.015 0 0];          % volume fraction of each component
alamo = [20.4 0.6 0.7 0.0615 0.176];    % thermal conductivity of each component (mCal cm-1 sec-1)
geomq = [0.125 0.125 .75];              % g-factor quartz
geomom = [0.5 0.5 0];                   % g-factor organic matter 
geomos = [0.125 0.125 0.75];            % g-factor other solids
heatcapsolids = 840;                    % specific heat capacity solids (J/kg/C)
heatcapwater = 4186.8;                  % specific heat capacity water (J/kg/C)


%%End of inputsection, all input parameters are converted to the correct
%%unit and/or put in a matrix form

part = repmat(part,ncomp+1,1);
alamo = repmat(alamo, ncomp+1,1);%thermal conductivity of each component (mCal cm-1 sec-1)
%(Conversion to J/L/T is done later)
geomq = repmat(geomq,ncomp+1,1);
geomos = repmat(geomos,ncomp+1,1);
geomom = repmat(geomom,ncomp+1,1);

%Conversion of heat capacity to J/M/C
if strcmp(units{1,3},'kg')==1 
    corr_M = 1;
elseif strcmp(units{1,3},'g')==1
    corr_M = 1/1000;
elseif strcmp(units{1,3},'mg')==1 
    corr_M = 1/1000000;
end
heatcapsolids_corr = heatcapsolids*corr_M; %in J/M/C
heatcapwater_corr = heatcapwater*corr_M; %in J/M/C

%Conversion bulk density water to M L^-3
if strcmp(units{1,1},'cm')==1 
    corr_L = 1/100^3;
elseif strcmp(units{1,1},'m')==1
    corr_L = 1;
elseif strcmp(units{1,1},'mm')==1 
    corr_L = 1/1000^3;
end
bd_water = 1000*corr_L/corr_M;
%1000 is bulk density in kg/m^3

%Convert thermal conductivity of to (mCal cm-1 sec-1 °C-1) to J L-1 T-1 C-1
%Calculate the conversion factor needed for the length
if strcmp(units{1,1},'cm')==1 
    corr_L = 1;
elseif strcmp(units{1,1},'m')==1
    corr_L = 100;
elseif strcmp(units{1,1},'mm')==1 
    corr_L = 0.1;
end

%Calculate the conversion factor needed for the time
if strcmp(units{1,2},'day' )==1
    corr_T = 60*60*24;
elseif strcmp(units{1,2},'hour' )==1
    corr_T = 60*60;
elseif strcmp(units{1,2},'min' )==1 
    corr_T = 60;
elseif strcmp(units{1,2},'s' )==1
    corr_T = 1;    
end

%use the correction factors
alamo = alamo*corr_T*corr_L*4.1868/1000;
%4.1868 is conversion factor for Calorie to Joule
%1000 is to correct for mili

