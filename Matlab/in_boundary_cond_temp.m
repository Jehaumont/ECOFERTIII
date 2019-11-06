function  [temp_BC, current_atmospheric_BC,p] = in_boundary_cond_temp(temp_BC,t,tmax,temp_table)

%Inputfile for temperature boundary conditions 
%Also retrieves the current BC.

%temp_BC is an array containing 4 elements:
%column 1:  type of upper boundary condition
%       if column 1 = 1: fixed temperature at the top 
%       if column 1 = 2: atmospheric boundary condition    
%column 2 : value of top boundary condition
%       if column 1 = 1: value of fixed temperature at top (°C)
%       if column 1 = 2: non applicable
%column 3 : type of bottom boundary condition
%   if column 3 = 1 : fixed temperature at bottom (°C)
%   if column 3 = 2 : zero flux bc at bottom
%column 4: 
%   if column 3 = 1 : value of fixed temperature at bottom 
%   if column 3 = 2 : value not used in calculation.   

%When atmospheric BC are applied, user should define the time to complete one cycle and an input a matrix containing
%maximum and minimum temperature for each cycle. 

%%%%%%%%%%%%
%Joachim Vansteenkiste 1/10/09

%Read in climatic data if climatic boundary condition is applied
if temp_BC(1) ==2 
        %Enter the time to complete one periodic cycle (T)                            
        climate_temperature =temp_table;
        p = climate_temperature(2,1) - climate_temperature(1,1); %assuming that this interval is 
                %equal to the time to complete one cycle (typically one day)   
        %temp_atmosperic_BC = [Time MaximumTemp MinimumTemp]
        %=ones(size(temp_atmospheric_BC,1),1)*[30 10]; %(°C) (For example
        %purposes)
        index = max(find(t>=climate_temperature(:,1)));
        current_atmospheric_BC = [climate_temperature(index,:)];
else %fixed temperature at the top of the profile
    current_atmospheric_BC = 0;
    p=0;
end

