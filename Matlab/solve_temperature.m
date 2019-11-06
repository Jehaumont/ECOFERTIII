function soilStateCilinderParams = solve_temperature(climateState,simulationSettings,soilConsParams,...
    soilStateCilinderParams)

%% DOCUMENTATION
%Solve the 1-D heat flow equation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IN:
%dt: time step
%dx: spacing between the nodes
%ncomp= number of compartments
%t: current time
%WC: water content at each node at the current time step
%wcp: water content at each node of the previous time step
%temp: temperature profile of the beginning of the time step
%soil_parameters(1,2) = wcs: saturated water content
%temptopp: temperature at soil surface of the previous time step
%tmax: maximum time of the simulation
%
%OUT:
%temp = temperature profile at the current time step
%temptopp = temperature at soil surface at current time step
%
%CALLS: in_boundary_cond_temp
%CALLED BY: wavemat105_sol.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This m-file is based on soiltemp.for. Some changes were made with respect
%to the Fortran version. In this case, the original code from fortran is
%also mentioned as a comment. 
%Joachim V. 18/12/09
%

%% FUNCTION INPUT

% climate data
temp_table = climateState.MinMaxTemp;

% simulation settings
dx = simulationSettings.dx;
ncomp = simulationSettings.nComp;
t = simulationSettings.t;
tmax = simulationSettings.tmax;
units = simulationSettings.units;

% soil constant parameters
soil_parameters = [soilConsParams.wcr,soilConsParams.wcs,soilConsParams.alfa,...
                  soilConsParams.N,soilConsParams.ks,soilConsParams.lambda,soilConsParams.alfa_r];
solute_param1 = [soilConsParams.bulkDensity,soilConsParams.lm,soilConsParams.ratio,...
                soilConsParams.alf,soilConsParams.f];
temp_BC = soilConsParams.temp_BC;

% soil cilinder state parameters
dt = soilStateCilinderParams.dt;
temp = soilStateCilinderParams.temp;
temptopp = soilStateCilinderParams.temptopp;
WC = soilStateCilinderParams.WC;
wcp = soilStateCilinderParams.wcp;

%% FUNCTION MAIN BODY
%Read Boundary Conditions
[temp_BC, current_atmospheric_BC,p]= in_boundary_cond_temp(temp_BC,t,tmax,temp_table);

%Set top and bottom temperature from data about boundary conditions
if temp_BC(1) ==1 %Fixed temperature at the top
    temptop = temp_BC(2);
elseif temp_BC(1) ==2 %Atmospheric Boundary Condition
    tmax = current_atmospheric_BC(3);
    tmin = current_atmospheric_BC(2);
    temptop = (tmax + tmin)/2 + (tmax -tmin)/2 * sin(2*pi*t/p - 7*pi/12);
end
if temp_BC(3) ==1
    temp(ncomp) = temp_BC(4); %Set last node to a fixed temperature
end

%Set current value to be the value of the precedent time step
tempp=temp;

%Calculte thermal conductivity
[allam, beta] = calc_thermal_conduc(units,ncomp,soil_parameters, solute_param1,WC,wcp);
allam(:,:)=900;
beta(:,:)=2;
%[allam,beta] = calc_thermal_conduc2(units,ncomp,soil_parameters, solute_param1,WC,wcp);
%[allam, beta] = calc_thermal_conduc3(units,ncomp,soil_parameters, solute_param1,WC,wcp);
%[allam, beta] = calc_thermal_conduc4(units,ncomp,soil_parameters, solute_param1,WC,wcp);

% %for benchmarking purposes:
% allam = ones(ncomp+1,1).*1;
% beta = ones(ncomp+1,1).*1;

%Old fortran code could only simulate fixed temperature BC, matlab code can
% handle zero flux bottom BC as well.

% if temp_BC(3) ==1 %Fixed temperature at the bottom ==> Use the old Fortran code (works better than matlab, bottom temperature not completely fixed)
% %Fortran code shows large oscillations in the first nodes under climatic
% %top BC
% %FORTRAN VERSION %Calculate elements of tridiagonal matrix
%    for i=1:ncomp-1
%        AA = (dt*allam(i))/(beta(i)*2*dx*dx);
%        BB = (dt*allam(i+1))/(beta(i+1)*2*dx*dx);
%        A1(i) = -AA;
%        B1(i) = AA+BB+1;
%        C1(i) =-BB;
%        if i==1
%            D1(i) = AA*temptop+(1-AA-BB)*tempp(1) + BB*tempp(2);
%            F1(1) = C1(i)/B1(i);
%            G1(1) = (D1(i)-A1(i)*temptop)/B1(i);
%        elseif i ==ncomp
%            D1(i) = AA*tempp(i-1)+(1-AA-BB)*tempp(i) + BB*tempp(i+1);
%            F1(i) = C1(i)/(B1(i) -F1(i-1)*A1(i));
%            G1(i) = (D1(i)-A1(i)*G1(i-1))/(B1(i) -A1(i)*F1(i-1));
%        else
%            D1(i) = AA*tempp(i-1)+(1-AA-BB)*tempp(i) + BB*tempp(i+1);
%            F1(i) = C1(i)/(B1(i) -F1(i-1)*A1(i));
%            G1(i) = (D1(i)-A1(i)*G1(i-1))/(B1(i) -A1(i)*F1(i-1));
%        end
%    end
% % %Solve via backward substitution (adapted from Fortran version)
%  for i=ncomp-1:-1:1
%      temp(i)=G1(i)-F1(i)*temp(i+1);
%  end


%elseif temp_BC(3)==2 %Zero flux bottom boundary condition ==> Use new matlab code
%Code here presented can also work with the fixed temperature BC, but does
%not work properly( last node is not fixed completely, small variation
%occur)
% %Calculate elements of tridiagonal matrix - Matlab version
i=[1:ncomp];
X = (dt*allam(i))./(beta(i)*2*dx*dx);
Y = (dt*allam(i+1))./(beta(i+1)*2*dx*dx);
B1(1)= X(1)+Y(1)+1;
C1(1) = -Y(1);
D1(1) = X(1)*temptopp+(1-X(1)-Y(1))*tempp(1) + Y(1)*tempp(2) + X(1)*temptop;
if temp_BC(3) ==1 %fixed value at the bottom of the profile
    A1(ncomp-1) = -X(ncomp);
    B1(ncomp) = X(ncomp)+Y(ncomp)+1;
    D1(ncomp) = X(ncomp)*tempp(ncomp-1)+(1-X(ncomp)-Y(ncomp))*tempp(ncomp) + 2*Y(ncomp)*temp_BC(4);
elseif temp_BC(3) ==2 %zero flux condition
    A1(ncomp-1) = -X(ncomp);
    B1(ncomp) = X(ncomp)+Y(ncomp)+1-Y(ncomp);
    D1(ncomp) = X(ncomp)*tempp(ncomp-1)+(1-X(ncomp)-Y(ncomp))*tempp(ncomp) + Y(ncomp)*tempp(ncomp);
end
i=[2:ncomp-1];
A1(i-1) = -X(i);
B1(i) = X(i)+Y(i)+1;
C1(i) =-Y(i);
D1(i) = X(i)'.*tempp(i-1)+(1-X(i)'-Y(i)').*tempp(i) + Y(i)'.*tempp(i+1);

%Solve the system
    temp =(inv(diag(B1)+diag(C1,1)+diag(A1,-1))*D1').';
%end

%Save calculated top temperature
temptopp = temptop;

%% FUNCTION OUTPUT
soilStateCilinderParams.temp = temp;
soilStateCilinderParams.temptopp = temptopp;