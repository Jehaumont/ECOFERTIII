function [Smax_given,Smax_param] = in_root_uptake()

% Uptake of water by plants is described via Smax,the maximum uptake rate. 
% In this model there are 2 options: 
% (0) linear function to describe the uptake: root uptake = arer - brer * x where x is the depth of the node
% (1) manual input of Smax for each compartment
Smax_given =0; %Number as in options described above

if Smax_given ==0      %define the parameters arer and brer
    arer = 0.04;          % T-1
    brer = 0.004/6;         % T-1 L-1
else                  %else Smax_given == 1 then manual input the matrix (one value for each compartment)
    error('input matrix in the file in_root_uptake')
end

%Actual uptake of water is calculated by mulitplying a factor alpha
%(function of pressure head) with the maximum uptake rate.
%Is there a linear relationschip between reduction factor of the root water uptake
%(alpha) and the pressure head? (yes = 1, no = 0=hyperbolic); 
ilinr = 1;
%Specify critical pressure heads
p0 =-15;        %below this pressure head start extraction (above = anaerobis)
p1 = -30 ;      %below this pressure head optimal extraction
p2h = -320;     %below this pressure head no longer optimal at high evaporative demand (=epa_high)
p2l = -600 ;    %below this pressure head no longer optimal at low evaporative demand (=epa_low)
p3 =-8000 ;     %below this pressure head no extraction (= wilting point)
epa_high = 0.5 ; %(L T-1)
epa_low = 0.1 ; %(L T-1)

%% Store Smax parameters in one variable
if Smax_given
    Smax_param=[p0 p1 p2l p2h p3 epa_high epa_low ilinr Smax];
else
    Smax_param = [p0 p1 p2l p2h p3 epa_high epa_low ilinr arer brer ];
end

