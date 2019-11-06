function  [rtex,EPRA] = RER(t,ph,simplant,dx,units,plant_date,harvest_date,drz,epa,fraction_plant)

% Calculation of the extraction profile
%
%IN:
%   t: current time
%   ph: ph distribution
%   simplant: plant simulation (logical)
%   dx: spatial discretization
%OUT:
%   rtex: root extraction rate (1/T)
%   EPRA: total root extraction (1/T)
%CALLS:
%   In_ETsplit.m
%CALLED BY:
%   state_var.m
%--------------------
% F. Hupet 15/10/2000
% Javaux, 2006
%Modified by M. Sall 21/11/2009

%Read input from the inputfile in_root_uptake
[Smax_given,Smax_param] = in_root_uptake();

ncs=length(ph);%nber of comp.

%Extract from these parameters the parameters to determine Smax
if Smax_given
    Smax=Smax_param(9:end);
else
    arer= Smax_param(9);
    brer= Smax_param(10);
end

boolean_plant = (t>= plant_date);
boolean_harvest = (t<= harvest_date);
boolean_plant_harvest = boolean_plant + boolean_harvest;
index = find(boolean_plant_harvest == 2);

if simplant && isempty(index)==0
    irz=ceil(drz/dx); %Number of compartments in the root zone
    
    %Calculate Smax if not given
    %calculate maximal root extraction rate (Smax) if Smax is not input
    if Smax_given==0
        Smax(1:ncs,1)=0;
        i=1:irz-1;
        heigth=dx*i;
        Smax(i,1)=arer-(brer.*heigth);           %%%% (T-1)
        if irz>drz/dx
            Smax(irz,1)=0;
        end
    end
    
    %Upscale Smax
    % if fraction_plant == 0
    %     Smax = Smax/0.005;
    % else
    Smax = Smax/fraction_plant;
    % end
    %IN WAVE.FOR 
    %Calculate the sum of maximal root extraction rate
    sum_smax=sum(Smax);
    %Calculate fraction of root distribution for each compartment
    fraction=Smax/sum_smax;
    % recalculate smax to avoid the depassment of epa
    if sum_smax < epa
        %    disp('potential evaporation exceeds maximum root extraction')
    else
        Smax=fraction*epa/dx;     %%%ESSAI A ENLEVER
    end
    
    
    %if irz>drz/dx;
    %    Smax(irz,1)=fraction*epa/(irz*dx-drz);
    %end
    
    %Set negative values of Smax to zero
    dummy = Smax<0;
    Smax = (1-dummy).*Smax + (dummy).*0;
    
    %Reduction of potential root water uptake
    alpha_h(1:ncs,1)=0;
    rtex(1:ncs,1)=0;
    
    p0 = Smax_param(1);
    p1 = Smax_param(2);
    p2l = Smax_param(3);
    p2h = Smax_param(4);
    p3 = Smax_param(5);
    epa_high = Smax_param(6);
    epa_low = Smax_param(7);
    ilinr = Smax_param(8);
    
    p2 = p2h;
    if epa < epa_low
        p2 = p2l;
    elseif epa < epa_high
        p2 = p2h +((epa_high - epa)/(epa_high - epa_low))*(p2l - p2h); %Based on fortran
    end
    
    
    %Determination of alpha_h for each compartiment
    dummy1 = ph >=p0;
    dummy2 = ph>=p1 & ph < p0;
    dummy3 =  ph >= p2 & ph < p1;
    dummy4 = ph >= p3 & ph < p2 ;
    dummy5 = ph < p3;
    
    if ilinr
        alpha_h =...
            dummy1.*0 +...
            dummy2.*(p0-ph)./(p0-p1) +...
            dummy3*1+...
            dummy4.*ilinr.*((p3-ph)./(p3-p2))+...
            dummy5.*0;
    else
        dummy1.*0 +...
            dummy2.*(p0-ph)./(p0-p1) +...
            dummy3*1+...
            dummy4.*((1-ilinr)*(10.^((p2-ph)./p3)))+...
            dummy5.*((1-ilinr)*(10.^((p2-ph)./p3)))
        %Theory: dummy5.*0; but no smooth curve as in figure 2-9
    end
    
    alpha_h = alpha_h';
    %Set all reduction factors outside the rootzone tot 0 (no uptake)
    alpha_h (irz+1:end) = 0;
    
    %Determination of root extractio for each compartiment
    
    rtex=alpha_h.*Smax;       %(T-1)
    EPRA=cumsum(rtex*dx);     % L T-1
    
else
    rtex=zeros(ncs,1);
    EPRA=0;
end

if t>=684
    disp('')
end