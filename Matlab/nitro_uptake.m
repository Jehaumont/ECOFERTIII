function[uptakem,uptakei,tot_upt, first_time, uptake_matrix,rdens,und,unc]=nitro_uptake(csol,dt,dx,harvest_date,...
        ncs,nsol,soil_parameters,plant_uptake_param,solute_param1,solute_param2,...
        plant_date,simplant,t,rtex,tot_upt,wco,first_time,drz,NDemand,PLM2,PropFact,RLengthLa,rdens,ncrop,fraction_plant,...
		rcsolo)


% Calculate the plant N uptake
%IN:
%t,dt,tmax= time, time increment and maximum time(min)
%nsol, ncs=number of solute spieces and number of compartiments
%WC:water content at time t, wcp: water content at time t-dt
%cmp:solute mobile concentration 
%wco= water content
%soil_parameters= parameters of the soil
%solute_param1,solute_param2 = parameters of the solute (from scrip%in_soliue_parameters)
%nitro_param = nitrogen parameters (from scrip in_nitro_om_param)
%rtex= root water extraction calculated in the water flow module

%:
%tot_upt= cumulative N uptake
%uptakem and uptakei = N uptake in the mobile and immobile region
%CALLED BY solve_solute
%--------------------------------------------------------------------------
%M.SALL 17/03/09
ratio=solute_param1(:,3);
dif=solute_param2(:,2*nsol+1:3*nsol);
ar=solute_param2(:,3*nsol+1:4*nsol);
br=solute_param2(:,4*nsol+1:5*nsol);

%%%%%%%%%%%%%%%%%%%%%%%%%%        

volume=dx;

%Find if there is a crop and which crop.
boolean_plant = (t>= plant_date);
boolean_harvest = (t<= harvest_date);
boolean_plant_harvest = boolean_plant + boolean_harvest;
index = 1;

%Read nitrogen uptake parameters
rorad=plant_uptake_param(index,1);       %%rorad=nitro_param(4);
rdo=plant_uptake_param(index,2);         %%rdo=nitro_param(5);
g = plant_uptake_param(index,3);           %%%g= nitro_param(13);  
%rnmaxp=plant_uptake_param(index,4);  %%rnmaxp=nitro_param(1);
rdens0=plant_uptake_param(index,5);   %%rdens0=nitro_param(2);
alfa_rdens=plant_uptake_param(index,6);  %%alfa_rdens=nitro_param(3);

rnmaxp = NDemand/(fraction_plant)*PLM2/10^4;

if simplant ==0 || isempty(index)==1 %No plant simulation or no crop on field
     uptakem(1:ncs,1:nsol)=0;
     uptakei(1:ncs,1:nsol)=0;
     uptake_matrix = [t, 0, 0, 0, 0,0];
     tot_upt = 0;
     rdens = zeros(ncs,1);
     drza = 0;
     unc=zeros(ncs,3);
     und=zeros(ncs,3);
     tconn = 0;
     tund = 0;
     totvr = 0;
     totdem = 0;
else %There is a cropmodeling and the crop is on the field
    if first_time == 1
        first_time =0;
    end


%Calculate the potential uptake during the (variable) timestep
if(rnmaxp-tot_upt)<= 0 % Total uptake of that day already exceeds demand
   uptakem(1:ncs,1:nsol)=0;
   uptakei(1:ncs,1:nsol)=0;   
   totdem = 0; %No extra uptake required
   totvr = 0;
else
    totdem = (rnmaxp-tot_upt);
     %Calculate the maximal total convective nitrogen uptake (tconn) (M L-2)
    tconn=0;
    for sp=2:3
        for i=1:ncs
            hcsolo(sp) = max(0,(csol(i,sp)-rcsolo(sp)));
            unc(i,sp) = (rtex(i)*hcsolo(sp).*wco(i)*dt*dx);
            tconn=tconn+(unc(i,sp));
        end
    end
    
    %Check it toral convective is enough to fulfill demand
    if tconn > totdem %=>Demand can ce fulfilled with convective uptake
        unc(1:ncs,1:nsol)=unc(1:ncs,1:nsol)*totdem/tconn;
        und(1:ncs,1:nsol) =0;
        tconn=totdem;
        tund=0;
    else %==>Active transport is needed to complete demand
        tpdnup = totdem - tconn;
        tund=0;
        %Calculate the total active transport possible
        for sp=2:3
            for i=1:ncs
                hcsolo(sp) = max(0,(csol(i,sp)-rcsolo(sp)));
                diffus_rm(i,sp) = dif(i,sp).*ar(i,sp).*exp(br(i,sp).*wco(i))/wco(i);
        	    und(i,sp)=rdens(i)*rorad*2*pi*diffus_rm(i,sp)*hcsolo(sp)*wco(i)*dx*dt/(rdo*1e05);
        	    tund=tund+und(i,sp);
            end
        end
        if tund > tpdnup  %==>Active transport can complete the demand
            for sp = 2:3
                und(1:ncs,sp)=und(1:ncs,sp)*tpdnup/tund;
            end
            tund = tpdnup;

        end
    end
    
    %calculate the uptake rates
    for sp=2:3
        for i=1:ncs
            uptakem(i,sp)= (unc(i,sp)+und(i,sp)).*ratio(i)/(dx*dt);
            uptakei(i,sp)= (unc(i,sp)+und(i,sp)).*(1-ratio(i))/(dx*dt);
        end
    end
    
end


for sp = 2:3
    for i=1:ncs-1
        sol_sinkm(i,sp)= - uptakem(i,sp);
        supply = csol(i,sp)*(wco(i))/dt;
        if sol_sinkm(i,sp)<=0 %solute has disappeared from the profile by uptake
            demand=-sol_sinkm(i,sp); %becomes a positive number
            if demand >= supply
                deficit=supply-demand; %deficit is a negative number
                sol_sinkm(i,sp)= min(sol_sinkm(i,sp)- deficit, 0);
                corr_uptakem(i,sp)= max(uptakem(i,sp)+ deficit, 0);
            else
                corr_uptakem(i,sp) = uptakem(i,sp);
            end
            if uptakem(i,sp)~=0
                disp('')
                unc(i,sp) = unc(i,sp)/uptakem(i,sp)*corr_uptakem(i,sp);
                und(i,sp) = und(i,sp)/uptakem(i,sp)*corr_uptakem(i,sp);
            else
                unc(i,sp) = 0;
                und(i,sp) = 0;
            end
            
        end
    end
end

tund = sum(sum(und));
tconn = sum(sum(unc));
totvr = tund + tconn;
tot_upt = tot_upt + tund + tconn; %Total uptake is reset every day!

uptakem = corr_uptakem;
uptakem(ncs,2:3) = 0;

end
unc(ncs,1:3) = 0;
und(ncs,1:3) = 0;
uptake_matrix = [t, tconn, tund, totvr, tot_upt, totdem];