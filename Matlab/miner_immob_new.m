function [minerm,mineri,rcarbmin,soil_om,om_appl,carbman, rnitman, carblit, rnitlit]=...
     miner_immob_new(t,dt,temp,ncs,WC,csol,ph,soil_om,wco,dx,cm,cim,soil_parameters,om_param, miner_param, plant_uptake_param,om_appl,initsol);

      
%Calculate the nitrogen mineralisation or immobilisation

%IN:
%t= time(min)
%dt = the time increment (min)
%temp=temperature (°C)
%WC:water content at time t, wcp: water content at time t-dt
%cmp:solute mobile concentration at time t-dt
%soil parameters for water flow
%nitro_param =nitrogen parameters (from script in_nitro_om_param)
%om_param = organic matter parameters (from script in_nitro_om_param)
%c_manp,n_manp,c_litp,n_litp,c_hump,n_hump: carbone and nitrogen in manure,
%litter and humus (from in initial om)
%OUT:
%minerm= amount of mineralisation or immobilisation in the mobile zone
%minerm= amount of mineralisation or immobilisation in the immobile zone
%rcarbmin= organic carbon

%CALL: OM_APPLIC
%CALLED BY: SOLVE_SOLUTE

%-------------------------------------------------------------------------
% Mamadou SALL 17/03/09
csolo = csol;
%% SOIL PARAMETER
 wcr=soil_parameters(:,1);
 wcs=soil_parameters(:,2);
 alfa=soil_parameters(:,3);
 n=soil_parameters(:,4);
 m=1-1./n;

 %%%%% NITROGEN AND ORGANIC MATER PARAMETER
 ro=miner_param(1);  
 fe=miner_param(2);  
 fh=miner_param(3);  
 k_lito=om_param(:,5); 
 k_mano=om_param(:,6);  
 k_humo=om_param(:,7); 
 
%%%%%%%%%%%%%%%%INITIAL OM
 c_manp=soil_om(:,1);
 n_manp=soil_om(:,2);
 c_litp=soil_om(:,3);
 n_litp=soil_om(:,4);
 c_hump=soil_om(:,5);
 n_hump=soil_om(:,6);
 cco2p = soil_om(:,7);
 
 
volume=dx;
% Find the addition of organic carbon and organic nitrogen due to organic fertilisation
%--------------------------------------------------------------------------
[om_applic,om_appl]=om_app(t,dt,om_appl);
%--------------------------------------------------------------------------
% Consider that organic matter is filled in the first compartiment 
%organic matter application(om_applic) is added in the manure pool and litter pool (see Nitboco.FOR) 
carbman = 0.3*om_applic(1,1); % organic carbon added in the manure pool of the first compariment
rnitman = carbman*(1/ro);   % organic nitrogen added in the manure pool of the first compariment
carblit = om_applic(1,1)-carbman; % organic carbon added in the litter pool of the first compariment
rnitlit = om_applic(1,2)-rnitman; % organic nitrogen added in the litter pool of the first compariment

% new values of organic carbone and organic nitrogen in the first compartiment
c_manp(1)=c_manp(1)+carbman;
c_litp(1)= c_litp(1)+carblit; 
n_manp(1)= n_manp(1)+rnitman;
n_litp(1)= n_litp(1)+rnitlit;

cnlito=n_litp;
cnmano=n_manp;
cchumo=c_hump;
cnhumo=n_hump;
ccmano=c_manp;
cclito=c_litp;
cco2o = cco2p;


% CALCULATE reduction factor FOR ACTUAL WATER CONTENT and temperature
% temp is the soil temperature (°C)
se=(wco'-wcr)./(wcs-wcr);
wc1=(wcr+(wcs-wcr)./(1+(alfa*100).^n).^m)';       %%WC for ph=100
wc2=(wcr+(wcs-wcr)./(1+(alfa.*1000).^n).^m)';    %%WC for ph=1000
wc3= (wcr+(wcs-wcr)./(1+(alfa.*15000).^n).^m)';  %%% WC f
%% reduction factor for temperature
redtemp=3.^((temp-16)/10);
for i=1:ncs
 %%   reduction factor for moisture
   if wco(i)>wcs(i)                     %%%0-wcs
       redmoist(i)=0.5;
   elseif wco(i)<wcs(i) & wco(i)>wc1(i);         %%%wcs-100
        redmoist(i)=0.5+0.5*(wcs(i)-wco(i))/(wcs(i)-wc1(i));
   elseif wc3(i)<wco(i) & wco(i)<wc2(i);        %%%% 1000-15000
        redmoist(i)=(WC(i)-wc3(i))/(wc2(i)-wc3(i)); 
   elseif wco(i)<wc3(i)                       %%%%%%%%%%%%>15000
        redmoist(i)=0;
   else                                     %%%%%%%%%%%%%%%100-1000
       redmoist(i)=1;                   
   end
end
%%%reduction factor
red_fact=redtemp'.*redmoist';
red_fact =red_fact./red_fact;
%%%effective coefficient for the actual time
rmano= n_manp/c_manp-fe/ro;     %%pxm in Wave.FOR
rlito=n_litp/c_litp -fe/ro;    %%pxl in Wave.FOR

rklit=k_lito;
rkman=k_mano;
rkhum=k_humo;
redmoi=redmoist;

 %CALCULATE AMOUNT OF MINERALISATION

for i=1:ncs
    red_fact = redmoi(i) * redtemp(i);
    pxl=(cnlito(i)/cclito(i)-(fe/ro));
    pxm=(cnmano(i)/ccmano(i)-(fe/ro));
    eff_rate_lit = rklit(i)*red_fact;
    eff_rate_man = rkman(i)*red_fact;
    eff_rate_hum = rkhum(i)*red_fact;
    check_mim = ( eff_rate_hum*cnhumo(i)+ pxm*eff_rate_man*ccmano(i)+...
    pxl*eff_rate_lit*cclito(i))*dt;
    if check_mim >= 0 %Mineralisation
        change_nh4 = check_mim;
        change_no3 = 0.;
        eff_rate_lit = rklit(i)*red_fact;
        eff_rate_man = rkman(i)*red_fact;
    else  %Immobilisation
        if((csolo(i,2)-1.8)*wco(i)*dx >= (-check_mim)) %Enough ammonium to immobilise
            change_nh4 = check_mim;
            change_no3 = 0;
            eff_rate_lit = rklit(i)*red_fact;
            eff_rate_man = rkman(i)*red_fact;
        else %Not enough ammonium, use as much as possible and try with nitrate
            change_nh4 = -(csolo(i,2)-1.8)*wco(i)*dx;
            chno3 = (csolo(i,2)-1.8)*wco(i)*dx+check_mim;
            if((csolo(i,3)-6.2)*wco(i)*dx >= (-chno3)) %Enough nitrate
                change_no3 = chno3;
                eff_rate_lit = rklit(i)*red_fact;
                eff_rate_man = rkman(i)*red_fact;
            else %Not enough nitrate, use as much as possible
                change_no3 = -(csolo(i,3)-6.2d0)*wco(i)*dx;
                shortage_min_n = (csolo(i,3)-6.2d0)*wco(i)*dx+chno3;
                if pxm < 0 & pxl >= 0 
                    vimmb = pxm*eff_rate_man*ccmano(i)*dt;
                    vimma = vimmb-shortage_min_n;
                    eff_rate_lit = rklit(i)*red_fact;
                    eff_rate_man = vimma/(pxm*ccmano(i)*dt)
                elseif pxl < 0 & pxm >= 0
                    vimmb = pxl*eff_rate_lit*cclito(i)*dt;
                    vimma = vimmb-shortage_min_n
                    eff_rate_lit = vimma/(pxl*cclito(i)*dt);
                    eff_rate_man = rkman(i)*red_fact;
                elseif pxl< 0 & pxm < 0 
                    vimmb = pxl*eff_rate_lit*cclito(i)*dt;
                    vimma = vimmb-shortage_min_n;
                    if vimma > 0 
                        vimmb = pxm*eff_rate_man*ccmano(i)*dt;
                        vimma = vimmb+vimma;
                        eff_rate_lit = 0;
                        eff_rate_man = vimma/(pxm*ccmano(i)*dt);
                    else
                        error ('troubles in rminimm')
                        %call stop_simulation ('programme stopped : check err_file')
                    end
                else
                    error('troubles in rminimm')
                    % stop_simulation ('programme stopped : check err_file')
                end
            end
        end
    end
 
% the gain/losses of the nitrogen litter pool
    % through mineralisation/immobilisation and nitrogen humification
    cnlito(i)=cnlito(i)+(-pxl*eff_rate_lit-fe*fh/ro*eff_rate_lit)*cclito(i)*dt;
    %c the gain/losses of the nitrogen litter pool
    %c through mineralisation/immobilisation and nitrogen humification
    cnmano(i)=cnmano(i)+(-pxm*eff_rate_man-fe*fh/ro*eff_rate_man)*ccmano(i)*dt;
    % the production of co2
    cco2o(i)=cco2o(i)+(1-fe)*(eff_rate_man*ccmano(i)+eff_rate_lit*cclito(i))*dt+eff_rate_hum*cchumo(i)*dt;
    
   % the production of carbon humus
    cchumo(i)=cchumo(i)+fe*fh*(eff_rate_man*ccmano(i)+eff_rate_lit* cclito(i))*dt-eff_rate_hum*cchumo(i)*dt;
    % the changes in the n-humus pool through mineralisation and n-litter humification
    cnhumo(i)=cnhumo(i)+(fe*fh/ro*(eff_rate_lit*cclito(i)+eff_rate_man*ccmano(i))-eff_rate_hum*cnhumo(i))*dt;
    % the change in the manure pool
    ccmano(i)=ccmano(i)+(-fe*fh-(1-fe))*eff_rate_man*ccmano(i)*dt;
    % the change in the litter pool
    cclito(i)=cclito(i)+(-fe*fh-(1-fe))*eff_rate_lit*cclito(i)*dt;
    % retain the rate of change of the ammonia and nitrate species
    % to be incorporated in the sequential transformation process (mg/(day*liter))
    minerm(i,2)=change_nh4/(dt*volume);%%%rminm(i,2)=change_nh4/(dt*volume)
    minerm(i,3)=change_no3/(dt*volume);
    mineri(i,2)=0.0;
    mineri(i,3)=0.0;
    rcarbmin(i)=(1-fe)*(eff_rate_man*ccmano(i)+...
        eff_rate_lit*cclito(i))/volume + eff_rate_hum*cchumo(i)/volume;
end

n_litp=cnlito;
n_manp=cnmano;
c_hump=cchumo;
n_hump=cnhumo;
c_manp=ccmano;
c_litp=cclito;

soil_om=[c_manp n_manp c_litp n_litp c_hump n_hump cco2o];   



