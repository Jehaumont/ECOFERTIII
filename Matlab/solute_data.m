function [cm,cim,csol,acsolmo,acsolio,initsol]=solute_data (cm,cim,csol,acsolmo,acsolio,t,dx,tcsolo_ini,wco,wcmo,wcma,wcio,immobile,solute_param1,solute_param2,nsol,ncs,pvela,initsol);

% Calculation of initial solute data (only first timestep)

%IN: tcsolo_ini: initial mass of solute in the profile (ML-²)
%wco: total water content
%wcmo: water content for the mobile region
%wcma: water content between nodes for the mobile region
%wcio: water content for the immobile region
%pvela:pvela: pore velocity between the nodes
%solute_param1,solute_param2: parameters of the solutes
%OUT:
%cm: concentration in the profile for the mobile region
% cim: concentration in the profile for the immobile region
%CALL: none
%CALLED BY: solve_solute
%------------------------------------------------------------------
%M. Sall 10/04/09

%Read solute parameters
bd=solute_param1(:,1);
f=solute_param1(:,5);
kd=solute_param2(:,nsol+1:2*nsol);
ar=solute_param2(:,3*nsol+1:4*nsol);
br=solute_param2(:,4*nsol+1:5*nsol);
dif=solute_param2(:,2*nsol+1:3*nsol);
lm=solute_param1(:,2);
lm(ncs+1)=lm(ncs);
volume=dx;

if initsol==0
    for j=1:nsol
        for i=1:ncs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            if immobile
                tcsolio(i,j)=tcsolo_ini(i,j)*wcio(i)/wco(i);
                tcsolmo(i,j)=tcsolo_ini(i,j)*wcmo(i)/wco(i);
                cim(i,j)=tcsolio(i,j)/(volume*(kd((i),j)*(1-f((i)))*bd((i))+wcio(i)));
                cm(i,j)=tcsolmo(i,j)/(volume*(kd((i),j)*f((i))*bd((i))+wcmo(i)));
                csol(i,j)=(cim(i,j)*wcio(i)+cm(i,j)*wcmo(i))/wco(i);
                acsolio(i,j)=kd((i),j)*(1-f((i)))*cim(i,j);
                acsolmo(i,j)=kd((i),j)*	f((i))*cm(i,j);
                acsolo(i,j)=acsolio(i,j)+acsolmo(i,j);
            else
                tcsolio(i,j)=tcsolo_ini(i,j)*wcio(i)/wco(i);
                tcsolmo(i,j)=tcsolo_ini(i,j)*wcmo(i)/wco(i);
                csol(i,j)=tcsolo_ini(i,j)/(volume*(kd((i),j)*bd((i)) +wco(i)));
                cim(i,j)=0.0;
                cm(i,j)=csol(i,j);
                acsolio(i,j)=0.0;
                acsolo(i,j)=kd((i),j)*csol(i,j);
                acsolmo(i,j)=acsolo(i,j);
            end
            %%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end

if initsol==2
    for j=1:nsol
        for i=1:ncs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            if immobile
                tcsolio(i,j)=tcsolo_ini(i,j)*wcio(i)/wco(i);
                tcsolmo(i,j)=tcsolo_ini(i,j)*wcmo(i)/wco(i);
                cim(i,j)=tcsolio(i,j)/(volume*(kd((i),j)*(1-f((i)))*bd((i))+wcio(i)));
                cm(i,j)=tcsolmo(i,j)/(volume*(kd((i),j)*f((i))*bd((i))+wcmo(i)));
                csol(i,j)=(cim(i,j)*wcio(i)+cm(i,j)*wcmo(i))/wco(i);
                acsolio(i,j)=kd((i),j)*(1-f((i)))*cim(i,j);
                acsolmo(i,j)=kd((i),j)*	f((i))*cm(i,j);
                acsolo(i,j)=acsolio(i,j)+acsolmo(i,j);
            else
                tcsolio(i,j)=tcsolo_ini(i,j);
                tcsolmo(i,j)=tcsolo_ini(i,j);
                csol(i,j)=csol(i,j);
                cim(i,j)=0.0;
                cm(i,j)=csol(i,j);
                acsolio(i,j)=0.0;
                acsolo(i,j)=kd((i),j)*csol(i,j);
                acsolmo(i,j)=acsolo(i,j);
            end
            %%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end

 