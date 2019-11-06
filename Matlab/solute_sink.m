function[decsolm,decsoli,sol_sinkm,sol_sinki]=solute_sink(dt,nsol,ncs,cm,cim,wcmo,...
    wcio,solute_param1,solute_param2,immobile);

% Calculate the sink term without without considering nitrogen module
%IN:
%dt
%nsol and ncs= number of solutes and number of compartiments
%cm and cim= solute concentration in the mobile and immobile region
%wcmo and wcio= water content in the mobile and immobile region (from script%wat_sol)
%solute_param1 and solute_param2= solute parameters  (from script%in_solute_parameters)
%immobile= use or not use immobile concept(from script%in_solute_parameters)
%OUT:
%decsolm and decsoli= amount of decayed solute in the mobile and immobile %region
%sol_sinkm and sol_sinki= sink term for the mobile and immobile region
%CALL: none
%CALLED BY: solve solute
%----------------------------------------------------------------
%M. Sall 17/03/09

%solute parameters 
bd=solute_param1(:,1);
alf=solute_param1(:,4);
f=solute_param1(:,5);
rates=solute_param2(:,1:nsol);
kd=solute_param2(:,nsol+1:2*nsol);

%calculate current decayed solute without acconting the mobile/immobile %transfer
for sp=1:nsol
    for i=1:ncs
        decsolm(i,sp) = -rates(i,sp).*cm(i,sp).*(wcmo(i));
        decsoli(i,sp) = -rates(i,sp) .*cim(i,sp) .*wcio(i) ;
    end
end
%calculate the sink term
for sp =1:nsol
    for i=1:ncs-1
 %%the mobile sink term
        sol_sinkm(i,sp)=decsolm(i,sp)-alf(i)*(cm(i,sp)-cim(i,sp)) ;   
		defsolm=cm(i,sp)*(wcmo(i)+f((i))*bd((i))*kd((i)))/dt+sol_sinkm(i,sp);                        
		totsolm=-sol_sinkm(i,sp);      
 %% reduction of the mobile sink term
		if defsolm <0 & totsolm > 0                 
			facsolm= max(0, 1+defsolm/totsolm);       
			sol_sinkm(i,sp)= sol_sinkm(i,sp)*facsolm;             
			decsolm(i,sp)=decsolm(i,sp)*facsolm;                      
        end
%%	the immobile sinkterm           
        if immobile
            sol_sinki(i,sp)=decsoli(i,sp);  
            defsoli=cim(i,sp)*(wcio(i)+f((i))*bd((i))*kd((i)))/dt+sol_sinki(i,sp);                        
            totsoli=-sol_sinki(i,sp);       
 %% reduction of the immobile sink term
                if defsoli <0 & totsoli > 0                 
                    facsoli= max(0, 1+defsoli/totsoli);       
                    sol_sinki(i,sp)= sol_sinki(i,sp)*facsoli;             
                    decsoli(i,sp)=decsoli(i,sp)*facsoli;                      
                end
        else
             decsoli(i,sp)=0;
             solsinki(i,sp)=0;
        end
    end
end
% for the botom compartiment 
sol_sinkm(ncs,:)=zeros(1,nsol);
sol_sinki(ncs,:)=zeros(1,nsol);
decsolm(ncs,:)=zeros(1,nsol);
decsoli(ncs,:)=zeros(1,nsol);




