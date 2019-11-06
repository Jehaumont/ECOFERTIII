function cm = inv_sol_matrix(immobile,ncs,cs,f,dt, pveloh, wat_flxsah, wcmo, wcmob, dx, diffush, wcmah, pvelah,betha1,betha2, betha3,betha4, csp, cmp, bd, kd, sol_sinkm,nsol)

%Solve solute matrix

  for j=[1:nsol];                     
%TOP VALUES
     I=1;
ab(1,1)=  dt*pveloh(1)*wat_flxsah(1)*(wcmo(1)-wcmob(1))/(16*dx^2)+...        
    (diffush(1,j))*((wcmah(1)))/(2*dx^2);
bb(1,1)= dt*pveloh(1)*pvelah(2)*(wcmo(I)-wcmob(I))/(16*dx^2)+...
    (diffush((2),j))*(wcmah(2))/(2*dx^2);

cb(1,1)=wat_flxsah(1)/(2*dx);                           

db(1,1)=wat_flxsah(2)/(2*dx);

bl(1,1)= wcmo(1)/dt +ab(1,1)+bb(1,1)+(bd(1)*f(1)*kd(1,j))/dt-betha2(1)*db(1,1)+betha4(1)*cb(1,1);

cl(1,1)= -bb(1,1)-betha3(I)*db(1,1);
dl(1,1)=csp(j)*(ab(1,1)- betha1(I)*cb(1,1))-...
    cs(j)*(ab(1,1)+betha1(1)*cb(1,1))+...
    cmp(1,j)*(wcmob(1)/dt -ab(1,1)-bb(1,1)+(bd(I)*f(1)*kd(1,j))/dt+betha2(1)*db(1,1)-betha4(1)*cb(1,1))+...
    cmp((2),j)*(bb(1,1)+betha3(1)*db(1,1))+ sol_sinkm(1,j); % 


%INTERMEDIATE VALUES
   I= [2:ncs-1]; % for
       
ab(I,1)= (diffush(I,j))'.*(wcmah(I))./(2*dx^2) +...
dt*pveloh(I).*((pvelah(I))).*(wcmo(I)-wcmob(I))./(16*dx^2);

bb(I,1)= (diffush((I+1),j))'.*wcmah(I+1)./(2*dx^2) +... 
dt*pveloh(I).*((pvelah(I+1))).*(wcmo(I)-wcmob(I))./(16*dx^2);

cb(I,1)=wat_flxsah(I)./(2*dx);

db(I,1)=wat_flxsah(I+1)./(2*dx);

al(I-1,1)= -ab(I,1)+betha1(I).*cb(I,1);

bl(I,1)= (wcmo(I))'./dt +ab(I,1)+bb(I,1)+(bd(I).*f(I).*kd(I,j))/dt-betha2(I).*db(I,1)+betha4(I).*cb(I,1);

cl(I,1)= -bb(I,1)-betha3(I).*db(I,1);

dl(I,1)= cmp((I-1),j).*(ab(I,1) - betha1(I).*cb(I,1))+ ...
    cmp(I,j).*((wcmob(I))'/dt -ab(I,1)-bb(I,1)+(bd(I).*f(I).*kd(I,j))./dt+betha2(I).*db(I,1)-betha4(I).*cb(I,1))+...
    cmp((I+1),j).*(bb(I,1)+betha3(I).*db(I,1))+ sol_sinkm(I,j); 
     
    
 %BOTTOM VALUES
     I=ncs;
ab(ncs,1)= ((diffush(ncs,j)))*(wcmah(ncs))/(2*dx^2) +...
dt*pveloh(ncs)*((pvelah(ncs)))*(wcmo(ncs)-wcmob(ncs))/(16*dx^2);

bb(ncs,1)= 0;  %%%%%%%vmm(ncs+1/2)?

cb(ncs,1)=(wat_flxsah(ncs))/(2*dx);

db(ncs,1)=(wat_flxsah(ncs+1))/(2*dx);         

al(ncs-1,1)= -ab(ncs,1)+betha1(ncs)*cb(ncs,1); 

bl(ncs,1)= wcmo(ncs)/dt +ab(ncs,1)+bb(ncs,1)+(bd(ncs)*f(ncs)*kd(ncs,j))/dt -...
    betha2(ncs)*db(ncs,1)+betha4(ncs)*cb(ncs,1);

dl(ncs,1)=cmp((ncs-1),j)*(ab(ncs,1)-betha1(ncs)*cb(ncs,1)) +...
cmp(ncs,j)*(wcmob(ncs)/dt-ab(ncs,1)+(bd(ncs)*f(ncs)*kd(ncs,j)/dt)+...
    betha2(ncs)*db(ncs,1)-bb(ncs,1)-betha4(ncs)*cb(ncs,1))+...
    cmp(ncs+1,j)*(bb(ncs,1) +betha3(ncs)*db(ncs,1));

 %solution of the tridiagonal matrix

 cm(:,j)=(inv(diag(bl)+diag(cl,1)+diag(al,-1))*dl).'; %%
   

end 