function  [csol, cm, cim, acsolmo, acsolio,ppdepth,tcsolo] = plowing(t,solute_param1, solute_param2,csol, cm, cim, acsolmo,...
     acsolio, wcio, wco, nsol, immobile,wcmo,tcsolo,dx,ppdepth,solsur)
 
%Function describing redistribution of solute due to plowing

%Vansteenkiste J. 16/03/2010

%Read in parameters
bd=(solute_param1(:,1));
f=(solute_param1(:,5));
kd=(solute_param2(:,nsol+1:2*nsol));

%Find if there has been plowing
index = max(find(t>=ppdepth(:,1)));
current_plow = ppdepth(index,:);

%Set to zero
ppdepth(index,2:end)=0;


if current_plow(end)~=0 %There is plowing in the current time step
%Calculate the total amount of solute in (reservoir + plow layer)
%Redistribute all solute over the plow layer
pltsol = cumsum(tcsolo, 1);
pltsol = pltsol(current_plow(3),:);
pltsol = pltsol + solsur;

i = [1:current_plow(3)];
tcsolo(i,:) = ones(current_plow(3),1)*(pltsol./current_plow(3));
tcsolio(i,:) = tcsolo(i,:).*repmat(wcio(i)',1,nsol)./repmat(wco(i)',1,nsol);
tcsolmo(i,:) = tcsolo(i,:).*repmat(wcmo(i)',1,nsol)./repmat(wco(i)',1,nsol);
if immobile
    disp('still to  be programmed')
else
csolio(i,1:nsol) =0;
end
cim(i,:) = csolio(i,:);
csolmo(i,:) = tcsolmo(i,:)./(dx* kd(i,:).*repmat(bd(i),1,nsol).*repmat(f(i),1,nsol) + repmat(wcmo(i)',1,nsol));
cm(i,:) = csolmo(i,:);
csolo = csolio(i,:).*repmat(wcio(i)',1,nsol)+csolmo(i,:).*repmat(wcmo(i)',1,nsol)./repmat(wco(i)',1,nsol);
if immobile
    disp('still to be programmed')
else
    acsolio(i,:) = 0;
end
acsolmo(i,:) = kd(i,:).*repmat(f(i),1,nsol).* csolmo(i,:);
acsolo = acsolio(i,:)+ acsolmo(i,:);
 
end