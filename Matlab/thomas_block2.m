function ph = thomas_block2(ph,WC,WCB,kh,kh_in,CH,rtex,dt,DX,dx_inter,...
   phsurf,phbot,flxsbot,flxar,boco_top_type,boco_bot_type,ncbot);

%Numerical inversion of the 3diagonal matrix to find ph for t+dt%
%IN:
%   ph,WC,WCB,kh,kh_in,CH,rtex,dt,DX
%   phsurf,phbot,flxsbot,flxar,compartiments_number
%   boco_top_type,boco_top,boco_bot_type,boco_bot
%OUT:
%   ph
%CALL:
%   none
%CALLED BY: 
%   solveflow.m
%
%------------------------------
%modified by Javaux 11/04/00


if boco_top_type == 1 
   h0=dt/DX(1);
	h1=h0/dx_inter(1);
	h2=h0/dx_inter(2);
	a=h2*kh_in(2);
   c=h1*kh_in(1);
   A(1,1)=-a;
   B(1,1)=CH(1)+a+c;
   E(1,1)=c*phsurf(1)+CH(1)*ph(1)+h0*(kh_in(1)-kh_in(2))-dt*rtex(1)-WC(1)+WCB(1);
end 
     
if boco_top_type == 2 
	h0=dt/DX(1);
 	h1=h0/dx_inter(1);
	h2=h0/dx_inter(2);
   a=h2*kh_in(2);
   A(1,1)=-a;
   B(1,1)=CH(1)+a;
   E(1,1)=CH(1).*ph(1)-h0.*(flxar+kh_in(2))-dt*rtex(1)-WC(1)+WCB(1);
end

% intermediate values
 I=[2:ncbot-1];%AJOUT
h0=dt./DX(I);
h1=h0./dx_inter(I);
h2=h0./dx_inter(I+1);
a=h2.*kh_in(I+1);
c=h1.*kh_in(I);
A(I,1)=-a.';
B(I,1)=(CH(I)+a+c).';
D(I-1,1)=-c.';
E(I,1)=(CH(I).*ph(I)+h0.*(kh_in(I)-kh_in(I+1))-dt*rtex(I)-WC(I)+WCB(I)).';

%  bottom values
if boco_bot_type == 1 | boco_bot_type==5;
   h0=dt/DX(ncbot);
   h1=h0/dx_inter(ncbot);
	a=h0/dx_inter(ncbot+1)*kh_in(ncbot+1);%
	c=h1*kh_in(ncbot);
   B(ncbot,1)=CH(ncbot)+a+c;
   D(ncbot-1,1)=-c;
	E(ncbot,1)=CH(ncbot)*ph(ncbot)-WC(ncbot)+WCB(ncbot)+a*phbot+h0*(kh_in(ncbot)-kh_in(ncbot+1))-dt*rtex(ncbot);
end

if boco_bot_type == 2 | boco_bot_type == 4 | boco_bot_type == 6 | boco_bot_type == 7
	h0=dt/DX(ncbot);
   h1=h0/dx_inter(ncbot);
	c=h1*kh_in(ncbot);
   B(ncbot,1)=CH(ncbot)+c;
   D(ncbot-1,1)=-c;
	E(ncbot,1)=CH(ncbot)*ph(ncbot)-WC(ncbot)+WCB(ncbot)+h0*(kh_in(ncbot)+flxsbot)-dt*rtex(ncbot);
end   

% ph Matrix: inversion
ph(1:ncbot)=(inv(diag(B)+diag(A,1)+diag(D,-1))*E).';


