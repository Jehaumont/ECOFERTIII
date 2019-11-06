function balance_error = check_balance(wat_flxs, WC, WCB, rtex, dt, dx, ncs, err_tol)
												
%Calculation of the mass balance error
%
%IN: 
%	wat_flxs
%	WC, WCB
%	rtex
%   dt,dx,compartiments_number,err_tol
%OUT:
%   balance_error
%CALL:
%	none
%CALLED BY:
%   solve_flow.m
%--------------------------
% M. Vanclooster 2/2/2000

%initialization

diff_th=(WC-WCB)/dt;
diff_flux=(-wat_flxs(1:ncs)+wat_flxs(2:ncs+1))/dx;
err= abs(diff_th-diff_flux+rtex);
err_max = max(err);

if err_max > err_tol
   balance_error = 1;
else
   balance_error = 0;
end
