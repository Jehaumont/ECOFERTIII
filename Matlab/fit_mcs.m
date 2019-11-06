%function fitting_glob()
%FITTING_GLOB Resolve the inverse problem for unsaturated flow using global optimisation.
%	Soil parameters are estimated by minimizing an objective function.
%	(general least square problem).
	
%	INPUTS: function to edit before running FITTING_GLOB
%		.DATA_INV_SOL define observed data
%		.PARAM_DATA 	define fixed and fitted parameters
%							define parametric space in wich the optimization is carried out
%		.see inputs of function MAIN to define the direct problem
%		.FITTING_GLOB 	define the arguments of mcs to improve the inversion

%	OUTPUTS:
%		.xbest: 		estimated parameters
%		.IC: 			+- 95% confidence interval
%		.C: 			parameters covariance matrix
%		.A: 			parameters correlation matrix
%		.ssq: 		error sum of square (SSQ)
%		.rmse: 		root mean square error (RMSE)
%		.corr_coef: linear correlation coefficient (r²)

%	REFERENCES
%		Kool, J.B., and J.C. Parker. 1988.
% 		Analysis of the inverse problem for transient unsaturated flow.
% 		Water Resources Research, 24(6) : 817-830.
%
%		Huyer, W., A. Neumaier. 1998.	 
%		Global minimization by multilevel coordinate search.
%  	Institut für Mathematik, Universität Wien, Austria.

%	S. Lambot (October 2000)


% Global variables shared with OBJ_FUN
global DATA_INV t_obs %weigth_ph weigth_wc weigth_fl
global PARAM
global Ymodel Ymeas t_obs ERR WEIGHT
global cpt nf

% Read observed data
disp('Loading of measured data ...');
load('nl30')
t_obs=nl30(:,1);
wc_node_obs=nl30(:,2);
[DATA_INV] = wc_node_obs';
Ymeas = DATA_INV(:,1);
%[weigth_ph,weigth_wc,weigth_fl] = calc_weigth(DATA_INV);

% Read parameters data
[PARAM,U,V,DIS1,Qp,NM,smax,nf,stop,iinit,local,gamma] = param_data;
np = length(U);

% Preliminary calculations
cpt = 0;
warning off;

% Global optimization using Genetic Algorithm
%disp('Global optimization (Genetic Algorithm) ...');
%bounds = [U V];
%[xbest1,endPop,bestSols,trace] = ga(bounds,'obj_fun_ga');

% Global optimization by multilevel coordinate search (mcs) : Scan 1
disp('Global optimization (MCS) ...');
[xbest1,fbest,xmin,fmi,ncall,ncloc,flag] = ...
   mcs('feval','obj_fun',U,V,0,smax,nf,stop,iinit,local,gamma,ones(np,np));

% Local optimization using Nedler-Mead algorithm
disp('Local optimization (Nelder-Mead) ... ');
xbest2 = fminsearch('obj_fun',xbest1);

% STATISTICS
disp('Statistics calculation ...');

% Preliminary calculation
[OF] = obj_fun(xbest2);
ssq = OF;
Ymodelfit = Ymodel;
ERRfit = ERR;
p = size(xbest2 ,1);
n = size(Ymodel,1);
WEIGHT = WEIGHT*diag(ones(size(Ymodel))); % only if simplification in obj_fun is done

% Calculate the jacobian (J)
J = zeros(n,p);
for k = 1:p
   delta = zeros(p,1);
   delta(k) = xbest2(k)*0.01;
   [OF] = obj_fun(xbest2+delta);
   if OF == 1
      [OF] = obj_fun(xbest2-delta);
   end
   J(:,k) = (Ymodel-Ymodelfit)/delta(k);
end

% Calculate parameters covariance matrix (C)
S = chol(WEIGHT);	%Cholesky decomposition
Jw = S*J;
ERRw = S*ERR;
ERR_var = (ERRw'*ERRw)/(n-p);
C = ERR_var*inv(Jw'*Jw);		
Cii = diag(C);

% Calculate parameters correlation matrix (A)
Ci = Cii*ones(1,p);
Cj = Ci';
A = C./(sqrt(Ci).*sqrt(Cj));

% Calculate +- 95% confidence interval (IC)
IC = [sqrt(Cii)*qt(0.975,n-p)];

% Calculate root mean square error (rmse)
rmse = sqrt((ERRw'*ERRw)/n);

% Calculate correlation coefficient
corr_coef = corrcoef(Ymodelfit,Ymeas);

% Display results
disp('Estimated parameters [MCS  Nelder-Mead]');	disp([xbest1 xbest2]);
disp('+- 95% confidence interval');						disp(IC);
disp('Parameters covariance matrix'); 					disp(C);
disp('Parameters correlation matrix'); 				disp(A);
disp(['SSQ = ' num2str(ssq)]);
disp(['RMSE = ' num2str(rmse)]);
disp(['r² = ' num2str(corr_coef(2,1)^2)]);

% Delete memory consuming variables
clear S;
clear WEIGHT;