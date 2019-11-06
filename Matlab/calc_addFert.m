%%%%%---- CALCULATE AMOUNT OF FERTILIZATION NEEDED -----%%%%%
% CALLED by: wavemat105_sol
% CALLS: 
%
% J. Haumont 06/03/2019

function [soilInnerStateParams,soilOuterStateParams] = calc_addFert(cropStateParams,...
                                                                    managementSettings,...
                                                                    simulationSettings,...
                                                                    soilCommonStateParams,...
                                                                    soilInnerStateParams,...
                                                                    soilOuterStateParams)
%%%%%%% DOCUMENTATION %%%%%%%
%%%%% INPUT %%%%%
% reaminNmin: remaining N concentration in soil [g/cm²] 
% method: method for calculating TV
%         "N_pre","KNS","KNS60","FAD","BD"
% cropType: which crop is planted 
%           1: cauliflower
%           2: leek
% timeNextFert: date which a next fertilization moment is planned [Julinan day[]
%               = (i+1)-th fertilization (if another moment in cycle planned)
%               = harvest date           (if no other moment in cycle planned)
% strategy: fertilization method
%           "Broad","Fert","Row"

%%%%% OUTPUT %%%%%
% addFert: amount of fetilization needs to be applied during the next day
%          [g/cm²]

%% FUNCTION INPUT

% management settings
KNS_miner_correction = managementSettings.KNS_miner_corr;

fertMoment = [reshape(managementSettings.fertMomentCauli,...
                      numel(managementSettings.fertMomentCauli),1) + ...
                      simulationSettings.tstart,...
              ones(numel(managementSettings.fertMomentCauli),1);...
              reshape(managementSettings.fertMomentLeek,...
                      numel(managementSettings.fertMomentLeek),1) + ...
                      simulationSettings.tstart,...
              ones(numel(managementSettings.fertMomentLeek),1).*2];
          
[~,idx] = sort(fertMoment(:,1));
fertMoment = fertMoment(idx,:);
harvest_date = sort([managementSettings.hDateCauli;managementSettings.hDateLeek]);
method = managementSettings.fertDecisionRule;

plant_date = sort([managementSettings.pDateCauli; managementSettings.pDateLeek]);

t = simulationSettings.t;

tcsolo = soilCommonStateParams.tcsolo; %[g/cm²]


%% FUNCTION MAIN BODY
% FIND THE END OF THE TIME WINDOW THAT SHOULD BE REGARDED FOR THE
% CALCULATION OF THE FERTILIZER ADDITION
if any(fertMoment(:,1)-t == 1)
    
    idx = find(fertMoment(:,1)-t ==1); % the current time should be only 1 day before the fertilization moment
    cropType = fertMoment(idx,2);
    timeNextHarvest = harvest_date(find(t+1>=plant_date & t<=harvest_date,1,'last')); % find the date of the next harvest
    
    try
        timeNextFert = fertMoment(idx+1,1); % find the date of the next moment to apply fertilizer
    catch
        timeNextFert =  harvest_date(find(t+1 >= plant_date & t<=harvest_date,1,'last'))+1;
    end
    
    if timeNextFert > timeNextHarvest
        
        endTimeWindow = timeNextHarvest; % end of the time window is the end of the interval for which the calculation of the fertilizer application should account
        
    elseif timeNextFert < timeNextHarvest
        
        endTimeWindow = timeNextFert;
    end
end
    
% DETERMINE THE TARGET VALUE FOR FERTILIZATION
if cropType == 1 % if cauliflower is on the field    
    switch method
        case "FIX3"
            % fixed target values no consideration of current soil nitrogen
            % content
            soilSampleFlag = 0; 
            TV = [130 25];
            if endTimeWindow < timeNextHarvest
                index = 1;
                
            elseif endTimeWindow == timeNextHarvest
                index = 2;
            end
            TV = TV(index);
        
        case "FIX3_50min"
            soilSampleFlag = 0;
            TV = [65, 12.5];
            
            if endTimeWindow < timeNextHarvest
                index = 1;
            elseif endTimeWindow == timeNextHarvest
                index = 2;
            end
            TV = TV(index);
        
        case "FIX3_50plus"
            soilSampleFlag = 0;
            TV = [195, 37.5];
            
            if endTimeWindow < timeNextHarvest
                index = 1;
            
            elseif endTimeWindow == timeNextHarvest
                index = 2;
            end
            
            TV = TV(index);
            
        case "KNS"
            soilSampleFlag = 1;
            
            TV = [219, 218, 214, 308, 287, 245, 181, 124, 92];
            soilSampleDepth = [30, 30, 30, 60, 60, 60, 60, 60, 60];
            sampleTime = 1:1:8;
            
            idx_pDate = find( t+1-plant_date >= 0, 1, 'last');
            DAP = t - plant_date(idx_pDate)+1;
            
            if DAP~=0
                week_num = ceil(DAP/7);
            elseif DAP == 0
                week_num = 1;
            end
            
            if week_num > sampleTime(end)
                error('The crop should be fertilized earlier')
            end         
            
            TV = TV(week_num);
            soilSampleDepth = soilSampleDepth(week_num);
            
            % Correct TV for mineralisation (0.8kgN/day)
            TV = TV-KNS_miner_correction*(timeNextHarvest-t);
            
        case "N_pre"
            NUptakePredicted = PBLiski(endTimeWindow,"gompInt",strategy,t);
            TV = NUptakePredicted(end,1);
        
        otherwise
            error(['Error: Method for calulating the needed fertilizer'...
                   'is unknown. Check calc_addFert for possible methods'])
    end
    
    
elseif cropType == 2 % if leek is on the field
    
    switch method
        case "FIX3"
            % fixed target values no consideration of current soil nitrogen
            % content
            soilSampleFlag = 0; 
            TV = 100;
        
        case "FIX3_50min"
            soilSampleFlag = 0;
            TV = 50;
            
        case "FIX3_50plus"
            soilSampleFlag = 0;
            TV = 150;
        
        case "KNS"
            soilSampleFlag = 1;
            TV = [159, 159, 158, 156, 153, 148, 256, 243, 224, 196, 157,...
                  117, 87, 70];
            soilSampleDepth = [30, 30, 30, 30, 30, 30, 60, 60, 60, 60, 60,...
                               60, 60, 60];
            
            sampleTime = 1:1:14;
            
            idx_pDate = find( t+1 - plant_date >= 0, 1, 'last');
            DAP = t - plant_date(idx_pDate) + 1;
            
            if DAP~=0
                week_num = ceil(DAP/7);
            elseif DAP == 0
                week_num = 1;
            end
            
            if week_num > sampleTime(end)
                error('The crop should be fertilized earlier')
            end
            
            TV = TV(week_num);
            soilSampleDepth = soilSampleDepth(week_num);           
            
            % Correct TV for mineralisation (0.8kgN/day)
            TV = TV-KNS_miner_correction*(timeNextHarvest-t);
            
        case "N_pre"
            NUptakePredicted = PBLiski(endTimeWindow,"gompInt",strategy,t);
            TV = NUptakePredicted(end,1);
            
        otherwise
            error(['Error: Method for calulating the needed fertilizer'...
                   'is unknown. Check calc_addFert for possible methods'])
    end
end

% CALCULATE THE REMAINING AMONIUM AND NITRATE CONCENTRATION IN THE PROFILE
if soilSampleFlag == 1
    
    NH4conc = sum(tcsolo(1:soilSampleDepth,2))*10^5; % convert g/cm² to kg/ha
    NO3conc = sum(tcsolo(1:soilSampleDepth,3))*10^5; % convert g/cm² to kg/ha
    Nremaining = NH4conc + NO3conc; % total remaining nitrogen concentration available for uptake

elseif soilSampleFlag == 0
    % no soil sample is taken and thus apparantly there's no remaining nitrogen
    Nremaining = 0; 
end

addFert = (TV-Nremaining)/10^5; % convert kg/ha to g/cm²

if addFert < 0
    addFert = 0;
end

if isnan(addFert)
    disp('')
end

switch managementSettings.fertStrategy
    
    case "Broad"
        soilInnerStateParams.fsol(end+1,:) = [simulationSettings.t+1, 0, addFert/2, addFert/2]; % divide the fertilizer amount equally over ammonium and nitrate
        soilOuterStateParams.fsol(end+1,:) = [simulationSettings.t+1, 0, addFert/2, addFert/2]; % also apply fertilizer outside rooted area
        disp([num2str(addFert),' g/cm² is added to the soil'])
    
    case "Row"
        
        disp([num2str(addFert),' g/cm² is added to the rooted soil'])        
        soilInnerStateParams.fsol(end+1,:) = [simulationSettings.t+1, 0, addFert/2, addFert/2]; % only apply fertilizer to the rooted area
        soilOuterStateParams.fsol(end+1,:) = [simulationSettings.t+1, 0, 0, 0]; % set fertilizer application outside rooted area to 0, should be specified for correct boundary conditions
    otherwise
        error("The specified fertilization strategy is not implmented, please choose between 'Broad' or 'Row'")
end

end

%%%%%----- HELPER FUNCTIONS -----%%%%%
function NUptakePredicted = PBLiski(endTimeWindow,growthCurve,strategy,cropStateVars,t)
%%%%%% DOCUMENTATION %%%%%%
%%%% INPUT %%%%%
% DAP: days after planting [-]
% endTimeWindow: date which a next fertilization moment is planned [Julinan day]
%              = (i+1)-th fertilization (if another moment in cycle
%                                        planned)
%              = harvest date           (if no other moment in cycle
%                                        planned)
% growthCurve: which type of growth curve should be fitted to the data
%              gomp: gompertz curve with fixed intercept 0
%              gompInt: gompertz curve with flexible intercept
%              fplr: four parameter logistic
%              tplr: tree parameter logistic 
%%%%% OUTPUT %%%%%
% NUptakePredicted: predicted amount of N that will be taken up until the
%                   next fertilization moment or harvest moment [g/cm²]

%% define the right growth curve function and the derivatives
switch growthCurve
    case "gomp"
        fun = @(param,t) param(1).*exp(-param(2).*param(3).^t);
        %  dgdi is partial derivative of fun for parameter i 
        deriv{1,1} = "dfp1";
        deriv{2,1} = "dfp2";
        deriv{3,1} = "dfp3";
        
        deriv{1,2} = @(param,t) exp(-param(2).*param(3).^t);
        deriv{2,2} = @(param,t) -param(1).*param(3).^t.*exp(-param(2).*param(3).^t);
        deriv{3,2} = @(param,t) -param(2).*param(2).*t.*param(1).*param(3).^(t-1) .*...
                                exp(-param(2).*param(3).^t);
        nParams = 3;
    case "gompInt"
        fun = @(param,t) param(4) + (param(1)-param(4)).*exp(-param(2).*...
                         param(3).^t);
        %  dgdi is partial derivative of fun for parameter i
        deriv{1,1} = "dfp1";
        deriv{2,1} = "dfp2";
        deriv{3,1} = "dfp3";
        deriv{4,1} = "dfp4";
        
        deriv{1,2} = @(param,t) exp(-param(2).*param(3).^t);
        deriv{2,2} = @(param,t) (param(4)-param(1)).*param(3).^t.*exp(-param(2).*...
                                 param(3).^t);
        deriv{3,2} = @(param,t) -param(2).*t.*(param(4)-param(1)).*param(3).^(t-1)...
                                 .* exp(-param(2).*param(3).^t);
        deriv{4,2} = @(param,t) 1-exp(-param(2).*param(3).^t);
        nParams = 4;
    case "fplr"
        error('Error: This growth curve is not yet available')
    case "tplr"
        error('Error: This growth curve is not yet available')
    otherwise
        error('Error: Cannot make predictions with this growth curve')
end

%% Get initial information for mixed effect model
% load the database with all N-uptake curves for 43 seasons based on the
% chosen strategy. Define starting values for the mixed model (beta0) based
% on analysis in R 
switch strategy
    case "Broad"
        
        dataBase = readtable(['C:\Users\U0123281\Box Sync\ECOFERT_eddie\',...
                              'EcofertI_JVL\Chapter 6\6.1\Results\',...
                              'test_NfacTiming\Leek_output_tAd.csv']);
        beta0 = [197.203,20.432,0.946,17.393]; 
    case "Fert"
        error('Error: This strategy is not yet available')
    case "Row"
        error('Error: This strategy is not yet available')
    otherwise
        error('Error: The defined strategy is not included in the model')
end

if t<max(dataBase.DATE)  % check wether the simulation year is already in the database
                         % In this case the database should be split in
                         % 'observed' and 'unobserved' data
                        
    yearNewObs = dataBase.Year(dataBase.DATE==t);
    DAP = dataBase.DAP(dataBase.DATE==t);
    
    % apparant observed individuals 
    datObs = dataBase(dataBase.Year~=yearNewObs,["Year","DATE","DAP","TotNUptakeperHa"]);
    
    % apparant new observed individuals
    datNew = dataBase(dataBase.Year==yearNewObs,["Year","DATE","DAP","TotNUptakeperHa"]);
    nObsNew = DAP;
    datPast = datNew(1:nObsNew,:);
    
elseif t> max(dataBase.DATE)
    
    datObs = dataBase(:,["Year","DATE","DAP","TotNUptakeperHa"]);
    yearNewObs = unique(cropStateVars.Year); 
    datPast = table(cropStateVars.Year,cropStateVars.DATE,cropStateVars.DAP,...
                    cropStateVars.TotNUptakeperHa);
    datPast.Properties.VariableNames = {'Year','DATE','DAP','TotNUptakeperHa'};
    nObsNew = length(cropStateVars.TotNUptakeperHa);
end

years = unique([datObs.Year;datPast.Year]);

%% Fit non linear mixed model with fixed effects coFe and random effects coRa
options = statset('FunValCheck','off','Display','iter');
CovPattern = ones(nParams);
LogLikMethod = "gq";
NIterations = [150 150 100];

[coFe,PSI,stats,coRa] = nlmefitsa([datObs.DAP;datPast.DAP],...
    [datObs.TotNUptakeperHa;datPast.TotNUptakeperHa],...
    [datObs.Year;datPast.Year],[],fun,beta0,'Options',options,...
    'CovPattern',CovPattern,'LogLikMethod',LogLikMethod,'NIterations',NIterations);

%% Apply PBLiski prediction method 
xiBar = cell(1,nParams); % pseudo data of the fully observed individuals
xoBar = cell(1,nParams); % pseudo data of the new observed individual
xiBarYear = [];
xoBarYear = [];

% Calculate the pseudo data xiBar xoBar of the explenatory variable
for i=1:length(years)
    tmpParams = coFe + coRa(:,i);
    if years(i) ~= yearNewObs
        DAPobs = datObs.DAP(datObs.Year==years(i));
        nObs = length(DAPobs);
        for j=1:size(deriv,1)
            xiBar{1,j} = [xiBar{1,j};deriv{j,2}(tmpParams,DAPobs)];
        end
        xiBarYear = [xiBarYear,repmat(years(i),1,nObs)];
    else
        DAPobs = 1:endTimeWindow;
        for j=1:size(deriv,1)
            xoBar{1,j} = [xoBar{1,j};deriv{j,2}(tmpParams,DAPobs')];
        end     
    end
end

xiBar = cell2mat(xiBar);
xoBar = cell2mat(xoBar);
xiBarYear = xiBarYear';

Yi = datObs.TotNUptakeperHa;
xi = datObs.DAP;
Yo = datPast.TotNUptakeperHa;
xo = 1:endTimeWindow;
xoP = xo(1:nObsNew);
xoF = xo(nObsNew+1:end);

ai = coRa(:,1:end-1) + coFe;
ao = coFe + coRa(:,end);

zi = [];
j = 0;
% calculate the pseudo data Zi of the reponse variable
for i=1:length(years)
    if years(i) ~= yearNewObs
        
        j=j+1;
        DAP = datObs.DAP(datObs.Year==years(i));
        tmpYi = datObs.TotNUptakeperHa(datObs.Year==years(i));
        xiBarSubset = xiBar(xiBarYear==years(i),:);
        a = ai(:,j);
        tmpZi = tmpYi - fun(a,DAP) + xiBarSubset*a;
        zi = [zi;tmpZi];
    end
end

xoFbar = xoBar((nObsNew+1):size(xoBar,1),:);
xoPbar = xoBar(1:nObsNew,:);
zo = Yo - fun(ao,xoP) + xoPbar*ao;


YoFmax = []; % upper limit of the prediction interval 
YoFmin = []; % lower limit of the prediction interval
YoF = [];    % predicted response variable

VoP = xoPbar*PSI*xoPbar'+ stats.rmse^2*eye(nObsNew);
j=0;
covFe = zeros(nParams);
for i=1:length(years)
    if(years(i)~=yearNewObs)
        j=j+1;
        xiBarSubset =  xiBar(xiBarYear==years(i),:);
        Vi = xiBarSubset*PSI*xiBarSubset'+stats.rmse^2*eye(size(xiBarSubset,1));
        covFe = covFe + xiBarSubset'*(Vi\xiBarSubset) + xoPbar'*(VoP\xoPbar);
    end
end

% calculate future values YoF of the response variables and the prediction
% intervals YoFmin, YoFmax

for k = 1:(endTimeWindow-nObsNew)
    VoFP = xoFbar(k,:)*PSI*xoPbar';
    ZoF = xoFbar(k,:)*coFe + VoFP*(VoP\(zo(1:nObsNew,:)-xoPbar*coFe));
    YoF(k) = ZoF + fun(ao,xoF(k)) - xoFbar(k,:)*ao;
    
    VoF = xoFbar(k,:)*PSI*xoFbar(k,:)' + stats.rmse^2;
    OmegaF = VoF -VoFP*(VoP\VoFP');
    Mo = xoFbar(k,:)- VoFP*(VoP\xoPbar);
    VarYoF = OmegaF + Mo*(covFe\Mo');
    critT = tinv(0.975,500);
    YoFmin = [YoFmin,YoF(k) - critT*sqrt(VarYoF)];
    YoFmax = [YoFmax,YoF(k) + critT*sqrt(VarYoF)];
end

YoP = datPast.TotNUptakeperHa;
Y = [YoP',YoF];
YoFmin = [YoP',YoFmin];
YoFmax = [YoP',YoFmax];

NUptakePredicted = [Y,YoFmin,YoFmax];
end
