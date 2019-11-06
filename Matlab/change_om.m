function [soilInnerStateParams, soilOuterStateParams] = change_om(fractSoil_new, fractPlant_new,...
    soilOuterStateParams,soilInnerStateParams)

%% FUNCTION INPUT
om_node1_new = soilOuterStateParams.om_node_new;
om_node2_new = soilInnerStateParams.om_node_new;

soil_om1_new = soilOuterStateParams.soil_om;
soil_om2_new = soilInnerStateParams.soil_om;

%% FUNCTION MAIN BODY
for i=1:7
    if fractSoil_new==1
        om_node1_new{1,i}(:,end+1) = soil_om1_new(:,i);
        om_node2_new{1,i}(:,end+1) = soil_om1_new(:,i)*0;
        
    elseif fractPlant_new==1
        om_node2_new{1,i}(:,end+1) = soil_om2_new(:,i);
        om_node1_new{1,i}(:,end+1) = soil_om2_new(:,i)'*0;
    else
        om_node2_new{1,i}(:,end+1) = soil_om2_new(:,i);
        om_node1_new{1,i}(:,end+1) = soil_om1_new(:,i);
    end
end

%% FUNCTION OUTPUT
soilOuterStateParams.om_node_new = om_node1_new;
soilInnerStateParams.om_node_new = om_node2_new;