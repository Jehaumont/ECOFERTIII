function [soilInnerStateParams,soilOuterStateParams] = change_temp(fract1_new, fract2_new,...
            soilInnerStateParams,soilOuterStateParams)

%% FUNCTION INPUT
temp_node1_new = soilOuterStateParams.temp_node_new;
temp_node2_new = soilInnerStateParams.temp_node_new;

temp1 = soilOuterStateParams.temp;
temp2 = soilInnerStateParams.temp;

%% FUNCTION MAIN BODY

if fract1_new==1
temp_node1_new(:,end+1) = temp1';
temp_node2_new(:,end+1) = temp1'*0;

elseif fract2_new==1
temp_node2_new(:,end+1) = temp2';
temp_node1_new(:,end+1) = temp2'*0;
else
temp_node2_new(:,end+1) = temp2';
temp_node1_new(:,end+1) = temp1';
end

%% FUNCTION OUTPUT
soilOuterStateParams.temp_node_new = temp_node1_new;
soilInnerStateParams.temp_node_new = temp_node2_new;