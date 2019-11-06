function [soilOuterStateParams,soilInnerStateParams] = change_tc(fract1_new, fract2_new,...
    soilOuterStateParams,soilInnerStateParams)


%% FUNCTION INPUT
tc_node1_new = soilOuterStateParams.tc_node_new;
tc_node2_new = soilInnerStateParams.tc_node_new;

tcsolo1 = soilOuterStateParams.tcsolo;
tcsolo2 = soilInnerStateParams.tcsolo;

%% FUNCTION MAIN BODY
for i=1:3
    if fract1_new==1
        tc_node1_new{1,i}(:,end+1) = tcsolo1(:,i);
        tc_node2_new{1,i}(:,end+1) = tcsolo1(:,i)*0;
        
    elseif fract2_new==1
        tc_node2_new{1,i}(:,end+1) = tcsolo2(:,i);
        tc_node1_new{1,i}(:,end+1) = tcsolo2(:,i)*0;
        
    else
        tc_node2_new{1,i}(:,end+1) = tcsolo2(:,i);
        tc_node1_new{1,i}(:,end+1) = tcsolo1(:,i);
    end
end

%% FUNCTION OUTPUT

soilOuterStateParams.tc_node_new = tc_node1_new;
soilInnerStateParams.tc_node_new = tc_node2_new;
