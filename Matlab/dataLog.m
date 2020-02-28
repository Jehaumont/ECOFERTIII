function structure = dataLog(structure , type, fileSettings, force, fieldNames)

if nargin == 3
    force = false;
    fieldNames = 'all';

elseif nargin == 4 
    fieldNames = 'all';
end

% set field names that should be saved
if strcmp(fieldNames,'all')
    fields = fieldnames(structure);
else
    if ~(class(fieldNames) == "cell")
        if class(fieldNames) == "string"
            fields = arrayfun(@(x)char(fieldNames(x)), 1:numel(varNames),...
                            'uni', false);
        else
            error('The specified field names should entered as a 1xN cell or string array');
        end
    else
        fields = fieldNames;
    end
end                          
    
resultsFolder = string(fileSettings.ResultsPath);
DataLog_size = fileSettings.DataLog_size;

for i=1:length(fields)

    % save if data size exceeds datalog block size
    if size(structure.(fields{i}),2) >= DataLog_size & ~strcmpi(force,"force") 
        data = structure.(fields{i})(:,1:end-2);
        
        % try to save data to file
        try
            path_to_file = resultsFolder + type + "\" + fields(i) + ".txt";
            dlmwrite(path_to_file, data','-append','newline','pc','precision',6);
        
        % if error make directory and save
        catch
            saveFolder = resultsFolder + type + "\";
            mkdir(saveFolder)
            dlmwrite(path_to_file, data','-append','newline','pc','precision',6);
        end
        
        % keep only the last 2 data points
        structure.(fields{i}) = structure.(fields{i})(:,end-1:end);
    
    % save data when forced
    elseif strcmpi(force,"force")
        path_to_file = resultsFolder + type + "\" + fields(i) + ".txt";
        data = structure.(fields{i});
        dlmwrite(path_to_file, data','-append','newline','pc','precision',6);
    end    
end


