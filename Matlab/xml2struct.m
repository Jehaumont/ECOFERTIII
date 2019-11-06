function param = xml2struct(filePath,variables)

%% DOCUMENTATION
% This function reads in the input parameters from the stored .xml files in
% the filePath. The input 'variables' can be used if you only want to load
% 1 or more specific variables from the .xml file. If more than 1 specific
% variable should be loaded the name should be stored in a string array
% filePath = folder path were all input files are stored as .xml files
% variables = string array of the specific variables names you want to load

% default value for variables is all --> all variables are loaded 
if nargin<2
    variables = 'all';
end

str = parseXML(filePath);
param = struct();
count = 0;
for i = {str.Children.Name}
    count = count+1;
    switch string(i)
        case 'Parameter' 
            data = str.Children(count).Children;
            name = data(2).Children.Data;
            dType = data(8).Children.Data;
            switch dType
                case 'numeric'
                    isEmpty = strcmp(data(4).Attributes.Value,'empty');
                    if isEmpty
                        value = [];
                    else
                        value = str2num(data(4).Children.Data);
                    end
                case 'cell'
                    if isempty(str2num(data(4).Children.Data))
                        value = strsplit(data(4).Children.Data,', ');
                    else
                        value = str2num(data(4).Children.Data);
                        cellSize = str2num(data(10).Children.Data);
                        arraySize = strsplit(data(4).Attributes.Value,'x');
                        arraySize{1} = str2num(arraySize{1});
                        arraySize{2} = str2num(arraySize{2});
                        arraySize = cell2mat(arraySize);
                        cols = repmat(arraySize(2)/cellSize(2),1,cellSize(2));
                        rows = repmat(arraySize(1)/cellSize(1),1,cellSize(1));
                        value = mat2cell(value,rows,cols);
                    end
                case 'string'
                    value = data(4).Children.Data;
                    value = string(strsplit(value,', '));
                case 'boolean'
                    value = int8(str2num(data(4).Children.Data));

                otherwise
                    error(['The data type is not correctly specified in the input file',...
                          ' correct the data type of variable ',name])      
            end
            unit = data(6).Children.Data;      
        case '#text'
            continue
    end
    eval(['param.',name,'=value;'])
end

if ~strcmp(variables,'all') % if the variables is not equal to all a specific set of variables is stored given as output
    for i = 1:length(variables)
        fieldName = variables(1);
        expression = ['tmp.',fieldName,'= param.',fieldName,';'];
        eval(expression);
    end
    param = tmp;
end



%% helper function
function theStruct = parseXML(filename)
% PARSEXML Convert XML file to a MATLAB structure.
try
    tree = xmlread(filename);
catch
    error('Failed to read XML file %s.',filename);
end

% Recurse over child nodes. This could run into problems
% with very deeply nested trees.
try
    theStruct = parseChildNodes(tree);
catch
    error('Unable to parse XML file %s.',filename);
end
% ----- Local function PARSECHILDNODES -----
function children = parseChildNodes(theNode)
% Recurse over node children.
children = [];
if theNode.hasChildNodes
    childNodes = theNode.getChildNodes;
    numChildNodes = childNodes.getLength;
    allocCell = cell(1, numChildNodes);
    
    children = struct(             ...
        'Name', allocCell, 'Attributes', allocCell,    ...
        'Data', allocCell, 'Children', allocCell);
    
    for count = 1:numChildNodes
        theChild = childNodes.item(count-1);
        children(count) = makeStructFromNode(theChild);
    end
end
% ----- Local function MAKESTRUCTFROMNODE -----
function nodeStruct = makeStructFromNode(theNode)
% Create structure of node info.

nodeStruct = struct(                        ...
    'Name', char(theNode.getNodeName),       ...
    'Attributes', parseAttributes(theNode),  ...
    'Data', '',                              ...
    'Children', parseChildNodes(theNode));

if any(strcmp(methods(theNode), 'getData'))
    nodeStruct.Data = char(theNode.getData);
else
    nodeStruct.Data = '';
end
% ----- Local function PARSEATTRIBUTES -----
function attributes = parseAttributes(theNode)
% Create attributes structure.

attributes = [];
if theNode.hasAttributes
    theAttributes = theNode.getAttributes;
    numAttributes = theAttributes.getLength;
    allocCell = cell(1, numAttributes);
    attributes = struct('Name', allocCell, 'Value', ...
        allocCell);
    
    for count = 1:numAttributes
        attrib = theAttributes.item(count-1);
        attributes(count).Name = char(attrib.getName);
        attributes(count).Value = char(attrib.getValue);
    end
end

