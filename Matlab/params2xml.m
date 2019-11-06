function params2xml(params,filePath)


if isstruct(params)
   % this case is for when params is an output from
   % wavemat105_sol_jeremie.m
   filePath = char(join(filePath,""));
   paramsCellVal = struct2cell(params);
   paramsCellNames = fieldnames(params);
   paramsCellType = repmat({'output'},length(paramsCellVal),1);
   paramsCellDtype = structfun(@class,params,'UniformOutput',false);
   paramsCellDtype = struct2cell(paramsCellDtype);
   idx = cell2mat(cellfun(@(dtype)strcmp(dtype,'double'),paramsCellDtype,'UniformOutput',false));
   paramsCellDtype(idx) = {'numeric'};
   paramsCellUnit = repmat({'--'},length(paramsCellVal),1);
   paramsCellExp = repmat({'--'},length(paramsCellVal),1);
   paramsCell = [paramsCellNames,paramsCellType,paramsCellVal,paramsCellUnit,...
                paramsCellDtype,paramsCellExp];
else
   % This case is for when params is generated in initialParameters.m
   paramsCell = params;
end
idx = strfind(filePath,'\');
fileName = filePath((idx(end)+1):end);

docNode = com.mathworks.xml.XMLUtils.createDocument(fileName);
nVar = size(paramsCell,1);

for iVar = 1:nVar
    a = paramsCell(iVar,:);
    
    param_node = docNode.createElement('Parameter');
    docNode.getDocumentElement.appendChild(param_node);
        
    name_node = docNode.createElement('name');
    name_text = docNode.createTextNode(a{1});
    name_node.appendChild(name_text);
    param_node.appendChild(name_node);
    
    value_node = docNode.createElement('value');
    if strcmp(a{5},"string")
        formatSpec = repmat('%s, ',1,length(a{3}));
        formatSpec = formatSpec(1:end-2);
        value_text = docNode.createTextNode(sprintf(formatSpec,a{3}));
        value_node.setAttribute('size',[num2str(1),'x',num2str(length(a{3}))]);
    elseif strcmp(a{5},"boolean")
        value_text = docNode.createTextNode(sprintf('%f',a{3}));
        value_node.setAttribute('size','1');
    elseif strcmp(a{5},"cell")
        classType = class(a{3});
        cellSize = size(a{3});
        switch classType
            case 'cell'
                a{3} = cell2mat(a{3});
                [r,c] = size(a{3});
                formatspec = [];
                
                % detect small numbers that would be written in scientific
                % notation
                if any(a{3}(:)<0.0001)
                    numFormat = '%.10f,';
                else
                    numFormat = '%f,';
                end
                    
                for i=1:r
                    tmp = repmat(numFormat,1,c);
                    tmp = [tmp(1:end-1),';\n'];
                    formatspec = [formatspec,tmp];
                end
                formatspec = formatspec(1:end-3);
                value_text = docNode.createTextNode(sprintf(formatspec,a{3}));
                value_node.setAttribute('size',[num2str(r),'x',num2str(c)]);
                cellSize_node = docNode.createElement('cell_size');
                cellSize_text = docNode.createTextNode(sprintf('%f,%f',cellSize));
                
            case 'string'
                formatSpec = repmat('%s, ',1,length(a{3}));
                formatSpec = formatSpec(1:end-2);
                value_text = docNode.createTextNode(sprintf(formatSpec,a{3}));
                value_node.setAttribute('size',[num2str(1),'x',num2str(length(a{3}))]);
        end
        
    else % the numeric case
        [r,c] = size(a{3});
        formatspec = [];
        % detect small numbers that would be written in scientific
        % notation
        if any(a{3}(:)<0.0001)
            numFormat = '%.10f,';
        else
            numFormat = '%f,';
        end
        if c == 0 && r ==0 % if you want to store an empty array
            formatspec = numFormat(1:end-1);
            value_text = docNode.createTextNode(sprintf(formatspec,[]));
            value_node.setAttribute('size','empty');
        else
            for i=1:r
                tmp = repmat(numFormat,1,c);
                tmp = [tmp(1:end-1),';\n'];
                formatspec = [formatspec,tmp];
            end
            formatspec = formatspec(1:end-3);
            value_text = docNode.createTextNode(sprintf(formatspec,a{3}'));
            value_node.setAttribute('size',[num2str(r),'x',num2str(c)]);
        end
    end
    
    value_node.appendChild(value_text);
    param_node.appendChild(value_node);
    
    unit_node = docNode.createElement('unit');
    unit_text = docNode.createTextNode(a{4});
    unit_node.appendChild(unit_text);
    param_node.appendChild(unit_node);
    
    dtype_node = docNode.createElement('data_type');
    dtype_text = docNode.createTextNode(a{5});
    dtype_node.appendChild(dtype_text);
    param_node.appendChild(dtype_node);
    
    cellExist = exist('cellSize_node');
    if cellExist == 1
        cellSize_node.appendChild(cellSize_text);
        param_node.appendChild(cellSize_node);
    end
    
    source_node = docNode.createElement('source');
    source_text = docNode.createTextNode(a{2});
    source_node.appendChild(source_text);
    param_node.appendChild(source_node);
    
    comment_node = docNode.createElement('Meaning');
    comment_text = docNode.createTextNode(a{6});
    comment_node.appendChild(comment_text);
    param_node.appendChild(comment_node);
    
    clear param_node value_node unit_node dtype_node source_node cellSize_node
end

xmlFileName = [filePath,'.xml'];
xmlwrite(xmlFileName,docNode);
