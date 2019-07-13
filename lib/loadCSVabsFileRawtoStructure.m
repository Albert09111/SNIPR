function array = loadCSVabsFileRawtoStructure(FileName, PathName)
%function array = loadCSVabsFileRawtoStructure(FileName, PathName)

%[FileName,PathName] = uigetfile('*.csv','Select an .csv file corresponding to the Wavelength; Absorbance matrix...');
TotalName = [PathName,FileName];
FileID = fopen(TotalName);

%load all lines of file
array = [];
counter = 1;
temp_string = fgetl(FileID);
while length(temp_string) ~= 0
    if isnumeric(temp_string) == 1
        break;
    end
    line1 = parse_csv_text_string(temp_string);
    for c1 = 1:length(line1)
        array{counter,c1} = line1{c1};
    end
    temp_string = fgetl(FileID);
    counter = counter + 1;
end
fclose(FileID);
end

function text_array = parse_csv_text_string(temp_string)
%function text_array = parse_csv_text_string(temp_string)

index_list = [];
for c1 = 1:length(temp_string)
    if temp_string(c1) == ','
        index_list(end+1) = c1;
    end
end

index1 = 1;
for c1 = 1:length(index_list)
    if index1 == index_list(c1)
        text_array{1,c1} = '';
    else
        text_array{1,c1} = temp_string(index1:index_list(c1)-1);
    end
    index1 = index_list(c1)+1;
end
if index1 <= length(temp_string)
    text_array{1,c1+1} = temp_string(index1:end);
end
end