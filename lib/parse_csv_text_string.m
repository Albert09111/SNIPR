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