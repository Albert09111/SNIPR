function fid = writeCSVFilefromStructure(FileName, array)
%function fid = writeCSVFilefromStructure(FileName, array)
fid = fopen(FileName,'w');
for c1 = 1:size(array,1)
    for c2 = 1:size(array,2)
        if size(array{c1,c2},1) == 0
            fprintf(fid, ',');
        else
            fprintf(fid, '%s,', array{c1,c2});
        end
    end
    fprintf(fid, '\n');
end
fclose(fid);