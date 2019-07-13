function out_seq = loadNupackDefect(FileName,index)
%function out_seq = loadNupackDefect(FileName,index)
%index = 1: takes first number (non-normalized defect)
%index = 2: takes second nuber (normalized defect)

out_seq = 1;
fid = fopen(FileName,'r');
if fid == -1
    return;
end
temp_string = fgetl(fid);
out_seq = 1;
while ischar(temp_string) && isempty(temp_string) ~= 1
    if temp_string(1) ~= '%'
        out_seq = temp_string;
        index = index - 1;
        if index == 0
            out_seq = str2num(out_seq);
            fclose(fid);
            return;
        end
    end
    temp_string = fgetl(fid);
end
fclose(fid);
