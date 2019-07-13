function energy = loadNupackEnergy(FileName)
%function energy = loadNupackEnergy(FileName)

out_seq = 1;
fid = fopen(FileName,'r');
if fid == -1
    return;
end
temp_string = fgetl(fid);
energy = [];
while ischar(temp_string)
    if temp_string(1) ~= '%'
        energy = str2num(temp_string);
        fclose(fid);
        return;
    end
    temp_string = fgetl(fid);
end
fclose(fid);
