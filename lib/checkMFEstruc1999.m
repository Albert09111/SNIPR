function struc_out = checkMFEstruc1999(rseq)
%function struc_out = checkMFEstruc1999(rseq)

if length(rseq) == 0
    struc_out = char(zeros(1,10)+'(');
    return 
end
filename = randseq(8);
fid = fopen([filename,'.in'],'w');
fprintf(fid,'%s',rseq);
fclose(fid);
a = 1;
counter = 1;
while a ~= 0 && counter < 5
    a = system(sprintf('mfe -material rna1999 -T 37 %s',filename));
    counter = counter + 1;
end
fid = fopen([filename,'.mfe'],'r');
if fid == -1
    struc_out = rseq;
    return;
end
temp_string = fgetl(fid);
while length(temp_string) ~= 0
    temp_string = fgetl(fid);
end
temp_string = fgetl(fid); % % line
temp_string = fgetl(fid); % length
temp_string = fgetl(fid); % deltaG
struc_out = fgetl(fid);
fclose(fid);
system(sprintf('rm %s.in %s.mfe',filename,filename));