function defect_level = checkDefectArb1999(rseq,rstruc)
%function defect_level = checkDefectArb1999(rseq,rstruc)

filename = randseq(16);
fid = fopen([filename,'.in'],'w');
fprintf(fid,'%s\n%s',rseq,rstruc);
fclose(fid);
a = 1;
counter = 1;
while a ~= 0 && counter < 50
    a = system(sprintf('complexdefect -T 37 -material rna1999 %s > %s.out',filename,filename));
    counter = counter + 1;
end
defect_level = loadNupackDefect([filename,'.out'],2);
system(sprintf('rm %s.in %s.out',filename,filename));