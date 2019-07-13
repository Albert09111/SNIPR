function deltaG = checkStructureFreeEnergyRNA1999(rseq,rstruc)
%function deltaG = checkStructureFreeEnergyRNA1999(rseq,rstruc)
%can handle multi strand constructs
%may have errors if no base pairing between strands

if length(rseq) == 0
    fprintf('No sequences specified in checkStructureFreeEnergyRNA.\n');
    deltaG = -1e20;
    return 
end
filename = randseq(16);
fid = fopen([filename,'.in'],'w');
if iscell(rseq)
    fprintf(fid,'%d\n',length(rseq));
    for c1 = 1:length(rseq)
        fprintf(fid,'%s\n',rseq{c1});
    end
    for c1 = 1:length(rseq)
        fprintf(fid,'%d ',c1);
    end
    fprintf(fid,'\n%s',rstruc);
    fclose(fid);
    a = 1;
    counter = 1;
    while a ~= 0 && counter < 5
        a = system(sprintf('energy -material rna1999 -T 37 -multi %s > %s.out',filename,filename));
        counter = counter + 1;
    end    
else
    fprintf(fid,'%s\n%s',rseq,rstruc);
    fclose(fid);
    a = 1;
    counter = 1;
    while a ~= 0 && counter < 5
        a = system(sprintf('energy -material rna1999 -T 37 %s > %s.out',filename,filename));
        counter = counter + 1;
    end
end
%disp([filename,'.out']);
deltaG = loadNupackEnergy([filename,'.out']);

system(sprintf('rm %s.in %s.out',filename,filename));