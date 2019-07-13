function dna_comp_string = d_revcomp(dna_string)
%function dna_comp_string = d_revcomp(dna_string)

if (length(findstr(dna_string,'u')) + length(findstr(dna_string,'U'))) > 0
    dna_string = rna2dna(dna_string);
end
dna_string = fliplr(dna_string);
for c1 = 1:length(dna_string)
    if dna_string(c1) == 'G'
        dna_string(c1) = 'C';
    elseif dna_string(c1) == 'C'
        dna_string(c1) = 'G';
    elseif dna_string(c1) == 'A'
        dna_string(c1) = 'T';
    elseif dna_string(c1) == 'T'
        dna_string(c1) = 'A';
    end
end
dna_comp_string = dna_string;