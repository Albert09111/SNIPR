function trig_seq1 = mutate_seq(type, trig_seq,SNP_index1)
%function trig_seq1 = mutation(trig_seq,SNP_index1)
if type =='RNA'
    base_seq1 = 'GUAC';
    
    if trig_seq(SNP_index1)=='A'
        base_seq1 = 'UAC';
        trig_seq1=trig_seq;
        index1 = find(base_seq1 == trig_seq1(SNP_index1));
        new_base_seq1 = base_seq1;
        new_base_seq1(index1) = []; %%mutation base pool
        trig_seq0 = trig_seq1;

        mutation_pool_index = randi(2);
        trig_seq1(SNP_index1) = new_base_seq1(mutation_pool_index);
    elseif trig_seq(SNP_index1) == 'C'
        base_seq1 = 'GAC';
        trig_seq1=trig_seq;
        index1 = find(base_seq1 == trig_seq1(SNP_index1));
        new_base_seq1 = base_seq1;
        new_base_seq1(index1) = []; %%mutation base pool
        trig_seq0 = trig_seq1;

        mutation_pool_index = randi(2);
        trig_seq1(SNP_index1) = new_base_seq1(mutation_pool_index);
    else
        trig_seq1=trig_seq;
        index1 = find(base_seq1 == trig_seq1(SNP_index1));
        new_base_seq1 = base_seq1;
        new_base_seq1(index1) = []; %%mutation base pool
        trig_seq0 = trig_seq1;

        mutation_pool_index = randi(2);
        trig_seq1(SNP_index1) = new_base_seq1(mutation_pool_index);
        
    end
      
elseif type == 'DNA'
    base_seq1 = 'GTAC';
    trig_seq1=trig_seq;
    index1 = find(base_seq1 == trig_seq1(SNP_index1));
    new_base_seq1 = base_seq1;
    new_base_seq1(index1) = []; %%mutation base pool
    trig_seq0 = trig_seq1;
    
    mutation_pool_index = randi(3);
    trig_seq1(SNP_index1) = new_base_seq1(mutation_pool_index);
else
end


end
