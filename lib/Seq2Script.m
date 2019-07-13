function Seq2Script(design_name,WT_seq,SNP_seq,protein_seq,SNP_TARGET)

fprintf('Starting design for %s SNIPRs...\n',design_name);
fprintf('### Generating NUPACK design scripts...\n')

toehold_rev = [4];

design_structure_info = {'name','struct_complex_correct_on','struct_complex_correct_off','struct_complex_mut_on','struct_complex_mut_off'};

for c2 = 1:size(toehold_rev,2)
    toehold_rev_length = toehold_rev(c2);
    toehold_fwd = [toehold_rev_length-1,toehold_rev_length,toehold_rev_length+1];
    for c3 = 1:size(toehold_fwd,2)
        toehold_fwd_length = toehold_fwd(c3);
        sample_name = design_name;
        cur_targ_name = sprintf('%s_R%dF%d',sample_name,toehold_rev_length,toehold_fwd_length);

        %Generate the NUPACK design file:
        folder_name = sprintf('DesignFiles/%s',sample_name);

        mkdir(folder_name);
        cd(folder_name);
        number_of_trials = 10;
        material_type = 'rna';
        design_temperature = 37;

        if SNP_TARGET
            wt_target_seq = dna2rna(WT_seq);
            mut_target_seq = dna2rna(SNP_seq);
        else
            wt_target_seq = dna2rna(SNP_seq);
            mut_target_seq = dna2rna(WT_seq);
        end 
        %wt_target_seq = mutate_seq('RNA',mut_target_seq,31);

        struct1 = localalign(wt_target_seq, mut_target_seq, 'alpha', 'nt');
        sequence_align1 = struct1.Alignment{1};

        snp_target_leng = length(mut_target_seq);
        wt_target_leng = length(wt_target_seq);
        SNP_position = find(sequence_align1(2,:) == ' ');

        SNP_length = length(SNP_position);
        
        if ~isempty(find(sequence_align1(1,:)=='-'))

            mutation_type = 'insertion';

        elseif ~isempty(find(sequence_align1(3,:)=='-'))

            mutation_type = 'deletion';
            
        else 
            
            mutation_type = 'substitution';

        end

        %%%%set the design parameter
        bm_length = 21-toehold_rev_length;

        % specify the mismatch position 
        mismatch_position = 10;

        if mismatch_position > SNP_position(1) || mismatch_position == SNP_position(1)

            mismatch_position = mismatch_position - 2;

        else

        end
        
        SNP_position = SNP_position(1);
        bm_start = SNP_position-mismatch_position;
        toe_fwd_start = bm_start+bm_length;
        bulge_start = toe_fwd_start+toehold_fwd_length;
        docking_start = bulge_start+10;

        over_hang = mut_target_seq(1:bm_start-1);
        bm_reigon = mut_target_seq(bm_start:toe_fwd_start-1);% mutant target can open the hairpin
        toehold_fwd_seq = mut_target_seq(toe_fwd_start:toe_fwd_start+toehold_fwd_length-1);
        toehold_rev_seq = mut_target_seq(toe_fwd_start:toe_fwd_start+toehold_rev_length-1);
        bulge_region = mut_target_seq(bulge_start:docking_start-1);
        docking_region = mut_target_seq(docking_start:end);

        overhang_length = length(over_hang);

        docking_length = length(docking_region);
        
        prehpin_length = docking_length+10+toehold_fwd_length;
             
        if toehold_fwd_seq(end)=='G'&& toehold_fwd_seq(end-1) ~= 'G'
           toehold_rev_seq = [toehold_rev_seq(1:end-3),toehold_fwd_seq(end),toehold_fwd_seq(end-2),toehold_fwd_seq(end-1)];
        else
        end

        before_mutation_length = SNP_position-bm_start;
        after_mutation_length = toe_fwd_start-SNP_position-1+toehold_fwd_length;

        struct_target1 = sprintf('U3 U%d',wt_target_leng);
        struct_target2 = sprintf('U3 U%d',snp_target_leng);
        struct_hairpin = sprintf('U3 U%d D21( U20 ) U50',prehpin_length); % with AUG at loop site, 20 nt loop
        struct_complex_correct_on = sprintf('U3 U%d D%d( U10 D%d( + U3 ) U10 ) U%d',overhang_length,bm_length+toehold_fwd_length,docking_length,91+toehold_rev_length); % Open state 
        struct_complex_correct_off = sprintf('U3 U%d U%d U10 D%d( + U3 ) U%d D21( U20) U50',overhang_length,bm_length+toehold_fwd_length,docking_length,toehold_fwd_length+10); % Close state
        %check stop codon in stem region.
        

        
        stem_seq = [toehold_rev_seq,bm_reigon];
        aa_stem = nt2aa(stem_seq);

        shift_length = 1;
        while contains(aa_stem,'*')   

            fprintf('design:%s \n warning: shifted %d \n',cur_targ_name,shift_length);
            cur_targ_name = sprintf('%s_shift%d',cur_targ_name,shift_length); % move SNP position one nucleotide ahead.

            bm_length = 21-toehold_rev_length;

            bm_start = SNP_position-mismatch_position-shift_length;
            toe_fwd_start = bm_start+bm_length;
            bulge_start = toe_fwd_start+toehold_fwd_length;
            docking_start = bulge_start+10;

            over_hang = mut_target_seq(1:bm_start-1);
            bm_reigon = mut_target_seq(bm_start:toe_fwd_start-1);% mutant target can open the hairpin
            toehold_fwd_seq = mut_target_seq(toe_fwd_start:toe_fwd_start+toehold_fwd_length-1);
            toehold_rev_seq = mut_target_seq(toe_fwd_start:toe_fwd_start+toehold_rev_length-1);
            bulge_region = mut_target_seq(bulge_start:docking_start-1);
            docking_region = mut_target_seq(docking_start:end);

            docking_length = length(docking_region);
            overhang_length = length(over_hang);
            prehpin_length = docking_length+10+toehold_fwd_length;

            stem_seq = [toehold_rev_seq,bm_reigon];
            aa_stem = nt2aa(stem_seq);

            if toehold_rev_seq(end)=='G'&& toehold_rev_seq(end-1) ~= 'G'
               toehold_rev_seq = [toehold_rev_seq(1:end-3),toehold_rev_seq(end),toehold_rev_seq(end-2),toehold_rev_seq(end-1)];
            else 
            end
            before_mutation_length = SNP_position-bm_start;
            after_mutation_length = toe_fwd_start-SNP_position-1+toehold_fwd_length;

            struct_target1 = sprintf('U3 U%d',wt_target_leng);
            struct_target2 = sprintf('U3 U%d',snp_target_leng);
            struct_hairpin = sprintf('U3 U%d D21( U20 ) U50',prehpin_length); % with AUG at loop site, 20 nt loop
            struct_complex = sprintf('U3 U%d D%d( U10 D%d( + U3 ) U10 ) U%d',overhang_length,bm_length+toehold_fwd_length,docking_length,91+toehold_rev_length); % Open state
            struct_complex_mut_on = sprintf('U3 U%d D%d( U1 D%d( U10 D%d( + U3 ) U10 ) U1 ) U%d',overhang_length,before_mutation_length,after_mutation_length,docking_length,91+toehold_rev_length); % Open state 
            struct_complex_correct_off = sprintf('U3 U%d U%d U10 D%d( + U3 ) U%d D21( U20) U50',overhang_length,bm_length+toehold_fwd_length,docking_length,toehold_fwd_length+10); % Close state
            %check stop codon in stem region.
            shift_length = shift_length+1;  
        end

        if contains(mutation_type,'insertion')
            
            struct_complex_mut_on = sprintf('U3 U%d D%d(  D%d( U10 D%d( + U3 ) U10 ) U1 ) U%d',overhang_length,before_mutation_length,after_mutation_length,docking_length,91+toehold_rev_length); % Open state 
            struct_complex_mut_off = sprintf('U3 U%d U%d U10 D%d( + U3 ) U%d D21( U20) U50',overhang_length,bm_length+toehold_fwd_length-1,docking_length,toehold_fwd_length+10); % Close state

        elseif contains(mutation_type,'deletion')
            
            struct_complex_mut_on = sprintf('U3 U%d D%d( U%d D1( D%d( U10 D%d( + U3 ) U10 ) ) ) U%d',overhang_length,before_mutation_length,SNP_length,after_mutation_length,docking_length,91+toehold_rev_length); % Open state 
            struct_complex_mut_off = sprintf('U3 U%d U%d U10 D%d( + U3 ) U%d D21( U20) U50',overhang_length,bm_length+toehold_fwd_length+SNP_length,docking_length,toehold_fwd_length+10); % Close state
        
        elseif contains(mutation_type,'substitution')
            
            struct_complex_mut_on = sprintf('U3 U%d D%d( U1 D%d( U10 D%d( + U3 ) U10 ) U1 ) U%d',overhang_length,before_mutation_length,after_mutation_length,docking_length,91+toehold_rev_length); % Open state 
            struct_complex_mut_off = sprintf('U3 U%d U%d U10 D%d( + U3 ) U%d D21( U20) U50',overhang_length,bm_length+toehold_fwd_length,docking_length,toehold_fwd_length+10); % Close state
        else  
        end     
       
       %fprintf('mutation on structure length = %d \n',length(DUnotation2DotParens(struct_complex_mut_on)));
       %fprintf('mutation off structure length = %d \n',length(DUnotation2DotParens(struct_complex_mut_off)));
        
       domain_hpin_loop = 'NNNAGAGGAGANNNNNNAUG';
       domain_hpin_bulge = repelem('N',10);
       domain_post_hpin = [repelem('N',21),protein_seq];

       hpin_seq = ['GGG',docking_region,bulge_region,toehold_fwd_seq,bm_reigon,toehold_rev_seq,domain_hpin_loop,...
                    d_revcomp(toehold_rev_seq),d_revcomp(bm_reigon),domain_post_hpin];     
       hpin_parthe = DUnotation2DotParens(struct_hairpin);

       length(hpin_seq);
       length(hpin_parthe);

       file_name = sprintf('%s.txt',cur_targ_name);
       fid = fopen(file_name,'w');

       fprintf(fid,'material = %s\n',material_type);
       fprintf(fid,'temperature = %d \n\n',design_temperature);

       fprintf(fid,'domain post_hpin= %s\n',domain_post_hpin);
       fprintf(fid,'domain WT_target_seq = %s\n',wt_target_seq);
       fprintf(fid,'domain SNP_target_seq = %s\n',mut_target_seq);
       fprintf(fid,'domain preGGG = GGG \n');
       fprintf(fid,'domain hpin_loop = %s\n',domain_hpin_loop);
       fprintf(fid,'domain hpin_bulge = %s\n',domain_hpin_bulge);
       fprintf(fid,'domain docking_region = %s\n',docking_region);
       fprintf(fid,'domain toehold1 = %s\n',toehold_fwd_seq);
       fprintf(fid,'domain toehold2 = %s\n',toehold_rev_seq);
       fprintf(fid,'domain bm = %s\n\n',bm_reigon);

       fprintf(fid,'strand target1 = preGGG WT_target_seq \n'); % target1 is wt
       fprintf(fid,'strand target2 = preGGG SNP_target_seq \n');% target2 is snp           
       fprintf(fid,'strand hairpin0 = preGGG docking_region* hpin_bulge toehold1* bm* toehold2* hpin_loop toehold2 bm post_hpin \n\n');

       fprintf(fid,'complex target_wt = target1 \n');
       fprintf(fid,'complex target_snp = target2 \n');           
       fprintf(fid,'complex hairpin = hairpin0 \n');
       fprintf(fid,'complex assemble_on = target2 hairpin0 \n');
       fprintf(fid,'complex assemble_off = target1 hairpin0 \n\n');

       fprintf(fid,'target_wt.structure = %s\n',struct_target1);
       fprintf(fid,'target_snp.structure = %s\n',struct_target2);
       fprintf(fid,'hairpin.structure = %s\n',struct_hairpin);
       fprintf(fid,'assemble_on.structure = %s\n',struct_complex_correct_on);%OPEN state
       fprintf(fid,'assemble_off.structure = %s\n\n',struct_complex_mut_off);%OFF state

       fprintf(fid,'stop[%%] = 50 \n');
       fprintf(fid,'prevent = AAAAAAAA, CCCCCCCC, GGGGGGGG, UUUUUUUU\n');

       fclose(fid);
       cd ../..

       design_structure_info(end+1,:) = {cur_targ_name,struct_complex_correct_on,struct_complex_correct_off,struct_complex_mut_on,struct_complex_mut_off};
    end
    fclose('all');
    writeCSVFilefromStructure(sprintf('%s/design_structure_info.csv',folder_name),design_structure_info);
end