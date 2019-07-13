function EnergyAnalysis(name)

fprintf('### Analyzing reaction energies for the %s designs...\n',name);

file=dir(sprintf('DesignFiles/%s/%s*.csv',name,name));
hpin_cancer = {'design_name','hpin_seq','WT_target','SNP_target'};
structure_info_set = loadCSVabsFileRawtoStructure(sprintf('DesignFiles/%s/design_structure_info.csv',name),'');
structure_info_set(1,:) = [];

folder_name = sprintf('DesignFiles/%s',name);
cd(folder_name);

energy_order_table = {'order','hpin_seq','wt_targ_seq','mut_targ_seq',...
                      'Overall_Defects', 'reacE','reacE_mut1'};

for c1 = 1:length(file)
 cur_name = file(c1).name;
 input_table = loadCSVabsFileRawtoStructure(cur_name,'');
 input_table(1,:) = [];
 
 % find the cooresponding structure
 design_name = cur_name(1:end-4);
 
 for c4 = 1:size(structure_info_set)
     structure_name = structure_info_set{c4,1};
     mutation_notation_index = strfind(structure_name,'>');
     if mutation_notation_index
         structure_name(mutation_notation_index) = 'm';
     else
     end
     
     if ~isempty(strfind(structure_name,design_name))
         structure_index = c4;
     else
     end
 end
 
 structure_name = structure_info_set{structure_index,1};
 structure_correct_on = DUnotation2DotParens(structure_info_set{structure_index,2});
 on_break_index = strfind(structure_correct_on,'+');
 structure_complex_correct_on =[structure_correct_on(1:on_break_index-1),'.....',structure_correct_on(on_break_index+1:end)];
 
 structure_mut_on = DUnotation2DotParens(structure_info_set{structure_index,4});
 mut_on_index = strfind(structure_mut_on,'+');
 structure_complex_mutON =[structure_mut_on(1:mut_on_index-1),'.....',structure_mut_on(mut_on_index+1:end)];
 
 structure_correct_off = DUnotation2DotParens(structure_info_set{structure_index,3});
 off_index = strfind(structure_correct_off,'+');
 structure_complex_correct_OFF = [structure_correct_off(1:off_index-1),'.....',structure_correct_off(off_index+1:end)];

 structure_mut_off = DUnotation2DotParens(structure_info_set{structure_index,5});
 off_index = strfind(structure_mut_off,'+');
 structure_complex_mut_OFF = [structure_mut_off(1:off_index-1),'.....',structure_mut_off(off_index+1:end)];
 
 loop_seq = 'AAAAAAAA';
 deltaG_table = [];
 ddg1_table = [];

 %fprintf('start_computing_design_%s \n',cur_name);
 for c2 = 1:size(input_table,1)
     
         %fprintf('Computing design %d of %d...\n',c2,size(input_table,1));
         hpin = input_table{c2,3};
         snp_targ = input_table{c2,4};
         wt_targ = input_table{c2,5};

         wt_targ0 = wt_targ;
         snp_targ0 = snp_targ;
         hpin0=hpin(4:end);
     
        long_seq = [snp_targ0,loop_seq,hpin0];
        long_seq_mut = [wt_targ0,loop_seq,hpin0];

        dg1 = checkStructureFreeEnergyRNA1999(long_seq,structure_complex_correct_OFF);
        dg2 = checkStructureFreeEnergyRNA1999(long_seq,structure_complex_correct_on);
        dg3 = checkStructureFreeEnergyRNA1999(long_seq_mut,structure_complex_mut_OFF);
        dg4 = checkStructureFreeEnergyRNA1999(long_seq_mut,structure_complex_mutON);
        
        ddg1 = dg2-dg1;
        ddg2 = dg4-dg3;
        deltaG_table(end+1,:) = [dg1, dg2, dg3, dg4, ddg1, ddg2];
        ddg1_table(end+1,:) = ddg1;
 end
 
% put the design in order based on the ddg1.
[sorted_ddg1, ddg1_indices] = sort(ddg1_table, 'descend');
ddg1_order = [ddg1_indices,sorted_ddg1];

for c4 = 1:size(input_table,1)
    index = ddg1_indices(c4,1);
    rectE = sorted_ddg1(c4,1);
    rectE_mut = deltaG_table(index,6);
    
    hairpin_seq = input_table{index,3};
    WT_target = input_table{index,5};
    SNP_target = input_table{index,4};
    
    
    energy_order_table(end+1,:) = {num2str(index),input_table{index,3},WT_target(4:end),SNP_target(4:end),...
                                   input_table{index,21},num2str(rectE),num2str(rectE_mut)};
end

end

%eval(sprintf('delete %s_R*.csv',name))

filename = sprintf('energy_%s.csv',name);
writeCSVFilefromStructure(filename,energy_order_table);
cd ../..;