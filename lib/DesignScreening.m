function DesignScreening(name,number_of_selected_desigs)

fprintf('### Screening the %s designs...\n',name);

energy_rank_list = dir(sprintf('DesignFiles/%s/energy*.csv',name));
seq_library = {};

for c1 = 1:length(energy_rank_list)
    file_name = energy_rank_list(c1).name;
    seq_info = loadCSVabsFileRawtoStructure(sprintf('DesignFiles/%s/%s',name,file_name),'');
    seq_info(1,:) = [];
    seq_library = [seq_library;seq_info];
end

% Define the penalty scores
% dG 
% RBS secodary structure
% overall secondary structure

sorted_score_library = seq_library;
sorted_score_library{1,8} = 'penalty_score';
penalty_score_set = [];

for c2 = 1:size(seq_library,1)
    hpin_seq = seq_library{c2,2};
    wt_seq = seq_library{c2,3};
    snp_seq = seq_library{c2,4};
    
    reactE = str2num(seq_library{c2,6});  
    mut_reacE = str2num(seq_library{c2,7});
    
    energy_penalty = abs(reactE-(-1.0))*1000;
    
    if mut_reacE < 0
        
       energy_penalty = 1000000;
        
    else
        
    end
    
    overall_def = str2num(seq_library{c2,7});
    overall_def_penalty = overall_def*500;
    
    RBS_loop_index = regexp(hpin_seq,'...AGAGGAGA......AUG');
    RBS_loop_length = length('...AGAGGAGA......AUG');
    bm_length = 17;
    
    RBS_region = hpin_seq(RBS_loop_index-3:RBS_loop_index+RBS_loop_length+bm_length);
    RBS_def = checkDefectArb1999(RBS_region,repelem('.',length(RBS_region)));
    
    RBS_penalty = RBS_def*1000;
    
    penalty_score = energy_penalty+overall_def_penalty+RBS_penalty;
    
    penalty_score_set(end+1,:) = penalty_score; 
    sorted_score_library{c2,8} = num2str(penalty_score);
    
end

% sort the sequence based on penalty score
[sorted_score, sorted_index] = sort(penalty_score_set);
sorted_design_library = sorted_score_library(sorted_index,:);

top_design = {};
top_design(1,:) = {'order','SNIPR sequence','WT sequence','Mutant sequence','Defects','SNP ReactE','WT ReactE','Penalty Score'};
top_design(2:number_of_selected_desigs+1,:) = sorted_design_library(1:number_of_selected_desigs,:);
fprintf('\n');
fprintf(' The best %d designs are as follows:\n',number_of_selected_desigs);
for c3 = 2:size(top_design,1)
    top_design{c3,1} = num2str(c3-1);
    fprintf('   Design %d:...\n',c3-1);
    fprintf('     SNIPR Seq:%s\n',top_design{c3,2});
    fprintf('     Reaction energy:%s\n',top_design{c3,6});
    fprintf('     Penalty score: %s\n',top_design{c3,8});

end
 fprintf('\n\n');
writeCSVFilefromStructure(sprintf('output/%s_design.csv',name),top_design);








