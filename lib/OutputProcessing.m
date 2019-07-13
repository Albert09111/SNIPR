function OutputProcessing(name)

fprintf('### Extracting designs from the %s SNIPR library...\n',name);

main_dir_name = sprintf('DesignFiles/%s',name);
dir_list = dir(main_dir_name);
dir_list(1:2) = [];

not_dir_list = [];
for n1 = 1:length(dir_list)
    dir_name = [main_dir_name,'/',dir_list(n1).name];
    if ~isdir(dir_name)
        not_dir_list(end+1)=n1;    
    end
end

dir_list(not_dir_list) = [];

for c = 1:length(dir_list)
    dir_name = dir_list(c).name;
    %fprintf('processing %s...\n',dir_name);
    folder_name = sprintf('%s/%s',main_dir_name,dir_name);
    npofile  = dir(sprintf('%s/*.npo',folder_name));

    design_seq_info = {};
    FLAG = 1;
    strand_num = 0;

    for c1 = 1:length(npofile)

        cur_name = sprintf('%s/%s',folder_name,npofile(c1).name);

        seq_output = parseNUPACKfile2SeqInfo(cur_name);
        seq_list = seq_output(:,2:end)';
        if FLAG
            strand_names = seq_output(:,1)';
            design_seq_info(end+1,:) = strand_names;
            FLAG = 0;
        end
        
        strand_info = {};
        for c5 = 1:size(seq_output,1)
            strand_name = seq_output{c5,1};
            strand_seq = seq_output{c5,2};
            strand_struct = seq_output{c5,3};
            strand_def = seq_output{c5,4};
            strand_dG = seq_output{c5,5};
            
            strand_info{c5}= {strand_seq,strand_struct,strand_def,strand_dG};
            % store the strand info into a cell;
        end
        design_seq_info(end+1,:) = strand_info;
    end

    remove_list = [];
    
    for c2 = 1:length(strand_names)
        if issame(strand_names{c2},'hairpin') || issame(strand_names{c2},'rep_hairpin') || issame(strand_names{c2},'switch')||issame(strand_names{c2},'hairpin0')||issame(strand_names{c2},'hairpin_seq')
            
            for c3 = 2:size(design_seq_info,1)
                cur_hairpin_seq = design_seq_info{c3,c2}{1};
                
                indices = strfind(cur_hairpin_seq,'AGAGGAGA');
                
                %fprintf('checking stop codon for desigm...')
                
                for c4 = 1:length(indices)
                        shift_index = indices(c4)+length('AGAGGAGA')+1;
                        sub_seq = cur_hairpin_seq(shift_index:end);
                        temp = strfind(sub_seq,'AUG');
                        %assume that all AGAGGAGA are actually intended RBS
                        while ~isempty(temp)
                            if temp(1) < 4
                                temp(1) = [];
                            else
                                break;
                            end
                        end
                        
                        if temp(1) > 10
                           continue;
                        end
                        
                        sub_seq = sub_seq(temp(1):end);
                        %fprintf('%s\n',nt2aa(sub_seq));
                        if ~isempty(strfind(nt2aa(sub_seq),'*'))
                            remove_list(end+1,1) = c3;
                        end                        
                end                
            end
            design_seq_info(remove_list,:) = [];  
        end                       
    end
    
    if ~isempty(design_seq_info)
        def_num = [];
        design_seq_info_new = design_seq_info;
        design_seq_info_new(1,:) = [];
        
        for c5 = 1:size(design_seq_info_new,1)
            sum_def_num = 0;
            for c6 = 1:size(design_seq_info_new,2)
                sum_def_num = sum_def_num+str2num(design_seq_info_new{c5,c6}{3});              
            end
            def_num(end+1) = sum_def_num/size(design_seq_info_new,2);
        end
        
        [sorted_def, def_index] = sort(def_num);
        
        seq_info_output = {};
        output_list = {};
        
        seq_list = {};
        struct_list = {};
        def_list = {};
        dG_list = {};
        
        for c8 = 1:length(strand_names)
            seq_list{end+1} = sprintf('%s_seq',strand_names{c8});
            struct_list{end+1} = sprintf('%s_struct',strand_names{c8});
            def_list{end+1} = sprintf('%s_def',strand_names{c8});
            dG_list{end+1} = sprintf('%s_dG',strand_names{c8});
        end
        
        output_list = [seq_list,struct_list,def_list,dG_list,'overall def'];
        seq_info_output(end+1,:) = output_list;
        
        
        for c6 = 1:length(def_index)
            index = def_index(c6);
            strand_seq_list = {};
            strand_struct_list ={};
            strand_def_list = {};
            strand_dG_list = {};
            for c9 = 1:length(strand_names)
                strand_seq_list{end+1} = design_seq_info_new{index,c9}{1};
                strand_struct_list{end+1} = design_seq_info_new{index,c9}{2};
                strand_def_list{end+1} = design_seq_info_new{index,c9}{3};
                strand_dG_list{end+1} = design_seq_info_new{index,c9}{4};
            end
            
            strand_info_output_list = [strand_seq_list,strand_struct_list,strand_def_list,strand_dG_list,num2str(def_num(index))];
            
            seq_info_output(end+1,:) = strand_info_output_list;
        end       
    end
     
    writeCSVFilefromStructure([folder_name,'.csv'],seq_info_output);
end
