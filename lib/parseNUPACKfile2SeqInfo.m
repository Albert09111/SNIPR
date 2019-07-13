function seq_info = parseNUPACKfile2SeqInfo(file_name)
    seq_info = {};

    fid = fopen(file_name,'r');        
    while ~feof(fid)        
        tline = fgetl(fid);
        name_index1 = strfind(tline,'structures:');
        
        if ~isempty(name_index1)
            FLAG = 1;            
            count = 0;
            while FLAG
                tline = fgetl(fid);
                count = count+1;
                index_stop = strfind(tline,'tubes:');
                
                if ~isempty(index_stop)
                    FLAG = 0;
                else                
                        index_name = strfind(tline,'name:')+length('name:');
                        index_structure = strfind(tline,'structure:')+length('structure:');
                        index_seq = strfind(tline,'sequence:')+length('sequence:');
                        index_defect = strfind(tline,'normalized defect:')+length('normalized defect:');
                        index_dG = strfind(tline,'complex free energy[kcal/(mol K)]:')+length('complex free energy[kcal/(mol K)]:');

                        if ~isempty(index_name)
                            complex_name = tline(index_name:end);
                            %fprintf('name: %s \n',complex_name);
                            complex_name(find(isspace(complex_name))) = [];
                        end

                        
                        if ~isempty(index_seq)
                            complex_seq = tline(index_seq:end);
                            %fprintf('sequence: %s \n',complex_seq);
                            complex_seq(find(isspace(complex_seq))) = [];
                        end

                        if ~isempty(index_structure)
                            complex_structure = tline(index_structure:end);
                            %fprintf('strucure: %s \n',complex_structure);
                            complex_structure(find(isspace(complex_structure))) = [];
                        end

                        if ~isempty(index_defect)
                            complex_defect = tline(index_defect:end);
                            %fprintf('defect: %s \n',complex_defect);
                            complex_defect(find(isspace(complex_defect))) = [];
                        end

                        if ~isempty(index_dG)
                            complex_dG = tline(index_dG:end);
                            %fprintf('free energy: %s \n',complex_dG);
                            complex_dG(find(isspace(complex_dG))) = [];
                        end
                end
                
                if ~rem(count,8)
                    seq_info(end+1,:) = {complex_name,complex_seq,complex_structure,complex_defect,complex_dG};                   
                end
                
            end

        end
        
    end
    fclose(fid);
end
