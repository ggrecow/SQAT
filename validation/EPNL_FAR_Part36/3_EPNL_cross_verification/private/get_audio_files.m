function [NAF, AURAFONE, FLAURA] = get_audio_files(main_path, nFiles, tool_tag)

for j = 1:length(tool_tag)
    
    tag_files = { ['VR-app-case5-' tool_tag{j} '.wav'], ...
                      ['VR-dep-case6-' tool_tag{j} '.wav'], ...
                      ['V2-app-case7-' tool_tag{j} '.wav'], ...
                      ['V2-dep-case8-' tool_tag{j} '.wav'], ...
                      ['V2-dep-case8a-' tool_tag{j} '.wav'] }; % name of the input signals
    
    for i = 1:nFiles
        
        switch tool_tag{j}
            case 'NAF'
                [ NAF{i}.audio, NAF{i}.fs] = audioread( [ main_path tag_files{i} ] );
            case 'AURAFONE'
                [ AURAFONE{i}.audio, AURAFONE{i}.fs ] = audioread( [ main_path tag_files{i} ] );
            case 'FLAURA'
                [ FLAURA{i}.audio, FLAURA{i}.fs ] = audioread( [ main_path tag_files{i} ] );
        end
        
    end
end

end