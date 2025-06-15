%% Align comments at a specified column
function function_align_comments(file_name, comment_column)
    % Default column for comments
    if nargin < 2
        comment_column = 40; 
    end
    
    % Try to open the file for reading
    file_identifier = fopen(file_name, 'r');
    if file_identifier == -1
        error('Could not open the file: %s. Please check the file path.', file_name);
    end
    
    % Read the file line by line
    lines = {};
    while ~feof(file_identifier) % feof: Test for end of file
        line = fgetl(file_identifier); % fgetl: Read line from file, removing newline characters
        lines{end+1} = line; %#ok<AGROW>
    end
    fclose(file_identifier);
    
    % Process each line to align comments, but skip headings
    for i = 1:length(lines)
        line = lines{i};
        comment_index = strfind(line, '%'); % strfind: Find strings within other strings
        
        % Skip lines that are full comments or section headers (%%)
        if isempty(comment_index) || startsWith(strtrim(line), '%') || startsWith(strtrim(line), '%%') 
            % startsWith: Determine if string starts with substring in Requirements Table block. strtrim: Remove leading and trailing whitespace from strings
            % Skip heading or section header without aligning
            continue;
        end
        
        % Align comments only if there's code before the comment
        code_part = strtrim(line(1:comment_index(1)-1)); % Code before comment
        comment_part = strtrim(line(comment_index(1):end)); % Comment part
        if length(code_part) + 1 > comment_column
            newLine = [code_part ' ' comment_part]; % Keep as is if too long
        else
            spaces = repmat(' ', 1, comment_column - length(code_part)); % Create spaces. repmat: Repeat copies of array
            newLine = [code_part spaces comment_part]; % Align comment
        end
        lines{i} = newLine; %#ok<AGROW>
    end
    
    % Write the updated file back
    file_identifier = fopen(file_name, 'w');
    if file_identifier == -1
        error('Could not write to the file: %s.', file_name); % "%s" is a format specifier used in functions like sprintf, fprintf, and input to represent a string. It tells MATLAB to insert or display a string at that position in the output.
    end
    for i = 1:length(lines)
        fprintf(file_identifier, '%s\n', lines{i});
    end
    fclose(file_identifier);
end
