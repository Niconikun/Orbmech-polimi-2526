function [M] = readHorizonsCar(filename)

    % 1. Read file
    txt = fileread(filename);
    
    % 2. Extract the text between $$SOE and $$EOE
    dataSection = regexp(txt, '\$\$SOE(.*?)\$\$EOE', 'tokens', 'dotall');
    if isempty(dataSection)
        error('The section was not found');
    end
    dataSection = dataSection{1}{1};
    
    % 3. Pattern for the extraction of data
    pattern = ['X\s*=\s*([-\d.E+]+)\s*Y\s*=\s*([-\d.E+]+)\s*Z\s*=\s*([-\d.E+]+)\s*' ...
               'VX\s*=\s*([-\d.E+]+)\s*VY\s*=\s*([-\d.E+]+)\s*VZ\s*=\s*([-\d.E+]+)\s*' ...
               'LT\s*=\s*([-\d.E+]+)\s*RG\s*=\s*([-\d.E+]+)\s*RR\s*=\s*([-\d.E+]+)'];
    
    % 4. Extraction of data
    tokens = regexp(dataSection, pattern, 'tokens');
    if isempty(tokens)
        error('The section was not found');
    end
    
    % 5. Conversion in matrix
    M = cellfun(@(x) str2double(x), tokens, 'UniformOutput', false);
    M = vertcat(M{:});

end