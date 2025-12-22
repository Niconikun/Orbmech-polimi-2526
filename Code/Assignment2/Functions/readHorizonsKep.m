function [M] = readHorizonsKep(filename)

    % 1. Read file
    txt = fileread(filename);
    
    % 2. Extract the text between $$SOE and $$EOE
    dataSection = regexp(txt, '\$\$SOE(.*?)\$\$EOE', 'tokens', 'dotall');
    if isempty(dataSection)
        error('The section was not found');
    end
    dataSection = dataSection{1}{1};
    
    % 3. Pattern for the extraction of data
    pattern = ['EC=\s*([\dE\+\-\.]+).*?' ...
               'QR=\s*([\dE\+\-\.]+).*?' ...
               'IN=\s*([\dE\+\-\.]+).*?' ...
               'OM=\s*([\dE\+\-\.]+).*?' ...
               'W\s*=\s*([\dE\+\-\.]+).*?' ...
               'Tp=\s*([\dE\+\-\.]+).*?' ...
               'N\s*=\s*([\dE\+\-\.]+).*?' ...
               'MA=\s*([\dE\+\-\.]+).*?' ...
               'TA=\s*([\dE\+\-\.]+).*?' ...
               'A\s*=\s*([\dE\+\-\.]+).*?' ...
               'AD=\s*([\dE\+\-\.]+).*?' ...
               'PR=\s*([\dE\+\-\.]+)'];
    
    % 4. Extraction of data
    tokens = regexp(dataSection, pattern, 'tokens');
    if isempty(tokens)
        error('The section was not found');
    end
    
    % 5. Conversion in matrix
    M = cellfun(@(x) str2double(x), tokens, 'UniformOutput', false);
    M = vertcat(M{:});

end