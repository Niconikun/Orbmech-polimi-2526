function [ephem_data, summary] = parse_horizons_ephemeris(filename)
% PARSE_HORIZONS_EPHEMERIS Parse NASA Horizons ephemerides data file
%
% INPUT:
%   filename - Path to the Horizons ephemerides text file
%
% OUTPUT:
%   ephem_data - Structure array with ephemerides data
%   summary - Structure with mission summary information

    % Initialize outputs
    ephem_data = [];
    summary = struct();
    
    % Read the entire file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    
    file_content = textscan(fid, '%s', 'Delimiter', '\n', 'WhiteSpace', '');
    fclose(fid);
    lines = file_content{1};
    
    % Parse header information
    in_header = true;
    in_data = false;
    data_lines = {};
    
    for i = 1:length(lines)
        line = strtrim(lines{i});
        
        % Skip empty lines
        if isempty(line)
            continue;
        end
        
        % Check for data section markers
        if strcmp(line, '$$SOE')
            in_data = true;
            in_header = false;
            continue;
        elseif strcmp(line, '$$EOE')
            in_data = false;
            break;
        end
        
        % Parse header information
        if in_header
            % Extract mission parameters with more flexible patterns
            if contains(line, 'Apogee altitude')
                tokens = regexp(line, 'Apogee altitude\s*=\s*~?\s*([\d,]+)\s*km', 'tokens');
                if ~isempty(tokens)
                    summary.apogee_alt = str2double(strrep(tokens{1}{1}, ',', ''));
                end
            end
            
            if contains(line, 'Perigee altitude')
                tokens = regexp(line, 'Perigee altitude\s*=\s*~?\s*([\d,]+)\s*km', 'tokens');
                if ~isempty(tokens)
                    summary.perigee_alt = str2double(strrep(tokens{1}{1}, ',', ''));
                end
            end
            
            if contains(line, 'Inclination')
                tokens = regexp(line, 'Inclination\s*=\s*([\d.]+)\s*degrees', 'tokens');
                if ~isempty(tokens)
                    summary.inclination = str2double(tokens{1}{1});
                end
            end
            
            % More flexible orbital period pattern
            if contains(line, 'Orbital period')
                % Try multiple patterns
                tokens = regexp(line, 'Orbital period\s*=\s*~?\s*([\d.]+)\s*hours', 'tokens');
                if isempty(tokens)
                    tokens = regexp(line, 'Orbital period\s*=\s*~?\s*([\d.]+)', 'tokens');
                end
                if ~isempty(tokens)
                    summary.period_hours = str2double(tokens{1}{1});
                end
            end
        end
        
        % Collect data lines
        if in_data
            data_lines{end+1} = line;
        end
    end
    
    % Parse ephemerides data
    if ~isempty(data_lines)
        ephem_data = parse_ephemerides_data(data_lines);
    end
    
    % Calculate additional orbital parameters
    if isfield(summary, 'apogee_alt') && isfield(summary, 'perigee_alt')
        R_E = 6378.137; % Earth radius in km
        summary.a = (summary.apogee_alt + summary.perigee_alt + 2*R_E) / 2; % semi-major axis
        summary.e = (summary.apogee_alt - summary.perigee_alt) / (summary.apogee_alt + summary.perigee_alt + 2*R_E); % eccentricity
    end
    
    % If period wasn't found in header, calculate it from semi-major axis
    if ~isfield(summary, 'period_hours') && isfield(summary, 'a')
        mu_E = 398600.4418; % Earth's gravitational parameter [km^3/s^2]
        period_seconds = 2*pi*sqrt(summary.a^3/mu_E);
        summary.period_hours = period_seconds / 3600;
    end
end

function data = parse_ephemerides_data(data_lines)
% Parse the actual ephemerides data lines
    
    data = struct('time', [], 'ra', [], 'dec', [], 'range', [], 'range_rate', [], ...
                  'solar_elongation', [], 'phase_angle', [], 'ra_deg', [], 'dec_deg', []);
    
    for i = 1:length(data_lines)
        line = data_lines{i};
        
        % Skip non-data lines
        if isempty(line) || contains(line, 'Date__') || contains(line, '$$')
            continue;
        end
        
        % Split by commas
        tokens = strsplit(line, ',');
        if length(tokens) < 8
            continue;
        end
        
        % Parse time
        time_str = strtrim(tokens{1});
        try
            data(i).time = datetime(time_str, 'InputFormat', 'yyyy-MMM-dd HH:mm');
        catch
            try
                data(i).time = datetime(time_str, 'InputFormat', 'yyyy-MM-dd HH:mm');
            catch
                data(i).time = NaT;
            end
        end
        
        % Parse Right Ascension (HH MM SS.ss)
        ra_str = strtrim(tokens{4});
        data(i).ra = ra_str;
        data(i).ra_deg = hms2degrees(ra_str);
        
        % Parse Declination (±DD MM SS.s)
        dec_str = strtrim(tokens{5});
        data(i).dec = dec_str;
        data(i).dec_deg = dms2degrees(dec_str);
        
        % Parse range (km)
        range_str = strtrim(tokens{8});
        data(i).range = str2double(range_str);
        
        % Parse range rate (km/s)
        range_rate_str = strtrim(tokens{9});
        data(i).range_rate = str2double(range_rate_str);
        
        % Parse solar elongation and phase angle
        if length(tokens) >= 11
            solar_elong_str = strtrim(tokens{10});
            if contains(solar_elong_str, '/')
                solar_elong_str = strsplit(solar_elong_str, '/');
                data(i).solar_elongation = str2double(solar_elong_str{1});
            else
                data(i).solar_elongation = str2double(solar_elong_str);
            end
            
            phase_angle_str = strtrim(tokens{11});
            data(i).phase_angle = str2double(phase_angle_str);
        end
    end
end

function deg = hms2degrees(hms_str)
% Convert HH MM SS.ss to degrees (1 hour = 15 degrees)
    tokens = regexp(hms_str, '(\d+)\s+(\d+)\s+([\d.]+)', 'tokens');
    if ~isempty(tokens)
        hours = str2double(tokens{1}{1});
        minutes = str2double(tokens{1}{2});
        seconds = str2double(tokens{1}{3});
        deg = 15 * (hours + minutes/60 + seconds/3600); % 1 hour = 15 degrees
    else
        deg = NaN;
    end
end

function deg = dms2degrees(dms_str)
% Convert ±DD MM SS.s to decimal degrees
    tokens = regexp(dms_str, '([+-]?)(\d+)\s+(\d+)\s+([\d.]+)', 'tokens');
    if ~isempty(tokens)
        sign_char = tokens{1}{1};
        degrees = str2double(tokens{1}{2});
        minutes = str2double(tokens{1}{3});
        seconds = str2double(tokens{1}{4});
        
        deg = degrees + minutes/60 + seconds/3600;
        if strcmp(sign_char, '-')
            deg = -deg;
        end
    else
        deg = NaN;
    end
end