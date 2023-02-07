function [data] = cea_detonation_read(file_name, data)
% OUTREAD: Reads CEA rocket output files.
%
% Input
%   file_name: The name of the '.out' data file that CEA generates
%   data: 'cea_detonation_read' data file to append the results of processing
%       the contents of 'file_name' to. Only works if the reactants part of
%       the input file is the same between runs. Has a datatype of
%       containers.Map. This is an OPTIONAL variable.
%
% Output
%   data: Contains the contents of the output file
%   data(type): The type of CEA analysis; 'det'
%   data(type)('fuel*'): Various fuel properties including the fuel name,
%       temperature, weight percentage, and energy (for various *)
%   data(type)('ox*'): Various oxidizer properties including the oxidizer
%       name, temperature, weight percentage, and energy (for various *)
%   data(type)(property): All of the properties output by CEA including
%       performance, thermal, and species (for equilibrium only) in a 3D
%       matrix. Unburned gas will have a suffix of "1", while burned gas
%       will not have any suffix. Frozen reactions will have the suffix
%       "_fr", while equilibrium reactions will have the suffix "_eq".
%
% Notes
%   All properties are in SI units (kg, m, s, K, Pa, J).
%   
%   Does not verify that the same reactants are used when appending to a
%       given data file. Ask the current maintainer to implement this
%       feature if you need it.
%   Does not allow lower CEA trace values. Ask the current maintainer to
%       implement this feature if you need it.
%   Reorders the pressure, mixture, and ratio inputs such that they are
%       sorted in ascending order.
%
% Author: Thomas Beimrohr
% Maintainer: Thomas Beimrohr
% Contact: (812) 946-2715

% Loop through every line of the file
fileID = fopen(file_name, 'r');
if fileID == -1
    error('Can''t open file %s', file_name)
end
% initialize variables from CEA
map = containers.Map;
map('P1') = [];
map('T1') = [];
map('H1') = [];
map('M1') = [];
map('GAM1') = [];
map('SON1') = [];
map('P') = [];
map('T') = [];
map('RHO') = [];
map('H') = [];
map('U') = [];
map('G') = [];
map('S') = [];
map('M') = [];
map('MW') = [];
map('dLV_dLT_p') = [];
map('dLV_dLP_t') = [];
map('CP') = [];
map('GAMs') = [];
map('SON') = [];
map('VISC') = [];
map('CP_eq') = [];
map('COND_eq') = [];
map('PR_eq') = [];
map('CP_fr') = [];
map('COND_fr') = [];
map('PR_fr') = [];
map('P_P1') = [];
map('T_T1') = [];
map('M_M1') = [];
map('RHO_RHO1') = [];
map('DET_MACH') = [];
map('DET_VEL') = [];

line = strtrim(fgetl(fileID));
found_a_section = false;
ratio_type = '';
while ischar(line)
    line = strtrim(line);

    if startsWith2(line, 'DETONATION PROPERTIES')
        % Instantiate the map and fuel/ox arrays
        found_a_section = true;
        fuel = {};
        fuel_wt = [];
        fuel_energy = [];
        fuel_T = [];
        ox = {};
        ox_wt = [];
        ox_energy = [];
        ox_T = [];
        map('type') = 'det';
        
    elseif startsWith2(line, '(SEE NOTE)')
        line2 = fgetl(fileID); 
        if startsWith2(line2, ' FUEL')
            % Get the fuels
            line_next = line2;
            while startsWith2(line_next, ' FUEL')
            fuel_index = size(fuel, 1) + 1;
            fuel_split = strsplit(strtrim(line_next(length(' FUEL '):end)));
            fuel{fuel_index} = fuel_split{1};
            fuel_wt(fuel_index) = str2num(fuel_split{2});
            fuel_energy(fuel_index) = str2num(fuel_split{3}) * 1000;
            fuel_T(fuel_index) = str2num(fuel_split{4});
            fuel_index = fuel_index + 1;
            line_next = fgetl(fileID);
            end

            % Get the oxidants
            while startsWith2(line_next, ' OXIDANT')
            ox_index = size(ox, 1) + 1;
            ox_split = strsplit(strtrim(line_next(length(' OXIDANT '):end)));
            ox{ox_index} = ox_split{1};
            ox_wt(ox_index) = str2num(ox_split{2});
            ox_energy(ox_index) = str2num(ox_split{3}) * 1000;
            ox_T(ox_index) = str2num(ox_split{4});
            line_next = fgetl(fileID);
            end
        elseif startsWith2(line2, ' OXIDANT')
            % Get the oxidants
            line_next = line2;
            while startsWith2(line_next, ' OXIDANT')
            ox_index = size(ox, 1) + 1;
            ox_split = strsplit(strtrim(line_next(length(' OXIDANT '):end)));
            ox{ox_index} = ox_split{1};
            ox_wt(ox_index) = str2num(ox_split{2});
            ox_energy(ox_index) = str2num(ox_split{3}) * 1000;
            ox_T(ox_index) = str2num(ox_split{4});
            line_next = fgetl(fileID);
            end
        else
            ox_index = size(ox, 1) + 1;
            ox_split = strsplit(strtrim(line2));
            ox{ox_index} = ox_split{1};
            ox_wt(ox_index) = str2num(ox_split{2});
            ox_energy(ox_index) = str2num(ox_split{3}) * 1000;
            ox_T(ox_index) = str2num(ox_split{4});
        end


    elseif startsWith2(line, 'P1')
        P1_start = strfind(line, 'R') + 1;
        map('P1') = [map('P1'); 1e5.* str2double(strtrim(split(strtrim(line(P1_start:end)),' ')))'];
        assert(~isempty(map('P1')));
    elseif startsWith2(line, 'T1')
        T1_start = strfind(line, 'K') + 1;
        map('T1') = [map('T1'); str2double(strtrim(split(strtrim(line(T1_start:end)),' ')))'];
        assert(~isempty(map('T1')));
    elseif startsWith2(line, 'H1')
        H1_start = strfind(line, 'G') + 1;
        map('H1') = [map('H1'); 1e3.* str2double(strtrim(split(strtrim(line(H1_start:end)),' ')))'];
        assert(~isempty(map('H1')));
    elseif startsWith2(line, 'M1')
        M1_start = strfind(line, ')') + 1;
        map('M1') = [map('M1'); str2double(strtrim(split(strtrim(line(M1_start:end)),' ')))'];
        assert(~isempty(map('M1')));
    elseif startsWith2(line, 'GAMMA1')
        GAM1_start = strfind(line, 'A1') + 2;
        map('GAM1') = [map('GAM1'); str2double(strtrim(split(strtrim(line(GAM1_start:end)),' ')))'];
        assert(~isempty(map('GAM1')));
    elseif startsWith2(line, 'SON VEL1')
        SV1_start = strfind(line, 'C') + 1;
        map('SON1') = [map('SON1'); str2double(strtrim(split(strtrim(line(SV1_start:end)),' ')))'];
        assert(~isempty(map('SON1')));

    elseif startsWith2(line, 'P, ')
        P_start = strfind(line, 'R') + 1;
        map('P') = [map('P'); 1e5.* str2double(strtrim(split(strtrim(line(P_start:end)),' ')))'];
        assert(~isempty(map('P')));
    elseif startsWith2(line, 'T,')
        T_start = strfind(line, 'K') + 1;
        map('T') = [map('T'); str2double(strtrim(split(strtrim(line(T_start:end)),' ')))'];
        assert(~isempty(map('T')));
    elseif startsWith2(line, 'RHO,')
        RHO_start = strfind(line, 'M') + 1;
        int_var = str2double(strtrim(split(strtrim(line(RHO_start:end)),' ')))';
        for i = 1:(length(int_var))/2
            int_var2(i) = int_var(2*i-1)*10^int_var(2*i);
        end
        map('RHO') = [map('RHO'); int_var2];
        clear('int_var2');
        assert(~isempty(map('RHO')));
    elseif startsWith2(line, 'H,')
        H_start = strfind(line, 'G') + 1;
        map('H') = [map('H'); 1e3.* str2double(strtrim(split(strtrim(line(H_start:end)),' ')))'];
        assert(~isempty(map('H')));
    elseif startsWith2(line, 'U,')
        U_start = strfind(line, 'G') + 1;
        map('U') = [map('U'); 1e3.* str2double(strtrim(split(strtrim(line(U_start:end)),' ')))'];
        assert(~isempty(map('U')));
    elseif startsWith2(line, 'G, K')
        G_start = strfind(line, 'KG') + 2;
        map('G') = [map('G'); 1e3.* str2double(strtrim(split(strtrim(line(G_start:end)),' ')))'];
        if isempty(map('G')) || any(isnan(map('G')) == 1) && length(map('G')) < 2
            map('G') = [map('G'); 1e3.* str2double(strtrim(split(strtrim(line(G_start:end)),'.')))'];
        end
        assert(~isempty(map('G')));
    elseif startsWith2(line, 'S, K')
        S_start = strfind(line, 'K)') + 2;
        map('S') = [map('S'); 1e3.* str2double(strtrim(split(strtrim(line(S_start:end)),' ')))'];
        assert(~isempty(map('S')));
    elseif startsWith2(line, 'M,')
        M_start = strfind(line, ')') + 1;
        map('M') = [map('M'); str2double(strtrim(split(strtrim(line(M_start:end)),' ')))'];
        assert(~isempty(map('M')));
    elseif startsWith2(line, 'MW,')
        MW_start = strfind(line, 'T') + 1;
        map('MW') = [map('MW'); str2double(strtrim(split(strtrim(line(MW_start:end)),' ')))'];
        assert(~isempty(map('MW')));
    elseif startsWith2(line, '(dLV/dLP')
        dlv_dlp_t_start = strfind(line, 't') + 1;
        map('dLV_dLP_t') = [map('dLV_dLP_t'); str2double(strtrim(split(strtrim(line(dlv_dlp_t_start:end)),' ')))'];
        assert(~isempty(map('dLV_dLP_t')));
    elseif startsWith2(line, '(dLV/dLT')
        dlv_dlt_p_start = strfind(line, 'p') + 1;
        map('dLV_dLT_p') = [map('dLV_dLT_p'); str2double(strtrim(split(strtrim(line(dlv_dlt_p_start:end)),' ')))'];
        assert(~isempty(map('dLV_dLT_p')));
    elseif startsWith2(line, 'Cp')
        CP_start = strfind(line, 'K)') + 2;
        map('CP') = [map('CP'); 1e3.* str2double(strtrim(split(strtrim(line(CP_start:end)),' ')))'];
        assert(~isempty(map('CP')));
    elseif startsWith2(line, 'GAMMAs')
        GAMs_start = strfind(line, 's') + 1;
        map('GAMs') = [map('GAMs'); str2double(strtrim(split(strtrim(line(GAMs_start:end)),' ')))'];
        assert(~isempty(map('GAMs')));
    elseif startsWith2(line, 'SON VEL,')
        SON_start = strfind(line, 'C') + 1;
        map('SON') = [map('SON'); str2double(strtrim(split(strtrim(line(SON_start:end)),' ')))'];
        assert(~isempty(map('SON')));
    elseif startsWith2(line, 'VISC')
        VISC_start = strfind(line, 'E') + 1;
        map('VISC') = [map('VISC'); 1e-4.* str2double(strtrim(split(strtrim(line(VISC_start:end)),' ')))'];
        assert(~isempty(map('VISC')));

    elseif startsWith2(line, 'WITH EQ')
        line2 = fgetl(fileID);
        line2 = fgetl(fileID);
        CPEQ_start = strfind(line2, 'K)') + 2;
        map('CP_eq') = [map('CP_eq'); 1e3.* str2double(strtrim(split(strtrim(line2(CPEQ_start:end)),' ')))'];
        assert(~isempty(map('CP_eq')));
        line2 = fgetl(fileID);
        CONDEQ_start = strfind(line2, 'Y') + 1;
        map('COND_eq') = [map('COND_eq'); 1e-1.* str2double(strtrim(split(strtrim(line2(CONDEQ_start:end)),' ')))'];
        assert(~isempty(map('COND_eq')));
        line2 = fgetl(fileID);
        PREQ_start = strfind(line2, 'ER') + 2;
        map('PR_eq') = [map('PR_eq'); str2double(strtrim(split(strtrim(line2(PREQ_start:end)),' ')))'];
        assert(~isempty(map('PR_eq')));
    elseif startsWith2(line, 'WITH FR')
        line2 = fgetl(fileID);
        line2 = fgetl(fileID);
        CPFR_start = strfind(line2, 'K)') + 2;
        map('CP_fr') = [map('CP_fr'); 1e3.* str2double(strtrim(split(strtrim(line2(CPFR_start:end)),' ')))'];
        assert(~isempty(map('CP_fr')));
        line2 = fgetl(fileID);
        CONDFR_start = strfind(line2, 'Y') + 1;
        map('COND_fr') = [map('COND_fr'); 1e-1.* str2double(strtrim(split(strtrim(line2(CONDFR_start:end)),' ')))'];
        assert(~isempty(map('COND_fr')));
        line2 = fgetl(fileID);
        PRFR_start = strfind(line2, 'ER') + 2;
        map('PR_fr') = [map('PR_fr'); str2double(strtrim(split(strtrim(line2(PRFR_start:end)),' ')))'];
        assert(~isempty(map('PR_fr')));

    elseif startsWith2(line, 'P/P1')
        PP1_start = strfind(line, 'P1') + 2;
        map('P_P1') = [map('P_P1'); str2double(strtrim(split(strtrim(line(PP1_start:end)),' ')))'];
        assert(~isempty(map('P_P1')));
    elseif startsWith2(line, 'T/T1')
        TT1_start = strfind(line, 'T1') + 2;
        map('T_T1') = [map('T_T1'); str2double(strtrim(split(strtrim(line(TT1_start:end)),' ')))'];
        assert(~isempty(map('T_T1')));
    elseif startsWith2(line, 'M/M1')
        MM1_start = strfind(line, 'M1') + 2;
        map('M_M1') = [map('M_M1'); str2double(strtrim(split(strtrim(line(MM1_start:end)),' ')))'];
        assert(~isempty(map('M_M1')));
    elseif startsWith2(line, 'RHO/RHO1')
        RR1_start = strfind(line, 'O1') + 2;
        map('RHO_RHO1') = [map('RHO_RHO1'); str2double(strtrim(split(strtrim(line(RR1_start:end)),' ')))'];
        assert(~isempty(map('RHO_RHO1')));
    elseif startsWith2(line, 'DET MACH')
        DM_start = strfind(line, 'R') + 1;
        map('DET_MACH') = [map('DET_MACH'); str2double(strtrim(split(strtrim(line(DM_start:end)),' ')))'];
        assert(~isempty(map('DET_MACH')));
    elseif startsWith2(line, 'DET VEL')
        DV_start = strfind(line, 'C') + 1;
        map('DET_VEL') = [map('DET_VEL'); str2double(strtrim(split(strtrim(line(DV_start:end)),' ')))'];
        assert(~isempty(map('DET_VEL')));
    elseif startsWith2(line, 'NOTE.')
        % End of a section
        map_append = containers.Map({'fuel', 'fuel_wt', 'fuel_energy', ...
            'fuel_t', 'ox', 'ox_wt', 'ox_energy', 'ox_t'}, {fuel, ...
            fuel_wt, fuel_energy, fuel_T, ox, ox_wt, ox_energy, ox_T});
        map = [map; map_append];

    elseif startsWith2(line, 'MOLE FRACTIONS')
        line2 = fgetl(fileID);
        while ~startsWith2(line2, '* THERMODYNAMIC')
            line2 = fgetl(fileID);
            key_cutoff = 16;
            if length(line2) < key_cutoff
                line2 = strtrim(fgetl(fileID));
                continue;
            end
            key = line2(1:key_cutoff);
            key = cea_rocket_get_key_name(key);
            map(key) = str2double(strtrim(split(strtrim(line2(DV_start:end)),' ')))';
            assert(~isempty(map(key)));
        end

    end
    % Get a new line without trimming just in case there is no line to get.
    % The new line will be trimmed at the top of the loop.
    line = fgetl(fileID);
end

% Create the type data map if it doesn't exist
if ~isKey(data, map('type'))
    data(map('type')) = containers.Map();
end

% Adjust and insert the new values into the data_type variable
int_val = map('T1');
map('T1') = int_val(~isnan(map('T1')));
data_type = data(map('type'));
map_keys = union(keys(map), keys(data_type));
key_blacklist = {'%fuel', 'fuel_wt', 'fuel_energy', 'fuel_t', 'ox_wt','ox_energy', 'ox_t', 'o/f', 'type','fuel','ox', 'phi', 'r'};

new_P1_row = false;
if ~isKey(data_type, 'T1')
    data_type('T1') = [];
end
P1_index = length(map('T1'))+1;
data_type_P1 = data_type('T1');
if isempty(data_type_P1)
    data_type_P1 = map('T1');
elseif P1_index > length(data_type_P1)
    data_type_P1 = [data_type_P1(1:P1_index-1), map('T1'), data_type_P1(P1_index:end)];
    new_P1_row = true;
end
data_type('T1') = data_type_P1;

for i = 1:length(map_keys)
    key = map_keys{i};
    if ~isempty(strmatch(key, key_blacklist, 'exact'))
        % Skip blacklisted keys
        continue
    end
    if ~isKey(map, key)
        % data_type has the key but map does not. This can end up
        % happening with molar concentrations. Fix this by adding
        % an array of zeros according to the size of the
        % ratio_type.
        map(key) = zeros(1, length(map('T1')));
    end
    if ~isvector(map(key)) || iscell(map(key)) ...
            || length(map(key)) ~= length(map('T1'))
        if isnumeric(map(key))
            if isempty(map(key))
                map(key) = zeros(1, length(map('T1')));
            else
                int_val = map(key);
                map(key) = int_val(~isnan(map(key)));
            end
        else
            continue;
        end
    end
    if ~isKey(data_type, key) && ~isKey(data_type, 'visc')
        % This key doesn't exist, so create the 3D matrix. 'visc'
        % is the last key, so only new keys for the first run go
        % here. If a molar concentration isn't in the first run,
        % but is in subsequent runs, it gets dealt with in the else
        % statement.

        vals3d = zeros(1, 1, length(map('T1')));
        vals3d(1, 1, :) = map(key);
    else
        if ~isKey(data_type, key)
            % Deal with molar concentration problems as explained
            % in the if statement above. Note that MATLAB can't
            % create a 1 x 1 x 1 3D array, it automatically
            % squeezes it down to a 1 x 1 2D array. Therefore, we
            % need some logic to account for MATLAB's shortcomings.
            ref_size = size(data_type('visc'));
            if length(ref_size) < 3
                ref_size(3) = 1;
            end
            vals3d = zeros(ref_size(1), ref_size(2), ref_size(3));
        else
            % Otherwise start vals3d as the old data
            vals3d = data_type(key);
        end

        % Create the pin row if applicable
        if new_P1_row
            vals3d_begin = vals3d(1:P1_index-1, :, :);
            if isempty(vals3d_begin)
                vals3d_begin = [];
            end
            vals3d_end = vals3d(P1_index:end, :, :);
            if isempty(vals3d_end)
                vals3d_end = [];
            end
            vals3d = cat(1, vals3d_begin, ...
                zeros(1, size(vals3d, 2), size(vals3d, 3)), ...
                vals3d_end);
        end

        % Insert the ratio_type rows
        vals3d_original_size = size(vals3d, 3);
        vals3d = cat(3, vals3d, zeros(size(vals3d, 1), size(vals3d, 2), length(map('T1'))));
        vals3d(1, 1, vals3d_original_size+1:end) = map(key);
    end
    data_type(key) = vals3d;
end

% data(map('type')) = map;
if ~found_a_section
    % Notify the user if no sections were found
    error('No sections were found in "%s". Check your input file.', ...
        file_name);
end
% Close the file so it can be deleted if needed
fclose(fileID);