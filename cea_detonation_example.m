% CEA_ROCKET_EXAMPLE: Example file for the MATLAB CEA wrapper. For in-depth
% documentation read the headers of cea_rocket_run.m,
% cea_rocket_run_single.m, and cea_rocket_read.m
clear
clc

% Change this variable to true to rerun CEA instead of using saved values
%addpath('C:\path\to\cea\', '-end');
addpath('C:\Users\Thomas\Desktop\CEA_DET_RKT', '-end');
savepath();
fclose all;

CEA_RUN = true;
CEA_SAVE_FILE = 'cea.mat';

pressure = 25:100:2025;

% The CEA MATLAB code takes a MATLAB map (called a dictionary in Python or
% hash in C) as input. The dictionary uses MATLAB character arrays as the
% keys, and the value data type varies by which key is used. Details of
% each key are listed in cea_rocket_run.m
% For example: inp('key') = value.
inp = containers.Map;
inp('prob') = 'detonation';
inp('type') = 'det';              % Sets the type of CEA calculation
inp('p') = pressure;                % Chamber pressure
inp('p_unit') = 'psi';              % Chamber pressure units
inp('o/f') = 2;               % Mixture ratio
inp('fuel') = 'CH4';             % Fuel name from thermo.inp
inp('fuel_t') = 298;                % Fuel inlet temperature
inp('fuel_wt%') = 100;
inp('ox') = 'H2O2';              % Ox name from thermo.inp
inp('ox_t') = 298;                  % Ox inlet temperature
inp('ox_wt%') = 100;
inp('file_name') = 'DetTest.inp';    % Input/output file name
if CEA_RUN
    data = cea_rocket_run(inp);     % Call the CEA MATLAB code
    data_det = squeeze(data('det'));
    save(CEA_SAVE_FILE, 'data_det');
else
    load(CEA_SAVE_FILE);
end

% Use keys(data_det) to see the contents of each map
% respectively. Every output of CEA is contained in these keys, including
% molar concentrations. Most keys contain a 3D array with columns
% corresponding to the pressure, O/F, and area/pressure ratio inputs
% respectively. If only a single value is given for one of these inputs,
% the output will still be a 3D array. The squeeze() MATLAB function must
% be used to reduce the number of dimensions appropriately. Read the notes
% at the top of cea_rocket_read.m for more details.
temp = squeeze(data_det('T'))';

figure(1)
plot(pressure,temp)
title({'Detonation Test Run','H2O2 & CH4'})
xlabel('Initial Pressure [psi]')
ylabel('Temperature of Burned Gas [K]')
grid on
