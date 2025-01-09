% Initialize parameters
clear; 
clc;
close all;

% Constants and Parameters
units = 'mm';
fc = 11.5;  % Frequency in GHz
M = 40;     % Number of rows
N = M;      % Number of columns
p = 12;     % Unit period

% OAM and angle parameters
l1 = -1;    % First OAM mode
l2 = 2;     % Second OAM mode
theta = 20*pi/180;   % First source tilt angle
pha = 0*pi/180;      % First source azimuth angle
theta2 = 20*pi/180;  % Second source tilt angle
pha2 = 180*pi/180;   % Second source azimuth angle
Z3 = 350;           % Source distance

% Calculate source positions
feed1_x = Z3*sin(theta)*cos(pha);
feed1_y = Z3*sin(theta)*sin(pha);
feed1_z = Z3*cos(theta);

feed2_x = Z3*sin(theta2)*cos(pha2);
feed2_y = Z3*sin(theta2)*sin(pha2);
feed2_z = Z3*cos(theta2);

% Hexagonal grid parameters
hex_width = 3*p/2/sqrt(3);
hex_height = p;
num_rows = 40;
num_cols = 40;
array_height_mm = num_rows*hex_height;
array_width_mm = num_cols*hex_width;

% Generate coordinates using meshgrid
[row_idx, col_idx] = meshgrid(1:num_rows, 1:num_cols);
x_coords = -array_width_mm/2 + (col_idx(:) - 1) * hex_width;
y_coords = array_height_mm/2 - (row_idx(:) - 1) * hex_height;
y_coords(mod(col_idx(:), 2) == 0) = y_coords(mod(col_idx(:), 2) == 0) + hex_height/2;

% Calculate radius mask and valid indices
radius_mask = sqrt(x_coords.^2 + y_coords.^2) < 400/2;
valid_indices = find(radius_mask);

% Initialize phase arrays
num_elements = length(x_coords);
[phase_l1, phase_l2] = deal(ones(num_elements, 1) * -1);          % OAM phases
[sphere_phase1, sphere_phase2] = deal(ones(num_elements, 1) * -1); % Spherical wave phases
[steer_phase1, steer_phase2] = deal(ones(num_elements, 1) * -1);  % Steering phases
[total_phase1, total_phase2] = deal(ones(num_elements, 1) * -1);  % Total phases for each beam
apzxy = ones(num_elements, 1) * -1;                               % Final combined phase

% Calculate phases for valid elements
for idx = valid_indices'
    % 1. Spherical wave phases calculation
    r1 = sqrt((x_coords(idx)-feed1_x)^2 + (y_coords(idx)-feed1_y)^2 + (feed1_z)^2)-feed1_z;
    r2 = sqrt((x_coords(idx)-feed2_x)^2 + (y_coords(idx)-feed2_y)^2 + (feed2_z)^2)-feed2_z;
    
    sphere_phase1(idx) = mod(fc*360/300*r1, 360);
    sphere_phase2(idx) = mod(fc*360/300*r2, 360);
    
    % 2. OAM phases calculation
    phase_l1(idx) = -l1 * atan2(x_coords(idx), y_coords(idx));
    phase_l1(idx) = mod(rad2deg(phase_l1(idx)), 360);
    
    phase_l2(idx) = -l2 * atan2(y_coords(idx), x_coords(idx));
    phase_l2(idx) = mod(rad2deg(phase_l2(idx)), 360);
    
    % 3. Steering phases calculation
    steer_phase1(idx) = fc*360/300*(x_coords(idx)*sin(theta)*cos(pha) + ...
                                   y_coords(idx)*sin(theta)*sin(pha));
    steer_phase1(idx) = mod(steer_phase1(idx), 360);
    
    steer_phase2(idx) = fc*360/300*(x_coords(idx)*sin(theta2)*cos(pha2) + ...
                                   y_coords(idx)*sin(theta2)*sin(pha2));
    steer_phase2(idx) = mod(steer_phase2(idx), 360);
    
    % 4. Total phase for each beam (OAM + spherical + steering)
    total_phase1(idx) = mod(phase_l1(idx), 360);% + sphere_phase1(idx) + steer_phase1(idx)
    total_phase2(idx) = mod(phase_l2(idx), 360);
    
    % 5. Final combined phase
    apzxy(idx) = total_phase1(idx);%mod(angle(exp(1i*deg2rad(total_phase1(idx))) +exp(1i*deg2rad(total_phase2(idx))))/pi*180, 360);
end

% Plot the phase distribution
OAMTools.HexArrayGenerator.plotHexPhaseDistribution(x_coords, y_coords, apzxy, p);

% Calculate and plot far-field patterns
options.freq = fc;
options.dynRange = 40;
options.colormap = 'jet';
options.plotType = {'3d', 'cuts', 'contour', 'phase3d', '2d'};

% Calculate far-field pattern
OAMTools.OAMFunctions.calculateFarFieldPattern(x_coords, y_coords, apzxy, options);

% Save the final phase data
final_data = [x_coords(valid_indices), y_coords(valid_indices), apzxy(valid_indices)];
writematrix(final_data, 'phase_data.txt', 'Delimiter', '\t');
disp('Phase data has been saved to phase_data.txt');