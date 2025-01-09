classdef HexArrayGenerator
    methods(Static)
        function [x_coords, y_coords, phases] = generateHexGrid(M, N, p, array_radius, ...
                feed1_x, feed1_y, feed1_z, feed2_x, feed2_y, feed2_z, ...
                l1, l2, fc, theta, pha, theta2, pha2)
            
            % Hexagonal grid parameters
            hex_width = 3*p/2/sqrt(3);
            hex_height = p;
            num_rows = M;
            num_cols = N;
            array_height_mm = num_rows*hex_height;
            array_width_mm = num_cols*hex_width;

            % Generate coordinates using meshgrid
            [row_idx, col_idx] = meshgrid(1:num_rows, 1:num_cols);
            x_coords = -array_width_mm/2 + (col_idx(:) - 1) * hex_width;
            y_coords = array_height_mm/2 - (row_idx(:) - 1) * hex_height;
            y_coords(mod(col_idx(:), 2) == 0) = y_coords(mod(col_idx(:), 2) == 0) + hex_height/2;

            % Calculate radius mask and valid indices
            radius_mask = sqrt(x_coords.^2 + y_coords.^2) < array_radius;
            valid_indices = find(radius_mask);

            % Initialize phase arrays
            num_elements = length(x_coords);
            [phase_l1, phase_l2] = deal(ones(num_elements, 1) * -1);
            [sphere_phase1, sphere_phase2] = deal(ones(num_elements, 1) * -1);
            [steer_phase1, steer_phase2] = deal(ones(num_elements, 1) * -1);
            [total_phase1, total_phase2] = deal(ones(num_elements, 1) * -1);
            phases = ones(num_elements, 1) * -1;

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
                
                % 4. Total phase for each beam
                total_phase1(idx) = mod(phase_l1(idx), 360);  % + sphere_phase1(idx) + steer_phase1(idx)
                total_phase2(idx) = mod(phase_l2(idx), 360);
                
                % 5. Final combined phase
                phases(idx) = mod(angle(exp(1i*deg2rad(total_phase1(idx))) + ...
                              exp(1i*deg2rad(total_phase2(idx))))/pi*180, 360);
            end
        end
        
        function plotHexPhaseDistribution(x_coords, y_coords, phases, p)
            figure('Position', [100, 100, 600, 500]);
            hold on;
            
            % Find valid elements
            valid_mask = phases ~= -1;
            x_valid = x_coords(valid_mask);
            y_valid = y_coords(valid_mask);
            phases_valid = phases(valid_mask);
            
            % Calculate hexagon vertices
            angles = (0:60:360) * pi/180;
            hex_x = p/sqrt(3) * cos(angles);
            hex_y = p/sqrt(3) * sin(angles);
            
            % Plot each hexagonal element
            for i = 1:length(x_valid)
                patch(hex_x + x_valid(i), hex_y + y_valid(i), phases_valid(i), ...
                    'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.1);
            end
            
            % Customize plot appearance
            colormap('jet');
            cb = colorbar;
            cb.Label.String = 'Phase (degrees)';
            cb.Label.FontSize = 12;
            title('Hexagonal Array Phase Distribution', 'FontSize', 14);
            xlabel('X Position (mm)', 'FontSize', 12);
            ylabel('Y Position (mm)', 'FontSize', 12);
            axis equal;
            grid on;
            
            % Set axis limits with some padding
            max_dim = max(abs([x_valid; y_valid])) * 1.1;
            xlim([-max_dim max_dim]);
            ylim([-max_dim max_dim]);
            
            % Set color limits
            caxis([0 360]);
            colorbar('Ticks', 0:60:360, 'TickLabels', {'0°','60°','120°','180°','240°','300°','360°'});
            
            hold off;
        end
    end
end