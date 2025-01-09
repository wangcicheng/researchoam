classdef OAMFunctions
    methods(Static)
        function calculateFarFieldPattern(x_coords, y_coords, apzxy, options)
            % Set default options if not provided
            if nargin < 4
                options = struct();
            end
            if ~isfield(options, 'freq')
                options.freq = 11.5;
            end
            if ~isfield(options, 'dynRange')
                options.dynRange = 40;
            end
            if ~isfield(options, 'colormap')
                options.colormap = 'jet';
            end
            if ~isfield(options, 'plotType')
                options.plotType = {'3d', 'cuts', 'contour', 'phase3d', '2d'};
            end
            
            % Calculate basic parameters
            theta_scan = -90:1:90;
            phi_scan = 0:2:360;
            [THETA, PHI] = meshgrid(theta_scan, phi_scan);
            theta_rad = deg2rad(THETA);
            phi_rad = deg2rad(PHI);
            
            % Wave number calculation
            c = 3e8;
            fc = options.freq * 1e9;
            lambda = c/fc;
            k = 2*pi/lambda;
            
            % Initialize array factor
            AF = zeros(size(theta_rad));
            
            % Calculate array factor
            valid_indices = find(sqrt(x_coords.^2 + y_coords.^2) < 400/2);
            for idx = valid_indices'
                x = x_coords(idx)/1000;  % Convert to meters
                y = y_coords(idx)/1000;
                phase = apzxy(idx);
                
                phase_term = k * (x * sin(theta_rad) .* cos(phi_rad) + ...
                                y * sin(theta_rad) .* sin(phi_rad));
                AF = AF + exp(1i * (phase_term + deg2rad(phase)));
            end
            
            % Normalize array factor
            AF = AF / max(abs(AF(:)));
            
            % Calculate element pattern
            E_theta = cos(theta_rad);
            E_phi = zeros(size(theta_rad));
            element_pattern = sqrt(E_theta.^2 + E_phi.^2);
            
            % Calculate total field
            E_total = AF .* element_pattern;
            E_total_norm = abs(E_total);
            E_total_db = 20*log10(E_total_norm);
            E_total_db = E_total_db - max(E_total_db(:));
            E_total_phase = mod(rad2deg(angle(E_total)), 360);

            % Generate 2D plots
            if any(strcmp(options.plotType, '2d'))
                % Create new figure for 2D plots
                figure('Position', [100, 100, 800, 400]);
                
                % 2D Amplitude Plot
                subplot(1,2,1);
                surf(PHI, THETA, E_total_db');
                shading interp;
                view(2);
                ylabel('\theta (degree)');
                xlabel('\phi (degree)');
                title('2D Far-field Amplitude Pattern (dB)');
                colorbar;
                caxis([-60 0]);
                axis tight;
                try
                    colormap(gca, flipud(othercolor('RdYlBu11')));
                catch
                    colormap(gca, flipud(jet));
                end
                
                % 2D Phase Plot
                subplot(1,2,2);
                surf(PHI, THETA, E_total_phase');
                shading interp;
                view(2);
                ylabel('\theta (degree)');
                xlabel('\phi (degree)');
                title('2D Far-field Phase Pattern (degree)');
                colorbar;
                caxis([0 360]);
                axis tight;
                colormap(gca, 'hsv');
            end

            % Generate other plots based on options
            if any(strcmp(options.plotType, '3d'))
                % 3D Amplitude Pattern
                [X, Y, Z] = sph2cart(phi_rad, pi/2-theta_rad, E_total_norm);
                
                figure('Position', [100 100 600 600]);
                surf(X, Y, Z, E_total_db);
                shading interp;
                colorbar;
                colormap(options.colormap);
                caxis([-options.dynRange 0]);
                title('3D Far-Field Amplitude Pattern', 'FontSize', 14);
                xlabel('x', 'FontSize', 12);
                ylabel('y', 'FontSize', 12);
                zlabel('z', 'FontSize', 12);
                axis equal;
                grid on;
                view(45, 30);
            end
            
            if any(strcmp(options.plotType, 'phase3d'))
                % 3D Phase Pattern
                [X, Y, Z] = sph2cart(phi_rad, pi/2-theta_rad, E_total_norm);
                
                figure('Position', [100 100 600 600]);
                surf(X, Y, Z, E_total_phase);
                shading interp;
                colorbar;
                colormap('hsv');
                caxis([0 360]);
                title('3D Far-Field Phase Pattern', 'FontSize', 14);
                xlabel('x', 'FontSize', 12);
                ylabel('y', 'FontSize', 12);
                zlabel('z', 'FontSize', 12);
                axis equal;
                grid on;
                view(45, 30);
            end
            
            if any(strcmp(options.plotType, 'cuts'))
                % Cut Planes
                figure('Position', [100 100 1200 400]);
                
                % E-plane cut (φ = 0°)
                subplot(1,2,1);
                phi_idx = find(phi_scan == 0);
                plot(theta_scan, E_total_db(phi_idx,:), 'LineWidth', 2);
                grid on;
                title('E-plane Pattern (φ = 0°)', 'FontSize', 14);
                xlabel('\theta (degrees)', 'FontSize', 12);
                ylabel('Normalized Amplitude (dB)', 'FontSize', 12);
                ylim([-40 0]);
                xlim([-90 90]);
                
                % H-plane cut (φ = 90°)
                subplot(1,2,2);
                phi_idx = find(phi_scan == 90);
                plot(theta_scan, E_total_db(phi_idx,:), 'LineWidth', 2);
                grid on;
                title('H-plane Pattern (φ = 90°)', 'FontSize', 14);
                xlabel('\theta (degrees)', 'FontSize', 12);
                ylabel('Normalized Amplitude (dB)', 'FontSize', 12);
                ylim([-40 0]);
                xlim([-90 90]);
            end
            
            if any(strcmp(options.plotType, 'contour'))
                % Amplitude Contour Plot
                figure('Position', [100 100 600 600]);
                contourf(THETA, PHI, E_total_db, 20);
                colorbar;
                colormap(options.colormap);
                caxis([-options.dynRange 0]);
                title('Far-Field Amplitude Pattern Contour (dB)', 'FontSize', 14);
                xlabel('\theta (degrees)', 'FontSize', 12);
                ylabel('\phi (degrees)', 'FontSize', 12);
                
                % Phase Contour Plot
                figure('Position', [100 100 600 600]);
                contourf(THETA, PHI, E_total_phase, 20);
                colorbar;
                colormap('hsv');
                caxis([0 360]);
                title('Far-Field Phase Pattern Contour', 'FontSize', 14);
                xlabel('\theta (degrees)', 'FontSize', 12);
                ylabel('\phi (degrees)', 'FontSize', 12);
            end
        end
    end
end