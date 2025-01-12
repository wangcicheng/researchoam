% 读取阵列数据
data = readmatrix('data6sides.txt');
x_coords = data(:,1);
y_coords = data(:,2);
phase = data(:,3);

% 去除无效单元（phase=-1的点）
valid_indices = phase ~= -1;
x_coords = x_coords(valid_indices);
y_coords = y_coords(valid_indices);
phase = phase(valid_indices);

% 设置参数
fc = 11.5; % 频率 GHz
c = 3e8;   % 光速
lambda = c/(fc*1e9); % 波长
k = 2*pi/lambda;    % 波数

% 计算六边形单元的面积和等效半径
p = 12; % 单元周期(mm)
hex_area = 2*sqrt(3)*p^2/4;  % 六边形面积
r_eff = sqrt(hex_area/pi);    % 等效圆半径

% 首先绘制相位分布图
figure(1);
% subplot(1,2,1);
scatter(x_coords, y_coords, 50, phase, 'filled');
colorbar;
colormap('jet');
title('Phase Distribution');
xlabel('x (mm)');
ylabel('y (mm)');
axis equal;
grid on;

% % 添加极坐标相位分布
% subplot(1,2,2);
% theta_element = atan2(y_coords, x_coords);
% r_element = sqrt(x_coords.^2 + y_coords.^2);
% polarscatter(theta_element, r_element, 50, phase, 'filled');
title('Phase Distribution (Polar)');
colorbar;

% 计算远场的角度范围
theta = -90:0.5:90;  % theta范围，度
phi = 0:2:360;      % phi范围，度
[THETA, PHI] = meshgrid(theta, phi);

% 转换为弧度
theta_rad = theta*pi/180;
phi_rad = phi*pi/180;
[THETA_RAD, PHI_RAD] = meshgrid(theta_rad, phi_rad);

% 初始化远场
E_total = zeros(length(phi), length(theta));

% 计算每个单元的远场贡献
for i = 1:length(x_coords)
    % 计算单元的位置向量（转换为米）
    r = [x_coords(i)*1e-3, y_coords(i)*1e-3, 0];
    
    % 计算方向余弦
    u = sin(THETA_RAD).*cos(PHI_RAD);
    v = sin(THETA_RAD).*sin(PHI_RAD);
    w = cos(THETA_RAD);
    
    % 计算相位项
    phase_term = k*(u*r(1) + v*r(2));
    
    % % 六边形单元的方向性因子（使用Bessel函数近似）
    % kr = k*r_eff*1e-3*sqrt(u.^2 + v.^2);
    % element_pattern = 2*besselj(1, kr)./(kr + eps);
    % element_pattern(kr == 0) = 1;  % 处理零点
    

    % 偶极子方向图 (cos(theta))
    element_pattern = cos(THETA_RAD);
    % 单元的幅度和相位
    amplitude = 1;  % 统一幅度
    element_phase = phase(i)*pi/180;
    
    % 累加每个单元的贡献
    E_total = E_total + amplitude * element_pattern .* exp(1j*(phase_term + element_phase));
end

% 归一化远场
E_total_norm = abs(E_total)/max(max(abs(E_total)));
E_total_dB = 20*log10(E_total_norm + eps);

% 绘制phi=0的1D方向图
figure(2);
plot(theta, E_total_dB(1,:), 'LineWidth', 1.5);
grid on;
xlabel('\theta (degrees)');
ylabel('Normalized Pattern (dB)');
title('\phi = 0° Cut');
ylim([-40 0]);

% % 绘制2D远场图
% figure(3);
% [X, Y] = pol2cart(PHI_RAD, THETA_RAD);
% surf(X, Y, E_total_dB);
% colorbar;
% shading interp;
% view(2);
% title('2D Far-field Pattern (dB)');
% xlabel('u = sin\theta cos\phi');
% ylabel('v = sin\theta sin\phi');

% 绘制3D远场图
figure(3);
[X, Y, Z] = sph2cart(PHI_RAD, pi/2-THETA_RAD, abs(E_total_norm));
surf(X, Y, Z, E_total_dB);
colorbar;
shading interp;
axis equal;
title('3D Far-field Pattern');
xlabel('x');
ylabel('y');
zlabel('z');
colormap('jet');

view(45, 45);

% 极坐标方向图
figure(4);
polarplot(theta_rad, E_total_dB(1,:));
title('Polar Plot at \phi = 0°');
rlim([-40 0]);
thetalim([-90 90]);
grid on;


% 保存所有图形
for i = 1:4
    figure(i);
    saveas(gcf, sprintf('figure_%d.png', i));
end