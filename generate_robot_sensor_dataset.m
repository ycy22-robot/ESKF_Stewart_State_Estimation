%% ========== 0. Custom Frequency and General Settings ==========
clear; clc; close all;
rng(0);

% ========== User Customization ==========
% Multipliers
freq_multi = 1;
freq_noise = 1;

% Encoder noise
encoder_resolution = 0.06e-3 * freq_noise;     % 0.05 mm
encoder_noise_std  = encoder_resolution / 3;   % Standard deviation

% IMU parameters (example)
acc_noise_std  = 8e-2 * freq_noise;         % Accelerometer noise std (m/s^2)
gyro_noise_std = 0.1 * freq_noise;          % Gyro noise std (deg/s)
acc_bias       =  0.015 * randn(1,3) - 0.015;   % Accelerometer bias (m/s^2)
gyro_bias      =  10 * rand(1,3) - 5;           % Gyro bias (deg/s)
acc_bias_noise_std = 2e-5;
gyro_bias_noise_std = 1e-5;

% Save ground truth biases (Note: gyro bias converted to rad/s)
gyro_bias_rad = gyro_bias * (pi/180);  % Convert to rad/s

% Simulation duration
T = 20;    % Total time (s)

% 1) Reference sampling frequency for ground truth trajectory
freq_real = 250;             % Frequency (Hz)
dt_real   = 1 / freq_real;   % Time resolution (s)

% 2) IMU sampling frequency
freq_imu = 250;              % IMU sampling frequency (Hz)
dt_imu   = 1 / freq_imu;     % IMU timestep (s)

% 3) Rod length sensor sampling frequency
freq_rod = 250;              % Rod sampling frequency (Hz)
dt_rod   = 1 / freq_rod;     % Rod sampling timestep (s)

% Generate time sequence (reference)
t_real = 0 : dt_real : T; 
n_real = length(t_real);

% Bias random walk noise std (per step drift, units match bias)
acc_bias_walk_std  = acc_bias_noise_std;   
gyro_bias_walk_std = gyro_bias_noise_std;

% Preallocation
acc_bias_traj  = zeros(n_real, 3);
gyro_bias_traj = zeros(n_real, 3);

% Initialize
acc_bias_traj(1,:)  = acc_bias;
gyro_bias_traj(1,:) = gyro_bias;

% Generate random walk bias trajectory
for i = 2:n_real
    acc_bias_traj(i,:)  = acc_bias_traj(i-1,:)  + acc_bias_walk_std  * randn(1,3);
    gyro_bias_traj(i,:) = gyro_bias_traj(i-1,:) + gyro_bias_walk_std * randn(1,3);
end

% Save bias ground truth
true_biases = struct();
true_biases.acce_bias = acc_bias_traj;    
true_biases.gyro_bias = gyro_bias_traj;   
save('true_bias_values.mat', 'true_biases');

%% ========== 1. Generate End-Effector Trajectory (XYZ + RPY) ==========

% Preallocate
pose   = zeros(n_real, 6);  % [x, y, z, roll, pitch, yaw] (deg)
vel    = zeros(n_real, 6);  % [vx, vy, vz, roll_dot, pitch_dot, yaw_dot]
accel  = zeros(n_real, 6);  % [ax, ay, az, roll_ddot, pitch_ddot, yaw_ddot]
omega  = zeros(n_real, 3);  % [wx, wy, wz] angular velocity (deg/s)

% Position time functions
amplitude_pos = [0.1, 0.1, 0.2];    
frequency_pos = [0.5, 0.7, 0.3]*freq_multi;   
offset_pos    = [0,   0,   1];     

% Orientation time functions
amplitude_ori = [10, 8, 20];       
frequency_ori = [0.5, 0.6, 0.3]*freq_multi;   
offset_ori    = [0,   0,   0];     

for i = 1 : n_real
    ti = t_real(i);

    % Position - XYZ
    pose(i, 1) = offset_pos(1) + amplitude_pos(1)*cos(2*pi*frequency_pos(1)*ti);
    pose(i, 2) = offset_pos(2) + amplitude_pos(2)*cos(2*pi*frequency_pos(2)*ti);
    pose(i, 3) = offset_pos(3) + amplitude_pos(3)*cos(2*pi*frequency_pos(3)*ti);

    % Orientation - RPY (deg)
    pose(i, 4) = offset_ori(1) + amplitude_ori(1)*cos(2*pi*frequency_ori(1)*ti);
    pose(i, 5) = offset_ori(2) + amplitude_ori(2)*cos(2*pi*frequency_ori(2)*ti);
    pose(i, 6) = offset_ori(3) + amplitude_ori(3)*cos(2*pi*frequency_ori(3)*ti);

    % Velocity 
    vel(i, 1) = - amplitude_pos(1) * 2*pi*frequency_pos(1) * sin(2*pi*frequency_pos(1)*ti);
    vel(i, 2) = - amplitude_pos(2) * 2*pi*frequency_pos(2) * sin(2*pi*frequency_pos(2)*ti);
    vel(i, 3) = - amplitude_pos(3) * 2*pi*frequency_pos(3) * sin(2*pi*frequency_pos(3)*ti);

    % Euler angle derivatives (deg/s)
    vel(i, 4) = - amplitude_ori(1) * 2*pi*frequency_ori(1) * sin(2*pi*frequency_ori(1)*ti);
    vel(i, 5) = - amplitude_ori(2) * 2*pi*frequency_ori(2) * sin(2*pi*frequency_ori(2)*ti);
    vel(i, 6) = - amplitude_ori(3) * 2*pi*frequency_ori(3) * sin(2*pi*frequency_ori(3)*ti);

    % Acceleration 
    accel(i, 1) = - amplitude_pos(1) * (2*pi*frequency_pos(1))^2 * cos(2*pi*frequency_pos(1)*ti);
    accel(i, 2) = - amplitude_pos(2) * (2*pi*frequency_pos(2))^2 * cos(2*pi*frequency_pos(2)*ti);
    accel(i, 3) = - amplitude_pos(3) * (2*pi*frequency_pos(3))^2 * cos(2*pi*frequency_pos(3)*ti);

    % Euler angle accelerations (deg/s^2)
    accel(i, 4) = - amplitude_ori(1) * (2*pi*frequency_ori(1))^2 * cos(2*pi*frequency_ori(1)*ti);
    accel(i, 5) = - amplitude_ori(2) * (2*pi*frequency_ori(2))^2 * cos(2*pi*frequency_ori(2)*ti);
    accel(i, 6) = - amplitude_ori(3) * (2*pi*frequency_ori(3))^2 * cos(2*pi*frequency_ori(3)*ti);

    % Angular velocity (convert Euler angle rates to body angular velocity)
    phi       = deg2rad(pose(i, 4));      
    theta     = deg2rad(pose(i, 5));      
    psi       = deg2rad(pose(i, 6));      
    phi_dot   = deg2rad(vel(i, 4));       
    theta_dot = deg2rad(vel(i, 5));       
    psi_dot   = deg2rad(vel(i, 6));       

    transform_matrix = [
        1,          0,         -sin(theta);
        0,    cos(phi),  sin(phi)*cos(theta);
        0,   -sin(phi),  cos(phi)*cos(theta)
    ];
    euler_rates = [phi_dot; theta_dot; psi_dot];

    omega_rad = transform_matrix * euler_rates;
    omega(i, :) = rad2deg(omega_rad);  % (deg/s)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ========== 2. Inverse Kinematics for Rod Lengths (Noiseless & Noisy) ==========

% InvKine(pose) returns [L, s],
% L: 1x6 rod lengths, s: 3x6 rod direction vectors (implement as needed)
num_samples_real = length(t_real);
L_real = zeros(num_samples_real, 6);         
s_real = zeros(3, 6, num_samples_real);      
for i = 1 : num_samples_real
    [L_temp, s_temp] = InvKine(pose(i, :));  % User-defined
    L_real(i, :)     = L_temp;
    s_real(:, :, i)  = s_temp;
end

% Rod length measurement noise (encoder resolution + noise)
L_measured = L_real + encoder_noise_std * randn(size(L_real));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ========== 3. Generate IMU Ground Truth & Noisy Data ==========

% Gravity
g = [0; 0; 9.81];

% Preallocate (correspond to t_real)
imu_acc_true  = zeros(num_samples_real, 3);  % Without noise/bias
imu_gyro_true = zeros(num_samples_real, 3);  % Without noise/bias
imu_acc  = zeros(num_samples_real, 3);       % With noise/bias
imu_gyro = zeros(num_samples_real, 3);       % With noise/bias

for i = 1 : num_samples_real
    % Current orientation (radians)
    roll  = deg2rad(pose(i, 4));
    pitch = deg2rad(pose(i, 5));
    yaw   = deg2rad(pose(i, 6));

    % Rotation matrix (world to body)
    Rx = [1 0 0; 
          0 cos(roll) -sin(roll);
          0 sin(roll)  cos(roll)];
    Ry = [ cos(pitch) 0 sin(pitch);
           0          1 0;
          -sin(pitch) 0 cos(pitch)];
    Rz = [ cos(yaw) -sin(yaw) 0;
           sin(yaw)  cos(yaw) 0;
           0         0        1];
    R = Rz * Ry * Rx;

    % True linear acceleration (world frame, no gravity)
    linear_acc = [accel(i,1); accel(i,2); accel(i,3)];

    % IMU acceleration (body frame): R'*(linear_acc - g)
    imu_acc_true(i,:) = (R' * (linear_acc - g))';

    % Noisy, biased
    imu_acc(i,:) = imu_acc_true(i,:) + acc_bias_traj(i,:) + acc_noise_std * randn(1,3);

    % True angular velocity (deg/s, body frame)
    imu_gyro_true(i,:) = omega(i,:);

    % Noisy, biased
    imu_gyro(i,:) = imu_gyro_true(i,:) + gyro_bias_traj(i,:) + gyro_noise_std * randn(1,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ========== 4. Downsampling/Interpolation for IMU & Rod Data ==========

% 4.1 IMU sampling times
t_imu = 0 : dt_imu : T;  
imu_acc_imu_true  = interp1(t_real, imu_acc_true,  t_imu, 'linear');
imu_gyro_imu_true = interp1(t_real, imu_gyro_true, t_imu, 'linear');
imu_acc_imu       = interp1(t_real, imu_acc,       t_imu, 'linear');
imu_gyro_imu      = interp1(t_real, imu_gyro,      t_imu, 'linear');

% 4.2 Rod length sampling times
t_rod = 0 : dt_rod : T;
L_measured_rod = interp1(t_real, L_measured, t_rod, 'linear');

% 4.3 Ground truth at IMU frequency
pos_imu   = interp1(t_real, pose(:,1:3), t_imu, 'linear');   % Position (m)
ori_imu   = interp1(t_real, pose(:,4:6), t_imu, 'linear');   % Orientation (deg)
vel_imu   = interp1(t_real, vel(:,1:3),  t_imu, 'linear');   % Linear velocity (m/s)
omega_imu = interp1(t_real, omega,       t_imu, 'linear');   % Angular velocity (deg/s)
acc_imu   = interp1(t_real, accel(:,1:3), t_imu, 'linear');  % Acceleration (m/s^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ========== 5. Output Data to Text Files ==========

% 5.1 IMU data: [timestamp, acce(3), gyro(3)]
imu_data = [t_imu', imu_acc_imu, imu_gyro_imu];
writematrix(imu_data, 'imu_data.txt', 'Delimiter', 'tab');

% 5.2 Rod length data: [timestamp, rod_lengths(6)]
rod_data = [t_rod', L_measured_rod];
writematrix(rod_data, 'rod_data.txt', 'Delimiter', 'tab');

% 5.3 Ground truth: [timestamp, position(3), orientation(3), velocity(3), angular_velocity(3)]
ground_truth = [ ...
    t_imu', ...
    pos_imu, ...
    ori_imu, ...
    vel_imu, ...
    omega_imu, ...
    acc_imu ...
];
writematrix(ground_truth, 'ground_truth.txt', 'Delimiter', 'tab');

disp('========== Simulation and sampled data generated and written to text files ==========');

function [L, s] = InvKine(X)
% Input: 
%   X - End-effector pose, [x, y, z, roll, pitch, yaw] 
%       (translation in meters, orientation in degrees, RPY)
%
% Output:
%   L - Vector of rod (leg) lengths (1x6)
%   s - Direction vectors of each rod (3x6)
%
% Structure parameters: 
%   A - Base platform joint coordinates
%   B - Moving platform joint coordinates

p = [X(1); X(2); X(3)];    % Position of the moving platform origin in the base frame

% Structure radii and initial angles (in degrees)
R_a = 800e-3;              % Base platform radius (meters)
theta_a = 18/2;            % Base platform angle (degrees)
R_b = 500e-3;              % Moving platform radius (meters)
theta_b = 10/2;            % Moving platform angle (degrees)

% Orientation (RPY Euler angles, degrees)
R = Eular(deg2rad(X(4:6)));   % Rotation matrix from Euler angles (implement Eular() function)

% Joint positions (angles in degrees)
phi_a_0 = [-theta_a; theta_a; -theta_a+120; theta_a+120; -theta_a+240; theta_a+240];
phi_b_0 = [theta_b; -theta_b+120; theta_b+120; -theta_b+240; theta_b+240; -theta_b+360];
phi_a = phi_a_0;
phi_b = phi_b_0 - 60;

% Coordinates of joints on the base and moving platform
a_A = [R_a*cosd(phi_a)'; R_a*sind(phi_a)'; zeros(1,6)];
b_B = [R_b*cosd(phi_b)'; R_b*sind(phi_b)'; zeros(1,6)];

% Kinematic model: vector from base joint to moving platform joint
d = p + R*b_B - a_A;

% Direction vectors of each rod (unit vectors)
s = d ./ sqrt(d(1,:).^2 + d(2,:).^2 + d(3,:).^2);

% Length of each rod
L = sqrt(d(1,:).^2 + d(2,:).^2 + d(3,:).^2);

end

function R = Eular(RPY)
    a = RPY(1); 
    b = RPY(2);
    c = RPY(3);
    R = [cos(b)*cos(c), -cos(a)*sin(c)+sin(a)*sin(b)*cos(c),  sin(a)*sin(c)+cos(a)*sin(b)*cos(c); 
         cos(b)*sin(c),  cos(a)*cos(c)+sin(a)*sin(b)*sin(c), -sin(a)*cos(c)+cos(a)*sin(b)*sin(c); 
         -sin(b),         sin(a)*cos(b),                       cos(a)*cos(b)                   ];
end