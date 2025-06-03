%% Read simulation data and run ESKF algorithm based on rod length observation

clear; clc; close all;
rng(0)

%% 1. Load simulation data files
disp('Reading simulation data files...');

try
    imu_data = readmatrix('imu_data.txt');        % [timestamp, acc(3), gyro(3)]
    rod_data = readmatrix('rod_data.txt');        % [timestamp, rod_lengths(6)]
    ground_truth = readmatrix('ground_truth.txt');% [timestamp, position(3), orientation(3), velocity(3), angular_velocity(3)]
    imu_bias_init = readmatrix('imu_bias_init.txt');
    disp('All data loaded successfully!');
catch ME
    error('Cannot read data files: %s', ME.message);
end

%% 2. Extract variables
disp(['IMU data points: ', num2str(size(imu_data,1))]);
disp(['Rod data points: ', num2str(size(rod_data,1))]);
disp(['Ground truth points: ', num2str(size(ground_truth,1))]);

imu_timestamps = imu_data(:, 1);
imu_acc = imu_data(:, 2:4);
imu_gyr = deg2rad(imu_data(:, 5:7)); % deg/s to rad/s

rod_timestamps = rod_data(:, 1);
L_measured = rod_data(:, 2:7);

gt_timestamps = ground_truth(:, 1);
gt_position = ground_truth(:, 2:4);
gt_orientation = ground_truth(:, 5:7);
gt_velocity = ground_truth(:, 8:10);
gt_angular_velocity = ground_truth(:, 11:13);
gt_acceleration = ground_truth(:, 14:16);

imu_bias_gyro = deg2rad(imu_bias_init(1, :)');
imu_bias_acc = imu_bias_init(2, :)';

%% 3. Print sensor frequency info
imu_dt = mean(diff(imu_timestamps));
rod_dt = mean(diff(rod_timestamps));
disp(['IMU frequency: ', num2str(1/imu_dt), ' Hz']);
disp(['Rod sensor frequency: ', num2str(1/rod_dt), ' Hz']);

%% 4. Visualization
figure('Name', 'Rod Length Data');
plot(rod_timestamps, L_measured);
xlabel('Time (s)'); ylabel('Rod Length (m)');
title('Rod Length Sensor Data');
legend('Leg 1', 'Leg 2', 'Leg 3', 'Leg 4', 'Leg 5', 'Leg 6');
grid on;

figure('Name', 'IMU Data');
subplot(2,1,1);
plot(imu_timestamps, imu_acc);
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
title('IMU Acceleration Data');
legend('a_x', 'a_y', 'a_z');
grid on;
subplot(2,1,2);
plot(imu_timestamps, rad2deg(imu_gyr));
xlabel('Time (s)'); ylabel('Angular Velocity (deg/s)');
title('IMU Angular Velocity Data');
legend('\omega_x', '\omega_y', '\omega_z');
grid on;

figure('Name', 'Ground Truth Trajectory');
subplot(2,1,1);
plot(gt_timestamps, gt_position);
xlabel('Time (s)'); ylabel('Position (m)');
title('True Position Trajectory');
legend('X', 'Y', 'Z');
grid on;
subplot(2,1,2);
plot(gt_timestamps, gt_orientation);
xlabel('Time (s)'); ylabel('Orientation (deg)');
title('True Orientation Trajectory (RPY)');
legend('Roll', 'Pitch', 'Yaw');
grid on;

%% Parameter configuration for simulation
acc_needling_amplitude  = 0.08 * 50;      % Acceleration spike amplitude (m/s^2)
gyro_needling_amplitude = 0.1 * 50;       % Gyro spike amplitude (deg/s)
rod_needling_amplitude  = 2e-5 * 50;      % Rod length spike amplitude (m)
needling_duration       = 0.01;           % Duration per spike (s)
acc_needling_period     = 1.0 * 0.2;      % Spike period for acceleration
gyro_needling_period    = 1.0 * 0.2;      % Spike period for gyro
rod_needling_period     = 2.0 * 0.3;      % Spike period for rod

imu_data_with_noise = imu_data;
rod_data_with_noise = rod_data;

total_time = imu_timestamps(end) - imu_timestamps(1);
num_acc_needlings  = max(1, floor(total_time / acc_needling_period));
num_gyro_needlings = max(1, floor(total_time / gyro_needling_period));
num_rod_needlings  = max(1, floor(total_time / rod_needling_period));

acc_needling_time_start  = imu_timestamps(1);
acc_needling_time_end    = imu_timestamps(end);
gyro_needling_time_start = imu_timestamps(1);
gyro_needling_time_end   = imu_timestamps(end);
rod_needling_time_start  = rod_timestamps(1);
rod_needling_time_end    = rod_timestamps(end);

acc_time_span  = acc_needling_time_end  - acc_needling_time_start;
gyro_time_span = gyro_needling_time_end - gyro_needling_time_start;
rod_time_span  = rod_needling_time_end  - rod_needling_time_start;

%% Add acceleration spike noise
acc_needling_start_times  = sort(acc_needling_time_start  + (acc_time_span  - needling_duration) * rand(num_acc_needlings, 1));
acc_needling_axes = randi(3, num_acc_needlings, 1);

for i = 1:num_acc_needlings
    start_time = acc_needling_start_times(i);
    end_time   = start_time + needling_duration;
    axis_idx   = acc_needling_axes(i);
    needling_indices = find(imu_timestamps >= start_time & imu_timestamps <= end_time);
    if isempty(needling_indices), continue; end
    t_local = linspace(0, pi, length(needling_indices));
    sign_flip = randsample([-1, +1], 1);
    needling_shape = sign_flip * acc_needling_amplitude * sin(t_local);
    imu_data_with_noise(needling_indices, axis_idx+1) = ...
        imu_data_with_noise(needling_indices, axis_idx+1) + needling_shape';
end

%% Add gyro spike noise
gyro_needling_start_times = sort(gyro_needling_time_start + (gyro_time_span - needling_duration) * rand(num_gyro_needlings, 1));
gyro_needling_axes        = randi(3, num_gyro_needlings, 1);

for i = 1:num_gyro_needlings
    start_time = gyro_needling_start_times(i);
    end_time   = start_time + needling_duration;
    axis_idx   = gyro_needling_axes(i);
    needling_indices = find(imu_timestamps >= start_time & imu_timestamps <= end_time);
    if isempty(needling_indices), continue; end
    t_local = linspace(0, pi, length(needling_indices));
    sign_flip = randsample([-1, +1], 1);
    needling_shape = sign_flip * gyro_needling_amplitude * sin(t_local);
    imu_data_with_noise(needling_indices, axis_idx+4) = ...
        imu_data_with_noise(needling_indices, axis_idx+4) + needling_shape';
end

%% Add rod length spike noise
rod_needling_start_times  = sort(rod_needling_time_start  + (rod_time_span  - needling_duration) * rand(num_rod_needlings, 1));
rod_needling_legs        = randi(6, num_rod_needlings, 1);

for i = 1:num_rod_needlings
    start_time = rod_needling_start_times(i);
    end_time   = start_time + needling_duration;
    leg_idx    = rod_needling_legs(i);
    needling_indices = find(rod_timestamps >= start_time & rod_timestamps <= end_time);
    if isempty(needling_indices), continue; end
    t_local = linspace(0, pi, length(needling_indices));
    sign_flip = randsample([-1, +1], 1);
    needling_shape = sign_flip * rod_needling_amplitude * sin(t_local);
    rod_data_with_noise(needling_indices, leg_idx+1) = ...
        rod_data_with_noise(needling_indices, leg_idx+1) + needling_shape';
end

disp('Spike noise added to IMU and rod data.');

%% Visualize noisy data
figure('Position', [100, 100, 1200, 800]);
subplot(2,2,1);
hold on;
plot(imu_timestamps, imu_data(:,2:4), 'LineWidth', 1);
plot(imu_timestamps, imu_data_with_noise(:,2:4), 'LineWidth', 1, 'LineStyle', '--');
title('Acceleration Data (with Spikes)');
xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
legend('Original X', 'Original Y', 'Original Z', 'Noisy X', 'Noisy Y', 'Noisy Z');
grid on;

subplot(2,2,2);
hold on;
plot(imu_timestamps, imu_data(:,5:7), 'LineWidth', 1);
plot(imu_timestamps, imu_data_with_noise(:,5:7), 'LineWidth', 1, 'LineStyle', '--');
title('Gyro Data (with Spikes)');
xlabel('Time (s)'); ylabel('Angular Velocity (rad/s)');
legend('Original X', 'Original Y', 'Original Z', 'Noisy X', 'Noisy Y', 'Noisy Z');
grid on;

subplot(2,2,3);
hold on;
plot(rod_timestamps, rod_data(:,2:4), 'LineWidth', 1);
plot(rod_timestamps, rod_data_with_noise(:,2:4), 'LineWidth', 1, 'LineStyle', '--');
title('Rod Length 1-3 (with Spikes)');
xlabel('Time (s)'); ylabel('Rod Length (m)');
legend('Original L1', 'Original L2', 'Original L3', 'Noisy L1', 'Noisy L2', 'Noisy L3');
grid on;

subplot(2,2,4);
hold on;
plot(rod_timestamps, rod_data(:,5:7), 'LineWidth', 1);
plot(rod_timestamps, rod_data_with_noise(:,5:7), 'LineWidth', 1, 'LineStyle', '--');
title('Rod Length 4-6 (with Spikes)');
xlabel('Time (s)'); ylabel('Rod Length (m)');
legend('Original L4', 'Original L5', 'Original L6', 'Noisy L4', 'Noisy L5', 'Noisy L6');
grid on;

%% Prepare data for algorithms
imu_data = [imu_timestamps, imu_gyr, imu_acc];
rod_data = [rod_timestamps, L_measured];

%% Run ESKF (Rod observation)
options_ESKF = struct(); 
options_ESKF.imu_dt = imu_dt;
options_ESKF.gyro_var = 1.1 * deg2rad(0.1)^2;
options_ESKF.acce_var = 1.1 * (8e-2)^2;
options_ESKF.bias_gyro_var = 1e-6;
options_ESKF.bias_acce_var = 1e-5;
options_ESKF.rod_length_noise = 1.1 * (0.02e-3)^2;
options_ESKF.imu_bias_gyro = imu_bias_gyro;
options_ESKF.imu_bias_acc = imu_bias_acc;

tic
est_states_ESKF = run_eskf_rod(imu_data, rod_data, options_ESKF);
toc
disp('ESKF finished.');

%% Run EKF (Rod observation, nominal params)
options_EKF = struct();
options_EKF.imu_dt = imu_dt;
options_EKF.gyro_var = 20 * deg2rad(0.1)^2;
options_EKF.acce_var = 1.2 * (8e-2)^2;
options_EKF.bias_gyro_var = 1e-6;
options_EKF.bias_acce_var = 1e-5;
options_EKF.rod_length_noise = 50 * (0.02e-3)^2;
options_EKF.imu_bias_gyro = imu_bias_gyro;
options_EKF.imu_bias_acc = imu_bias_acc;

tic
est_states_EKF = run_ekf_rod(imu_data, rod_data, options_EKF);
toc
disp('EKF finished.');

%% Run EBFK (Forward kinematics based numeric method)
tic
est_states_FK = run_fk(rod_data);
toc
disp('EBFK finished.');

%% Save simulation results to .mat file
disp('Saving simulation results...');

ground_truth_data = struct();
ground_truth_data.timestamps = gt_timestamps;
ground_truth_data.position = gt_position;
ground_truth_data.orientation = gt_orientation;
ground_truth_data.velocity = gt_velocity;
ground_truth_data.angular_velocity = gt_angular_velocity;
ground_truth_data.acceleration = gt_acceleration;

imu_measure_data = struct();
imu_measure_data.timestamps = imu_timestamps;
imu_measure_data.gyro = rad2deg(imu_gyr);
imu_measure_data.acce = imu_acc;

leg_measure_data = struct();
leg_measure_data.timestamps = rod_timestamps;
leg_measure_data.L = L_measured;

save_filename = 'simulation_results_nominal.mat';

try
    save(save_filename, ...
         'ground_truth_data', ...
         'imu_measure_data', ...
         'leg_measure_data', ...
         'est_states_ESKF', ...
         'est_states_EKF', ...
         'est_states_FK', ...
         '-v7.3');
    disp(['Simulation results saved to: ', save_filename]);
catch ME_save
    error('Failed to save results: %s', ME_save.message);
end