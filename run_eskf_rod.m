function [estimated_states] = run_eskf_rod(imu_data, rod_data, options)
% RUN_ESKF_ROD  ESKF state estimation based on rod length observation.
%   imu_data: [timestamp, gyro(3), acc(3)], each row
%   rod_data: [timestamp, rod1, ..., rod6], each row
%   options:  config struct
%   Returns:  estimated_states with fields:
%     timestamp, position, velocity, orientation(cell), euler, rod_lengths,
%     gyro_bias, acce_bias, acceleration, angular_velocity, filtered fields

if nargin < 3, options = struct(); end
eskf = eskf_initialize(options);

n_samples = size(imu_data, 1);
estimated_states = struct(...
    'timestamp',         zeros(n_samples, 1), ...
    'position',          zeros(n_samples, 3), ...
    'velocity',          zeros(n_samples, 3), ...
    'orientation',       {cell(n_samples, 1)}, ...
    'euler',             zeros(n_samples, 3), ...
    'rod_lengths',       zeros(n_samples, 6), ...
    'gyro_bias',         zeros(n_samples, 3), ...
    'acce_bias',         zeros(n_samples, 3), ...
    'acceleration',      zeros(n_samples, 3), ...
    'angular_velocity',  zeros(n_samples, 3), ...
    'P_history',         {cell(n_samples, 1)});

rod_inited = false;
result_count = 0;
rod_idx = 1;

for i = 1:n_samples
    imu = struct('timestamp', imu_data(i, 1), ...
                'gyro', imu_data(i, 2:4)', ...
                'acce', imu_data(i, 5:7)');
    % Initialization by rod observation
    if ~rod_inited
        while rod_idx <= size(rod_data, 1) && rod_data(rod_idx, 1) <= imu.timestamp
            rod_obs = struct('timestamp', rod_data(rod_idx, 1), ...
                             'lengths', rod_data(rod_idx, 2:7)');
            eskf = eskf_observe_rod_lengths(eskf, rod_obs);
            rod_inited = true;
            result_count = result_count + 1;
            estimated_states.timestamp(result_count) = rod_obs.timestamp;
            estimated_states.position(result_count, :) = eskf.p';
            estimated_states.velocity(result_count, :) = eskf.v';
            estimated_states.orientation{result_count} = eskf.R;
            estimated_states.euler(result_count, :) = rotm2euler(eskf.R)' * 180/pi;
            estimated_states.rod_lengths(result_count, :) = rod_obs.lengths';
            estimated_states.gyro_bias(result_count, :) = eskf.bg';
            estimated_states.acce_bias(result_count, :) = eskf.ba';
            estimated_states.P_history{result_count} = eskf.cov;
            acc_est = eskf.R * (imu.acce - eskf.ba) + eskf.g;
            estimated_states.acceleration(result_count, :) = acc_est';
            ang_vel = imu.gyro - eskf.bg;
            estimated_states.angular_velocity(result_count, :) = rad2deg(ang_vel');
            rod_idx = rod_idx + 1;
        end
        if ~rod_inited, continue; end
        continue;
    end

    eskf = eskf_predict(eskf, imu);

    while rod_idx <= size(rod_data, 1) && rod_data(rod_idx, 1) <= imu.timestamp
        rod_obs = struct('timestamp', rod_data(rod_idx, 1), ...
                         'lengths', rod_data(rod_idx, 2:7)');
        eskf = eskf_observe_rod_lengths(eskf, rod_obs);
        rod_idx = rod_idx + 1;
    end

    result_count = result_count + 1;
    estimated_states.timestamp(result_count) = imu.timestamp;
    estimated_states.position(result_count, :) = eskf.p';
    estimated_states.velocity(result_count, :) = eskf.v';
    estimated_states.orientation{result_count} = eskf.R;
    estimated_states.euler(result_count, :) = rotm2euler(eskf.R)' * 180/pi;
    [current_lengths, ~] = calculate_rod_lengths(eskf.p, eskf.R, ...
                                         eskf.platform.fixed_points, ...
                                         eskf.platform.moving_points);
    estimated_states.rod_lengths(result_count, :) = current_lengths';
    estimated_states.gyro_bias(result_count, :) = eskf.bg';
    estimated_states.acce_bias(result_count, :) = eskf.ba';
    est_acc = eskf.R * (imu.acce - eskf.ba) + eskf.g;
    estimated_states.acceleration(result_count, :) = est_acc';
    est_ang_vel = imu.gyro - eskf.bg;
    estimated_states.angular_velocity(result_count, :) = rad2deg(est_ang_vel');
    estimated_states.P_history{result_count} = eskf.cov;
end

% Filter result fields
estimated_states.filtered_velocity = zeros(result_count, 3);
estimated_states.filtered_acceleration = zeros(result_count, 3);
estimated_states.filtered_angular_velocity = zeros(result_count, 3);

median_window = 5;     % Must be odd
sg_window = 11;        % Must be odd
sg_order = 2;
initial_points = max(median_window, sg_window);

for i = 1:min(initial_points, result_count)
    estimated_states.filtered_velocity(i,:) = estimated_states.velocity(i,:);
    estimated_states.filtered_acceleration(i,:) = estimated_states.acceleration(i,:);
    estimated_states.filtered_angular_velocity(i,:) = estimated_states.angular_velocity(i,:);
end

if result_count > initial_points
    for dim = 1:3
        v = estimated_states.velocity(:,dim);
        v_med = medfilt1(v, median_window);
        estimated_states.filtered_velocity(:,dim) = sgolayfilt(v_med, sg_order, sg_window);

        a = estimated_states.acceleration(:,dim);
        a_med = medfilt1(a, median_window);
        estimated_states.filtered_acceleration(:,dim) = sgolayfilt(a_med, sg_order, sg_window);

        w = estimated_states.angular_velocity(:,dim);
        w_med = medfilt1(w, median_window);
        estimated_states.filtered_angular_velocity(:,dim) = sgolayfilt(w_med, sg_order, sg_window);
    end
end
estimated_states.velocity  = estimated_states.filtered_velocity;
estimated_states.acceleration = estimated_states.filtered_acceleration;
estimated_states.angular_velocity = estimated_states.filtered_angular_velocity;
end

function euler = rotm2euler(R)
% Convert 3x3 rotation matrix to Euler angles (ZYX, rad)
pitch = asin(-R(3,1));
if abs(cos(pitch)) > 1e-10
    roll = atan2(R(3,2), R(3,3));
    yaw  = atan2(R(2,1), R(1,1));
else
    roll = 0;
    yaw = atan2(-R(1,2), R(2,2));
end
euler = [roll; pitch; yaw];
end

function [lengths, directions] = calculate_rod_lengths(p, R, fixed_points, moving_points)
% Calculate rod lengths and directions given position and rotation
lengths = zeros(6, 1);
directions = zeros(3, 6);
for i = 1:6
    p_moving_world = R * moving_points(i,:)' + p;
    rod_vector = p_moving_world - fixed_points(i,:)';
    lengths(i) = norm(rod_vector);
    directions(:, i) = rod_vector / lengths(i);
end
end

function eskf = eskf_initialize(options)
% Initialize ESKF estimator
if nargin < 1, options = struct(); end
if ~isfield(options, 'imu_dt'), options.imu_dt = 0.01; end
if ~isfield(options, 'gyro_var'), options.gyro_var = 1e-3; end
if ~isfield(options, 'acce_var'), options.acce_var = 1e-2; end
if ~isfield(options, 'bias_gyro_var'), options.bias_gyro_var = 1e-6; end
if ~isfield(options, 'bias_acce_var'), options.bias_acce_var = 1e-4; end
if ~isfield(options, 'imu_bias_gyro'), options.imu_bias_gyro = zeros(3, 1); end
if ~isfield(options, 'imu_bias_acc'), options.imu_bias_acc = zeros(3, 1); end
if ~isfield(options, 'rod_length_noise'), options.rod_length_noise = 1e-3; end
if ~isfield(options, 'gnss_ang_noise'), options.gnss_ang_noise = 1.0 * pi/180; end
if ~isfield(options, 'update_bias_gyro'), options.update_bias_gyro = true; end
if ~isfield(options, 'update_bias_acce'), options.update_bias_acce = true; end

eskf = struct();
eskf.options = options;
eskf.current_time = 0.0;
eskf.p = zeros(3, 1);
eskf.v = zeros(3, 1);
eskf.R = eye(3);
eskf.bg = eskf.options.imu_bias_gyro;
eskf.ba = eskf.options.imu_bias_acc;
eskf.g = [0; 0; 9.81];
eskf.dx = zeros(18, 1);
eskf.cov = eye(18) * 1e-3;
eskf.cov(4:6, 4:6) = eye(3) * 1e-4;
eskf = build_noise(eskf);
eskf.first_observation = true;
eskf.first_gnss = true;

% Stewart platform geometry parameters
eskf.platform = struct();
eskf.platform.R_b = 800e-3;
eskf.platform.theta_b = 18/2;
eskf.platform.R_a = 500e-3;
eskf.platform.theta_a = 10/2;
eskf.platform.h_mid = 900e-3;
phi_b_0 = [-eskf.platform.theta_b; eskf.platform.theta_b; -eskf.platform.theta_b+120; eskf.platform.theta_b+120; -eskf.platform.theta_b+240; eskf.platform.theta_b+240];
phi_a_0 = [eskf.platform.theta_a; -eskf.platform.theta_a+120; eskf.platform.theta_a+120; -eskf.platform.theta_a+240; eskf.platform.theta_a+240; -eskf.platform.theta_a+360];
phi_b = phi_b_0;
phi_a = phi_a_0 - 60;
eskf.platform.a_A = [eskf.platform.R_a*cosd(phi_a)'; eskf.platform.R_a*sind(phi_a)'; zeros(1,6)];
eskf.platform.b_B = [eskf.platform.R_b*cosd(phi_b)'; eskf.platform.R_b*sind(phi_b)'; zeros(1,6)];
eskf.platform.moving_points = eskf.platform.a_A';
eskf.platform.fixed_points = eskf.platform.b_B';
end

function eskf = build_noise(eskf)
% Build process and observation noise matrices
ev = eskf.options.acce_var;
et = eskf.options.gyro_var;
eg = eskf.options.bias_gyro_var;
ea = eskf.options.bias_acce_var;
Q_diag = [0,0,0, ev,ev,ev, et,et,et, eg,eg,eg, ea,ea,ea, 0,0,0];
eskf.Q = diag(Q_diag);
rod_n = eskf.options.rod_length_noise;
eskf.V_rod = diag(rod_n * ones(1, 6));
end

function eskf = eskf_observe_rod_lengths(eskf, rod_lengths)
% Rod length update for ESKF
if rod_lengths.timestamp < eskf.current_time, return; end
fixed_points = eskf.platform.fixed_points;
moving_points = eskf.platform.moving_points;
if eskf.first_observation
    initial_pose = FwKineQ(rod_lengths.lengths); % User should implement this
    eskf.p = initial_pose(1:3);
    roll = initial_pose(4); pitch = initial_pose(5); yaw = initial_pose(6);
    Rx = [1 0 0; 0 cos(roll) -sin(roll); 0 sin(roll) cos(roll)];
    Ry = [cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)];
    Rz = [cos(yaw) -sin(yaw) 0; sin(yaw) cos(yaw) 0; 0 0 1];
    eskf.R = Rz * Ry * Rx;
    eskf.first_observation = false;
    eskf.current_time = rod_lengths.timestamp;
    return;
end
[predicted_lengths, directions] = calculate_rod_lengths(eskf.p, eskf.R, fixed_points, moving_points);
H = zeros(6, 18);
for i = 1:6
    a_i = moving_points(i,:)';
    s_i_unit = directions(:,i);
    a_i_skew = skew_symmetric(a_i);
    H(i, 1:3) = s_i_unit';
    H(i, 7:9) = -s_i_unit' * eskf.R * a_i_skew;
end
V = eskf.V_rod;
innovation = rod_lengths.lengths - predicted_lengths;
K = eskf.cov * H' / (H * eskf.cov * H' + V);
eskf.dx = eskf.dx + K * innovation;
eskf.cov = (eye(size(eskf.cov)) - K * H) * eskf.cov;
eskf = eskf_update_and_reset(eskf);
eskf.current_time = rod_lengths.timestamp;
end

function eskf = eskf_predict(eskf, imu)
% ESKF prediction step using IMU
dt = imu.timestamp - eskf.current_time;
if dt > (5 * eskf.options.imu_dt) || dt < 0
    eskf.current_time = imu.timestamp;
    return;
end
p = eskf.p; v = eskf.v; R = eskf.R; bg = eskf.bg; ba = eskf.ba; g = eskf.g;
acc_unbias = imu.acce - ba;
gyro_unbias = imu.gyro - bg;
acc_world = R * acc_unbias;
new_p = p + v * dt + 0.5 * (acc_world + g) * dt^2;
new_v = v + (acc_world + g) * dt;
new_R = R * so3_exp(gyro_unbias * dt);
F = eye(18);
F(1:3, 4:6) = eye(3) * dt;
F(4:6, 7:9) = -R * skew_symmetric(acc_unbias) * dt;
F(4:6, 13:15) = -R * dt;
F(4:6, 16:18) = eye(3) * dt;
F(7:9, 7:9) = so3_exp(-gyro_unbias * dt);
F(7:9, 10:12) = -eye(3) * dt;
dx = F * eskf.dx;
cov = F * eskf.cov * F' + eskf.Q;
eskf.p = new_p;
eskf.v = new_v;
eskf.R = new_R;
eskf.dx = dx;
eskf.cov = cov;
eskf.current_time = imu.timestamp;
end

function eskf = eskf_update_and_reset(eskf)
% Inject error state and reset
dp = eskf.dx(1:3);
dv = eskf.dx(4:6);
dtheta = eskf.dx(7:9);
dbg = eskf.dx(10:12);
dba = eskf.dx(13:15);
dg = eskf.dx(16:18);
eskf.p = eskf.p + dp;
eskf.v = eskf.v + dv;
eskf.R = eskf.R * so3_exp(dtheta);
if eskf.options.update_bias_gyro, eskf.bg = eskf.bg + dbg; end
if eskf.options.update_bias_acce, eskf.ba = eskf.ba + dba; end
eskf.g = eskf.g + dg;
eskf = eskf_project_cov(eskf);
eskf.dx = zeros(18, 1);
end

function eskf = eskf_project_cov(eskf)
% Project covariance after injection
J = eye(18);
dtheta = eskf.dx(7:9);
J(7:9, 7:9) = eye(3) - 0.5 * skew_symmetric(dtheta);
eskf.cov = J * eskf.cov * J';
end

function S = skew_symmetric(v)
% Skew-symmetric matrix from vector
S = [    0, -v(3),  v(2);
       v(3),     0, -v(1);
      -v(2),  v(1),     0];
end

function R = so3_exp(w)
% SO(3) exponential map: rotation vector to rotation matrix
theta = norm(w);
if theta < 1e-10, R = eye(3); return; end
a = w / theta;
a_skew = skew_symmetric(a);
R = eye(3) + sin(theta) * a_skew + (1 - cos(theta)) * a_skew^2;
end