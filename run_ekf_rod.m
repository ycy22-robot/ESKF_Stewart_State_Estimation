function [estimated_states] = run_ekf_rod(imu_data, rod_data, options)
% RUN_EKF_ROD EKF state estimation with rod length measurements.
%   imu_data: Nx7 [timestamp, gyro(3), acc(3)]
%   rod_data: Mx7 [timestamp, rod1,...,rod6]
%   options : config struct
%   Returns: estimated_states with fields:
%     timestamp, position, velocity, quaternion, orientation(cell),
%     euler, acceleration, angular_velocity, rod_lengths, gyro_bias, acce_bias

if nargin < 3, options = struct(); end
ekf = ekf_initialize(options);

n_samples = size(imu_data, 1);
estimated_states = struct(...
    'timestamp',         zeros(n_samples, 1), ...
    'position',          zeros(n_samples, 3), ...
    'velocity',          zeros(n_samples, 3), ...
    'quaternion',        zeros(n_samples, 4), ...
    'orientation',       {cell(n_samples, 1)}, ...
    'euler',             zeros(n_samples, 3), ...
    'acceleration',      zeros(n_samples, 3), ...
    'angular_velocity',  zeros(n_samples, 3), ...
    'rod_lengths',       zeros(n_samples, 6), ...
    'gyro_bias',         zeros(n_samples, 3), ...
    'acce_bias',         zeros(n_samples, 3), ...
    'P_history',         {cell(n_samples, 1)});

rod_inited = false;
result_count = 0;
rod_idx = 1;

for i = 1:n_samples
    imu = struct('timestamp', imu_data(i, 1), ...
                 'gyro', imu_data(i, 2:4)', ...
                 'acce', imu_data(i, 5:7)');
    if ~rod_inited
        while rod_idx <= size(rod_data, 1) && rod_data(rod_idx, 1) <= imu.timestamp
            rod_obs = struct('timestamp', rod_data(rod_idx, 1), ...
                             'lengths', rod_data(rod_idx, 2:7)');
            ekf = ekf_observe_rod(ekf, rod_obs);
            ekf.R = quat2rotm(ekf.q');
            rod_inited = true;
            result_count = result_count + 1;
            estimated_states.timestamp(result_count) = rod_obs.timestamp;
            estimated_states.position(result_count, :) = ekf.p';
            estimated_states.velocity(result_count, :) = ekf.v';
            estimated_states.quaternion(result_count, :) = ekf.q';
            estimated_states.orientation{result_count} = ekf.R;
            estimated_states.P_history{result_count} = ekf.P;
            euler = rotm2euler(ekf.R);
            estimated_states.euler(result_count, :) = rad2deg(euler');
            acc_est = ekf.R * (imu.acce - ekf.ba) + ekf.g;
            estimated_states.acceleration(result_count, :) = acc_est';
            ang_vel = imu.gyro - ekf.bg;
            estimated_states.angular_velocity(result_count, :) = rad2deg(ang_vel');
            estimated_states.rod_lengths(result_count, :) = rod_obs.lengths';
            estimated_states.gyro_bias(result_count, :) = ekf.bg';
            estimated_states.acce_bias(result_count, :) = ekf.ba';
            rod_idx = rod_idx + 1;
        end
        if ~rod_inited, continue; end
        continue;
    end

    ekf = ekf_predict(ekf, imu);

    while rod_idx <= size(rod_data, 1) && rod_data(rod_idx, 1) <= imu.timestamp
        rod_obs = struct('timestamp', rod_data(rod_idx, 1), ...
                         'lengths', rod_data(rod_idx, 2:7)');
        ekf = ekf_observe_rod(ekf, rod_obs);
        rod_idx = rod_idx + 1;
    end

    ekf.R = quat2rotm(ekf.q');
    result_count = result_count + 1;
    estimated_states.timestamp(result_count) = imu.timestamp;
    estimated_states.position(result_count, :) = ekf.p';
    estimated_states.velocity(result_count, :) = ekf.v';
    estimated_states.quaternion(result_count, :) = ekf.q';
    estimated_states.orientation{result_count} = ekf.R;
    estimated_states.P_history{result_count} = ekf.P;
    euler = rotm2euler(ekf.R);
    estimated_states.euler(result_count, :) = rad2deg(euler');
    [current_lengths, ~] = calculate_rod_lengths(ekf.p, ekf.R, ...
                                ekf.platform.fixed_points, ekf.platform.moving_points);
    estimated_states.rod_lengths(result_count, :) = current_lengths';
    estimated_states.gyro_bias(result_count, :) = ekf.bg';
    estimated_states.acce_bias(result_count, :) = ekf.ba';
    estimated_acc = ekf.R * (imu.acce - ekf.ba) + ekf.g;
    estimated_states.acceleration(result_count, :) = estimated_acc';
    estimated_states.angular_velocity(result_count, :) = rad2deg((imu.gyro - ekf.bg)');
end

% Optional: Filtering (median + SG)
estimated_states.filtered_velocity = zeros(result_count, 3);
estimated_states.filtered_acceleration = zeros(result_count, 3);
estimated_states.filtered_angular_velocity = zeros(result_count, 3);

median_window = 5;
sg_window = 11;
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
pitch = asin(-R(3,1));
if abs(cos(pitch)) > 1e-10
    roll = atan2(R(3,2), R(3,3));
    yaw = atan2(R(2,1), R(1,1));
else
    roll = 0;
    yaw = atan2(-R(1,2), R(2,2));
end
euler = [roll; pitch; yaw];
end

function [lengths, directions] = calculate_rod_lengths(p, R, fixed_points, moving_points)
lengths = zeros(6, 1);
directions = zeros(3, 6);
for i = 1:6
    p_moving_world = R * moving_points(i,:)' + p;
    rod_vector = p_moving_world - fixed_points(i,:)';
    lengths(i) = norm(rod_vector);
    directions(:, i) = rod_vector / lengths(i);
end
end

function R = quat2rotm(q)
qw = q(1); qx = q(2); qy = q(3); qz = q(4);
R = [1-2*(qy^2+qz^2),   2*(qx*qy - qw*qz), 2*(qx*qz+qw*qy);
     2*(qx*qy+qw*qz), 1-2*(qx^2+qz^2),   2*(qy*qz - qw*qx);
     2*(qx*qz-qw*qy), 2*(qy*qz+qw*qx),   1-2*(qx^2+qy^2)];
end

function ekf = ekf_initialize(options)
if nargin < 1, options = struct(); end
if ~isfield(options, 'imu_dt'), options.imu_dt = 0.01; end
if ~isfield(options, 'gyro_var'), options.gyro_var = 1e-3; end
if ~isfield(options, 'acce_var'), options.acce_var = 1e-2; end
if ~isfield(options, 'bias_gyro_var'), options.bias_gyro_var = 1e-6; end
if ~isfield(options, 'bias_acce_var'), options.bias_acce_var = 1e-4; end
if ~isfield(options, 'imu_bias_gyro'), options.imu_bias_gyro = zeros(3, 1); end
if ~isfield(options, 'imu_bias_acc'), options.imu_bias_acc = zeros(3, 1); end
if ~isfield(options, 'rod_length_noise'), options.rod_length_noise = 1e-3; end
if ~isfield(options, 'update_bias_gyro'), options.update_bias_gyro = true; end
if ~isfield(options, 'update_bias_acce'), options.update_bias_acce = true; end

ekf = struct();
ekf.options = options;
ekf.current_time = 0.0;
ekf.p = zeros(3, 1);
ekf.v = zeros(3, 1);
ekf.q = [1; 0; 0; 0];
ekf.bg = ekf.options.imu_bias_gyro;
ekf.ba = ekf.options.imu_bias_acc;
ekf.g = [0; 0; 9.81];
ekf.P = eye(16) * 1e-2;
ekf.P(7:10, 7:10) = eye(4) * 1e-4;
ekf = build_noise_ekf(ekf);
ekf.first_observation = true;
ekf.platform = struct();
ekf.platform.R_b = 800e-3;
ekf.platform.theta_b = 18/2;
ekf.platform.R_a = 500e-3;
ekf.platform.theta_a = 10/2;
ekf.platform.h_mid = 900e-3;
phi_b_0 = [-ekf.platform.theta_b; ekf.platform.theta_b; -ekf.platform.theta_b+120; ekf.platform.theta_b+120; -ekf.platform.theta_b+240; ekf.platform.theta_b+240];
phi_a_0 = [ekf.platform.theta_a; -ekf.platform.theta_a+120; ekf.platform.theta_a+120; -ekf.platform.theta_a+240; ekf.platform.theta_a+240; -ekf.platform.theta_a+360];
phi_b = phi_b_0;
phi_a = phi_a_0 - 60;
ekf.platform.a_A = [ekf.platform.R_a*cosd(phi_a)'; ekf.platform.R_a*sind(phi_a)'; zeros(1,6)];
ekf.platform.b_B = [ekf.platform.R_b*cosd(phi_b)'; ekf.platform.R_b*sind(phi_b)'; zeros(1,6)];
ekf.platform.moving_points = ekf.platform.a_A';
ekf.platform.fixed_points = ekf.platform.b_B';
end

function ekf = build_noise_ekf(ekf)
options = ekf.options;
gyro = options.gyro_var;
acce = options.acce_var;
bg = options.bias_gyro_var;
ba = options.bias_acce_var;
ekf.Qc = zeros(12, 12);
ekf.Qc(1:3, 1:3) = eye(3) * gyro;
ekf.Qc(4:6, 4:6) = eye(3) * acce;
ekf.Qc(7:9, 7:9) = eye(3) * bg;
ekf.Qc(10:12, 10:12) = eye(3) * ba;
rod = options.rod_length_noise;
ekf.R_rod = eye(6) * rod;
end

function ekf = ekf_observe_rod(ekf, rod)
if rod.timestamp < ekf.current_time, return; end
fixed_points = ekf.platform.fixed_points;
moving_points = ekf.platform.moving_points;
if ekf.first_observation
    initial_pose = FwKineQ(rod.lengths); % Implement this
    ekf.p = initial_pose(1:3);
    roll = initial_pose(4);
    pitch = initial_pose(5);
    yaw = initial_pose(6);
    ekf.q = eul2quat(roll,pitch,yaw);
    ekf.first_observation = false;
    ekf.current_time = rod.timestamp;
    return;
end
p = ekf.p; q = ekf.q; R = quat2R(q);
z = rod.lengths;
z_pred = zeros(6,1);
H = zeros(6, 16);
for i = 1:6
    a_i = moving_points(i,:)';
    b_i = fixed_points(i,:)';
    d_i = p + R*a_i - b_i;
    lambda = norm(d_i);
    if lambda < 1e-8, lambda = 1e-8; end
    z_pred(i) = lambda;
    dldp = d_i' / lambda;
    J_Ra = compute_Ra_dq(q, a_i);
    dldq = (d_i' / lambda) * J_Ra;
    H(i,:) = [ dldp, zeros(1,3), dldq, zeros(1,6) ];
end
residual = z - z_pred;
R_rod = ekf.R_rod;
S = H * ekf.P * H' + R_rod;
K = ekf.P * H' / S;
dx = K * residual;
ekf.p  = ekf.p + dx(1:3);
ekf.v  = ekf.v + dx(4:6);
ekf.q  = ekf.q + dx(7:10); ekf.q = ekf.q / norm(ekf.q);
ekf.bg = ekf.bg + dx(11:13);
ekf.ba = ekf.ba + dx(14:16);
I_KH = eye(16) - K * H;
ekf.P = I_KH * ekf.P * I_KH' + K * R_rod * K';
ekf.P = (ekf.P + ekf.P') / 2;
ekf.current_time = rod.timestamp;
end

function R = quat2R(q)
qw = q(1); qx = q(2); qy = q(3); qz = q(4);
R = [ 1-2*(qy^2+qz^2), 2*(qx*qy - qw*qz), 2*(qx*qz + qw*qy);
      2*(qx*qy + qw*qz), 1-2*(qx^2+qz^2), 2*(qy*qz - qw*qx);
      2*(qx*qz - qw*qy), 2*(qy*qz + qw*qx), 1-2*(qx^2+qy^2) ];
end

function J = compute_Ra_dq(q, a)
qw = q(1); qv = q(2:4);
qv_cross_a = cross(qv, a);
J = zeros(3,4);
J(:,1) = 2*qw * a + 2*qv_cross_a;
for i = 1:3
    e = zeros(3,1); e(i) = 1;
    J(:,i+1) = -2*qv(i)*a + 2*a(i)*qv + 2*(qv'*a) * e + 2*qw*cross(e, a);
end
end

function q = eul2quat(roll, pitch, yaw)
cr = cos(roll/2); sr = sin(roll/2);
cp = cos(pitch/2); sp = sin(pitch/2);
cy = cos(yaw/2); sy = sin(yaw/2);
w = cr * cp * cy + sr * sp * sy;
x = sr * cp * cy - cr * sp * sy;
y = cr * sp * cy + sr * cp * sy;
z = cr * cp * sy - sr * sp * cy;
q = [w; x; y; z];
end

function ekf = ekf_predict(ekf, imu)
dt = imu.timestamp - ekf.current_time;
if dt > (5 * ekf.options.imu_dt) || dt < 0
    ekf.current_time = imu.timestamp;
    return;
end
p = ekf.p; v = ekf.v; q = ekf.q; bg = ekf.bg; ba = ekf.ba; g = ekf.g;
gyro_unbias = imu.gyro - bg;
acc_unbias = imu.acce - ba;
qw = q(1); qx = q(2); qy = q(3); qz = q(4);
R = [1-2*(qy^2+qz^2), 2*(qx*qy-qw*qz), 2*(qx*qz+qw*qy);
     2*(qx*qy+qw*qz), 1-2*(qx^2+qz^2), 2*(qy*qz-qw*qx);
     2*(qx*qz-qw*qy), 2*(qy*qz+qw*qx), 1-2*(qx^2+qy^2)];
acc_world = R * acc_unbias;
new_p = p + v * dt + 0.5 * (acc_world + g) * dt^2;
new_v = v + (acc_world + g) * dt;
new_q = update_quaternion(q, gyro_unbias, dt);
new_bg = bg;
new_ba = ba;
[F, G] = compute_jacobians(q, R, gyro_unbias, acc_unbias);
Fk = eye(16) + F * dt;
Gk = G * dt;
new_P = Fk * ekf.P * Fk' + Gk * ekf.Qc * Gk';
new_P = (new_P + new_P') / 2;
ekf.p = new_p;
ekf.v = new_v;
ekf.q = new_q;
ekf.bg = new_bg;
ekf.ba = new_ba;
ekf.P = new_P;
ekf.current_time = imu.timestamp;
end

function [F, G] = compute_jacobians(q, R, gyro, acc)
F = zeros(16,16);
F(1:3,4:6) = eye(3);
qw = q(1); qx = q(2); qy = q(3); qz = q(4);
qv = [qx; qy; qz];
qv_cross_acc = cross(qv, acc);
acc_skew = [0 -acc(3) acc(2); acc(3) 0 -acc(1); -acc(2) acc(1) 0];
block1 = qw * acc + qv_cross_acc;
term1 = - acc * qv';
term2 = (qv' * acc) * eye(3);
term3 = qv * acc';
term4 = - qw * acc_skew;
block2 = term1 + term2 + term3 + term4;
Fvq = 2 * [block1, block2];
F(4:6,7:10) = Fvq;
F(4:6,14:16) = -R;
wx = gyro(1); wy = gyro(2); wz = gyro(3);
Fqq = 0.5 * [0, -wx, -wy, -wz;
             wx, 0, wz, -wy;
             wy, -wz, 0, wx;
             wz, wy, -wx, 0];
F(7:10,7:10) = Fqq;
Fqbg = -0.5 * [ -qx, -qy, -qz;
                  qw, -qz,  qy;
                  qz,  qw, -qx;
                 -qy,  qx,  qw];
F(7:10,11:13) = Fqbg;
G = zeros(16,12);
G(7:10, 1:3) = Fqbg;
G(4:6, 4:6) = -R;
G(11:13, 7:9) = eye(3);
G(14:16,10:12) = eye(3);
end

function new_q = update_quaternion(q, gyro, dt)
omega_norm = norm(gyro);
if omega_norm > 1e-8
    axis = gyro / omega_norm;
    angle = omega_norm * dt;
    delta_q = [cos(angle/2);
               axis(1) * sin(angle/2);
               axis(2) * sin(angle/2);
               axis(3) * sin(angle/2)];
    new_q = quaternion_multiply(q, delta_q);
    new_q = new_q / norm(new_q);
else
    new_q = q;
end
end

function q_out = quaternion_multiply(q1, q2)
w1 = q1(1); x1 = q1(2); y1 = q1(3); z1 = q1(4);
w2 = q2(1); x2 = q2(2); y2 = q2(3); z2 = q2(4);
w = w1*w2 - x1*x2 - y1*y2 - z1*z2;
x = w1*x2 + x1*w2 + y1*z2 - z1*y2;
y = w1*y2 - x1*z2 + y1*w2 + z1*x2;
z = w1*z2 + x1*y2 - y1*x2 + z1*w2;
q_out = [w; x; y; z];
end