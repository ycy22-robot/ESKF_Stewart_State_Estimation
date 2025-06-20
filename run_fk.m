function estimated_states = run_fk(rod_data)
% RUN_FK Calculates Stewart platform end-effector states and derivatives using forward kinematics.
%        Uses causal (online) differentiation and 2nd-order Butterworth filter to estimate velocity,
%        acceleration, and angular velocity.
%
% Inputs:
%   rod_data: n_samples×7 matrix, each row in the format:
%             [timestamp, L1, L2, L3, L4, L5, L6]
% Output:
%   estimated_states: struct with the following fields:
%      timestamp        - Timestamps (Nx1)
%      position         - Position (Nx3) - from FwKine (m)
%      velocity         - Velocity (Nx3) - causal diff + Butterworth filtered (m/s)
%      acceleration     - Acceleration (Nx3) - causal diff + Butterworth filtered (m/s^2)
%      euler            - RPY Euler angles (Nx3), in degrees - from FwKine
%      orientation      - Rotation matrices (cell array, each is 3×3)
%      angular_velocity - Angular velocity (Nx3), in deg/s - causal diff + Butterworth filtered
%
% Dependencies: Signal Processing Toolbox (for butter, filter)
%               FwKineQ function
%               eul2rotm_custom function

% --- Input Validation ---
if ~exist('butter', 'file') || ~exist('filter', 'file')
    error('Signal Processing Toolbox is required for Butterworth filtering (butter, filter).');
end
if size(rod_data, 2) ~= 7
    error('Input rod_data must have 7 columns: [timestamp, L1, L2, L3, L4, L5, L6]');
end
if size(rod_data, 1) < 2 % Need at least 2 points for backward difference
    error('Need at least 2 data samples for calculations.');
end

% --- Filter Parameters ---
original_cutoff_freq = 18; % Original cutoff frequency (Hz)
filter_order = 3; % Butterworth filter order

% --- Calculate Sampling Frequency ---
timestamp = rod_data(:,1);
dt_vec = diff(timestamp);
if any(dt_vec <= 1e-9)
   warning('Non-increasing or very small timestamps detected. Using mean of positive differences for sampling frequency.');
   positive_dt = dt_vec(dt_vec > 1e-9);
   if isempty(positive_dt)
       error('Cannot determine a valid sampling interval from timestamps.');
   end
   dt_mean = mean(positive_dt);
else
    dt_mean = mean(dt_vec);
end
fs = 1 / dt_mean; % Average sampling frequency (Hz)
nyquist = fs / 2;
fprintf('Estimated sampling frequency: %.2f Hz, Nyquist frequency: %.2f Hz\n', fs, nyquist);

% --- Dynamically Adjust Cutoff Frequency ---
% Cutoff should not exceed Nyquist, and is recommended to be well below it
adjusted_cutoff_freq = min(original_cutoff_freq, 0.45 * nyquist);
if adjusted_cutoff_freq < original_cutoff_freq
    warning('Sampling frequency is low, cutoff freq adjusted from %.2f Hz to %.2f Hz (not exceeding 0.45x Nyquist)', ...
            original_cutoff_freq, adjusted_cutoff_freq);
end
cutoff_freq = adjusted_cutoff_freq;

% --- Design Butterworth Low-pass Filter ---
[b_lp, a_lp] = butter(filter_order, cutoff_freq / nyquist, 'low');
fprintf('Designed %d-order Butterworth low-pass filter with cutoff %.2f Hz (fs = %.2f Hz).\n', ...
         filter_order, cutoff_freq, fs);

% Number of samples
n_samples = size(rod_data, 1);

% --- Pre-allocate Output Variables ---
position    = zeros(n_samples, 3);
euler       = zeros(n_samples, 3); % Store original Euler angles in radians
orientation = cell(n_samples,1);
velocity    = zeros(n_samples, 3);
acceleration= zeros(n_samples, 3);
euler_rate  = zeros(n_samples, 3);
angular_velocity = zeros(n_samples, 3); % rad/s initially

% --- Initialize Filter States ---
% For causal 'filter', we need initial conditions. Zero is standard.
% The state vector length is max(length(a), length(b)) - 1
zi_vel = zeros(max(length(a_lp), length(b_lp)) - 1, 3); % States for velocity filter (x,y,z)
zi_acc = zeros(max(length(a_lp), length(b_lp)) - 1, 3); % States for acceleration filter (x,y,z)
zi_eul = zeros(max(length(a_lp), length(b_lp)) - 1, 3); % States for euler_rate filter (r,p,y)

% --- Process Data Sample by Sample ---
disp('Processing data online (Forward Kinematics, Differentiation, Filtering)...');
q_prev = []; % Initialize previous q for FwKineQ
position_prev = [];
euler_prev = [];
timestamp_prev = [];
velocity_prev = zeros(1,3); % Store previous filtered velocity for accel calculation

for i = 1:n_samples
    % --- Step 1: Forward Kinematics ---
    d = rod_data(i, 2:7)';
    if i == 1 || isempty(q_prev)
        q = FwKineQ(d);
        % Initialize previous values for first iteration's difference calculation
        position_prev = q(1:3)';
        euler_prev = q(4:6)';
        timestamp_prev = timestamp(i) - dt_mean; % Estimate a previous timestamp
        velocity(i,:) = 0; % Assume zero initial velocity
        acceleration(i,:) = 0; % Assume zero initial acceleration
        euler_rate(i,:) = 0; % Assume zero initial rates
    else
        q = FwKineQ(d, q_prev);
    end

    position(i,:) = q(1:3)';
    euler(i,:)    = q(4:6)'; % Store in radians
    orientation{i} = eul2rotm_custom(q(4:6));

    % --- Step 2: Differentiation (Backward Difference) & Filtering ---
    if i > 1
        dt = timestamp(i) - timestamp_prev;
        if dt <= 1e-9 % Avoid division by zero or instability with tiny dt
            warning('run_FK:TimeStepIssue', 'Time step is zero or negative at index %d. Using previous values.', i);
            velocity(i,:) = velocity(i-1,:);
            euler_rate(i,:) = euler_rate(i-1,:);
            acceleration(i,:) = acceleration(i-1,:);
            % Note: Filter states zi_* are NOT updated when dt is invalid
        else
            % 2a: Velocity
            vel_noisy = (position(i,:) - position_prev) / dt;
            [velocity(i,:), zi_vel] = filter(b_lp, a_lp, vel_noisy, zi_vel, 1); % dim=1 for row vector input

            % 2b: Euler Rates
            % Handle wrapping: calculate difference, wrap it, then divide by dt
            euler_diff = euler(i,:) - euler_prev;
            euler_diff_wrapped = atan2(sin(euler_diff), cos(euler_diff)); % Wrap difference to [-pi, pi]
            eul_rate_noisy = euler_diff_wrapped / dt;
            [euler_rate(i,:), zi_eul] = filter(b_lp, a_lp, eul_rate_noisy, zi_eul, 1);

            % 2c: Acceleration (differentiate FILTERED velocity)
            accel_noisy = (velocity(i,:) - velocity_prev) / dt;
            [acceleration(i,:), zi_acc] = filter(b_lp, a_lp, accel_noisy, zi_acc, 1);
        end
    end % end if i > 1

    % --- Step 3: Angular Velocity ---
    phi   = euler(i,1);   % Current roll (rad)
    theta = euler(i,2);   % Current pitch (rad)
    T = [1, 0, -sin(theta);
         0, cos(phi), sin(phi)*cos(theta);
         0, -sin(phi), cos(phi)*cos(theta)];
    angular_velocity(i,:) = (T * euler_rate(i,:)')'; % Use filtered euler_rate

    % --- Update previous values for next iteration ---
    q_prev = q;
    position_prev = position(i,:);
    euler_prev = euler(i,:);
    timestamp_prev = timestamp(i);
    velocity_prev = velocity(i,:); % Store current filtered velocity

    % % Progress indicator (optional)
    % if mod(i, 100) == 0
    %     fprintf('Progress: Processed %d / %d samples.\n', i, n_samples);
    % end

end
disp('Online processing complete.');

%% --- Step 4: Organize Output Structure ---
disp('Finalizing output structure...');
estimated_states = struct();
estimated_states.timestamp = timestamp;
estimated_states.position  = position;       % Original position from FK (m)
estimated_states.velocity  = velocity;       % Filtered velocity (m/s)
estimated_states.acceleration = acceleration; % Filtered acceleration (m/s^2)
estimated_states.euler     = rad2deg(euler); % Original Euler angles (degrees)
estimated_states.orientation = orientation;  % Rotation matrices
estimated_states.angular_velocity = rad2deg(angular_velocity); % Filtered angular velocity (degrees/s)

disp('run_FK finished.');

end

%% Auxiliary function: ZYX Euler angles to Rotation Matrix
function R = eul2rotm_custom(eul)
% EUL2ROTM_CUSTOM Converts ZYX Euler angles to rotation matrix.
% Input: eul = [roll; pitch; yaw] (phi, theta, psi) in radians.
% Output: R (3x3 rotation matrix)
phi   = eul(1); % roll (X-axis rotation)
theta = eul(2); % pitch (Y-axis rotation)
psi   = eul(3); % yaw (Z-axis rotation)

Rx = [1, 0, 0; 0, cos(phi), -sin(phi); 0, sin(phi), cos(phi)];
Ry = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
Rz = [cos(psi), -sin(psi), 0; sin(psi), cos(psi), 0; 0, 0, 1];

% Apply rotations in ZYX order: R = Rz * Ry * Rx
R = Rz * Ry * Rx;
end