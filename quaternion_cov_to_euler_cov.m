function P_euler = quaternion_cov_to_euler_cov(P_quaternion, q)
    % Inputs:
    % P_quaternion: 4x4 covariance matrix of quaternion from EKF
    % q: current estimated quaternion [qw, qx, qy, qz]
    
    % Extract quaternion components
    qw = q(1);
    qx = q(2);
    qy = q(3);
    qz = q(4);
    
    % Compute Euler angles
    % roll (φ)
    roll_numerator = 2 * (qw * qx + qy * qz);
    roll_denominator = 1 - 2 * (qx^2 + qy^2);
    roll = atan2(roll_numerator, roll_denominator);
    
    % pitch (θ)
    pitch_argument = 2 * (qw * qy - qz * qx);
    % Ensure the value is within [-1, 1] to avoid numerical issues
    pitch_argument = max(min(pitch_argument, 1), -1);
    pitch = asin(pitch_argument);
    
    % yaw (ψ)
    yaw_numerator = 2 * (qw * qz + qx * qy);
    yaw_denominator = 1 - 2 * (qy^2 + qz^2);
    yaw = atan2(yaw_numerator, yaw_denominator);
    
    % Compute the Jacobian matrix from quaternion to Euler angles
    J = zeros(3, 4);
    
    % For roll (φ)
    droll_dqw = 2 * qx / (roll_numerator^2 + roll_denominator^2);
    droll_dqx = 2 * qw / (roll_numerator^2 + roll_denominator^2);
    droll_dqy = 2 * qz / (roll_numerator^2 + roll_denominator^2) - 4 * qy * roll_numerator / (roll_numerator^2 + roll_denominator^2);
    droll_dqz = 2 * qy / (roll_numerator^2 + roll_denominator^2);
    
    % For pitch (θ)
    dpitch_dqw = 2 * qy / sqrt(1 - pitch_argument^2);
    dpitch_dqx = -2 * qz / sqrt(1 - pitch_argument^2);
    dpitch_dqy = 2 * qw / sqrt(1 - pitch_argument^2);
    dpitch_dqz = -2 * qx / sqrt(1 - pitch_argument^2);
    
    % For yaw (ψ)
    dyaw_dqw = 2 * qz / (yaw_numerator^2 + yaw_denominator^2);
    dyaw_dqx = 2 * qy / (yaw_numerator^2 + yaw_denominator^2);
    dyaw_dqy = 2 * qx / (yaw_numerator^2 + yaw_denominator^2) - 4 * qy * yaw_denominator / (yaw_numerator^2 + yaw_denominator^2);
    dyaw_dqz = 2 * qw / (yaw_numerator^2 + yaw_denominator^2) - 4 * qz * yaw_denominator / (yaw_numerator^2 + yaw_denominator^2);
    
    % Assemble the Jacobian matrix
    J(1,:) = [droll_dqw, droll_dqx, droll_dqy, droll_dqz];
    J(2,:) = [dpitch_dqw, dpitch_dqx, dpitch_dqy, dpitch_dqz];
    J(3,:) = [dyaw_dqw, dyaw_dqx, dyaw_dqy, dyaw_dqz];
    
    % Apply the error propagation formula
    P_euler = J * P_quaternion * J';
end