function P_euler = rotation_matrix_cov_to_euler_cov(P_rotation, R)
    % Inputs:
    % P_rotation: 3x3 covariance matrix of the rotation error vector in ESKF
    % R: current estimated rotation matrix

    % Extract Euler angles from the rotation matrix
    roll = atan2(R(3,2), R(3,3));
    pitch = -asin(R(3,1));
    yaw = atan2(R(2,1), R(1,1));

    % Compute the Jacobian matrix from the rotation error vector to Euler angle variation
    % This transformation consists of two steps:
    % 1. Mapping from the rotation error vector to the rotation matrix perturbation
    % 2. Mapping from the rotation matrix perturbation to the Euler angle perturbation

    % Compute partial derivatives of Euler angles with respect to rotation matrix elements
    droll_dR = zeros(1,9);
    dpitch_dR = zeros(1,9);
    dyaw_dR = zeros(1,9);

    % For roll = atan2(R(3,2), R(3,3))
    denom_roll = R(3,2)^2 + R(3,3)^2;
    droll_dR(8) = R(3,3) / denom_roll;  % d(roll)/dR(3,2)
    droll_dR(9) = -R(3,2) / denom_roll; % d(roll)/dR(3,3)

    % For pitch = -asin(R(3,1))
    dpitch_dR(7) = -1 / sqrt(1 - R(3,1)^2); % d(pitch)/dR(3,1)

    % For yaw = atan2(R(2,1), R(1,1))
    denom_yaw = R(2,1)^2 + R(1,1)^2;
    dyaw_dR(4) = R(1,1) / denom_yaw;  % d(yaw)/dR(2,1)
    dyaw_dR(1) = -R(2,1) / denom_yaw; % d(yaw)/dR(1,1)

    % Mapping from rotation matrix perturbation to the rotation error vector
    % For small perturbations, the rotation error vector δθ approximates the change in the rotation matrix as "R * skew(δθ)"
    % We need to compute the Jacobian of each element of R with respect to δθ

    % Construct the mapping from rotation error vector to rotation matrix element perturbation
    dR_dtheta = zeros(9, 3);

    % The structure of the perturbed rotation matrix: R * (I + skew(δθ))
    % where skew(δθ) is the skew-symmetric matrix of δθ

    % Compute each column of R * skew(e_i), where e_i is the standard basis vector
    for i = 1:3
        delta = zeros(3, 1);
        delta(i) = 1;

        % Construct the skew-symmetric matrix
        skew_delta = [0, -delta(3), delta(2); 
                       delta(3), 0, -delta(1); 
                       -delta(2), delta(1), 0];

        % Compute R * skew(delta) and flatten it
        R_skew = R * skew_delta;
        dR_dtheta(:, i) = R_skew(:);
    end

    % Combine the two Jacobian matrices
    J_euler_R = [droll_dR; dpitch_dR; dyaw_dR]; % 3x9
    J_R_theta = dR_dtheta;                      % 9x3

    % Compute the total Jacobian matrix
    J = J_euler_R * J_R_theta;  % 3x3

    % Apply the error propagation formula
    P_euler = J * P_rotation * J';
end