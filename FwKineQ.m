function p = FwKineQ(lg, q0)
    % tic

    %====================
    % Parameter settings
    %====================
    maxIter = 20000;      % Maximum number of iterations, adjust as needed
    tol     = 1e-4;       % Convergence criterion

    % If q0 is not provided, use default initial value
    if nargin < 2 || isempty(q0)
        P0 = [0; 0; 1; 0; 0; 0];  % Default initial value
    else
        P0 = q0(:);
    end

    % Leg length must be a column vector
    lg = lg(:);

    % Pre-allocate space for iteration variables
    P = zeros(6, maxIter+1);
    P(:,1) = P0;

    dl = Inf;   % To store the norm of the rod length error
    i  = 1;     % Iteration counter

    %====================
    % Precompute hinge point coordinates
    %====================
    [a_A, b_B_local] = precomputeHingePoints();

    %====================
    % Main loop
    %====================
    for k = 2:maxIter

        %--- Step #3: Construct T matrix ---
        a = P(4,k-1) * pi/180;
        b = P(5,k-1) * pi/180;

        B = [1,       0,      sin(b);
             0,  cos(a), -sin(a)*cos(b);
             0,  sin(a),  cos(a)*cos(b)];

        T = [eye(3),       zeros(3,3);
             zeros(3,3),   B        ];

        %--- Step #4: Inverse kinematics, get current rod lengths ---
        [l, n, R] = mInv(P(:,k-1), a_A, b_B_local);

        %--- Step #2: Compute Jacobian ---
        J = jaco(n, R, b_B_local);

        % Compute length error
        Dl = lg - l';
        dl = norm(Dl, 2);  % 2-norm

        %--- Check convergence ---
        if dl < tol
            break;
        end

        %--- Step #5: Update P ---
        P(:,k) = P(:,k-1) + (J*T) \ Dl;
        i = i + 1;  % Iteration count
    end

    % If iteration ends or convergence is reached, use the last result
    p = [P(1:3,i); deg2rad(P(4:6,i))];
    % toc
end

%=======================================================================
% Inverse kinematics helper function
%=======================================================================
function [L, n, R] = mInv(X, a_A, b_B_local)
    p = [X(1); X(2); X(3)];
    R = Eular(deg2rad(X(4:6))); % Pose transformation: Euler angles to rotation matrix
    d = p + R*b_B_local - a_A;  % Rod vectors
    L = sqrt(sum(d.^2, 1));     % Norm for each column, get lengths
    n = d ./ L;                 % Rod direction vectors
    % b_B = R*b_B_local;        % End platform hinge points in base coordinates (not used)
end

%=======================================================================
% Jacobian computation
%=======================================================================
function J = jaco(n, R, b_B)
    % Use mInv to compute needed n and R
    % [~, n, R] = mInv(P, a_A, b_B);

    % Precompute R * b_B
    R_s = R * b_B; % Result is 3x6 matrix

    % Manually compute cross product, operate column-wise
    cross_terms = [
        R_s(2, :) .* n(3, :) - R_s(3, :) .* n(2, :);
        R_s(3, :) .* n(1, :) - R_s(1, :) .* n(3, :);
        R_s(1, :) .* n(2, :) - R_s(2, :) .* n(1, :)
    ];

    % Concatenate to get 6x6 Jacobian matrix
    J = [n', cross_terms'];
end

%=======================================================================
% Precompute hinge point coordinates
%=======================================================================
function [a_A, b_B_local] = precomputeHingePoints()
    R_a     = 800e-3;
    theta_a = 18/2;
    R_b     = 500e-3;
    theta_b = 10/2;

    phi_a_0 = [-theta_a; theta_a; -theta_a+120; ...
               theta_a+120; -theta_a+240; theta_a+240];
    phi_b_0 = [theta_b; -theta_b+120; theta_b+120; ...
               -theta_b+240; theta_b+240; -theta_b+360];

    phi_a = phi_a_0;            % Base platform
    phi_b = phi_b_0 - 60;       % End platform

    a_A = [R_a*cosd(phi_a)'; R_a*sind(phi_a)'; zeros(1,6)];
    b_B_local = [R_b*cosd(phi_b)'; R_b*sind(phi_b)'; zeros(1,6)];
end

%=======================================================================
% Euler angles → Rotation matrix
%=======================================================================
function R = Eular(rpy)
    roll  = rpy(1);
    pitch = rpy(2);
    yaw   = rpy(3);
    cr = cos(roll);  sr = sin(roll);
    cp = cos(pitch); sp = sin(pitch);
    cy = cos(yaw);   sy = sin(yaw);
    R = [ cy*cp, cy*sp*sr - sy*cr, cy*sp*cr + sy*sr;
          sy*cp, sy*sp*sr + cy*cr, sy*sp*cr - cy*sr;
          -sp,   cp*sr,            cp*cr ];
end