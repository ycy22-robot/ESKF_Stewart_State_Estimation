function p = FwKineQ(lg, q0)
    % tic
    
    %====================
    % 参数预设
    %====================
    maxIter = 20000;      % 最大迭代次数，可根据需要调整
    tol     = 1e-4;      % 收敛判据
    
    % 如果未传入 q0，则使用默认初值
    if nargin < 2 || isempty(q0)
        P0 = [0; 0; 1; 0; 0; 0];  % 默认初值
    else
        P0 = q0(:);
    end
    
    % leg长度必须是列向量
    lg = lg(:);
    
    % 预先分配迭代变量空间
    P = zeros(6, maxIter+1);
    P(:,1) = P0;
    
    dl = Inf;  % 用于存储支链长度误差的范数
    i  = 1;    % 迭代计数器
    
    %====================
    % 铰链点坐标预计算
    %====================
    [a_A, b_B_local] = precomputeHingePoints();
    
    %====================
    % 主循环
    %====================
    for k = 2:maxIter
        
        %--- Step #3：构造 T 矩阵 ---
        a = P(4,k-1) * pi/180;
        b = P(5,k-1) * pi/180;
        
        B = [1,       0,      sin(b);
             0,  cos(a), -sin(a)*cos(b);
             0,  sin(a),  cos(a)*cos(b)];
        
        T = [eye(3),       zeros(3,3);
             zeros(3,3),   B        ];
         
        %--- Step #4：逆运动学，得到当前各连杆长度 ---
        [l, n, R] = mInv(P(:,k-1), a_A, b_B_local);

        %--- Step #2：计算雅可比 ---
        J = jaco(n, R, b_B_local);  
        
        % 计算长度误差
        Dl = lg - l';  
        dl = norm(Dl, 2);  % 2-范数
        
        %--- 判断收敛 ---
        if dl < tol
            break;
        end
        
        %--- Step #5：更新 P ---
        P(:,k) = P(:,k-1) + (J*T) \ Dl;
        i = i + 1;  % 迭代计数
    end
    
    % 若迭代结束或满足收敛，则取最后一次迭代结果
    p = [P(1:3,i); deg2rad(P(4:6,i))];
    % toc
end

%=======================================================================
% 逆运动学辅助函数
%=======================================================================
function [L, n, R] = mInv(X, a_A, b_B_local)
    p = [X(1); X(2); X(3)];
    R = Eular(deg2rad(X(4:6))); % 姿态变换：欧拉角转旋转矩阵
    d = p + R*b_B_local - a_A;  % 支链向量
    L = sqrt(sum(d.^2, 1));     % 各列求范数，得到长度
    n = d ./ L;                 % 支链方向向量
    % b_B = R*b_B_local;          % 动平台铰链点在静平台坐标系中的位置
end

%=======================================================================
% 雅可比矩阵计算
%=======================================================================
function J = jaco(n, R, b_B)
    % 使用 mInv 计算所需的 n 和 R
    % [~, n, R] = mInv(P, a_A, b_B);

    % 提前计算 R * b_B
    R_s = R * b_B; % 结果是 3x6 矩阵
    
    % 手动计算交叉乘积，逐列操作展开
    cross_terms = [
        R_s(2, :) .* n(3, :) - R_s(3, :) .* n(2, :);
        R_s(3, :) .* n(1, :) - R_s(1, :) .* n(3, :);
        R_s(1, :) .* n(2, :) - R_s(2, :) .* n(1, :)
    ];
    
    % 拼接，得到 6x6 雅可比矩阵
    J = [n', cross_terms'];
end

%=======================================================================
% 铰链点坐标预计算
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

    phi_a = phi_a_0;            % 静平台
    phi_b = phi_b_0 - 60;       % 动平台
    
    a_A = [R_a*cosd(phi_a)'; R_a*sind(phi_a)'; zeros(1,6)];
    b_B_local = [R_b*cosd(phi_b)'; R_b*sind(phi_b)'; zeros(1,6)];
end

%=======================================================================
% 欧拉角→旋转矩阵
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