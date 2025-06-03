function P_euler = rotation_matrix_cov_to_euler_cov(P_rotation, R)
    % 输入:
    % P_rotation: ESKF中旋转误差向量的3x3协方差矩阵
    % R: 当前估计的旋转矩阵
    
    % 从旋转矩阵中提取欧拉角
    roll = atan2(R(3,2), R(3,3));
    pitch = -asin(R(3,1));
    yaw = atan2(R(2,1), R(1,1));
    
    % 计算从旋转误差向量到欧拉角变化的雅可比矩阵
    % 这个转换分两步：
    % 1. 从旋转误差向量到旋转矩阵扰动的映射
    % 2. 从旋转矩阵扰动到欧拉角扰动的映射
    
    % 计算欧拉角相对于旋转矩阵元素的偏导数
    droll_dR = zeros(1,9);
    dpitch_dR = zeros(1,9);
    dyaw_dR = zeros(1,9);
    
    % 对于roll = atan2(R(3,2), R(3,3))
    denom_roll = R(3,2)^2 + R(3,3)^2;
    droll_dR(8) = R(3,3) / denom_roll;  % d(roll)/dR(3,2)
    droll_dR(9) = -R(3,2) / denom_roll; % d(roll)/dR(3,3)
    
    % 对于pitch = -asin(R(3,1))
    dpitch_dR(7) = -1 / sqrt(1 - R(3,1)^2); % d(pitch)/dR(3,1)
    
    % 对于yaw = atan2(R(2,1), R(1,1))
    denom_yaw = R(2,1)^2 + R(1,1)^2;
    dyaw_dR(4) = R(1,1) / denom_yaw;  % d(yaw)/dR(2,1)
    dyaw_dR(1) = -R(2,1) / denom_yaw; % d(yaw)/dR(1,1)
    
    % 从旋转矩阵扰动到旋转误差向量的映射
    % 对于小扰动，旋转误差向量δθ可以通过"R * skew(δθ)"近似旋转矩阵的变化
    % 我们需要计算R的每个元素相对于δθ的雅可比
    
    % 构建从旋转误差向量到旋转矩阵元素扰动的映射
    dR_dtheta = zeros(9, 3);
    
    % 扰动旋转矩阵的结构: R * (I + skew(δθ))
    % 其中skew(δθ)是δθ的反对称矩阵
    
    % 计算R * skew(e_i)的每一列，其中e_i是标准基向量
    for i = 1:3
        delta = zeros(3, 1);
        delta(i) = 1;
        
        % 构建反对称矩阵
        skew_delta = [0, -delta(3), delta(2); 
                       delta(3), 0, -delta(1); 
                       -delta(2), delta(1), 0];
        
        % 计算R * skew(delta)并将其展平
        R_skew = R * skew_delta;
        dR_dtheta(:, i) = R_skew(:);
    end
    
    % 组合两个雅可比矩阵
    J_euler_R = [droll_dR; dpitch_dR; dyaw_dR]; % 3x9
    J_R_theta = dR_dtheta;                      % 9x3
    
    % 计算总的雅可比矩阵
    J = J_euler_R * J_R_theta;  % 3x3
    
    % 应用误差传播公式
    P_euler = J * P_rotation * J';
end