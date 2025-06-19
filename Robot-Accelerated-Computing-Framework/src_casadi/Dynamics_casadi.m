function dyn = Dynamics_casadi(rbt_df, geom, gravity_vec)

import casadi.*

dyn = struct();
gravity_vec = gravity_vec(:);
dyn.gravity_vec = gravity_vec;

tic;
fprintf('\n');
fprintf('========================================\n');
fprintf('         CasADi 动力学计算开始\n');
fprintf('========================================\n');

if geom.static_model
    fprintf('⚙️  模式: 静态模型 (仅重力项)\n');
else
    fprintf('⚙️  模式: 完整动力学模型\n');
end

fprintf('📊 系统信息: %d个连杆, %d个自由度\n', rbt_df.frame_num-1, rbt_df.dof);

%% 初始化能量项
p_e = SX(0);
k_e = SX(0);

q_vec = vertcat(rbt_df.coordinates{:});
dq_vec = vertcat(rbt_df.d_coordinates{:});
ddq_vec = vertcat(rbt_df.dd_coordinates{:});

%% 计算每个连杆的能量
fprintf('\n');
fprintf('🔋 计算连杆能量:\n');
for num = 2:rbt_df.frame_num
    k_e_n = SX(0);
    if rbt_df.use_inertia{num}
        fprintf('   连杆 %d: ✓\n', num-1);
        
        p_e = p_e - rbt_df.m{num} * (gravity_vec' * geom.p_c{num});
        
        if ~geom.static_model
            K_trans = rbt_df.m{num} * geom.v_cw{num}' * (geom.v_cw{num}) / 2;
            K_rot = geom.w_b{num}' * rbt_df.I_by_Llm{num} * geom.w_b{num} / 2;
            k_e_n = K_trans + K_rot;
        end
    else
        fprintf('   连杆 %d: ✗ (跳过) \n', num-1);
    end
    k_e = k_e + k_e_n;
end

dyn.Lagrange.K = k_e;
dyn.Lagrange.P = p_e;
L = k_e - p_e;        
dyn.Lagrange.L = L;

%% 计算关节力矩（拉格朗日方程）
fprintf('\n');
fprintf('🔧 计算关节力矩 (拉格朗日方程):\n');
dyn.tau = SX.zeros(rbt_df.dof, 1);

for i = 1:rbt_df.dof
    fprintf('   关节 %d: 计算中...', i);
    q_i = rbt_df.coordinates{i};
    dq_i = rbt_df.d_coordinates{i};
    dK_ddq = jacobian(k_e, dq_i);
    d_dK_ddq_dt = jacobian(dK_ddq, q_vec) * dq_vec + jacobian(dK_ddq, dq_vec) * ddq_vec;
    dL_dq = jacobian(L, q_i);
    dyn.tau(i) = d_dK_ddq_dt - dL_dq;
    
    fprintf(' ✓\n');
end
fprintf('   关节力矩计算完成 ✓\n');

M = jacobian(jacobian(k_e, dq_vec), dq_vec);
G = -jacobian(p_e, q_vec);
C_times_dq = dyn.tau - M * ddq_vec - G';
dyn.M = M;
dyn.G = G;
dyn.C_times_dq = C_times_dq;

dyn.tau_MCG = dyn.tau;

%% 添加弹簧力（如果有）
fprintf('\n');
fprintf('🌸 检查弹簧力:\n');
dyn.tau_spring = SX.zeros(rbt_df.dof, 1);
spring_count = 0;

for i = 1:rbt_df.frame_num
    if isfield(rbt_df, 'spring_dl') && i <= length(rbt_df.spring_dl)
        if ~(isnumeric(rbt_df.spring_dl{i}) && all(isnan(rbt_df.spring_dl{i})))
            spring_count = spring_count + 1;
            if isfield(rbt_df, 'spring_formula') && i <= length(rbt_df.spring_formula)
                tau_s = rbt_df.spring_formula{i};
                
                if isfield(rbt_df, 'q_for_frame') && i <= length(rbt_df.q_for_frame) && ~isempty(rbt_df.q_for_frame{i})
                    q_i = rbt_df.q_for_frame{i};
                    
                    for j = 1:rbt_df.dof
                        partial_q = jacobian(q_i, rbt_df.coordinates{j});
                        if ~is_zero(partial_q)
                            dyn.tau_spring(j) = dyn.tau_spring(j) + partial_q * tau_s;
                        end
                    end
                end
            end
        end
    end
end

if spring_count > 0
    fprintf('   发现 %d 个弹簧元素 ✓\n', spring_count);
else
    fprintf('   无弹簧元素 ✓\n');
end

dyn.tau = dyn.tau_MCG + dyn.tau_spring;

dynamics_time = toc;
fprintf('\n');
fprintf('⏱️  动力学计算耗时: %.4f 秒\n', dynamics_time);

%% 计算回归矩阵
dyn = Dynamics_calc_regressor_casadi(rbt_df, dyn);

fprintf('\n');
fprintf('========================================\n');
fprintf('         CasADi 动力学计算完成\n');
fprintf('========================================\n');
fprintf('\n');

end

function dyn = Dynamics_calc_regressor_casadi(rbt_df, dyn)

import casadi.*
tic;
fprintf('🔍 计算回归矩阵:\n');

if isfield(rbt_df, 'bary_params') && ~isempty(rbt_df.bary_params)
    param_vec = rbt_df.bary_params(:);
    param_count = length(param_vec);
    fprintf('   参数个数: %d\n', param_count);
    
    H_symbolic = jacobian(dyn.tau, param_vec);
    dyn.H_symbolic = H_symbolic;
    dyn.inertial_param = param_vec;
    
    fprintf('   回归矩阵维度: %dx%d ✓\n', size(H_symbolic, 1), size(H_symbolic, 2));
else
    fprintf('   ⚠️  未找到参数，跳过回归矩阵计算\n');
    dyn.H_symbolic = [];
    dyn.inertial_param = [];
end

regressor_time = toc;
fprintf('⏱️  回归矩阵计算耗时: %.4f 秒\n', regressor_time);

end