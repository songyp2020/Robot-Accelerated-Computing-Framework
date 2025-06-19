function geom = GeometryCalculation_casadi(rbt_df, static_model)

import casadi.*

if nargin < 2
    static_model = false;
end
geom.static_model = static_model;

tic;
geom = calc_geom_casadi(rbt_df, geom);
disp(['GeometryCalculation_casadi: time ' num2str(toc) ' sec']);

end

function geom = calc_geom_casadi(rbt_df, geom)
import casadi.*
f = utils_casadi();

geom.T_0n = cell(1, rbt_df.frame_num);
geom.R = cell(1, rbt_df.frame_num);
geom.p_n = cell(1, rbt_df.frame_num);
geom.T_0nc = cell(1, rbt_df.frame_num);
geom.p_c = cell(1, rbt_df.frame_num);
geom.v_cw = cell(1, rbt_df.frame_num);
geom.w_b = cell(1, rbt_df.frame_num);

Tsb = cell(1, rbt_df.frame_num);

if geom.static_model
    disp('static_model = true, 忽略动力学部分（速度/加速度）');
else
    disp('(default) static_model = false, 计算完整运动学');
end

for num = 1:rbt_df.frame_num
    disp(['计算坐标系 ' num2str(rbt_df.link_nums{num}) ' 的运动学']);
    
    if num == 1
        geom.T_0n{num} = rbt_df.M{num};
        Tsb{num} = SX.eye(4);
        geom.R{num} = geom.T_0n{num}(1:3, 1:3);
        geom.p_n{num} = geom.T_0n{num}(1:3, 4);
        
        geom.T_0nc{num} = geom.T_0n{num};
        geom.p_c{num} = geom.p_n{num};
        geom.v_cw{num} = SX.zeros(3, 1);
        geom.w_b{num} = SX.zeros(3, 1);
        continue
    end
    
    prev_idx = rbt_df.prev_link_num{num} + 1;
    
    exp_S_theta = f.matrix_exp_6(rbt_df.screw{num}, rbt_df.theta{num});
    
    Tsb{num} = Tsb{prev_idx} * exp_S_theta;
    
    geom.T_0n{num} = Tsb{num} * rbt_df.M{num};
    
    geom.R{num} = geom.T_0n{num}(1:3, 1:3);
    geom.p_n{num} = geom.T_0n{num}(1:3, 4);
    
    if rbt_df.use_inertia{num}
        T_com = f.translation_transfmat(rbt_df.r_by_ml{num});
        geom.T_0nc{num} = geom.T_0n{num} * T_com;
        geom.p_c{num} = geom.T_0nc{num}(1:3, 4);
    else
        geom.T_0nc{num} = geom.T_0n{num};
        geom.p_c{num} = geom.p_n{num};
    end
    
    % varsStruct = struct('screw', rbt_df.screw{num}, ...
    %     'theta', rbt_df.theta{num}, ...
    %     'exp_S_theta', exp_S_theta, ...
    %     'Tsb', Tsb{num}, ...
    %     'M', rbt_df.M{num}, ...
    %     'T_0n', geom.T_0n{num}, ...
    %     'R', geom.R{num}, ...
    %     'p_n', geom.p_n{num}, ...
    %     'T_0nc', geom.T_0nc{num}, ...
    %     'p_c', geom.p_c{num});
    
    % f.debug_vars(rbt_df, varsStruct, num);
    
    if ~geom.static_model
        J_c = jacobian(geom.p_c{num}, vertcat(rbt_df.coordinates{:}));
        geom.v_cw{num} = J_c * vertcat(rbt_df.d_coordinates{:});
        R_mat = geom.R{num};
        dR_vec = jacobian(R_mat(:), vertcat(rbt_df.coordinates{:})) * vertcat(rbt_df.d_coordinates{:});
        dR = reshape(dR_vec, 3, 3);
        S = R_mat' * dR;
        w_b = SX([S(3,2); S(1,3); S(2,1)]);
        geom.w_b{num} = simplify(w_b);
    else
        geom.v_cw{num} = SX.zeros(3, 1);
        geom.w_b{num} = SX.zeros(3, 1);
    end
end

disp('运动学计算完成');
if ~geom.static_model
    disp('- 包含速度运动学（v_cw, w_b）');
end

end