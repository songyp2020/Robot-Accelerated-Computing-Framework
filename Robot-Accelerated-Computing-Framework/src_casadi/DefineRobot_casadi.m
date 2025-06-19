function rbt = DefineRobot_casadi(model_name, params)

import casadi.*
rbt = struct();
rbt.name = model_name;
rbt.frame_num = height(params);         
rbt.link_nums = params(:, 1);           
rbt.prev_link_num = params(:, 2);       
rbt.succ_link_num = params(:, 3);       
rbt.shat = params(:, 4);                
rbt.rs = params(:, 5);                  
rbt.M_R = params(:, 6);                 
rbt.theta = params(:, 7);               
rbt.h = params(:, 8);                   
rbt.use_inertia = params(:, 9);         
rbt.spring_dl = params(:, 10);       

rbt = convert_to_casadi(rbt);

rbt = DefineRobot_gentransfm_casadi(rbt);

rbt = DefineRobot_genparams_casadi(rbt);

rbt = gen_coordinates_casadi(rbt);

end

function rbt = convert_to_casadi(rbt)
import casadi.*
all_syms = [];
for i = 1:length(rbt.theta)
    if isa(rbt.theta{i}, 'sym')
        all_syms = [all_syms, symvar(rbt.theta{i})];
    end
end
all_syms = unique(all_syms);
rbt.sym_to_casadi = containers.Map();
for i = 1:length(all_syms)
    var_name = char(all_syms(i));
    casadi_var = SX.sym(var_name);
    rbt.sym_to_casadi(var_name) = casadi_var;
end


for i = 1:length(rbt.theta)
    if isa(rbt.theta{i}, 'sym')
        rbt.theta{i} = convert_sym_expr(rbt.theta{i}, rbt.sym_to_casadi);
    elseif isnumeric(rbt.theta{i})
        rbt.theta{i} = SX(rbt.theta{i});
    end
end


for i = 1:rbt.frame_num
    if isnumeric(rbt.shat{i})
        rbt.shat{i} = SX(rbt.shat{i});
    end
    if isnumeric(rbt.rs{i})
        rbt.rs{i} = SX(rbt.rs{i});
    end
    
    if isnumeric(rbt.M_R{i})
        rbt.M_R{i} = SX(rbt.M_R{i});
    end
end
for i = 1:length(rbt.spring_dl)
    sdl = rbt.spring_dl{i};
    if ~(isnumeric(sdl) && all(isnan(sdl)))
        if isnumeric(sdl)
            rbt.spring_dl{i} = SX(sdl);
        elseif isa(sdl, 'sym')
            vec = SX.zeros(0, 1);
            for k = 1:numel(sdl)
                expr = convert_sym_expr(sdl(k), rbt.sym_to_casadi);
                vec = [vec; expr];
            end
            rbt.spring_dl{i} = vec;
        end
    end
end
end

function casadi_expr = convert_sym_expr(sym_expr, sym_to_casadi)
import casadi.*

if isnumeric(sym_expr)
    casadi_expr = SX(sym_expr);
    return;
end

vars = symvar(sym_expr);
if isempty(vars)
    casadi_expr = SX(double(sym_expr));
    return;
end

sym_vars = [];
casadi_vars = [];
for i = 1:length(vars)
    var_name = char(vars(i));
    if isKey(sym_to_casadi, var_name)
        sym_vars = [sym_vars, vars(i)];
        casadi_vars = [casadi_vars; sym_to_casadi(var_name)];
    end
end

expr_str = char(sym_expr);

for i = 1:length(sym_vars)
    var_name = char(sym_vars(i));
    expr_str = regexprep(expr_str, ['\<' var_name '\>'], ['sym_to_casadi(''' var_name ''')']);
end

expr_str = regexprep(expr_str, '\<pi\>', num2str(pi));

try
    casadi_expr = eval(expr_str);
catch ME
    warning(['转换表达式失败: ' expr_str]);
    casadi_expr = SX(0);
end
end

function rbt = DefineRobot_gentransfm_casadi(rbt)
import casadi.*

f = utils_casadi();

M = cell(1, rbt.frame_num);
ws = cell(1, rbt.frame_num);
vs = cell(1, rbt.frame_num);
screw = cell(1, rbt.frame_num);
rbt.joint_type = [];

for num = 1:rbt.frame_num
    if isinf(rbt.h{num})
        ws{num} = SX.zeros(3, 1);
        h_temp = 1;
    else
        ws{num} = rbt.shat{num};
        h_temp = rbt.h{num};
    end
    
    vs{num} = -cross(ws{num}, rbt.rs{num}) + h_temp*rbt.shat{num};
    
    screw{num} = [ws{num}; vs{num}];
    
    M{num} = SX.zeros(4, 4);
    M{num}(1:3, 1:3) = rbt.M_R{num};
    M{num}(1:3, 4) = rbt.rs{num};
    M{num}(4, 4) = 1;
    
    screw_norm = norm(screw{num});
    if isa(screw_norm, 'casadi.SX')
        if isinf(rbt.h{num})
            rbt.joint_type = [rbt.joint_type, "P"];
        elseif rbt.h{num} == 0
            rbt.joint_type = [rbt.joint_type, "R"];
        else
            rbt.joint_type = [rbt.joint_type, "A"];
        end
    else
        if double(screw_norm) < eps
            rbt.joint_type = [rbt.joint_type, "F"];
        elseif isinf(rbt.h{num})
            rbt.joint_type = [rbt.joint_type, "P"];
        elseif rbt.h{num} == 0
            rbt.joint_type = [rbt.joint_type, "R"];
        else
            rbt.joint_type = [rbt.joint_type, "A"];
        end
    end
end

rbt.M = M;
rbt.ws = ws;
rbt.vs = vs;
rbt.screw = screw;
end

function rbt = DefineRobot_genparams_casadi(rbt)
import casadi.*

f = utils_casadi();
rbt.m = cell(1, rbt.frame_num);              
rbt.l = cell(1, rbt.frame_num);              
rbt.r = cell(1, rbt.frame_num);              
rbt.r_by_ml = cell(1, rbt.frame_num);       
rbt.L_vec = cell(1, rbt.frame_num);         
rbt.I_vec = cell(1, rbt.frame_num);        
rbt.L_mat = cell(1, rbt.frame_num);         
rbt.I_mat = cell(1, rbt.frame_num);         
rbt.I_by_Llm = cell(1, rbt.frame_num);       
rbt.spring_formula = cell(1, rbt.frame_num); 
rbt.spring_param = cell(1, rbt.frame_num);   

rbt.std_params = [];    
rbt.bary_params = [];   

spring_num = 0;

for num = 1:(rbt.frame_num-1)
    rbt.m{num+1} = SX.sym(['m' num2str(num)]);
    
    rbt.l{num+1} = SX.sym(['l' num2str(num)], 3, 1);
    
    rbt.r{num+1} = SX.sym(['r' num2str(num)], 3, 1);
    
    rbt.I_vec{num+1} = SX.sym(['I' num2str(num)], 6, 1);
    
    rbt.L_vec{num+1} = SX.sym(['L' num2str(num)], 6, 1);
    
    rbt.I_mat{num+1} = f.inertia_vec2tensor(rbt.I_vec{num+1});
    rbt.L_mat{num+1} = f.inertia_vec2tensor(rbt.L_vec{num+1});
    
    rbt.r_by_ml{num+1} = f.ml2r(rbt.m{num+1}, rbt.l{num+1});
    
    rbt.I_by_Llm{num+1} = f.Lmr2I(rbt.L_mat{num+1}, rbt.m{num+1}, rbt.r_by_ml{num+1});
    if ~(isnumeric(rbt.spring_dl{num+1}) && all(isnan(rbt.spring_dl{num+1})))
        spring_param_number = width(rbt.spring_dl{num+1});
        rbt.spring_param{num+1} = SX.sym(['K' num2str(num)], spring_param_number, 1);
        rbt.spring_formula{num+1} = 0;
        for j = 1:spring_param_number
            rbt.spring_formula{num+1} = rbt.spring_formula{num+1} + ...
                rbt.spring_dl{num+1}(j) * rbt.spring_param{num+1}(j);
        end
        spring_num = spring_num + spring_param_number;
    end
end
rbt.m{1} = SX(0);
rbt.l{1} = SX.zeros(3, 1);
rbt.r{1} = SX.zeros(3, 1);
rbt.r_by_ml{1} = SX.zeros(3, 1);
rbt.L_vec{1} = SX.zeros(6, 1);
rbt.I_vec{1} = SX.zeros(6, 1);
rbt.L_mat{1} = SX.zeros(3, 3);
rbt.I_mat{1} = SX.zeros(3, 3);
rbt.I_by_Llm{1} = SX.zeros(3, 3);

for num = 1:(rbt.frame_num-1)
    if rbt.use_inertia{num+1}
        rbt.bary_params = [rbt.bary_params; rbt.L_vec{num+1}; rbt.l{num+1}; rbt.m{num+1}];
        
        rbt.std_params = [rbt.std_params; rbt.I_vec{num+1}; rbt.r{num+1}; rbt.m{num+1}];
    end
end

rbt.inertial_param_number = length(rbt.bary_params);

for num = 1:(rbt.frame_num-1)
    if ~(isnumeric(rbt.spring_dl{num+1}) && all(isnan(rbt.spring_dl{num+1})))
        rbt.bary_params = [rbt.bary_params; rbt.spring_param{num+1}];
        rbt.std_params = [rbt.std_params; rbt.spring_param{num+1}];
    end
end

rbt.spring_param_number = spring_num;
end

function rbt = gen_coordinates_casadi(rbt)
import casadi.*

rbt.coordinates = {};
rbt.coordinates_joint_type = [];

for num = 1:length(rbt.link_nums)
    theta_vars = extract_casadi_variables(rbt.theta{num});
    
    for i = 1:length(theta_vars)
        s = theta_vars{i};
        if ~is_variable_in_list(s, rbt.coordinates)
            rbt.coordinates{end+1} = s;
            rbt.coordinates_joint_type = [rbt.coordinates_joint_type, rbt.joint_type(num)];
        end
    end
end
rbt.dof = length(rbt.coordinates);

rbt.d_coordinates = {};
rbt.dd_coordinates = {};
rbt.coordinates_t = {};

t = SX.sym('t');

for i = 1:length(rbt.coordinates)
    co = rbt.coordinates{i};
    
    var_name = get_casadi_variable_name(co, i);
    
    rbt.d_coordinates{end+1} = SX.sym(['d' var_name]);
    rbt.dd_coordinates{end+1} = SX.sym(['dd' var_name]);
    
    rbt.coordinates_t{end+1} = SX.sym([var_name 't']);
end

rbt.d_coordinates_t = {};
rbt.dd_coordinates_t = {};

for i = 1:length(rbt.coordinates_t)
    co_t = rbt.coordinates_t{i};
    
    var_name = get_casadi_variable_name(co_t, i);
    if contains(var_name, 't')
        base_name = strrep(var_name, 't', '');
    else
        base_name = var_name;
    end
    
    rbt.d_coordinates_t{end+1} = SX.sym(['d' base_name 't']);
    rbt.dd_coordinates_t{end+1} = SX.sym(['dd' base_name 't']);
end
rbt.subs_q2qt = create_substitution_pair(rbt.coordinates, rbt.coordinates_t);
rbt.subs_dq2dqt = create_substitution_pair(rbt.d_coordinates, rbt.d_coordinates_t);
rbt.subs_ddq2ddqt = create_substitution_pair(rbt.dd_coordinates, rbt.dd_coordinates_t);
rbt.subs_qt2q = create_substitution_pair(rbt.coordinates_t, rbt.coordinates);
rbt.subs_dqt2dq = create_substitution_pair(rbt.d_coordinates_t, rbt.d_coordinates);
rbt.subs_ddqt2ddq = create_substitution_pair(rbt.dd_coordinates_t, rbt.dd_coordinates);

rbt.q_for_frame = cell(1, rbt.frame_num);
rbt.dq_for_frame = cell(1, rbt.frame_num);
rbt.ddq_for_frame = cell(1, rbt.frame_num);

for i = 1:rbt.frame_num
    if rbt.joint_type(i) == "P" || rbt.joint_type(i) == "R" || rbt.joint_type(i) == "A"
        q = rbt.theta{i};
        rbt.q_for_frame{i} = q;
        
        qt = casadi_substitute(q, rbt.subs_q2qt.from, rbt.subs_q2qt.to);
        
        if ~isempty(rbt.coordinates)
            coord_idx = find_coordinate_index(q, rbt.coordinates);
            if ~isempty(coord_idx)
                rbt.dq_for_frame{i} = rbt.d_coordinates{coord_idx};
                rbt.ddq_for_frame{i} = rbt.dd_coordinates{coord_idx};
            else
                rbt.dq_for_frame{i} = SX.sym(['dq_frame_' num2str(i)]);
                rbt.ddq_for_frame{i} = SX.sym(['ddq_frame_' num2str(i)]);
            end
        end
    else
        rbt.q_for_frame{i} = [];
        rbt.dq_for_frame{i} = [];
        rbt.ddq_for_frame{i} = [];
    end
end

end

function vars = extract_casadi_variables(expr)
import casadi.*
vars = {};

if isa(expr, 'casadi.SX')
    try
        sym_vars = symvar(expr);
        if iscell(sym_vars)
            vars = sym_vars;
        else
            for i = 1:length(sym_vars)
                vars{end+1} = sym_vars(i);
            end
        end
    catch
        vars = {expr};
    end
elseif isa(expr, 'casadi.MX')
    try
        sym_vars = symvar(expr);
        if iscell(sym_vars)
            vars = sym_vars;
        else
            for i = 1:length(sym_vars)
                vars{end+1} = sym_vars(i);
            end
        end
    catch
        vars = {expr};
    end
else
    try
        if ischar(expr) || isstring(expr)
            vars = {SX.sym(char(expr))};
        else
            vars = {expr};
        end
    catch
        vars = {};
    end
end
end

function var_name = get_casadi_variable_name(var, fallback_idx)
var_name = '';

if isa(var, 'casadi.SX') || isa(var, 'casadi.MX')
    try
        var_name = var.name();
        if isempty(var_name)
            var_name = ['q' num2str(fallback_idx)];
        end
    catch
        var_name = ['q' num2str(fallback_idx)];
    end
elseif ischar(var) || isstring(var)
    var_name = char(var);
else
    var_name = ['q' num2str(fallback_idx)];
end
end

function result = is_variable_in_list(var, var_list)
result = false;
if isempty(var_list)
    return;
end

var_name = get_casadi_variable_name(var, 0);

for i = 1:length(var_list)
    existing_var = var_list{i};
    existing_name = get_casadi_variable_name(existing_var, i);
    
    if strcmp(var_name, existing_name)
        result = true;
        return;
    end
end
end

function subs_pair = create_substitution_pair(from_list, to_list)
subs_pair = struct();
subs_pair.from = from_list;
subs_pair.to = to_list;
end

function result = casadi_substitute(expr, from_vars, to_vars)
import casadi.*
if isempty(from_vars) || isempty(to_vars)
    result = expr;
    return;
end

try
    if iscell(from_vars) && ~isempty(from_vars)
        from_vec = vertcat(from_vars{:});
    else
        from_vec = from_vars;
    end
    
    if iscell(to_vars) && ~isempty(to_vars)
        to_vec = vertcat(to_vars{:});
    else
        to_vec = to_vars;
    end
    
    result = substitute(expr, from_vec, to_vec);
catch ME
    warning(ME.identifier, 'CasADi substitution failed: %s', ME.message);
    result = expr;
end
end

function idx = find_coordinate_index(q, coordinates)
idx = [];
if isempty(coordinates)
    return;
end

q_name = get_casadi_variable_name(q, 0);

for i = 1:length(coordinates)
    coord = coordinates{i};
    coord_name = get_casadi_variable_name(coord, i);
    
    if strcmp(q_name, coord_name)
        idx = i;
        return;
    end
end
end

