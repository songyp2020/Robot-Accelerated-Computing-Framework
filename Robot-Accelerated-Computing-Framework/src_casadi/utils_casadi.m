function result = utils_casadi

import casadi.*

result.vec2so3 = @vec2so3_casadi;
result.so32vec = @so32vec_casadi;
result.vec2se3 = @vec2se3_casadi;
result.se32vec = @se32vec_casadi;
result.matrix_exp_3 = @matrix_exp_3_casadi;
result.matrix_exp_6 = @matrix_exp_6_casadi;
result.inertia_vec2tensor = @inertia_vec2tensor_casadi;
result.inertia_tensor2vec = @inertia_tensor2vec_casadi;
result.translation_transfmat = @translation_transfmat_casadi;
result.ml2r = @ml2r_casadi;
result.Lmr2I = @Lmr2I_casadi;
result.debug_p_c = @debug_p_c;
result.debug_vars = @debug_vars;
result.display_casadi_symbols = @display_casadi_symbols;
end

function matrix = vec2so3_casadi(vector)
import casadi.*

if ~isa(vector, 'casadi.SX') && ~isa(vector, 'casadi.MX') && ~isa(vector, 'casadi.DM')
    vector = SX(vector);
end

matrix = SX.zeros(3, 3);
matrix(1,2) = -vector(3); matrix(1,3) = vector(2);
matrix(2,1) = vector(3);  matrix(2,3) = -vector(1);
matrix(3,1) = -vector(2); matrix(3,2) = vector(1);
end

function vector = so32vec_casadi(matrix)
import casadi.*

vector = SX.zeros(3, 1);
vector(1) = matrix(3, 2);
vector(2) = matrix(1, 3);
vector(3) = matrix(2, 1);
end

function matrix = vec2se3_casadi(vector)
import casadi.*

if ~isa(vector, 'casadi.SX') && ~isa(vector, 'casadi.MX') && ~isa(vector, 'casadi.DM')
    vector = SX(vector);
end

matrix = SX.zeros(4, 4);
matrix(1:3, 1:3) = vec2so3_casadi(vector(1:3));
matrix(1:3, 4) = vector(4:6);
end

function vector = se32vec_casadi(matrix)
import casadi.*

vector = SX.zeros(6, 1);
vector(1:3) = so32vec_casadi(matrix(1:3, 1:3));
vector(4:6) = matrix(1:3, 4);
end

function R = matrix_exp_3_casadi(omghat, theta)
import casadi.*

if ~isa(omghat, 'casadi.SX') && ~isa(omghat, 'casadi.MX') && ~isa(omghat, 'casadi.DM')
    omghat = SX(omghat);
end
if ~isa(theta, 'casadi.SX') && ~isa(theta, 'casadi.MX') && ~isa(theta, 'casadi.DM')
    theta = SX(theta);
end

omg_mat = vec2so3_casadi(omghat);
R = SX.eye(3) + sin(theta)*omg_mat + (1-cos(theta))*omg_mat*omg_mat;
end

function T = matrix_exp_6_casadi(screw, theta)
import casadi.*

if ~isa(screw, 'casadi.SX') && ~isa(screw, 'casadi.MX') && ~isa(screw, 'casadi.DM')
    screw = SX(screw);
end
if ~isa(theta, 'casadi.SX') && ~isa(theta, 'casadi.MX') && ~isa(theta, 'casadi.DM')
    theta = SX(theta);
end

omghat = screw(1:3);
screw_mat = vec2se3_casadi(screw);
omg_mat = screw_mat(1:3, 1:3);

T = SX.zeros(4, 4);
T(1:3, 1:3) = matrix_exp_3_casadi(omghat, theta);
T(1:3, 4) = (SX.eye(3)*theta + (1-cos(theta))*omg_mat + (theta-sin(theta))*omg_mat*omg_mat)*screw_mat(1:3, 4);
T(4, 4) = 1;
end

function I = inertia_vec2tensor_casadi(vector)
import casadi.*

if ~isa(vector, 'casadi.SX') && ~isa(vector, 'casadi.MX') && ~isa(vector, 'casadi.DM')
    vector = SX(vector);
end

I = SX.zeros(3, 3);
I(1,1) = vector(1); I(1,2) = vector(2); I(1,3) = vector(3);
I(2,1) = vector(2); I(2,2) = vector(4); I(2,3) = vector(5);
I(3,1) = vector(3); I(3,2) = vector(5); I(3,3) = vector(6);
end

function vector = inertia_tensor2vec_casadi(I)
import casadi.*

vector = SX.zeros(6, 1);
vector(1) = I(1, 1); vector(2) = I(1, 2); vector(3) = I(1, 3);
vector(4) = I(2, 2); vector(5) = I(2, 3); vector(6) = I(3, 3);
end

function matrix = translation_transfmat_casadi(v)
import casadi.*

if ~isa(v, 'casadi.SX') && ~isa(v, 'casadi.MX') && ~isa(v, 'casadi.DM')
    v = SX(v);
end

matrix = SX.eye(4);
matrix(1:3, 4) = v;
end

function r = ml2r_casadi(m, l)
import casadi.*

r = l / m;
end

function I = Lmr2I_casadi(L, m, r)
import casadi.*

r_skew = vec2so3_casadi(r);
I = L - m * r_skew' * r_skew;
end

function debug_p_c(rbt_df, p_c_expr, frame_idx, varargin)
import casadi.*

q_vec = vertcat(rbt_df.coordinates{:});
dq_vec = vertcat(rbt_df.d_coordinates{:});
ddq_vec = vertcat(rbt_df.dd_coordinates{:});
param_vec = rbt_df.bary_params(:);

if numel(varargin) == 4
    q_test     = varargin{1};
    dq_test    = varargin{2};
    ddq_test   = varargin{3};
    param_test = varargin{4};
else
    dof       = numel(rbt_df.coordinates);
    param_num = numel(rbt_df.bary_params);
    q_test     = (1:dof)';
    dq_test    = 0.1*(1:dof)';
    ddq_test   = 0.01*(1:dof)';
    param_test = ones(param_num,1);
end

p_c_func = Function('p_c_func', {q_vec, dq_vec, ddq_vec, param_vec}, {p_c_expr});

p_c_num = full(p_c_func(q_test, dq_test, ddq_test, param_test));
disp(['DEBUG p_c for frame ' num2str(frame_idx) ':']);
disp(p_c_num);
end

function debug_vars(rbt_df, varsStruct, frame_idx, varargin)
import casadi.*

q_vec = vertcat(rbt_df.coordinates{:});
dq_vec = vertcat(rbt_df.d_coordinates{:});
ddq_vec = vertcat(rbt_df.dd_coordinates{:});
param_vec = rbt_df.bary_params(:);

if numel(varargin) == 4
    q_test     = varargin{1};
    dq_test    = varargin{2};
    ddq_test   = varargin{3};
    param_test = varargin{4};
else
    dof       = numel(rbt_df.coordinates);
    param_num = numel(rbt_df.bary_params);
    q_test     = (1:dof)';
    dq_test    = 0.1*(1:dof)';
    ddq_test   = 0.01*(1:dof)';
    param_test = ones(param_num, 1);
end
disp(['DEBUG Geom Variables for Frame ', num2str(frame_idx), ' (numeric values):']);
fields = fieldnames(varsStruct);
for i = 1:length(fields)
    name = fields{i};
    expr = varsStruct.(name);
    try
        F = Function('F', {q_vec, dq_vec, ddq_vec, param_vec}, {expr});
        dm = full(F(q_test, dq_test, ddq_test, param_test));
        disp(['--- ', name, ':']);
        disp(dm);
    catch
        disp(['--- ', name, ' (symbolic):']);
        disp(expr);
    end
end
end

function display_casadi_symbols(expr)
try
    symbols = symvar(expr);
    
    fprintf('表达式: ');
    disp(expr);
    fprintf('符号变量数量: %d\n', length(symbols));
    
    if length(symbols) == 0
        fprintf('表达式中没有符号变量\n');
        return;
    end
    
    fprintf('\n符号变量列表:\n');
    fprintf('================\n');
    
    for i = 1:length(symbols)
        try
            symbol_name = symbols(i).name();
            fprintf('%d. %s\n', i, symbol_name);
        catch
            fprintf('%d. (unnamed_symbol)\n', i);
        end
    end
    
    fprintf('================\n\n');
    
catch ME
    fprintf('Error displaying symbols: %s\n', ME.message);
end
end