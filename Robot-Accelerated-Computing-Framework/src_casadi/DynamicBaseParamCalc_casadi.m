function base = DynamicBaseParamCalc_casadi(rbt_df, dyn)
import casadi.*;

disp('Calculating base dynamic parameters...');

[r, P_X, P, Kd] = find_dyn_param_deps_casadi(rbt_df, dyn);

n = length(rbt_df.bary_params);
s = n - r;
KG = [eye(r), Kd; zeros(s, r), eye(s)];
eye_n = eye(n);
inertia2beta = KG * eye_n(:, P)';
prime_param = rbt_df.bary_params;
beta_param = inertia2beta * prime_param;
base_param = beta_param(1:r);
P_b = P(1:r);

if isfield(dyn, 'H_symbolic')
    H_b = simplify(dyn.H_symbolic(:, P_b));
    H_full = dyn.H_symbolic;
else
    error('dyn.H_symbolic is required for computing H_b');
end

base.base_num = r;
base.prime_num = n;
base.redundant_num = s;
base.prime_param = prime_param;
base.beta_param = beta_param;
base.base_param = base_param;
base.prime2beta = inertia2beta;
base.H_b = H_b;
base.H = H_full;
end

function [r, P_X, P, Kd] = find_dyn_param_deps_casadi(rbt_df, dyn)
import casadi.*;
dof = rbt_df.dof;
param_num = length(rbt_df.bary_params);
sample_num = param_num * 2;
Z = zeros(dof * sample_num, param_num);

if isfield(dyn, 'H_func') && isa(dyn.H_func, 'casadi.Function')
    func_H = dyn.H_func;
else
    q_vec = vertcat(rbt_df.coordinates{:});
    dq_vec = vertcat(rbt_df.d_coordinates{:});
    ddq_vec = vertcat(rbt_df.dd_coordinates{:});
    param_vec = rbt_df.bary_params(:);
    func_H = Function('H_temp', {q_vec, dq_vec, ddq_vec, param_vec}, {dyn.H_symbolic});
    disp('Warning: dyn.H_func not found, using temporary H_temp from H_symbolic');
end

epsilon = 1e-6;
param_dummy = ones(param_num, 1) * epsilon;

rng(123);
for i = 1:sample_num
    q = rand(dof,1) * 2*pi - pi;
    dq = rand(dof,1) * 2*pi - pi;
    ddq = rand(dof,1) * 2*pi - pi;
    H_val = full(func_H(q, dq, ddq, param_dummy));
    Z((i-1)*dof+1 : i*dof, :) = H_val;
end

r = rank(Z);
[~, ~, P] = qr(Z, 'vector');
[~, R] = qr(Z(:, P));
R1 = R(1:r, 1:r);
R2 = R(1:r, r+1:end);

eye_n = eye(param_num);
P_X = [eye(r), R1\R2] * eye_n(:, P)';
Kd = R1\R2;

tol = 1e-10;
P_X(abs(P_X) < tol) = 0;
Kd(abs(Kd) < tol) = 0;
end