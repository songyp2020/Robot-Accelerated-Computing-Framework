function [H_func, H_b_func] = make_function_casadi(rbt)
import casadi.*
rbt_df = rbt.rbt_df;
dyn = rbt.dyn;

param_vec = rbt_df.bary_params(:);
q_vec = vertcat(rbt_df.coordinates{:});
dq_vec = vertcat(rbt_df.d_coordinates{:});
ddq_vec = vertcat(rbt_df.dd_coordinates{:});

H_symbolic = simplify(jacobian(dyn.tau, param_vec));

dyn.H_symbolic = H_symbolic;

H_func = Function('H_regressor', {q_vec, dq_vec, ddq_vec, param_vec}, {H_symbolic});
dyn.H_func = H_func;
disp('✓ Created full regressor H_func');

base = DynamicBaseParamCalc_casadi(rbt_df, dyn);
H_b_symbolic = base.H_b;
H_b_func = Function('H_b_regressor', {q_vec, dq_vec, ddq_vec, param_vec}, {H_b_symbolic});
dyn.H_b_func = H_b_func;
disp('✓ Created minimal regressor H_b_func');
end