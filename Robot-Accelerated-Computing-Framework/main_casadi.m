clear
clc

import casadi.*

addpath('C:\Users\syp\Desktop\MTM-casADI\model');
addpath('C:\Users\syp\Desktop\MTM-casADI\src_casadi');
addpath('C:\Install\casADI')

rbt_config = MTM_casadi;
gravity_vec = [0, 0, -9.8015];

rbt = struct();
rbt.name = 'MTM_casadi';

disp('==========================================');
disp('步骤1: 使用CasADi定义机器人模型...');
tic;
rbt.rbt_df = DefineRobot_casadi('MTM_casadi', rbt_config);
toc;

disp('==========================================');
disp('步骤2: 计算运动学...');
tic;
rbt.geom = GeometryCalculation_casadi(rbt.rbt_df, false);
toc;

disp('==========================================');
disp('步骤3: 计算动力学...');
tic;
rbt.dyn = Dynamics_casadi(rbt.rbt_df, rbt.geom, gravity_vec);
toc;

disp('==========================================');
disp('步骤4: 计算基惯性参数...');
tic;
base = DynamicBaseParamCalc_casadi(rbt.rbt_df, rbt.dyn);
disp(base);
toc;

disp('==========================================');
disp('步骤4: 生成回归矩阵函数...');
tic;
[H_func, H_b_func] = make_function_casadi(rbt);
toc;

rmpath('C:\Users\syp\Desktop\MTM-casADI\model');
rmpath('C:\Users\syp\Desktop\MTM-casADI\src_casadi');



