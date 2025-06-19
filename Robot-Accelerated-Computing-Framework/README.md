# Robot-Accelerated-Computing-Framework(RACF)

基于CasADi的机器人建模与分析加速框架

## 项目简介

Robot-Accelerated-Computing-Framework是一个用于机器人动力学建模、分析和参数识别的MATLAB工具包。该项目利用CasADi符号计算框架，充分利用CasAdi的自动微分、稀疏矩阵、梯度优化等特性，实现了机器人的运动学计算、动力学建模、回归矩阵生成以及基惯性参数计算等功能。

主要特性：

- 支持串、并联机器人的符号动力学建模
- 基于旋量理论的运动学计算
- 基于拉格朗日方程的动力学计算
- 支持带约束的牛顿欧拉递归动力学推导（待上传）
- 自动生成动力学回归矩阵
- 基惯性参数计算和优化
- 支持非线性建模
- 支持弹簧关节建模
- 高效的符号计算和代码生成

## 支持的机器人模型

该项目基于腔镜手术机器人带有并联关节的7自由度主操作臂、远心不动点结构的从操作臂及工业6轴机器人实验平台验证

- **MTM**: 7自由度主操作臂
- **RCM**: 4自由度远程中心运动机器人
- 常见业6轴串联机器人、冗余7自由度机器人

## 系统要求

### 必需软件

- MATLAB R2018b或更高版本
- CasADi v3.5+（MATLAB接口）

### 推荐硬件

- 内存：8GB或更大
- 处理器：Intel i5或更高

## 安装说明

1. **克隆项目**

   ```bash
   git clone https://github.com/songyp2020/Robot-Accelerated-Computing-Framework
   cd Robot-Accelerated-Computing-Framework
   ```
2. **安装CasADi**

   - 从[CasADi官网](https://web.casadi.org/)下载MATLAB版本
   - 将CasADi解压到合适目录
3. **配置路径**

   - 修改 `main_casadi.m`中的CasADi路径：

   ```matlab
   addpath('casADI Install Path')  % 修改为你的CasADi安装路径
   ```

## 使用方法

### 快速开始

1. **运行主程序**

   ```matlab
   main_casadi
   ```
2. **程序执行流程**

   - 步骤1：定义机器人模型
   - 步骤2：计算运动学
   - 步骤3：计算动力学
   - 步骤4：计算基惯性参数
   - 步骤5：生成回归矩阵函数
   - 步骤6：回归矩阵及最小回归矩阵高性能代码生成

### 自定义机器人模型

要添加新的机器人模型，请参考 `model/Demo7DOF.m`的格式：

```matlab
function rbt_config = YourRobot_casadi()
import casadi.*

% 定义关节变量
q1 = SX.sym('q1');
% ... 更多关节

% 定义几何参数
% ...

% 配置机器人结构
rbt_config = {
    % [连杆号, 前驱, 后继, 轴向量, 位置, 姿态, 角度, 节距, 惯性, 弹簧]
    % ...
};
end
```

## 文件结构

```
MTM-casADI/
├── main_casadi.m                    # 主程序入口
├── model/                           # 机器人模型定义
│   ├── MTM_casadi.m                # MTM机器人模型
│   └── RCM14_casadi.m              # RCM14机器人模型
├── src_casadi/                      # 核心算法实现
│   ├── DefineRobot_casadi.m        # 机器人定义和参数处理
│   ├── GeometryCalculation_casadi.m # 运动学计算
│   ├── Dynamics_casadi.m           # 动力学计算
│   ├── DynamicBaseParamCalc_casadi.m # 基参数计算
│   ├── make_function_casadi.m      # 函数生成
│   └── utils_casadi.m              # 工具函数库
└── README.md                        # 项目说明文档
```

### 核心模块说明

#### DefineRobot_casadi.m

- 机器人模型的符号定义
- 旋量轴和变换矩阵计算
- 惯性参数符号化
- 坐标变量生成

#### GeometryCalculation_casadi.m

- 前向运动学计算
- 雅可比矩阵生成
- 质心位置和速度计算

#### Dynamics_casadi.m

- 基于拉格朗日方程的动力学建模
- 质量矩阵、重力项、科里奥利力计算
- 弹簧力矩处理
- 动力学回归矩阵生成

#### DynamicBaseParamCalc_casadi.m

- 基惯性参数识别
- 参数线性相关性分析
- 最小参数集计算

#### make_function_casadi.m

- 高性能C/C++代码生成

## 理论基础

### 旋量理论

项目基于旋量理论进行运动学建模，使用指数映射表示关节运动：

$$
T = e^{[S]\theta} \cdot M
$$

其中：

- `[S]`是旋量的反对称矩阵表示
- `θ`是关节角度
- `M`是初始位姿矩阵

### 动力学建模

采用拉格朗日方程进行动力学建模：

$$
\tau = \frac{d}{dt}\left(\frac{\partial L}{\partial \dot{q}}\right) - \frac{\partial L}{\partial q}
$$

其中拉格朗日函数 `L = T - V`（动能减势能）

### 参数线性化

动力学方程可表示为参数线性形式：

$$
\tau = H(q,\dot{q},\ddot{q}) \cdot \pi
$$

其中：

- `H`是回归矩阵
- `π`是惯性参数向量

## API参考

### 主要函数

#### DefineRobot_casadi(model_name, params)

**参数**:

- `model_name`: 机器人名称
- `params`: 机器人配置参数

**返回**: 机器人结构体

#### GeometryCalculation_casadi(rbt_df, static_model)

**参数**:

- `rbt_df`: 机器人定义结构体
- `static_model`: 是否为静态模型

**返回**: 几何计算结果

#### Dynamics_casadi(rbt_df, geom, gravity_vec)

**参数**:

- `rbt_df`: 机器人定义
- `geom`: 几何计算结果
- `gravity_vec`: 重力向量

**返回**: 动力学计算结果

## 输出结果

程序运行后会生成：

1. **回归矩阵函数** (`H_func`, `H_b_func`)
2. **基惯性参数集**
3. **动力学模型**（质量矩阵、重力项等）
4. **运动学模型**（变换矩阵、雅可比矩阵等）

## 示例应用

- 机器人动力学参数识别
- 控制器设计与仿真
- 机器人性能分析
- 运动规划算法验证

## 常见问题

**Q: CasADi路径配置错误怎么办？**
A: 检查 `main_casadi.m`中的 `addpath`路径是否正确指向CasADi安装目录。

**Q: 内存不足怎么办？**
A: 对于复杂机器人模型，建议增加系统内存或简化模型复杂度。

**Q: 如何添加新的机器人模型？**
A: 参考 `model/Demo7DOF.m`的格式，创建新的基于旋量的机器人模型文件并在主程序中调用。

## 贡献指南

欢迎提交问题报告和改进建议。请确保：

1. 代码风格一致
2. 添加适当的文档说明
3. 测试新功能的正确性

## 许可证

[]

## 联系方式

ssongyp@163.com

---

*最后更新: 2025年*
