# 精确对角化

本模块提供全 Fock 空间的暴力精确对角化（ED）：显式构造费米子 `2^L` 维、玻色子 `d^L`
维空间中的算符矩阵与多体密度矩阵，可计算**任意**观测量（不限于二次型）。维度随模式数
指数增长，因此仅适用于小系统，或作为 [二次型系统](@ref) 高效 cdm 路线的交叉验证基准
（测试套件即按此方式对标）。

!!! warning
    二次型哈密顿量的二次观测量请优先使用 [二次型系统](@ref) 的系数矩阵方法；
    本页方法只为小系统、非二次观测量和 benchmark 而设计。

## Fock 空间算符

```@docs
fermionadagoperator
fermionaoperator
fermiondensityoperator
fermionoccupationoperator
fermionoperator
bosonadagoperator
bosonaoperator
bosondensityoperator
bosonoccupationoperator
bosonoperator
```

费米子算符自动包含 Jordan–Wigner 弦；玻色子算符需要指定截断维数 `d`。

## 真密度矩阵（dm）

```@docs
thermodm
fermionicthermodm
ImpurityModelBase.bosonicthermodm
```

!!! warning
    dm 的维度随模式数指数增长（费米 `2^L`、玻色 `d^L`），仅适用于小系统或作为
    cdm 路线（见 [二次型系统](@ref)）的交叉验证基准。

## 关联函数与开放系统

```@docs
correlation_2op_1t
correlation_2op_1τ
LindbladOperator
lindbladoperator
steady_state
```

`correlation_2op_1t` 计算双算符单时刻关联 ``\langle \hat a(t) \hat b\rangle``
（`reverse=true` 时为 ``\langle \hat a \hat b(t)\rangle``），在算符矩阵所在的
全空间（通常为 Fock 空间）中演化。

`lindbladoperator(H, jumpops)` 构造 Markov 主方程生成元

```math
\frac{\mathrm{d}\rho}{\mathrm{d}t} = -i[H, \rho] + \sum_k \left( 2 J_k \rho J_k^\dagger - \{J_k^\dagger J_k, \rho\} \right),
```

配合 `timeevo(ρ, L, t)` 与 `steady_state(L)` 使用。
