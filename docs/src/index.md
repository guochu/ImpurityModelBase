# ImpurityModelBase

**ImpurityModelBase** 是一个用于量子杂质模型（quantum impurity model）研究的 Julia 基础库。
它围绕「一个（或少数几个）杂质自由度耦合到巨大甚至无限大的粒子浴」这一物理图像，提供从
模型定义、精确对角化（ED）数值求解、到解析解基准的一整套工具。

## 功能概览

| 模块 | 功能 |
|:---|:---|
| 谱密度函数 | 连续谱（含预置半圆谱、Leggett 谱）、离散谱、δ 峰谱的统一表示与积分 |
| 粒子浴 | 玻色/费米浴、真空（零温）浴、BCS 超导浴；连续谱的自动离散化 |
| 精确对角化 | 二次哈密顿量 → 系数矩阵 → 热态、时间演化、格林函数、流；Lindblad 开放系统动力学 |
| 解析解 | 自由费米/玻色子、Toulouse 模型、独立玻色子模型、自旋-玻色退相干（含动力学去耦） |
| 工具函数 | 实/虚时间格林函数的傅里叶变换、线性预测外推 |

## 安装

本包尚未注册到 General registry，请在包管理器中直接通过 URL 安装：

```julia
using Pkg
Pkg.add(url="https://github.com/guochu/ImpurityModelBase.git")
```

或对于本地开发：

```julia
Pkg.develop(path="/path/to/ImpurityModelBase")
```

## 快速上手

下面用一个完整示例展示典型工作流：定义半圆谱费米浴 → 离散化 → 构建 Toulouse 模型
（单个杂质能级耦合浴）→ 计算热态与实时间格林函数。

```julia
using ImpurityModelBase

# 1. 半圆谱密度（带宽 2t），零温费米浴，化学势 μ = 0.5
f = semicircular(1.0)
fb = fermionicvacuum(f, μ=0.5)

# 2. 以步长 δw = 0.1 离散化为 L 个浴模式
db = discretevacuum(fb, δw=0.1)

# 3. Toulouse 模型：杂质能级 ϵ_d = -0.3 耦合到浴
m = Toulouse(db, ϵ_d=-0.3)

# 4. 系数密度矩阵（单粒子关联矩阵）与杂质占据
cdm = thermocdm(m)
n_imp = real(cdm[1, 1])

# 5. 杂质实时间格林函数 G(t)（t = 0.5 处取值）
h = cmatrix(m)
cache = eigencache(h)
Gt = freefermions_Gt(h, 1, 1, cache; β=Inf, μ=0.5)
Gt(0.5)
```

物理上，Toulouse 模型的哈密顿量为

```math
H = \epsilon_d\, d^\dagger d + \sum_k \epsilon_k\, c_k^\dagger c_k + \sum_k V_k \left( d^\dagger c_k + c_k^\dagger d \right),
```

浴谱密度 ``J(\omega) = 2\pi \sum_k |V_k|^2 \delta(\omega - \epsilon_k)`` 由 `semicircular` 给出。
更多细节见 [精确对角化](@ref) 与 [解析解](@ref)。

## 设计理念：cdm 与 dm

本包的核心设计之一是区分两种"密度矩阵"：

- **dm（密度矩阵）**：Fock 空间中的多体密度矩阵 ``e^{-\beta H}/Z``，维度随模式数指数增长
  （费米子 ``2^L``，玻色子 ``d^L``）。适用于包含相互作用的一般哈密顿量。
- **cdm（系数密度矩阵）**：单粒子关联矩阵 ``\rho_{ij} = \langle c_j^\dagger c_i \rangle``，仅为
  ``L \times L``。对二次型（无相互作用）哈密顿量，任何二次观测量的期望都可通过
  ``\langle A \rangle = \mathrm{tr}(\rho A)`` 以 ``\mathcal{O}(L^2)`` 的代价获得，
  无需构造指数维的多体态。

相应地，流算符（粒子流、热流）等也都有 `*_cmatrix`（基于 cdm）与 `*_hamiltonian`
（基于 dm）两套实现，前者高效，后者通用。

## 文档结构

- [谱密度函数](@ref)：谱密度的类型体系、预置谱、δ 峰与积分
- [粒子浴](@ref)：浴的构造、层级与离散化
- [精确对角化](@ref)：哈密顿量、热态、时间演化、格林函数、关联函数、开放系统
- [解析解](@ref)：各模型的精确解析结果，可用作 ED 的基准
- [工具函数](@ref)：傅里叶变换与线性预测

## 引用

若本库对你的研究有帮助，请引用：

```bibtex
@software{ImpurityModelBase,
  author = {Guo, Chu},
  title  = {ImpurityModelBase: A Julia library for quantum impurity models},
  url    = {https://github.com/guochu/ImpurityModelBase}
}
```
