# ImpurityModelBase

**ImpurityModelBase** 是一个用于量子杂质模型（quantum impurity model）研究的 Julia 基础库。
它围绕「一个（或少数几个）杂质自由度耦合到巨大甚至无限大的粒子浴」这一物理图像，提供从
模型定义、精确对角化（ED）数值求解、到解析解基准的一整套工具。

## 功能概览

| 模块 | 功能 |
|:---|:---|
| [谱密度函数](#谱密度函数) | 连续谱（含预置半圆谱、Leggett 谱）、离散谱、δ 峰谱的统一表示与积分 |
| [粒子浴](#粒子浴) | 玻色/费米浴、真空（零温）浴、BCS 超导浴；连续谱的自动离散化 |
| [精确对角化](#精确对角化) | 二次哈密顿量 → 系数矩阵 → 热态、时间演化、格林函数、流；Lindblad 开放系统动力学 |
| [解析解](#解析解) | 自由费米/玻色子、Toulouse 模型（含半圆浴的实/虚时间闭式解）、独立玻色子模型、自旋-玻色退相干（含动力学去耦）、Holstein 模型 CFE 解与 Bethe 晶格 DMFT |
| [工具函数](#工具函数) | 实/虚时间格林函数的傅里叶变换、线性预测外推 |

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

# 1. 半圆谱密度（半带宽 t = 1），零温费米浴，化学势 μ = 0.5
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

### 解析解示例

对常见的半圆谱浴，Toulouse 模型的实频率、实时间、Matsubara 频率与虚时间格林函数
**全部有闭式解**（无需数值积分与展宽参数），可直接用作 ED 数值结果的基准：

```julia
using ImpurityModelBase

# 半圆浴（半带宽 t = 1）上的推迟 / 实时间 / Matsubara / 虚时间格林函数
Gw = fermionic_toulouse_Gw_semicircular(0.5; ϵ_d=0.3, μ=0.0, t=1.0)
Gt = fermionic_toulouse_Gt_semicircular(0.7; ϵ_d=0.3, μ=0.0, t=1.0)
Giw = fermionic_toulouse_Giw_semicircular(2.1; ϵ_d=0.3, μ=0.0, t=1.0)
Gτ = fermionic_toulouse_Gτ_semicircular(0.4; β=2.0, ϵ_d=0.3, μ=0.0, t=1.0)
```

对通用谱密度，相应的数值实现为 `toulouse_Gw` / `toulouse_Gt` / `toulouse_Giw` / `toulouse_Gτ`。

Holstein 模型（电子-声子耦合）提供连续分数展开（CFE）解析解与 Bethe 晶格上的 DMFT 自洽求解
（Ciuchi, de Pasquale, Fratini & Feinberg, PRB 56, 4494 (1997)），标度参数 λ = g²/(ω₀t)、γ = ω₀/t：

```julia
ws = collect(range(-3.0, 3.0; length=801))
r = holstein_dmft_bethe(ws; λ=0.75, γ=2.0)   # r.A 即收敛后的谱函数 A(ω)
```

复现 Ciuchi et al. 论文各图表的完整教程见
[`docs/tutorials/holstein/holstein_tutorial.ipynb`](docs/tutorials/holstein/holstein_tutorial.ipynb)。

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

## 模块导览

### 谱密度函数

谱密度函数 `AbstractBoundedFunction` 是浴的统一抽象：连续谱（`BoundedFunction`）、
预置谱（`semicircular`、`Leggett`）、离散谱（`DiscreteSpectrum`）与 δ 峰（`DiracDelta`）
共用同一套积分（`quadgkwrapper`）与变换（`spectrumshift`）接口。

### 粒子浴

浴对象把谱密度与统计（`Fermion`/`Boson`）、温度（`β`）、化学势（`μ`）打包：
热浴 `fermionicbath`/`bosonicbath`，零温真空 `fermionicvacuum`/`bosonicvacuum`，
BCS 超导浴 `bcsbath`。连续谱可通过 `discretebath`/`discretevacuum` 自动离散化为
有限模式浴，供精确对角化使用。

### 精确对角化

二次哈密顿量（`QuadraticTerm` 组合、`quadratichamiltonian`）通过对角化得到
`eigencache`，之后热态关联矩阵（`thermocdm`）、时间演化（`timeevo`）、实时/虚时
格林函数（`freefermions_Gt`、`freefermions_Gτ` 等）都只做代数运算。含相互作用时
退回多体 `thermodm` + `gf_real`。BCS（`bcs_cmatrix`）与 Bogoliubov 变换、
边界驱动开放系统（Lindblad）亦有支持。

### 解析解

`src/analyticsolutions/` 收录存在精确解的模型：自由费米子/玻色子、相互作用单轨道
费米子、Toulouse 模型（任意谱密度的数值实现 + 半圆浴的全闭式解，含束缚态分析与
严格的因果性/求和规则）、独立玻色子模型（支持任意初始杂质密度矩阵）、自旋-玻色
退相干与动力学去耦、Holstein 模型的 CFE 解与 Bethe 晶格 DMFT。它们同时是
测试套件中 ED 数值结果的交叉验证基准。

### 工具函数

实/虚时间与实/虚频率格林函数之间的傅里叶变换（`Gt_to_Gw`、`Gτ_to_Giw`、
`Aw_to_Gτ` 等）以及基于线性预测（`LinearPrediction`）的长时间外推。

## 文档

完整的 API 文档（含各模型解析解的显式公式）由 Documenter 生成：

```julia
julia --project=docs docs/make.jl
```

构建结果位于 `docs/build/`，入口为 `docs/src/index.md`。

## 引用

若本库对你的研究有帮助，请引用：

```bibtex
@software{ImpurityModelBase,
  author = {Guo, Chu},
  title  = {ImpurityModelBase: A Julia library for quantum impurity models},
  url    = {https://github.com/guochu/ImpurityModelBase}
}
```
