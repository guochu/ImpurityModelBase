# 精确对角化

本模块提供杂质模型的精确对角化（ED）数值方法。核心思路是：

1. 用**项（term）**组合出哈密顿量（`AdagATerm`、`AdagAdagTerm`、`AATerm`、`QuarticTerm`）；
2. 二次型哈密顿量映射为**系数矩阵** `cmatrix`（BCS 型为 2L×2L BdG 形式）；
3. 基于特征分解缓存 `EigenCache` 高效计算热态、时间演化、格林函数与关联函数。

对二次型哈密顿量，本包区分两条路线（见[简介](index.md)）：

- **cdm 路线**：单粒子系数密度矩阵，`L×L`，任何二次观测量 ``\langle A \rangle = \mathrm{tr}(\mathrm{cdm}\cdot A)``；
- **dm 路线**：Fock 空间中的真密度矩阵（费米子 `2^L` 维、玻色子 `d^L` 维），可计算任意观测量。

哈密顿量 ``\hat h = \sum_{ij} h_{ij}\, \hat c_i^\dagger \hat c_j`` 的系数矩阵 ``h`` 满足
自由演化方程 ``\mathrm{d}\rho/\mathrm{d}t = -i[h^\mathsf{T}, \rho]``。

## 哈密顿量项

```@docs
AbstractTerm
QuadraticTerm
AdagATerm
tunneling
AdagAdagTerm
adagadag
AATerm
aa
QuarticTerm
interaction
NormalTerm
```

### 哈密顿量容器

```@docs
NormalHamiltonian
NormalQuadraticHamiltonian
GenericQuadraticHamiltonian
quadratichamiltonian
hamiltonian
```

示例——三格点紧束缚 + 相互作用：

```julia
using ImpurityModelBase

# h = ε (n₁+n₂) + t (c₁†c₂ + h.c.) + U n₁↑ n₁↓ 型的项级构造
ham = NormalQuadraticHamiltonian(3)
push!(ham, adaga(1, 1, coeff=0.5))    # ε n₁
push!(ham, adaga(2, 2, coeff=0.5))    # ε n₂
push!(ham, tunneling(1, 2, coeff=1.0))
push!(ham, tunneling(2, 1, coeff=1.0))
push!(ham, tunneling(2, 3, coeff=1.0))
push!(ham, tunneling(3, 2, coeff=1.0))
```

`hamiltonian` 也可直接从离散浴或模型构造（见 [Toulouse 模型](@ref) 与
[边界驱动模型](@ref)）。

### 系数矩阵

```@docs
cmatrix
```

- `NormalQuadraticHamiltonian` → `L×L` 矩阵；
- `GenericQuadraticHamiltonian`（含 `â†â†`/`ââ` 项，即 BCS 型）→ 反对称化的 `2L×2L`
  Bogoliubov–de Gennes 矩阵（由 [`ImpurityModelBase.bcs_cmatrix`](@ref) 从粒子-空穴块
  `h` 与配对块 `g` 组装）。

```@docs
ImpurityModelBase.bcs_cmatrix
```

## 特征分解与时间演化

```@docs
ImpurityModelBase.EigenCache
eigencache
timeevo
```

`timeevo(ρ₀, h, t, cache)` 返回 ``e^{t h}\rho_0 e^{t h^\dagger}``；
实时间演化取 `t = -im * time`，虚时间演化取 `t = -τ`（直接用于 `G(τ)`）。

## 热态

### 系数密度矩阵（cdm）

```@docs
thermocdm
fermionicthermocdm
```

cdm 是热高斯态的单粒子关联矩阵 ``\rho_{ij} = \langle \hat c_j^\dagger \hat c_i \rangle``，
由本征模式的占据数 `thermaloccupation(P, β, μ, λ)` 组装：

```julia
julia> cache = eigencache(h);

julia> ρ = fermionicthermocdm(cache; β=2.0, μ=0.1);  # L×L

julia> n₁ = real(ρ[1, 1])  # ⟨n̂₁⟩
```

### 真密度矩阵（dm）

```@docs
thermodm
fermionicthermodm
ImpurityModelBase.bosonicthermodm
```

!!! warning
    dm 的维度随模式数指数增长（费米 `2^L`、玻色 `d^L`），仅适用于小系统或作为
    cdm 路线的交叉验证基准。

## 格林函数

### 自由费米子（ED）

```@docs
freefermions_Gt
freefermions_greater_lesser
freefermions_Gτ
```

平衡态（给 `β`、`μ`）与非平衡（给初始 cdm `ρ₀`）均支持；给时间向量返回数组，
否则返回可调用函数：

```julia
Gt  = freefermions_Gt(h, 1, 1, cache; β=Inf, μ=0.5)  # 函数
vals = freefermions_Gt(h, 0:0.01:1, 1, 1, cache; β=Inf, μ=0.5)  # 数组
```

### 自由玻色子（ED）

```@docs
freebosons_Gt
freebosons_greater_lesser
freebosons_Gτ
```

## 关联函数与开放系统

```@docs
correlation_2op_1t
correlation_2op_1τ
LindbladOperator
lindbladoperator
steady_state
```

`correlation_2op_1t` 计算双算符单时刻关联 ``\langle \hat a(t) \hat b\rangle``
（`reverse=true` 时为 ``\langle \hat a \hat b(t)\rangle``）。

`lindbladoperator(H, jumpops)` 构造 Markov 主方程生成元

```math
\frac{\mathrm{d}\rho}{\mathrm{d}t} = -i[H, \rho] + \sum_k \left( 2 J_k \rho J_k^\dagger - \{J_k^\dagger J_k, \rho\} \right),
```

配合 `timeevo(ρ, L, t)` 与 `steady_state(L)` 使用。

## Fock 空间算符

以上均为系数矩阵层面的操作。若需要 Fock 空间中的显式矩阵表示（用于 dm 路线或
非二次观测量），可用：

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

## Toulouse 模型

Toulouse 模型即单杂质能级耦合离散浴：

```math
H = \epsilon_d\, d^\dagger d + \sum_k \epsilon_k\, c_k^\dagger c_k
+ \sum_k V_k \left( d^\dagger c_k + c_k^\dagger d \right).
```

```@docs
Toulouse
```

```julia
# 半圆谱零温费米浴 → 离散化 → Toulouse 模型
fb = fermionicvacuum(semicircular(1.0), μ=0.5)
m = Toulouse(discretevacuum(fb, δw=0.1), ϵ_d=-0.3)
```

模型级 API：

```@docs
separablecdm
separabledm
toulouse_greater_lesser
toulouse_neq_greater_lesser
particlecurrent_cmatrix
particlecurrent_hamiltonian
ImpurityModelBase.freehamiltonian
ImpurityModelBase.toulouse_cmatrix
```

其中 `thermocdm(m)` / `thermodm(m)` 给出 Toulouse 模型的热平衡 cdm / dm
（后者对玻色子浴需指定截断维数 `d`）；`toulouse_Gτ(m)` / `toulouse_Gt(m)`
为模型的 Matsubara 与实时间格林函数（其连续谱解析版本见 [解析解](@ref)）；
`heatcurrent_cmatrix` / `heatcurrent_hamiltonian` 与 `particlecurrent_*`
同理，给出热流算符的两种表示。

`separablecdm`/`separabledm` 构造杂质与浴可分离的初态（杂质占据 `nsys`，
浴处于热平衡），用于非平衡淬火问题。

!!! note
    BCS 浴的 Toulouse 模式顺序为 ``\hat a_1^\dagger \hat a_{-1}^\dagger
    \hat a_2^\dagger \hat a_{-2}^\dagger \cdots \hat a_1 \hat a_{-1} \cdots``，
    且要求 `μ = 0`。

## 边界驱动模型

边界驱动（boundary driving）输运模型：中间系统（`hsys` 为 L×L 矩阵）两端耦合
左、右两个离散浴，用于非平衡稳态输运模拟。不支持 BCS/BEC 浴。

```@docs
BoundaryDriving
```

```julia
# 泵-探测型输运：左浴高温高 μ，右浴零温
left = discretebath(Fermion, -5:0.1:5, semicircular(5.0); β=1.0, μ=1.0)
right = discretebath(Fermion, -5:0.1:5, semicircular(5.0); β=Inf, μ=0.0)
hsys = [0.0 1.0; 1.0 0.0]
m = BoundaryDriving(hsys, left, right)
```

`thermocdm(m)`（要求两浴同 `β`、`μ`）、`fermionicthermodm(m)`、
`separablecdm(m, ρ_sys)` 与左右粒子/热流算符与 Toulouse 模型的接口一致：

```@docs
ImpurityModelBase.leftparticlecurrent_cmatrix
ImpurityModelBase.leftheatcurrent_cmatrix
ImpurityModelBase.leftparticlecurrent_hamiltonian
ImpurityModelBase.leftheatcurrent_hamiltonian
ImpurityModelBase.fermionicseparabledm
```
