# 粒子浴

粒子浴（bath）是杂质模型的主体。一个浴由三部分定义：

1. **粒子种类** `P`（`Boson` 或 `Fermion`），
2. **谱密度**（见 [谱密度函数](@ref)），
3. **热力学状态**：逆温度 `β` 与化学势 `μ`（零温用 `Vacuum`，`β = Inf`）。

## 粒子种类与热分布

```@docs
AbstractParticle
Boson
Fermion
boseeinstein
fermidirac
thermaloccupation
```

两种分布函数分别为

```math
n_B(\epsilon) = \frac{1}{e^{\beta\epsilon} - 1}, \qquad
n_F(\epsilon) = \frac{1}{e^{\beta\epsilon} + 1}.
```

`thermaloccupation` 按粒子种类自动分发：

```julia
julia> thermaloccupation(Boson, 1.0, 0.0, 1.0) ≈ boseeinstein(1.0, 1.0)
true

julia> thermaloccupation(Fermion, 1.0, 0.5, 0.2) ≈ fermidirac(1.0, 0.5, 0.2)
true
```

!!! note
    `boseeinstein` 要求 `ϵ ≥ 0`（即能量在化学势之上）；`ϵ - μ < 0` 时玻色分布发散，
    会抛出 `ArgumentError`。

## 类型层级

```text
AbstractBath{P}
├── AbstractNormalBath{P}          # 正常（未配对）浴：谱 + β + μ
│   ├── Bath{P, F}                 #   有限温度
│   ├── Vacuum{P, F}               #   零温（β = Inf）
│   ├── DiscreteBath{P}            #   = Bath{P, DiscreteSpectrum}
│   └── DiscreteVacuum{P}          #   = Vacuum{P, DiscreteSpectrum}
├── AbstractBCSBath                # BCS 超导费米浴：谱 + β + μ + Δ
│   ├── BCSBath{F, T}
│   ├── BCSVacuum{F, T}
│   ├── DiscreteBCSBath{T}
│   └── DiscreteBCSVacuum{T}
└── AbstractBECBath                # BEC 玻色浴（尚未实现）
```

```@docs
AbstractBath
ImpurityModelBase.AbstractNormalBath
AbstractBCSBath
AbstractBECBath
AbstractFermionicNormalBath
AbstractBosonicNormalBath
particletype
```

## 正常浴（有限温度）

```@docs
Bath
BosonicBath
FermionicBath
bosonicbath
fermionicbath
bath
```

示例——构造有限温费米浴并查询其热力学量：

```julia
julia> f = semicircular(1.0);

julia> fb = fermionicbath(f; β=0.25, μ=0.5);

julia> fb.β, fb.T, fb.μ
(0.25, 4.0, 0.5)

julia> thermaloccupation(fb, 0.3) ≈ 1 / (exp(0.25 * (0.3 - 0.5)) + 1)
true
```

`Bath` 的谱密度可通过属性 `spectrum`（即字段 `f`）访问。

## 真空浴（零温）

```@docs
Vacuum
BosonicVacuum
FermionicVacuum
bosonicvacuum
fermionicvacuum
vacuum
```

`Vacuum` 中 `β` 恒为 `Inf`、`T` 恒为 `0`，只需指定 `μ`。

## BCS 超导浴

BCS 浴在正常浴之上增加配对参数 `Δ`（可为复数），描述超导库：

```math
H_{\mathrm{bath}} = \sum_k \epsilon_k c_k^\dagger c_k
- \sum_k \left( \Delta\, c_k^\dagger c_{-k}^\dagger + \Delta^*\, c_{-k} c_k \right).
```

```@docs
BCSBath
BCSVacuum
bcsbath
bcsvacuum
DiscreteBCSBath
DiscreteBCSVacuum
discretebcsbath
discretebcsvacuum
```

`bcsbath` 也可以从已有的正常费米浴加 `Δ` 升级得到：

```julia
julia> fb = fermionicbath(semicircular(1.0); β=3.0, μ=0.1);

julia> bcsbath(fb; Δ=0.4);  # 保留原浴的谱、β、μ，加上配对 Δ
```

!!! warning
    离散 BCS 浴的模式数是频率数的两倍（``k`` 与 ``-k`` 两个模式配对），
    故 `num_sites` 返回 `2 * length(ws)`。

## 离散浴与离散化

精确对角化只能处理有限个浴模式，因此连续谱必须先离散化。
本包实现了 *How to discretize a quantum bath for real-time evolution*
一文中 Eqs. (10a, 10b) 的方案：在频率区间 ``[w_n, w_{n+1}]`` 上，

```math
|V_n|^2 = \int_{w_n}^{w_{n+1}} dx\, f(x), \qquad
X_n = \frac{\int_{w_n}^{w_{n+1}} dx\, x\, f(x)}{|V_n|^2},
```

即每个离散模式携带区间内的总谱权重，且位于该区间的谱重心处。
这样保证了关联函数的长时行为（实时间演化）被正确捕捉。

```@docs
ImpurityModelBase.spectrum_couplings
```

```@docs
DiscreteBath
DiscreteVacuum
discretebosonicbath
discretefermionicbath
discretebosonicvacuum
discretefermionicvacuum
discretebath
discretevacuum
num_sites
```

三种离散化方式：

```julia
# (a) 直接给出频率与谱值
db = discretebath(Fermion, [0.5, 1.0, 1.5], [0.3, 0.5, 0.2]; β=5.0, μ=0.2)

# (b) 给定频率端点，对连续谱 f 做区间离散化（谱重心方案）
db = discretebath(Fermion, 0:0.1:10, f; β=5.0, μ=0.2)

# (c) 直接离散化一个连续浴，步长 δw
db = discretebath(fb; δw=0.1)
```

方式 (b)/(c) 会丢弃谱权重小于 `atol`（默认 `1e-6`）的区间，
因此返回的模式数可能少于 `length(freqs) - 1`。

## BEC 浴

`AbstractBECBath`（玻色-爱因斯坦凝聚浴）的类型已定义，但具体实现尚未完成。
