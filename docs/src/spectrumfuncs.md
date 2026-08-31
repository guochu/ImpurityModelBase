# 谱密度函数

杂质模型中，浴的物理信息完全由其**谱密度**（spectral density）刻画：

```math
J(\omega) = 2\pi \sum_k |V_k|^2 \delta(\omega - \epsilon_k),
```

其中 ``V_k`` 是杂质与第 ``k`` 个浴模式的耦合强度，``\epsilon_k`` 是模式能量。
本模块提供连续谱、离散谱与 δ 峰谱的统一抽象 `AbstractBoundedFunction`，
使浴的定义与后续 ED 离散化、解析解积分共享同一套接口。

## 有界谱函数

所有谱函数都实现 `AbstractBoundedFunction` 接口：函数只在区间 ``[\mathrm{lb}, \mathrm{ub}]``
内有定义，区间外取零。

```@docs
AbstractBoundedFunction
BoundedFunction
bounded
lowerbound
upperbound
```

`spectrum` 与 `bounded` 的区别在于前者会校验函数值非负——谱密度作为耦合强度的模方，
物理上必须非负：

```julia
julia> f = spectrum(ϵ -> sqrt(1 - ϵ^2), -1, 1)
```

!!! warning
    `BoundedFunction` 不做非负性检查；若将其用于浴的构造，请自行保证谱值非负，
    否则后续 `spectrumcouplings`（取平方根）等操作会失败。

## 预置谱

### 半圆谱（费米子）

```math
J(\omega) = \frac{2}{\pi t^2}\sqrt{t^2 - \omega^2}, \quad |\omega| \le t,
```

对应于无限维 Bethe 格子上的紧束缚模型（带宽 ``2t``），是最常用的费米浴谱。

```@docs
semicircular
```

### Leggett 谱（玻色子）

```math
J(\omega) = \frac{\alpha}{2} \frac{\omega^d}{\omega_c^{d-1}} e^{-\omega/\omega_c},
```

唯象玻色谱（Ohmic, ``d=1``; sub-Ohmic, ``d<1``; super-Ohmic, ``d>1``），
截断频率 ``\omega_c`` 处指数收敛。实现中以 `100ωc` 为上界。

```@docs
Leggett
```

## 离散谱

```@docs
DiscreteSpectrum
frequencies
spectrumvalues
spectrumcouplings
```

离散谱由频率向量 `ws`（升序）与谱值向量 `fs` 定义，`spectrumcouplings` 返回
``\sqrt{f_n}``，即每个离散模式的有效耦合强度。

## δ 峰谱

```@docs
DiracDelta
ImpurityModelBase.DeltaMultF
ImpurityModelBase.DeltasMultF
```

`DiracDelta(ω; α)` 表示位于 ``\omega``、强度为 ``\alpha`` 的 δ 峰。
它本身不可调用，但可以与其他有界函数组合（内部类型 `DeltaMultF` / `DeltasMultF`
表示 ``\alpha\,\delta(\omega-\omega_0) f(\omega)`` 型乘积），并被积分器正确处理。

## 积分

谱函数的标准积分通过 `quadgk` 进行，δ 峰部分被解析处理：

```@docs
quadgkwrapper
```

例如：

```julia
julia> f = spectrum(ϵ -> 1.0, 0, 1);  # [0,1] 上的常数谱

julia> quadgkwrapper(f)
(1.0, 0.0)

julia> d = DiracDelta(0.5, α=2.0);   # δ 峰谱的积分为其强度

julia> quadgkwrapper(d)
(2.0, 0.0)
```

## 频率轴平移

```@docs
spectrumshift
```

化学势 ``\mu`` 的作用等价于把谱的频率轴整体平移 ``-\mu``：在固定 ``\epsilon``
的谱上加上化学势时（`hamiltonian(b; include_chemical=true)` 的另一种实现方式），
`DiscreteBath` 会在能量 ``\epsilon_k - \mu`` 处放置模式。
