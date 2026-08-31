# 工具函数

## 傅里叶变换

杂质模型中的格林函数通常在四个轴之间转换：

```math
G(t) \xleftrightarrow{\ \omega\ } G(\omega),
\qquad
G(\tau) \xleftrightarrow{\ i\omega_n\ } G(i\omega_n).
```

实时间与实频率之间由数值傅里叶变换连接；虚时间与 Matsubara 虚频率之间
（``\omega_n = (2n-1)\pi/\beta``，费米子）由级数求和连接。

### 谱函数

谱函数与混合函数由格林函数的虚部给出：

```math
A(\omega) = -\frac{1}{\pi} \mathrm{Im}\, G(\omega), \qquad
J(\omega) = -\frac{1}{\pi} \mathrm{Im}\, \Delta(\omega).
```

```@docs
Gw_to_Aw
Δw_to_Jw
```

### 实时间 ↔ 实频率

```@docs
Gt_to_Gw
Gw_to_Gt
Aw_to_Gτ
```

注意 `Gt_to_Gw` 要求 `G(t)` 序列已收敛到 0（配合下文的线性预测使用）；
`Gw_to_Gt` 在区间 `[wmin, wmin+δw*length]` 外使用 ``1/\omega`` 渐近展开。

### 虚时间 ↔ 虚频率

```@docs
Gτ_to_Giw
Giw_to_Gτ
Δτ_to_Δiw
ifrequency
```

`ifrequencies(β, nmax)` 生成 `nmax` 到 `nmax+1` 范围内的全部 Matsubara 频点
（共 `2nmax+2` 个，与 `Gτ_to_Giw` 的输出约定一致）。

约定：`Gτ` 输入长 `N+1` 的序列（末点冗余，因 ``G(\beta) = 1 - G(0)``），
`Giw` 输出偶数个（`2nmax+2`）频点，首频点为 ``\omega_1 = \pi/\beta``。

### 典型工作流

```julia
using ImpurityModelBase

# ED 得到虚时间 G(τ)，变换到 Matsubara 轴
gτ = ...                    # Nτ+1 个点，间隔 δτ = β/Nτ
giw = Gτ_to_Giw(gτ; β=β)

# 实时间 G(t) → G(ω)，需要先线性预测延拓
gt = ...
gt_ext = linear_predict(gt, δt)
gw = Gt_to_Gw(gt_ext, ws; δt=δt)

# 谱函数
aw = Gw_to_Aw(gw)
```

## 线性预测

实时间数值演化只能覆盖有限时长 ``T``，直接做傅里叶变换会产生频谱展宽。
线性预测（linear prediction）用自回归模型把观测序列外推到更长的时间，
显著改善变换后的频率分辨率，尤其适合关联函数呈现缓慢衰减振荡的情形。
实现依据 Barthel & White, *PRB* **79**, 245101 (2009)。

```@docs
AbstractPredictionScheme
LinearPrediction
linear_predict
```

`LinearPrediction(obs, ws; stepsize, nfit, p)` 拟合 `nfit` 个数据点的
`p` 阶自回归模型（默认 `p = div(nfit, 2)`），返回可调用对象 `x(t)`
（线性插值）。`linear_predict` 迭代外推序列直至收敛（或达到 `maxiter`），
并按 `δt` 重新采样。
