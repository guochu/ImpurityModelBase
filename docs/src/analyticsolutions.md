# 解析解

本模块收录若干存在精确解析解的杂质/输运模型。它们有两个用途：

1. **物理研究**：直接研究这些模型的平衡/非平衡性质；
2. **数值基准**：为 [精确对角化](@ref) 的数值结果提供交叉验证（测试套件即以此为
   主要目的）。

所有谱相关的函数都接受 `AbstractBoundedFunction`（或直接接受浴对象，
自动提取其谱与热力学参数）。以下各节给出每个模型的哈密顿量与解析解的显式公式；
所有公式均与实现逐一对应。

## 自由费米子

单能级 ``\epsilon_d`` 上的自由费米子，热平衡态。本模块约定 ``\mu = -\epsilon_d``
（化学势从杂质能级处起算），平均占据为

```math
\bar n = \frac{1}{1 + e^{-\beta\mu}}.
```

实时间 greater / lesser 格林函数（``G^>(t) = -i\langle \hat a(t)\hat a^\dagger\rangle``、
``G^<(t) = i\langle \hat a^\dagger \hat a(t)\rangle``）：

```math
G^>(t) = -i\,(1-\bar n)\,e^{i\mu t}
= \frac{-i\,e^{i\mu t}}{1+e^{\beta\mu}},
\qquad
G^<(t) = i\,\bar n\,e^{i\mu t}
= \frac{i\,e^{i\mu t}}{1+e^{-\beta\mu}}.
```

```@docs
freefermion_greater
freefermion_lesser
freefermion_Gt
```

`freefermion_Gt` 实现为两分量之和：

```math
G(t) = G^>(t) + G^<(t)
= i\,\tanh\!\frac{\beta\mu}{2}\; e^{i\mu t}.
```

!!! note
    推迟（retarded）组合是 ``G^> - G^< = -i e^{i\mu t}``（自由粒子传播子）；
    若需要推迟格林函数请自行作差。

虚时间（Matsubara）格林函数，``0 \le \tau \le \beta``：

```math
G(\tau) = e^{\mu\tau}\,(1-\bar n),
```

满足反周期性质的包内约定 ``G(\beta) = 1 - G(0)``（与 `Giw_to_Gτ` 的末点约定一致）。

```@docs
freefermion_Gτ
freefermion_occupation
```

## 相互作用费米子

单轨道双自旋费米子 + Hubbard 型相互作用（精确可解，因为轨道只有一个）：

```math
H = \epsilon_d \sum_\sigma \hat n_\sigma
+ U\, \hat n_\uparrow \hat n_\downarrow
= -\mu \sum_\sigma \hat n_\sigma + U\, \hat n_\uparrow \hat n_\downarrow .
```

四个 Fock 态的玻尔兹曼权重为 ``1,\ e^{\beta\mu},\ e^{\beta\mu},\ e^{\beta(2\mu-U)}``，
配分函数 ``Z = 1 + 2e^{\beta\mu} + e^{\beta(2\mu-U)}``。实现中把分子分母同乘
``e^{-\beta\mu}``，即分母取 ``e^{-\beta\mu} + 2 + e^{\beta(\mu-U)}`` 的形式：

```math
G(\tau) = \frac{e^{\mu\tau-\beta\mu} + e^{(\mu-U)\tau}}
{e^{-\beta\mu} + 2 + e^{\beta(\mu-U)}},
```

```math
G^>(t) = -i\,\frac{e^{i\mu t-\beta\mu} + e^{i(\mu-U)t}}
{e^{-\beta\mu} + 2 + e^{\beta(\mu-U)}},
\qquad
G^<(t) = i\,\frac{e^{i\mu t} + e^{\beta(\mu-U)}\,e^{i(\mu-U)t}}
{e^{-\beta\mu} + 2 + e^{\beta(\mu-U)}}.
```

```@docs
fermion_greater
fermion_lesser
fermion_Gt
fermion_Gτ
```

`fermion_Gt` 同样实现为 ``G^>(t) + G^<(t)``。

## 自由玻色子

频率 ``\omega`` 的单模自由玻色子，平均占据

```math
\bar n_B = \frac{1}{e^{\beta\omega} - 1}.
```

```math
G^>(t) = -i\,(1+\bar n_B)\,e^{-i\omega t}
= \frac{-i\,e^{-i\omega t}}{1-e^{-\beta\omega}},
\qquad
G^<(t) = i\,\bar n_B\,e^{-i\omega t}.
```

二者之和恰为推迟格林函数 ``G^>(t)+G^<(t) = -i e^{-i\omega t}``。
虚时间轴（``0 \le \tau \le \beta``）：

```math
G(\tau) = -\,\frac{e^{-\omega\tau}}{1-e^{-\beta\omega}}
= -(1+\bar n_B)\,e^{-\omega\tau}.
```

```@docs
freeboson_greater
freeboson_lesser
freeboson_Gt
freeboson_Gτ
freeboson_occupation
```

## Toulouse 模型（解析解）

Toulouse 模型（单费米杂质 + 平衡自由费米浴）严格可解：

```math
H = \epsilon_d\, d^\dagger d + \sum_k \epsilon_k\, c_k^\dagger c_k
+ \sum_k V_k \left( d^\dagger c_k + c_k^\dagger d \right).
```

杂质格林函数由 Dyson 方程给出（谱密度 ``J(\omega)=2\pi\sum_k|V_k|^2\delta(\omega-\epsilon_k)``）：

```math
G(\omega) = \frac{1}{\omega - \epsilon_d - \Delta(\omega)},
\qquad
\Delta(\omega) = \int \mathrm{d}\epsilon\, \frac{J(\epsilon)}{\omega - \epsilon + i\delta},
```

其中 ``\Delta`` 为混合函数（hybridization function），与谱密度的关系为
``J(\omega) = -\mathrm{Im}\,\Delta(\omega)/\pi``。retarded 格林函数与温度无关。

### 实频率 / 实时间

实现中 ``\Delta`` 的频率轴以化学势为零点（``\omega + \mu - \epsilon``），而
``G`` 的分母以 ``\omega - \epsilon_d`` 计：

```math
G(\omega) = \frac{1}{\omega + i\delta - \epsilon_d - \Delta(\omega)},
\qquad
\Delta(\omega) = \int \mathrm{d}\epsilon\, \frac{J(\epsilon)}{\omega + \mu - \epsilon + i\delta}.
```

`toulouse_Δw`（无 `μ` 版本）则直接计算 ``\Delta(\omega) = \int \mathrm{d}\epsilon\,
J(\epsilon)/(\omega - \epsilon + i\delta)``。

实时间格林函数由数值傅里叶变换给出——先减去自由粒子渐近项
``1/(\omega+i\delta)`` 使被积函数收敛，其解析贡献 ``-i`` 单独补回：

```math
G(t) = \frac{1}{2\pi}\int_{\omega_{\min}}^{\omega_{\max}} \mathrm{d}\omega\,
\left[ G(\omega) - \frac{1}{\omega+i\delta} \right] e^{-i\omega t} \;-\; i .
```

```@docs
toulouse_Gw
toulouse_Δw
toulouse_Gt
```

### 虚频率 / 虚时间（Matsubara）

Matsubara 频点 ``\omega_n = (2n-1)\pi/\beta``（费米子）：

```math
G(i\omega_n) = \frac{1}{i\omega_n - \epsilon_d - \Delta(i\omega_n)},
\qquad
\Delta(i\omega_n) = \int \mathrm{d}\epsilon\, \frac{J(\epsilon)}{i\omega_n - \epsilon}.
```

虚时间格林函数由 Matsubara 级数求和得到（同样减去 ``1/(i\omega_n)`` 的自由项）：

```math
G(\tau) = \frac{1}{2} - \frac{1}{\beta}\sum_{n=-n_{\max}}^{n_{\max}+1}
\left[ G(i\omega_n) - \frac{1}{i\omega_n} \right] e^{-i\omega_n \tau}.
```

混合函数的虚时间表示（实现与测试的约定）：

```math
\Delta(\tau) = -\int \mathrm{d}\epsilon\, J(\epsilon)\,
\left(1 + e^{-\beta\epsilon}\right) e^{\epsilon\tau}.
```

```@docs
toulouse_Giw
toulouse_Gτ
toulouse_Δiw
toulouse_Δτ
```

示例——计算 Toulouse 模型在 Matsubara 轴上的格林函数并与 ED 对比：

```julia
using ImpurityModelBase

f = semicircular(1.0)
β = 10.0
giw = toulouse_Giw(f; β=β, ϵ_d=-0.5, μ=0.0)   # Matsubara 频点上的 G(iω_n)

# ED 参考：离散化浴后对角化
db = discretebath(Fermion, -1:0.02:1, f; β=β, μ=0.0)
m = Toulouse(db, ϵ_d=-0.5)
gτ_ed = toulouse_Gτ(m, range(0, β; length=100))
```

## 独立玻色子模型

费米杂质耦合声子（独立谐振子）浴：

```math
H = \epsilon_d d^\dagger d + \sum_k \omega_k a_k^\dagger a_k
+ d^\dagger d \sum_k g_k (a_k + a_k^\dagger),
```

`bands = 2` 时为自旋双带版本，含相互作用 ``U``：

```math
H = \epsilon_d \sum_\sigma d_\sigma^\dagger d_\sigma
+ U\, d_\uparrow^\dagger d_\uparrow d_\downarrow^\dagger d_\downarrow
+ \sum_k \omega_k a_k^\dagger a_k
+ \sum_\sigma d_\sigma^\dagger d_\sigma \sum_k g_k (a_k + a_k^\dagger).
```

声子浴对电子的 dressing 归结为一个乘性指数因子，其中重整化能移

```math
\Delta = \int \mathrm{d}\omega\, \frac{J(\omega)}{\omega}
```

由谱密度自动计算（关键字 `Δ`，默认按上式求积）；极化子指数

```math
\Phi(\tau) = \int \mathrm{d}\omega\, J(\omega)\,
\frac{1 - e^{-\omega\tau} + e^{-\beta\omega} - e^{-(\beta-\tau)\omega}}
{(1 - e^{-\beta\omega})\,\omega^2}.
```

格林函数为（重整化后的）自由费米子因子与 dressing 因子之积：

```math
G(\tau) = G^0(\tau;\, \mu' = -\epsilon_d + \Delta)\; e^{-\Phi(\tau)}
\qquad (\text{bands}=1),
```

```math
G(\tau) = G^0_U(\tau;\, \mu' = -\epsilon_d + \Delta,\, U' = U - 2\Delta)\; e^{-\Phi(\tau)}
\qquad (\text{bands}=2),
```

其中 ``G^0``、``G^0_U`` 分别为上文自由/相互作用费米子的 ``G(\tau)``。
实时间 greater / lesser 同理，指数取复时间 ``\tau = \pm it``：

```math
G^>(t) = G^>_0(t;\, \mu')\, e^{-\Phi(it)},
\qquad
G^<(t) = G^<_0(t;\, \mu')\, e^{-\Phi(-it)}.
```

```@docs
independentbosons_Gτ
independentbosons_greater
independentbosons_lesser
```

## 自旋-玻色退相干

自旋-玻色模型去掉隧穿项（无 ``\sigma_x`` 耦合）后的精确退相干动力学：

```math
H = \Delta \sigma_z + \sigma_z \sum_k V_k (a_k + a_k^\dagger) + \sum_k \omega_k a_k^\dagger a_k.
```

由于只有 ``\sigma_z`` 与浴耦合，布居（对角元）不随时间变化，只有相干（非对角元）
获得重整化相位并衰减：

```math
\rho_{12}(t) = \rho_{12}(0)\, e^{-i\Delta t}
\exp\left[-\int \mathrm{d}\omega\, J(\omega)\,
\frac{\coth(\beta\omega/2)}{\omega^2}\,(1 - \cos\omega t)\right].
```

```@docs
spinboson_dephasingdynamics
```

### 动力学去耦

`ddxx_spinboson_dephasingdynamics` 在 XX 脉冲序列（`N` 为偶数步、每步间隔 `δt`）
下演化退相干模型，可用于演示动力学去耦对相干性的保护。脉冲序列把衰减指数中的
``(1-\cos\omega t)`` 替换为滤波函数 ``F_N(\omega)``：

```math
\rho_{12}(N\delta t) = \rho_{12}(0)\,
\exp\left[-\int \mathrm{d}\omega\, J(\omega)\,
\frac{\coth(\beta\omega/2)}{\omega^2}\,\big(1-\cos\omega\delta t\big)\,
F_N(\omega)\right],
```

其中滤波函数为如下的双重求和（与实现逐项对应，``(-1)^{j+k}`` 为 XX 脉冲的符号调制）：

```math
F_N(\omega) = \sum_{j=1}^{N}\left[1 + 2\sum_{k=1}^{j-1} (-1)^{j+k}
\cos\big(\omega (j-k)\delta t\big)\right].
```

```@docs
ddxx_spinboson_dephasingdynamics
```

## Holstein 模型（实验性）

Holstein 模型（电子-声子耦合格点模型）的零温精确解（连分数）与有限温近似解
当前已实现于 `src/analyticsolutions/holstein/`，但导出暂被注释，
API 尚未稳定，此处不做详细介绍。相关无量纲参数转换
（`λ = g²/(ωt)`、`γ = ω/t`）：

```@docs
ImpurityModelBase.holstein_scaleless_parameters
ImpurityModelBase.holstein_bare_parameters
ImpurityModelBase.GreenFunction
ImpurityModelBase.holstein_G0w_to_Gw_zeroT
ImpurityModelBase.holstein_G0w_to_Gw_finiteT
ImpurityModelBase.holstein_G0w_to_Σw_zeroT
```
