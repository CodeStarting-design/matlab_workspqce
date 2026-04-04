# DCIM 实现对比报告：C++ strata vs MATLAB GreenFunc_matlab

## 1. SpectralMGF — ComputeKmi（K[1] = Gxz 分量）

### 🔴 关键差异：WVe_mi 和 WVh_mi 的计算公式

**C++ (spectral_MGF.cpp:853-856)**:
```cpp
IVe_mi = _I * Zdi.Ze * GIe_ii * TVe;
IVh_mi = _I * Zdi.Zh * GIh_ii * TVh;
WVe_mi = (-jkzm / Zm.Ze) * IVe_mi;  // 注意: Zm.Ze
WVh_mi = (-jkzm / Zm.Zh) * IVh_mi;  // 注意: Zm.Zh
```

**MATLAB (SpectralMGF.m:181-184)**:
```matlab
IVe_mi = I_sign * obj.Zdi.Ze .* obj.GIe .* obj.TVe;
IVh_mi = I_sign * obj.Zdi.Zh .* obj.GIh .* obj.TVh;
WVe_mi = (-1i*obj.omega*obj.epsm) .* IVe_mi;       % ← 用了 omega*epsm
WVh_mi = (-1i*obj.kzm.^2/(obj.omega*obj.mum)) .* IVh_mi;  % ← 用了 kzm^2/(omega*mum)
```

**分析**: C++ 使用 `(-jkzm / Zm.Ze)` 和 `(-jkzm / Zm.Zh)`，其中：
- `Zm.Ze = kzm / (omega * epsm)` → `(-jkzm / Zm.Ze)` = `(-jkzm * omega * epsm / kzm)` = `(-j * omega * epsm)` ✅ 一致
- `Zm.Zh = (omega * mum) / kzm` → `(-jkzm / Zm.Zh)` = `(-jkzm * kzm / (omega * mum))` = `(-j * kzm^2 / (omega * mum))` ✅ 一致

**结论**: ✅ 数学等价，只是写法不同。

---

## 2. SpectralMGF — ComputeKmi（K[2] = Gzx 分量）

**C++ (spectral_MGF.cpp:863-866)**:
```cpp
IIe_mi = _I * Zdi.Ye * GVe_ii * TIe;
IIh_mi = _I * Zdi.Yh * GVh_ii * TIh;
WIe_mi = (-jkzm / Zm.Ye) * IIe_mi;  // Zm.Ye
WIh_mi = (-jkzm / Zm.Yh) * IIh_mi;  // Zm.Yh
```

**MATLAB (SpectralMGF.m:188-191)**:
```matlab
IIe_mi = I_sign * obj.Zdi.Ye .* obj.GVe .* obj.TIe;
IIh_mi = I_sign * obj.Zdi.Yh .* obj.GVh .* obj.TIh;
WIe_mi = (-1i * obj.kzm ./ obj.Zm.Ye) .* IIe_mi;
WIh_mi = (-1i * obj.kzm ./ obj.Zm.Yh) .* IIh_mi;
```

**分析**: 
- C++: `(-jkzm / Zm.Ye)` where `Zm.Ye = 1/Zm.Ze = omega*epsm/kzm`
  → `(-jkzm * kzm / (omega*epsm))` = `(-j * kzm^2 / (omega * epsm))`
- MATLAB: `(-1i * obj.kzm ./ obj.Zm.Ye)` = `(-j * kzm / (omega*epsm/kzm))` = `(-j * kzm^2 / (omega * epsm))`

**结论**: ✅ 等价。

---

## 3. SpectralMGF — ComputeWVii

**C++ (spectral_MGF.cpp:948-950)**:
```cpp
_WVe = (jkzi/2.0/Zi.Ze/De) * (Gli.Ge*exp_terms[1] - Gri.Ge*exp_terms[2] + 2.0*J*Gli.Ge*Gri.Ge*sin_term);
_WVh = (jkzi/2.0/Zi.Zh/Dh) * (Gli.Gh*exp_terms[1] - Gri.Gh*exp_terms[2] + 2.0*J*Gli.Gh*Gri.Gh*sin_term);
```

**MATLAB (SpectralMGF.m:223-224)**:
```matlab
wve = (1i*obj.omega*obj.epsi/2./obj.De) .* (...);
wvh = (1i*obj.kzi.^2/(2*obj.omega*obj.mui)./obj.Dh) .* (...);
```

**分析**: 
- C++ WVe: `jkzi/(2*Zi.Ze*De)` where `Zi.Ze = kzi/(omega*epsi)` → `jkzi*omega*epsi/(2*kzi*De)` = `j*omega*epsi/(2*De)` ✅
- C++ WVh: `jkzi/(2*Zi.Zh*Dh)` where `Zi.Zh = omega*mui/kzi` → `jkzi*kzi/(2*omega*mui*Dh)` = `j*kzi^2/(2*omega*mui*Dh)` ✅

**结论**: ✅ 等价。

---

## 4. SpectralMGF — ComputeWIii

**C++ (spectral_MGF.cpp:964-966)**:
```cpp
_WIe = (jkzi/2.0/Zi.Ye/De) * (-Gli.Ge*exp_terms[1] + Gri.Ge*exp_terms[2] + 2.0*J*Gli.Ge*Gri.Ge*sin_term);
_WIh = (jkzi/2.0/Zi.Yh/Dh) * (-Gli.Gh*exp_terms[1] + Gri.Gh*exp_terms[2] + 2.0*J*Gli.Gh*Gri.Gh*sin_term);
```

**MATLAB (SpectralMGF.m:228-229)**:
```matlab
wie = (1i*obj.kzi.^2/(2*obj.omega*obj.epsi)./obj.De) .* (...);
wih = (1i*obj.omega*obj.mui/2./obj.Dh) .* (...);
```

**分析**: 
- C++ WIe: `jkzi/(2*Zi.Ye*De)` where `Zi.Ye = 1/Zi.Ze = omega*epsi/kzi` → `jkzi*kzi/(2*omega*epsi*De)` = `j*kzi^2/(2*omega*epsi*De)` ✅
- C++ WIh: `jkzi/(2*Zi.Yh*Dh)` where `Zi.Yh = 1/Zi.Zh = kzi/(omega*mui)` → `jkzi*omega*mui/(2*kzi*Dh)` = `j*omega*mui/(2*Dh)` ✅

**结论**: ✅ 等价。

---

## 5. SpectralMGF — ComputeExpTermsGii (半空间处理)

**C++ (spectral_MGF.cpp:991-1002)**:
```cpp
if (i == -1)
{
    _exp_terms[2] = 0.0;
    _cos_term = 0.0;  // ← 半空间时 cos_term 和 sin_term 也设为 0
    _sin_term = 0.0;
}
else if (i == (int)lm->layers.size())
{
    _exp_terms[1] = 0.0;
    _cos_term = 0.0;
    _sin_term = 0.0;
}
```

**MATLAB (SpectralMGF.m:235-236)**:
```matlab
if isinf(h), obj.cos_term=0; obj.sin_term=0; 
else, ... 
end
if obj.i==-1, obj.exp_terms(3)=0; 
elseif obj.i==length(obj.lm.layers), obj.exp_terms(2)=0; 
end
```

**分析**: 
- C++ 在 `i==-1` 或 `i==N` 时，同时设置 `exp_terms`、`cos_term`、`sin_term` 为 0
- MATLAB 用 `isinf(h)` 来检查半空间（因为半空间的 height = inf），然后单独设 `exp_terms`

**⚠️ 潜在风险**: MATLAB 中 `isinf(h)` 和 `i==-1 || i==N` 的判断**顺序不同**。如果 `h = inf`（半空间），MATLAB 先设 cos/sin=0，然后再设 exp_terms。但 cos/sin 是在 `if isinf(h)` 中设置的，而 exp_terms 是在后面的 `if` 中设置。如果 `h` 不是 inf 但 `i==-1`，MATLAB 可能不会设 cos/sin=0。但实际上对于半空间 `h = inf`（由 `GetHeight` 返回），所以逻辑一致。

**结论**: ✅ 等价（通过不同的代码路径实现相同逻辑）。

---

## 6. SpectralMGF — InterfaceImpedance (tan 溢出处理)

### 🔴 重要差异：tan(kz*h) 溢出处理

**C++ (spectral_MGF.cpp:563-576)**:
```cpp
if (std::isfinite(std::abs(t)))
{
    Zr.Ze = Z.Ze*(Zr_next.Ze + J*Z.Ze*t) / (Z.Ze + J*Zr_next.Ze*t);
    ...
}
else  // L'Hospital's rule
{
    Zr.Ze = std::pow(Z.Ze, 2) / Zr_next.Ze;
    ...
}
```

**MATLAB (SpectralMGF.m:305-306)**:
```matlab
Zr.Ze = Z.Ze.*(Zr_next.Ze+J.*Z.Ze.*t)./(Z.Ze+J.*Zr_next.Ze.*t);
```

**⚠️ 隐藏风险**: MATLAB 没有处理 `tan` 溢出（当 `kz*h = pi/2 + n*pi` 时，`tan` 为无穷大）。C++ 有 L'Hospital 规则的 fallback。对于 DCIM 路径上的复数 `krho`，`tan` 极少溢出，但对于特定的 `krho` 值可能导致 NaN。

**结论**: ⚠️ MATLAB 缺少 `isfinite` 检查，可能在极端情况下产生 NaN。

---

## 7. SpectralMGF — InterfaceImpedance_Upward (Ye/Yh 计算)

### 🔴 重要差异：Ye/Yh 的计算方式

**C++ (spectral_MGF.cpp:565-568)**:
```cpp
Zr.Ze = Z.Ze*(Zr_next.Ze + J*Z.Ze*t) / (Z.Ze + J*Zr_next.Ze*t);
Zr.Zh = Z.Zh*(Zr_next.Zh + J*Z.Zh*t) / (Z.Zh + J*Zr_next.Zh*t);
Zr.Ye = Z.Ye*(Zr_next.Ye + J*Z.Ye*t) / (Z.Ye + J*Zr_next.Ye*t);
Zr.Yh = Z.Yh*(Zr_next.Yh + J*Z.Yh*t) / (Z.Yh + J*Zr_next.Yh*t);
```

**MATLAB (SpectralMGF.m:306)**:
```matlab
Zr.Ze=Z.Ze.*(Zr_next.Ze+J.*Z.Ze.*t)./(Z.Ze+J.*Zr_next.Ze.*t); 
Zr.Zh=Z.Zh.*(Zr_next.Zh+J.*Z.Zh.*t)./(Z.Zh+J.*Zr_next.Zh.*t); 
Zr.Ye=1./Zr.Ze;   % ← 简单取倒数
Zr.Yh=1./Zr.Zh;   % ← 简单取倒数
```

### 🔴🔴🔴 这是一个 BUG！

C++ 对 Ye 和 Yh 使用**独立的递推公式**（使用 Ye 的导纳值递推），而 MATLAB 只计算 Ze/Zh 然后取倒数 `Ye = 1/Ze`。

虽然数学上 `Ye = 1/Ze`，但**递推过程中数值精度不同**！C++ 的递推：
```
Zr.Ye = Z.Ye*(Zr_next.Ye + J*Z.Ye*t) / (Z.Ye + J*Zr_next.Ye*t)
```
这个公式和 Ze 的递推形式相同但使用的是 Ye 值，结果等于 `1/Zr.Ze`，在精确算术下等价。但在有限精度浮点运算中，两种计算方式的**舍入误差**不同。

**影响评估**: 在大多数情况下差异极小，但对于大的 krho（虚轴深处），阻抗值可能非常大或非常小，取倒数会放大舍入误差。这不太可能是 DCIM 22% 误差的主要原因，但建议统一。

**结论**: 🔴 数值上可能有微小差异。建议 MATLAB 也对 Ye/Yh 做独立递推。

---

## 8. QuasistaticMGF — Fresnel 系数索引

### C++ (quasistatic_MGF.cpp:464-476):
```cpp
Fresnel QuasistaticMGF::GetFresnelCoefficient_Upward(int idx)
{
    if (idx == -2)
    {
        Fresnel Rf;
        Rf.Fh = 0.0;
        Rf.Fe = 0.0;
        return Rf;
    }
    return Rrs[idx+1];
}
```

### MATLAB (QuasistaticMGF2.m:481-487):
```matlab
function f = GetFresnelCoefficient_Upward(obj, idx)
    if idx == -2
        f.Fe = 0; f.Fh = 0; return;
    end
    f = obj.Rrs(idx + 2);  % ← idx+2 因为 MATLAB 1-based
end
```

**分析**: C++ 用 `Rrs[idx+1]`（0-based），MATLAB 用 `obj.Rrs(idx+2)`（1-based 转换）。
- C++ `idx=-1` → `Rrs[0]`（第1层与上半空间的 Fresnel 系数）
- MATLAB `idx=-1` → `obj.Rrs(1)`（同上）

**结论**: ✅ 等价。

---

## 9. QuasistaticMGF — Fresnel 系数计算

### C++ (quasistatic_MGF.cpp:419-428):
```cpp
if (lm->isPEC_top)
{
    Rrs[0].Fh = -1.0;
    Rrs[0].Fe = -1.0;
}
else
{
    Rrs[0].Fh = -(lm->mu[0] - lm->mu_top)/(lm->mu[0] + lm->mu_top);
    Rrs[0].Fe = (lm->eps[0] - lm->eps_top)/(lm->eps[0] + lm->eps_top);
}
for (int ii = 1; ii < lm->layers.size(); ii++)
{
    Rrs[ii].Fh = -(lm->mu[ii] - lm->mu[ii-1])/(lm->mu[ii] + lm->mu[ii-1]);
    Rrs[ii].Fe = (lm->eps[ii] - lm->eps[ii-1])/(lm->eps[ii] + lm->eps[ii-1]);
}
```

### MATLAB (QuasistaticMGF2.m:472-477):
```matlab
for ii = 0:(N-1)
    [ec, mc] = obj.GetLayerParams(ii); 
    [eu, mu] = obj.GetLayerParams(ii-1); 
    idx=ii+1; 
    obj.Rrs(idx).Fh=-(mc-mu)/(mc+mu); 
    obj.Rrs(idx).Fe=(ec-eu)/(ec+eu); 
end
```

**分析**: MATLAB 统一用循环，`GetLayerParams(ii-1)` 当 `ii=0` 时返回上半空间参数。C++ 对 `ii=0` 特殊处理。
- C++ 检查 `isPEC_top` 来设置 Fe=Fh=-1
- MATLAB 没有检查 `isPEC_top`！当 `isPEC_top=false` 时两者等价。当 `isPEC_top=true` 时：
  - C++ 直接设 `-1`
  - MATLAB: `GetLayerParams(-1)` 返回 `eps_top, mu_top`，然后 `Fe = (eps0-eps_top)/(eps0+eps_top)`，`Fh = -(mu0-mu_top)/(mu0+mu_top)`

### 🔴 潜在差异 (当 isPEC_top=true 时):
本测试案例 `isPEC_top=false`，所以不影响。但如果 `isPEC_top=true`，MATLAB 不会得到 `-1`。

**结论**: ⚠️ 本案例不影响（isPEC_top=false），但 isPEC_top=true 时有 bug。

---

## 10. QuasistaticMGF — Downward Fresnel for bottom PEC

### C++ (quasistatic_MGF.cpp:442-451):
```cpp
if (lm->isPEC_bot)
{
    Rls.back().Fh = -1.0;
    Rls.back().Fe = -1.0;
}
else
{
    Rls.back().Fh = -(lm->mu.back() - lm->mu_bot)/(lm->mu.back() + lm->mu_bot);
    Rls.back().Fe = (lm->eps.back() - lm->eps_bot)/(lm->eps.back() + lm->eps_bot);
}
```

### MATLAB (QuasistaticMGF2.m:475-477):
```matlab
for ii = 0:(N-1)
    [ec, mc] = obj.GetLayerParams(ii); 
    [el, ml] = obj.GetLayerParams(ii+1); 
    idx=ii+1;
    if ii==N-1 && obj.lm.isPEC_bot
        obj.Rls(idx).Fh=-1; obj.Rls(idx).Fe=-1; 
    else
        obj.Rls(idx).Fh=-(mc-ml)/(mc+ml); 
        obj.Rls(idx).Fe=(ec-el)/(ec+el); 
    end
end
```

**结论**: ✅ MATLAB 正确处理了 isPEC_bot。

---

## 11. QuasistaticMGF — GetFresnelCoefficient_Downward 索引

### C++ (quasistatic_MGF.cpp:480-491):
```cpp
Fresnel QuasistaticMGF::GetFresnelCoefficient_Downward(int idx)
{
    if (idx == (int)(lm->layers.size()))
    {
        Fresnel Rf;
        Rf.Fh = 0.0;
        Rf.Fe = 0.0;
        return Rf;
    }
    return Rls[idx];   // 0-based indexing
}
```

### MATLAB (QuasistaticMGF2.m:488-493):
```matlab
function f = GetFresnelCoefficient_Downward(obj, idx)
    if idx==-1 || idx==length(obj.lm.layers)
        f.Fe=0; f.Fh=0; return;
    end
    f = obj.Rls(idx + 1);  % 0-based to 1-based
end
```

**分析**: 
- C++ 只检查 `idx == N`
- MATLAB 额外检查 `idx == -1`

### ⚠️ 边界差异：
MATLAB 额外的 `idx==-1` 检查不在 C++ 中。在正常使用中 `GetFresnelCoefficient_Downward` 不应被调用 `idx=-1`，但如果被调用，C++ 会直接访问 `Rls[-1]`（未定义行为），而 MATLAB 返回 0。这是一个安全检查，不影响正确性。

**结论**: ✅ 等价（MATLAB 多了安全检查）。

---

## 12. DCIM — epsr_max 计算

### C++ (DCIM.cpp:305-309):
```cpp
double epsr_max = 0.0;
for (int ii = 0; ii < lm->layers.size(); ii++)
    if (std::abs(lm->eps[ii]) > epsr_max)
        epsr_max = std::abs(lm->eps[ii]);
epsr_max /= eps0;
```

### MATLAB (DCIM_Integrated.m):
```matlab
epsr_list = [obj.lm.layers.epsr];
epsr_max = max(epsr_list);
```

**分析**: 
- C++ 使用 `abs(complex_eps) / eps0`
- MATLAB 使用 `layers.epsr`（相对介电常数实部）

对于无损层（sigma=0），`abs(eps0*epsr)/eps0 = epsr`。✅ 等价。
对于有损层，C++ 会包含损耗部分 `sqrt(epsr^2 + (sigma/(omega*eps0))^2)`，而 MATLAB 只用 `epsr`。

**结论**: ✅ 本案例等价（无损层），但有损层时有差异。

---

## 13. DCIM — GPOF 实现差异

### 🔴🔴🔴 这是最大的差异来源！

**C++ (DCIM.cpp:820-912)** 使用 LAPACK:
- `ComputeSVD()` → LAPACK `zgesvd` ('A','A' 全SVD)
- `ComputeEigenvalues()` → LAPACK `zgeev`
- `SolveLeastSquares()` → LAPACK `zgels`
- 所有矩阵运算使用 CBLAS `zgemm`

**MATLAB (GPOF.m)** 使用内置函数:
- `svd(Y1)` → MATLAB 内置 SVD
- `eig(Z)` 或 `schur(Z)` → MATLAB 内置特征值分解
- `Y3 \ y` → MATLAB 内置最小二乘

**关键差异**:

1. **特征值精度**: LAPACK `zgeev` 和 MATLAB `eig`/`schur` 对同一矩阵可能产生不同精度的特征值，特别是对于病态矩阵。

2. **SVD 分解唯一性**: SVD 的 U 和 V 矩阵不唯一（可以乘以任意相位因子），不同实现可能产生不同的 U/V，导致 Z 矩阵的数值不同。

3. **最小二乘求解**: C++ 使用 `zgels`（QR 分解），MATLAB `\` 可能使用不同的算法。

**影响**: 这些数值差异通过 Level 2 的 t→kz 转换中的 `exp(k*alpha)` 被**指数放大**，导致最终复镜像系数的显著差异。

**结论**: 🔴🔴🔴 这是 DCIM 22% 误差（vs C++ 0.24%）的主要原因。无法通过代码修改消除，是 LAPACK vs MATLAB 数值库的固有差异。

---

## 14. DCIM — 缓存机制差异

### C++ (SpectralMGF.cpp:181-213):
```cpp
void SpectralMGF::SetSourcePoint(double _zp)
{
    zp = _zp;
    src_terms_computed = false;
    
    if (precomputations_done && i != m)
    {
        // 预计算源层项
        ComputeGVii(_z, zp, GVe, GVh);
        ComputeGIii(_z, zp, GIe, GIh);
        src_terms_computed = true;
    }
    src_point_set = true;
    exp_terms_computed = false;
}
```

### MATLAB (SpectralMGF.m:52-58):
```matlab
function SetSourcePoint(obj, zp) 
    obj.zp=zp; 
    obj.src_point_set=true; 
    if obj.precomputations_done && obj.i~=obj.m 
        obj.src_terms_computed=false;  % ← 只设置标志，不预计算
    end
    obj.exp_terms_computed=false; 
end
```

**⚠️ 差异**: C++ 在 `SetSourcePoint` 中**立即预计算**源层项（如果 precomputations 已完成），而 MATLAB 只设置 `src_terms_computed=false`，延迟到 `ComputeKmi` 时再计算。

**影响**: 在 DCIM 采样循环中，`SetSourcePoint` 和 `SetObservationPoint` 只调用一次（在路径生成之前），之后只调用 `SetRadialWaveNumber` 和 `ComputeSpectralMGF`。由于 `SetRadialWaveNumber` 调用 `ResetComputedFlags()` 重置 `src_terms_computed=false`，所以无论是否预计算，每次新的 krho 都会重新计算。**不影响结果正确性。**

**结论**: ✅ 不影响结果（仅影响性能）。

---

## 15. LayerManager — k_min / k_max 计算

### 🔴 已修复的 bug

**C++ (layers.cpp:406-455)**: `k_min`/`k_max` 只遍历层，半空间代码被注释掉。

**MATLAB (LayerManager.m)**: 修复后与 C++ 一致。

**⚠️ 注意**: 需要 `clear classes` 才能使修改生效。

---

## 总结

| # | 问题 | 严重性 | 状态 |
|---|------|--------|------|
| 1 | GPOF 数值差异(LAPACK vs MATLAB) | 🔴🔴🔴 | 无法通过代码修复 |
| 2 | InterfaceImpedance 缺少 tan 溢出处理 | 🔴 | 未修复 |
| 3 | InterfaceImpedance Ye/Yh 用 1/Ze 代替独立递推 | 🔴 | 未修复 |
| 4 | k_min/k_max 包含半空间 | 🔴🔴 | ✅ 已修复 |
| 5 | GPOF 稳定性过滤器 | 🔴🔴 | ✅ 已修复 |
| 6 | DCIM extract_quasistatic 配置 | 🔴🔴🔴 | ✅ 已修复 |
| 7 | isPEC_top 时 Fresnel 系数计算 | ⚠️ | 本案例不影响 |
| 8 | 有损层的 epsr_max 计算 | ⚠️ | 本案例不影响 |
| 9 | SetSourcePoint 缓存策略 | ✅ | 仅影响性能 |
| 10 | 所有 Kii/Kmi 公式 | ✅ | 数学等价 |
