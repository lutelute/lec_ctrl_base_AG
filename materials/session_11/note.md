# 第11章：制御系の過渡特性 (詳細版)

本章では、2次遅れ系におけるステップ応答を解析的に解き、「最大行き過ぎ量」や「整定時間」などの評価指標を数式から導出します。
120分の講義を想定し、減衰係数 $\zeta$ による挙動の違い（不足減衰、臨界減衰、過減衰）を詳しく解説します。

---

## 11.1 2次遅れ系のステップ応答

標準的な2次遅れ閉ループ系を考えます。
$$ G(s) = \frac{\omega_n^2}{s^2 + 2\zeta\omega_n s + \omega_n^2} $$
ここで $\omega_n > 0$、$\zeta > 0$ とします。
単位ステップ入力 $U(s) = 1/s$ に対する出力 $Y(s) = G(s) \frac{1}{s}$ を逆ラプラス変換して $y(t)$ を求めます。

特性方程式 $s^2 + 2\zeta\omega_n s + \omega_n^2 = 0$ の根（極）は、
$$ s = -\zeta\omega_n \pm \omega_n \sqrt{\zeta^2 - 1} $$
$\zeta$ の値によって3つのケースに分類されます。

### 11.1.1 ケース1: 不足減衰 (Underdamped, $0 < \zeta < 1$)
極は複素共役根になります。
$$ s = -\zeta\omega_n \pm j\omega_n \sqrt{1 - \zeta^2} = -\sigma \pm j\omega_d $$
ここで $\sigma = \zeta\omega_n$（減衰係数）、$\omega_d = \omega_n \sqrt{1-\zeta^2}$（減衰固有角周波数）。

部分分数展開と逆変換を行うと（詳細は第2章参照）、以下の応答が得られます。
> $$ y(t) = 1 - \frac{e^{-\zeta\omega_n t}}{\sqrt{1-\zeta^2}} \sin(\omega_d t + \phi), \quad \phi = \cos^{-1}\zeta $$

- 特徴: 振動しながら目標値 (1) に収束します。
- $\zeta$ が小さいほど減衰が遅く、激しく振動します。

### 11.1.2 ケース2: 臨界減衰 (Critically Damped, $\zeta = 1$)
極は重根 $s = -\omega_n$ です。
$$ Y(s) = \frac{\omega_n^2}{s(s+\omega_n)^2} = \frac{1}{s} - \frac{1}{s+\omega_n} - \frac{\omega_n}{(s+\omega_n)^2} $$
逆変換すると、
> $$ y(t) = 1 - e^{-\omega_n t} (1 + \omega_n t) $$

- 特徴: 振動せずに最も速く収束する境界条件です。

### 11.1.3 ケース3: 過減衰 (Overdamped, $\zeta > 1$)
極は異なる2つの実根 $s = -\omega_n (\zeta \pm \sqrt{\zeta^2-1})$ です。
これを $p_1, p_2$ とすると、
> $$ y(t) = 1 + \frac{p_2}{p_1 - p_2} e^{p_1 t} - \frac{p_1}{p_1 - p_2} e^{p_2 t} $$

- 特徴: 振動しませんが、収束は遅くなります。

---

## 11.2 過渡特性の評価指標（不足減衰の場合）

制御系設計では、振動的な応答（ケース1）の形状を評価するために以下の指標を用います。

### 11.2.1 最大行き過ぎ量 (Max Overshoot) $M_p$
出力が最大値に達する時刻（行き過ぎ時間 $T_p$）における、目標値からの超過分です。

まず $T_p$ を求めます。$\frac{dy}{dt} = 0$ となる最小の $t > 0$ です。
$$ \frac{dy}{dt} = \frac{\omega_n}{\sqrt{1-\zeta^2}} e^{-\zeta\omega_n t} \sin(\omega_d t) $$
$\sin(\omega_d t) = 0 \implies \omega_d t = \pi$ （最初の極大値）。
よって、
> $$ T_p = \frac{\pi}{\omega_d} = \frac{\pi}{\omega_n \sqrt{1-\zeta^2}} $$

これを $y(t)$ に代入して最大値を求めます。
$$ y(T_p) = 1 + e^{-\frac{\zeta\pi}{\sqrt{1-\zeta^2}}} $$
よってオーバーシュート $M_p$ は、
> $$ M_p = e^{-\frac{\pi \zeta}{\sqrt{1-\zeta^2}}} \times 100 [\%] $$

$M_p$ は $\zeta$ だけで決まり、$\omega_n$ には依存しません。
- $\zeta = 0.5$ のとき $M_p \approx 16.3\%$
- $\zeta = 0.7$ のとき $M_p \approx 4.6\%$

### 11.2.2 整定時間 (Settling Time) $T_s$
出力が目標値の $\pm 5\%$ 以内に収まるまでの時間です。
包絡線 $1 \pm \frac{e^{-\zeta\omega_n t}}{\sqrt{1-\zeta^2}}$ が $1.05$ (または $0.95$) になる時間を考えます。
近似的に時定数 $\tau = 1/(\zeta\omega_n)$ の3倍の時間とされます（$e^{-3} \approx 0.05$）。

> $$ T_s \approx \frac{3}{\zeta \omega_n} $$

- 応答を速くするには、実部 $\zeta\omega_n$ を大きくする必要があります。

### 11.2.3 立ち上がり時間 (Rise Time) $T_r$
出力が初めて目標値 (1.0) に到達する時間、または 10% から 90% に上昇する時間です。
定義によりますが、0から100%への到達時間は概算で
> $$ T_r \approx \frac{1.8}{\omega_n} $$

---

## 11.3 演習問題

**問題**:
2次遅れ系において、最大行き過ぎ量を 5% 以下に抑え、かつ整定時間（5%基準）を 1秒以内にしたい。
必要な $\zeta$ と $\omega_n$ の範囲を求めよ。

**解答**:
1.  **行き過ぎ量条件**:
    $M_p < 0.05 \implies e^{-\frac{\pi \zeta}{\sqrt{1-\zeta^2}}} < 0.05$
    両辺の自然対数をとる: $-\frac{\pi \zeta}{\sqrt{1-\zeta^2}} < \ln(0.05) \approx -3$
    $\frac{\pi \zeta}{\sqrt{1-\zeta^2}} > 3 \implies \pi^2 \zeta^2 > 9(1-\zeta^2) \implies (\pi^2 + 9)\zeta^2 > 9$
    $\zeta^2 > \frac{9}{9.86 + 9} \approx 0.47$
    $\zeta > 0.69$
    余裕を見て $\mathbf{\zeta \ge 0.7}$ とすれば十分である。

2.  **整定時間条件**:
    $T_s \approx \frac{3}{\zeta \omega_n} < 1$
    これより $\zeta \omega_n > 3$。
    もし $\zeta = 0.7$ を選ぶなら、
    $0.7 \omega_n > 3 \implies \omega_n > 4.3$ [rad/s]

**設計例**: $\zeta = 0.7, \omega_n = 5$ とすれば要求を満たす。
伝達関数は $G(s) = \frac{25}{s^2 + 7s + 25}$ となる。
