# 第5章：周波数伝達関数と周波数応答 (詳細版)

本章では、制御系の特性を知る上で重要となる「周波数入力に対する定常的な応答」について学びます。
120分の講義を想定し、理論的背景（なぜ $s=j\omega$ なのか）から、ベクトル軌跡、そして実用的なボード線図の描き方まで詳細に解説します。

---

## 5.1 周波数応答の基礎理論

### 5.1.1 定義
安定なLTIシステムに正弦波入力 $u(t) = A \sin(\omega t)$ を加えたとき、十分時間が経過した後の定常的な出力 $y_{ss}(t)$ は、同じ周波数の正弦波になりますが、**振幅** と **位相** が変化します。

$$ y_{ss}(t) = M A \sin(\omega t + \phi) $$

この「振幅倍率 $M(\omega)$」と「位相差 $\phi(\omega)$」の特性を **周波数応答** と呼びます。

### 5.1.2 ラプラス変換による導出（$s=j\omega$ の意味）
なぜ伝達関数 $G(s)$ に $s=j\omega$ を代入すると周波数特性が得られるのか、数学的に証明します。

**【証明】**
安定なシステム $G(s)$ に正弦波 $u(t) = \sin\omega t$ を加えます。
$$ U(s) = \frac{\omega}{s^2 + \omega^2} = \frac{\omega}{(s-j\omega)(s+j\omega)} $$
出力 $Y(s)$ は、
$$ Y(s) = G(s) U(s) = G(s) \frac{\omega}{(s-j\omega)(s+j\omega)} $$

これを部分分数展開します。$G(s)$ の極を $p_1, \dots, p_n$ とすると、
$$ Y(s) = \frac{a}{s-j\omega} + \frac{\bar{a}}{s+j\omega} + \sum_{k=1}^{n} \frac{C_k}{s-p_k} $$
ここで、システムは安定（全ての極 $p_k$ の実部が負）であるため、逆ラプラス変換したとき、$\sum$ の項は時間とともに減衰してゼロになります（過渡応答）。
定常状態（$t \to \infty$）で残るのは最初の2項だけです。

係数 $a$ を求めます（留数定理）。
$$ a = \left. (s-j\omega) Y(s) \right|_{s=j\omega} = \left. G(s) \frac{\omega}{s+j\omega} \right|_{s=j\omega} = G(j\omega) \frac{\omega}{2j\omega} = \frac{1}{2j} G(j\omega) $$
$\bar{a}$ は $a$ の共役なので、
$$ \bar{a} = -\frac{1}{2j} G(-j\omega) $$
ただし実係数システムでは $G(-j\omega) = \overline{G(j\omega)}$（共役）です。

$G(j\omega)$ を極形式で表します。
$$ G(j\omega) = |G(j\omega)| e^{j\angle G(j\omega)} = M e^{j\phi} $$
すると、
$$ a = \frac{1}{2j} M e^{j\phi}, \quad \bar{a} = -\frac{1}{2j} M e^{-j\phi} $$
定常出力 $y_{ss}(t)$ は、
$$
\begin{aligned}
y_{ss}(t) &= \mathcal{L}^{-1}\left[ \frac{a}{s-j\omega} + \frac{\bar{a}}{s+j\omega} \right] \\
&= a e^{j\omega t} + \bar{a} e^{-j\omega t} \\
&= \frac{M}{2j} e^{j(\omega t + \phi)} - \frac{M}{2j} e^{-j(\omega t + \phi)} \\
&= M \frac{e^{j(\omega t + \phi)} - e^{-j(\omega t + \phi)}}{2j} \\
&= M \sin(\omega t + \phi) \quad \blacksquare
\end{aligned}
$$

**結論**: 伝達関数 $G(s)$ の $s$ を $j\omega$ に置き換えた複素数 $G(j\omega)$ が、そのまま周波数応答を表します。
- **ゲイン (振幅倍率)**: $M = |G(j\omega)|$
- **位相差**: $\phi = \angle G(j\omega)$

---

## 5.2 ベクトル軌跡 (Nyquist Plot)

複素平面上で、$\omega$ を $0$ から $\infty$ まで変化させたときの $G(j\omega)$ の軌跡です。

### 5.2.1 1次遅れ系
$$ G(j\omega) = \frac{K}{1 + j\omega T} $$
分母・分子に共役複素数 $1-j\omega T$ を掛けて有理化します。
$$ G(j\omega) = \frac{K(1 - j\omega T)}{1 + (\omega T)^2} = \frac{K}{1+(\omega T)^2} - j \frac{K\omega T}{1+(\omega T)^2} $$
実部 $u$、虚部 $v$ とすると、
$$ u^2 + v^2 = \frac{K^2}{(1+(\omega T)^2)^2} (1 + (\omega T)^2) = \frac{K^2}{1+(\omega T)^2} = K u $$
$$ (u - K/2)^2 + v^2 = (K/2)^2 $$
これは中心 $(K/2, 0)$、半径 $K/2$ の円（の下半分）を描きます。

### 5.2.2 積分系
$$ G(j\omega) = \frac{1}{j\omega} = -j \frac{1}{\omega} $$
実部は常に0、虚部は $-\infty$ から $0$ へ上昇します（負の虚軸上）。

---

## 5.3 ボード線図 (Bode Diagram)

ベクトル軌跡は全体の形状把握に優れますが、周波数ごとの特性を読み取るには不便です。
そこで、**対数グラフ** を用いたボード線図が実務で多用されます。

### 5.3.1 定義
2つのグラフを並べて描きます。
1.  **ゲイン線図**: 横軸 $\log \omega$、縦軸 **ゲイン [dB]**
    $$ g_{dB} = 20 \log_{10} |G(j\omega)| $$
2.  **位相線図**: 横軸 $\log \omega$、縦軸 **位相 [deg]**
    $$ \phi = \angle G(j\omega) \times \frac{180}{\pi} $$

### 5.3.2 1次遅れ系のボード線図
$$ G(j\omega) = \frac{1}{1 + j\omega T} $$

#### ゲイン特性
$$ |G(j\omega)| = \frac{1}{\sqrt{1 + (\omega T)^2}} $$
$$ g_{dB} = -20 \log_{10} \sqrt{1 + (\omega T)^2} = -10 \log_{10} (1 + (\omega T)^2) $$

- **低周波域 ($\omega T \ll 1$)**: $1 + (\omega T)^2 \approx 1$ なので、$g_{dB} \approx 0$ dB。
- **高周波域 ($\omega T \gg 1$)**: $1 + (\omega T)^2 \approx (\omega T)^2$ なので、
  $$ g_{dB} \approx -10 \log_{10} (\omega T)^2 = -20 \log_{10} \omega - 20 \log_{10} T $$
  これは周波数が10倍（1 decade）になるごとに -20 dB 減少する直線です（傾き **-20 dB/dec**）。

2つの直線の交点 $\omega = 1/T$ を **折れ点周波数 (Corner Frequency)** と呼びます。

#### 位相特性
$$ \angle G(j\omega) = - \tan^{-1}(\omega T) $$
- $\omega \to 0$ で $0^\circ$
- $\omega = 1/T$ で $-45^\circ$
- $\omega \to \infty$ で $-90^\circ$

### 5.3.3 デシベル (dB) の感覚
- $20$ dB: 10倍
- $40$ dB: 100倍
- $6$ dB: 約2倍 ($20 \log 2 \approx 6.02$)
- $-3$ dB: 約 $1/\sqrt{2} \approx 0.707$ 倍（エネルギーが半分になる点）

---

## 5.4 演習問題

**問題**: 伝達関数 $G(s) = \frac{10}{s + 10}$ のボード線図の概略を描け。

**解答**:
まず標準形 $\frac{K}{Ts+1}$ に変形します。分母・分子を10で割ります。
$$ G(s) = \frac{1}{0.1s + 1} $$
これより、
- DCゲイン: $K=1 \implies 0$ dB
- 時定数: $T = 0.1$ [s]
- 折れ点周波数: $\omega_c = 1/T = 10$ [rad/s]

**描画手順**:
1.  ゲイン線図: $\omega < 10$ では $0$ dB の直線。$\omega = 10$ から $-20$ dB/dec で下がる直線を引く。（$\omega=100$ で $-20$ dB）
2.  位相線図: $\omega < 1$ で $0^\circ$、$\omega=10$ で $-45^\circ$、$\omega > 100$ で $-90^\circ$ を通る滑らかな曲線を描く。
