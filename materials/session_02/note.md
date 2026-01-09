# 第2章：複素数とラプラス変換 (詳細版)

本章では、制御工学の基礎言語である「ラプラス変換」について、その定義から性質、具体的な計算方法まで厳密に解説します。

---

## 2.1 複素数の基礎

制御工学では、信号やシステムの特性を表現するために複素数を多用します。

### 2.1.1 虚数単位と複素平面
虚数単位 $j$ を $j = \sqrt{-1}$ （すなわち $j^2 = -1$）と定義します。
複素数 $s$ は、実部 $\sigma$ と虚部 $\omega$ を用いて以下のように表されます（直交形式）。
$$ s = \sigma + j\omega $$

これを **複素平面（$s$平面）** 上の点 $(\sigma, \omega)$ として図示できます。
- 横軸：実軸 (Real Axis, Re)
- 縦軸：虚軸 (Imaginary Axis, Im)

### 2.1.2 極形式とオイラーの公式
ある複素数 $s$ の原点からの距離を $r = |s|$、実軸となす角を $\theta = \angle s$ とすると、
$$ s = r(\cos\theta + j\sin\theta) $$
と表せます（極形式）。

ここで、以下の **オイラーの公式 (Euler's Formula)** が極めて重要です。

> **定理：オイラーの公式**
> $$ e^{j\theta} = \cos\theta + j\sin\theta $$

**【証明】マクローリン展開による導出**
指数関数、三角関数のマクローリン展開（$x=0$ 近傍でのテイラー展開）は以下の通りです。
$$
\begin{aligned}
e^x &= 1 + x + \frac{x^2}{2!} + \frac{x^3}{3!} + \frac{x^4}{4!} + \dots \\
\cos x &= 1 - \frac{x^2}{2!} + \frac{x^4}{4!} - \dots \\
\sin x &= x - \frac{x^3}{3!} + \frac{x^5}{5!} - \dots
\end{aligned}
$$
ここで $x = j\theta$ を $e^x$ の式に代入します。$j^2=-1, j^3=-j, j^4=1, \dots$ に注意すると、
$$
\begin{aligned}
e^{j\theta} &= 1 + (j\theta) + \frac{(j\theta)^2}{2!} + \frac{(j\theta)^3}{3!} + \frac{(j\theta)^4}{4!} + \dots \\
&= 1 + j\theta - \frac{\theta^2}{2!} - j\frac{\theta^3}{3!} + \frac{\theta^4}{4!} + \dots \\
&= \left( 1 - \frac{\theta^2}{2!} + \frac{\theta^4}{4!} - \dots \right) + j \left( \theta - \frac{\theta^3}{3!} + \dots \right) \\
&= \cos\theta + j\sin\theta \quad \blacksquare
\end{aligned}
$$

これにより、極形式は $s = r e^{j\theta}$ と簡潔に記述できます。

---

## 2.2 ラプラス変換

### 2.2.1 定義
時間関数 $f(t)$ ($t \ge 0$) に対して、以下の積分変換を **ラプラス変換 (Laplace Transform)** と定義し、$\mathcal{L}[f(t)]$ または $F(s)$ と表記します。

> **定義：ラプラス変換**
> $$ F(s) = \mathcal{L}[f(t)] = \int_{0}^{\infty} f(t) e^{-st} dt $$

ここで $s$ は複素変数 ($s = \sigma + j\omega$) です。この積分が収束する範囲内で $F(s)$ は定義されます。

### 2.2.2 基本関数のラプラス変換（導出付き）

#### (1) 単位ステップ関数 (Unit Step Function) $u(t)$
$$
u(t) = \begin{cases} 1 & (t \ge 0) \\ 0 & (t < 0) \end{cases}
$$
**【導出】**
$$
\mathcal{L}[1] = \int_{0}^{\infty} 1 \cdot e^{-st} dt = \left[ -\frac{1}{s} e^{-st} \right]_{0}^{\infty}
$$
$\text{Re}(s) > 0$ であれば $\lim_{t\to\infty} e^{-st} = 0$ なので、
$$
= 0 - \left( -\frac{1}{s} \cdot 1 \right) = \frac{1}{s}
$$

#### (2) 指数関数 (Exponential Function) $e^{-at}$
**【導出】**
$$
\begin{aligned}
\mathcal{L}[e^{-at}] &= \int_{0}^{\infty} e^{-at} e^{-st} dt = \int_{0}^{\infty} e^{-(s+a)t} dt \\
&= \left[ -\frac{1}{s+a} e^{-(s+a)t} \right]_{0}^{\infty} = \frac{1}{s+a} \quad (\text{Re}(s+a) > 0)
\end{aligned}
$$

#### (3) ランプ関数 (Ramp Function) $t$
**【導出】** 部分積分 $\int u v' dt = uv - \int u' v dt$ を用います。
$$
\begin{aligned}
\mathcal{L}[t] &= \int_{0}^{\infty} t e^{-st} dt \\
&= \left[ t \cdot \left(-\frac{1}{s}e^{-st}\right) \right]_{0}^{\infty} - \int_{0}^{\infty} 1 \cdot \left(-\frac{1}{s}e^{-st}\right) dt
\end{aligned}
$$
第1項は $t\to\infty$ で $0$ に収束します（ロピタルの定理的振る舞い）。
$$
= 0 + \frac{1}{s} \int_{0}^{\infty} e^{-st} dt = \frac{1}{s} \cdot \frac{1}{s} = \frac{1}{s^2}
$$

#### (4) 正弦波・余弦波 ($\sin\omega t, \cos\omega t$)
オイラーの公式より $\cos\omega t = \frac{e^{j\omega t} + e^{-j\omega t}}{2}$, $\sin\omega t = \frac{e^{j\omega t} - e^{-j\omega t}}{2j}$ です。指数関数のラプラス変換の線形性を利用します。

**【導出: $\cos\omega t$】**
$$
\begin{aligned}
\mathcal{L}[\cos\omega t] &= \frac{1}{2} \left( \mathcal{L}[e^{j\omega t}] + \mathcal{L}[e^{-j\omega t}] \right) \\
&= \frac{1}{2} \left( \frac{1}{s - j\omega} + \frac{1}{s + j\omega} \right) \\
&= \frac{1}{2} \frac{(s + j\omega) + (s - j\omega)}{(s - j\omega)(s + j\omega)} = \frac{1}{2} \frac{2s}{s^2 + \omega^2} = \frac{s}{s^2 + \omega^2}
\end{aligned}
$$

同様に、
$$
\mathcal{L}[\sin\omega t] = \frac{\omega}{s^2 + \omega^2}
$$

---

## 2.3 ラプラス変換の重要な性質

これらは制御系の計算で頻繁に使われる「公式」です。証明を理解しておきましょう。

### 2.3.1 線形性 (Linearity)
$$ \mathcal{L}[\alpha f(t) + \beta g(t)] = \alpha F(s) + \beta G(s) $$
（積分の線形性より自明）

### 2.3.2 微分定理 (Differentiation Property)
制御工学で最も重要な性質です。時間領域での「微分」が、$s$領域での「$s$倍」に対応します。

> $$ \mathcal{L}\left[ \frac{df(t)}{dt} \right] = s F(s) - f(0) $$

**【証明】** 部分積分を用います。
$$
\begin{aligned}
\mathcal{L}[f'(t)] &= \int_{0}^{\infty} f'(t) e^{-st} dt \\
&= \left[ f(t) e^{-st} \right]_{0}^{\infty} - \int_{0}^{\infty} f(t) (-s) e^{-st} dt \\
&= (0 - f(0)) + s \int_{0}^{\infty} f(t) e^{-st} dt \\
&= s F(s) - f(0) \quad \blacksquare
\end{aligned}
$$
※ 2階微分の場合：
$$ \mathcal{L}[f''(t)] = s^2 F(s) - s f(0) - f'(0) $$

### 2.3.3 積分定理 (Integration Property)
> $$ \mathcal{L}\left[ \int_{0}^{t} f(\tau) d\tau \right] = \frac{1}{s} F(s) $$

**【証明】**
$g(t) = \int_{0}^{t} f(\tau) d\tau$ とおくと、$g'(t) = f(t)$ かつ $g(0)=0$ です。
微分定理より $\mathcal{L}[g'(t)] = s G(s) - g(0)$ なので、
$\mathcal{L}[f(t)] = s G(s) \implies G(s) = \frac{1}{s} F(s)$

### 2.3.4 推移定理（時間遅れ）
> $$ \mathcal{L}[f(t - L)] = e^{-Ls} F(s) \quad (L > 0) $$

### 2.3.5 最終値定理 (Final Value Theorem)
定常偏差の計算で使います。
> $$ \lim_{t \to \infty} f(t) = \lim_{s \to 0} s F(s) $$
※ ただし、極が左半平面にある（安定である）場合に限る。

**【証明の概略】**
微分定理 $\mathcal{L}[f'(t)] = sF(s) - f(0)$ において $s \to 0$ の極限をとります。
$$ \lim_{s \to 0} \int_{0}^{\infty} f'(t) e^{-st} dt = \lim_{s \to 0} (sF(s) - f(0)) $$
左辺は $e^0=1$ より $\int_{0}^{\infty} f'(t) dt = [f(t)]_{0}^{\infty} = f(\infty) - f(0)$ となります。
よって $f(\infty) - f(0) = \lim_{s \to 0} sF(s) - f(0)$ より成立。

---

## 2.4 逆ラプラス変換と部分分数展開

$F(s)$ から $f(t)$ を求めるには、通常は積分計算ではなく、**部分分数展開**をして基本的な変換対（表）を利用します。

### ケース1: 実数極のみ（単根）
例: $F(s) = \frac{s+3}{(s+1)(s+2)}$
$$ F(s) = \frac{A}{s+1} + \frac{B}{s+2} $$
とおきます。**ヘヴィサイドの展開定理（留数計算）** を用いると便利です。
$$ A = \left. (s+1)F(s) \right|_{s=-1} = \frac{-1+3}{-1+2} = 2 $$
$$ B = \left. (s+2)F(s) \right|_{s=-2} = \frac{-2+3}{-2+1} = -1 $$
よって $F(s) = \frac{2}{s+1} - \frac{1}{s+2}$。
逆変換すると、
$$ f(t) = 2e^{-t} - e^{-2t} $$

### ケース2: 重根がある場合
例: $F(s) = \frac{1}{s(s+1)^2}$
$$ F(s) = \frac{A}{s} + \frac{B_1}{s+1} + \frac{B_2}{(s+1)^2} $$
$A$ は通常通り求まります ($A=1$)。
$B_2$ も通常通り求まります ($B_2 = \left. (s+1)^2 F(s) \right|_{s=-1} = -1$)。
$B_1$ は微分を用いて求めます：
$$ B_1 = \frac{d}{ds} \left[ (s+1)^2 F(s) \right]_{s=-1} = \frac{d}{ds}\left[\frac{1}{s}\right]_{s=-1} = \left[-\frac{1}{s^2}\right]_{-1} = -1 $$
よって $F(s) = \frac{1}{s} - \frac{1}{s+1} - \frac{1}{(s+1)^2}$。
$\mathcal{L}^{-1}[\frac{1}{(s+a)^n}] = \frac{t^{n-1}}{(n-1)!}e^{-at}$ を利用して、
$$ f(t) = 1 - e^{-t} - t e^{-t} $$

### ケース3: 複素共役根がある場合（平方完成）
例: $F(s) = \frac{1}{s^2 + 2s + 5}$
分母の解が $s = -1 \pm j2$（複素数）になる場合です。平方完成を用います。
$$ s^2 + 2s + 5 = (s+1)^2 + 2^2 $$
これに合わせて分子を調整します。基本形は $\frac{\omega}{(s+a)^2 + \omega^2} \leftrightarrow e^{-at}\sin\omega t$ です。
$$ F(s) = \frac{1}{(s+1)^2 + 2^2} = \frac{1}{2} \cdot \frac{2}{(s+1)^2 + 2^2} $$
よって、
$$ f(t) = \frac{1}{2} e^{-t} \sin 2t $$
