# 第6章：状態空間表現 (詳細版)

本章では、古典制御（伝達関数）から現代制御（状態空間）への架け橋となる「状態空間表現」について学びます。
120分の講義を想定し、高階微分方程式からの変換手順、伝達関数との相互変換、座標変換（対角化）によるモード分解まで詳しく解説します。

---

## 6.1 状態空間表現の導出

### 6.1.1 状態変数とは？
システムのエネルギー状態や過去の履歴を一意に決定するために必要な最小限の変数の組を **状態変数 (State Variables)** と呼び、ベクトル $\boldsymbol{x}(t)$ で表します。
例えば、ばね・マス系では「位置」と「速度」があれば未来の挙動が決定できるため、これらが状態変数になります。

### 6.1.2 運動方程式からの導出
質量 $m$、減衰係数 $c$、ばね定数 $k$ の系に外力 $u(t)$ が加わるときの運動方程式：
$$ m \ddot{y} + c \dot{y} + k y = u $$
ここで、状態変数を以下のように定義します。
$$ x_1 = y \quad (\text{位置}) $$
$$ x_2 = \dot{y} \quad (\text{速度}) $$

それぞれの時間微分を考えます。
1.  $\dot{x}_1 = \dot{y} = x_2$
2.  $\dot{x}_2 = \ddot{y} = \frac{1}{m}(u - c\dot{y} - ky) = -\frac{k}{m}x_1 - \frac{c}{m}x_2 + \frac{1}{m}u$

これを行列形式にまとめると、**状態方程式 (State Equation)** になります。
$$
\begin{bmatrix} \dot{x}_1 \\ \dot{x}_2 \end{bmatrix} =
\begin{bmatrix} 0 & 1 \\ -\frac{k}{m} & -\frac{c}{m} \end{bmatrix}
\begin{bmatrix} x_1 \\ x_2 \end{bmatrix} +
\begin{bmatrix} 0 \\ \frac{1}{m} \end{bmatrix} u
$$
出力 $y$ は位置なので、**出力方程式 (Output Equation)** は、
$$ y = \begin{bmatrix} 1 & 0 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} $$

一般形：
> $$ \dot{\boldsymbol{x}} = A \boldsymbol{x} + B u $$
> $$ y = C \boldsymbol{x} + D u $$

ここで、$A$: システム行列、$B$: 入力行列、$C$: 出力行列、$D$: 直達行列 です。

---

## 6.2 伝達関数との関係

状態空間表現から伝達関数を一意に求めることができます（逆は一意ではありません）。

### 6.2.1 導出
状態方程式をラプラス変換します（初期値 $\boldsymbol{x}(0)=\boldsymbol{0}$）。
$$ s \boldsymbol{X}(s) = A \boldsymbol{X}(s) + B U(s) $$
$$ (sI - A) \boldsymbol{X}(s) = B U(s) $$
$$ \boldsymbol{X}(s) = (sI - A)^{-1} B U(s) $$

これを出力方程式 $Y(s) = C \boldsymbol{X}(s) + D U(s)$ に代入します。
$$ Y(s) = C (sI - A)^{-1} B U(s) + D U(s) = \{ C (sI - A)^{-1} B + D \} U(s) $$

よって、伝達関数 $G(s)$ は行列演算で求まります。
> $$ G(s) = C (sI - A)^{-1} B + D $$

ここで **$(sI - A)^{-1}$** をリゾルベント行列と呼びます。

### 6.2.2 具体例での計算
先ほどのばね・マス系（$m=1, c=2, k=1$ とする）で計算してみます。
$$ A = \begin{bmatrix} 0 & 1 \\ -1 & -2 \end{bmatrix}, \quad B = \begin{bmatrix} 0 \\ 1 \end{bmatrix}, \quad C = \begin{bmatrix} 1 & 0 \end{bmatrix}, \quad D = 0 $$

1.  $sI - A$ を計算：
    $$ sI - A = \begin{bmatrix} s & 0 \\ 0 & s \end{bmatrix} - \begin{bmatrix} 0 & 1 \\ -1 & -2 \end{bmatrix} = \begin{bmatrix} s & -1 \\ 1 & s+2 \end{bmatrix} $$

2.  逆行列 $(sI - A)^{-1}$ を計算：
    行列式 $\det(sI - A) = s(s+2) - (-1)(1) = s^2 + 2s + 1$
    $$ (sI - A)^{-1} = \frac{1}{s^2 + 2s + 1} \begin{bmatrix} s+2 & 1 \\ -1 & s \end{bmatrix} $$
    ※ $\det(sI-A)=0$ は特性方程式そのものです。つまり $A$ の固有値はシステムの極と一致します。

3.  $G(s)$ を計算：
    $$
    \begin{aligned}
    G(s) &= \begin{bmatrix} 1 & 0 \end{bmatrix} \frac{1}{(s+1)^2} \begin{bmatrix} s+2 & 1 \\ -1 & s \end{bmatrix} \begin{bmatrix} 0 \\ 1 \end{bmatrix} \\
    &= \frac{1}{(s+1)^2} \begin{bmatrix} 1 & 0 \end{bmatrix} \begin{bmatrix} 1 \\ s \end{bmatrix} \\
    &= \frac{1}{(s+1)^2} \cdot 1 = \frac{1}{s^2 + 2s + 1}
    \end{aligned}
    $$
    これは運動方程式 $\ddot{y} + 2\dot{y} + y = u$ を直接ラプラス変換した結果と一致します。

---

## 6.3 座標変換と対角化

状態変数の選び方は無数にあります。正則行列 $T$ を用いて $\boldsymbol{x} = T \boldsymbol{z}$ と変数変換すると、システムの見通しが良くなることがあります。

### 6.3.1 変換後のシステム行列
$$ \dot{\boldsymbol{x}} = A\boldsymbol{x} + Bu \implies T\dot{\boldsymbol{z}} = AT\boldsymbol{z} + Bu \implies \dot{\boldsymbol{z}} = (T^{-1}AT)\boldsymbol{z} + (T^{-1}B)u $$
新しいシステム行列は $\tilde{A} = T^{-1}AT$ となり、これは相似変換です（固有値は保存されます）。

### 6.3.2 対角化（モード分解）
$A$ が対角化可能であれば、$T$ として $A$ の固有ベクトルを並べた行列を選ぶことで、$\tilde{A}$ を対角行列 $\Lambda$ にできます。
$$ \Lambda = \begin{bmatrix} \lambda_1 & 0 \\ 0 & \lambda_2 \end{bmatrix} $$
システムは以下のように分解（非干渉化）されます。
$$ \dot{z}_1 = \lambda_1 z_1 + \tilde{b}_1 u $$
$$ \dot{z}_2 = \lambda_2 z_2 + \tilde{b}_2 u $$
それぞれの状態 $z_i$ は独立して動き、固有値 $\lambda_i$（極）に応じた挙動（モード）を示します。これを **モード分解** と呼びます。
制御系設計において、どのモードが支配的か、あるいはどのモードが可制御・可観測かを判断するのに役立ちます。

---

## 6.4 演習問題

**問題**:
システム行列 $A = \begin{bmatrix} 0 & 1 \\ -2 & -3 \end{bmatrix}$ の固有値を求め、安定性を判別せよ。

**解答**:
特性方程式 $|sI - A| = 0$ を解く。
$$ \det \begin{bmatrix} s & -1 \\ 2 & s+3 \end{bmatrix} = s(s+3) - (-2) = s^2 + 3s + 2 = (s+1)(s+2) = 0 $$
よって固有値は $\lambda_1 = -1, \lambda_2 = -2$。
すべての固有値の実部が負（$-1 < 0, -2 < 0$）であるため、このシステムは**安定**である。
