---
title: '常微分方程组的数值解法'
toc: true
excerpt_separator: <!--more-->
tags: 		
  - 数值算法
---



考虑用数值方法求解一阶非线性常微分方程的初值问题



<!--more-->



问题形式如下：

$$
\begin{align*}
y'(x) = f(x,y) ,\quad y(x_0) = y_0.
\end{align*}
$$



## Single-Step Method



考虑采用显式单步法求解，



$$
\begin{align*}
y_{n+1} = y_n + h \varphi(x_n,y_n;h).
\end{align*}
$$

定义其局部截断误差为


$$
\begin{align*}
T(x,y;h) = y(x+h ) - y(x) - h \varphi(x,y(x);h).
\end{align*}
$$



这等价于假设前面误差为 $0$ 也即 $y_n = y(x_n)$ 的时候 $x_n$ 到 $x_{n+1}$ 所产生的误差。我们称单步法为 $p$ 阶相容方法，也简称 $p$ 阶方法，若



$$
\begin{align*}
\vert T(x,y;h) \vert \le C h^{p+1}, \quad \exists C.
\end{align*}
$$



我们称单步法为相容的，若



$$
\begin{align*}
\lim_{h \rightarrow \infty} \frac{1}{h} T(x,y;h) = 0.
\end{align*}
$$



我们假设函数在所给的区域内是一致连续和一致有界的。由于



$$
\begin{align*}
&\quad \frac{1}{h} T(x,y;h) + \varphi(x,y;h)  - f(x,y) \\
&= \frac{1}{h}[ y(x+h) - y(x)] - y'(x) \\
&= \frac{1}{h} \int_0^h [ y '(x+h) - y'(x)] {\rm d}t \\
&\le \max_{0 \le t \le h} \vert y'(x+t) - y'(x) \vert \\
&= \max_{0 \le t \le h} \vert f(x+t, y(x+t)) - f(x,y(x)) \vert \rightarrow 0.
\end{align*}
$$



最后一步使用了 $f$ 的一致连续性，因此上面的推导说明单步法相容的充分必要条件是  

$$
\begin{align*}
\lim_{h \rightarrow 0} \varphi(x,y;h) = f(x,y)$. 
\end{align*}
$$


下面我们讨论算法的收敛性。假设 $\{ x_n \} $ 为区间 $[a,b]$ 上的等距剖分: $x_{n} = a + n h$. 

令误差函数 $e_n(h) = y(x_n) - y_n$ 以及整体误差为 $E(h) = \max_n \vert e_n(h) \vert$. 我们称单步法为 $p$ 阶收敛的，若




$$
\begin{align*}
E(h) \le C h^p, \quad \exists C \ge 0.
\end{align*}
$$



我们称单步法收敛，若



$$
\begin{align*}
\lim_{h \rightarrow 0 } E(h ) \rightarrow 0.
\end{align*}
$$




下面我们证明，若 $\varphi(x,y;h)$ 关于 $x,y,h$ 都为连续函数，并且关于 $y$ 满足Lipschitz条件，那么单步法的收敛性和相容性等价。

首先我们由相容性推出收敛性，首先我们需要如下引理，: 若



$$
\begin{align*}
\vert \xi_{j+1} \vert \le (1+A) \vert \xi_j \vert + B, \quad \forall j.  
\end{align*}
$$



那么



$$
\begin{align*}
\vert \xi_j \vert \le \vert \xi_0 \vert \exp(jA) + \frac{B}{A} ( \exp(jA) -1).
\end{align*}
$$



对应于 [1] 中的引理 3.1，证明使用归纳法即可，细节可以参考 [1]. 注意到



$$
\begin{align*}
&\quad e_{n+1} - e_n \\
&= [y(x_{n+1} ) - y(x_n)] - [ y_{n+1} - y_n ] \\
&= [y(x_{n+1} ) - y(x_n)]  - h \varphi(x_n,y_n;h) \\
&= T(x_n;y(x_n);h)  + h [\varphi(x_n,y(x_n);h) - \varphi(x_n,y_n;h)].
\end{align*}
$$



应用Lipschitz条件，可以得到



$$
\begin{align*}
\vert e_{n+1} - e_n \vert \le h L \vert e_n \vert + h \alpha(h) , \quad \text{ for some } \alpha(h) \rightarrow 0.
\end{align*}
$$



也即



$$
\begin{align*}
\vert e_{n+1} \vert  \le (1 + h L)\vert e_n \vert + h \alpha (h).
\end{align*}
$$



应用所给的引理，我们知道



$$
\begin{align*}
\vert e_n \vert \le \frac{\alpha(h)}{L} \left( \exp( L(x_n - x_1) )  - 1\right), \quad \forall n . 
\end{align*}
$$



这就说明了收敛性，并且可以注意算法 $p$ 阶相容可以推出于算法 $p$ 阶收敛。

收敛性需要考虑整体误差，但相容性仅需要考虑局部截断误差，上面的分析说明了为了分析算法的收敛性，只需要分析算法的相容性即可，这就将比较复杂的整体误差的计算归约到了比较简单的局部截断误差的计算上面，从而简化了问题。

另一方面，假定算法具有收敛性，下面我们推出算法的相容性。首先考虑下面的初值问题



$$
\begin{align*}
y'(x) =  \varphi(x,y;0), \quad y(x_0) = y_0.
\end{align*}
$$



给定的单步法对于该问题一定是相容的，因此算法一定收敛到该初值问题的解

由于我们假设算法是相容的，该算法也收敛到给定的初值问题 $y'(x) = f(x,y)$ 的解，这意味着



$$
\begin{align*}
\lim_{h \rightarrow \infty} \varphi(x,y;h) = \varphi(x,y;0) = f(x,y).
\end{align*}
$$



这正好就是相容性的等价条件。下面我们分析算法的稳定性和绝对稳定性。

我们称算法是稳定的，如果对应于初值为 $y_0, z_0$ 的单步法的解 $y_n, z_n$ 满足



$$
\begin{align*}
\vert y_n - z_n \vert \le C \vert y_0 - z_0 \vert, \quad \exists C \ge 0.
\end{align*}
$$



经过简单的分析可以证明，如果单步法的增量函数 $\varphi(x,y;h)$ 关于 $y$ 满足Lipschitz条件，那么该方法是稳定的。

证明同样参考自 [1], 我们知道



$$
\begin{align*}
y_{n+1} &= y_n + h \varphi(x_n,y_n;h) \\
z_{n+1} &= z_n +  h \varphi(x_n,z_n;h).
\end{align*}
$$



因此



$$
\begin{align*}
\vert y_{n} - z_{n} \vert \le (1+ h L )\vert y_{n-1} - z_{n-1} \vert \le (1+hL)^{n} \vert y_0  -z_0 \vert \le \exp(L(b-a)) \vert y_0 - z_0 \vert.
\end{align*}
$$



这就证明了方法的稳定性。而算法的绝对稳定性针对于试验方程 



$$
\begin{align*}
y'(x) = \lambda y(x), \quad {\rm Re}(\lambda) <0.
\end{align*}
$$



对于该试验方程，单步法可以写成如下的形式



$$
\begin{align*}
y_{n+1} = E(\lambda h) y_n.
\end{align*}
$$



如果满足 $\vert E(\lambda h) \vert <1$, 那么我们称该单步法为绝对稳定的。在复平面上，复变量 $\lambda h$ 满足  $\vert E(\lambda h) \vert <1$ 的区域称为绝对稳定性区域。该区域和实轴的交集，称为绝对稳定性区间。因此在分析绝对稳定性区间时，可以不妨设 $\lambda$ 为实数。



下面我们以几个常见的方法举例，对于Euler法， 


$$
\begin{align*}
y_{n+1} &= y_n + h f(x_n,y_n).
\end{align*}
$$




令 $\bar h = \lambda h$. 其绝对稳定性区域为


$$
\begin{align*}
R = \{ \bar h : \vert 1 + \bar h \vert < 1 \}.
\end{align*}
$$


这是复平面上面以 $(-1,0)$ 为中心，以 $1$ 为半径的圆的内部。考虑隐式Euler法，


$$
\begin{align*}
y_{n+1} = y_n + h f(x_{n+1},y_{n+1}). 
\end{align*}
$$


其绝对稳定性区域为


$$
\begin{align*}
R = \left \{ \bar h : \left \vert \frac{1}{1 - \bar h}  \right \vert < 1 \right \} = \{\bar h: \vert 1 - \bar h \vert >1 \}.
\end{align*}
$$


这是一个以 $(1,0)$ 为中心，以 $1$ 为半径的圆的外部，这包含了整个左半平面。考虑梯度公式


$$
\begin{align*}
y_{n+1} = y_n + \frac{h}{2} \left( f(x_n,y_n) + f(x_{n+1},y_{n+1}) \right).
\end{align*}
$$


其绝对稳定性区域为


$$
\begin{align*}
R = \left \{\bar h: \left \vert \frac{2 + \bar h}{2 - \bar h} \right \vert <1  \right\} = \{\bar h: \vert \bar h+2 \vert < \vert \bar h -2 \vert \}.
\end{align*}
$$


这是点 $(-2,0)$ 以及点 $(2,0)$ 的垂直平分线的左侧区域，正好是整个左半平面。

像隐式Euler法以及梯度方法这样绝对稳定性区域包含整个左半平面的方法被称为$A$-稳定的。



## Runge-Kutta Method



Runge-Kutta (RK) 方法的增量函数具有如下的形式



$$
\begin{align*}
\varphi(x_n,y_n; h ) &=  \sum_{i=1}^s b_i k_i \\
\text{where  }k_i &= f \left(x_n + c_i h, y_n + h \sum_{j=1}^s a_{ij} k_j \right),  \quad c_i  = \sum_{j=1}^s a_{ij}.
\end{align*}
$$


对于显式的 RK方法，有 $c_1 = 0$. 此时有



$$
\begin{align*}
k_1 &= f(x_n,y_n) \\
k_2&= f(x_n + c_2 h , y_n+ a_{21} k_1 h  ) \\
k_3 &= f(x_n +c_3 h , y_n+ a_{31} k_1 h + a_{32} k_2 h) \\
&\cdots 
\end{align*}
$$



RK法的目标是构造高阶的数值方法，为此我们计算 $y(x)$ 的高阶导数



$$
\begin{align*}
y'(x) &= f(x,y(x)) \\
y''(x) &= \frac{\partial f}{\partial x} (x,y(x)) + y'(x) \frac{\partial f}{\partial y} (x,y(x)) \\
&= \frac{\partial f}{\partial x} (x,y(x)) + f(x,y(x))\frac{\partial f}{\partial y} (x,y(x)) \\
y'''(x) &= \left(\frac{\partial^2 f}{\partial x^2} + f \frac{\partial^2 f}{\partial x \partial y} \right)  (x,y(x))   \\
&\quad+ \left( \frac{\partial f}{\partial x}  + f \frac{\partial f}{\partial y} \right) \frac{\partial f}{\partial y} (x,y(x)) \\
&\quad + f \left( \frac{\partial^2 f}{\partial x \partial y}  + f \frac{\partial^2 f}{\partial y^2} \right) (x,y(x)) \\
&=  \left(\frac{\partial^2 f}{\partial x^2} + 2 f \frac{\partial^2 f}{\partial x \partial y}  +\frac{\partial f}{\partial x}\frac{\partial f}{\partial y} + f \left(  \frac{\partial f}{\partial y}\right)^2 + f^2 \frac{\partial^2 f}{\partial y^2} \right)  (x,y(x))
\end{align*}
$$



容易验证，1阶RK方法就是显式Euler算法 $\varphi(x,y;h) = f(x,y)$.

下面我们考虑2阶RK方法，简称 RK2，对 $k_2$ 进行Taylor展开，得到



$$
\begin{align*}
k_2 = f(x_n,y_n) +c_2 h \frac{\partial f}{\partial x} (x_n,y_n)  + a_{21}  h f(x_n,y_n) \frac{\partial f}{\partial y} (x_n,y_n) + o(h).
\end{align*}
$$

因此



$$
\begin{align*}
\varphi(x,y;h) = (b_1 + b_2) f(x,y) + b_2 c_2 h \frac{\partial f}{\partial x} (x,y) + b_2 a_{21 } h \left(f \frac{\partial f}{\partial y} \right)(x,y).
\end{align*}
$$



为了使其为2阶算法，需要满足条件



$$
\begin{align*}
b_1 + b_2 =1, \quad b_2 c_2 = \frac{1}{2}, \quad b_2 a_{21} = \frac{1}{2}.
\end{align*}
$$



类似地，我们考虑 RK3，对于 $k_2$ 我们在之前的基础上多展开一项，得到



$$
\begin{align*}
k_2 &= f(x_n,y_n) +c_2 h \frac{\partial f}{\partial x} (x_n,y_n)  + a_{21}  h f(x_n,y_n) \frac{\partial f}{\partial y} (x_n,y_n) \\
&\quad + \frac{(c_2h)^2 }{2} \frac{\partial^2 f}{\partial x^2} (x_n,y_n) + (c_2 h )(a_{21} h) f(x_n,y_n)\frac{\partial^2 f}{\partial x\partial y} (x_n,y_n)  \\
&\quad+ \frac{(a_{21} h)^2}{2} \left( [f(x_n,y_n)]^2 \frac{\partial^2 f}{\partial y^2} (x_n,y_n) \right)   + o(h^2).
\end{align*}
$$



对 $k_3$ 进行Taylor展开，



$$
\begin{align*}
k_3 &= f(x_n,y_n) + c_3 h \frac{\partial f}{\partial x} (x_n,y_n) +  (a_{31} k_1 h + a_{32} k_2 h) \frac{\partial f}{\partial y} (x_n,y_n) \\
&\quad + \frac{1}{2} (c_3h)^2 \frac{\partial^2 f}{\partial x^2} (x_n,y_n)  + (c_3h ) (a_{31} k_1  h + a_{32} k_2 h) \frac{\partial^2 f }{\partial x \partial y}   (x_n,y_n) \\
&\quad + \frac{1}{2} (a_{31} k_1 h + a_{32} k_2 h)^2 \frac{\partial^2}{\partial y^2} (x_n,y_n) \\
&= f(x_n,y_n) + c_3  h \frac{\partial f}{\partial x} (x_n,y_n) \\
&\quad +\left(a_{31} f(x_n,y_n) h  + a_{32} f(x_n,y_n) h + a_{32} c_2 h^2 \frac{\partial f}{\partial x} (x_n,y_n)  + a_{32} a_{21}  h^2 f(x_n,y_n) \frac{\partial f}{\partial y} (x_n,y_n)\right) \frac{\partial f}{\partial y} (x_n,y_n) \\
&\quad + \frac{c_3^2 h^2}{2} \frac{\partial^2 f}{\partial x^2} (x_n,y_n) +c_3 (a_{31} + a_{32}) h^2 f(x_n,y_n )\frac{\partial^2 f }{\partial x \partial y}   (x_n,y_n) \\
&\quad + \frac{1}{2} (a_{31} + a_{32})^2  h^2 [f(x_n,y_n)]^2 \frac{\partial^2 f }{\partial y^2}   (x_n,y_n) + o(h^2).
\end{align*}
$$



额外使用行和条件



$$
\begin{align*}
c_2 = a_{21}, \quad c_3 = a_{31} + a_{32}. 
\end{align*}
$$



可以得到RK3 方法对应的增量函数为



$$
\begin{align*}
\varphi(x,y;h) &= (b_1 + b_2 +b_3) f(x,y) +  (b_2 c_2 + b_3c_3) h \frac{\partial f}{\partial x} (x,y) \\
&\quad + (b_2 c_2 + b_3 c_3) h \left(f \frac{\partial f}{\partial y} \right)(x,y) \\
&\quad + \frac{1}{2} (b_2 c_2^2 + b_3 c_3^2) h^2 \frac{\partial^2 f}{\partial x^2} (x,y)\\
&\quad + (b_2 c_2^2 + b_3 c_3^2) h^2 \left( f \frac{\partial^2 f}{\partial x \partial y} \right) (x,y) + b_3 a_{32} c_2 h^2 \left( \frac{\partial f}{\partial x} \frac{\partial f}{\partial y} \right)(x,y) \\
&\quad + b_3 a_{32} c_2 h^2 \left( f  \frac{\partial^2 f}{\partial x \partial y}  \right)(x,y) + \frac{1}{2} (b_2 c_2^2 +b_3 c_3^2) h^2 \left( f^2 \frac{\partial^2 f}{\partial y^2}\right)(x,y).
\end{align*}
$$



为了使其为3阶算法，需要满足条件



$$
\begin{align*}
\begin{cases}
b_1 + b_2 + b_3 = 1 \\
b_2 c_2 + b_3 c_3 = \dfrac{1}{2} \\
b_2 c_2^2 + b_3 c_3^2 = \dfrac{1}{3} \\
b_3 a_{32} c_2 = \dfrac{1}{6}.
\end{cases}
\end{align*}
$$



## Linear Multi-Step Method

本节考虑线性多步法，其具有如下的一般形式



$$
\begin{align*}
\sum_{j=0}^k \alpha_j y_{n+ j} = h \sum_{j=0}^k \beta_j f(x_{n+j}, y_{n+j}).
\end{align*}
$$


当 $k=1$ 时退化为之前所介绍过的单步法，类似地，定义局部截断误差



$$
\begin{align*}
T_{n+k}  = \sum_{j=0}^k \alpha_j y(x_{n+j}) - h \sum_{j=0}^k \beta_j f(x_{n+j}, y(x_{n+j})).
\end{align*}
$$



类似于单步法，我们可以通过Taylor展开确定局部截断误差主项，进而确定方法的阶数。我们称方法相容，若



$$
\begin{align*}
\lim_{h \rightarrow 0} \frac{1}{h}T_{n+k} = 0.
\end{align*}
$$



为了刻画相容性，我们引入如下的第一和第二特征多项式



$$
\begin{align*}
\rho(\lambda) = \sum_{j=0}^{k} \alpha_j \lambda^j, \quad \sigma(\lambda) = \sum_{j=0}^k \beta_j \lambda^j.
\end{align*}
$$

注意到



$$
\begin{align*}
T_{n+k} &= \sum_{j=0}^k \alpha_j y(x_{n+ j}) - h \sum_{j=0}^k \beta_j y'(x_{n+j}) \\
&= \sum_{j=0}^{k} \alpha_j \left[ y(x_n) + jh y'(x_n) + o(h^2) \right]  -h \sum_{j=0}^k \beta_j y'(x_n) + o(h^2) \\
&= y(x_n) \sum_{j=0}^k \alpha_j  + y'(x_n) h  \sum_{j=0}^k (j \alpha_j - \beta_j) +o(h^2).
\end{align*}
$$



因此多步长法相容当且仅当 



$$
\begin{align*}
\sum_{j=0}^k \alpha_j = 0, \quad \sum_{j=0}^k (j \alpha_j  - \beta_j)  = 0
\end{align*}
$$



这就等价于条件 $\rho(1) = 0$ ,$\rho'(1) = \sigma(1)$.





类似于，单步法的稳定性，我们可以定义如下的 (零) 稳定性。我们称方法为 (零) 稳定的， 若任意两解满足



$$
\begin{align*}
\vert u_n - v_n \vert \le C \vert u_0 - v_0 \vert, \quad \exists C \ge 0.
\end{align*}
$$



我们下面不加证明地给出如下结论，其证明较为繁琐，感兴趣的读者可以参考 [2].

我们称多步法满足根条件，若其第一特征多项式的零点都在单位圆内或在单位圆周上，且在单位圆周上的根为单根。

而线性多步法稳定的充分必要条件是其满足根条件。

进一步，若线性多步法相容且稳定，那么算法收敛。



记 $Y_n = [y_n; y_{n+1};\cdots; y_{n+k-1}]$. 由于 $\alpha_k=1$, 多步法的矩阵等价形式为



$$
\begin{align*}
Y_{n+1} = 
\begin{pmatrix}
0 & 1 & 0 & \cdots \\
0 & 0 & 1 & \cdots \\
\cdots & \cdots & \cdots & \cdots \\
-\alpha_0 & - \alpha_1 &\cdots & -\alpha_{k-1}
\end{pmatrix} + 
h 
\begin{pmatrix}
0 \\
0 \\
\vdots \\
\varphi
\end{pmatrix}.
\end{align*}
$$



经过简单的推导可以证明，上述的系数矩阵的特征多项式就是第一特征多项式。

因此检验多步法的根条件时，也可以仅仅通过上面的系数矩阵的谱来判断。



对于多步法的绝对稳定性，可以考虑如下的稳定性多项式


$$
\begin{align*}
\pi(\lambda, \bar h) = \rho(\lambda) - \bar h\sigma(\lambda).
\end{align*}
$$


若对于给定的 $\bar h$ 满足 ${\rm Re}(\bar h)\le 0$, 若上述的稳定性多项式对应的根小于1，称方法关于此 $\bar h$ 为绝对稳定的。

同样的，这样可以给出绝对稳定性区域和区间的概念。其中，绘制绝对稳定性区域时，可以用单位根描绘其边界，即考虑轨迹


$$
\begin{align*}
\bar h = \frac{\rho (\exp(i \theta))}{\sigma (\exp(i \theta))}.
\end{align*}
$$


这称为边界轨迹技巧。



## BDF2 Method



本节介绍Backward Differentiation Formula (BDF) 算法类中的BDF2方法，并且介绍其性质。

首先，我们从Newton插值多项式出发推导对应的公式，我们考虑变步长版本。

考虑三点 $(0,y_n)$, $(\alpha h,y_{n+1})$, $((1+\alpha)h, y_{n+2})$ 构成的Newton插值多项式，也即


$$
\begin{align*}
P_2(t) = y_n + \frac{y_{n+1} - y_n}{\alpha h} t + \dfrac{\dfrac{y_{n+2} - y_{n+1}}{h}- \dfrac{y_{n+1} - y_n}{\alpha h}}{(1+\alpha )h} t ( t - \alpha h).
\end{align*}
$$


其导数为


$$
\begin{align*}
P_2'(t) = \frac{y_{n+1} - y_n}{\alpha h} + \dfrac{\dfrac{y_{n+2} - y_{n+1}}{h}- \dfrac{y_{n+1} - y_n}{\alpha h}}{(1+\alpha )h}  (2t  -\alpha h).
\end{align*}
$$


我们希望 $P_2'((1+\alpha )h) \approx y_{n+2}' = f(x_{n+2},y_{n+2})$, 这就导出了BDF2的迭代公式，形如


$$
\begin{align*}
\frac{2 +\alpha}{1+ \alpha} y_{n+2} - \frac{1+\alpha}{\alpha} y_{n+1} + \frac{1}{\alpha (1+\alpha)} = h f(x_{n+2},y_{n+2}).
\end{align*}
$$


将系数归一化得到线性多步法的标准形式


$$
\begin{align*}
y_{n+2} - \frac{(1+\alpha)^2}{\alpha (2+\alpha )} y_{n+1} + \frac{1}{\alpha (2+ \alpha)} y_n = \frac{1+\alpha}{2+\alpha} h f(x_{n+2},y_{n+2}).
\end{align*}
$$




该方法的系数矩阵为


$$
\begin{align*}
A = 
\begin{pmatrix}
0 & 1 \\
\dfrac{1}{\alpha (2+\alpha)}  &\dfrac{(1+\alpha)^2}{\alpha (2+\alpha)} & 
\end{pmatrix}
= 
\begin{pmatrix}
0 & 1 \\
-\beta & \beta +1 
\end{pmatrix}
,\quad \beta = \frac{1}{\alpha (2+\alpha)}.
\end{align*}
$$


下面我们介绍强稳定性的概念，并且对于BDF2方法进行分析。

我们称方法为强稳定的，若存在 $h_0$ 使得对任意的 $0 \le h \le h_0$ ，都存在一个范数使得上述系数矩阵的范数不大于1.

我们不加证明地给出如下两个结论：

首先，强稳定性是零稳定性的充分条件。

其次，考虑变步长方法，若对应的等步长方法的第一特征多项式满足根条件，并且 $1$ 为单位圆周上的单根，那么存在一个步长的比例系数，使得对应的变步长方法是强稳定的。

我们不给出具体的证明，但以BDF2方法为例进行解释。

计算得到对应的等步长方法 $(\alpha=1)$ 的第一特征多项式的根为 $1$ 以及 $1/3$, 满足根条件并且 $1$ 为单位圆周上的单根。

我们下面构造一个矩阵范数，证明变步长的BDF2方法为强稳定的。

考虑 $\alpha=1$ 时系数矩阵对应的特征向量，这些特征向量所构成的矩阵为


$$
\begin{align*}
P = 
\begin{pmatrix}
1 & 3 \\
1 & 1
\end{pmatrix}
\end{align*}
$$


对于这个矩阵，


$$
\begin{align*}
 P^{-1} A P 
= -\frac{1}{2} \begin{pmatrix}
1 & -3 \\
-1 & 1 
\end{pmatrix}
\begin{pmatrix}
0 & 1 \\
-\beta & \beta +1 
\end{pmatrix}
\begin{pmatrix}
1 & 3 \\
1 & 1
\end{pmatrix} 
= 
\begin{pmatrix}
1 & 1 - 3 \beta\\
0 & \beta
\end{pmatrix}
\end{align*}
$$




我们希望对该矩阵取列和范数 $\Vert \cdot \Vert_1$, 使得矩阵的该范数不大于1. 可以得到条件 $0 \le \beta \le 1/2$.

进一步还可以获得更广的步长范围，考虑矩阵范数


$$
\begin{align*}
\Vert A \Vert = \Vert D^{-1} P^{-1} AP D \Vert_1, \quad D = 
\begin{pmatrix}
1 & 0 \\
0 & \epsilon
\end{pmatrix}
\end{align*}
$$


计算得到


$$
\begin{align*}
D^{-1} P^{-1} AP D = 
\begin{pmatrix}
1 & (1 - 3 \beta) \epsilon \\
0 & \beta
\end{pmatrix}.
\end{align*}
$$


此时选取一个更小的 $\epsilon$ 可以得到 $\beta$ 更广的范围，从而反过来推出 $\alpha$ 的范围。



下面我们接着以 BDF2方法为例，分析其其他形式，为了简单起见，我们仅仅考虑等步长版本


$$
\begin{align*}
y_{n+2} - \frac{4}{3} y_{n+1} + \frac{1}{3} y_n = \frac{2}{3} h  f(x_{n+2}, y_{n+2}).
\end{align*}
$$




回忆


$$
\begin{align*}
y'(x) &= f(x,y(x)) \\
y''(x) &=  \left(\frac{\partial f}{\partial x}+ f \frac{\partial f}{\partial y} \right)(x,y(x)).
\end{align*}
$$


我们计算方法的截断误差如


$$
\begin{align*}
T_{n+2} & = y(x_{n+2}) - \frac{4}{3} y(x_{n+1}) + \frac{1}{3} y(x_n) - \frac{2}{3} h  f(x_{n+2}, y(x_{n+2})) \\
&= y(x_n) + 2h f(x_n,y_n)+ \frac{(2h)^2}{2} \left(\frac{\partial f}{\partial x}+ f \frac{\partial f}{\partial y} \right) (x_n,y_n) \\
&\quad - \frac{4}{3} \left( y(x_n) + h f(x_n,y_n) + \frac{h^2}{2} \left(\frac{\partial f}{\partial x}+ f \frac{\partial f}{\partial y} \right) (x_n,y_n) \right) + \frac{1}{3} y(x_n) \\
&\quad - \frac{2}{3} h \left( f(x_n,y_n) + 2h \left( \frac{\partial f}{\partial x} + f\frac{\partial f}{\partial y} \right) (x_n,y_n)   \right) + o(h^2) \\
&= o(h^2).
\end{align*}
$$


因此方法为2阶算法。这说明算法一定是相容的，如果仅要去判断其相容性，可以简单查看其第一和第二特征多项式，


$$
\begin{align*}
\rho(\lambda) = \lambda^2 - \frac{4}{3} \lambda + \frac{1}{3}, \quad \sigma(\lambda) = \frac{2}{3}.
\end{align*}
$$


发现其满足条件


$$
\begin{align*}
\rho(1) = 0, \quad \rho'(1) = \sigma(1).
\end{align*}
$$


说明方法是相容的。由于其第一特征多项式的根为 $1$ 和 $1/3$，满足根条件，方法也是稳定的。

根据之前的分析，如下的轨迹将给出其绝对稳定性区域的边界


$$
\begin{align*}
\bar h = \dfrac{3 \exp(2 i \theta) - 4\exp(i \theta) +1  }{2}.
\end{align*}
$$


而图像与实轴的交集就是方法的绝对稳定性区间。



## Reference

[1] 关冶、陆金甫，数值分析基础 (第三版)

[2] Quarteroni, Alfio, Riccardo Sacco, and Fausto Saleri. Numerical mathematics. Vol. 37. Springer Science & Business Media, 2006.
