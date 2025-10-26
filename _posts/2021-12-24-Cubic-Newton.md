---
title: 'Cubic Regularized Newton Method'
toc: true
excerpt_separator: <!--more-->
tags: 
  - 优化
  - 牛顿法
---





三次正则牛顿方法或者称为Cubic Newton方法在非凸函数和凸函数上的全局收敛性分析。



<!--more-->



以BFGS等为代表的 [拟牛顿法](https://truenobility303.github.io/BFGS/) 可以达到超线性收敛，但是不能给出具体的收敛阶，而三次正则化Newton方法通过最小化函数的一个三次近似上界，在Newton方法的二阶梯度上加上正则项来达到全局收敛性。



## Cubic Regularized Newton



三次正则化Newton方法的可以使得Newton方法达到全局收敛，其假设二阶梯度 Hessian矩阵Lipschitz连续，也即：


$$
\Vert \nabla^2 f(y) - \nabla^2 f(x) \Vert \le L \Vert y - x \Vert
$$



在上述假设下，利用带积分余项的Taylor展开可以得到，

$$
\begin{align*}
\Vert \nabla f(y) - \nabla f(x)  - \nabla^2 f(x) (y-x) \Vert &\le \frac{L}{2} \Vert y - x \Vert^2 \\
\vert f(y) - f(x) - \langle \nabla f(x) , y-x \rangle - \frac{1}{2} \langle \nabla^2 f(x) (y-x) , y-x \rangle &\le \frac{L}{6} \Vert y - x \Vert^3
\end{align*}
$$


下面依次证明两式，对于第一个式子有

$$
\begin{align*}
&\quad \Vert \nabla f(y) - \nabla f(x)  - \nabla^2 f(x) (y-x) \Vert  \\
&= \Vert \int_0^1\nabla^2 f(x+t(y-x))(y-x) dt  - \nabla^2 f(x) (y-x) \Vert \\
& \le \Vert \int_0^1 (\nabla^2 f(x+t(y-x) - \nabla^2f(x)) dt\Vert \Vert y-x \Vert \\
&\le L  \int_0^1 t dt\Vert y-x \Vert^2  \\
&= \frac{L}{2} \Vert y - x \Vert^2  \\
\end{align*}
$$

对于第二个式子，有

$$
\begin{align*}
&\quad \vert f(y) - f(x) - \langle \nabla f(x) , y-x \rangle - \frac{1}{2} \langle \nabla^2 f(x) (y-x) , y-x \rangle \\
&= \vert \int_0^1 \langle\nabla^2(x +t(y-x)) (1-t) (y-x)dt,y-x \rangle  - \frac{1}{2} \langle \nabla^2 f(x) (y-x) , y-x \rangle \\
&= \vert \int_0^1 \langle\nabla^2(x +t(y-x)) (1-t) (y-x)dt - \frac{1}{2}\nabla^2 f(x) (y-x),y-x \rangle \vert \\
&\le \vert \int_0^1 \nabla^2(x +t(y-x)) (1-t) dt - \frac{1}{2}\nabla^2 f(x) \vert \Vert y-x \Vert^2 \\
&\le L \int_0^1 t(1-t) dt \Vert y-x \Vert^3 \\
&=\frac{L}{6} \Vert y - x \Vert^3 
\end{align*}
$$


三次正则化Newton方法的核心就在于利用第二个式子，希望最小化三次正则后的函数，

$$
\begin{align*}
f(y) &\le f(x) + \langle \nabla f(x) , y-x \rangle +\frac{1}{2} \langle \nabla^2 f(x) (y-x) , y-x \rangle + \frac{M}{6} \Vert y - x \Vert^3 \\
\bar f_M(y) &= \min_y [f(x) + \langle \nabla f(x) , y-x \rangle +\frac{1}{2} \langle \nabla^2 f(x) (y-x) , y-x \rangle + \frac{M}{6} \Vert y - x \Vert^3] \\
T_M(x) &= \text{argmin}_y [ \langle \nabla f(x) , y-x \rangle +\frac{1}{2} \langle \nabla^2 f(x) (y-x) , y-x \rangle + \frac{M}{6} \Vert y - x \Vert^3 ] 
\end{align*}
$$



算法使用 $T_M(x)$ 生成迭代序列即可，但需要选择的是 $M$, 该部分会在收敛性证明的部分详细介绍，


$$
x_{k+1} = T_M(x_k)
$$


利用驻点的条件可以得到，

$$
\begin{align*}
\nabla f(x) + \nabla^2 f(x) (T_M(x) -x)  + \frac{M}{2} \Vert T_M(x) - x \Vert (T_M(x) - x) &= 0 \\
\nabla f(x) + \nabla^2 f(x) (T_M(x) -x)  + \frac{M}{2} r_M(x) (T_M(x) - x) &= 0,\text{Let } r_M(x) = \Vert T_M(x) - x \Vert \\
\langle \nabla f(x) , T_M(x) -x \rangle + \langle \nabla^2 f(x) (T_M(x) -x),T_M(x) -x \rangle+ \frac{M}{2} r_M(x)^3 &= 0
\end{align*}
$$



据此可以得到  $\Vert \nabla f (T_M(x)) \Vert $  的界，


$$
\begin{align*}
\Vert \nabla f(x) + \nabla^2 f(x) (T_M(x) -x) \Vert &= \frac{M}{2} r_M(x)^2 \\
\Vert \nabla f(T_M(x)) - \nabla f(x)  - \nabla^2 f(x) (T_M(x)-x) \Vert &\le \frac{L}{2} r_M(x)^2  \\
\Vert \nabla f(T_M(x)) \Vert &\le \frac{L+M}{2} r_M(x)^2
\end{align*}
$$



利用三次正则化的对偶问题，可以得到一个非常重要的结论，下面的证明并不完全严谨，因为下面的证明直接假设了强对偶性满足，实际上可以通过分类讨论并且验证Slater条件说明强队偶性成立，


$$
\begin{align*}
\min_h v(h) &= \min_h[g^T h+ \frac{1}{2} h^T H h+ \frac{M}{6} \Vert h \Vert^3 ], \text{With } g = \nabla f(x),h = T_M(x)-x ,H = \nabla^2 f(x)\\
&= \min_{h,\tau} [g^T h+ \frac{1}{2} h^T H h + \frac{M}{6} \tau^{\frac{3}{2}}], \text{ s.t. } \Vert h \Vert^2 \le \tau \\
&= \min_{h,\tau} \max_{\lambda} [g^T h + \frac{1}{2} h^T Hh + \frac{M}{6} \tau ^{\frac{3}{2}} + \lambda(\frac{1}{2} \Vert h \Vert^2 -  \frac{1}{2} \tau)] \\
&= \max_{\lambda} \min_{h,\tau} [g^T h + \frac{1}{2} h^T Hh + \frac{M}{6} \tau ^{\frac{3}{2}} + \lambda(\frac{1}{2} \Vert h \Vert^2 -  \frac{1}{2} \tau)] \\
&= \max_{\lambda} \min_{h} [g^T h + \frac{1}{2} h^T Hh + \frac{\lambda}{2} \Vert h \Vert^2 - \frac{2\lambda^3}{3M^2} ], \text{With } \tau = (\frac{2\lambda}{M})^2 \\ 
\end{align*}
$$


根据最优值的条件，


$$
\begin{align*}
\frac{\partial v(h)}{\partial h} &= g + Hh+ \lambda h =0 \\
\frac{\partial v(h)}{\partial \lambda} &= \frac{1}{2} \Vert h \Vert^2 - \frac{2\lambda^2}{M^2} =0 \\
\end{align*}
$$



经过分类讨论可以知道，此处较为简略地给出结论,



$$
\begin{align*}
\nabla^2 f(x)+ \lambda I &\succeq 0 \\
\nabla^2 f(x) + \frac{M}{2} r_M(x) I &\succeq 0 \\
\end{align*}
$$



从而得到更新公式的显式表达式，可以看到等价于Newton方法加上了正则项系数，


$$
\begin{align*}
x_{k+1} &= x_k - (\nabla^2f(x_k) + \lambda_k I)^{-1}\nabla f(x_k) \\
&=x_k - (\nabla^2f(x_k) + \frac{M}{2}r_M(x) I)^{-1}\nabla f(x_k) 
\end{align*}
$$


进一步可以得到，

$$
\begin{align*}
\nabla^2f(x) + \frac{M}{2} r_M(x) I &\succeq 0 \\
(T_M(x) -x) ^T \nabla^2 f(x) (T_M(x) -x) + \frac{M}{2} r_M(x)^3 &\ge 0 \\
\langle \nabla f(x), T_M(x) -x \rangle &\ge 0
\end{align*}
$$


另外可以证明的是，迭代后新的函数值，可以被三次正则化的函数控制，也即：



$$
f(T_M(x)) \le f(x)
$$



证明较为简单，利用到之前的结论即可



$$
\begin{align*}
\bar f_M(x) &= \min_y [f(x) + \langle \nabla f(x) , y-x \rangle +\frac{1}{2} \langle \nabla^2 f(x) (y-x) , y-x \rangle + \frac{M}{6} \Vert y - x \Vert^3] \\ 
&\le \min_y [f(y) + \frac{L+M}{6} \Vert y-x \Vert^3] \\
\bar f_M(x)-f(x) &= \langle \nabla f(x) , T_M(x)-x \rangle  +\frac{1}{2} \langle \nabla^2 f(x) (T_M(x)-x) , T_M(x)-x \rangle + \frac{M}{6} r_M(x)^3  \\
&= \frac{1}{2}\langle \nabla f(x) , T_M(x)-x \rangle - \frac{M}{12} r_M(x)^3 \\
&\le -\frac{M}{12} r_M(x)^3 \\
\end{align*}
$$



最终如果取得 $M \ge L$ 可以得到单步下降的结论，

$$
\begin{align*}
f(T_M(x)) \le \bar f_M(x) \le f(x) - \frac{M}{12}r_M(x)^3 \le f(x)
\end{align*}
$$


## Non-Convex Convergence



本节证明三次正则Newton的收敛性，首先定义下面的极小值条件为，对于极小值点，我们关心其一阶条件和二阶条件能否得到满足，我们希望一阶和二阶条件可以被步长 $r_M(x)$ 所控制住，



对于一阶条件，使用关于 $\Vert \nabla f (T_M(x)) \Vert $  的界，
$$
\begin{align*}
\Vert \nabla f (T_M(x)) \Vert &\le \frac{L+M}{2} r_M(x)^2 \\
\sqrt{\frac{2}{L+M}  \Vert \nabla f (T_M(x)) \Vert} &\le r_M(x) 
\end{align*}
$$




对于二阶条件，使用关于Hessian矩阵的结论以及Hessian矩阵Lipschitz连续的条件，


$$
\begin{align*}
\nabla^2 T_M(x) &\succeq \nabla^2 f(x) - L r_M(x) I \\
&\succeq -(\frac{M}{2}+L) r_M(x),\text{With }\nabla^2f(x) + \frac{M}{2} r_M(x) I \succeq 0 \\
\end{align*}
$$


根据上述的一阶条件和二阶条件，定义如下的准则量，其被步长所控制住，


$$
\begin{align*}
\mu_M(x) &= \max[\sqrt{\frac{2}{L+M} \Vert \nabla f(x) \Vert} , -\frac{2}{2L+M} \lambda_{\min}(\nabla^2 f(x))] \\
\mu_M(T_M(x)) &\le r_M(x)
\end{align*}
$$



据此可以借助 $\mu_L(x)$ 为桥梁得到算法的收敛性，基于选取合适的 $M \in [L,2L]$, 若函数的Lipschitz系数已知，简单选择 $M=L$ 即可，若该系数未知，则需要选用搜索方法得到三次正则化中的 $M$ 

$$
\begin{align*}
f(x_0) - f(x_{\star}) &= \sum_{i=0}^{k-1} f(x_i) - f(x_{i+1}) \\
&\ge \frac{L}{12}\sum_{i=0}^{k-1} r_M(x_i)^3 \\
&\ge \frac{Lk}{12} \min_k \mu_M(x_{k})^3 \\
\end{align*}
$$


据此证明了算法的收敛性，并且观察其梯度条件，可以发现其明显优于普通的梯度法的收敛性，该结论应当也是在意料之中的，因为三次正则Newton方法利用了二阶梯度条件，


$$
\begin{align*}
\min_k \mu_L(x_k) &\le (\frac{12(f(x_0) - f(x_{\star}))}{kL})^{\frac{1}{3}}  \\
\min_k \Vert \nabla f(x_k) \Vert &\le (\frac{6(L+M)(f(x_0) - f(x_{\star}))}{kL})^{\frac{2}{3}} 
\end{align*}
$$


在 [梯度法](https://truenobility303.github.io/CG/) 一文中，我们实际上已经可以得到了普通的梯度法在非凸函数上的收敛率为 $O(k^{-\frac{1}{2}})$,  而三次正则Newton方法的收敛率为 $O(k^{-\frac{2}{3}})$ ，足以说明三次正则化的Newton方法具有更优的收敛率。



## Convex Convergence

本节分析三次正则Newton方法在凸函数上面的收敛率，但值得注意的是，本节证明中仅仅用到了星-凸函数的性质，定义如下：


$$
\begin{align*}
f(\alpha x_{\star} +(1-\alpha)x_k ) \le \alpha f(x_{\star}) + (1-\alpha) f(x_k), \alpha \in (0,1)
\end{align*}
$$


该性质对于凸函数显然是成立的，下面首先利用算法的迭代步得到，证明的过程中需要假设 $\Vert x_k - x_{\star} \Vert < D$


$$
\begin{align*}
f(x_{k+1}) &= f(T_M(x_k)) \\
&\le \min_y [f(y) + \langle \nabla f(x_k) , y-x_k \rangle +\frac{1}{2} \langle \nabla^2 f(x_k) (y-x_k) , y-x_k \rangle + \frac{M}{6} \Vert y - x_k \Vert^3] \\
&\le \min_y [f(y) + \frac{M+L}{6} \Vert y - x_k \Vert^3] \\
&\le \min_y  [f(y) + \frac{M+L}{6} \Vert y - x_k \Vert^3]:y = \alpha f(x_{\star}) +(1-\alpha) f(x_k) \\
f(x_{k+1}) - f(x_{\star}) &\le \min_y  [f(y)-f(x_{\star}) + \frac{M+L}{6} \Vert y - x_k \Vert^3]:y = \alpha x_{\star} +(1-\alpha) x_k \\
&= \min_{\alpha} [f(y) - f(x_{\star}) + \frac{(M+L)\alpha^3}{6} \Vert x_k -x_{\star} \Vert^3 ]\\
&\le \min_{\alpha} [(1-\alpha)(f(x_k) - f(x_{\star}) + \frac{(M+L)\alpha^3}{6} \Vert x_k -x_{\star} \Vert^3 ] \\
&\le \min_{\alpha} [(1-\alpha)(f(x_k) - f(x_{\star})) + \frac{L\alpha^3}{2} \Vert x_k -x_{\star} \Vert^3 ] ,\text{By } M \le 2L\\
&\le \min_{\alpha} [(1-\alpha)(f(x_k) - f(x_{\star})) + \frac{L\alpha^3D^3}{2}  ],\text{By Assumption} \\
&=(1-\frac{2}{3}\alpha_k) (f(x_k) -f(x_{\star})), \text{Let } \alpha_k^2 = \frac{2(f(x_k) - f(x_{\star}))}{3LD^3} \\
\end{align*}
$$


利用 $\alpha_k$ 的界可以得到算法的收敛性，当 $\alpha_k \ge 1$ 的时候只需要取步长 $\alpha=1$, 此时简单代入可以得到结果，考虑 $\alpha_k \le 1$的情况。


$$
\begin{align*}
f(x_{k+1}) - f(x_{\star}) &\le f(x_k) - f(x_{\star}) - \frac{2}{3} \alpha_k (f(x_k) - f(x_{\star})) \\
\alpha_{k+1}^2 &\le \alpha_k^2 - \frac{2}{3} \alpha_k^3  \le \alpha_k^2\\
\frac{1}{\alpha_{k+1}} - \frac{1}{\alpha_k} &= \frac{\alpha_k - \alpha_{k+1}}{\alpha_{k+1}\alpha_k} =\frac{\alpha_k^2 - \alpha_{k+1}^2}{\alpha_{k+1}\alpha_k(\alpha_k + \alpha_{k+1})} \ge \frac{1}{3} \\
\frac{1}{\alpha_k} &\ge \frac{1}{\alpha_0}+\frac{k}{3} =1+ \frac{k}{3} ,\text{Let } \alpha_0 = 1\\
\frac{2(f(x_k) - f(x_{\star}))}{3LD^3}=\alpha_k^2  &\le \frac{1}{(1+\frac{k}{3})^2} \\
f(x_k) - f(x_{\star}) &\le \frac{3LD^3}{2(1+\frac{k}{3})^2} \\
\end{align*}
$$


从最终的式子可以看见算法的收敛率为 $O(\frac{1}{k^2})$ , 这也与FISTA等算法的收敛率相匹配，作为一个二阶优化算法，达到一阶优化算法的最优收敛率也并不奇怪。





