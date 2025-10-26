---
title: 'Optimal Stochastic Nonsmooth Nonconvex Optimization'
toc: true
excerpt_separator: <!--more-->
tags:
  - 优化
---

Paper Reading: Optimal Stochastic Non-smooth Non-convex Optimization through Online-to-Non-convex Conversion

<!--more-->

文章提出如下基于online learning的算法，可以达到nonsmooth nonconvex optimization问题的最优复杂度。


![image-20230215153359260](/images/posts/NSNC/image-20230215153359260.png)


使用online 梯度下降算法作为子算法 $\mathcal{A}$, 并且每次向直径为 $D$ 的球内投影，可以得到如下的更新公式，

$$
\begin{align*}
x_n &= x_{n-1} +\Delta_n \\
g_n &= \nabla g(x_{n-1} + s_n \Delta_n; \xi_t) \\
\Delta_{n+1} &= \text{Clip}_D (\Delta_n + \eta g_n). 
\end{align*}
$$


其中 $s_n$ 从区间 $ [0,1] $ 之内均匀采样，首先观察到


$$
\begin{align*}
F(x_n) - F(x_{n-1}) = \int_{0}^1  \langle g_n, \Delta_n \rangle  {\rm d} s.
\end{align*}
$$



取期望并求和后可以得到，



$$
\begin{align*}
\mathbb{E} [ F(x_M)] = F(x_0) + \sum_{n=1}^M \mathbb{E} \langle g_n, \Delta_n - u_n \rangle + \sum_{n=1}^M \mathbb{E} \langle g_n,u_n \rangle. 
\end{align*}
$$



第一项即为online learning中regret的定义，而第二项选取合适的 $u_n$ 可以与Goldstein稳定点建立联系，将 $u_n$ 分成 $K$ 个阶段，每个阶段有 $T$ 份，并且每一段的取值如下给出，



$$
\begin{align*}
  u_k  = - D \frac{\sum_{t=1}^T \nabla F(w_t^k)}{\Vert \sum_{t=1}^T \nabla F(w_t^k) \Vert }.
\end{align*}
$$

据此可以定义如下的 $K$-偏移遗憾界:

$$
\begin{align*}
{\rm Reg}(u_1,\cdots,u_K) = \sum_{k=1}^K \sum_{n=(k-1)T+1}^{kT} \langle g_n, \Delta_n - u_k \rangle.
\end{align*}
$$

那么有，


$$
\begin{align*}
F^* &\le F(x_0) + \sum_{n=1}^M \mathbb{E} \langle g_n, \Delta_n - u_n \rangle + \sum_{n=1}^M \mathbb{E} \langle g_n,u_n \rangle \\
 &\le  F(x_0) + {\rm Reg}(u_1,\cdots,u_K) + \sum_{n=1}^M \mathbb{E} \langle g_n,u_n \rangle \\
 &= F(x_0) + {\rm Reg}(u_1,\cdots,u_K) + \sum_{k=1}^K  \mathbb{E} \left \langle \sum_{t=1}^T \nabla g(w_t^k;\xi_t^k),u_k \right\rangle \\
 &\le F(x_0) + {\rm Reg}(u_1,\cdots,u_K) + \sum_{k=1}^K  \mathbb{E} \left \langle \sum_{t=1}^T \nabla g(w_t^k),u_k \right\rangle + D G K \sqrt{T} \\
 &= F(x_0) + {\rm Reg}(u_1,\cdots,u_K) - D \sum_{k=1}^K  \mathbb{E} \left \Vert \sum_{t=1}^T \nabla g(w_t^k) \right\Vert + D G K \sqrt{T}.
\end{align*}
$$



其中 $G$ 为Lipschitz系数的上界，因此自然地为随机梯度方差的上界。



移项后并且代入online梯度下降算法的regret bound，可以得到 ${\rm Reg}(u_1,\cdots,u_K) = \mathcal{O}(K DG \sqrt{T})$, 那么


$$
\begin{align*}
\frac{1}{K} \sum_{k=1}^K \mathbb{E} \left \Vert \frac{1}{T} \sum_{t=1}^T \nabla F(w_t^k) \right \Vert \le \frac{F(x_0) - F^\ast}{D T K} + \frac{2G}{\sqrt{T}}.
\end{align*}
$$


令 $D T = \delta$ , 可以使得 $\Vert w_t^k - \bar w_t^k \Vert \le \delta$,  因此上面的算法可以在 $\mathcal{O}( G^2 \Delta\delta^{-1} \epsilon^{-3})$ 的时间内找到 $(\delta,\epsilon$) - 稳定点.



