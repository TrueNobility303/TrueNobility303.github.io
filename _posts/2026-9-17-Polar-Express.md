---
title: 'Muon优化器中最优的极表达'
toc: true
excerpt_separator: <!--more-->
tags:
  - 优化
  - 非凸优化
  - 神经网络优化
---

Paper Reading: The Polar Express: Optimal Matrix Sign Methods and Their Application to the Muon Algorithm. [ICLR 2026 杰出论文奖提名]

<!--more-->

在 [神经网络中的优化器](https://truenobility303.github.io/NN-Optimizer/), 我们介绍了Muon算法，其已经逐渐成为目前LLM优化的主流算法。在Muon算法的流程中，核心的步骤是对矩阵参数的动量进行正交化操作。令 $M$ 为矩阵动量，设 $M = U \Sigma V^\top$ 为其奇异值分解，矩阵正交化操作返回 $U V^\top$.


Muon采用的Newton-Schulz迭代，本质是选择一个多项式函数 $p$, 使得 $p(x) \approx 1$, 那么将这个函数运用于矩阵 $M$, 就可以得到 $p(M) = U p(\Sigma) V^\top \approx U V^\top$. 事实上，该问题等价于数值算法中常见的 [最佳逼近问题](https://truenobility303.github.io/Fitting/)，下面我们转化为经典问题求解。

在Muon算法的实际应用中，为了每一步迭代的高效性，希望每步迭代的函数是一个低次多项式。具体地，我们假设 $p = f_T \circ \cdots \circ f_2 \circ f_1$, 其中每个 $f_i$ 为一个奇多项式，例如三次或者五次多项式。 我们希望在上述复合函数空间上求解如下的最佳一致逼近问题：

$$
\begin{align*}
\min_{f} \max_{x \in [l,u]} | f(x) - 1|.
\end{align*}
$$

在实际中，区间 $[l,u]$ 可以简单选取为 $[\epsilon, 1]$, 其中 $\epsilon$ 是一个大于0的小量，例如0.001. 文章的核心贡献在于说明，当我们只想逼近常值函数1时，上述复合函数的最佳一致逼近问题可以用如下的简单贪婪算法求解：迭代地计算

$$
\begin{align*}
f_t^\ast&\in \arg\min_{f_t} \max_{x\in[l_t,u_t]}|f_t(x)-1|,\\
l_{t+1} &=\min_{x\in[l_t,u_t]}f_t^\ast(x),\\
u_{t+1} &=\max_{x\in[l_t,u_t]}f_t^\ast(x).
\end{align*}
$$

因此每一步化简为简单的 [最佳逼近问题](https://truenobility303.github.io/Fitting/)， 可以用我们之前介绍过的基于Chebyshev定理的经典Remez算法求解。因此，我们只需要证明上述的贪婪算法也可以找到全局最优解，我们下面的推导参考了 [苏剑林的博客](https://spaces.ac.cn/archives/10996):

## 贪婪算法的全局最优性

先记单步最优误差为

$$
E(a,b)=\min_{q\in\mathcal P_d^{\mathrm{odd}}}
\max_{x\in[a,b]}|q(x)-1|,\qquad 0<a\le b.
$$

我们需要用到缩放的不变性。由于 $q(x)\mapsto q(sx)$ 是奇多项式空间到自身的双射，对任意 $s>0$，

$$
E(sa,sb)=E(a,b).
$$

此外，最优多项式的值域一定以 $1$ 为中心。事实上，线性候选 $q(x)=2x/(a+b)$ 给出

$$
E(a,b)\le\frac{b-a}{b+a}<1.
$$

故最优多项式的值域 $[A,B]$ 满足 $A>0$。对固定的 $0<A\le B$，直接比较两个端点可得

$$
\min_{c>0}\max\{|cA-1|,|cB-1|\}
=\frac{B-A}{B+A},\qquad c^\ast=\frac{2}{A+B}.
$$

如果 $A+B\ne2$，将多项式乘以 $c^\ast$ 就能严格减小误差，与最优性矛盾。因此，令 $e_t=E(l_t,u_t)$，有

$$
[l_{t+1},u_{t+1}]=[1-e_t,1+e_t].
$$

多项式连续，区间的像仍是整个区间，所以这也是前 $t$ 步贪婪复合的准确值域，其最大误差正好为 $e_t<1$。

现在对步数 $T$ 作归纳。$T=1$ 时，结论由定义成立。假设 $T-1$ 步贪婪解已全局最优，考虑任意候选

$$
q\circ h,\qquad
h=\widetilde f_{T-1}\circ\cdots\circ\widetilde f_1,
\qquad q,\widetilde f_t\in\mathcal P_d^{\mathrm{odd}}.
$$

设 $h([l,u])=[a,b]$。先处理符号：若 $0\in[a,b]$，则存在 $x_0$ 使 $h(x_0)=0$，而 $q(0)=0$，故最终误差至少为 $1$，不可能优于贪婪解。若 $b<0$，将 $h$ 换为 $-h$，同时将 $q(y)$ 换为 $q(-y)$，复合函数不变，且各步仍为同次数上界的奇多项式。因此只需考虑 $0<a\le b$。 令

$$
c=\frac{2}{a+b},\qquad
\delta=\frac{b-a}{b+a}.
$$

函数 $ch$ 仍是可行的 $T-1$ 步复合，只需将缩放系数吸收到最后一步；其值域为 $[1-\delta,1+\delta]$。归纳假设于是给出

$$
\delta=\max_{x\in[l,u]}|ch(x)-1|\ge e_{T-1}.
$$

由于贪婪解的中间值域为 $[l_T,u_T]=[1-e_{T-1},1+e_{T-1}]$，上述不等式等价于

$$
[l_T,u_T]\subseteq[ca,cb].
$$

最后一步的多项式在更大的区间上逼近 $1$，最优误差不可能更小。结合连续性和缩放不变性，得到

$$
\begin{aligned}
\max_{x\in[l,u]}|q(h(x))-1|
&=\max_{y\in[a,b]}|q(y)-1|\\
&\ge E(a,b)\\
&=E(ca,cb)\\
&\ge E(l_T,u_T)\\
&=e_T.
\end{aligned}
$$

而 $T$ 步贪婪解恰好达到 $e_T$，故归纳成立。于是

$$
f_T^\ast\circ\cdots\circ f_1^\ast
\in\arg\min_{f_1,\ldots,f_T\in\mathcal P_d^{\mathrm{odd}}}
\max_{x\in[l,u]}|f_T\circ\cdots\circ f_1(x)-1|.
$$
