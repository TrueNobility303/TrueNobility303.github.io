---
title: '线性老虎机的遗憾下界'
toc: true
excerpt_separator: <!--more-->
tags: 		
  - 强化学习
  - 多臂老虎机
  - 线性老虎机
---

多臂老虎机的遗憾下界推导。

<!--more-->

在上一期的blog中，我们推导了线性老虎机，证明LinUCB算法可以达到 $\tilde{\mathcal{O}}(d \sqrt{T})$ 的遗憾上界。本期我们证明一个匹配的下界。证明参考自课本 [Bandit Algorithms](https://tor-lattimore.com/downloads/book/book.pdf).

在下面的证明中，我们都假设 $T \ge 2d$. 对于动作空间 $\Vert a \Vert \le 1$, 线性老虎机中的奖励为 $x = a^\top \theta + \eta$, 我们令 $\eta$ 服从标准高斯分布, 并且选取 $\theta = \{ \pm \Delta \}^d$. 对于算法的遗憾，我们有

$$
\begin{align*}
{\rm Reg}_T(\theta) =& \Delta \mathbb{E} \left[ \sum_{t=1}^T \sum_{i=1}^d \left( \frac{1}{\sqrt{d}} - a_{ti} {\rm sign}(\theta_i)  \right) \mid \theta \right] \\
\ge & \frac{\Delta \sqrt{d}}{2} \mathbb{E} \left[ \sum_{t=1}^T \sum_{i=1}^d \left( \frac{1}{\sqrt{d}} - a_{ti} {\rm sign}(\theta_i)  \right)^2 \mid \theta \right].
\end{align*}
$$

上述不等式可以直接通过展开平方项验证成立。从上式启发，我们现在定义 $U_i(x) = \sum_{i=1}^{\tau_i} (1/\sqrt{d} - A_{ti} x)^2$, 其中 $x \in \{\pm 1\}$, 而 $\tau_i$ 为一个停时。给定 $i$ 的前提下，对于 $\theta$, 我们定义它的一个对 $\theta'$ 满足 $\theta_i' = - \theta_i$ 而其他分量相同。如果将 $\tau_i$ 定义为第一次 $\sum_{s=1}^t A_{si}^2 \ge T/d$ 的时间，我们知道 $U_i(1) \le (4T/d +2)$. 那么使用下一节中的引理，我们知道当 $p,p'$ 分别为 $U_i(1)$ 在 $\theta, \theta'$ 下的分布时，成立着

$$
\begin{align*}
\vert \mathbb{E}[U_i(1) \mid \theta] - \mathbb{E}[U_i(1) \mid \theta'] \le \left(\frac{4T}{d} +2 \right) \delta(p,p'):= {\rm RHS}.
\end{align*}
$$

令$D(p \Vert p')$表示两个分布间的KL散度，进一步使用Pinksker不等式，我们知道

$$
\begin{align*}
{\rm RHS} \le \left(\frac{4T}{d} +2 \right) \sqrt{\frac{1}{2} D(p \Vert p')}.
\end{align*}
$$

由于 $U_i(1)$ 完全由 $(a_1,x_1,\cdots,a_\tau,x_{\tau_i})$ 所决定，利用数据处理不等式以及KL散度的链式法则，我们可以计算得到

$$
\begin{align*}
{\rm RHS} \le& \left(\frac{4T}{d} +2 \right) \sqrt{\frac{1}{2} \mathbb{E} \left[\sum_{t=1}^{\tau_i}  D(\mathcal{N}(a_t^\top \theta,1) \Vert \mathcal{N}(a_t^\top \theta',1) ) \mid \theta \right]} \\
=&\left(\frac{2T}{d} +1 \right) \sqrt{\mathbb{E} \left[  \sum_{t=1}^{\tau_i}  \langle a_t, \theta - \theta' \rangle^2  \mid \theta \right] } \\
=& \left(\frac{2T}{d} +1 \right) \Delta \sqrt{  \mathbb{E} \left[ \sum_{t=1}^{\tau_i} a_{ti}^2  \mid \theta \right] } \\
\le& \left(\frac{2T}{d} +1 \right) \Delta \sqrt{\frac{T}{d}} \le \frac{5 T \Delta}{2 d} \sqrt{\frac{T}{d}}.
\end{align*}
$$

代入关于RHS的界，我们可以知道

$$
\begin{align*}
& \mathbb{E} [ U_i(1) \mid \theta ] + \mathbb{E} [ U_i(-1) \mid \theta' ] \\
\ge & \mathbb{E} [ U_i(1) \mid \theta' ] + \mathbb{E} [ U_i(-1) \mid \theta' ] - \frac{5 T \Delta}{2 d} \sqrt{\frac{T}{d}} \\
=& 2 \mathbb{E} \left[ \sum_{t=1}^{\tau_i} \left( \frac{1}{d} + a_{ti}^2 \right)  \mid \theta' \right] - \frac{5 T \Delta}{2 d} \sqrt{\frac{T}{d}} \\
\ge& \frac{2T}{d} - \frac{5 T \Delta}{2 d} \sqrt{\frac{T}{d}} = \frac{T}{d},
\end{align*}
$$

其中最后一步依赖于设定 $\Delta =0.4 \sqrt{d/T}$. 最后，我们对所有的 $\theta$ 取平均，并且注意到上述每一对之和的下界，可以得到

$$
\begin{align*}
\sum_{\theta \in \{\pm \Delta \}^d} {\rm Reg}_T(\theta) 
\ge & \frac{\Delta \sqrt{d}}{2} \sum_{\theta \in \{\pm \Delta \}^d} \mathbb{E} \left[ \sum_{t=1}^T \sum_{i=1}^d U_i({\rm sign}(\theta_i)) \mid \theta \right] \\
\ge & \frac{2^d \Delta \sqrt{d} T}{4 } = \frac{2^d d \sqrt{T} }{10},
\end{align*}
$$

说明了至少存在一个 $\theta$ 使得遗憾界至少为 $\Omega(d \sqrt{T})$, 进而证明了上期介绍的LinUCB算法在忽略对数因子下为最优算法。

## A Useful Lemma for TV Distance

我们给出一个有用的概率引理, 来自课本的练习14.4. 令 $\delta(p,q)$ 为两个分布的TV距离，对于任意的 $f(x): X \rightarrow [a,b]$, 成立

$$
\begin{align*}
\left \vert  \int_X f(x) p(x)  {\rm d} x - \int_X f(x) q(x) {\rm d} x   \right \vert \le (b-a) \delta(p,q).
\end{align*}
$$

定义 $g(x) = f(x) - (a+b)/2$. 我们知道 $\vert g(x) \vert \le (b-a)/2$. 所以

$$
\begin{align*}
{\rm LHS} =& \left \vert  \int_X g(x)p(x) {\rm d}x  - \int_X g(x) q(x) {\rm d} x   \right \vert  \\
\le& \frac{b-a}{2}  \int_X \vert p(x) - q(x) \vert {\rm d} x = {\rm RHS}.
\end{align*} 
$$



