---
title: '神经网络中的优化器'
toc: true
excerpt_separator: <!--more-->
tags:
  - 优化
  - 非凸优化
  - 神经网络优化
---

本文整理神经网络中的著名优化器的推导过程，包括自适应算法Adagrad、Adam，以及为了节省内存所提出的Adafactor、GaLore及相关优化器，以及近期关注度很高的在特征空间上作用的新型优化器Shampoo、SOAP以及Muon, 试图用简单的方式总结其联系。

<!--more-->

## Adagrad (matrix version)



考虑在线学习的设定，给定一系列凸函数 $\{f_t\}$, 希望生成序列 $\{ x_t\}$ 最小化如下的遗憾界

$$
\begin{align*}
R_T = \sum_{t=1}^T f_t(x_t) - \inf_x \sum_{t=1}^T f_t(x).
\end{align*}
$$

在线梯度下降算法产生序列 $x_{t+1} = x_t - \eta g_t$，其中 $g_t \in \partial f_t(x_t)$, 在最坏情况下该算法的遗憾界不能被改进。 Adagrad 算法 [1] 希望寻找到一个正定矩阵 $S_t$ 作为 preconditioner, 使得预处理后的梯度迭代 $x_{t+1} = x_t - \eta S_t^{-1} g_t$, 在某些场景下具有更好的保证。我们从算法的证明中得到最优的 $S_t$, 我们定义自适应范数 $\Vert \,\cdot\, \Vert_{S_t}$, 这使得下面分析中正好出现凸性对应的内积项，那么

$$
\begin{align*}
&\quad \Vert x_{t+1} - x^\ast \Vert_{S_{t}}^2 \\ 
&= \Vert x_{t} - \eta S_t^{-1} g_t -  x^\ast \Vert_{S_{t}}^2 \\
&= \Vert x_t - x^\ast \Vert^2_{S_{t}} - \langle \eta g_t, x_t - x^\ast \rangle +  \eta^2 \langle g_t, S_t^{-1} g_t \rangle \\
&\le  \Vert x_t - x^\ast \Vert^2_{S_{t-1}} - \eta (f_t(x_t) -f_t(x^\ast))  + \eta^2 \langle g_t, S_t^{-1} g_t \rangle + \Vert x_t - x^\ast \Vert^2_{S_{t} - S_{t-1}} \\
&\le  \Vert x_t - x^\ast \Vert^2_{S_{t-1}} - \eta (f_t(x_t) -f_t(x^\ast))  + \eta^2 \langle g_t, S_t^{-1} g_t \rangle + \Vert x_t - x^\ast \Vert^2 \lambda_{\max} (S_t - S_{t-1}) \\
&\le \Vert x_t - x^\ast \Vert^2_{S_{t-1}} - \eta (f_t(x_t) -f_t(x^\ast))  + \eta^2 \langle g_t, S_t^{-1} g_t \rangle + \Vert x_t - x^\ast \Vert^2 {\rm tr} (S_t - S_{t-1}).
\end{align*}
$$

令 $D = \sup_t \Vert x_t - x^\ast \Vert$。或者使用投影梯度法可以使得所有迭代都位于一个紧集上，我们同样可以记紧集的大小为 $D$. 对所有的 $t = 1,\cdots,T$ 进行求和，可以得到遗憾界的上界为

$$
\begin{align*}
R_T \le \eta^{-1} D^2 ({\rm tr}(S_T) +1) + \eta \sum_{t=1}^T \langle g_t ,S_t^{-1} g_t \rangle.
\end{align*}
$$

从上式启发，我们限制 $S$ 有限得约束下进行求解：

$$
\begin{align*}
S_t = \arg \min_{S} \sum_{t=1}^T  \langle g_t, S^{-1} g_t \rangle , \quad {\rm s.t.} \quad S \succeq O, ~ {\rm tr}(S) \le c.
\end{align*}
$$


使用Largange乘子法求解该问题，令 $A = \sum_{t=1}^T g_t g_t^\top$, Largange函数为

$$
\begin{align*}
L(S,\theta, Z )= {\rm tr}(S^{-1} A) + \theta ({\rm tr}(S) -c) + tr(S Z).
\end{align*}
$$

一阶最优性条件为

$$
\begin{align*}
-S^{-1} A S^{-1} + \theta I + Z = O.
\end{align*}
$$ 

当 $S \succ O$ 此时根据互补松弛条件 $Z = O$, 那么上式可以推出 $S \propto A^{-1/2}$. 当 $S$ 奇异的时候可以通过取极限的方式得到，为了便于读者理解此处暂略。至此就推出了矩阵形式的Adagrad算法:

$$
\begin{align*}
x_{t+1} = x_t - \eta \left(\sum_{t=1}^T g_t g_t^\top \right)^{-1/2} g_t.
\end{align*}
$$

## Adagrad (diagonal version), RMSProp, and Adam

在大规模训练场景中，矩阵形式的Adagrad 算法由于需要维护 $d \times d$ 的矩阵，实际难以接受。因此原文也提出了使用对角近似的方法$S_t = {\rm diag}(s_t)$. 通过完全相同的推导后，最优的 $s_t$ 可以通过求解如下的优化问题给出：

$$
\begin{align*}
s_t = \arg \min_{s} \sum_{t=1}^T \langle g_t , {\rm diag}(s)^{-1} g_t \rangle, \quad {\rm s.t.} \quad S \succeq 0, ~ \langle 1, s \rangle \le c.
\end{align*}
$$

类似上面的矩阵优化问题，可以通过Lagrange函数的KKT条件求解得到 $s_t \propto \sqrt{ \sum_{t=1}^T g_t \odot g_t }$. RMSProp 算法 [2] (出现在Hinton的slides中）的改动在于上述的二阶矩使用带权的滑动平均计算，而不是uniform的权重，这可以给更近的梯度赋予更高的权重，更好地利用近期的信息：

$$
\begin{align*}
v_t = \beta_2 v_{t-1} + (1- \beta_2) g_t \odot g_t.
\end{align*}
$$


而著名的Adam算法 [3] 进一步加入动量更新，相当于带动量的SGD与RMSProp进行结合：

$$
\begin{align*}
m_t = \beta_1 m_{t-1} + (1- \beta_1) g_t.
\end{align*}
$$

然后更新为 $x_{t+1} = x_t - \eta m_t / (\sqrt{v_t} + \epsilon)$.

## Adafactor and Adalayer/Adam-Mini 

对于很大规模的训练，Adam由于需要同时存储一阶矩以及二阶矩，会造成比SGD更大的存储开销，Adafactor 希望利用矩阵分解的方式降低二阶矩得存储开销. 在神经网络中的参数通常为矩阵, 因此 $V_t \in \mathbb{R}^{m \times n}$, Adafactor 算法[4] 希望计算一个近似的秩1逼近，得到列向量 $c_t \in \mathbb{R}^{m \times 1}$ 和行向量 $r_t \in \mathbb{R}^{1 \times n}$ 使得 $V_t \approx c_t r_t$.  最直接的想法是根据Eckart-Young 定理求出最大的左右奇异向量作为Frobenius范数意义下的最佳秩一逼近，但此时得到的估计不能保证是非负的我们无法计算其平方根。为此，Adafactor原文考虑如下的矩阵散度意义下的非负矩阵分解问题，

$$
\begin{align*}
(c_t,r_t) = \arg \min_{c,r} d(V_t, c r), \quad {\rm s.t.} \quad c \succeq 0, ~r \succeq 0,
\end{align*}
$$

其中矩阵散度定义为 $d(P,Q) = \sum_{i=1}^m \sum_{j=1}^n P_{ij} \log ({P_{ij}}/{Q_{ij}}) - P_{ij}  +Q_{ij}$. 经过简单的计算，可以验证得到上述优化问题的 (一个）解为

$$
\begin{align*}
c_t = V_t 1_n, \quad r_t = 1_m^\top V_t / 1_m^\top V_t  1_n.
\end{align*}
$$

由上述分析启发，我们可以得到Adafactor算法对于矩阵梯度 $G_t = \nabla f(X_t)$的二阶矩更新：

$$
\begin{align*}
c_t &= \beta_2  c_{t-1} + (1- \beta_2) (G_t \odot G_t) 1_n, \\
r_t &= \beta_2 r_{t-1} + (1- \beta_2) 1_m^\top (G_t \odot G_t).
\end{align*}
$$

相当于仅在一行一列上面记录。近期工作（Adalayer [5] / Adam-mini [6]) 为了进一步节省内存，提出每个layer只存储一个参数 $s_t$ 的更粗的存储方法，直接记录每一块参数二阶矩的均值，更新为

$$
\begin{align*}
s_t = \beta_2 s_{t-1} + (1- \beta_2)  {\rm mean} (G_t \odot G_t).
\end{align*}
$$

## Shampoo2 and SOAP 

Shampoo2 (这里2表示square的含义) 与原始的Shampoo [8] 稍有不同，且提出较晚 [7], 但有实验表明似乎有更优的效果，而且 [7] 从秩一逼近Adagrad的角度推出该算法，具有较强的理论基础，因此这里先选择作为介绍。

对于矩阵参数，Adagrad希望近似的preconditioner为 $S_t = \sum_{t=1}^T G_t \otimes G_t$, 这里使用kronecker积表示矩阵的外积。利用恒等式 $A \otimes B = {\rm vec}(B) {\rm vec}(A)^\top$，最佳秩一逼近问题可以表示为

$$
\begin{align*}
(L_t, R_t) = \arg \min_{L, R} \Vert S_t - R \otimes L \Vert_F^2.
\end{align*}
$$

根据Eckart-Young 定理该问题的显式解由最大的奇异向量给出。但是实际上 $S_t$ 作为一个 $mn \times mn$ 大小的矩阵，我们不可能对其进行奇异值分解。我们采用如下幂法近似计算：

$$
\begin{align*}
{\rm vec}(L_k) = S_t ~ {\rm vec}(R_{k-1}), \quad {\rm vec}(R_k) = S_t^\top ~ {\rm vec}(L_{k-1}).
\end{align*}
$$

注意到我们这里没有进行归一化，因此归一化因子可以被吸收在更新的步长 $\eta$ 中。带入 $S_t = \sum_{t=1}^T G_t \otimes G_t$ 并且利用恒等式 $(B \otimes A) {\rm vec}(X) = {\rm vec}(AXB^\top$), 有

$$
\begin{align*}
L_k = \sum_{t=1}^T G_t R_{k-1} G_t^\top, \quad R_k = \sum_{t=1}^T G_t^\top L_{k-1} G_t.
\end{align*}
$$

[7] 中声称初始化为 $(L_0,R_0) = (I,I)$ 的单步幂法 $k=1$ 就已经足够，而这可以推出Shampoo2的迭代，注意到上式中的平均仍然使用滑动的方法更新，得到

$$
\begin{align*}
L_t &= \beta_2 L_{t-1} + (1 - \beta_2) G_t G_t^\top, \\
R_t &= \beta_2 R_{t-1} + (1- \beta_2) G_t^\top G_t.
\end{align*}
$$

得到 $(L_t,R_t)$ 之后我们用其近似Adagrad更新：

$$
\begin{align*}
{\rm vec}(X_t) = {\rm vec}(X_{t-1}) - \eta (R_t \otimes L_t)^{-1/2} {\rm vec}(G_t).
\end{align*}
$$

同样利用恒等式 $(B \otimes A) {\rm vec}(X) = AXB^\top$ 得到矩阵形式的更新为

$$
\begin{align*}
X_t = X_{t-1} - \eta  L_t^{-1/2} G_t R_t^{-1/2}.
\end{align*}
$$

这就从幂法的角度推导出了Shampoo2的迭代，其中使用的preconditioner为原Shampoo算法中的平方，[9] 从另一个角度理解上面的更新，证明Shampoo2完全等价于在特征空间上的 Adafactor 算法。定义 $U_t$ 和 $V_t$ 分别为 $L_t$ 以及 $R_t$对应的特征向量，然后定义旋转后的梯度 $G_t' = U_t ^\top G_t V_t$. 现在考虑对于 $G_t'$ 作用 $\beta_2=0$ 的Adafactor算法，得到的分解后的行向量和列向量为

$$
\begin{align*}
c_t' = (G_t' \odot  G_t') 1_n, \quad r_t' = 1_m^\top (G_t' \odot G_t').
\end{align*}
$$

用分解后的二阶矩计算得到在特征空间内的更新方向, 然后再旋转回原空间进行更新，得到

$$
\begin{align*}
X_t =  X_t - U_t (G_t' / \sqrt{c_t r_t}) V_t^\top.
\end{align*}
$$

上述过程就描述了 [9] 中所说的在特征空间上的Adafactor算法，下面我们验证该算法与Shampoo2的等价性。利用奇异值分解 $G_t = U_t \Sigma_t V_t^\top$, 可以得到

$$
\begin{align*}
&\quad L_t^{-1/2} G_t R_t^{-1/2} \\
&= (G_t G_t^\top)^{-1/2} G_t (G_t^\top G_t)^{-1/2} \\
&=(U_t \Sigma_t^{-1} U_t^\top) (U_t \Sigma_t V_t^\top) (V_t \Sigma_t^{-1} V_t^\top) = U_t \Sigma_t^{-1} V_t^\top.
\end{align*}
$$

容易知道旋转过后 $G_t' = \Sigma_t$, 因此 

$$
\begin{align*}
U_t (G_t' / \sqrt{c_t r_t}) V_t^\top = U_t \Sigma_t^{-1} V_t^\top.
\end{align*}
$$

上述推导十分简单，但却可以给我们一个启发，既然Shampoo2相当于在特征空间上使用Adafactor算法，是否可以考虑也在特征空间上使用Adam算法呢，这将导出最终的SOAP算法 [9]。事实上，上述想法并非原创，该想法在Galore算法中 [10] 事实上就已经出现，但两者的区别在于GaLore的目的是为了投影到特征空间后节省内存开销，而SOAP的目的是发明更高效的算法，因此GaLore中将动量也投影到对应的子空间内用Adam进行更新，而SOAP的动量并不进行投影。第二个区别是SOAP使用了Shampoo中相同的滑动平均计算 $(L_t,R_t)$ , 而GaLore更直接地直接投影到了当前的子空间上，同样处于节省内存的考量。SOAP的文章中也展现出了令人满意的实验结果。

## Shampoo and Muon 

从特征空间的角度出发，如果在特征空间上使用SignSGD算法，就可以自然推出Muon / Shampoo算法。Shampoo算法 [8] 的 $(L_t,R_t)$ 采用相同的方式更新，不同的是迭代为

$$
\begin{align*}
X_t = X_{t-1} - \eta L_t^{-1/4} G_t R_t^{-1/4}.
\end{align*}
$$

原文类似于Adagrad论文的推导，可以得到遗憾界，其大小依赖于 $L_T$ 以及 $R_T$. 在Shampoo的迭代中取 $\beta_2 = 0$, 并且利用奇异值分解 $G_t = U_t \Sigma_t V_t^\top$, 可以得到

$$
\begin{align*}
&\quad L_t^{-1/4} G_t R_t^{-1/4} \\
&= (G_t G_t^\top)^{-1/4} G_t (G_t^\top G_t)^{-1/4} \\
&=(U_t \Sigma_t^{-1/2} U_t^\top) (U_t \Sigma_t V_t^\top) (V_t \Sigma_t^{-1/2} V_t^\top) = U_t V_t^\top.
\end{align*}
$$

这相当于把奇异值都设置成1，对矩阵进行whitening操作。上式的另一个含义是取 $G_t$ 的最佳正交逼近，也即

$$
\begin{align*}
U_t V_t = \arg \min_{OO^\top = I } \Vert O - G_t \Vert_F^2.
\end{align*}
$$

Muon 算法 [11] 使用了一个不动点迭代算法求解上述的最佳正交逼近问题，实际表明这个迭代算法仅需要5步就可以取得了相较于Adam优化器显著更优的效果。更具体地，使用如下迭代

$$
\begin{align*}
\hat G_t^1 &= a G_t + b (G_t G_t^\top) G_t + c (G_t G_t^\top)^2 G_t \\
&= U_t (a \Sigma_t + b \Sigma_t^3 + c \Sigma_t^5 ) V_t^\top = U_t \phi(\Sigma_t) V_t^\top.
\end{align*}
$$

不断重复上述迭代 k次，得到

$$
\begin{align*}
\hat G_t^k = U_t \phi^k(\Sigma_t) V_t^\top.
\end{align*}
$$

选取系数 $(a,b,c)$ 使得 $\phi^k(x) \rightarrow 1$ 即可完成目的。

## Reference 

[1] Duchi, John, Elad Hazan, and Yoram Singer. "Adaptive subgradient methods for online learning and stochastic optimization." JMLR, 2011 (extention of COLT, 2010).

[2] Hinton, Geoffrey, Nitish Srivastava, and Kevin Swersky. "Neural networks for machine learning lecture 6a overview of mini-batch gradient descent." Technical report, 2012.

[3] Kingma, Diederik P., and Jimmy Ba. "Adam: A method for stochastic optimization." In ICLR, 2015.

[4] Shazeer, Noam, and Mitchell Stern. "Adafactor: Adaptive learning rates with sublinear memory cost." In ICML, 2018.

[5] Zhao, Rosie, Depen Morwani, David Brandfonbrener, Nikhil Vyas, and Sham M. Kakade. "Deconstructing What Makes a Good Optimizer for Autoregressive Language Models." In ICLR, 2025.

[6] Zhang, Yushun, Congliang Chen, Ziniu Li, Tian Ding, Chenwei Wu, Diederik P. Kingma, Yinyu Ye, Zhi-Quan Luo, and Ruoyu Sun. "Adam-mini: Use fewer learning rates to gain more." In ICLR, 2025. 

[7] Morwani, Depen, Itai Shapira, Nikhil Vyas, Eran Malach, Sham Kakade, and Lucas Janson. "A New Perspective on Shampoo's Preconditioner." In ICLR, 2025. 

[8] Gupta, Vineet, Tomer Koren, and Yoram Singer. "Shampoo: Preconditioned stochastic tensor optimization." In ICML, 2018.

[9] Vyas, Nikhil, Depen Morwani, Rosie Zhao, Mujin Kwun, Itai Shapira, David Brandfonbrener, Lucas Janson, and Sham Kakade. "SOAP: Improving and stabilizing shampoo using adam." In ICLR, 2025.

[10] Zhao, Jiawei, Zhenyu Zhang, Beidi Chen, Zhangyang Wang, Anima Anandkumar, and Yuandong Tian. "GaLore: Memory-Efficient LLM Training by Gradient Low-Rank Projection." In ICML, 2024.

[11] https://github.com/KellerJordan/Muon.
