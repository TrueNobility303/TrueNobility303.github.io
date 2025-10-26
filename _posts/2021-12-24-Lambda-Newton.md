---
title: 'Regularized Newton Method'
toc: true
excerpt_separator: <!--more-->
tags:
  - 优化
  - 牛顿法
---



Paper Reading: Regularized Newton Method with Global $\mathcal{O}(1/k^2)$ Convergence.

<!--more-->



## Regularized Newton Method



该论文的主要思想是研究三次正则Newton方法中的正则项系数，首先回顾三次正则Newton方法，或者称为Cubib Newton方法, 文章仅考虑优化函数为凸函数的情形，


$$
\begin{align}
x_{k+1} &= \text{argmin}_y [ \langle \nabla f(x_k) , y-x_k \rangle +\frac{1}{2} \langle \nabla^2 f(x_k) (y-x_k) , y-x_k \rangle + \frac{M}{6} \Vert y - x_k \Vert^3 ] \\
&=x_k - (\nabla^2f(x_k)+ Lr_k )^{-1} \nabla f(x_k), \text{Let }r_k = \Vert x_{k+1} - x_k \Vert, M = 2L 
\end{align}
$$


同样与Cubic Newton方法相同，假设Hessian矩阵Lipschitz连续的条件满足，


$$
\begin{align}
\Vert \nabla^2 f(x) - \nabla^2 f(y) \Vert \le 2L \Vert x - y\Vert
\end{align}
$$




同理在上述假设下，利用带积分余项的Taylor展开可以得到，


$$
\begin{align}
\Vert \nabla f(y) - \nabla f(x)  - \nabla^2 f(x) (y-x) \Vert &\le L\Vert y - x \Vert^2 \\
\vert f(y) - f(x) - \langle \nabla f(x) , y-x \rangle - \frac{1}{2} \langle \nabla^2 f(x) (y-x) , y-x \rangle &\le \frac{L}{3} \Vert y - x \Vert^3
\end{align}
$$






该论文提出了一种新的正则化方式，从而避免了上述的 $r_k$ 的计算，而且可以达到相同的收敛率。


$$
\begin{align}
x_{k+1} &= x_k - (\nabla^2f(x_k)+ \lambda_k I )^{-1} \nabla f(x_k) ,\text{Let } \lambda_k = \sqrt{L \Vert \nabla f_k(x_k)\Vert}\\
\end{align}
$$
 

该正则项的选取实际上为了满足和Cubic Newton方法近似的结果，


$$
\lambda_k \approx L r_k
$$


我们后面会看到该论文是如何将该想法变为现实，并且证明收敛率的, 我们首先希望对 $r_k$ 进行观察， 可以发现正则项系数 $\lambda_k$ 更大与我们的需求，


$$
\begin{align}
Lr_k &= L\Vert x_{k+1} - x_k \Vert \\
&=L\Vert (\nabla^2 f(x_k) + \lambda_k I)^{-1} \nabla f(x_k) \Vert \\
&\le \frac{L\Vert \nabla f(x_k) \Vert}{\lambda_k}  \\
&= \lambda_k 
\end{align}
$$


首先重写迭代更新公式, 进一步可以得到更新后的函数值的梯度的界, 


$$
\begin{align}
\lambda_k (x_{k+1} - x_{k}) &=  -\nabla f(x_k) - \nabla^2 f(x_k) (x_{k+1} -x_k) \\
\Vert \nabla f(x_{k+1}) \Vert &= \Vert \nabla f(x_{k+1}) -\nabla f(x_k) - \nabla^2 f(x_k) (x_{k+1} -x_k) - \lambda_k (x_{k+1}-x_k) \Vert\\
&\le \Vert \nabla f(x_{k+1}) -\nabla f(x_k) - \nabla^2 f(x_k) (x_{k+1} -x_k) \Vert+ \lambda_k \Vert (x_{k+1}-x_k) \Vert \\
&\le L \Vert x_{k+1}- x_k \Vert^2 +  \lambda_k \Vert (x_{k+1}-x_k) \Vert \\ 
&= L r_k^2 + \lambda_k r_k \\
&\le 2\lambda_k r_k \\
&\le \frac{2\lambda_k^2}{L} \\
&=2\Vert \nabla f(x_k) \Vert
\end{align}
$$


该结论虽然证明简单，但设计精巧，对于算法的收敛性证明具有关键性作用。



## Coverage

利用函数的三次上界，可以发现该算法可以使得函数值是一个下降序列，且可以出现Cubic Newton方法的下降序列类似的形式，


$$
\begin{align}
f(x_{k+1}) - f(x_k) &\le \nabla f(x_k)^T (x_{k+1}-x_k) + \frac{1}{2} (x_{k+1}-x_k)^T \nabla^2 f(x_k) (x_{k+1} -x_k) + \frac{L}{3} r_k^3 \\
&=(\nabla f(x_k)+ \nabla^2 f(x_k) (x_{k+1} - x_k))^T (x_{k+1}-x_k) - \frac{1}{2} (x_{k+1}-x_k)^T \nabla^2 f(x_k) (x_{k+1} -x_k) + \frac{L}{3} r_k^3 \\
&\le (\nabla f(x_k)+ \nabla^2 f(x_k) (x_{k+1} - x_k))^T (x_{k+1}-x_k) + \frac{L}{3} r_k^3 \\
&=-\lambda_k (x_{k+1} - x_k))^T (x_{k+1}-x_k) + \frac{L}{3} r_k^3 \\
&=-\lambda_k r_k^2 +\frac{L}{3} r_k^3 \\
&\le -\frac{2L \lambda_k r_k^2}{3} 
\end{align}
$$


在凸函数的前提下，我们像Cubic Newton方法相同地假设有 $\Vert x_k - x_{\star} \Vert \le D$, 可以得到函数值和梯度的联系，


$$
\begin{align}
f(x_k) - f(x_{\star}) &\le \nabla f(x_k)^T(x_k - x_{\star}) \\
&\le \Vert \nabla f(x_k) \Vert  \Vert x_k - x_{\star} \Vert \\
&\le D\Vert \nabla f(x_k) \Vert \\
\end{align}
$$


为了得到最后的结论，需要考虑下面的集合，该定义衡量了梯度的平稳性，


$$
\mathcal{I_k}: \Vert \nabla f(x_{i+1}) \Vert \le \frac{1}{4}\Vert \nabla f(x_i) \Vert, i\le k
$$


考虑简单的情况，如果超过半数的元素满足上述的平稳性，也即 $\vert \mathcal{I}_k \vert \ge \frac{k}{2} $, 可以显然地得到指数级的收敛率，


$$
\begin{align}
f(x_k) - f(x_{\star}) &\le D\Vert \nabla f(x_k) \Vert \\
&\le D (\frac{1}{4})^{\frac{k}{2}} 2^{\frac{k}{2}} \Vert \nabla f(x_0 ) \Vert \\
&=D (\frac{1}{2})^{\frac{k}{2}} \Vert \nabla f(x_0 ) \Vert
\end{align}
$$


据此我们看到了平稳性梯度对收敛性的重要贡献，但下面更神奇的是，我们将会看到，即便所有的元素都不满足梯度更新的平稳性假设，由于算法本身保证了 $\Vert \nabla f(x_{k+1}) \Vert  \le 2 \Vert \nabla f(x_{k}) \Vert$,  可以得到类似Cubic Newton法中的结论，从而保证类似的收敛性，


$$
\begin{align}
\frac{1}{4} \Vert \nabla f(x_{k}) \Vert &\le \Vert \nabla f(x_{k+1}) \Vert \le 2\lambda_k r_k  = 2r_k\sqrt{L \Vert \nabla f(x_k) \Vert} \\ 
\Vert \nabla f(x_{k}) \Vert &\le  64 r_k^2 L \\
f(x_{k+1}) - f(x_k) &\le -\frac{2L \lambda_kr_k^2}{3}  \\
&\le -\frac{L} {96} (\frac{f(x_k) - f(x_{\star})}{LD})^{\frac{3}{2}} \\
&= -\tau (f(x_k) - f(x_{\star}))^{\frac{3}{2}}, \text{Let } \tau = \frac{1}{96 L^{\frac{1}{2}} D^{\frac{3}{2}} }
\end{align}
$$


类似Cubic Newton中的证明思路，定义辅助量 $\alpha_k$ ,


$$
\begin{align}
\alpha_k^2 &:= \tau^2 (f(x_k) - f(x_{\star})) \\
\alpha_{k+1}^2 &=\tau^2 (f(x_{k+1}) - f(x_{\star})) \\
&=\tau^2 (f(x_{k}) - f(x_{\star}) -\tau (f(x_k) - f(x_{\star}))^{\frac{3}{2}}) \\
&= \alpha_k^2 - \alpha_k^3 \\
\end{align}
$$


据此得到了Cubic Newton中的序列的形式，回顾该序列的性质，


$$
\begin{align}
\alpha_{k+1}^2 &\le \alpha_k^2 - \frac{2}{3} \alpha_k^3  \le \alpha_k^2\\
\frac{1}{\alpha_{k+1}} - \frac{1}{\alpha_k} &= \frac{\alpha_k - \alpha_{k+1}}{\alpha_{k+1}\alpha_k} =\frac{\alpha_k^2 - \alpha_{k+1}^2}{\alpha_{k+1}\alpha_k(\alpha_k + \alpha_{k+1})} \ge \frac{1}{3} \\
\frac{1}{\alpha_k} &\ge \frac{1}{\alpha_0}+\frac{k}{3} =1+ \frac{k}{3} ,\text{Let } \alpha_0 = 1\\
\alpha_k^2  &\le \frac{1}{(1+\frac{k}{3})^2} \\
\tau^2 (f(x_k) - f(x_{\star})) &\le  \frac{1}{(1+\frac{k}{3})^2}\\
(f(x_k) - f(x_{\star})) &\le  \frac{1}{(1+\frac{k}{3})^2 \tau^2}\\
\end{align}
$$


值得注意的是，尽管上述的证明仅仅针对于所有元素都不满足平稳梯度的情况，但对于更好的情况，上述证明显然是成立的，只需要根据函数值序列的单调性，选取满足平稳梯度的一个子列进行证明即可，本质上并无区别。


最终得到了算法在凸函数上具有收敛率为 $\mathcal{O}(1/k^2)$ , 优于梯度法的  $\mathcal{O}(1/k)$，且该算法形式非常简单，并不需要像Cubic Newton一样解一个三次正则方程得到的子问题，只需要根据迭代公式计算正则项 $\lambda_k$ 即可。


