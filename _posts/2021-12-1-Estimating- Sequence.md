---
title: 'Estimate Sequence and AGD'
toc: true
excerpt_separator: <!--more-->
tags:
  - 优化
  - 加速
---



Nesterov加速方法中的更细公式非常精巧，此处记录Yurii Nesterov 在其 Lectures on Convex Optimization 一书中对该更新公式给出的解释，核心思想在于使用一个序列的收敛性替代函数值的收敛性。



<!--more-->

由于本文主要介绍 [Nesterov加速方法](https://truenobility303.github.io/Nesterov-Acceleration/)中公式来由，不熟悉的读者建议先阅读该方法。

## Estimating Sequence

估计序列的核心思想在于使用一个数列的收敛性替代函数的收敛性，由于在优化算法的收敛性证明时，需要一个函数值的收敛性，这是比较难以分析的。但如果我们得到一个数列，通过数列的收敛性替代函数值的收敛性，可能可以通过数列的构造，控制算法的收敛性。著名的Nesterov加速算法就是使用上述技术达到了某些问题上的最优算法。

估计序列需要一个估计对于收敛性的数列$\lambda_k$和一个控制函数值的函数序列$\phi_k(x)$ ，满足以下关系，


$$
\begin{align*}
\lambda_k &\rightarrow \infty \\
\phi_k(x) &\le (1-\lambda_k) f(x) + \lambda_k \phi_0(x) \\
f(x_k) &\le \phi^{\star}_k = \min \phi_k(x) 
\end{align*}
$$


此时函数值的收敛就转化为数列的收敛，


$$
\begin{align*}
f(x_k) & \le \min \phi_k(x) = \min (1-\lambda_k) f(x)+ \lambda_k \phi_0(x) \le (1-\lambda _k) f(x_{\star}) + \lambda_k \phi_0(x_{\star}) \\
f(x_k ) - f(x_{\star}) & \le \lambda_k ( \phi_o(x_{\star}) - f(x_{\star}))
\end{align*}
$$


在上述观察的基础上，我们使用迭代方法，构造上述的估计序列，


$$
\begin{align*}
\lambda_{k+1} &= (1-\theta_k) \lambda_k \\
\phi_{k+1}(x) &= (1-\theta_k) \phi_k(x) + \theta_k(f(y_k) + \nabla f(y_k)^T(x - y_k) + \frac{\mu}{2} \Vert x- y_k \Vert_2^2)
\end{align*}
$$


容易证明上述构造的序列满足估计序列的性质，利用数学归纳法即可，


$$
\begin{align*}
\phi_{k+1}(x) &= (1-\theta_k) \phi_k(x) + \theta_k(f(y_k) + \nabla f(y_k)^T(x - y_k) + \frac{\mu}{2} \Vert x- y_k \Vert_2^2) \\
& \le (1-\theta_k) \phi_k(x) + \theta_k f(x) \\
& \le (1- \theta_k)((1-\lambda_k) f(x)+ \lambda_{k} \phi_0(x)) + \theta_k f(x) \\
&=(1-\lambda_{k+1}) f(x) + \lambda_{k+1} \phi_0(x)
\end{align*}
$$


为了使得$\phi_k$ 成为$f_k$的一个上界估计，我们希望其最优值$\phi^{\star}$存在闭式解，使用简单的二次函数，我们希望有，


$$
\phi_k(x) = \phi_k^{\star} + \frac{\gamma_k}{2} \Vert x - v_k \Vert_2^2
$$


通过迭代式子，我们计算$\gamma_k, v_k$的更新公式，


$$
\begin{align*}
\phi_{k+1}(x) &= (1-\theta_k) \phi_k(x) + \theta_k(f(y_k) + \nabla f(y_k)^T(x - y_k) + \frac{\mu}{2} \Vert x- y_k \Vert_2^2) \\
\nabla^2 \phi_{k+1}(x) &=  \frac{\gamma_{k+1}}{2}  = \frac{\gamma_k}{2} + \frac{\mu}{2}  \\
\gamma_{k+1} &= (1-\theta_k) \gamma_k + \theta_k\mu \\
\nabla \phi_{k+1}(v_{k+1}) &= (1-\theta_k) \gamma_k(v_{k+1} - v_k)+\theta_k  (\nabla f(y_k) + \mu (v_{k+1} - y_k)) = 0   \\
((1-\theta_k) \gamma_k + \theta_k \mu)v_{k+1} &= (1-\theta_k) \gamma_k v_k + \theta_k \mu y_k - \theta_k \nabla f(y_k) \\
v_{k+1} &= \frac{(1-\theta_k) \gamma_k v_k + \theta_k \mu y_k - \theta_k \nabla f(y_k)}{\gamma_{k+1}} \\
\phi_{k+1}(x) &= \phi_{k+1}^{\star} + \frac{\gamma_{k+1}}{2}\Vert x - v_{k+1} \Vert_2^2 = (1-\theta_k) \phi_k(x) + \theta_k(f(y_k) + \nabla f(y_k)^T(x - y_k) + \frac{\mu}{2} \Vert x- y_k \Vert_2^2)
\end{align*}
$$

以及


$$
\begin{align*}
\phi_{k+1}^{\star} &= (1-\theta_k)(\phi_k^{\star} + \frac{\gamma_k}{2} \Vert v_k - y_k \Vert_2^2) + \theta_kf(y_k)  - \frac{\gamma_{k+1}}{2} \Vert  v_{k+1} - y_k \Vert_2^2\\
&= (1-\theta_k)(\phi_k^{\star} + \frac{\gamma_k}{2} \Vert v_k - y_k \Vert_2^2) + \theta_kf(y_k)  - \frac{\gamma_{k+1}}{2} \Vert  v_{k+1} - y_k \Vert_2^2\\
&= (1-\theta_k)(\phi_k^{\star} + \frac{\gamma_k}{2} \Vert v_k - y_k \Vert_2^2) + \theta_kf(y_k)  - \frac{1}{2 \gamma_{k+1}} \Vert (1-\theta_k) \gamma_k (v_{k} - y_k) - \theta_k \nabla f(y_k) \Vert_2^2\\
&= (1-\theta_k)\phi_k^{\star} + \theta_k f(y_k) - \frac{\theta_k^2}{2\gamma_{k+1} } \Vert \nabla f(y_k) \Vert_2^2+ \frac{(1-\theta_k) \theta_k \gamma_k}{\gamma_{k+1}} \nabla f(y_k)^T(v_k - y_k) + (\frac{(1-\theta_k)\gamma_k}{2} - \frac{(1-\theta_k)^2 \gamma_k^2}{2 \gamma_{k+1}}) \Vert v_k - y_k \Vert_2^2 \\
&= (1-\theta_k)\phi_k^{\star} + \theta_k f(y_k) - \frac{\theta_k^2}{2\gamma_{k+1} } \Vert \nabla f(y_k) \Vert_2^2+\frac{(1-\theta_k) \theta_k \gamma_k}{\gamma_{k+1}} (\nabla f(y_k)^T(v_k - y_k) + \frac{\mu}{2} \Vert y_k - v_k \Vert_2^2) \\
\end{align*}
$$

总结其更新公式，


$$
\begin{align*}
\gamma_{k+1} &= (1-\theta_k) \gamma_k + \theta_k\mu \\
v_{k+1} &= \frac{(1-\theta_k) \gamma_k v_k + \theta_k \mu y_k - \theta_k \nabla f(y_k)}{\gamma_{k+1}} \\
\phi_{k+1}^{\star} 
&= (1-\theta_k)\phi_k^{\star} + \theta_k f(y_k) - \frac{\theta_k^2}{2\gamma_{k+1} } \Vert \nabla f(y_k) \Vert_2^2+\frac{(1-\theta_k) \theta_k \gamma_k}{\gamma_{k+1}} (\nabla f(y_k)^T(v_k - y_k) + \frac{\mu}{2} \Vert y_k - v_k \Vert_2^2) \\
\end{align*}
$$

我们已经几乎得到了最终的答案，下面只需要选取合适的$y_k$即可，



仍然考虑递推地证明，如果我们已经有 $\phi_k^{\star}  \ge f(x_k)$ , 我们利用递推式,


$$
\begin{align*}
\phi_{k+1}^{\star} 
&= (1-\theta_k)\phi_k^{\star} + \theta_k f(y_k) - \frac{\theta_k^2}{2\gamma_{k+1} } \Vert \nabla f(y_k) \Vert_2^2+\frac{(1-\theta_k) \theta_k \gamma_k}{\gamma_{k+1}} (\nabla f(y_k)^T(v_k - y_k) + \frac{\mu}{2} \Vert y_k -v_k \Vert_2^2) \\ 
&\ge (1-\theta_k) f(x_k) + \theta_k f(y_k) - \frac{\theta_k^2}{2\gamma_{k+1} } \Vert \nabla f(y_k) \Vert_2^2+\frac{(1-\theta_k) \theta_k \gamma_k}{\gamma_{k+1}} \nabla f(y_k)^T(v_k - y_k)  \\  
&\ge  (1-\theta_k) (f(y_k) + \nabla f(y_k)^T (x_k - y_k)) + \theta_k f(y_k) - \frac{\theta_k^2}{2\gamma_{k+1} } \Vert \nabla f(y_k) \Vert_2^2+\frac{(1-\theta_k) \theta_k \gamma_k}{\gamma_{k+1}} \nabla f(y_k)^T(v_k - y_k) \\  
&= f(y_k) -\frac{\theta_k^2}{2\gamma_{k+1} } \Vert \nabla f(y_k) \Vert_2^2 + (1-\theta_k) \nabla f(y_k)^T (\frac{\theta_k \gamma_k}{\gamma_{k+1}}(v_k-y_k) +(x_k - y_k))  \\
&= f(y_k) -t_k \Vert \nabla f(y_k) \Vert_2^2 + (1-\theta_k) \nabla f(y_k)^T (\frac{\theta_k \gamma_k}{\gamma_{k+1}}(v_k-y_k) +(x_k - y_k)), \text{Let } t_k = \frac{\theta_k^2}{2\gamma_{k+1}} \\
&\le f(x_{k+1}) + (1-\theta_k) \nabla f(y_k)^T (\frac{\theta_k \gamma_k}{\gamma_{k+1}}(v_k-y_k) +(x_k - y_k)),\text{By Choosing} x_{k+1} = y_k - t_k \nabla y_k ,t_k \le \frac{1}{L} \\
& \le f(x_{k+1}) , \text{Let } \frac{\theta_k \gamma_k}{\gamma_{k+1}}(v_k-y_k) +(x_k - y_k) = 0
\end{align*}
$$


最后一步得到了要使得递推式成立$y_k$的条件，从而我们产生了所需要的估计序列，总结更新公式如下，


$$
\begin{align*}
\gamma_{k+1} &= (1-\theta_k) \gamma_k + \theta_k\mu ,\text{With }\gamma_{k+1} = \frac{\theta_k^2}{2t_k} \\
y_k &= \frac{\gamma_{k+1} x_k + \theta_k \gamma_k v_k}{\gamma_{k+1} + \theta_k \gamma_k}\\
v_{k+1} &= \frac{(1-\theta_k) \gamma_k v_k + \theta_k \mu y_k - \theta_k \nabla f(y_k)}{\gamma_{k+1}} \\
\end{align*}
$$


对其进行进一步化简，就可以Nesterov加速算法的更新公式，实际上仍有一定的计算量，此处暂略。

上述的推导对于近端梯度的场景仍然是成立的，只需要将上面的梯度 $\nabla f(y_k)$ 替换为近端梯度的更新量即可。


$$
\begin{align*}
\gamma_{k+1} &= (1-\theta_k) \gamma_k + \mu \theta_k , \text{Where }\gamma_k = \frac{\theta_{k-1}^2}{t_{k-1}} \\
y &=  x_k + \frac{\theta_k \gamma_k}{\gamma_k + \mu \theta_k} (v_k - x_k) = \frac{\gamma_{k+1}x_k +\theta_k \gamma_k v_k}{\gamma_k+ \mu \theta_k} \\
x_{k+1} &= \text{prox}_{th}(y - t_k \nabla g(y)) \\
v_{k+1} &= x_k + \frac{1}{\theta_k}(x_{k+1}- x_k)
\end{align*}
$$


