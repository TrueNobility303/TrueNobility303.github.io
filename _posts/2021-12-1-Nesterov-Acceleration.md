---
title: 'AGD: Accelerated Gradient Descent'
toc: true
excerpt_separator: <!--more-->
tags:
  - 优化
  - 加速
---

Nesterov对于近端梯度下降的加速方法以及详细证明。



<!--more--> 

## Nesterov‘s Method



Nesterov加速方法采用如下精巧的更新公式，对普通的近端梯度下降算法进行加速，使其收敛速度达到该问题最优，


$$
\begin{align}
\gamma_{k+1} &= (1-\theta_k) \gamma_k + \mu \theta_k , \text{Where }\gamma_k = \frac{\theta_{k-1}^2}{t_{k-1}} \\
y &=  x_k + \frac{\theta_k \gamma_k}{\gamma_k + \mu \theta_k} (v_k - x_k) = \frac{\gamma_{k+1}x_k +\theta_k \gamma_k v_k}{\gamma_k+ \mu \theta_k} \\
x_{k+1} &= \text{prox}_{th}(y - t_k \nabla g(y)) \\
v_{k+1} &= x_k + \frac{1}{\theta_k}(x_{k+1}- x_k)
\end{align}
$$


类似于 [近端梯度下降方法](https://truenobility303.github.io/Sub-Proximal/)  中的证明， 或者因为FISTA为Nesterov加速方法的特例，也可以参见 [FISTA](https://truenobility303.github.io/FISTA/) 中的部分


$$
\begin{align}
f(x_{k+1}) &= f(y_k - t_k G_k) \\
&= g(y_k - t_k G_k) + h(y_k - t_k G_k) \\
&\le g(y_k) - t_k \nabla g(y_k)^T G_k + \frac{Lt_k^2}{2} \Vert G_k \Vert_2^2 + h(y_k - t_k G_k) ,\text{ By Lipschitz}\\ 
&\le  g(y_k) - t_k \nabla g(y_k)^T G_k + \frac{Lt_k^2}{2} \Vert G_k \Vert_2^2  + h(x_k) - (G_k - \nabla g(y_k))^T(x_k - y_k + t_k G_k ) ,\text{By Sub Gradient} \\
&\le g(y_k) +h(x_k)+( \nabla g(y_k)- G_k)^T(x_k - y_k) + (\frac{L t_k^2}{2}- t_k) \Vert G_k \Vert_2^2  \\
&\le g(x_k) +h(x_k) - G_k^T(x_k - y_k) + (\frac{L t_k^2}{L}- t_k) \Vert G_k \Vert_2^2 - \frac{\mu}{2} \Vert x_k - y_k \Vert_2^2 , \text{ By Stronly Convexity of } g \\
&=f(x_k )  + G_k^T(y_k - x_k) + (\frac{L t_k^2}{2}- t_k) \Vert G_k \Vert_2^2 -\frac{\mu}{2} \Vert x_k - y_k \Vert_2^2
\end{align}
$$



类似地对于最优值，



$$
\begin{align}
f(x_{k+1})&= f(y_k - t_k G_k) \\
&= g(y_k - t_k G_k) + h(y_k - t_k G_k) \\
&\le g(y_k) - t_k \nabla g(y_k)^T G_k + \frac{Lt_k^2}{2} \Vert G_k \Vert_2^2 + h(y_k - t_k G_k) ,\text{ By Lipschitz}\\ 
&\le  g(y_k) - t_k \nabla g(y_k)^T G_k + \frac{Lt_k^2}{2} \Vert G_k \Vert_2^2  + h(x_{\star}) - (G_k - \nabla g(y_k))^T(x_{\star} - y_k + t_k G_k ) ,\text{By Sub Gradient} \\
&\le g(y_k) +h(x_{\star})+( \nabla g(y_k)- G_k)^T(x_{\star} - y_k) + (\frac{L t_k^2}{2}- t_k) \Vert G_k \Vert_2^2  \\
&\le g(x_{\star}) +h(x_{\star}) - G_k^T(x_{\star} - y_k) + (\frac{L t_k^2}{L}- t_k) \Vert G_k \Vert_2^2 - \frac{\mu}{2} \Vert x_{\star} -y_k \Vert_2^2, \text{ By Strongly Convexity of } g \\
&=f(x_{\star})  + G_k^T( y_k-x_{\star}) + (\frac{L t_k^2}{2}- t_k) \Vert G_k \Vert_2^2 - \frac{\mu}{2} \Vert x_{\star} -y_k \Vert_2^2
\end{align}
$$



同样地进行凸组合可以得到，


$$
\begin{align}
f(x_{k+1}) \le (1-\theta_k)f(x_k) + \theta_k f(x_{\star}) + G_k^T(y_k - (1-\theta_k) x_k - \theta_k x_{\star}) - \frac{(1-\theta_k) \mu}{2} \Vert   x_k - y_k \Vert_2^2 - \frac{\theta_k \mu}{2} \Vert x_{\star} - y_k \Vert_2^2 + (\frac{L t_k^2}{2}-t_k) \Vert G_k\Vert_2^2\\
\end{align}
$$




进行一番精巧绝伦的计算，


$$
\begin{align}
v_{k+1} &= x_k + \frac{1}{\theta_k} (y - t_k G_k - x_k) \\
&= \frac{1}{\theta_k}(y - (1-\theta_k) x_k) - \frac{1}{\theta_k} t_kG_k \\
&= \frac{1}{\theta_k}(y - (1-\theta_k) x_k) - \frac{\theta_k}{\gamma_{k+1}}G_k \\
\gamma_{k+1} v_{k+1} &= \frac{\gamma_{k+1}}{\theta_k}(y - (1-\theta_k) x_k) - \theta_k G_k \\
&= \frac{1}{\theta_k}(\gamma_{k+1} y - (1-\theta_k) \gamma_{k+1} x_k) - \theta_k G_k \\
&= \frac{1}{\theta_k}((\gamma_k + \mu \theta_k - \theta_k \gamma_k) y - (1-\theta_k) \gamma_{k+1} x_k) - \theta_k G_k \\
&= \frac{1-\theta_k}{\theta_k}((\gamma_k + \mu \theta_k)y - \gamma_{k+1}x_k) + \theta_k \mu y - \theta_k G_k\\
&= \frac{1-\theta_k}{\theta_k}(\theta_k \gamma_k v_k) + \theta_k \mu y - \theta_k G_k\\
&=(1-\theta_k) \gamma_k v_k + \theta_k \mu y - \theta_k G_k\\
&=(\gamma_{k+1} - \mu \theta_k) v_k+ \mu \theta_k y - \theta_k G_k\\
&= \gamma_{k+1}v_k + \mu \theta_k(y-v_k) - \theta_k G_k \\
\end{align}
$$



更进一步，



$$
\begin{align}
\frac{\gamma_{k+1}}{2} \Vert v_{k+1} - x_{\star} \Vert_2^2 &= \frac{\gamma_{k+1}}{2} \Vert \frac{\gamma_{k+1}v_k + \mu \theta_k(y-v_k) - \theta_k G_k}{\gamma_{k+1}} - x_{\star} \Vert_2^2 \\
&= \frac{\gamma_{k+1}}{2} \Vert v_k + \frac{\mu \theta_k}{\gamma_{k+1}}(y-v_k) - \theta_k G_k - x_{\star} \Vert_2^2 \\
&= \frac{\gamma_{k+1}}{2} [\Vert v_k - x_{\star} \Vert_2^2 +\Vert \theta_k G_k \Vert_2^2 - \theta_k G_k^T(v_k - x_{\star}+ \frac{\mu\theta_k}{\gamma_{k+1}}(y-v_k) )\\ 
&=\frac{\gamma_{k+1}}{2} \Vert v_k - x_{\star}+\frac{\mu \theta_k}{\gamma_{k+1}}(y-v_k) \Vert_2^2 +\frac{t_k}{2} \Vert G_k \Vert_2^2 - \theta_k G_k^T(v_k - x_{\star}+ \frac{\mu\theta_k}{\gamma_{k+1}}(y-v_k) )\\
&= \frac{\gamma_{k+1}}{2} \Vert v_k - x_{\star} \Vert_2^2 +\frac{\mu^2 \theta_k^2}{2 \gamma_{k+1}} \Vert y - v_k \Vert_2^2 +\mu \theta_k(y - v_k)^T(v_k - x_{\star}) + \frac{t_k}{2} \Vert G_k \Vert_2^2 - \theta_k G_k^T(v_k - x_{\star}+ \frac{\mu\theta_k}{\gamma_{k+1}}(y-v_k) )\\
&=\frac{\gamma_{k+1}}{2} \Vert v_k - x_{\star} \Vert_2^2 +\frac{\mu^2 \theta_k^2}{2 \gamma_{k+1}} \Vert y - v_k \Vert_2^2 + \frac{t_k}{2} \Vert G_k \Vert_2^2 + \frac{\mu \theta_k}{2}[ \Vert y -x_{\star} \Vert_2^2 - \Vert y - v_k \Vert_2^2 - \Vert v_k - x_{\star} \Vert_2^2] -  \theta_k G_k^T(v_k - x_{\star}+ \frac{\mu\theta_k}{\gamma_{k+1}}(y-v_k) )\\
&=\frac{\gamma_{k+1} - \mu \theta_k}{2} \Vert v_k - x_{\star} \Vert_2^2 +\frac{\mu \theta_k(\mu \theta_k - \gamma_{k+1})}{2 \gamma_{k+1}} \Vert y - v_k \Vert_2^2 + \frac{t_k}{2} \Vert G_k \Vert_2^2 + \frac{\mu \theta_k}{2} \Vert y - x_{\star} \Vert_2^2- \theta_k G_k^T(v_k - x_{\star}+ \frac{\mu\theta_k}{\gamma_{k+1}}(y-v_k) ) \\
& =\frac{\gamma_{k+1} - \mu \theta_k}{2} \Vert v_k - x_{\star} \Vert_2^2 + \frac{t_k}{2} \Vert G_k \Vert_2^2 + \frac{\mu \theta_k}{2} \Vert y - x_{\star} \Vert_2^2 - \frac{\mu \theta_k(1-\theta_k) \gamma_k}{2 \gamma_{k+1}}- \theta_k G_k^T(v_k - x_{\star}+ \frac{\mu\theta_k}{\gamma_{k+1}}(y-v_k) )\\
& \le\frac{\gamma_{k+1} - \mu \theta_k}{2} \Vert v_k - x_{\star} \Vert_2^2 + \frac{t_k}{2} \Vert G_k \Vert_2^2 + \frac{\mu \theta_k}{2} \Vert y - x_{\star} \Vert_2^2 - \theta_k G_k^T(v_k - x_{\star}+ \frac{\mu\theta_k}{\gamma_{k+1}}(y-v_k) )\\
&=\frac{\gamma_{k+1} - \mu \theta_k}{2} \Vert v_k - x_{\star} \Vert_2^2 + \frac{t_k}{2} \Vert G_k \Vert_2^2 + \frac{\mu \theta_k}{2} \Vert y - x_{\star} \Vert_2^2 +G_k^T(\theta_k x_{\star} - \theta_kv_k-\frac{\mu\theta_k^2}{\gamma_{k+1}}(y-v_k) )\\
\end{align}
$$



又根据更新公式可以知道，


$$
\begin{align}
\gamma_{k+1} &= (1-\theta_k) \gamma_k + \mu \theta_k \\
y &=  x_k + \frac{\theta_k \gamma_k}{\gamma_k + \mu \theta_k} (v_k - x_k) \\
x_k - y &=\frac{\theta_k \gamma_k}{\gamma_k + \mu \theta_k} (x_k - v_k)  \\
y -v_k  &= (1- \frac{\theta_k \gamma_k}{\gamma_k + \mu \theta_k})(x_k - v_k)\\
&=  \frac{(1-\theta_k)\gamma_k+ \mu \theta_k}{\gamma_k + \mu \theta_k}(x_k -v_k) \\
&= \frac{\gamma_{k+1}}{\gamma_k + \mu \theta_k}(x_k -v_k) \\
&=  \frac{\gamma_{k+1}}{\gamma_k\theta_k}(x_k -y) \\
G_k^T(\theta_k x_{\star} - \theta_kv_k-\frac{\mu\theta_k^2}{\gamma_{k+1}}(y-v_k) ) &= G_k^T(\theta_k x_{\star} - \theta_kv_k+\theta_ky - \theta_k y-\frac{\mu\theta_k^2}{\gamma_{k+1}}(y-v_k) ) \\
&=G_k^T(\theta_k x_{\star}+\theta_k (y-v_k) -\frac{\mu\theta_k^2}{\gamma_{k+1}}(y-v_k) - \theta_k y) \\
&=G_k^T(\theta_k x_{\star}+\theta_k(1-\frac{\mu\theta_k}{\gamma_{k+1}})(y-v_k) - \theta_k y) \\
&=G_k^T(\theta_k x_{\star}+\theta_k(\frac{\gamma_{k+1} -\mu\theta_k}{\gamma_{k+1}})(y-v_k) - \theta_k y) \\
&=G_k^T(\theta_k x_{\star}+\frac{\theta_k(1-\theta_k)\gamma_k}{\gamma_{k+1}}(y-v_k) - \theta_k y) \\
&=G_k^T(\theta_k x_{\star}+(1-\theta_k)(x_k - y) - \theta_k y) \\
&=G_k^T(\theta_k x_{\star}+(1-\theta_k)x_k -  y) \\
\end{align}
$$




因此可以得到另外一个关键的式子，


$$
\begin{align}
\frac{\gamma_{k+1}}{2} \Vert v_{k+1} - x_{\star} \Vert_2^2  &\le \frac{\gamma_{k+1} - \mu \theta_k}{2} \Vert v_k - x_{\star} \Vert_2^2 + \frac{t_k}{2} \Vert G_k \Vert_2^2 + \frac{\mu \theta_k}{2} \Vert y - x_{\star} \Vert_2^2 +G_k^T(\theta_k x_{\star}+(1-\theta_k)x_k -  y)\\
\end{align}
$$


继续进行化简，


$$
\begin{align}
f(x_{k+1}) &\le (1-\theta_k)f(x_k) + \theta_k f(x_{\star}) + G_k^T(y_k - (1-\theta_k) x_k - \theta_k x_{\star}) - \frac{(1-\theta_k) \mu}{2} \Vert   x_k - y_k \Vert_2^2 - \frac{\theta_k \mu}{2} \Vert x_{\star} - y_k \Vert_2^2 + (\frac{L t_k^2}{2}-t_k) \Vert G_k\Vert_2^2\\ 
&\le (1-\theta_k)f(x_k) + \theta_k f(x_{\star}) + G_k^T(y_k - (1-\theta_k) x_k - \theta_k x_{\star}) - \frac{\theta_k \mu}{2} \Vert x_{\star} - y_k \Vert_2^2 + (\frac{L t_k^2}{2}-t_k) \Vert G_k\Vert_2^2\\ 
& \le (1-\theta_k)f(x_k) + \theta_k f(x_{\star}) + G_k^T(y_k - (1-\theta_k) x_k - \theta_k x_{\star}) - \frac{\theta_k \mu}{2} \Vert x_{\star} - y_k \Vert_2^2 - \frac{t_k}{2} \Vert G_k\Vert_2^2, \text{Choosing } t_k \le \frac{1}{L}\\ 
&\le  (1-\theta_k)f(x_k) + \theta_k f(x_{\star}) - \frac{\theta_k \mu}{2} \Vert x_{\star} - y_k \Vert_2^2 - \frac{t_k}{2} \Vert G_k\Vert_2^2 +  \frac{\gamma_{k+1} - \mu \theta_k}{2} \Vert v_k - x_{\star} \Vert_2^2 + \frac{t_k}{2} \Vert G_k \Vert_2^2 + \frac{\mu \theta_k}{2} \Vert y - x_{\star} \Vert_2^2 -\frac{\gamma_{k+1}}{2} \Vert v_{k+1} - x_{\star} \Vert_2^2\\
&=(1-\theta_k)f(x_k) + \theta_k f(x_{\star})  +  \frac{\gamma_{k+1} - \mu \theta_k}{2} \Vert v_k - x_{\star} \Vert_2^2 -\frac{\gamma_{k+1}}{2} \Vert v_{k+1} - x_{\star} \Vert_2^2\\
&=(1-\theta_k)f(x_k) + \theta_k f(x_{\star})  +  \frac{(1-\theta_k) \gamma_k}{2} \Vert v_k - x_{\star} \Vert_2^2 -\frac{\gamma_{k+1}}{2} \Vert v_{k+1} - x_{\star} \Vert_2^2\\
\end{align}
$$




得到了最为关键的式子，


$$
\begin{align}
f(x_{k+1}) - f(x_{\star}) + \frac{\gamma_{k+1}}{2} \Vert v_{k+1} - x_{\star} \Vert_2^2 &\le (1-\theta_k)(f(x_k) - f(x_{\star})+ \frac{\gamma_k}{2} \Vert v_k - x_{\star} \Vert_2^2) \\
f(x_{k+1}) - f(x_{\star}) &\le \lambda_k(f(x)- f(x_{\star})+ \frac{\gamma_0}{2} \Vert x_0  - x_{\star} \Vert_2^2)  \\
&= \lambda_k(1-\theta_0)(f(x_0)- f(x_{\star})+ \frac{\gamma_0}{2} \Vert x_0  - x_{\star} \Vert_2^2) ,
\text{With } \lambda_k = \prod_{i=1}^k (1- \theta_i)
\end{align}
$$



该方法的本质上是将函数值的收敛性转化为序列$\lambda_k$的收敛性，而对于我们选定的序列$\lambda_k$, 

$$
\begin{align}
\lambda_{k+1} &= (1-\theta_k) \lambda_k = \frac{\gamma_{k+1} - \mu \theta_k}{\gamma_k} \lambda_k \le \frac{\gamma_{k+1}}{\gamma_k} \lambda_k \\
\frac{1}{\sqrt{\lambda_{k+1}}} -  \frac{1}{\sqrt{\lambda_{k}}} &= \frac{\sqrt{1-\theta_{k+1}}-1}{\sqrt{\lambda_{k+1}}}  \le \frac{\sqrt{1-\theta_{k+1}+ \theta_{k+1}^2 / 4} -1}{\sqrt{\lambda_{k+1}}}  = \frac{\theta_{k+1}}{2\sqrt{\lambda_{k+1}}}  \\
\frac{1}{\sqrt{\lambda_{k+1}}} -  \frac{1}{\sqrt{\lambda_{k}}} &\ge \frac{\theta_{k+1}}{2\sqrt{\gamma_{k+1} / \gamma_1}} = \frac{\sqrt{t_{k} \gamma_1}}{2} \\
\frac{2}{\sqrt{\lambda_{k+1}}} -  \frac{2}{\sqrt{\lambda_{1}}} &\ge  \sum_{i=1}^k \sqrt {t_i \gamma_1}\\
\lambda_{k+1} &\le \frac{4}{(2 + \sum_{i=1}^k \sqrt{t_i \gamma_1})^2}, \text{With } \lambda_1 = 1
\end{align}
$$





设定恒定步长，可以得到，


$$
\begin{align}
\lambda_{k+1} &\le \frac{4}{k+1} \\
f(x_{k+1}) - f(x_{\star}) & \le \lambda_k(1-\theta_0)(f(x_0)- f(x_{\star})+ \frac{\gamma_0}{2} \Vert x_0  - x_{\star} \Vert_2^2) \\
&=\frac{\lambda_k \gamma_0}{2} \Vert x_0 - x_{\star} \Vert_2^2 ,\text{By Choosing } \theta_0 = 1 \\
&\le \frac{2}{(k+1)^2} \Vert x_0 - x_{\star} \Vert_2^2
\end{align}
$$


这与 [FISTA](https://truenobility303.github.io/FISTA/)  的结果相同，本质上FISTA就是该算法在$\mu= 0 $的特例，代入可以得到，


$$
\begin{align}
y_k &= (1- \theta_k) x_{k} + \theta_k v_{k} \\
x_{k+1} &= \text{prox}_{th}(y_k - t_k \nabla g(y_k)) \\
v_{k+1} &= x_{k} + \frac{1}{\theta_{k+1}} (x_{k+1} - x_{k})
\end{align}
$$


正好就是 FISTA的更新公式完全一致。



但该算法的神奇之处在于可以针对$\mu \ne  0$的 情况，也即在强凸函数上面，此时可以达到线性收敛，


$$
\begin{align}
\text{If Choose } \gamma_k &\ge \mu \\
\gamma_{k+1} &= (1-\theta_k) \gamma_k + \mu \theta_k\ge (1-\theta_k) \mu + \theta_k\mu = \mu \\ 
\lambda_k &= \prod_{i=1}^k (1- \theta_i) \le (1-\sqrt{\mu t })^k
\end{align}
$$


这个收敛速度相比于原本的近端梯度下降方法是质的飞跃，可以将收敛速度进行如下对比，为了达到$\epsilon$- 最优解：

* 对于一般凸函数，近端梯度方法需要$O(\frac{1}{\epsilon})$次迭代，而加速后的算法或者FISTA算法仅需要 $O(\frac{1}{\sqrt{\epsilon}})$次迭代。

* 对于具有强凸性质的函数，近端梯度方法仍然需要$O(\frac{1}{\epsilon})$次迭代，而加速后的算法仅仅需要$O(\log \frac{1}{\epsilon})$次迭代。

