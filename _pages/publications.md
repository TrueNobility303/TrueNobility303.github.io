---
layout: archive
title: "Research"
permalink: /publications/
author_profile: true
---

<style>
  .detail-toggle {
    color: #0366d6;
    cursor: pointer;
    text-decoration: underline;
    font-size: 0.9em;
    margin-left: 6px;
  }
  .detail-toggle:hover {
    color: #0056b3;
  }
  .detail-image-container {
    display: none;
    margin-top: 8px;
    margin-bottom: 12px;
  }
  .detail-image-container img {
    max-width: 100%;
    height: auto;
    border: 1px solid #ddd;
    border-radius: 6px;
    box-shadow: 0 2px 5px rgba(0,0,0,0.1);
  }
</style>

<script>
  function toggleDetail(id) {
    var container = document.getElementById(id);
    var btn = document.getElementById('btn-' + id);
    if (container.style.display === "block") {
      container.style.display = "none";
      btn.innerText = "[Expand]";
    } else {
      container.style.display = "block";
      btn.innerText = "[Hide]";
    }
  }
</script>

I have broad interests in modern optimization theory, with a particular focus on game-structure optimization. In my representative work ⭐, I have collaborated with my excellent coauthors to solve many fundamental problems in this area, such as the oracle complexity of convex, nonconvex, minimax, and bilevel optimization.

<h2> Research Highlights </h2>
<ol class="custom-ol">
<li> Second-Order Minimax Optimization <a href="https://drive.google.com/file/d/18vEccWx-tONoFAaHDexAVlu-QW2tqRO_/view?usp=sharing">[Slides in SIAM OP 26]</a> 
 Upper Bound: $\tilde{\mathcal{O}}(1/T^{1.75})$ <a href="http://arxiv.org/abs/2506.08362">[COLT 2025]</a> $\to$ $\tilde{\mathcal{O}}(1/T^{2})$  <a href="https://arxiv.org/abs/2608.08463">[arXiv 2026]</a> <br>
 Lower Bound: $\Omega(1/T^{2.5})$ <a href="https://arxiv.org/pdf/2604.19462">[arXiv 2026]</a>
<span class="detail-toggle" id="btn-detail-1" onclick="toggleDetail('detail-1')">[Expand]</span> <br>
 <div id="detail-1" class="detail-image-container">
   <img src="/images/research/Minimax.png" alt="Second-Order Minimax Optimization Detail">
 </div>
</li>

<li> Bilevel Optimization <a href="https://drive.google.com/file/d/1lghkhCXdDP9IQWwBK_E_sfCjokgFlrGu/view?usp=sharing">[Slides in SIAM OP 26]</a> <span class="detail-toggle" id="btn-detail-2" onclick="toggleDetail('detail-2')">[Expand]</span> <br>
 Upper Bound: $\tilde{\mathcal{O}}(\kappa^{3.5} \epsilon^{-2} + \sigma^2 \kappa^{11} \epsilon^{-6})$ <a href="https://arxiv.org/abs/2306.14853">[JMLR 2025]</a> <br>
 Lower Bound: $\Omega(\kappa^{2.5} \epsilon^{-2} + \sigma^2 \kappa^{4.5} \epsilon^{-4})$ <a href="https://arxiv.org/abs/2511.22331">[arXiv 2025]</a> 
 
 <div id="detail-2" class="detail-image-container">
   <img src="/images/research/Bilevel.png" alt="Bilevel Optimization Detail">
 </div>
</li>
</ol>

<!-- A list of open problems is maintained <a href="https://truenobility303.github.io/openproblems/">here</a>. -->

<h2> Working Papers </h2>
<ol class="custom-ol">
<font size="3">  
<li><p> <b>Lesi Chen<sup>1</sup></b>, Xinliang Zhang <sup>1</sup>, Junru Li, Chengchang Liu, Luo Luo, and Jingzhao Zhang,  <i> Solving Convex-Concave Problems with $\tilde{\mathcal{O}}(\epsilon^{-4/(3p+1)})$ pth-Order Oracle Complexity</i>,  arXiv preprint. <a href="https://arxiv.org/pdf/2604.19462">[arXiv 2026]</a> <br>
Preliminary version in Conference on Learning Theory. <b>(Best Student Paper, 2/556) </b> <a href="http://arxiv.org/abs/2506.08362">[COLT 2025]</a> ⭐ 

</p></li> 
<li><p> 
<b>Lesi Chen<sup>1</sup></b>, Xinliang Zhang <sup>1</sup>, Hengyu Wang <sup>1</sup>, Chengchang Liu, Yongchao Chen, and Jingzhao Zhang, <i>Halpern Iteration Achieves $\tilde{\mathcal{O}}(\epsilon^{-1/p})$ pth-Order Oracle Complexity for Monotone Variational Inequalities</i>, arXiv preprint. <a href="https://arxiv.org/abs/2608.08463">[arXiv 2026]</a> ⭐
</p></li>
 <li><p>
<b>Lesi Chen<sup>1</sup> </b>, Kaiyi Ji <sup>1</sup>, and Jingzhao Zhang, <i>On the Condition Number Dependency in Bilevel Optimization</i>, arXiv preprint. <a href="https://arxiv.org/abs/2511.22331">[arXiv 2025]</a> ⭐
</p></li>  
 <li><p>  <b>Lesi Chen</b>, Chengchang Liu, Luo Luo, John C.S. Lui, and Jingzhao Zhang,
<i> Optimal Convex Optimization with Inexact Second-Order Oracles </i>, arXiv preprint. <a href="https://arxiv.org/abs/2607.24520">[arXiv 2026]</a>
</p></li>
</font>
</ol>

 
<h2> Featured Publications </h2>

<!-- <sup>1</sup> indicates co-first-authors. -->

The conference and journal publications that have overlapped are grouped into one item.
<ol class="custom-ol">
<font size="3">      
 <li><p>  <b>Lesi Chen</b>, Chengchang Liu, Luo Luo, and Jingzhao Zhang,
<i> Faster Newton Methods for Convex and Nonconvex Optimization in Gradient Complexity </i>, in Conference on Learning Theory.
 <a href="https://arxiv.org/abs/2501.17488">[COLT 2026]</a> ⭐
</p></li>
<li><p> 
<b>Lesi Chen</b>, Junru Li, El Mahdi Chayti and Jingzhao Zhang, <i> Faster Gradient Methods for Highly-Smooth Stochastic Bilevel Optimization</i>, in International Conference on Learning Representations. <a href="https://arxiv.org/abs/2509.02937">[ICLR 2026]</a> 
</p>
</li>
<li><p> <b>Lesi Chen<sup>1</sup></b>, Yaohua Ma<sup>1</sup>, and Jingzhao Zhang,
 <i> Near-Optimal Nonconvex-Strongly-Convex Bilevel Optimization with Fully First-Order Oracles </i>, Journal of Machine Learning Research, 1-56. 
 <a href="https://arxiv.org/abs/2306.14853">[JMLR 2025]</a> ⭐ <br>
 </p></li>
<li><p> <b>Lesi Chen<sup>1</sup></b>, Chengchang Liu<sup>1</sup>, and Jingzhao Zhang,  <i> Second-Order Min-Max Optimization with Lazy Hessians</i>, in International Conference on  Learning Representations. <b>(Oral, <2%) </b>   <a href="https://arxiv.org/pdf/2410.09568">[ICLR 2025]</a> 
</p></li>
<li><p> <b>Lesi Chen</b> and Luo Luo, <i> Near-Optimal Algorithms for Making the Gradient Small in Stochastic Minimax Optimization</i>, Journal of Machine Learning Research, 1-44. 
 <a href="https://arxiv.org/abs/2208.05925">[JMLR 2024]</a>
</p></li> 
<li><p> Huaqing Zhang<sup>1</sup>, <b>Lesi Chen<sup>1</sup></b>, Jing Xu, and Jingzhao Zhang, <i>
 Functionally Constrained Algorithm Solves Convex Simple Bilevel Problems</i>, in Conference on Neural Information Processing Systems. <br>
 <a href="https://arxiv.org/abs/2409.06530">[NeurIPS 2024]</a>
 </p></li>
 <li> <p>
<b>Lesi Chen<sup>1</sup></b>, Jing Xu<sup>1</sup>, and Jingzhao Zhang, <i> On Finding Small Hyper-Gradients in Bilevel Optimization: Hardness Results and Improved Analysis</i>, in Conference on Learning Theory. 
<a href="https://arxiv.org/abs/2301.00712">[COLT 2024]</a> 
  </p> </li>
<li><p> <b>Lesi Chen</b>, Jing Xu, and Luo Luo, <i> Faster Gradient-Free Algorithms for Nonsmooth Nonconvex Stochastic Optimization</i>,
 in International Conference on Machine Learning. 
  <a href="https://arxiv.org/abs/2301.06428"> [ICML 2023]</a>
 </p> </li>
<li><p>  <b>Lesi Chen</b>, Boyuan Yao, and Luo Luo, <i> Faster Stochastic Algorithms for Minimax Optimization under Polyak-Łojasiewicz Condition</i>, 
 in Conference on Neural Information Processing Systems.
  <a href="https://arxiv.org/abs/2307.15868"> [NeurIPS 2022]</a>
 </p> </li>
</font>
</ol>

See [Google Scholar](https://scholar.google.com/citations?user=ynGzhugAAAAJ&hl=en&oi=ao) for a complete list.
