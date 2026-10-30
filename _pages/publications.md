---
layout: archive
title: "Research"
permalink: /publications/
author_profile: true
---

<style>
  /* Publication detail disclosure */
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

  .detail-disclosure {
    display: inline;
  }

  .detail-disclosure summary {
    display: inline;
    list-style: none;
    cursor: pointer;
  }

  .detail-disclosure summary::-webkit-details-marker {
    display: none;
  }

  .detail-disclosure .hide-label {
    display: none;
  }

  .detail-disclosure[open] .expand-label {
    display: none;
  }

  .detail-disclosure[open] .hide-label {
    display: inline;
  }

  .detail-image-container {
    display: none;
    margin-top: 8px;
    margin-bottom: 12px;
  }

  .detail-disclosure[open] .detail-image-container {
    display: block;
  }

  .detail-image-container img {
    max-width: 100%;
    height: auto;
    border: 1px solid #ddd;
    border-radius: 6px;
    box-shadow: 0 2px 5px rgba(0,0,0,0.1);
  }

  .highlight-box {
    background-color: #f6f8fa;
    border: 1px solid #d1d5da;
    border-radius: 6px;
    padding: 12px 16px;
    margin-top: 10px;
    margin-bottom: 10px;
    font-family: ui-monospace, SFMono-Regular, SF Mono, Menlo, Consolas, Liberation Mono, monospace;
    font-size: 0.9em;
    color: #24292e;
    line-height: 1.5;
  }

  /* Tab wrapper */
  .pub-tabs-wrapper {
    position: relative;
  }

  /* Visually hidden radio inputs */
  .pub-tab-radio {
    position: absolute;
    opacity: 0;
    width: 1px;
    height: 1px;
    overflow: hidden;
    pointer-events: none;
  }

  /* Tab buttons */
  .pub-tabs {
    display: flex;
    gap: 8px;
    border-bottom: 1px solid #d1d5da;
    margin-bottom: 20px;
    margin-top: 15px;
  }

  .pub-tab-btn {
    display: inline-block;
    background-color: #f6f8fa;
    border: 1px solid #d1d5da;
    border-bottom: none;
    padding: 8px 18px;
    font-size: 0.85em;
    font-weight: 600;
    letter-spacing: 0.5px;
    color: #586069;
    cursor: pointer;
    border-top-left-radius: 6px;
    border-top-right-radius: 6px;
    transition: all 0.2s ease;
    margin-bottom: -1px;
  }

  .pub-tab-btn:hover {
    color: #24292e;
    background-color: #e1e4e8;
  }

  /* Active tab styling */
  #pub-tab-highlighted-radio:checked ~ .pub-tabs label[for="pub-tab-highlighted-radio"],
  #pub-tab-all-radio:checked ~ .pub-tabs label[for="pub-tab-all-radio"] {
    background-color: #ffffff;
    color: #24292e;
    border-bottom: 2px solid #ffffff;
  }

  /* Keyboard focus */
  #pub-tab-highlighted-radio:focus-visible ~ .pub-tabs label[for="pub-tab-highlighted-radio"],
  #pub-tab-all-radio:focus-visible ~ .pub-tabs label[for="pub-tab-all-radio"] {
    outline: 2px solid #0366d6;
    outline-offset: 2px;
  }

  /* Hide all tab contents by default */
  .pub-tab-content {
    display: none;
  }

  /* Show selected tab content */
  #pub-tab-highlighted-radio:checked ~ #tab-highlighted {
    display: block;
  }

  #pub-tab-all-radio:checked ~ #tab-all {
    display: block;
  }

</style>


<!-- ===================================================== -->
<!-- TAB WRAPPER -->
<!-- ===================================================== -->

<div class="pub-tabs-wrapper">
  <!-- Radio controls -->
  <input
    type="radio"
    class="pub-tab-radio"
    name="publication-tabs"
    id="pub-tab-highlighted-radio"
    checked>

  <input
    type="radio"
    class="pub-tab-radio"
    name="publication-tabs"
    id="pub-tab-all-radio">

  <!-- Tab buttons -->
  <div class="pub-tabs">

    <label
      class="pub-tab-btn"
      for="pub-tab-highlighted-radio">
      HIGHLIGHTED
    </label>

    <label
      class="pub-tab-btn"
      for="pub-tab-all-radio">
      ALL
    </label>

  </div>


  <!-- ===================================================== -->
  <!-- HIGHLIGHTED TAB -->
  <!-- ===================================================== -->

  <div id="tab-highlighted" class="pub-tab-content">

    <h2>Theory of Game-Structured Optimization</h2>

    My series of work provides (near-)optimal algorithms and matching lower bounds on target precision $\epsilon$ for solving game-structured optimization such as minimax and bilevel optimization.

    <ol class="custom-ol">

      <li>
        Second-Order Minimax Optimization
        <a href="https://drive.google.com/file/d/18vEccWx-tONoFAaHDexAVlu-QW2tqRO_/view?usp=sharing">[Slides in SIAM OP 26]</a>

        <details class="detail-disclosure">
          <summary class="detail-toggle">
            <span class="expand-label">[Expand]</span>
            <span class="hide-label">[Hide]</span>
          </summary>

          <div class="detail-image-container">
            <img src="/images/research/Minimax.jpg" alt="Second-Order Minimax Optimization Detail">
          </div>
        </details>

        <br>

        Upper Bound:
        <a href="http://arxiv.org/abs/2506.08362">[COLT 2025]</a>
        <a href="https://arxiv.org/abs/2609.23557">[arXiv 2026]</a>

        <br>

        Lower Bound:
        <a href="https://arxiv.org/pdf/2604.19462">[arXiv 2026]</a>
      </li>


      <li>
        Bilevel Optimization
        <a href="https://drive.google.com/file/d/1lghkhCXdDP9IQWwBK_E_sfCjokgFlrGu/view?usp=sharing">[Slides in SIAM OP 26]</a>

        <details class="detail-disclosure">
          <summary class="detail-toggle">
            <span class="expand-label">[Expand]</span>
            <span class="hide-label">[Hide]</span>
          </summary>

          <div class="detail-image-container">
            <img src="/images/research/Bilevel.jpg" alt="Bilevel Optimization Detail">
          </div>
        </details>

        <br>

        Upper Bound:
        <a href="https://arxiv.org/abs/2306.14853">[JMLR 2025]</a>

        <br>

        Lower Bound:
        <a href="https://arxiv.org/abs/2511.22331">[arXiv 2025]</a>
      </li>

    </ol>


    <h2>AI for Optimization</h2>

    I've been exploring these questions for more than one year, and I'm excited to see them fully resolved under the assistance of GPT.

    <ol class="custom-ol">

      <font size="3">

        <li>
          <p>
            Qihao Zhou <sup>1</sup>, <b>Lesi Chen <sup>1</sup></b>, Luo Luo, Chengchang Liu,
            <i>Envelope Cubic Regularized Newton Method for Nonconvex-Strongly-Concave Minimax Optimization</i>,
            arXiv preprint.
            <a href="https://arxiv.org/abs/2610.20345">[arXiv 2026]</a>
          </p>
        </li>

        <li>
          <p>
            Xinliang Zhang <sup>1</sup>, <b>Lesi Chen<sup>1</sup></b>, Chengchang Liu, Jingzhao Zhang,
            <i>Near-Optimal Convex Optimization with Lazy Second-Order Oracles</i>,
            arXiv preprint.
            <a href="https://arxiv.org/abs/2610.03222">[arXiv 2026]</a>
          </p>
        </li>

        <li>
          <p>
            Xinliang Zhang <sup>1</sup>, <b>Lesi Chen<sup>1</sup></b>, Linxuan Pan<sup>1</sup>, Chengchang Liu, Junchi Yang, Jingzhao Zhang,
            <i>Optimal High-Order Methods for Solving Monotone Variational Inequalities</i>,
            arXiv preprint.
            <a href="https://arxiv.org/abs/2609.23557">[arXiv 2026]</a>
          </p>
        </li>

      </font>

    </ol>

  </div>


  <!-- ===================================================== -->
  <!-- ALL TAB -->
  <!-- ===================================================== -->

<div id="tab-all" class="pub-tab-content">

<h3>Featured Publications</h3>

<p>
See
<a href="https://scholar.google.com/citations?user=ynGzhugAAAAJ&hl=en&oi=ao" target="_blank" rel="noopener noreferrer">Google Scholar</a>
for a complete list.
</p>

<ol class="custom-ol">

<font size="3">

<!-- Working Paper 1 -->
<li id="extension-paper">
<p>
<b>Lesi Chen<sup>1</sup></b>, Xinliang Zhang <sup>1</sup>, Junru Li, Chengchang Liu, Luo Luo, and Jingzhao Zhang,
<i>Solving Convex-Concave Problems with $\tilde{\mathcal{O}}(\epsilon^{-4/(3p+1)})$ pth-Order Oracle Complexity</i>,
arXiv preprint.
<a href="https://arxiv.org/pdf/2604.19462">[arXiv 2026]</a>
<br>
[Extension of <a href="#colt2025-paper">COLT 2025 below</a>]
</p>
</li>


<li>
<p>
<b>Lesi Chen</b>, Chengchang Liu, Luo Luo, John C.S. Lui, and Jingzhao Zhang,
<i>Optimal Convex Optimization with Inexact Second-Order Oracles</i>,
arXiv preprint.
<a href="https://arxiv.org/abs/2607.24520">[arXiv 2026]</a>
</p>
</li>


<li>
<p>
<b>Lesi Chen<sup>1</sup></b>, Kaiyi Ji <sup>1</sup>, and Jingzhao Zhang,
<i>On the Condition Number Dependency in Bilevel Optimization</i>,
arXiv preprint.
<a href="https://arxiv.org/abs/2511.22331">[arXiv 2026]</a> ⭐
</p>
</li>


<li>
<p>
<b>Lesi Chen</b>, Chengchang Liu, Luo Luo, and Jingzhao Zhang,
<i>Faster Newton Methods for Convex and Nonconvex Optimization in Gradient Complexity</i>,
in Conference on Learning Theory.
<a href="https://arxiv.org/abs/2501.17488">[COLT 2026]</a> ⭐
</p>
</li>


<li>
<p>
<b>Lesi Chen</b>, Junru Li, El Mahdi Chayti and Jingzhao Zhang,
<i>Faster Gradient Methods for Highly-Smooth Stochastic Bilevel Optimization</i>,
in International Conference on Learning Representations.
<a href="https://arxiv.org/abs/2509.02937">[ICLR 2026]</a>
</p>
</li>


<li>
<p>
<b>Lesi Chen<sup>1</sup></b>, Yaohua Ma<sup>1</sup>, and Jingzhao Zhang,
<i>Near-Optimal Nonconvex-Strongly-Convex Bilevel Optimization with Fully First-Order Oracles</i>,
Journal of Machine Learning Research, 1-56.
<a href="https://arxiv.org/abs/2306.14853">[JMLR 2025]</a> ⭐
<br>
</p>
</li>


<li id="colt2025-paper">
<p>
<b>Lesi Chen</b>, Chengchang Liu, Luo Luo, and Jingzhao Zhang,
<i>Solving Convex-Concave Problems with $\tilde{\mathcal{O}}(\epsilon^{-4/7})$ pth-Order Oracle Complexity</i>,
in Conference on Learning Theory.
<b>(Best Student Paper, 2/556)</b>
<a href="http://arxiv.org/abs/2506.08362">[COLT 2025]</a>
[Extension <a href="#extension-paper">available above</a>] ⭐
</p>
</li>


<li>
<p>
<b>Lesi Chen<sup>1</sup></b>, Chengchang Liu<sup>1</sup>, and Jingzhao Zhang,
<i>Second-Order Min-Max Optimization with Lazy Hessians</i>,
in International Conference on Learning Representations.
<b>(Oral, &lt;2%)</b>
<a href="https://arxiv.org/abs/2410.09568">[ICLR 2025]</a>
</p>
</li>


<li>
<p>
<b>Lesi Chen</b> and Luo Luo,
<i>Near-Optimal Algorithms for Making the Gradient Small in Stochastic Minimax Optimization</i>,
Journal of Machine Learning Research, 1-44.
<a href="https://arxiv.org/abs/2208.05925">[JMLR 2024]</a>
</p>
</li>


<li>
<p>
Huaqing Zhang<sup>1</sup>, <b>Lesi Chen<sup>1</sup></b>, Jing Xu, and Jingzhao Zhang,
<i>Functionally Constrained Algorithm Solves Convex Simple Bilevel Problems</i>,
in Conference on Neural Information Processing Systems.
<br>
<a href="https://arxiv.org/abs/2409.06530">[NeurIPS 2024]</a>
</p>
</li>


<li>
<p>
<b>Lesi Chen<sup>1</sup></b>, Jing Xu<sup>1</sup>, and Jingzhao Zhang,
<i>On Finding Small Hyper-Gradients in Bilevel Optimization: Hardness Results and Improved Analysis</i>,
in Conference on Learning Theory.
<a href="https://arxiv.org/abs/2301.00712">[COLT 2024]</a>
</p>
</li>


<li>
<p>
<b>Lesi Chen</b>, Jing Xu, and Luo Luo,
<i>Faster Gradient-Free Algorithms for Nonsmooth Nonconvex Stochastic Optimization</i>,
in International Conference on Machine Learning.
<a href="https://arxiv.org/abs/2301.06428">[ICML 2023]</a>
</p>
</li>


<li>
<p>
<b>Lesi Chen</b>, Boyuan Yao, and Luo Luo,
<i>Faster Stochastic Algorithms for Minimax Optimization under Polyak-Łojasiewicz Condition</i>,
in Conference on Neural Information Processing Systems.
<a href="https://arxiv.org/abs/2307.15868">[NeurIPS 2022]</a>
</p>
</li>

</font>
</ol>

</div>
