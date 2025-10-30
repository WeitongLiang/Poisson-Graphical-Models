# Poisson Graphical Models and Their Extensions

This repository contains an independent replication and evaluation of **Poisson Graphical Models (PGM)** and their two extensions: the **Quadratic PGM (QPGM)** and **Sublinear PGM (SPGM)**, originally proposed in Yang et al. (2013). These models are used to estimate sparse graphical structures in count-valued data, with applications in biological and biomedical domains.

The project includes:
- Custom implementation of simulation and model fitting routines in R
- Structural recovery evaluation using precision, recall, F1 score, and AUC
- Application to breast cancer miRNA expression data from TCGA

**Full write-up**: [Poisson Graphical Models & Their Extensions.pdf](./Poisson%20Graphical%20Models%20&%20Their%20Extentions.pdf)

---

## Project Overview

Graphical models are powerful tools for learning the dependency structure among multivariate variables. While Gaussian Graphical Models are widely used, they are not well-suited for count data. This motivates the use of **Poisson Graphical Models**, which directly model counts via exponential family distributions.

However, standard PGMs impose non-positivity constraints on edge weights. The QPGM and SPGM address this by:
- QPGM: Introducing a **quadratic base measure**, enabling both positive and negative interactions.
- SPGM: Using a **sublinear transformation** on the sufficient statistics to allow more flexible dependency modeling and heavier tails.

---

## Simulation Studies

We conducted large-scale simulations to evaluate model recovery:

- **PGM**: Replicated success-rate curves under 2D lattice graphs, using ℓ₁-penalized Poisson regression (`glmnet`).  
- **QPGM & SPGM**: Simulated data from scale-free graphs using positive edge weights. Recovery metrics were computed after fitting via `optim()` with ℓ₁ penalty.

Results show that:
- PGM recovers structure well under large sample sizes.
- SPGM achieves higher accuracy than QPGM under positive dependencies.
- Our models required larger samples than originally reported to reach similar performance.

---

## Real Data Application

We applied QPGM and SPGM to breast cancer miRNA expression data from **TCGA** (n = 445, p = 353). Preprocessing followed the pipeline in the `XMRF` package. Model fitting was done via node-wise likelihood maximization using custom optimizers. The resulting networks are sparse and biologically interpretable.

---

## ⚙️ Code and Implementation

All code is written in base R and includes:

- Simulation scripts for PGM, QPGM, SPGM
- Penalized likelihood estimation routines (using `glmnet` or `optim`)
- Graph recovery metrics and plotting utilities

We avoided relying on the `XMRF` package due to its numerical instability and reproducibility issues, and instead re-implemented all methods from scratch.

---

---

## 📚 References

1. Baxter, M. (1990). Generalised linear models, by P. McCullagh and J.A. Nelder. *The Mathematical Gazette*, 74(469), 320–321.  
2. Besag, J. (1974). Spatial interaction and the statistical analysis of lattice systems. *J. R. Stat. Soc. Series B (Methodological)*, 36(2), 192–225.  
3. Griffith, D. A. (2002). A spatial filtering specification for the auto-poisson model. *Statistics & Probability Letters*, 58(3), 245–251.  
4. Jalali, A., Ravikumar, P., Vasuki, V., & Sanghavi, S. (2011). On learning discrete graphical models using group-sparse regularization. In *AISTATS* (pp. 378–387).  
5. Kaiser, M. S., & Cressie, N. (1997). Modeling Poisson variables with positive spatial dependence. *Statistics & Probability Letters*, 35(4), 423–432.  
6. Krishnamoorthy, A. (1951). Multivariate binomial and Poisson distributions. *Sankhyā: The Indian Journal of Statistics*, 117–124.  
7. Lauritzen, S. L. (1996). *Graphical Models* (Vol. 17). Oxford University Press.  
8. Meinshausen, N., & Bühlmann, P. (2006). High-dimensional graphs and variable selection with the lasso.  
9. Tay, J. K., Narasimhan, B., & Hastie, T. (2023). Elastic net regularization paths for all generalized linear models. *Journal of Statistical Software*, 106, 1–31.  
10. Wainwright, M. J., Lafferty, J., & Ravikumar, P. (2006). High-dimensional graphical model selection using ℓ₁-regularized logistic regression. *NeurIPS*, 19.  
11. Yahav, I., & Shmueli, G. (2007). An elegant method for generating multivariate Poisson random variables. arXiv:0710.5670.  
12. Yang, E., Allen, G. I., Liu, Z., & Ravikumar, P. (2012). Graphical models via generalized linear models. *NeurIPS*, 25.  

---

## 👥 Authors

- Weitong Liang  
- Xueyan Hu
