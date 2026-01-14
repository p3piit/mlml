# mlml  
*Multilevel Machine Learning utilities for clustered data*

## Overview

**mlml** is an R package developed as part of my master’s thesis on **multilevel machine learning**.  
It provides implementations of:

- **Generalized Mixed-Effects Regression Trees (GMERT)**
- **Generalized Mixed-Effects Random Forests (GMERF)**
- Optimized variants using matrix precomputations
- Simulation utilities for clustered binary data
- Prediction and evaluation tools (e.g. F1-score)

The package is designed for **clustered / hierarchical data**, where observations are grouped  
(e.g. students within schools, patients within hospitals) and both **fixed** and **random effects**  
must be modeled.

---

## Purpose of the Package (and How I Use It)

This package serves three main goals:

1. **Support my thesis analysis** with reusable, tested code  
2. **Implement methods from the literature** in a transparent way  
3. **Provide a single entry point** that loads all required tools for multilevel ML  

I use **mlml** as a *personal research toolbox* during my thesis work.

This:

- Loads all required dependencies  
- Gives access to all model-fitting functions  
- Ensures reproducibility through `renv`  

Instead of loading many packages manually, I rely on **mlml** to centralize everything. I know that usually is bad 
practice to import whole libraries and not single functions in a package, but in this way I'll get all the packages
that I need with only one line of code

---

## Implemented Methods

The package includes:

### GMERT (Tree-based)

- `fit_gmert()`  
- `predict_gmert()`  
- `fit_gmert_small()` (optimized version)

### GMERF (Random Forest-based)

- `fit_gmerf()`  
- `predict_gmerf()`  
- `fit_gmerf_small()` (optimized version)

### Utilities

- `sim_data_gmert()` – simulate clustered binary data  
- `split_gmert_data()` – cluster-aware train/test split  
- `f1_fun()` – F1 score for majority/minority class  

---

## Dependencies

The package imports:

- `tidyverse`  
- `rpart`  
- `ranger`  
- `MASS`  
- `lme4`  
- `parallel`  
- `ggplot2`  
- `performance`  

---

## License

This project is licensed under the **MIT License**.
---

## References

The implemented methods are based on the following papers:

**Hajjem, Larocque & Bellavance (2011)**  
*Mixed Effects Regression Trees for Clustered Data*  
Statistics & Probability Letters, 81(4), 451–459  
https://doi.org/10.1016/j.spl.2010.12.003  

**Hajjem, Larocque & Bellavance (2017)**  
*Generalized Mixed Effects Regression Trees*  
Statistics & Probability Letters, 126, 114–118  
https://doi.org/10.1016/j.spl.2017.02.033  

**Pellagatti et al. (2021)**  
*Generalized Mixed-Effects Random Forest*  
Statistical Analysis and Data Mining, 14(3), 241–257  
https://doi.org/10.1002/sam.11505  

---

## 👤 Author

**Paolo Colussi**  
Master’s student  
Email: p.colussi@uu.nl  
