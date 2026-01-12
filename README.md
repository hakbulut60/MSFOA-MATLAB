# MSFOA – Modified Starfish Optimization Algorithm (MATLAB)

This repository provides a MATLAB implementation of the **Modified Starfish Optimization Algorithm (MSFOA)**, an enhanced variant of the original Starfish Optimization Algorithm (SFOA).

## 🔧 Key Features
- Exponentially decreasing adaptive step size
- Momentum vector integration
- Multi-candidate local search mechanism
- Dynamic exploration–exploitation balance

## 📌 Reference (Mandatory Citation)

If you use this code in your research, please cite the following paper:

> H. Akbulut,  
> *A modified starfish optimization algorithm (M-SFOA) for global optimization problems and its application to heart disease risk prediction*,  
> **Expert Systems with Applications**, Volume 307, 2026, Article 131088.  
> DOI: 10.1016/j.eswa.2026.131088

BibTeX is available in `docs/citation.bib`.

## ▶ Usage Example

```matlab
Npop = 30;
Max_it = 500;
nD = 30;
lb = -100;
ub = 100;

[xbest, fbest, curve] = MSFOA(Npop, Max_it, lb, ub, nD, @sphere);
