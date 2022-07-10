# CMS
Performs physics-based model reduction using component mode synthesis (CMS) based on the Craig-Bampton method. 
Works with 3D and 2D structural FE models created with COMSOL Multiphysics. The formulations are based on ["Sub-structure Coupling for Dynamic Analysis"](https://link.springer.com/book/10.1007/978-3-030-12819-7) by H.Jensen & C.Papadimitriou.

The code is written in MATLAB and communicates with COMSOL Multiphysics through LiveLink for MATLAB. It was developed for my diploma thesis ["Advanced Model Reduction Techniques for Structural Dynamics Simulations"](https://ir.lib.uth.gr/xmlui/handle/11615/57557?locale-attribute=en) at the Department of Mechanical Engineering of the University of Thessaly. MATLAB must be started with a COMSOL server (using LiveLink).
In the future, a user guide might be written.

# Main Features
Non-parametrized CMS: The non-parametrized reduced-order matrices are constant.
Three variants:
1) Without interface reduction
2) With global-level interface reduction
3) With local-level interface reduction

Parametrized CMS: The parametrized reduced-order matrices depend on model parameters. Model parametrization is used in structural dynamics simulations.
Three variants:
1) Without interface reduction
2) With global-level interface reduction
3) With local-level interface reduction

Additional Features: Consideration of residual normal modes (static correction), meta-model for interpolation of interface modes using minimum number of support points, definition of interfaces by user, optimization of number of kept modes, efficient RAM usage, highly parallelized

# Some results

## Non-parametrized CMS
Non-parametrized CMS was applied to a high-fidelity FE model of a highway bridge consisting of 944,613 DOFs. The model is divided in 22 substructures and results in 928,200 internall DOFs (not shared between two or more components) and 16,413 interface DOFs (existing at the interface of two or more components).

### Division in substructures (in COMSOL)
![](https://github.com/FK-MAD/CMS/blob/main/Metsovo%20bridge%20results/metsovo%2022%20parameters%20iso%20-%20numbered.png?raw=true)


### Reduction of internal DOF per component
![](https://github.com/FK-MAD/CMS/blob/main/Metsovo%20bridge%20results/internall%20dofs%20full%20vs%20reduced.svg?raw=true)


### Fractional modal frequency error between full and reduced-order models
![](https://github.com/FK-MAD/CMS/blob/main/Metsovo%20bridge%20results/no%20vs%20global%20vs%20local.svg?raw=true)

## Parametrized CMS
Parametrized CMS was applied on the same bridge model. Each sub-structure is associated with one parameter (22 parameters) affecting its modulus of elasticity. Three parametrized reduced-order models (pROMs) were developed:
- NIR-pROM: reduced-order model where reduction occurs only on internal DOFs of each component (No Interface Reduction). It consists of 46 internal DOFs and 16,413 interface DOFs (16,459 DOFs in total).
- GIR-pROM/SX1: reduced-order model where reduction occurs both on internal DOFs of each component and on the DOFs of the interface at the global level (Global Interface Reduction). It consists of 46 internal DOFs and 36 interface DOFs (82 DOFs in total).
- LIR-pROM/C: reduced-order model where reduction occurs both on internal DOFs of each component and on the DOFs of the interface at the local level (Local Interface Reduction). It consists of 46 internal DOFs and 291 interface DOFs (337 DOFs in total).

In this application, the pROMs were used to predict the first few modal frequencies of the bridge for different values of the 22 parameters affecting the modulus of elasticity of each sub-structure. To test the speed and accuracy of each pROM compared to the full model, 100 runs (predictions) were performed with parameters sampled from a 22-dim Gaussian distribution at every run.

### Accuracy of pROMs

![](https://github.com/FK-MAD/CMS/blob/main/Metsovo%20bridge%20results/errors_all_v3.svg?raw=true)

*Median fractional modal frequency error (from 100 runs) - as a function of eigenmode number - between the predictions of the full model and the three pROMs*

### Speed of pROMs

| Model  | Full FE Model | NIR-pROM | GIR-pROM/SX1 | LIR-pROM/C |
| ------------- | :-------------: | :-------------: | :-------------: | :-------------: |
| Mean total iteration time [sec] | 116.5 | 32.6 | 6.9 | 0.02  |
| Speedup over the full FE model | 1x | 3.6x | 16.9x | 5825x  |

*Mean computational times (from 100 runs) of a single simulation step and speedup for the full FE model and the three pROMs*

# License
This work is licensed under a
[Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International Public License][cc-by-nc-nd].

[![CC BY-NC-ND 4.0][cc-by-nc-nd-image]][cc-by-nc-nd]

[cc-by-nc-nd]: http://creativecommons.org/licenses/by-nc-nd/4.0/
[cc-by-nc-nd-image]: https://licensebuttons.net/l/by-nc-nd/4.0/88x31.png
[cc-by-nc-nd-shield]: https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg
