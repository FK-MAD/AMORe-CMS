# CMS
Performs physics-based model reduction using component mode synthesis (CMS) based on the Craig-Bampton method. 
Works with 3D and 2D structural FE models created with COMSOL Multiphysics. The formulations are based on ["Sub-structure Coupling for Dynamic Analysis"](https://link.springer.com/book/10.1007/978-3-030-12819-7) by H.Jensen & C.Papadimitriou.

The code is written in MATLAB and communicates with COMSOL Multiphysics through LiveLink for MATLAB. I developed it for my diploma thesis ["Advanced Model Reduction Techniques for Structural Dynamics Simulations"](https://ir.lib.uth.gr/xmlui/handle/11615/57557?locale-attribute=en) at the Department of Mechanical Engineering of the University of Thessaly. MATLAB must be started with a COMSOL server (using LiveLink).
In the future, a user guide might be written.

# Main Features
## Non-parametrized CMS
The non-parametrized reduced-order matrices are constant, independent of model parameters. Three variants of reduced-order models (ROMs) can be created based on treatment of the degrees of freedom (DOFs) at the interface between two or more components.

- NIR-ROM: ROM where reduction occurs only on internal DOFs of each component (No Interface Reduction)
- GIR-ROM: ROM where reduction occurs both on internal DOFs of each component and on the DOFs of the interface at the global level (Global Interface Reduction)
- LIR-ROM: ROM where reduction occurs both on internal DOFs of each component and on the DOFs of the interface at the local level (Local Interface Reduction)

## Parametrized CMS
The parametrized reduced-order matrices depend on model parameters. Model parametrization is used in structural dynamics simulations. Again, three variants of parametrized ROMs (pROMs) can be created based on treatment of interface DOFs.

- NIR-pROM: pROM where reduction occurs only on internal DOFs of each component (No Interface Reduction)
- GIR-pROM: pROM where reduction occurs both on internal DOFs of each component and on the DOFs of the interface at the global level (Global Interface Reduction)
- LIR-pROM: pROM where reduction occurs both on internal DOFs of each component and on the DOFs of the interface at the local level (Local Interface Reduction)

## Additional Features
Some additional features that the user can set through the input file are:
- Consideration of residual normal modes (static correction)
- Meta-model for interpolation of interface modes (using minimum number of support points) without explicitly solving the interface eigenproblem
- Automatic definition of interfaces
- Optimization of number of kept modes to achieve a set accuracy
- Ability to run in parallel on systems with multi-core CPUs

# Case Study from my Diploma Thesis
To test the computational efficiency and accuracy of ROMs and pROMs created using this technique, a high-fidelity FE model of a highway bridge consisting of 944,613 DOFs was used.

## Non-parametrized CMS
Non-parametrized CMS was applied to the bridge model. The model is divided in 22 substructures and results in 928,200 internall DOFs (not shared between two or more components) and 16,413 interface DOFs (existing at the interface of two or more components).

Three ROMs were created each one with different treatment of interface modes. The complexity of ROMs in term of number of kept DOFs can be seen below.
| Model  | Full FE Model | NIR-ROM | GIR-ROM | LIR-ROM |
| ------------- | :-------------: | :-------------: | :-------------: | :-------------: |
| Internal DOF | 928,200	| 46 | 46	| 46 |
| Interface DOF | 16,413 | 16,413 | 36 | 291 |
| Total DOF | 944,613 | 16,459 | 82 |	337 |

### Division in substructures (in COMSOL)
![](https://github.com/FK-MAD/CMS/blob/main/Metsovo%20bridge%20results/metsovo%2022%20parameters%20iso%20-%20numbered.png?raw=true)

### Number of DOF per component of the full-order model and NIR-ROM
![](https://github.com/FK-MAD/CMS/blob/main/Metsovo%20bridge%20results/internall%20dofs%20full%20vs%20reduced.svg?raw=true)

### Fractional modal frequency error - as a function of eigenmode number - between the predictions of the full-order model and the three ROMs
![](https://github.com/FK-MAD/CMS/blob/main/Metsovo%20bridge%20results/no%20vs%20global%20vs%20local.svg?raw=true)

It can be seen that the maximum error for all ROMs is **$\sim10^{-2}$ or approximately 1%**.


## Parametrized CMS
Parametrized CMS was applied on the same bridge model. Each sub-structure is associated with one parameter (22 parameters) affecting its modulus of elasticity.

In this application, the pROMs were used to predict the first twenty modal frequencies of the bridge for different values of the 22 parameters affecting the modulus of elasticity of each sub-structure. To test the speed and accuracy of each pROM compared to the full model, 100 runs (predictions) were performed with parameters sampled from a 22-dim Gaussian distribution at every run.

Three pROMs were created with same number of internal and interface DOFs as the ROMs presented above. Apart from the difference in number of kept DOFs, there are also differences in the way interface modes are calculated at each step of the simulation for each pROM.
- NIR-pROM: No interface reduction is performed, interface modes are explicitly calculated at each step of the simulation
- GIR-pROM/SX1: At each step of the simulation, interface modes are interpolated (not explicitly solved) using a proposed meta-model
- LIR-pROM/C: Interface modes are calculated once for the reference model and kept constant at each step of the simulation

### Accuracy of pROMs
#### Median fractional modal frequency error (from 100 runs) - as a function of eigenmode number - between the predictions of the full-order model and the three pROMs
![](https://github.com/FK-MAD/CMS/blob/main/Metsovo%20bridge%20results/errors_all_v3.svg?raw=true)

It can be seen that the maximum error for all pROMs is **$\sim10^{-2}$ or approximately 1%**.

#### Fractional modal frequency error between the predictions of the full-order model and the three pROMs at each sample point for the first four modal frequencies
![](https://github.com/FK-MAD/CMS/blob/main/Metsovo%20bridge%20results/errors_all_multi_v3%20-%201.svg?raw=true)

### Speed of pROMs
The mean computational times (from 100 runs) of a single simulation step for the full-order model and the three pROMs can be seen below.

| Model  | Full FE Model | NIR-pROM | GIR-pROM/SX1 | LIR-pROM/C |
| ------------- | :-------------: | :-------------: | :-------------: | :-------------: |
| Mean total iteration time [sec] | 116.5 | 32.6 | 6.9 | 0.02  |
| Speedup over the full FE model | 1x | 3.6x | 16.9x | 5825x  |

It can be seen that the fastest performing model achieves a speedup of **three orders of magnitude** over the full-order model.

# License
This work is licensed under a
[Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International Public License][cc-by-nc-nd].

[![CC BY-NC-ND 4.0][cc-by-nc-nd-image]][cc-by-nc-nd]

[cc-by-nc-nd]: http://creativecommons.org/licenses/by-nc-nd/4.0/
[cc-by-nc-nd-image]: https://licensebuttons.net/l/by-nc-nd/4.0/88x31.png
[cc-by-nc-nd-shield]: https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg
