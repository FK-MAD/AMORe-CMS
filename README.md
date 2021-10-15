# CMS

Performs model reduction using component mode synthesis (CMS) based on the Craig-Bampton method. 
Works with 3D and 2D structural FE models created with COMSOL Multiphysics.

The code is written in MATLAB. It was developed for my diploma thesis "Advanced Model Reduction Techniques for Structural Dynamics Simulations" at the Department of Mechanical Engineering of the University of Thessaly. MATLAB must be started with a COMSOL server (using LiveLink).
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

# License

This work is licensed under a
[Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International Public License][cc-by-nc-nd].

[![CC BY-NC-ND 4.0][cc-by-nc-nd-image]][cc-by-nc-nd]

[cc-by-nc-nd]: http://creativecommons.org/licenses/by-nc-nd/4.0/
[cc-by-nc-nd-image]: https://licensebuttons.net/l/by-nc-nd/4.0/88x31.png
[cc-by-nc-nd-shield]: https://img.shields.io/badge/License-CC%20BY--NC--ND%204.0-lightgrey.svg
