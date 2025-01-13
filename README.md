# An Improved Suspension Balance Model

This project has the code for simulating the viscosimetric flow of non Brownian suspensions in the Stokes regime, implemented in [OpenFOAM-v7&reg;](https://openfoam.org/release/7/). It corresponds to the Suspension Balance Model presented in Schlatter et al. (2025), which features a frame-independent formulation of the particles’ normal stress with an improved momentum interpolation scheme that prevents numerical oscillations. The particle’s stress model also includes a local formulation for the microscopically generated extra stress, ensuring grid convergence.

# Usage

1 - Download or clone the project.

2 - Source the OpenFOAM bashrc file. For an installation directory OF_PATH:
```
OF_PATH$ source OpenFOAM-7/etc/bashrc
```

3 - Compile the code:

```
An-Improved-Suspension-Balance-Model$ ./Allwmake
```

Run example:

```
An-Improved-Suspension-Balance-Model$ cd tutorials/Poiseuille/Circular
An-Improved-Suspension-Balance-Model/tutorials/Poiseuille/Circular$ blockMesh
An-Improved-Suspension-Balance-Model/tutorials/Poiseuille/Circular$ suspensionBalanceFoam
```


# Reference

L. Schlatter, G. G. da Silva Ferreira, and P. L. da Cunha Lage. An improved suspension balance model applied to shear-induced phase segregation. International Journal of Multiphase Flow, 184:105120, 2025. https://doi.org/10.1016/j.ijmultiphaseflow.2024.105120.
