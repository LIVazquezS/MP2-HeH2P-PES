# MP2 potential energy surface of HeH2p

Authors: L.I. Vazquez-Salazar and M. Meuwly

This repository contains the data to construct the potential energy surface at MP2 level with the basis set aug-cc-pVTZ for the HeH_{2}^{+} system.

## Ab initio calculations
The data is saved in the file `pes_mp2.csv`. The first column is the value of the angle (as $angle = 1-cos(\theta)/2$, second is the coordinate $R$, third is the coordinate $r$ and the last column is the value of the energy in eV.
To fit the PES, you need the RKHS package that can be found in [RKHS](https://github.com/MeuwlyGroup/RKHS). To do the fitting of the PES follow instructions there.

## Evaluation of the PES
Complementary, a `pes_mp2.kernel` file containing the coefficients of the kernel fitting is provided. This file can be used to evaluate the potential energy surface using the RKHS package.
The file `evaluate_kernel.f90` contains the code to evaluate the PES in a defined grid $\theta \in [0,180]$, $r \in [0,10]$ and $R \in [0,20]$. 
The file `evalute_kernel.f90` can be compiled with the adjunct make file. It requires the RKHS module in the same directory. 
The output will be `evaluate_kernel.x` that can be executed as follows

```
./evaluate_kernel.x pes_mp2.kernel new_pes.csv
```
Which will create a file `new_pes.csv` with the values of the PES in the defined grid.

## Fitting of 2-body potentials

The potential for HeH+ and H2+ are fitted following the analytical expression in: Phys.Chem.Chem.Phys. 2019, 21, 24976. 
In the folder Diatomic, the file `fitting.py` can be used to obtain the coefficients of the analytical expression. A plot with the correlation plot and the fitting of the curve will be showed.
The file `h2p.csv` contains the data for the H2+ system and `heh2p.csv` contains the data for the HeH+ system.

## Contact

For any questions, please contact  Markus Meuwly (m.meuwly@unibas.ch) or Luis Vazquez-Salazar (luisitza.vazquezsalazar@unibas.ch) 












