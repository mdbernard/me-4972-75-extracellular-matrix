# Design Requirements

## Objectives

1. The model shall be statistically equivalent to the known structure of ECM.
2. The model shall have statistically equivalent behavior to the known behavior of ECM.

## Metrics
3. The model shall have its structure validated as statistically equivalent to the known structure of ECM using two-point correlation.
4. The model shall have its structure validated as statistically equivalent to the known structure of ECM using lineal-path correlation.
5. The model shall have its behavior validated using LAMMPS and VMD by comparing the dynamics parameters describing the model and ECM in vivo. <**TODO:** define "dynamics parameters">

## Constraints
6. The model shall have a simulation runtime in LAMMPS less than <**TODO:** runtime upper limit>.
7. The model shall have a file size less than <**TODO:** file size upper limit>.
- <**TODO:** Create precision requirements based on parameters which will become apparent as we become more familiar with the LAMMPS software.>

## Criteria
8. The chosen model shall have better scores than alternative models using the mathematical methods defined in the Metrics section.