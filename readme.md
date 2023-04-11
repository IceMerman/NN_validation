# Neural network validation for power systems

This project aims to build and train neural networks (NN) for power system validation of operational points.

# Dataset

The dataset is generated in MatLab using the library [Matpower](https://matpower.org/)

## Methodology

The algorithm for [DB Generation](Demanda/PPflow_DBIA.m) consist of:  

    1. Run DC optimal power flow without the grid (infinite line capacity)
    2. Validate if the current DC solution is AC feasible and safe
       1. Check voltage maginite
       2. Check line capacity
       3. Check no load shedded
       4. Check no aditional reactive power usage
    3. Store the power output of the generators
    4. Compute AC power flow
    5. Store AC results
    6. Store load profile
    7. Store safety label (1 if the case is feasible an sure)

# Python exploratory analysis

In this repo you can find a series of python notebook where data exploration was made

# References

# TODO: 
 - [ ] Change algorith to change the generation in order to generate feasible scenarios