Here, we provide some supplementary files for the manuscript[1]. 

- Movie 1 shows the simulated deformation process of different samples. To download the movie, one can click the item, then click "View raw". 

- Each folder contains the ABAQUS input files and the UEL subroutines

[1] Coupled magneto-mechanical growth in hyperelastic materials: surface patterns modulation and shape control in bio-inspired structures, Zhanfeng Li, Yafei Wang, Zuodong Wang, Chennakesava Kadapa, Mokarram Hossain, Xiaohu Yao, Jiong Wang

[10.1016/j.jmps.2025.106089](https://www.sciencedirect.com/science/article/pii/S0022509625000651)

# Simulation(needs expertise)

Each folder also contains necessary files to run the simulation using ABAQUS. To run the simulation, please make sure 

- Visual Studio and Parallel Studio are linked to ABAQUS, such that it can be used for subroutine development. [Here](https://www.researchgate.net/publication/349991987_Linking_ABAQUS_20192020_and_Intel_oneAPI_Base_Toolkit_FORTRAN_Compiler) is a linking guide from the internet. 

- Submit the job through ABAQUS COMMAND window


# Simulated results

### Example 1: uniaxial loading of a tube

![Tube](https://github.com/Jeff97/Magneto-growth-MixedFEM/blob/main/Tube.jpg)

### Example 2: mesh convergence test

![Test](https://github.com/Jeff97/Magneto-growth-MixedFEM/blob/main/Convergence.jpg)

### Example 3: modulation of surface patterns

![Surface patterns](https://github.com/Jeff97/Magneto-growth-MixedFEM/blob/main/Wrinkle.jpg)

### Example 4: the inversion process of the algal genus Volvox

![Volvox](https://github.com/Jeff97/Magneto-growth-MixedFEM/blob/main/Volvox.jpg)