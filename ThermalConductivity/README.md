# Protocol to calculate the thermal conductivity ($\kappa$) of a matherial, applyed to a paraffine/graphene composite.


In this protocol we are going to compute the thermal conductivity $\kappa$ of a matherial using the open-source software [LAMMPS](https://www.lammps.org/).

As an example, we are going to use a paraffine/graphene nanocomposite, wich will be built from starting structures generated via the freely available [CHARMM-GUI](https://charmm-gui.org/) web service. The example structure will be built in section 1, and relaxed in step 2, so to apply the protocol to a general, already equilibrated structure, the user can skip to section 3 - onwards.

## 1. Create the starting structure data file

Via the [CHARMM-GUI](https://charmm-gui.org/) it is possible to obtain the LAMMPS data file for a pre-equilibrated melt of a desired polymer (in this case paraffine). Also, the website allows the creation of a LAMMPS data file of a pre-equilibrated graphene sheet matching the size of the setional area of the polymer melt (for our example, the graphene sheet will be aligned to the _xy_ plane of the polymer melt).

The script [01_stack_layers.py](00_scripts/01_stack_layers.py) (in the 00_scripts directory) will prepare a single LAMMPS data file starting with the CHARMM-GUI-generated files. It's possible to define in the script how many graphene sheets should be created, and the starting distance between the graphene and the paraffine melt. To simplify the tractation, in the script a simplified version of the GAFF forcefield will be used to parametrize the paraffine polymer, but hte user can easly modify the script to implement a different forcefield.

The user can find more informations within the comments of the script.

## 2. Equilibrate the starting structure

In the [02_equilibrate](02_equilibrate) directory there  are some LAMMPS scripts that can be used to minimize the starting structure, and than gradually equilibre its density, to prepare the system for the actual calculations. The simulations must be checked to find the proper lenght of the final equilibration run. 

The script [02_density_profile.py](00_scripts/02_density_profile.py) (in the 00_scripts directory) can be used to calculate the density profile of paraffine around the graphene sheets along the _z_ direction. 

## 3. Compute $\kappa$ with the Müller-Plathe method

In the [03_TC_MP](03_TC_MP) directory the user can find the LAMMPS imput scripts that can be used to compute $\kappa$ with the Müller-Plathe method. The data should be collected applying the heat flux along each dimension. The same LAMMPS input script can be used for each dimension, just changing the `dir` variable (wich can be changed drectly when invoking LAMMPS, being an _index_ variable). Other useful parameters can be easly changed in the first lines of the script, to find the best set for the system of interest.

The [01_lammps_run.sh](03_TC_MP/01_lammps_run.sh) script can be used to run the calculation along all the three directions sequentially, but the calculations could also be performed in parallel by adjusting the scripts accordingly. The user can find more informations within the comments of the LAMMPS input script.

Once the simulations are completed, the [03_check_MP.py](00_scripts/03_check_MP.py) script can be used to compute $\kappa$ (doing the units conversion to report its value in $W/(m*K)$ units). The script will also create some useful plots that can be used to check the simulation parameters. Specifically, the script will plot the final temperature gradient, and will compute $\kappa$ based on three different methods to compute the $\nabla T$. The rought $\nabla T$ is calculated as the difference in the temperature betweeh the hottest and coldest slab. The other two methods use the slope of the fitted line from the left (L) and right (R) side of the temperature gradient plot.

Care must be taken in the example system we used, as its $\kappa$ is anisotropic.

## 4. Compute $\kappa$ with the Green-Kubo approach

In the [04_TC_GK](04_TC_GK) the user can find the script to compute $\kappa$ with the Green-Kubo approach. With a single script, sequential simulations are performed in a loop, to obtain enough values to assess the raiability of the results.
The user can find more information within the LAMMPS input script in  the directory.


