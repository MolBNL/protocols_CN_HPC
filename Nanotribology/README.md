# Project description
In this protocol, we'll simulate the scratching of a silica tip on an alumina ($\mathrm{Al_2O_3}$) surface, representing an atomic force microscopy setup. The friction of the tip on the surface at different indentations of the tip will be computed.

The protocol is divided into different steps, and an illustrative bash script that can be used to run all the steps sequentially can be found in the [data](./data) directory.

In the [data](./data) directory are saved all the [LAMMPS](https://www.lammps.org/) scripts, which will also load some common settings that can be found in the three files, whose names begin with the string "00".
All scripts expect the following subdirectories to be present:
- Surf, where all the simulations pertaining to the surface equilibration will be saved
- Sys, where the equilibration of the whole system and the tip indentation simulations will be saved
- Scr, where the scratching simulations will be saved 
- FF, where the comb3 forcefield files can be found

## STEP 01
Firstly, the alumina ($\mathrm{Al_2O_3}$) surface is created, starting from the alumina lattice information. 
The user can modify the script [01_create_box_and_surf.lmpin](./data/01_create_box_and_surf.lmpin) to adjust the dimensions of the surface.

## STEP 02 to 05
The next steps will minimize and equilibrate sequentially the created surface. The user can easily adjust the duration of each simulation to make sure the system is properly equilibrated.

## STEP 06
Once the surface is properly equilibrated, the tip is created on top of the surface itself. The tip will be treated as a rigid body. The user can easily adjust the tip radius to properly model the system of interest. The script will automatically create a half-sphere of silica on top of the surface, placing it at the center of the surface on the x dimension. In the z dimension, at least a distance equal to the radius of the silica tip will be kept between the tip itself and the surface.

## STEP 07
Finally, the whole system is properly equilibrated; the user can easily adjust the time of the simulation.

## STEP 08 
In this step, the tip will be moved downwards toward the alumina surface. This will be done in two steps. Firstly, at a higher speed, the tip is lowered until it's located at 3 Å from the surface. Then the speed is decreased (the user can modify the speed if needed), and a "simulation loop" is started. 
In this loop, at each iteration, a reference distance will be selected from a list of values provided by the user (these values must strictly be in decreasing order). The tip will be lowered until its distance to the surface is really close to the reference distance. Then, a LAMMPS restart file will be written, the value of the reference distance will be updated to the next one, and the simulation will continue.

## STEP 09
At this step, the scratching of the surface along the y dimension will be performed. To select the indentation level for the simulation, the 'ind' variable can be defined directly in the command used to invoke LAMMPS (see the [lammps_run.sh](./data/lammps_run.sh)).
In the 'log' file produced by the simulation, the friction in the y directions is reported as the 'c_friction[2]' value.

