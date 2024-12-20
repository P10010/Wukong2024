# PBD (based on Wukong codebase)

This repository is forked from the [WuKong](https://github.com/liyuesolo/Wukong2024) codebase.

However, the only functionalities it actually uses are libigl, polyscope and the script to create a new project.
As a consequence, the other existing projects have been removed.

## Building from source in Linux

Clone the repository

> $ git clone https://github.com/P10010/Wukong2024  --recurse-submodules

Make and build the project:

> $ cd Wukong2024
> 
> $ mkdir build
> 
> $ cd build
> 
> $ cmake -DCMAKE_BUILD_TYPE=Release ..
> 
> $ make -j8

If building fails due to missing dependencies, check that all the packages listed in `dependencies.txt` are installed. 

Depending on your distribution, you might need to copy Eigen into `/usr/local/include/`

Some issues, like not finding `uint64_t` or a missing `cast`, can be fixed "manually" by changing the related source code from Wukong depending on the messages you get while you try to execute `make`.

`libmetis` is part of apt. If your distro, does not have a libmetis package, you can build it manually from here https://github.com/KarypisLab/METIS. First you need a library called GKlib (it's explained there) AND you should follow the steps they suggested in this issue https://github.com/KarypisLab/METIS/issues/83

Run the project with the boat scene we have set up (from the `build/Projects/PBD` directory):

> $ cd build/Projects/PBD
> 
> $ ./PBD ../../../Projects/PBD/data/sailboat/parts/sail.obj

Building has been tested with GCC 11.

