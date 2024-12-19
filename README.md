# PBD (based on Wukong codebase)

This repository is forked from the WuKong codebase.

However, the only functionalities it actually uses are libigl, polyscope and the script to create a new project.
As a consequence, the other existing projects have been removed.

[//]: # (Depending on what you need you may choose to build the specific project by changing the CMakeLists.txt in the Project folder)

[//]: # (Ideally basic simulation models such as,)

[//]: # (-DiscreteShell, FEM2D/3D, EoLRods-)

[//]: # (should have only the basic implementation such that they can be inherited whenever needed. )

## Building from source in Linux

Clone the repository

> $ git clone https://github.com/P10010/Wukong2024

Make and build the project:

> $ cd Wukong2024
> 
> $ mkdir build
> 
> $ cd build
> 
> $ cmake ..
> 
> $ make

If building fails due to missing dependencies, check that all the packages listed in .devcontainer/Dockerfile are installed. 

Check the Dockerfile for other steps that might be needed,
such as the step to copy eigen /usr/local/include/ if you get messages telling you they can't find some eigen things even if it is installed. 

If your distribution uses apt, you can probably just copy and paste those commands.

Some issues, like not finding uint64_t or a missing cast, can be fixed "manually" by changing the related source code from Wukong depending on the messages you get while you try to execute make.

libmetis is part of apt. If your distro, does not have a libmetis package, you can build it manually from here https://github.com/KarypisLab/METIS. First you need a library called GKlib (it's explained there) AND you should follow the steps they suggested in this issue https://github.com/KarypisLab/METIS/issues/83

Run the project with the boat scene we have set up (from the build folder):

> 
> $ ./Projects/PBD/PBD ../Projects/PBD/data/sailboat/parts/sail.obj

Building has been tested with GCC 11.

## Instructions for Windows

Install WSL2 from Powershell or Command Prompt:
> $ wsl --install

From the WSL terminal, follow the instructions given for Linux.

[//]: # (From the WSL terminal, clone the repository:)

[//]: # (> $ git clone https://github.com/P10010/Wukong2024 --recurse-submodules)

[//]: # ()
[//]: # (Navigate to the repository folder:)

[//]: # (> $ cd Wukong2024)

[//]: # (Open in VSCode:)

[//]: # (> $ code .)

[//]: # ()
[//]: # (Rename .devcontiner/devcontainer_windows.json to .devcontainer/devcontainer.json, replacing the existing one. Follow the instructions for _Run Docker in VSCode_.)

[//]: # ()
[//]: # (Give WSL access to the internet for building the docker image. Run the following to open a file editor and replace the existing IP address with 8.8.8.8)

[//]: # (> $ sudo nano /etc/resolv.conf)

[//]: # (## Docker)

[//]: # ()
[//]: # (Download the docker image. Change tag "linux" if needed.)

[//]: # (> $ docker pull wukongsim/wukong:linux)

[//]: # ()
[//]: # (If you wish to build the docker image from scratch from the Dockerfile &#40;not recommended&#41; or rebuild after modifying the dockerfile, run the following in the command line from the directory Wukong2024/.devcontainer.)

[//]: # (> $ docker build -t wukongsim/wukong:linux .)

[//]: # ()
[//]: # (If finished modifying the Dockerfile, push to dockerhub.)

[//]: # (> $ docker push wukongsim/wukong:linux)

[//]: # ()
[//]: # (### Install NVIDIA Docker)

[//]: # ()
[//]: # (> $ curl -fsSL https://nvidia.github.io/libnvidia-container/gpgkey | sudo gpg --dearmor -o /usr/share/keyrings/nvidia-container-toolkit-keyring.gpg)

[//]: # ()
[//]: # (> $ curl -s -L https://nvidia.github.io/libnvidia-container/stable/deb/nvidia-container-toolkit.list | \)

[//]: # (    sed 's#deb https://#deb [signed-by=/usr/share/keyrings/nvidia-container-toolkit-keyring.gpg] https://#g' | \)

[//]: # (    sudo tee /etc/apt/sources.list.d/nvidia-container-toolkit.list)

[//]: # (    )
[//]: # (> $ sudo sed -i -e '/experimental/ s/^#//g' /etc/apt/sources.list.d/nvidia-container-toolkit.list)

[//]: # ()
[//]: # (> $ sudo apt-get update)

[//]: # ()
[//]: # (> $ sudo apt-get install -y nvidia-container-toolkit)

[//]: # ()
[//]: # (> $ sudo nvidia-ctk cdi generate --output=/etc/cdi/nvidia.yaml)

[//]: # ()
[//]: # (> $ sudo nvidia-ctk runtime configure --runtime=docker)

[//]: # ()
[//]: # (> $ sudo systemctl restart docker)

[//]: # ()
[//]: # (### Enable Display)

[//]: # ()
[//]: # (Enable Docker to connect to host display to spawn GUI windows. Run from command line on host machine, repeat in case of display error.)

[//]: # (> $ xhost +)

[//]: # ()
[//]: # (### Run Docker in VSCode)

[//]: # ()
[//]: # (Open the repository folder in VSCode. Install Docker and Dev Containers extensions in VSCode.)

[//]: # ()
[//]: # (In VSCode, type `control + p`, then type `>Reopen in Container` &#40;with the '>'&#41;. This option will show up in the >< tab in the bottom left corner of vscode.)

[//]: # ()
[//]: # (The above command will open a dev container using the docker image we provided.)

[//]: # ()
[//]: # (There you go! )

[//]: # ()
[//]: # (### Run Docker from Command Line)

[//]: # ()
[//]: # (Navigate to repository home directory.)

[//]: # (> $ cd Wukong2024)

[//]: # ()
[//]: # (Open the Docker container. \)

[//]: # (`-v ./:/{directory-name-in-container}` maps current working directory &#40;.&#41; to /{directory-name-in-container} in the Docker container. \)

[//]: # (`--network=host -e DISPLAY=$DISPLAY --privileged` enables use of host display to spawn GUI windows.)

[//]: # (> $ docker run -v ./:/{directory-name-in-container} -it --network=host -e DISPLAY=$DISPLAY --privileged --rm wukongsim/wukong:linux bash)

[//]: # ()
[//]: # (Build the code.)

[//]: # (> $ cd {directory-name-in-container})

[//]: # ()
[//]: # (> $ ./build.sh)

[//]: # ()
[//]: # (### Docker Image Building Time)

[//]: # (Building this docker image can take a while, for downloading MKL libraries and compiling SuiteSparse from the source code &#40;just to remove a single print&#41;. )

[//]: # (In case you have a powerful workstation, considering changing all the `make -j8` to `make -j128`.)

[//]: # (### Projects Tested Compiling)

[//]: # (- Discrete Shell [x] Linux [] MacOs)

[//]: # (- FEM3D  [x] Linux [] MacOs)

[//]: # (- EoLRods  [x] Linux [] MacOs)

[//]: # (- Isohedral Tiling  [x] Linux [] MacOs)

### Coding Convention

### Naming Convention

    NAMESPACE_EXAMPLE
    ClassExample
    functionExample
    variable_example
    TypenameExample



### More Info
If WuKong contributes to an academic publication, cite it as:
```bib
@misc{wukong,
  title = {WuKong},
  author = {Yue Li and others},
  note = {https://github.com/liyuesolo/Wukong2024},
  year = {2024}
}
```
