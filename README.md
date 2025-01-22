# TFMC
A command-line executable performing Flip Monte Carlo (FMC) simulations of glassy trimer mixtures. 

## Scientific overview 
Our project generalizes the particle-swap move ([A. Ninarello _et al._, Phys. Rev. X 7, 021039 (2019)](https://link.aps.org/doi/10.1103/PhysRevX.7.021039)) to 3-dimensional trimer-mixtures. The latter has proven to drastically speed-up the dynamics of glassy polydisperse systems in a range of dimensions. 

Our model consists of molecules composed of three particles A, B, C of equal mass $m=1$ and of diameters 
$\sigma_\mathrm{A} = 0.9,~\sigma_\mathrm{B} = 1.0,~\sigma_\mathrm{C} = 1.1$. Particles fill the simulation volume with density $\rho=1.2$. 
We employ periodic boundary conditions with the minimum image conventions in a cubic
box of size $L=(N/\rho)^{1/3}$.  

The interactions between particles (both intra-molecular and inter-molecular) are governed
by a purely repulsive Weeks-Chandler-Andersen (WCA) potential. For bounded particles within the same molecule, 
we use a purely attractive _Finitely Extensible Nonlinear Elastic_ (FENE) potential. 

For more informations about the scientific background and particularly our package dynamical speedup compared to standard simulation methods, the reader is invited to consult `docs/report.pdf`.

## Installation
Our project is guaranteed to run on linux-based systems as we personally use it 
on a cluster running a very old version of Debian. Moreover, we tested it for C++ 11, 14 and 17 on the latest version of Ubuntu using
github CI workflows and it works just as fine. The build tests for macos systems
did not pass. We are currently investigating the problem and we hope to address this issue very soon. The procedure to install our project is described next.

1. Clone the repository:
   ```bash
   git clone https://github.com/NikitaY69/trimer-flip-montecarlo.git
   cd trimer-flip-montecarlo/
   ```
2. Make sure to possess the following dependencies: 
    - Boost (latest version)
    - CMake (latest version)

    All other third_party libraries are directly fetched using CMake. We decided not to fetch Boost because of its long installation time. The building phase (see below) might work with older versions of CMake but we do not guarantee it. 

    To install Boost on __linux-based systems__, you can use 
    ```bash
    sudo apt-get install libboost-all-dev
    ```
    On __MacOS__, you can execute
    ```bash
    brew install boost
    ```
3. Build the project
    ```bash
    mkdir build
    cd build
    cmake ..
    make
    make install
    ```
    In case you are working on a remote machine and you don't have access to root files, you can specify the installation root by replacing the third command with 
    ```bash
    cmake -DCMAKE_INSTALL_PREFIX=/path/to/local ..
    ```
4. Run tests
   ```bash
   ctest --output-on-failure
   ```
## Usage
### Executable
TFMC provides a command-line executable that is parsed as follows
```bash
TFMC --init ${INPUT_FILE} --params ${PARAMS_JSON_FILE} [--observables U MSD Fs]
```
- `INPUT_FILE`: path to starting configuration with data structure `MOL_INDEX TYPE X Y Z` (see `tests/config/initconf.xyz` for a reference).

    The `MOL_INDEX` column is not mandatory as our script automatically supposes that particles are sorted in the order of the order of their molecule index. 

    `TYPE` particle type as in A B C (1 2 3)

    `X Y Z` coordinates

- `PARAMS_JSON_FILE`: simulation parameters (see `tests/params/params_template.json` for a reference). Next is explained the meaning of each parameter:
    - `N`: number of particles (must be provided)
    - `T`: temperature in reduced-units (default 2.0)
    - `cycles`: number of decorrelation cycles for dynamical observables (default 1)
    - `linPoints`: number of linearly-spaced configuration snapshots inside a single cycle (default 50)
    - `logPoints`: number of logarithmically-spaced configuration snapshots AND observables calculations inside a single cycle (default 50)
    - `p_flip`: flip-move probability attempt (default 0.2)
    - `rootdir`: path to output rootdir (must be provided)
    - `tau`: number of MC sweeps inside a single cycle (default 100000)
    - `tw`: waiting-time between 2 subsequent cycles (default 1)

- `U MSD Fs`: list of observables than can be computed: total potential energy, mean-squared displacement and self-part of the intermediate scattering function. If no `--observables` is provided, only configurations are written. Observables are computed in the order of their appearance, e.g `--observables Fs MSD` will output Fs before MSD.

We are currently developping an observables-only possibility where observables are computed ON TOP of configurations. The later has already been implemented but we did not massively test it.

### Outputs
TFMC outputs configurations in `rootdir/configs/` and observables in `rootdir/obs.txt`. Configs are written just like the `INPUT_FILE` but without the `MOL_INDEX` column. Observables are written with the data structure `t cycle obs1 obs2 ...`

## Documentation
The documentation of TFMC is available at [https://nikitay69.github.io/trimer-flip-montecarlo/](https://nikitay69.github.io/trimer-flip-montecarlo/).

## What TFMC is not ?
TFMC is NOT a post-processing package. The executable is solely designed to equilibrate configurations at desired temperature and compute observables on the fly. Because the `rootdir` automatically copies the `JSON_PARAMS_FILE`, deploying automatized post-processing workflows for massive datasets is extremely easy with Python. On top of that, the HEADER of the `obs.txt` eases the automatization of workflows. 