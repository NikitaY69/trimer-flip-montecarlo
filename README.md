# Flip Monte Carlo algorithm for 3-dimensional trimer-mixtures

This repository provides a C++ package to run Flip Monte Carlo (FMC) simulations
of 3-dimensional trimer-mixtures. 

FMC was written
as a continuation of our [Swap Monte Carlo package for 2-dimensional polydisperse liquids](https://github.com/NikitaY69/swap-montecarlo). The latter was originally created during a Master 1 internship done at the _Laboratoire de Physique de l'École Normale Supérieure_ under the supervision of Camille Scalliet. We decided to separate the two projects to highlight its evolution and novel contributions.

## Compilation and execution
On each branch, the module can be compiled with
```
cd src/
g++ *.cpp -o EXEC_NAME -lstdc++fs -O2 -mfma -mbmi2 -flto -lboost_program_options
```
The executable is then parsed as follows
```
EXEC_NAME --input INPUT_FILE --outdir OUT_DIRECTORY --N 5 --T 0.04 --tau 100000 --tw 1 --cycles 1 --lin 50 --log 50 --p_swap 0.2 [--MSD --Cb --Fs --U]
```
- `input`: path to starting configuration with data structure `X Y SIGMA` repeated `N` times (see `tutorial/cfg.xy` for example)
- `outdir`: path of output directory
- `N`: number of interacting particles
- `T`: temperature in reduced units
- `tau`: single-run time (as in aging plots, time after which observables are reset)
- `tw`: waiting-time between each cycle (as in aging plots, time after which a new cycle is measured)
- `cycles`: number of cycles of decorrelation tours
- `lin`: number of lin-spaced checkpoints (configurations)
- `log`: number of log-spaced checkpoints (configurations AND observables)
- `p_swap`: SWAP probability
- `MSD`, `Cb`, `Fs`, `U`: flags to compute the _Mean-Squared Displacement_, the _bond-breaking correlation function_, the _self-part of the intermediate scattering function_ and the _average potential energy_


Because the model has a density of one, the size of the system is simply taken as $L=\sqrt{N}$. Particles then interact with each other in a box $[-L/2,L/2]\times [-L/2,L/2]$ with the minimum image convention. <br>
The chosen unit of time is one monte-carlo **sweep** (that is $N$ consecutive displacement/swap trials). The total number of steps is calculated as in $n_\texttt{steps}=\texttt{tw}*(\texttt{cycles}-1)+\tau$.

## How to use it ?
Our SMC module is usable for different purposes as can testify the various branches avalailable. The `master` branch can be used for both thermalization and production runs. The `observables-only` branch is made to compute observables AFTER a simulation has finished. A `continue` branch finally manages to carry on finished runs. 

P.S. once you have compiled an executable under each relevent branch, you no longer need to switch between branches as you can just continuously re-use your executables.

### Thermalization
Thermalization is ensuring the fact that starting from any initial configuration, equilibrium is reached at desired temperature. <br>
To check on equilibration, one must use $\texttt{tw}>1$ and $\texttt{cycles}>1$. Usually, one launches a first run with a `tw`=`cycles`=1 to have an approximate of the relaxation time $\tau_\alpha$ with $C_B$ or $F_s$. A good value for `tw` is then slightly below $\tau_\alpha$ using `cycles`>1; finally one reaches equilibrium when 1- no memory-effects are visible between different cycles and 2- the potential energy is constant.

### Production
Calculating most physical observables relies on having good ensemble statistics. To do so, after reaching equilibrium at desired temperature, one must get observables evolution from different starting configurations. To proceed, you must use $\texttt{tw}=\texttt{cycles}=1$ and select the observables of interest.<br>
Most of the time, you can also use a slurm job with arrays each one running a specific starting configuration.

### Observables-only
If you have ran a simulation without any observable flag, it is still possible to
compute them *after* the simulation has ended. Using the corresponding executable
(compiled from `observables-only`), the command is simpler than before
```
EXEC_NAME --outdir OUT_DIRECTORY [--MSD --Cb --Fs --U]
```
- `outdir`: SMC output directory (the simulation must be ended: if there is no
`params.txt` file or `configs/` directory, it will return an error)
- `MSD`, `Cb`, `Fs`, `U`: check the master command

An `obs.txt` file will automatically be created in the existing directory with all demanded observables.