# tensor-networks-julia
Tensor Network Applications in Julia

Preliminary study of Tensor Network applications in Julia.

To run (with Julia installed):
```
git clone https://github.com/DennisReddyhoff/tensor-networks-julia.git
cd tensor-networks-julia
julia --project=.
```

## 1) DMRG for Transverse Field Ising Model

[Pluto Notebook (Hosted on pluto.land)](https://pluto.land/n/ir51nblx)

[Pluto Noteboook (Github)](notebooks/dmrg.jl)

Also see [DMRG.jl](src/DMRG.jl) and [DMRGPlots.jl](src/DMRGPlots.jl)

- Custom Julia DMRG implementation following https://tensornetwork.org/mps/algorithms/dmrg/#toc_3
- Manual environment building and canonical forms vs ITensor
- Comparison of DMRG vs exact diagonalization
- Bond-by-bond optimisation and full sweeps

## 2) Tensor Networks for Computational Fluid Dynamics

[Pluto Notebook (Hosted on pluto.land)](https://pluto.land/n/2q35hmix)
You will need to run with Binder to use the interactive features

[Pluto Noteboook (Github)](notebooks/1d_cfd.jl)

- Using Tensor Networks to solve 1D Advection-Diffusion Equation
- Encoding fields as quantum amplitudes using MPS
- Time Evolving Block Decimation for evolution in time
- Interactive simulation using PlutoUI

### In Progress:
On branch `hubbard`
- 1D Hubbard Problem using iTensor
- 1D Hubbard Problem using exact solution

### TODO:
- Benchmarking
