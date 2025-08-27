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

[Pluto Notebook (Hosted on pluto.land)](https://pluto.land/n/7zlqrbkw)

[Pluto Noteboook (Github)](notebooks/dmrg.jl)

Also see [DMRG.jl](src/DMRG.jl) and [DMRGPlots.jl](src/DMRGPlots.jl)

- Custom Julia DMRG implementation following https://tensornetwork.org/mps/algorithms/dmrg/#toc_3
- Manual environment building and canonical forms vs ITensor
- Comparison of DMRG vs exact diagonalization
- Bond-by-bond optimisation and full sweeps


### TODO:

- 1D Hubbard Problem using iTensor
- 1D Hubbard Problem using exact solution
- Benchmarking