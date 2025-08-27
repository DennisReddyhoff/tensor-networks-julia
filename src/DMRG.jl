using ITensors
using ITensorMPS
using KrylovKit
using LinearAlgebra

#======== Basic DMRG Sweep Using 1D TFIM =========

The TFIM is defined as

$H = -J * \sum_{i=1}^{L-1} Sz_i * Sz_{i+1}  - g * \sum_{i=1}^{L} Sx_i$

for a `spin-1/2` system, where 
	$S_z$ and $S_x$ are Pauli gates, 
	$L$ is the number of spins in the chain, 
	$J$ is the nearest neighbour in $z$ and 
	$g$ is the transverse field along $x$. 

==================================================#

# Build and return H and sites for given L
function get_H(L=3, J=1.0, g=1.2)

	sites = siteinds("S=1/2", L) # spin 1/2 indices
	
	# Hamiltonian using ITensors OpSum
	# H = -J * sum_{i=1}^{L-1} Sz[i] * Sz[i+1]  - g * sum_{i=1}^{L} Sx[i]

	ops = OpSum()
	# Ising interactions between pair (i, i+1)
	# -J * sum_{i=1}^{L-1} Sz[i] * Sz[i+1]
	# Sz is Pauli-Z
	for i in 1:L-1
		ops += (-J, "Sz", i, "Sz", i+1)
	end
	# Transverse field
	# - g * sum_{i=1}^{L} Sx[i]
	# Pauli-X 
	for i in 1:L
		ops += (-g, "Sx", i)
	end
	H = MPO(ops, sites)
	return H, sites
end

# DMRG sweep, left -> right -> left
function dmrg_sweep!(M, H; maxdim=4, cutoff=1e-12)
    
	L = length(M)
    energies = Float64[]
    PH = ProjMPO(H)

    # Sweep left→right
    for i in 1:L-1
        # Canonical center = left site of the two-site block
        position!(PH, M, i)

        # Two-site tensor
        psi = M[i] * M[i+1]

        # Solve projected Hamiltonian
        vals, vecs = eigsolve(PH, psi, 1, :SR, krylovdim=20)
        gs = vecs[1]
        energy = vals[1]
        push!(energies, real(energy))

        # SVD split
        U, S, V = svd(gs, siteind(M, i); maxdim=maxdim, cutoff=cutoff)
        M[i] = U
        M[i+1] = S * V
    end
	
    # Sweep right→left
    for i in (L-1):-1:1
        # Canonical center = left site of two-site block
        position!(PH, M, i)

        psi = M[i] * M[i+1]

        vals, vecs = eigsolve(PH, psi, 1, :SR, krylovdim=20)
        gs = vecs[1]
        energy = vals[1]
        push!(energies, real(energy))

        U, S, V = svd(gs, siteind(M, i); maxdim=maxdim, cutoff=cutoff)
        # S is still absorbed in to right site
        M[i] = U 
        M[i+1] = S * V
    end

    return M, energies
end

# Get exact solution of H, only good for small L
function get_E_exact(H, sites)
	H_dense = contract(H)
	rows = combiner(sites'; tags="row")
	cols = combiner(sites;  tags="col")
	H_mat = H_dense * rows * cols
	H_array = Array(H_mat, combinedind(rows), combinedind(cols))
	vals_exact, vecs_exact = eigen(Hermitian(H_array))
	return minimum(vals_exact)
end

# Run n_sweeps full sweeps
function run_dmrg(L, J=1.0, g=1.2; maxdim=4, cutoff=1e-12, n_sweeps=10)
    H, sites = get_H(L, J, g)
    M = randomMPS(sites, maxdim)
	E_now = Inf
    for sweep in 1:n_sweeps
        M, sweep_energies = dmrg_sweep!(M, H; maxdim=maxdim, cutoff=cutoff)
        E_now = minimum(sweep_energies)
    end
	
    return M, E_now, get_E_exact(H, sites)
end


