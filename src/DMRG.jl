using ITensors
using ITensorMPS

# Problem Definition

L = 4 # N spins in chain
J = 1.0 # Nearest neighbour in z
g = 1.2 # transverse field along x
sites = siteinds("S=1/2", L) # spin 1/2 indices

# Hamiltonian using ITensors OpSum
# H = -J * sum_{i=1}^{L-1} Sz[i] * Sz[i+1]  - g * sum_{i=1}^{L} Sx[i]
begin
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
	nothing
end

# Step 0

maxdim = 16
A = randomMPS(sites; linkdims=maxdim)
A = orthogonalize(A, 1) # Right-canonical with centre at 1


# Build left env
function build_L(A, H)
	Lenv = Vector{ITensor}(undef, L+1)
	Lenv[1] = ITensor(1.0)
	for n in 1:L-1
		Lenv[n+1] = Lenv[n] * A[n] * H[n] * dag(prime(A[n], "Site"))
	end
	return Lenv
end													  

# Build right env
function build_R(A, H)
	Renv = Vector{ITensor}(undef, L+2)
	Renv[L+1] = ITensor(1.0)
	for n in L:-1:2
		Renv[n] = H[n] * A[n] * dag(prime(A[n], "Site")) * Renv[n+1]
	end
	return Renv
end

# Step 1a

PH = ProjMPO(H)
position!(PH, A, 1)

phi = A[1] * A[2]

vals, vecs = eigsolve(PH, phi, 1, :SR)

ground_state = vecs[1]
energy = vals[1]

# Step 1b

U, S, V = svd(ground_state, siteind(A, 1); maxdim=maxdim, cutoff=1e-12)
A[1] = U * S
A[2] = V