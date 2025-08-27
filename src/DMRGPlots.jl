include("DMRG.jl")
using Plots

function get_error_mat(L_list, sweeps_list; dmrg_fun)
    errors = Dict()
	for L in L_list
	    for sweeps in sweeps_list
	        _, e_opt, e_exact = dmrg_fun(L; n_sweeps=sweeps)
	        errors[(L, sweeps)] = abs(e_opt - e_exact)
	    end
	end
    error_matrix = [errors[(L, sweeps)] for L in L_list, sweeps in sweeps_list]
    return error_matrix
end

function plot_heatmap(L_list, sweeps_list, error_matrix)
	plt_hm = heatmap(
				sweeps_list, L_list, error_matrix;
		        xlabel="Number of Sweeps",
		        ylabel="Chain Length L",
		        title="DMRG Energy Error vs Exact",
		        yflip=true,
		        color=:viridis,
		        clims=(minimum(error_matrix), maximum(error_matrix))
	)
    return plt_hm
end

function plot_v_exact(L_list, sweeps_list, error_matrix)
    plt_error = plot(xlabel="Number of Sweeps",
               ylabel="Energy Error",
               title="DMRG Energy Error vs Exact")
	
    for (i, L) in enumerate(L_list)
        y = error_matrix[i, :]  # ith row corresponds to L
        plot!(sweeps_list, y, label="L = $L", lw=2, marker=:o)
    end
    
    return plt_error
end

function save_plots(L_list, sweeps_list; maxdim=32, cutoff=1e-12, J=1.0, g=1.2)
    dmrg_fun = (L; n_sweeps) -> run_dmrg(L, J, g; maxdim=maxdim, cutoff=cutoff, n_sweeps=n_sweeps)
    error_mat = get_error_mat(L_list, sweeps_list; dmrg_fun=dmrg_fun)
    heatmap = plot_heatmap(L_list, sweeps_list, error_mat)
    error_plot = plot_v_exact(L_list, sweeps_list, error_mat)
    savefig(heatmap, "./plots/dmrg_heatmap.png")
    savefig(error_plot, "./plots/dmrg_v_exact.png")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    L_list = 3:6
    sweeps_list = [1, 10, 25, 50, 100]
    save_plots(L_list, sweeps_list; J=1.0, g=1.2)
    println("Plots saved to ./plots")
end