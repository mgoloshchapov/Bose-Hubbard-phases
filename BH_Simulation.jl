using ITensors, ITensorMPS
using NPZ
using Suppressor
using "BH_Simulation.jl"


#U is assumed to be 1.0
function PhaseDiagram(μ, t, D, nsweeps, max_occupation, base_filename="data/")
    navg = zeros((length(μ), length(t)))
    K = zeros((length(μ), length(t)))
    xi = zeros((length(μ), length(t)))
    name = "PhaseDiagram:N=50_maxoccupation=5_D=$(D)_nsweeps=$(nsweeps).npz"
    filename = base_filename * name

    @suppress Threads.@threads for i in eachindex(μ)
        @suppress Threads.@threads for j in eachindex(t)
            H, sites = BoseHubbard(t[j], μ[i], 1.0; max_occupation=max_occupation);
            energy, psi = Simulation(H, sites; maxdim=D, nsweeps=nsweeps)
            navg[i, j] = AverageOccupation(psi);
            K[i, j], xi[i, j] = CorrelationLengths(psi);
        end
        npzwrite(filename, Dict("mu" => μ, "t" => t, "U" => 1.0, "navg" => result, "K" => K, "xi" => xi))
    end
end
