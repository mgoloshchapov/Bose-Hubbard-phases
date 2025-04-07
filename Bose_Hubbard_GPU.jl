using ITensors, ITensorMPS
using NPZ
using Suppressor
using CUDA

function BoseHubbard(t, μ, U; N=50, max_occupation=5)
    sites = siteinds("Qudit", N; dim=max_occupation)
    os = OpSum()
    for j=1:N-1
        os += -μ,"n",j
        os += -U/2.0,"n",j
        os += U/2.0,"n",j,"n",j
        os += -t,"Adag",j,"A",j+1
        os += -t,"Adag",j+1,"A",j
    end

    H = MPO(os, sites)
    return H, sites
end;


function Simulation(H, sites; maxdim=30, nsweeps=30)
    psi0 = cu(random_mps(sites; linkdims=maxdim))

    cutoff = [1E-10]

    energy, psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)
    return energy, psi
end;


# Average occupation number computation
function AverageOccupation(psi)
    n = expect(psi, "n")
    return sum(n)/length(n)
end;


μ = 0.0:1.2:1.1
U = 1.0
t = 0.0:0.02:0.4
bonddim_list = [3]
nsweeps_list = [10]

base_filename = "data/"

@suppress for k in eachindex(bonddim_list)
    result = zeros((length(μ), length(t)))
    name = "PhaseDiagram_GPU:N=50_maxoccupation=5_D=$(bonddim_list[k])_nsweeps=$(nsweeps_list[k]).npz"
    filename = base_filename * name

    @suppress for i in eachindex(μ)
        @suppress for j in eachindex(t)
            H, sites = BoseHubbard(t[j], μ[i], U; max_occupation=5);
            energy, psi = Simulation(H, sites; maxdim=bonddim_list[k], nsweeps=nsweeps_list[k])
            result[i, j] = AverageOccupation(psi);
        end
        npzwrite(filename, Dict("mu" => μ, "t" => t, "U" => U, "navg" => result))
    end
end