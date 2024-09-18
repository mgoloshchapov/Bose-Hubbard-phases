using ITensors, ITensorMPS
using NPZ
using Suppressor
using CurveFit

function BoseHubbard(t, μ, U; N=51, max_occupation=5)
    sites = siteinds("Qudit", N; dim=max_occupation)
    os = OpSum()
    for j=1:N
        os += -μ,"n",j
        os += -U/2.0,"n",j
        os += U/2.0,"n",j,"n",j

        if j < N
            os += -t,"Adag",j,"A",j+1
            os += -t,"Adag",j+1,"A",j
        end
    end

    H = MPO(os, sites)
    return H, sites
end;


function ExtendedBoseHubbard(t, μ, U, V; N=51, max_occupation=5)
    sites = siteinds("Qudit", N; dim=max_occupation)
    os = OpSum()
    for j=1:N
        os += -μ,"n",j
        os += -U/2.0,"n",j
        os += U/2.0,"n",j,"n",j

        if j < N
            os += V, "n",j,"n",j+1
            os += -t,"Adag",j,"A",j+1
            os += -t,"Adag",j+1,"A",j
        end
    end

    H = MPO(os, sites)
    return H, sites
end;


function Simulation(H, sites; maxdim=30, nsweeps=30)
    psi0 = random_mps(sites; linkdims=maxdim)

    cutoff = [1E-10]

    energy, psi = dmrg(H,psi0;nsweeps,maxdim,cutoff)
    return energy, psi
end;


# Average occupation number computation
function Occupation(psi)
    n = expect(psi, "n")
    return n
end;


function AverageOccupation(psi)
    n = expect(psi, "n")
    return sum(n)/length(n)
end;


function Correlation(psi)
    C = correlation_matrix(psi, "Adag", "A")
    return C
end


function CorrelationCord(C; avg=false)
    L = size(C)[1]
    Cr = zeros(L)

    if avg
        for r in eachindex(Cr)
            for i in eachindex(Cr)
                Cr[r] += C[i, 1 + (i + r - 1)%(L)]/L
            end
        end
    else
        Cr = C[Int(ceil(L//2)),:]
    end

    return Cr
end


function CorrelationLengths(psi; avg=false, max_r=10)
    C = Correlation(psi)
    Cr = CorrelationCord(C; avg=avg)

    L = length(Cr)
    center=Int(ceil(L//2))
    r = 1:max_r
    a, b = linear_fit(r, log.(abs.(Cr[center+1:center+max_r])))
    c, d = linear_fit(log.(r), log.(abs.(Cr[center+1:center+max_r])))

    xi = -1.0/b
    K = -2.0*d
    return xi, K
end


function StructureFactor(psi)
    Cn = correlation_matrix(psi, "n", "n")
    L = size(Cn)[1]
    Spi = 0.0
    for i in 1:L
        for j in 1:L
            Spi += (-1)^abs(i-j)*Cn[i,j] / L^2
        end
    end

    return Spi
end


#U is assumed to be 1.0
function PhaseDiagram(μ, t, D, nsweeps, max_occupation; N=51, base_filename="data/")
    navg = zeros((length(μ), length(t)))
    xi = zeros((length(μ), length(t)))
    K = zeros((length(μ), length(t)))
    name = "PhaseDiagram:N=50_maxoccupation=5_D=$(D)_nsweeps=$(nsweeps).npz"
    filename = base_filename * name

    idx_list = []
    for i in eachindex(μ)
        for j in eachindex(t)
            push!(idx_list, (i,j))
        end
    end

    @suppress Threads.@threads for (i,j) in idx_list
        H, sites = BoseHubbard(t[j], μ[i], 1.0; N=N, max_occupation=max_occupation);
        energy, psi = Simulation(H, sites; maxdim=D, nsweeps=nsweeps)
        navg[i, j] = AverageOccupation(psi);
        xi[i, j], K[i, j] = CorrelationLengths(psi);
    end
    npzwrite(filename, Dict("mu" => μ, "t" => t, "U" => 1.0, "navg" => navg, "K" => K, "xi" => xi))
end


function PhaseDiagramEBHM(μ, t, V, D, nsweeps, max_occupation; N=51, base_filename="data/")
    navg = zeros((length(μ), length(t)))
    xi = zeros((length(μ), length(t)))
    K = zeros((length(μ), length(t)))
    Spi = zeros((length(μ), length(t)))
    ni = zeros((N, length(μ), length(t)))
    ninj = zeros((N, N, length(μ), length(t)))

    name = "EBHM_Correlators_N_$(N)_maxoccupation_5_D_$(D)_nsweeps_$(nsweeps).npz"
    filename = base_filename * name

    #Create list of indices, so not to think about double loop in threads
    idx_list = []
    for i in eachindex(μ)
        for j in eachindex(t)
            push!(idx_list, (i,j))
        end
    end

    @suppress Threads.@threads for (i,j) in idx_list
        H, sites = ExtendedBoseHubbard(t[j], μ[i], 1.0, V; N=N, max_occupation=max_occupation);
        energy, psi = Simulation(H, sites; maxdim=D, nsweeps=nsweeps)
        navg[i, j] = AverageOccupation(psi);
        xi[i, j], K[i, j] = CorrelationLengths(psi);
        Spi[i, j] = StructureFactor(psi)
        ni[:,i,j] = expect(psi, "n")
        ninj[:,:,i,j] = correlation_matrix(psi, "n", "n")
    end
    npzwrite(filename, Dict("mu" => μ, "t" => t, "U" => 1.0, "V" => V, "navg" => navg, "K" => K, "xi" => xi, "Spi" => Spi, "ninj" => ninj, "ni" => ni))
end