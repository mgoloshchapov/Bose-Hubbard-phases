### Bose-Hubbard and extended Bose-Hubbard models

Consider a lattice with bosons on each site, where operators are defined as the following:
- $\hat{a}_i, \hat{a}^{\dagger}_i$ - boson annihilation and creation operators for site $i$,
- $\hat{n}_{i}$ - boson occupation number operator.

_Bose-Hubbard model includes three terms:_
- chemical potential term $\mu$, which tells the price of adding more particles to the system, 
- hopping amplitude $t$, which defines how likely the bosons are hopping between the nearest neighbour lattice sites
- onsite interaction potential $U$, which prevents bosons from occupying the same lattice site

_Bose-Hubbard hamiltonian:_

$H = -t\sum_{\left<i, j\right>}\left(\hat{a}_i^\dagger\hat{a}_j + \hat{a}_i\hat{a}_j^{\dagger}\right) + \frac{U}{2}\sum_i\hat{n}_i(\hat{n}_i - 1) - \mu \sum_i \hat{n}_i$


If we include repulsion between nearest neighbours, we will end up with extended Bose-Hubbard model, which demonstrates several new phases even in one-dimensional case. Besides Mott insulator and superfluid phases present in the original Bose-Hubbard model, you can observe density wave, supersolid and Haldane insulator phases in 1D extended Bose-Hubbard model. 

_Extended Bose-Hubbard hamiltonian:_

$H = -t\sum_{\left<i, j\right>}\left(\hat{a}_i^\dagger\hat{a}_j + \hat{a}_i\hat{a}_j^{\dagger}\right) + \frac{U}{2}\sum_i\hat{n}_i(\hat{n}_i - 1) - V\sum n_i n_j - \mu \sum_i \hat{n}_i$



### DMRG Calculations and Correlators

Using DMRG and Julia ITensor package, we can compute groundstates of the $N$-site 1D Bose-Hubbard and extended Bose-Hubbard models in MPS form. After that we can calculate correlation matrices and different order parameters, including: 

- $\Gamma_{i,j} = \left< \hat{a}_i^\dagger \hat{a}_j\right>$ - correlation matrix for creation and annihilation operators

- $n_{i,j} = \left< \hat{n}_i \hat{n}_j\right>$ - correlation matrix for occuptaion numbers

- $n_{i} = \left< n_i\right>$ - occupation number of sites

- $\Gamma(r) = \Gamma_{N/2, N/2+r}$ - correlation function 

- $\Gamma(r) \sim e^{-r/\xi}, \\; \Gamma \sim r^{-K/2}$ - correlation lengths $\xi$ and $K$ from fits to exponential and polynomial decay correspondingly

- $S_{\pi} = \frac{1}{N^2}\sum_{i,j} (-1)^{|i-j|} \left< \hat{n}_i \hat{n}_j\right>$ - structure factor



