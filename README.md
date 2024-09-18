Consider a lattice with bosons on each site, where operators are defined as the following:
- $\hat{a}_{i}, \hat{a}^{\dagger}_{i}$ - boson annihilation and creation operators for site $i$,
- $\hat{n}_{i}$ - boson occupation number operator.

_Bose-Hubbard model includes three terms:_
- chemical potential term $\mu$, which tells the price of adding more particles to the system, 
- hopping amplitude $t$, which defines how likely the bosons are hopping between the nearest neighbour lattice sites
- onsite interaction potential $U$, which prevents bosons from occupying the same lattice site

_Bose-Hubbard hamiltonian:_

$\begin{equation}H = -t\sum_{\left<i, j\right>}\left(\hat{a}_{i}^\dagger\hat{a}_{j} + \hat{a}_{i}\hat{a}_{j}^\dagger\right) + \frac{U}{2}\sum_{i}\hat{n}_i(\hat{n}_{i} - 1) - \mu \sum_{i} \hat{n}_{i}\end{equation}$


If we include repulsion between nearest neighbours, we will end up with extended Bose-Hubbard model, which demonstrates several new phases even in one-dimensional case.

_Extended Bose-Hubbard hamiltonian:_

$\begin{equation}H = -t\sum_{\left<i, j\right>}\left(\hat{a}_i^\dagger\hat{a}_j + \hat{a}_i\hat{a}_j^\dagger\right) + \frac{U}{2}\sum_{i}\hat{n}_i(\hat{n}_i - 1) + \frac{V}{2}\sum_{i}\hat{n}_i\hat{n}_j - \mu \sum_i \hat{n}_i\end{equation}$

