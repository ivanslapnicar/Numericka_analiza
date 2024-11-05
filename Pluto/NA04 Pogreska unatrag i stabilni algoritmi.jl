### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 98b7578d-1bdb-4be9-a4de-e5cfef2df4db
md"""
# Teorija smetnje, pogreška unatrag i stabilni algoritmi

# Teorija smetnje

Odgovorimo na sljedeće pitanje:

> _Koliko se mijenja rezulat u ovisnosti o promjeni ulaznih podataka?_

Za funkciju $f(x)$ i za neki ulazni podatak $x$, želimo dobiti ocjenu za
__apsolutnu pogrešku__ u odnosu na __promjenu__ ulaznog podatka za $\delta x$,

$$
\|f(x+\delta x)-f(x)\|\leq \kappa \|\delta x\|,$$

 i __relativnu pogrešku__ u odnosu na __relativnu promjenu__ ulaznog podataka $\displaystyle\frac{\| \delta x\|}{\|x\|}$,

$$
\frac{\| f(x+\delta x)-f(x)\|}{\| f(x) \|}\leq \kappa \frac{\| \delta x \|}{ \|x\|},$$


Vrijedi

$$
\| f(x+\delta x)-f(x)\| = \frac{\| f(x+\delta x)-f(x)\|}{\| \delta x \|} \|\delta x\| \equiv \kappa \|\delta x\|.$$

Veličina $\kappa$ je __uvjetovanost__ ili __kondicija__. Ona podsjeća na derivaciju, a kazuje koliko se najviše uveća smetnja u ulaznim podacima.

Slično, u izrazu

$$
\frac{\| f(x+\delta x)-f(x)\|}{\| f(x) \|}= \frac{\| f(x+\delta x)-f(x)\|\cdot  \|x\| }{\|\delta x\| \cdot\| f(x)\|}
\cdot \frac{\|\delta x\|}{\|x\|} \equiv \kappa \frac{\|\delta x\|}{\|x\|}$$

 $\kappa$ nam kazuje koliko se najviše relativno uveća relativna smetnja u ulaznim podacima.
"""

# ╔═╡ 95835e70-25e3-4dc7-bbf9-8d5f23ce3e44
md"""
# Pogreška unatrag

Neka vrijednost funkcije $f(x)$ računamo pomoću algoritma $\mathrm{alg}(x)$.
__Pogreška algoritma__ je 

$$
\|\mathrm{alg}(x)-f(x)\|,$$

a __relativna pogreška algoritma__ je 

$$
\frac{\| \mathrm{alg}(x)-f(x)\|}{\| f(x) \|}.$$

Ove pogreške je teško ili čak nemoguće procijeniti direktno. Stoga se promatra __pogreška unatrag__,

$$
\mathrm{alg}(x)=f(x+\delta x),$$
odnosno

> izračunata vrijednost funkcije $f$ za ulazni podatak $x$ jednaka je točnoj vrijednosti funkcije $f$ u smetanom ulaznom podatku za neku (nepoznatu) smetnju. 
"""

# ╔═╡ 63a8eeef-913c-4b0c-8cbc-9db8a07892b3
md"""
# Stabilni algoritmi

Algoritam je __stabilan__ ako uvijek vrijedi 

$$
\mathrm{alg}(x)=f(x+\delta x)$$

za neki "mali" $\delta x$.
"""

# ╔═╡ Cell order:
# ╟─98b7578d-1bdb-4be9-a4de-e5cfef2df4db
# ╟─95835e70-25e3-4dc7-bbf9-8d5f23ce3e44
# ╟─63a8eeef-913c-4b0c-8cbc-9db8a07892b3
