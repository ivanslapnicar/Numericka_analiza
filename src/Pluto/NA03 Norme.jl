### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 1ad65234-9ba8-4eaf-af4d-b881f8d5724d
begin
	using LinearAlgebra
	import Random
	Random.seed!(1244)
	x=rand(-9:9,5)
end

# ╔═╡ 54201b21-914a-4152-b187-3ecfe9982e65
md"""
# Norme


Općenito, __norma__ na vektorskom prostoru $X$ je svaka funkcija $\| \phantom{x} \| : X\to \mathbb{R}$ sa sljedećim svojstvima:

1.  $\| x\|=0 \Leftrightarrow x=0$
2.  $\| \lambda x\|=|\lambda| \|x\|$
3.  $\| x+y\| \leq \|x\|+\|y\|  \qquad$ (nejednakost trokuta)

"""

# ╔═╡ 92d26690-7f76-4a73-ac98-38c42f9c46f4
md"""
## Vektorske norme

Za $X=\mathbb{R}^n$ imamo

$$\|x\|_p=\big(\sum_{i=1}^n |x_i|^p\big)^{1/p}$$

Posebno:

*  $\|x\|_1=\sum_{i=1}^n |x_i|$
*  $\|x\|_2=\sqrt{\sum_{i=1}^n x_i^2}= \sqrt{x\cdot x}$
*  $\|x\|_\infty = \max\limits_{i=1,\ldots,n} |x_i|$
"""

# ╔═╡ c1ffcca3-d537-4c4f-86b1-87215325c115
norm(x,1), norm(x), norm(x,Inf)

# ╔═╡ 31a749ca-adb3-4e0e-9484-4506e1586a74
md"""
## Matrične norme

Iz svake vektorske norme možemo izvesti matričnu normu (__inducirane norme__):

$$\|A\| = \max\limits_{x\neq 0} \frac{\|Ax\|}{\|x\|}=\max\limits_{\|x\|=1} \|Ax\|$$

Posebno:

*  $\|A\|_1=\max\limits_{j=1:n} \sum_{i=1}^n |a_{ij}|$  - najveća 1-norma stupca
*  $\|A\|_{\infty}=\max\limits_{i=1:n} \sum_{j=1}^n |a_{ij}|$ - najveća 1-norma retka
*  $\|A\|_2$ - najveća singularna vrijednost  matrice $A$

__Frobeniusova__ ili __Euklidska__ norma

$$\|A\|_F =\sqrt{\sum_{i,j=1}^n a_{ij}^2}$$

nije inducirana norma.

Matrične norme još imaju i svojstvo 

$$
\|A\cdot B\|\leq \|A\| \cdot \| B\|.$$
"""

# ╔═╡ 365a36f4-e332-40ea-9113-d322dc47d480
A=rand(-4:4,5,5)

# ╔═╡ 10bc88f1-3675-4446-b31a-033848b6160e
norm(A,1), norm(A), norm(A,2), norm(A,Inf), opnorm(A),
maximum(svdvals(A)), opnorm(A,1), opnorm(A,Inf)

# ╔═╡ 196d75ee-c304-4501-8ee3-b685d2fc2a65
md"""
## Skalarni produkt, norma i ortogonalnost funkcija


__Skalarni produkt__ na vektorskom prostoru $X$ je svako preslikavanje 
$\cdot : X\times X \to \mathbb{R}$ sa sljedećim svojstvima:

1.  $x\cdot x\geq 0$
1.  $x\cdot x=0 \Leftrightarrow x=0$
2.  $x\cdot y=y\cdot x$
3.  $(\alpha x)\cdot y =\alpha (x\cdot y)$
3.  $(x+y)\cdot z=x\cdot z+y \cdot z$

Ukoliko je na vektorskom prostoru definiran skalarni produkt, normu možemo definirati kao

$$\|x\|=\sqrt{x\cdot x}.$$

Također, ako je $x \cdot y=0$ kažemo da su vektori $x$ i $y$ __međusobno ortogonalni (okomiti)__.  

Na primjer, standardna vektorska norma

$$\|x\|_2=\sqrt{\sum_{i=1}^n x_i^2}= \sqrt{x\cdot x}$$

je definirana pomoću skalarnog produkta vektora, 

$$x\cdot y=\sum_{i=1}^n  x_i y_i,$$

a vektori $x$ i $y$ su ortogonalni, odnosno $x\perp y$, ako je 
$x\cdot y=0$.

Skalarni produkt funkcija definiramo pomoću određenog integrala:

$$f\cdot g = \int_a^b f(x)g(x) \, dx.$$

Ostale definicije ostaju iste:

$$\| f\|_2= \sqrt{f\cdot f} = \sqrt{\int_a^b [f(x)]^2 \, dx},$$

$$f\perp g \Longleftrightarrow f\cdot g =0.$$
"""

# ╔═╡ 0b47171d-ff6e-4e85-acc0-7f5aec8c11f4


# ╔═╡ Cell order:
# ╟─54201b21-914a-4152-b187-3ecfe9982e65
# ╟─92d26690-7f76-4a73-ac98-38c42f9c46f4
# ╠═1ad65234-9ba8-4eaf-af4d-b881f8d5724d
# ╠═c1ffcca3-d537-4c4f-86b1-87215325c115
# ╟─31a749ca-adb3-4e0e-9484-4506e1586a74
# ╠═365a36f4-e332-40ea-9113-d322dc47d480
# ╠═10bc88f1-3675-4446-b31a-033848b6160e
# ╟─196d75ee-c304-4501-8ee3-b685d2fc2a65
# ╠═0b47171d-ff6e-4e85-acc0-7f5aec8c11f4
