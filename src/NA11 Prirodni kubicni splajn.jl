### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 2d895d76-dcb9-4d43-b985-b7e9c97d0cd0
begin
	using Polynomials
	using Plots
	using SpecialMatrices
end

# ╔═╡ cf690c86-f1eb-418f-992f-e569937ff499
begin
	# Generiranje točaka
	using Random, LinearAlgebra
	Random.seed!(123)
	n=5
	x=sort(rand(n+1))
	y=rand(n+1)
end

# ╔═╡ 5ec7a38e-cecd-45c0-94a1-39d58e437add
md"""
# Prirodni kubični splajn


Neka je zadana funkcija $f(x)$ na intervalu $[a,b]$.

Odaberimo $n+1$ točku 

$$
a\equiv x_0<x_1<x_2<\cdots <x_n\equiv b$$ 

i izračunajmo vrijednosti 

$$
y_i=f(x_i), \quad i=0,1,\ldots,n.$$
 
Na intervalu $[x_{i-1},x_i]$ funkciju $f$ aproksimiramo kubičnim polinomom $C_i$,
tako da je na intervalu $[a,b]$ funkcija $f$ aproksimirana funkcijom 

$$
C(x)=C_i(x), \quad x\in[x_{i-1},x_i]$$

Od funkcije $C(x)$ tražimo 

* __neprekidnost__
* __neprekidnost prve derivacije__ i
* __neprekidnost druge derivacije__.

Dakle,

$$\begin{aligned}
C_i(x_{i-1})&=y_{i-1}, \quad &i=1,\ldots,n, \\
C_i(x_{i})&=y_{i} \quad &i=1,\ldots, n,\\
C'_i(x_i)&=C'_{i+1}(x_i), \quad &i=1,\ldots,n-1, \\
C'_i(x_i)&=C'_{i+1}(x_i), \quad &i=1,\ldots,n-1,
\end{aligned}$$

pa imamo sustav od $4n-2$ jednadžbe i $4n$ nepoznanica (svaki od $n$ polinoma ima 4 koeficijenta).

Vrijede sljedeće tvrdnje:

$$
C_i(x)=y_{i-1}-s_{i-1}\frac{h_i^2}{6}+b_i(x-x_{i-1})+\frac{s_{i-1}}{6h_i}(x_i-x)^3
+\frac{s_i}{6h_i}(x-x_{i-1})^3,$$

gdje je 

$$\begin{aligned}
b_i&=d_i-(s_i-s_{i-1})\frac{h_i}{6},\\
d_i&=\frac{y_i-y_{i-1}}{h_i},\\
h_i&=x_i-x_{i-1},
\end{aligned}$$

a brojevi $s_i$, $i=0,1,\ldots,n$, zadovoljavaju sustav jednadžbi 

$$
s_{i-1}h_i+2s_i(h_i+h_{i+1})+s_{i+1}h_{i+1}=6(d_{i+1}-d_i),\quad i=1,\ldots,n-1.$$

Ako zadamo $s_0$ i $s_n$, sustav će imati jedinstveno rješenje. 

Najčešće su zadani __prirodni uvjeti__:

$$
s_0=0, \quad s_n=0.$$ 

U tom slučaju, $s_1,\ldots,s_{n-1}$ su rješenja sustava



$$
{\small
\begin{pmatrix} 2(h_1+h_2) & h_2 & 0 & \cdots & 0 & 0 \\
h_2 & 2(h_2+h_3) & h_3 & \cdots & 0 & 0 \\
0 & h_3 & 2(h_3+h_4) & \cdots & 0 & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & 0 & \cdots & 2(h_{n-2}+h_{n-1}) & h_{n-1} \\
0 & 0 & 0 & \cdots & h_{n-1}  & 2(h_{n-1}+h_{n})\\
\end{pmatrix}
\begin{pmatrix}
s_1\\ s_2 \\ s_3 \\ \vdots \\ s_{n-2} \\ s_{n-1}
\end{pmatrix}
= 
\begin{pmatrix}
6(d_2-d_1)\\
6(d_3-d_2)\\
6(d_4-d_3) \\
\vdots \\
6(d_{n-1}-d_{n-2}\\
6(d_n-d_{n-1})
\end{pmatrix}.}$$


Dokaz se nalazi u udžbeniku [Numerička matematika, str. 29](http://www.mathos.unios.hr/pim/Materijali/Num.pdf).

Matrica sustava je __tridijagonalna__ i __pozitivno definitna__ pa se sustav može riješiti metodom Choleskog (bez pivotiranja) u $O(n)$ operacija. 



Vrijede __ocjene pogreške__:

$$\begin{aligned}
\max |f(x)-C(x)| &\leq \frac{5}{384} \max h_i^4 \\
\max |f'(x)-C'(x)| &\leq \frac{1}{24} \max h_i^3 \\
\max |f''(x)-C''(x)| &\leq \frac{3}{8} \max h_i^2. 
\end{aligned}$$

Ocjene se mogu promatrati i na svakom intervalu posebno.
"""

# ╔═╡ c12ee2d1-c5fa-44de-92e5-3c4e096e0187
md"""
## Primjer - Interpolacija slučajnih točaka
"""

# ╔═╡ 7b7c78a0-19fc-11eb-1d5f-65535f27c308
# include("Vandermonde.jl")

# ╔═╡ a990adad-3abe-4dc3-be83-02e3ce5c46d4
function myspline(x,y)
    h=x[2:end]-x[1:end-1]
    d=(y[2:end]-y[1:end-1])./h
    H=SymTridiagonal(2*(h[1:end-1]+h[2:end]),h[2:end-1])
    b₀=6*(d[2:end]-d[1:end-1])
    s=H\b₀
    s=[0;s;0]
    # Definirajmo polinome
    b=d-(s[2:end]-s[1:end-1]).*h/6
	n=length(x)-1
    C=Array{Any}(undef,n)
    C=[xx -> 
        y[i]-s[i]*h[i]^2/6+b[i]*(xx-x[i])+s[i]*(x[i+1]-xx)^3/(6*h[i])+s[i+1]*(xx-x[i])^3/(6*h[i]) 
        for i=1:n]
    return C
end 

# ╔═╡ e0f49f72-b3c9-40e8-94ba-6646edb0966d
function plotspline(C,x,xx)
	# Točke na splajnu
	ySpline=Array{Float64}(undef,length(xx))
	for i=1:length(xx)
		for k=1:length(C)
			if xx[i]<=x[k+1]
				ySpline[i]=C[k](xx[i])
	            break
	        end
	    end
	end
	return ySpline
end

# ╔═╡ 1e6ed7c3-a1bb-4d85-8c48-18410365aa47
C=myspline(x,y)

# ╔═╡ bca31410-9afb-4b45-bc46-d4a98d926301
begin
	# Crtanje
	lsize=200
	xx=range(x[1],x[end],length=lsize)
	scatter(x,y,label="Točke")
	ySpline=plotspline(C,x,xx)
	plot!(xx,ySpline,label="Splajn")
end

# ╔═╡ 4afd29e7-8e7c-44b6-8897-62dd7317b983
md"""
Usporedimo splajn s interpolacijskim polinomom:
"""

# ╔═╡ e667e8e4-79ca-43cb-b0c1-2e673e13cfe3
begin
	A=Vandermonde(x)
	p=Polynomial(A\y)
	yPoly=p.(xx)
	scatter(x,y,label="Točke")
	plot!(xx,[ySpline yPoly],label=["Splajn" "Polinom"])
end

# ╔═╡ 4e9a4998-2a8b-4476-8b6f-d086bb5ca6c3
md"""
## Interpolacija funkcija

Usporedimo splajn i interpolacijski polinom za funkcije 

$$f(x)=\sin(x), \quad x\in[0,\pi]$$

i 

$$f(x)=1-|x-1|,\quad  x\in[0,2].$$

"""

# ╔═╡ 20f11b0f-110e-465a-84d5-0f78e99b14ff
begin
	# n₁=5; a=0; b=pi; f(x)=sin.(x)
	n₁=10; a=0; b=2; f(x)=1 .-abs.(x .-1)
	
	x₁=collect(range(a,stop=b,length=n₁+1))
	y₁=f(x₁)
	xx₁=collect(range(a,stop=b,length=lsize))
	
	# Polinom
	A₁=Vandermonde(x₁)
	p₁=Polynomial(A₁\y₁)
	yPoly₁=p₁.(xx₁)
	
	# Splajn
	C₁=myspline(collect(x₁),y₁)
	ySpline₁=plotspline(C₁,x₁,xx₁)
	
	# Funkcija 
	yFun=f(xx₁)
	
	# Crtanje
	scatter(x₁,y₁,label="Točke")
	plot!(xx₁,[ySpline₁ yPoly₁ yFun],label=["Splajn" "Polinom" "Funkcija"])
end

# ╔═╡ Cell order:
# ╟─5ec7a38e-cecd-45c0-94a1-39d58e437add
# ╟─c12ee2d1-c5fa-44de-92e5-3c4e096e0187
# ╠═2d895d76-dcb9-4d43-b985-b7e9c97d0cd0
# ╠═7b7c78a0-19fc-11eb-1d5f-65535f27c308
# ╠═cf690c86-f1eb-418f-992f-e569937ff499
# ╠═a990adad-3abe-4dc3-be83-02e3ce5c46d4
# ╠═e0f49f72-b3c9-40e8-94ba-6646edb0966d
# ╠═1e6ed7c3-a1bb-4d85-8c48-18410365aa47
# ╠═bca31410-9afb-4b45-bc46-d4a98d926301
# ╟─4afd29e7-8e7c-44b6-8897-62dd7317b983
# ╠═e667e8e4-79ca-43cb-b0c1-2e673e13cfe3
# ╟─4e9a4998-2a8b-4476-8b6f-d086bb5ca6c3
# ╠═20f11b0f-110e-465a-84d5-0f78e99b14ff
