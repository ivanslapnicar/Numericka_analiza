### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ 66ac2c80-5d8d-4427-a002-14f14bdc2535
using PlutoUI, LinearAlgebra, Random

# ╔═╡ 52742d2f-73d7-4b5c-8f48-69f21d345133
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ fb66b990-4511-476f-8e77-87aa9c265041
md"""
# QR rastav


__QR rastav__ matrice $A$ tipa $m\times n$,  $m\geq n$,
glasi

$$
A=QR,$$

pri čemu je $Q$ __ortonormirana matrica__ dimenzije $m\times m$, odnosno

$$
Q^TQ=Q Q^T=I,$$

a $R$ je $m\times n$ gornje trokutasta matrica.

Ortonormiranu matricu kraće zovemo i __ortogonalna matrica__.

Na primjer,

$$\begin{aligned}
\begin{bmatrix} a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33} \\
a_{41} & a_{42} & a_{43} \\
a_{51} & a_{52} & a_{53}
\end{bmatrix}=
\begin{bmatrix}
q_{11} & q_{12} & q_{13} & q_{14} & q_{15} \\
q_{21} & q_{22} & q_{23} & q_{24} & q_{25} \\
q_{31} & q_{32} & q_{33} & q_{34} & q_{35} \\
q_{41} & q_{42} & q_{43} & q_{44} & q_{45} \\
q_{51} & q_{52} & q_{53} & q_{54} & q_{55}
\end{bmatrix}
\begin{bmatrix}
r_{11} & r_{12} & r_{13} \\
0 & r_{22} & r_{23} \\
0 & 0 & r_{33} \\
0 & 0 & 0 \\
0 & 0 & 0 
\end{bmatrix}. 
\end{aligned}$$

S ovim je definiran i __ekonomični QR rastav__

$$\begin{aligned}
\begin{bmatrix} a_{11} & a_{12} & a_{13} \\
a_{21} & a_{22} & a_{23} \\
a_{31} & a_{32} & a_{33} \\
a_{41} & a_{42} & a_{43} \\
a_{51} & a_{52} & a_{53}
\end{bmatrix}=
\begin{bmatrix}
q_{11} & q_{12} & q_{13} \\
q_{21} & q_{22} & q_{23} \\
q_{31} & q_{32} & q_{33} \\
q_{41} & q_{42} & q_{43} \\
q_{51} & q_{52} & q_{53}
\end{bmatrix}
\begin{bmatrix}
r_{11} & r_{12} & r_{13} \\
0 & r_{22} & r_{23} \\
0 & 0 & r_{33}
\end{bmatrix}.
\end{aligned}$$

# Gram-Schmidtova ortogonalizacija

Izjednačavanje stupaca počevši od prvog daje: 

$$
\begin{aligned}
t&=a_{:1}\\
r_{11}&=\|t\|_2 \\
q_{:1}&=t\frac{1}{r_{11}}\\
r_{12}&= q_{:1}^Ta_{:2} \\
t&=a_{:2}-q_{:1}r_{12} \\
r_{22}&=\|t\|_2 \\
q_{:2}&=t\frac{1}{r_{22}} \\
r_{13}&=q_{:1}^Ta_{:3} \\
r_{23}&=q_{:2}^Ta_{:3} \\
t&=a_{:3}-q_{:1}r_{13}-q_{:2}r_{23}\\
r_{33}&=\|t\|_2 \\
q_{:3}&=t\frac{1}{r_{33}}.
\end{aligned}$$

Indukcijom slijedi __Gram-Schmidtov postupak ortogonalizacije__.
"""

# ╔═╡ 24c9f0ba-6e18-4362-840b-1a3581032ff5
function GramSchmidtQR(A::Array)
	m,n=size(A)
	R=zeros(n,n)
	Q=Array{Float64}(undef,m,n)
	R[1,1]=norm(A[:,1])
	Q[:,1]=A[:,1]/R[1,1]
	for k=2:n
	    for i=1:k-1
	        R[i,k]=Q[:,i]⋅A[:,k]
	    end
	    t=A[:,k]-sum([R[i,k]*Q[:,i] for i=1:k-1])
	    R[k,k]=norm(t)
	    Q[:,k]=t/R[k,k]
	end
	Q,R
end 

# ╔═╡ b6ee2c00-78e7-4d22-a851-d22c4a4946bc
begin
	Random.seed!(123)
	A=rand(8,5)
end

# ╔═╡ b8994eec-d651-4921-b1a3-184a029e83b7
Q,R=GramSchmidtQR(A)

# ╔═╡ 19a0f0b8-dbab-4609-9f7f-17339af6bf63
Q

# ╔═╡ 54f456e7-133d-4f96-901f-5bcda000178f
Q'*Q

# ╔═╡ 1d11a929-212b-42d0-ab52-c793ed8905d2
R

# ╔═╡ 9038e960-55c0-437a-a9fd-73f48860a36a
# Rezidual
A-Q*R

# ╔═╡ b9bf7235-dcb8-47bf-8c3f-1179ba9e08e0
md"""
Algoritam `GramSchmidtQR()` je numerički nestabilan pa je bolje koristiti __modificirani Gram-Schmidtov algoritam__ ili __Householderove reflektore__ ili __Givensove rotacije__ (vidi [Matrix Computations, poglavlje 5](https://books.google.hr/books?id=X5YfsuCWpxMC&printsec=frontcover&hl=hr#v=onepage&q&f=false)).
"""

# ╔═╡ 3fd60d80-25d7-45e2-9776-7250841e8f74
md"""
# Householderovi reflektori

__QR rastav vektora__ $x$ jednak je

$$
H \begin{bmatrix} x_1 \\ x_2 \\ \vdots \\ x_m 
\end{bmatrix}  =r,$$

gdje je 

$$
H=I - \frac{2}{v^Tv}v v^T, \qquad  
v=\begin{bmatrix}
x_1\pm \|x\|_2 \\ x_2 \\ x_3 \\ \vdots \\ x_m
\end{bmatrix}.$$ 

__Householderov reflektor__ $H$ je __simetrična__ i __ortogonalna__ matrica (__dokažite!__). Ovisno o izboru predznaka u definicije vektora $v$ vrijedi

$$
r=\begin{bmatrix} \mp \|x\| \\ 0 \\ \vdots \\ 0
\end{bmatrix}.$$

Zbog numeričke stabilnost se najčešće uzima

$$
v_1=x_1+\mathop{\mathrm{sign}} (x_1) \|x\|_2.$$

Matrica $H$ se __ne računa eksplicitno__ već se produkt $Hx$ računa po formuli

$$
Hx=x-\frac{2(v^Tx)}{v^Tv}v=x-\frac{2 (v\cdot x)}{v\cdot v}v$$

za koju je potrebno $O(6m)$ operacija.
"""

# ╔═╡ 4a371793-1d5c-413f-8e89-dce986c972ef
function HouseholderVector(x::Array)
    # Računa v
    v=copy(x)
    v[1]=x[1]+sign(x[1])*norm(x)
    v
end

# ╔═╡ 814c6dd6-faa1-4211-8f64-bf349624dc37
begin
	x=rand(8)
	v=HouseholderVector(x)
	β=(2/(v⋅v))*(v⋅x)
	x-β*v
end

# ╔═╡ e282bbd3-6505-4afc-b304-3c327fa5db62
x

# ╔═╡ eb8f508a-9018-4896-bd82-e707e01157bc
norm(x)

# ╔═╡ df0bb186-4b2c-429a-b279-5843f97cd3f0
md"""
QR rastav matrice se računa rekurzivnim QR rastavom vektora pomoću Householderovih reflektora:
"""

# ╔═╡ 3bb4d086-a79a-44b1-92c3-5e29be4c36de
function HouseholderQR(A₁::Array)
    # Računa Q i R
    A=copy(A₁)
    m,n=size(A)
    Q=Matrix{Float64}(I,m,m) # eye
    for k=1:n
        v=HouseholderVector(A[k:m,k])
        β=(2/(v⋅v))*v
        A[k:m,k:n]=A[k:m,k:n]-β*(v'*A[k:m,k:n])
        Q[k:m,:]=Q[k:m,:]-β*(v'*Q[k:m,:])
    end
    R=triu(A)
    Q',R
end
    

# ╔═╡ d7ddbafd-8923-45c1-a0f2-55aa1dc86185
A

# ╔═╡ cb302477-773e-4fd9-a9d6-a9df50e45b18
Qₕ,Rₕ=HouseholderQR(A)

# ╔═╡ 12069de5-8064-437b-b96c-222baa7b4c7d
Qₕ'*A

# ╔═╡ 3ceae442-3a07-4400-933c-a04bdc8ee742
Rₕ

# ╔═╡ ef4702ee-44cc-49d8-85bb-39249c431036
md"""
Program `myHouseholderQR()` je ilustrativan. Profesionalni programi imaju sljedeća svojstva:

* računaju s blok matricama (uobičajena dimenzija bloka je 32 ili 64),
* izračuna se vektor $\hat v=v/v_1$. Vrijedi $\hat v_1=1$, dok se ostali elemenenti vektora $\hat v$ spremaju u strogi donji trokut matrice $A$,
* ako se traži matrica $Q$, akumulacija se vrši unatrag koristeći spremljene vektore $v$ (tako se smanjuje broj operacija),
* postoji opcija vraćanja ekonomičnog rastava,
* postoji opcija računanja s __pivotiranjem__ - u svakom koraku se na prvo mjesto dovede stupac s najvećom normom pa je 

$$AP=QR,\quad |R_{kk}|\geq |R_{k+1,k+1}|$$

pa se može utvrditi i __numerički rang__ matrice.
"""

# ╔═╡ 4587b182-80fd-4bd7-87af-17074ea7c9d9
# ?qr # Pogledajmo upute

# ╔═╡ f944cdbe-794c-405e-a215-687cbb12b22b
# Izračunajmo QR objekt
F=qr(A)

# ╔═╡ 64e052c2-b798-4fb8-80fd-c018df2c8625
F.Q'*A

# ╔═╡ ef2ebf49-e7b2-4300-9a5c-3a2e42251e74
F.Q*F.R

# ╔═╡ 6bf029d7-68e2-47a4-87b9-8f66154933d4
Fₚ=qr(A,Val(true))

# ╔═╡ 68d0e359-9fb9-4a9b-9f13-e8039597e948
# Vektor pivotiranja
Fₚ.p

# ╔═╡ dcc30d6f-3170-4e6d-ad51-741bce871a42
# Matrica pivotiranja
Fₚ.P

# ╔═╡ ad9107ee-2abd-4dd7-b821-2e0a4bae4374
# Provjera s matricom
Fₚ.Q*Fₚ.R-A*Fₚ.P

# ╔═╡ 17cebd37-73b8-4699-a119-4237ab85d084
# Provjera s vektorom
Fₚ.Q*Fₚ.R-A[:,Fₚ.p]

# ╔═╡ a199bc51-275d-44a6-bc48-e97cb533922f
F.factors

# ╔═╡ 550adbc6-118d-4dab-b13b-5cc2eb59b287
F.R

# ╔═╡ 4c803427-419f-42ee-8672-1af245889d33
v1=HouseholderVector(A[:,1])

# ╔═╡ 3c9bf5b9-3a35-46e8-bcad-a9b16b0e731d
v1/v1[1]

# ╔═╡ 60166ed0-dfe1-474a-91df-64c433867489
md"""
## Brzina

Broj računskih operacija potrebnih za računanje QR rastava matrice $n\times n$ je $O\big(\frac{4}{3}n^3\big)$ za računanje matrice $R$ i  $O\big(\frac{4}{3}n^3\big)$ za računanje matrice $Q$. 

"""

# ╔═╡ 7a2fb72c-2c5c-4e1a-895a-fb0eb5608876
begin
	n=512
	A₁=rand(n,n);
end

# ╔═╡ c7f654da-3ad0-4781-a1b0-f7f31a19bf74
# Ispis je u Julia terminalu
@time F₁=qr(A₁);

# ╔═╡ 8a76f428-69aa-4770-a5f5-766e47fa9f9c
@time lu(A₁);

# ╔═╡ c0eaba40-c919-49af-96f8-700173e3a72a
@time qr(A₁,Val(true));

# ╔═╡ f6b568ec-20e5-4575-a62c-70d470589ef1
@time HouseholderQR(A₁);

# ╔═╡ 7e40ad34-0497-410d-9669-20ad8a9ca74b
md"""
## Točnost

Za matrice $\hat Q$ i $\hat R$ izračunate Householderovom metodom vrijedi: 

$$
\begin{aligned}
\hat Q^T\hat Q& =I+E, \qquad \|E \|_2\approx \varepsilon,\\ 
\| A-\hat Q\hat R\|_2& \approx \varepsilon\|A\|_2.
\end{aligned}$$

Također, postoji egzaktna ortogonalna matrica $Q$ za koju je 

$$
\| A- Q\hat R\|_2\approx \varepsilon\|A\|_2.$$
"""

# ╔═╡ 652da960-3488-11eb-0414-0b6f313755a5
Q₁=Matrix(F₁.Q)

# ╔═╡ e6c3a470-3488-11eb-0387-053b5428d074
Q₁'*Q₁

# ╔═╡ 03251d60-3489-11eb-2754-57656ea9b972
norm(A₁)

# ╔═╡ 0cd0b360-3489-11eb-3cd8-b9f2d4831da2
norm(A₁-F₁.Q*F₁.R)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
PlutoUI = "~0.7.21"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "abb72771fd8895a7ebd83d5632dc4b989b022b5b"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "ae4bbcadb2906ccc085cf52ac286dc1377dceccc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "b68904528fd538f1cb6a3fbc44d2abdc498f9e8e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.21"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═66ac2c80-5d8d-4427-a002-14f14bdc2535
# ╠═52742d2f-73d7-4b5c-8f48-69f21d345133
# ╟─fb66b990-4511-476f-8e77-87aa9c265041
# ╠═24c9f0ba-6e18-4362-840b-1a3581032ff5
# ╠═b6ee2c00-78e7-4d22-a851-d22c4a4946bc
# ╠═b8994eec-d651-4921-b1a3-184a029e83b7
# ╠═19a0f0b8-dbab-4609-9f7f-17339af6bf63
# ╠═54f456e7-133d-4f96-901f-5bcda000178f
# ╠═1d11a929-212b-42d0-ab52-c793ed8905d2
# ╠═9038e960-55c0-437a-a9fd-73f48860a36a
# ╟─b9bf7235-dcb8-47bf-8c3f-1179ba9e08e0
# ╟─3fd60d80-25d7-45e2-9776-7250841e8f74
# ╠═4a371793-1d5c-413f-8e89-dce986c972ef
# ╠═814c6dd6-faa1-4211-8f64-bf349624dc37
# ╠═e282bbd3-6505-4afc-b304-3c327fa5db62
# ╠═eb8f508a-9018-4896-bd82-e707e01157bc
# ╟─df0bb186-4b2c-429a-b279-5843f97cd3f0
# ╠═3bb4d086-a79a-44b1-92c3-5e29be4c36de
# ╠═d7ddbafd-8923-45c1-a0f2-55aa1dc86185
# ╠═cb302477-773e-4fd9-a9d6-a9df50e45b18
# ╠═12069de5-8064-437b-b96c-222baa7b4c7d
# ╠═3ceae442-3a07-4400-933c-a04bdc8ee742
# ╟─ef4702ee-44cc-49d8-85bb-39249c431036
# ╠═4587b182-80fd-4bd7-87af-17074ea7c9d9
# ╠═f944cdbe-794c-405e-a215-687cbb12b22b
# ╠═64e052c2-b798-4fb8-80fd-c018df2c8625
# ╠═ef2ebf49-e7b2-4300-9a5c-3a2e42251e74
# ╠═6bf029d7-68e2-47a4-87b9-8f66154933d4
# ╠═68d0e359-9fb9-4a9b-9f13-e8039597e948
# ╠═dcc30d6f-3170-4e6d-ad51-741bce871a42
# ╠═ad9107ee-2abd-4dd7-b821-2e0a4bae4374
# ╠═17cebd37-73b8-4699-a119-4237ab85d084
# ╠═a199bc51-275d-44a6-bc48-e97cb533922f
# ╠═550adbc6-118d-4dab-b13b-5cc2eb59b287
# ╠═4c803427-419f-42ee-8672-1af245889d33
# ╠═3c9bf5b9-3a35-46e8-bcad-a9b16b0e731d
# ╟─60166ed0-dfe1-474a-91df-64c433867489
# ╠═7a2fb72c-2c5c-4e1a-895a-fb0eb5608876
# ╠═c7f654da-3ad0-4781-a1b0-f7f31a19bf74
# ╠═8a76f428-69aa-4770-a5f5-766e47fa9f9c
# ╠═c0eaba40-c919-49af-96f8-700173e3a72a
# ╠═f6b568ec-20e5-4575-a62c-70d470589ef1
# ╟─7e40ad34-0497-410d-9669-20ad8a9ca74b
# ╠═652da960-3488-11eb-0414-0b6f313755a5
# ╠═e6c3a470-3488-11eb-0387-053b5428d074
# ╠═03251d60-3489-11eb-2754-57656ea9b972
# ╠═0cd0b360-3489-11eb-3cd8-b9f2d4831da2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
