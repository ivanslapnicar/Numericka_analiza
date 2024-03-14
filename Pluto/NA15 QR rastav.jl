### A Pluto.jl notebook ###
# v0.19.40

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
	Random.seed!(125)
	A=randn(8,5)
	# A[:,3]=A[:,4]+1e-5*rand(8)
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
norm(A-Q*R)

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

# ╔═╡ 5a42e1e7-7704-42d1-bbca-2eb09cf80cbc
a=[5.0;1]

# ╔═╡ 7040f957-6000-4be5-8fa9-314a6140c0d8
ϕ=0.3

# ╔═╡ 0c8b4d0e-62c0-4ee5-999f-274ab410b305
g=[cos(ϕ) sin(ϕ);sin(ϕ) -cos(ϕ)]

# ╔═╡ 26ab8ae1-3118-41aa-bb55-4f122dd62150
g*a

# ╔═╡ 4a371793-1d5c-413f-8e89-dce986c972ef
function HouseholderVector(x::Vector)
    # Računa v
    v=copy(x)
    v[1]=x[1]+sign(x[1])*norm(x)
    v
end

# ╔═╡ ef6f8af2-4ab8-4958-be52-3302448ff887
vr=HouseholderVector(a)

# ╔═╡ d391ffc6-974f-4a27-b0a6-a556392f4365
Hr=I-2vr*vr'/(vr'*vr)

# ╔═╡ 3b6cd3d5-fbf4-468f-8bbe-1149df7629f2
Hr*a

# ╔═╡ 814c6dd6-faa1-4211-8f64-bf349624dc37
begin
	x=randn(8)
	v=HouseholderVector(x)
	β=(2/(v⋅v))*(v⋅x)
	x-β*v
end

# ╔═╡ 7a7aa5c9-2dcd-45cc-8ac3-569c7a55b6de
β

# ╔═╡ e282bbd3-6505-4afc-b304-3c327fa5db62
x

# ╔═╡ eb8f508a-9018-4896-bd82-e707e01157bc
norm(x)

# ╔═╡ df0bb186-4b2c-429a-b279-5843f97cd3f0
md"""
QR rastav matrice se računa rekurzivnim QR rastavom vektora pomoću Householderovih reflektora:
"""

# ╔═╡ 3bb4d086-a79a-44b1-92c3-5e29be4c36de
function HouseholderQR(A₁::Matrix{T}) where T
    # Računa Q i R
    A=copy(A₁)
    m,n=size(A)
    Q=Matrix{T}(I,m,m) # eye
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

# ╔═╡ ac488c01-8acb-4b5d-9251-1e6c6574997e
begin
	# Probajmo za kompleksne matrice
	B=randn(ComplexF64,6,3)
	Qb,Rb=HouseholderQR(B)
	norm(Qb'*Qb-I), norm(Qb*Rb-B)
end

# ╔═╡ 8347726d-721a-46ad-b7af-332764ddae11
Rb

# ╔═╡ ef4702ee-44cc-49d8-85bb-39249c431036
md"""
Program `HouseholderQR()` je ilustrativan. Profesionalni programi imaju sljedeća svojstva:

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

# ╔═╡ 90ad6140-fce2-457f-98fc-aecccab86589
Matrix(F.Q)

# ╔═╡ c756e15a-dab1-4846-b0f3-fb232956ff21
F.R

# ╔═╡ 64e052c2-b798-4fb8-80fd-c018df2c8625
F.Q'*A

# ╔═╡ ef2ebf49-e7b2-4300-9a5c-3a2e42251e74
norm(F.Q*F.R-A)

# ╔═╡ 6bf029d7-68e2-47a4-87b9-8f66154933d4
Fₚ=qr(A,Val(true))

# ╔═╡ 3406770c-943e-4d44-8e42-d3533621cf22
A

# ╔═╡ 68d0e359-9fb9-4a9b-9f13-e8039597e948
# Vektor pivotiranja
Fₚ.p

# ╔═╡ dcc30d6f-3170-4e6d-ad51-741bce871a42
# Matrica pivotiranja
Fₚ.P

# ╔═╡ ad9107ee-2abd-4dd7-b821-2e0a4bae4374
# Provjera s matricom
norm(Fₚ.Q*Fₚ.R-A*Fₚ.P)

# ╔═╡ 17cebd37-73b8-4699-a119-4237ab85d084
# Provjera s vektorom
norm(Fₚ.Q*Fₚ.R-A[:,Fₚ.p])

# ╔═╡ a199bc51-275d-44a6-bc48-e97cb533922f
F.factors

# ╔═╡ 4c803427-419f-42ee-8672-1af245889d33
v1=HouseholderVector(A[:,1])

# ╔═╡ 3c9bf5b9-3a35-46e8-bcad-a9b16b0e731d
v1/v1[1]

# ╔═╡ 550adbc6-118d-4dab-b13b-5cc2eb59b287
F.R

# ╔═╡ 60166ed0-dfe1-474a-91df-64c433867489
md"""
## Brzina

Broj računskih operacija potrebnih za računanje QR rastava matrice $n\times n$ je $O\big(\frac{4}{3}n^3\big)$ za računanje matrice $R$ i  $O\big(\frac{4}{3}n^3\big)$ za računanje matrice $Q$. 

"""

# ╔═╡ 7a2fb72c-2c5c-4e1a-895a-fb0eb5608876
begin
	n=512
	A₁=randn(n,n);
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
PlutoUI = "~0.7.58"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.2"
manifest_format = "2.0"
project_hash = "1867d9ce1bd88115b124f124b5d7cd866c186b11"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "0f748c81756f2e5e6854298f11ad8b2dfae6911a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "71a22244e352aa8c5f0f2adde4150f62368a3f2e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.58"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
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
# ╠═5a42e1e7-7704-42d1-bbca-2eb09cf80cbc
# ╠═7040f957-6000-4be5-8fa9-314a6140c0d8
# ╠═0c8b4d0e-62c0-4ee5-999f-274ab410b305
# ╠═26ab8ae1-3118-41aa-bb55-4f122dd62150
# ╠═ef6f8af2-4ab8-4958-be52-3302448ff887
# ╠═d391ffc6-974f-4a27-b0a6-a556392f4365
# ╠═3b6cd3d5-fbf4-468f-8bbe-1149df7629f2
# ╠═4a371793-1d5c-413f-8e89-dce986c972ef
# ╠═814c6dd6-faa1-4211-8f64-bf349624dc37
# ╠═7a7aa5c9-2dcd-45cc-8ac3-569c7a55b6de
# ╠═e282bbd3-6505-4afc-b304-3c327fa5db62
# ╠═eb8f508a-9018-4896-bd82-e707e01157bc
# ╟─df0bb186-4b2c-429a-b279-5843f97cd3f0
# ╠═3bb4d086-a79a-44b1-92c3-5e29be4c36de
# ╠═d7ddbafd-8923-45c1-a0f2-55aa1dc86185
# ╠═cb302477-773e-4fd9-a9d6-a9df50e45b18
# ╠═12069de5-8064-437b-b96c-222baa7b4c7d
# ╠═3ceae442-3a07-4400-933c-a04bdc8ee742
# ╠═ac488c01-8acb-4b5d-9251-1e6c6574997e
# ╠═8347726d-721a-46ad-b7af-332764ddae11
# ╟─ef4702ee-44cc-49d8-85bb-39249c431036
# ╠═4587b182-80fd-4bd7-87af-17074ea7c9d9
# ╠═90ad6140-fce2-457f-98fc-aecccab86589
# ╠═f944cdbe-794c-405e-a215-687cbb12b22b
# ╠═c756e15a-dab1-4846-b0f3-fb232956ff21
# ╠═64e052c2-b798-4fb8-80fd-c018df2c8625
# ╠═ef2ebf49-e7b2-4300-9a5c-3a2e42251e74
# ╠═6bf029d7-68e2-47a4-87b9-8f66154933d4
# ╠═3406770c-943e-4d44-8e42-d3533621cf22
# ╠═68d0e359-9fb9-4a9b-9f13-e8039597e948
# ╠═dcc30d6f-3170-4e6d-ad51-741bce871a42
# ╠═ad9107ee-2abd-4dd7-b821-2e0a4bae4374
# ╠═17cebd37-73b8-4699-a119-4237ab85d084
# ╠═a199bc51-275d-44a6-bc48-e97cb533922f
# ╠═4c803427-419f-42ee-8672-1af245889d33
# ╠═3c9bf5b9-3a35-46e8-bcad-a9b16b0e731d
# ╠═550adbc6-118d-4dab-b13b-5cc2eb59b287
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
