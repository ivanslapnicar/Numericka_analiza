### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ decffbb4-4f7d-4b49-9d2e-c91b4b00587c
using PlutoUI, LinearAlgebra, Random

# ╔═╡ 59014114-8bb5-4fa5-96c7-d51551d61ba1
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ fcc44b72-e162-4351-8601-f7402e2ed694
md"""
# Metoda najmanjih kvadrata


Neka je zadan sustav s više jednadžbi od nepoznanica:

$$Ax=b, \quad m>n.$$

Ako sustav ima rješenje, tada je je $Ax-b=0$, odnosno $\| Ax-b\|=0$ za svaku vektorsku  normu.

Ako sustav nema rješenje, tada je prirodno tražiti rješenje za koje je 

$$
\|Ax-b \|_{1,2,\infty}\to \min$$

za odabranu vektorsku normu.
"""

# ╔═╡ a062b872-1eaa-11eb-005f-9d66fad5ee28
md"
__Teorem.__ Ako je $\mathop{\mathrm{rang}} A=n$, tada se __jedinstveno__ rješenje $x$ za koje 

$$
\|Ax-b \|_{2}\to \min$$

dobije rješavanjem sustava __normalnih jednadžbi__:

$$
A^T A x=A^T b. \tag{*}$$
"

# ╔═╡ a6fc5380-1eaa-11eb-11a9-3544f456255c
md"
_Dokaz_ : Definirajmo

$$
Q(x)=\|Ax-b\|_2^2=(x^TA^T-b^T)(Ax-b)=x^TA^T A x -2x^T A^T b+b^Tb.$$

Vrijedi

$$\begin{aligned}
Q(x+h)&=(x^T+h^T)A^TA(x+h)-2(x^T+h^T)A^Tb+b^Tb \\
&=Q(x) +2h^T(A^TAx-A^Tb)+h^TA^TAh\\ &= Q(x)+\|Ah\|_2^2 \\
&\geq Q(x),
\end{aligned}$$

pa se minimum zaista postiže u $x$.

Rješenje je jedinstveno jer $Q(x)=Q(y)$ povlači $\|Ah\|_2=0$ pa je ili $h=0$ ili $\mathop{\mathrm{rang}} A<n$ što je kontradikcija i teorem je dokazan.
"

# ╔═╡ b4ddcb00-1eaa-11eb-23c6-d15643cd207a
md"
__Geometrijsko značenje.__ Vektori $Ax$ i $Ax -b$ su međusobno okomiti, 

$$
(Ax)^T\cdot (Ax - b)=x^T (A^TAx - A^Tb)=0.$$ 

Dakle, $Ax$ je ortogonalna projekcija vektora $b$ na skup $\{Ay:\ y \textrm{ proizvoljan}\}$.

Rješenje $x$ zove se __kvadratična prilagodba__ 
sustavu $A x=b$ u smislu najmanjih kvadrata. __Kvalitetu prilagodbe__ mjerimo s

$$
q=\sqrt{\frac{Q(x)}{Q(0)}}=\frac{\|A x - b\|_2}{\|b\|_2 }.$$
"

# ╔═╡ 5c65704d-666f-4f15-bc8f-7741457f9af0
md"""
## Primjer 1

Riješimo sustav 


$$\begin{aligned}
x+y&=0\\
y+z&=1\\
x+z&=0\\
-x+y+z&=1\\
-x-z&=0
\end{aligned}$$

u smislu najmanjih kvadrata.
"""

# ╔═╡ 365fc919-988d-4a1b-b42c-b8ab6931f860
A=[1//1 1 0;0 1 1;1 0 1;-1 1 1;-1 0 -1]

# ╔═╡ dea071d2-66b7-44f0-a75c-e0c67e574561
b=[0//1,1,0,1,0]

# ╔═╡ e6753e70-2402-11eb-3055-594fe9e7c50c
A'*A

# ╔═╡ eedd09d0-2402-11eb-219e-2dfe890d6b22
A'*b

# ╔═╡ d6e05d5d-6878-4865-934a-6d8846f1d157
x=(A'*A)\(A'*b)

# ╔═╡ 49874317-1bb6-4882-ae83-4644094bf87e
# Kvaliteta prilagodbe
q=√(norm(A*x-b)/norm(b))

# ╔═╡ 1052a0b7-748d-46d5-a5fb-9d3b1ba2b65e
md"
Ako je sustav predefiniran, standardna naredba `\` odmah računa kvadratičnu prilagodbu, pri čemu se koristi QR rastav.
"

# ╔═╡ a0dbeff0-8d5f-11eb-0c1a-bb82ced59a92
A\b

# ╔═╡ d8458b22-f8b3-4fda-8351-3abe5af8dc46
float(x)

# ╔═╡ c2f32713-6b34-4c34-9180-d759039891c5
md"""
## Primjer 2
"""

# ╔═╡ d92d5bed-8689-481f-b56f-d46fa4f835c1
begin
	Random.seed!(123)
	A₁=rand(20,10)
	b₁=rand(20);
end

# ╔═╡ 74feeef4-02b4-454c-9d86-c29d775c89c0
x₁=A₁\b₁

# ╔═╡ 73d14810-4471-40ad-b1ca-f1d210ad2eb2
q₁=√(norm(A₁*x₁-b₁)/norm(b₁))

# ╔═╡ f7f6d6b8-44cc-4337-a601-3c4c5ef3fb77
md"""
# Teorija smetnje

__Osjetljivost problema najmanjih kvadarata__ dana je sljedećim ocjenama (vidi [Matrix Computations, poglavlje 5, str. 266](https://books.google.hr/books?id=X5YfsuCWpxMC&printsec=frontcover&hl=hr#v=onepage&q&f=false)).

Za matricu $A$ __kondiciju__ definiramo na sljedeći način:

$$
\kappa_2(A)=\sqrt{\kappa_2(A^TA)}=\|A\|_2 \|(A^TA)^{-1} A^T\|_2.$$

Neka su $x$ i $\hat x$, kvadratične prilagodbe sustava $Ax=b$ i 
$(A+\delta A)\hat x=b+\delta b$. __Reziduali__ su definirani s

$$
\begin{aligned}
r&=Ax-b\\
\hat r&=(A+\delta A)\hat x-(b+\delta b).
\end{aligned}$$

Neka je 

$$
\epsilon=\max \bigg\{ \frac{\|\delta A\|_2}{\|A\|_2},\frac{\|\delta b\|_2}{\|b\|_2}\bigg\}$$

i neka je 

$$
q=\frac {\|r\|_2}{\|b\|_2}\equiv\sin\theta <1.$$

Vrijedi:

$$
\begin{aligned}
\frac{\|\hat x-x\|_2}{\|x\|_2}&\leq \epsilon \bigg[\frac{2\,\kappa_2(A)}{\cos \theta} +\tan\theta \,\kappa_2^2(A)\bigg]+O(\epsilon^2),\\
\frac{\|\hat r-r\|_2}{\|b\|_2}&\leq \epsilon\,[1+ 2\,\kappa_2(A)](m-n)+O(\epsilon^2).
\end{aligned}$$

Vidimo da je rezidual manje osjetljiv od samog mjesta na kojem se postiže.
"""

# ╔═╡ fbcad473-e073-4edf-a46f-f1d28a3b753d
cond(A₁)

# ╔═╡ 4dcccfd8-fdf6-4cdc-8314-c2a0d002baff
δA₁=1e-4*(rand(20,10).-0.5)

# ╔═╡ bb45e359-ca50-49cd-a23c-16552fac65ed
xp₁=(A₁+δA₁)\b₁

# ╔═╡ 4bd4192c-4d5e-4869-b598-9edd04128c7a
begin
	r₁=A₁*x₁-b₁
	rp₁=(A₁+δA₁)*xp₁-b₁
end

# ╔═╡ d7a8a461-17f6-4a9e-9343-f58314e6a9ee
norm(xp₁-x₁)/norm(x₁), norm(rp₁-r₁)/norm(b₁)

# ╔═╡ 314c1246-1772-415f-9f9a-38b221b6eb96
md"""
# Analiza greške i točnost

Ako je $\mathop{\mathrm{rang}}A =n$, matrica $A^TA$ je simetrična i pozitivno definitna pa se sustav (*) može riješiti metodom Choleskog.

Za izračunato rješenje $\hat x$ vrijedi

$$
(A^TA +E)\hat x=A^Tb,$$

gdje je 

$$
\|E\|_2\approx \varepsilon \| A^TA\|_2,$$

pa za relativnu pogrešku vrijedi ocjena

$$
\frac{\|\hat x -x\|_2}{\|x\|_2}\approx \varepsilon \kappa_2(A^TA) =\varepsilon \kappa^2_2(A).$$


Dakle, relativna pogreška rješenja dobivenog pomoću metode normalnih jednadžbi ovisi o __kvadratu kondicije__ pa je bolje koristiti QR rastav.
"""

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
# ╠═decffbb4-4f7d-4b49-9d2e-c91b4b00587c
# ╠═59014114-8bb5-4fa5-96c7-d51551d61ba1
# ╟─fcc44b72-e162-4351-8601-f7402e2ed694
# ╟─a062b872-1eaa-11eb-005f-9d66fad5ee28
# ╟─a6fc5380-1eaa-11eb-11a9-3544f456255c
# ╟─b4ddcb00-1eaa-11eb-23c6-d15643cd207a
# ╟─5c65704d-666f-4f15-bc8f-7741457f9af0
# ╠═365fc919-988d-4a1b-b42c-b8ab6931f860
# ╠═dea071d2-66b7-44f0-a75c-e0c67e574561
# ╠═e6753e70-2402-11eb-3055-594fe9e7c50c
# ╠═eedd09d0-2402-11eb-219e-2dfe890d6b22
# ╠═d6e05d5d-6878-4865-934a-6d8846f1d157
# ╠═49874317-1bb6-4882-ae83-4644094bf87e
# ╟─1052a0b7-748d-46d5-a5fb-9d3b1ba2b65e
# ╠═a0dbeff0-8d5f-11eb-0c1a-bb82ced59a92
# ╠═d8458b22-f8b3-4fda-8351-3abe5af8dc46
# ╟─c2f32713-6b34-4c34-9180-d759039891c5
# ╠═d92d5bed-8689-481f-b56f-d46fa4f835c1
# ╠═74feeef4-02b4-454c-9d86-c29d775c89c0
# ╠═73d14810-4471-40ad-b1ca-f1d210ad2eb2
# ╟─f7f6d6b8-44cc-4337-a601-3c4c5ef3fb77
# ╠═fbcad473-e073-4edf-a46f-f1d28a3b753d
# ╠═4dcccfd8-fdf6-4cdc-8314-c2a0d002baff
# ╠═bb45e359-ca50-49cd-a23c-16552fac65ed
# ╠═4bd4192c-4d5e-4869-b598-9edd04128c7a
# ╠═d7a8a461-17f6-4a9e-9343-f58314e6a9ee
# ╟─314c1246-1772-415f-9f9a-38b221b6eb96
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
