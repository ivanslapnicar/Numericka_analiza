### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 08b094da-0680-4727-bdc4-83f71de41b32
using Polynomials, Random

# ╔═╡ 70d60c9a-ff90-4003-a2e9-b3f0b5c88fb2
md"""
# Izvrednjavanje funkcija


Računalo može izvoditi samo četiri osnovne operacije, `+`, `-`, `*` i `/` pa se sve ostale funkcije 
računaju pomoću polinoma (npr. Taylorova formula uz ocjenu ostatka ili bolje formule).

Neka je zadan polinom __stupnja__ $n$:

$$
p_n(x)=a_0+a_1 x+a_2x^2+a_3 x^3+\cdots + a_{n-1}x^{n-1}+a_n x^n,\quad a_n\neq 0.$$

## Brzina

Direktno računanje vrijednosti $p_n(x)$ treba $O(n^2)$ operacija.

Uz __pamćenje potencija__ imamo sljedeći algoritam:
"""

# ╔═╡ 5c6e0ac4-9b1e-498b-a2f5-994e2a596783
function mypolyval(p::Polynomial,x::Number)
    s=p[0]
    t=one(p[0])
    for i=1:length(p)-1
        t*=x
        s+=p[i]*t
    end
    s
end 

# ╔═╡ d40a5115-1f38-4e4a-9bde-37408f90a5f5
p=Polynomial([1,2,3,4,5])

# ╔═╡ 58d8436a-740c-4801-aba3-92d01070a4e6
mypolyval(p,3)

# ╔═╡ 834e9612-5b50-4924-a158-1effe9046a74
mypolyval(p,π)

# ╔═╡ fda3e34e-1901-11eb-126a-795d469edfcc
p(3), p(π)

# ╔═╡ 2b402caa-7cbb-4402-8307-0a36af1410c5
md"""
Funkcija `mypolyval()` koristi $2n$ množenja i $n$ zbrajanja.
"""

# ╔═╡ 5c6e1617-6682-43ec-bec5-8245f1a55d34
Random.seed!(123)

# ╔═╡ 6da79660-153d-11eb-10f7-b9e079919759
pbig=Polynomial(randn(1000));

# ╔═╡ 923c54ff-bf4b-42df-96ae-96b22c728aa0
@time mypolyval(pbig,1.5)

# ╔═╡ e4ac8c0c-f8ec-4b95-9ca7-f8aa1b3dbdfe
@time pbig(1.5)

# ╔═╡ 29ac5eb4-89f4-4a9e-a548-b087e6b8da4f
md"
Let's see how is this done:
"

# ╔═╡ f85661e2-25c1-4174-ab29-cbc19ccf70d9
@which pbig(1.5)

# ╔═╡ 529715c6-48ed-445b-bc8e-3337b43e16c3
@which evalpoly(1.5,pbig)

# ╔═╡ 867f44a5-590d-4011-9133-73cb99f41481
@which Base.evalpoly(1.5,pbig.coeffs)

# ╔═╡ e3c1f285-fd03-47f5-b705-914e4cb35170
md"""
__Hornerova shema__ (Horner 1819, Newton 1669) treba $n$ množenja i $n$ zbrajanja:

$${\displaystyle {\begin{aligned}a_{0}&+a_{1}x+a_{2}x^{2}+a_{3}x^{3}+\cdots +a_{n}x^{n}\\&=a_{0}+x{\bigg (}a_{1}+x{\Big (}a_{2}+x{\big (}a_{3}+\cdots +x(a_{n-1}+x\,a_{n})\cdots {\big )}{\Big )}{\bigg )}\,,\end{aligned}}}$$

"""

# ╔═╡ b3718874-ea9f-4f5d-ad83-d22cb272b489
function myhorner(p::Polynomial,x::Number)
    s=p[end]
    for i=length(p)-2:-1:0
        # s*=x
        # s+=p[i]
        s=s*x+p[i]
    end
    s
end

# ╔═╡ 24d7653c-c02b-4239-ace6-b6b99ffea022
myhorner(p,3)

# ╔═╡ e3382c91-e26f-4636-b068-483dab80a3bc
@time myhorner(pbig,1.5)

# ╔═╡ cce27db4-32cd-4ff1-8531-880985527ff2
md"""
Hornerova shema je __optimalna__ u smislu da je općenito za izvrednjavanje polinoma $p_n(x)$ potrebno barem $n$ množenja. 

Mogući su, naravno, posebni slučajevi, kao $x^{100}$. 
"""

# ╔═╡ eb155720-1904-11eb-0006-a31743a0a7f2
log2(100)

# ╔═╡ a7ab8371-c193-431d-9b1d-d3269a44459d
md"""
## Točnost

Neka je $\hat q$ vrijednost $p_n(x)$ izračunata u aritmetici s točnošću stroja $\varepsilon$. Tada vrijedi ocjena
(vidi [Accuracy and Stability of Numerical Algorithms, str. 95](https://books.google.hr/books?id=5tv3HdF-0N8C&printsec=frontcover&hl=hr#v=onepage&q&f=false)):

$$
\big|\, p_n(x)-\hat q\,\big| \leq \frac{2n\varepsilon}{1-2n\varepsilon} \sum_{i=0}^n |a_i||x|^i.$$
"""

# ╔═╡ b5ffd702-3c09-4b37-a58f-b376826fbbdb
begin
	p₁=Polynomial([1,2,3,4,5])
	myhorner(p₁,sqrt(2))
end

# ╔═╡ 8410e079-7970-48a2-abd9-38b465547b90
BigFloat(sqrt(2))

# ╔═╡ 45b1021e-e34d-43fe-8b47-182ec35adaa9
sqrt(BigFloat(2))

# ╔═╡ 4700a424-1f5e-4124-b851-c196d45dd8cc
pb₁=Polynomial([BigInt(1),2,3,4,5])

# ╔═╡ 167c82a2-014b-4ff6-b5fe-036096f00007
myhorner(pb₁,sqrt(BigFloat(2)))

# ╔═╡ ca410bf6-18b3-4d97-93a6-2fa6b42de453
myhorner(p₁,sqrt(200000))

# ╔═╡ fa022c78-1d6a-4490-bee9-feff0f7dbba1
myhorner(pb₁,sqrt(BigFloat(200000)))

# ╔═╡ 7b15e5da-8318-4b1c-9e72-07746cd8f2fe
begin
	r=[1,sqrt(2),3,4,5,6,sqrt(50)]
	p₂=fromroots(r)
end

# ╔═╡ 412c6ac6-5fc3-42a7-97a8-21d074a3b6b6
pb₂=fromroots(map(BigFloat,r))

# ╔═╡ 3f4613bb-42f4-4df8-90bc-3b134ce9910c
myhorner(p₂,sqrt(2)+0.1)

# ╔═╡ ce0209f9-8f95-42a6-b0ad-8ad4025ba74d
myhorner(pb₂,sqrt(BigFloat(2))+0.1)

# ╔═╡ 9bc8612e-84ce-40e0-bb76-6802b2957918
myhorner(p₂,-sqrt(10000))

# ╔═╡ 0d802709-35ac-424e-82e1-6d355dde5503
myhorner(pb₂,-sqrt(BigFloat(10000)))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Polynomials = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
Polynomials = "~2.0.15"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "77b022fd16313d405285aeb9a0c1598467d99e2c"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "8a62af3e248a8c4bad6b32cbbe663ae02275e32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "ac0aaa807ed5eaf13f67afe188ebc07e828ff640"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.10.0"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "4cc0c5a83933648b615c36c2b956d94fda70641e"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "6985021d02ab8c509c841bb8b2becd3145a7b490"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.3.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "716e24b21538abc91f6205fd1d8363f39b442851"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Polynomials]]
deps = ["Intervals", "LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "a1f7f4e41404bed760213ca01d7f384319f717a5"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "2.0.25"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TZJData]]
deps = ["Artifacts"]
git-tree-sha1 = "d39314cdbaf5b90a047db33858626f8d1cc973e1"
uuid = "dc5dba14-91b3-4cab-a142-028a31da12f7"
version = "1.0.0+2023c"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TimeZones]]
deps = ["Artifacts", "Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "Scratch", "TZJData", "Unicode", "p7zip_jll"]
git-tree-sha1 = "89e64d61ef3cd9e80f7fc12b7d13db2d75a23c03"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.13.0"
weakdeps = ["RecipesBase"]

    [deps.TimeZones.extensions]
    TimeZonesRecipesBaseExt = "RecipesBase"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─70d60c9a-ff90-4003-a2e9-b3f0b5c88fb2
# ╠═08b094da-0680-4727-bdc4-83f71de41b32
# ╠═5c6e0ac4-9b1e-498b-a2f5-994e2a596783
# ╠═d40a5115-1f38-4e4a-9bde-37408f90a5f5
# ╠═58d8436a-740c-4801-aba3-92d01070a4e6
# ╠═834e9612-5b50-4924-a158-1effe9046a74
# ╠═fda3e34e-1901-11eb-126a-795d469edfcc
# ╟─2b402caa-7cbb-4402-8307-0a36af1410c5
# ╠═5c6e1617-6682-43ec-bec5-8245f1a55d34
# ╠═6da79660-153d-11eb-10f7-b9e079919759
# ╠═923c54ff-bf4b-42df-96ae-96b22c728aa0
# ╠═e4ac8c0c-f8ec-4b95-9ca7-f8aa1b3dbdfe
# ╟─29ac5eb4-89f4-4a9e-a548-b087e6b8da4f
# ╠═f85661e2-25c1-4174-ab29-cbc19ccf70d9
# ╠═529715c6-48ed-445b-bc8e-3337b43e16c3
# ╠═867f44a5-590d-4011-9133-73cb99f41481
# ╟─e3c1f285-fd03-47f5-b705-914e4cb35170
# ╠═b3718874-ea9f-4f5d-ad83-d22cb272b489
# ╠═24d7653c-c02b-4239-ace6-b6b99ffea022
# ╠═e3382c91-e26f-4636-b068-483dab80a3bc
# ╟─cce27db4-32cd-4ff1-8531-880985527ff2
# ╠═eb155720-1904-11eb-0006-a31743a0a7f2
# ╟─a7ab8371-c193-431d-9b1d-d3269a44459d
# ╠═b5ffd702-3c09-4b37-a58f-b376826fbbdb
# ╠═8410e079-7970-48a2-abd9-38b465547b90
# ╠═45b1021e-e34d-43fe-8b47-182ec35adaa9
# ╠═4700a424-1f5e-4124-b851-c196d45dd8cc
# ╠═167c82a2-014b-4ff6-b5fe-036096f00007
# ╠═ca410bf6-18b3-4d97-93a6-2fa6b42de453
# ╠═fa022c78-1d6a-4490-bee9-feff0f7dbba1
# ╠═7b15e5da-8318-4b1c-9e72-07746cd8f2fe
# ╠═412c6ac6-5fc3-42a7-97a8-21d074a3b6b6
# ╠═3f4613bb-42f4-4df8-90bc-3b134ce9910c
# ╠═ce0209f9-8f95-42a6-b0ad-8ad4025ba74d
# ╠═9bc8612e-84ce-40e0-bb76-6802b2957918
# ╠═0d802709-35ac-424e-82e1-6d355dde5503
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
