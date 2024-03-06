### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 379d9b77-5d2f-4477-81ab-1957b75f13e1
using LinearAlgebra, Random

# ╔═╡ f0e54010-5a21-11ed-12bc-45182a10006b
md"
# Procjena uvjetovanosti
"

# ╔═╡ ac6bf8ee-5d9f-4216-a71e-973b63fbc320
Random.seed!(431)

# ╔═╡ c58cf6e0-f332-44c1-945d-c162bb2465f7
n=1024

# ╔═╡ f9769a5c-6468-41ed-bd1b-0e1ae48b6a27
A=rand(n,n);

# ╔═╡ 64c7043a-9a38-41a6-b79c-a8604fa18b91
b=rand(n);

# ╔═╡ 2872a885-da17-4fd3-a3f4-e32a3871adc2
@time lu(A);

# ╔═╡ e1818076-86d5-48df-a448-b517e6cace30
@time A\b;

# ╔═╡ e846aa72-5b03-4af3-b8f9-54da7411ff98
# Ovo traje previše dugo
@time cond(A)

# ╔═╡ 908011e2-e506-496c-944e-3a346398f2f6
A[1:3,1:3]

# ╔═╡ 7cb17280-2003-4825-90b7-50ef837b54e5
b[1:3]

# ╔═╡ 4b3ce556-34b9-4308-8b25-b054995c1ba4
@which cond(A)

# ╔═╡ 9361923a-6019-44ee-811f-68d217b90aef
md"
`cond()` računa uvjetovanost pomoću rastava singularnih vrijednosti (SVD).

Umjesto toga, sustav se može riješiti pomoću LAPack-ove ekspertne funkcije [__gesvx()__](https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/lapack-routines/lapack-linear-equation-routines/lapack-linear-equation-driver-routines/gesvx.html#gesvx), koja uz malo dodatnog vremena računa procjenu inverzne kondicije u 1-normi. Za proračun kondicije koristi se izračunati LU rastav i funkcija [__gecon()__](https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/lapack-routines/lapack-linear-equation-routines/lapack-linear-equation-computational-routines/estimate-the-condition-number-lapack-computation/gecon.html#gecon).

Ekvivalentno, možemo računati kondiciju u 1-normi, `cond(A,1)`.
"

# ╔═╡ c83d84e4-4f4d-4a97-b50f-1a5805e54c32
@time x=LAPACK.gesvx!(A,b)

# ╔═╡ 65be9c7e-c4cc-4c03-8537-651f648c2f44
1/x[2]

# ╔═╡ d611d92b-872e-4cfe-8864-1e706f7c6fd6
@time cond(A,1)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.1"
manifest_format = "2.0"
project_hash = "f479011a250b0799ce99df8978eea5d0a8ab069c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.0+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"
"""

# ╔═╡ Cell order:
# ╟─f0e54010-5a21-11ed-12bc-45182a10006b
# ╠═379d9b77-5d2f-4477-81ab-1957b75f13e1
# ╠═ac6bf8ee-5d9f-4216-a71e-973b63fbc320
# ╠═c58cf6e0-f332-44c1-945d-c162bb2465f7
# ╠═f9769a5c-6468-41ed-bd1b-0e1ae48b6a27
# ╠═64c7043a-9a38-41a6-b79c-a8604fa18b91
# ╠═2872a885-da17-4fd3-a3f4-e32a3871adc2
# ╠═e1818076-86d5-48df-a448-b517e6cace30
# ╠═e846aa72-5b03-4af3-b8f9-54da7411ff98
# ╠═908011e2-e506-496c-944e-3a346398f2f6
# ╠═7cb17280-2003-4825-90b7-50ef837b54e5
# ╠═4b3ce556-34b9-4308-8b25-b054995c1ba4
# ╟─9361923a-6019-44ee-811f-68d217b90aef
# ╠═c83d84e4-4f4d-4a97-b50f-1a5805e54c32
# ╠═65be9c7e-c4cc-4c03-8537-651f648c2f44
# ╠═d611d92b-872e-4cfe-8864-1e706f7c6fd6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
