### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ ad5d9c87-f4f6-4f2a-8f69-1b4a72b98dec
using PlutoUI, Polynomials

# ╔═╡ 64255f65-cdd6-4288-8da1-b3a4ffb26536
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ a481f540-027d-11eb-10b5-5bd4e827e274
md"""
# Uvod

__Numerička analiza je znanost o računanju rješenja problema koji su matematički postavljeni u polju relanih ili kompleksnih brojeva.__  

Navedimo dva primjera.
"""

# ╔═╡ c2a8c4e0-027d-11eb-26fe-1f7bf12091f4
md"""
# Računanje broja $\pi \approx 3.14159265358979$

Broj $\pi$ je omjer opsega kružnice i njenog promjera. 
"""

# ╔═╡ eab6ace0-027d-11eb-0f82-db588f0a253b
md"""

## Arhimedova metoda

Upišite pravilni $n$-terokut u kružnicu radijusa $1$. Izračunajte opseg njegove gornje polovice.
To je najlakše napraviti kada je $n=2^k$. Koristeći geometriju,

$$
p(n) = 2n \sin \: \frac{\pi}{2n}.$$

Premda se ova formula poziva na $\pi$, vrijednosti $p(n)$ za $n=1,2,\cdots$, možemo izračunati bez poznavanja $\pi$ i bez korištenja funkcije $\sin$ na kalkulatoru (Računanje sinusa zahtijeva poznavanje broja $\pi$!).
"""

# ╔═╡ 3d725c70-027b-11eb-10cb-5b6e56311aa1
md"""
Znamo da je

$
\sin \: \frac{\pi}{2} =1, \quad \sin \: \frac{\pi}{4} = \frac{\sqrt{2}}{2}.$

Dakle, 

$
p(1)= 2 \cdot \sin \: \frac{\pi}{2} =2, \quad p(2) = 4 \cdot \sin \: \frac{\pi}{4} = 2\sqrt{2}= 2.828427125.$

Ali, što je $p(4) = 8 \cdot \sin \: \displaystyle\frac{\pi}{8}$. Formula za sinus polovice kuta daje

$$
\begin{aligned}
\sin \: \frac{\theta}{2} = \sqrt{\frac{1-\cos \: \theta }{2}},\\ \\
\cos \: \theta = \sqrt{1-\sin^2 \: \theta }.\end{aligned}$$


Vrijedi

$
\sin \: \frac{\pi}{8} = 0.382683432.$

Stoga je

$
p(4) = 8 \cdot \sin \: \frac{\pi}{8}=3.061467459.$

Nastavljajući postupak imamo

n | $\sin\displaystyle\frac{\pi}{2n}$ | $p(n)$
:---|:---|:---
1 | 1 | 2
2 | $\sqrt{0.5}$ | 2.82842712
4 |  0.382683432 | 3.061467459
8 | 0.195090322 | 3.121445152
16 | 0.09801714 | 3.136548491
32 | 0.049067674 | 3.140331157

Ova metoda je spora, ali "sigurna". 
"""

# ╔═╡ 3dea5041-85b9-4d8b-9dcf-2d39c8481314
steps=30

# ╔═╡ a47b53a0-027f-11eb-256d-f7d5ce0f97e3
let
	c=Float64(0.0)
	for n=1:steps
    	s=√((1-c)/2)
    	c=√(1-s^2)
    	println("2n = ", 2^(n+1), ", približna vrijednost od π = ",2^(n+1)*s)
	end
end

# ╔═╡ 4d68f850-027b-11eb-14a1-1741ea4ce7c2
# Točan π
π

# ╔═╡ fc4ecaf9-d735-4fc2-aa48-e7dfda2a5727
supertype(Float64)

# ╔═╡ f68459d5-b408-472b-b3f5-a68da07cf1ee
subtypes(Real)

# ╔═╡ 5d3240fd-7c59-4b6d-baaa-9a91c2a523a8
supertype(AbstractFloat)

# ╔═╡ a6664117-a273-4d8a-9f49-e5455d916f3c
subtypes(Number)

# ╔═╡ 6ac49527-05f0-4a0a-ad99-aba5d4b07c17
md"""
__Zadatak:__ Povečajte broj koraka i objasnite što se događa.
"""

# ╔═╡ 8443676d-70db-4f63-8b8a-a7e96175f813
md"""
Kasnije ćemo opisati modernije poboljšanje koje je i brže.
Uz $h=1/(2n)$ vrijedi

$$p(n) = \frac{\sin \: \pi h}{h} = \pi -a_2 h^2 + a_4 h^4 - \cdots,$$

gdje je $a_k = \displaystyle\frac{\pi^{k+1}}{(k+1)!}$. Dakle, ova formula teži k $\pi$ po stopi približno jednakoj
$O(h^2)$. To je _aproksimacijski problem_. Ne postoji konačan algoritam za računanje broja $\pi$, jer se radi o transcedentnom broju (iracionalan i nije korijen niti jednog polinoma s cjelobrojnim koeficijentima). 
Međutim, možemo ga aproksimirati po volji točno.
"""

# ╔═╡ bd7a0ad8-8b31-42cf-94e9-881711a3c40e
md"""
# Kvdratna jednadžba

Sljedeći problem ilustrira potpuno različite fenomene. 

Izračunajmo korijene $p(x) = ax^2 + bx +c =0$ za realne konstante $a$, $b$ i $c$.

Rješenja su

$$
x_{1,2} = \frac{ -b \pm \sqrt{b^2 - 4ac}}{2a}. \tag{1}$$

Kako ćemo pristupiti izradi programa za računanje ovih korijena?

Trebamo riješiti pet posebnih slučajeva.

__Slučaj I.__ $a=0, b \neq 0$.

Ovo više nije kvadratna jednadžba, nego linearna. Jedino rješenje je

$$
x_1 = -c/b.$$

__Slučaj II.__ $a=b=0$

Ako je $c \neq 0$, nema rješenja. Ako je $c=0$, bilo koji $x$ je rješenje.

Slučajevi koje smo rješavali u srednjoj školi koriste diskriminantu
$b^2 - 4ac$.

__Slučaj III.__ $b^2 -4ac < 0$. Dva kompleksna rješenje (nisu realna).

$$
x_{1,2} = -\frac{b}{2a} \pm \mathbf{i} \frac{\sqrt{4ac-b^2}}{2a}, \quad \mathbf{i}^2 = -1.$$

__Slučaj IV.__ $b^2 - 4ac =0$. Jedan dvostruki korijen (realni)

$$
x_1 = x_2 = -\frac{b}{2a}.$$

__Slučaj V.__ $b^2 -4ac > 0$. Dva različita realna rješenje.
Koristimo formulu (1)?
"""

# ╔═╡ b5fe78a7-5258-41ba-9544-c0cd64496ef5
begin
	a=1
	b=2
	c=10.0^(-17)
	x₁=(-b-√(b*b-4*a*c))/(2*a)
	x₂=(-b+√(b*b-4*a*c))/(2*a)
	
	x₁,x₂
end

# ╔═╡ d7a46520-824d-11eb-2090-ad945a312f56
x₁

# ╔═╡ dfabe810-824d-11eb-3659-cf6fd64c38af
x₂

# ╔═╡ bca1eb35-dd3d-41af-8130-a8035551ec49
md"""
Dva realna rješenja su (prvih 17 znamenki)

$$
x_1 = -2, \quad x_2 = -5 \cdot 10^{-18}.$$

Gornji algoritam izračuna $x_1$ točno, ali $x_2 = 0$. Standardni brojevi s plivajućim zarezom (floating-point) u dvostrukoj točnosti, `Float64`, spremaju približno $15$ decimalnih znamenki (54 binarne znamenke) i u tih 15 znamenki je $\sqrt{b^2 -4ac}-b =0$. To je zato što oduzimamo dva bliska broja, a jedan od njih je približno točan, pa je razlika isključivo "greške zaokruživanja". Jednostavno zapažanje rješava ovaj problem.

"Veliki" korijen kvadratne jednadžbe u __slučaju V__ je

$$
x_1 = \frac{ -b - \mathrm{sign}(b) \sqrt{b^2 - 4ac}}{2a},$$

a dva korijena zadovoljavaju 

$$
x_1 x_2 = \frac{c}{a}.$$

Primjetite da, osim unutar kvadratnog korijena, zbrajamo brojeve istog predznaka!
Nakon nekoliko transformacija imamo formulu za sitni korijen

$$
x_2 = \frac{c}{a x_1} = \frac{-2c}{ b + \mathrm{sign}(b) \sqrt{b^2 - 4ac}}.$$

Koristeći ovu formulu oba korijena računamo blizu točnosti stroja.

U ovom primjeru imali smo egzaktnu formulu, ali u aritmetici plivajućeg zareza standarna formula daje rezultate koji se značajno razlikuju od rezultata u realnoj aritmetici.
"""

# ╔═╡ 2f084ae2-1e38-4c0b-83b6-da9f025565b1
begin
	x₃=(-b-sign(b)*√(b*b-4*a*c))/(2*a)
	x₄=-2*c/(b+sign(b)*√(b*b-4*a*c))
	x₃,x₄
end

# ╔═╡ f8cf034f-b149-4ff9-8a7e-bc8143581d51
md"""
Sljedeća funkcija implementira svih pet sučajeva. Probajte funkciju na različitim ulazima tako da pokrijete sve slučajeve.
"""

# ╔═╡ 0738969f-7c4f-4f8e-92de-62cb0cd6510c
function korijeni(a,b,c)
    # Funkcija za računanje korijena kvadratne jednadžbe
    # sa zadanim koeficijentima a, b i c.
    # Ova funkcija ne vodi računa o skaliranju. 
    if a==0
        #  Provjerimo posebne slučajeve za a=0 
        if b==0
            if c==0
                return "svi brojevi su rješenja a=b=c=0"
            else
                return "nema rješenja a=b=0, c ne 0"
            end
        else
            x₁=-c/b
            x₂=x₁
            poruka="jedno rješenje a=0"
        end
    else
        Δ= b*b-4*a*c
        if Δ < 0
            # Dva kompleksna rješenja izračunata pomoću realne aritmetike
            ximaginarni=√(-Δ)/(2*a)
            xrealni=-b/(2*a)
            x₁=xrealni+im*ximaginarni
            # x₂ je kompleksno konjugirani x₁, 
            # x₂ = xrealni - im*ximaginarni
            x₂=conj(x₁)
        else
            if b==0
                # Julia lako računa s kompleksnim brojevima,
                # pa u ovom slučaju možemo koristiti formulu.
                x₁=√(-c)/a
                x₂=-x₁
            else
                # Slučaj s dva različita realna korijena.
                x₁=(-b-sign(b)*√(Δ))/(2*a)
                x₂=-2*c/(b+sign(b)*√(Δ))
            end
        end
        poruka="korijeni su dobri"
    end
    x₁, x₂, poruka
end

# ╔═╡ cd29c4bc-9cfa-40bc-bc80-74ee5a635a94
korijeni(1,0,9)

# ╔═╡ af57717b-4dca-4b54-8272-e98fe76fd845
korijeni(a,b,c)

# ╔═╡ 9315e7e0-0ee8-11eb-2eeb-4160ebee122d
md"""
# Julia je brza i otvorena

`Julia` je __brza__ i prva je rješila problem dva jezika, npr. Matlab ili Python za razvoj, C ili C++ ili FORTRAN za brzinu. Svaka funkcija se nakon prvog poziva kompajlira pomoću [JIT kompajlera](https://en.wikipedia.org/wiki/Just-in-time_compilation), u ovom slučaju [LLVM](https://en.wikipedia.org/wiki/LLVM). 

`Julia` je __otvorena__ jer
* Kompletan izvorni kod je uvijek dostupan na GitHub-u.
* MIT licenca
* Makro `@which` olakšava snalaženje.
* __Može se pogledati LLVM kod i asemblerski kod.__
"""

# ╔═╡ cb1cb5b0-0ee8-11eb-10c2-59d499c44695
# Sadržaj paketa
varinfo(Polynomials)

# ╔═╡ d63410b0-0ee8-11eb-31a3-e71de384fadb
# ?Polynomial

# ╔═╡ db5f4eb0-0ee8-11eb-3e4f-670b968e851d
p=Polynomial([c,b,a])

# ╔═╡ 298912b0-0ee9-11eb-2c06-11cd5e89feed
@which roots(p)

# ╔═╡ 97780d3c-80c4-4c7a-bb20-70c63bb12b5f
roots(p)

# ╔═╡ 32961ec2-0ee9-11eb-2a47-1568df8b95d7
md"
Može se pogledati i datoteku na GitHub-u za najnoviju verziju.
"

# ╔═╡ 3e6d474e-0ee9-11eb-278f-4d620b298db3
# Ispis je u Julia terminalu
# @code_llvm korijeni(1,0,9)

# ╔═╡ 4b09e2c0-0ee9-11eb-38c6-db7ffd4c00b3
# ispis je u Julia terminalu
# @code_native korijeni(1,0,7)

# ╔═╡ 004fe019-b74a-467a-9e41-919e8a95a13c
p₁=Polynomial([2//1,3//1,4//1])

# ╔═╡ 7aab7932-5ec7-4763-af8f-d9cbd4f94fce
roots(p₁)

# ╔═╡ 54a31c16-424c-49c8-92d6-dbba0d4443fa
korijeni(4//1,3//1,2//1)

# ╔═╡ aaedea4f-5ee0-4a23-b0e8-77d857fdcb97
Rational{BigInt}(√2)*Rational(√3)

# ╔═╡ e651042d-b285-4bbf-98b8-f1b6b0e66a77
Rational{BigInt}(√6)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Polynomials = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"

[compat]
PlutoUI = "~0.7.54"
Polynomials = "~4.0.6"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.1"
manifest_format = "2.0"
project_hash = "c49261613a668006e073727a005e4c9ff5fa2e0a"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "c278dfab760520b8bb7e9511b968bf4ba38b7acc"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.3"

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

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

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

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

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

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "b211c553c199c111d998ecdaf7623d1b89b69f93"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.12"

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
git-tree-sha1 = "bd7c69c7f7173097e7b5e1be07cee2b8b7447f51"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.54"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "a9c7a523d5ed375be3983db190f6a5874ae9286d"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.6"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

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

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

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
# ╠═ad5d9c87-f4f6-4f2a-8f69-1b4a72b98dec
# ╠═64255f65-cdd6-4288-8da1-b3a4ffb26536
# ╟─a481f540-027d-11eb-10b5-5bd4e827e274
# ╟─c2a8c4e0-027d-11eb-26fe-1f7bf12091f4
# ╟─eab6ace0-027d-11eb-0f82-db588f0a253b
# ╟─3d725c70-027b-11eb-10cb-5b6e56311aa1
# ╠═3dea5041-85b9-4d8b-9dcf-2d39c8481314
# ╠═a47b53a0-027f-11eb-256d-f7d5ce0f97e3
# ╠═4d68f850-027b-11eb-14a1-1741ea4ce7c2
# ╠═fc4ecaf9-d735-4fc2-aa48-e7dfda2a5727
# ╠═f68459d5-b408-472b-b3f5-a68da07cf1ee
# ╠═5d3240fd-7c59-4b6d-baaa-9a91c2a523a8
# ╠═a6664117-a273-4d8a-9f49-e5455d916f3c
# ╟─6ac49527-05f0-4a0a-ad99-aba5d4b07c17
# ╟─8443676d-70db-4f63-8b8a-a7e96175f813
# ╟─bd7a0ad8-8b31-42cf-94e9-881711a3c40e
# ╠═b5fe78a7-5258-41ba-9544-c0cd64496ef5
# ╠═d7a46520-824d-11eb-2090-ad945a312f56
# ╠═dfabe810-824d-11eb-3659-cf6fd64c38af
# ╟─bca1eb35-dd3d-41af-8130-a8035551ec49
# ╠═2f084ae2-1e38-4c0b-83b6-da9f025565b1
# ╟─f8cf034f-b149-4ff9-8a7e-bc8143581d51
# ╠═0738969f-7c4f-4f8e-92de-62cb0cd6510c
# ╠═cd29c4bc-9cfa-40bc-bc80-74ee5a635a94
# ╠═af57717b-4dca-4b54-8272-e98fe76fd845
# ╟─9315e7e0-0ee8-11eb-2eeb-4160ebee122d
# ╠═cb1cb5b0-0ee8-11eb-10c2-59d499c44695
# ╠═d63410b0-0ee8-11eb-31a3-e71de384fadb
# ╠═db5f4eb0-0ee8-11eb-3e4f-670b968e851d
# ╠═298912b0-0ee9-11eb-2c06-11cd5e89feed
# ╠═97780d3c-80c4-4c7a-bb20-70c63bb12b5f
# ╟─32961ec2-0ee9-11eb-2a47-1568df8b95d7
# ╠═3e6d474e-0ee9-11eb-278f-4d620b298db3
# ╠═4b09e2c0-0ee9-11eb-38c6-db7ffd4c00b3
# ╠═004fe019-b74a-467a-9e41-919e8a95a13c
# ╠═7aab7932-5ec7-4763-af8f-d9cbd4f94fce
# ╠═54a31c16-424c-49c8-92d6-dbba0d4443fa
# ╠═aaedea4f-5ee0-4a23-b0e8-77d857fdcb97
# ╠═e651042d-b285-4bbf-98b8-f1b6b0e66a77
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
