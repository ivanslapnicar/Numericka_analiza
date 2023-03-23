### A Pluto.jl notebook ###
# v0.19.22

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
	c=0.0
	for n=1:steps
    	s=√((1-c)/2)
    	c=√(1-s^2)
    	println("2n = ", 2^(n+1), ", približna vrijednost od π = ",2^(n+1)*s)
	end
end

# ╔═╡ 4d68f850-027b-11eb-14a1-1741ea4ce7c2
# Točan π
π

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
* Može se pogledati LLVM kod i asemblerski kod.
"""

# ╔═╡ cb1cb5b0-0ee8-11eb-10c2-59d499c44695
# Sadržaj paketa
varinfo(Polynomials)

# ╔═╡ d63410b0-0ee8-11eb-31a3-e71de384fadb
# ?Polynomial

# ╔═╡ db5f4eb0-0ee8-11eb-3e4f-670b968e851d
p=Polynomial([c,b,a])

# ╔═╡ 21ed22d0-0ee9-11eb-2407-3b2f0eb6ab8b
roots(p)

# ╔═╡ 298912b0-0ee9-11eb-2c06-11cd5e89feed
@which roots(p)

# ╔═╡ 32961ec2-0ee9-11eb-2a47-1568df8b95d7
md"
Može se pogledati i datoteku na GitHub-u za najnoviju verziju.
"

# ╔═╡ 3e6d474e-0ee9-11eb-278f-4d620b298db3
# Ispis je u Julia terminalu
@code_llvm korijeni(1,0,9)

# ╔═╡ 4b09e2c0-0ee9-11eb-38c6-db7ffd4c00b3
# ispis je u Julia terminalu
@code_native korijeni(1,0,7)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Polynomials = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"

[compat]
PlutoUI = "~0.7.14"
Polynomials = "~2.0.15"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "31d0151f5716b655421d9d75b7fa74cc4e744df2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.39.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[HypertextLiteral]]
git-tree-sha1 = "72053798e1be56026b81d4e2682dbe58922e5ec9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.0"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[Intervals]]
deps = ["Dates", "Printf", "RecipesBase", "Serialization", "TimeZones"]
git-tree-sha1 = "323a38ed1952d30586d0fe03412cde9399d3618b"
uuid = "d8418881-c3e1-53bb-8760-2df7ec849ed5"
version = "1.5.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "29714d0a7a8083bba8427a4fbfb00a540c681ce7"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.3"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "3927848ccebcc165952dc0d9ac9aa274a87bfe01"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "0.2.20"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "a8709b968a1ea6abc2dc1967cb1db6ac9a00dfb6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.5"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[PlutoUI]]
deps = ["Base64", "Dates", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "d1fb76655a95bf6ea4348d7197b22e889a4375f4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.14"

[[Polynomials]]
deps = ["Intervals", "LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "029d2a5d0e6c2b5d87ac690aa58dcf40c2e2acb1"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "2.0.15"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[TimeZones]]
deps = ["Dates", "Future", "LazyArtifacts", "Mocking", "Pkg", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "6c9040665b2da00d30143261aea22c7427aada1c"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.5.7"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
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
# ╠═21ed22d0-0ee9-11eb-2407-3b2f0eb6ab8b
# ╠═298912b0-0ee9-11eb-2c06-11cd5e89feed
# ╟─32961ec2-0ee9-11eb-2a47-1568df8b95d7
# ╠═3e6d474e-0ee9-11eb-278f-4d620b298db3
# ╠═4b09e2c0-0ee9-11eb-38c6-db7ffd4c00b3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
