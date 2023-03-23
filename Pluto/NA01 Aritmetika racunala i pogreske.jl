### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 62d67e4a-ae9f-43a2-bcc8-53ee5da9e268
using PlutoUI

# ╔═╡ 0425b895-a11b-49ee-b715-228384218624
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ 76d37869-e20b-4211-8227-1f0616e3d8f2
md"""
# Aritmetika računala i pogreške

# Apsolutna i relativna pogreška

Neka $\alpha$ aproksimira $a$. Onda je

$$
\begin{aligned}
err&=|a-\alpha| \\ \\ relerr&=\frac{err}{|a|}=\frac{|a-\alpha|}{|a|}.\end{aligned}$$
"""

# ╔═╡ 5c635357-8163-4954-949a-999dc48998f0
md"
α = $(@bind α Slider(0:0.1:6,show_value=true))
"

# ╔═╡ 9ccbb154-8618-4190-bfce-985ba66c8380
begin
	# Probajte α=a:0.01:2a
	a=5.0
	err=abs(a-α)
	relerr=err/abs(a)
	α, err, relerr
end

# ╔═╡ 91448849-f9f0-459b-87c0-b9fc5a386770
md"""
# Aritmetika plivajućeg zareza

Korisna knjiga o IEEE standardu za aritmetiku plivajućeg zareza:

M. Overton, Numerical Computing with IEEE Floating Point Arithmetic, SIAM Publications, Philadephia, 2001.

Koristan članak:

[David Goldberg, What Every Computer Scientist Should Know About Floating-Point Arithmetic](https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html).

## Sustav brojeva s plivajućim zarezom

 $x$ je __broj s plivajućim zarezom__ (_floating point number_) ako je oblika

$$x = \pm d \cdot \beta^e, \quad \beta \in \{ 2,10 \}$$

__Baza__ 2 je za standardna računala opće namjene, __baza__ 10 je za džepne kalkulatore.

 $e$ je __eksponent__ i zadovoljava

$$e_{\min} \leq e \leq e_{\max},\quad e_{\min} < 0 < e_{\max}$$

Pretpostavit ćemo da koristimo aritmetiku s bazom 2, ali ćemo primjere uglavnom davati u bazi 10.

__Mantisa__ $d$ je oblika

$$
\begin{aligned}
	d &= 0.d_1 \dots d_t = d_1 \beta^{-1} + d_2 \beta^{-2}
	+ \dots + d_t \beta^{-t}\\
d  &\in \{ 0,1\}\\
	d_1 &= 1 \qquad \mbox{ normalized }   \\
	d_1 &= 0 \qquad \mbox{ unnormalized }   \\
\end{aligned}$$

Standardni oblik brojeva s plivajućim zarezom je normalizirani oblik osim pri dnu raspona eksponenata.

Tijekom ulaza  / izlaza (_input/output_) brojevi se konvertiraju iz binarnog oblika u dekadski i natrag.

Aritmetika računala je standardizirana pomoću IEEE 754 standard for binary arithmetic.  Gotovo sva računala (procesori) rade po tom standardu (➡ `paranoia`).
"""

# ╔═╡ 1f8858f4-104d-4dd8-99de-b4bbf278e720
md"""
## Strojna jedinica i točnost stroja

Skup

$$
\{x \colon \lfloor \log_2 \: |x| \rfloor \in [e_{min},e_{max}] \}$$

je podskup skupa realnih brojeva koji se nalaze u normaliziranom rasponu brojeva s plivajućim zarezom. $fl(x)$ je realan broj $x$ zaokružen na najbliži broj s plivajućim zarezom.

__Strojna jedinica__ (_machine unit_) je najveća relativna udaljenost između realnog broja koji se nalazi u rasponu brojeva s plivajućim zarezom i najbližeg broja s plivajućim zarezom,

$$
\epsilon_M = \max_{\lfloor \log_2
\:|x|\rfloor \in
[e_{\min},e_{\max}]} \frac{|x - fl(x)|}{|x|}  = 2^{-t}$$

__Točnost stroja__ (_machine precision_) je relativna udaljenost između dva susjedna broja s plivajućim zarezom. Za $\beta=2$ očito vrijedi $\epsilon=2\epsilon_M$.

Važni primjeri su __jednostruka točnost__ (_single precision_) i 
__dvostruka točnost__ (_double precision_):

__IEEE Standard Single Precision (Float32)__  $\beta = 2$, $t = 24$

$$
\begin{aligned}
\epsilon_M  &= 2^{-24} \approx	5.9605 \times 10^{-8}\\
\epsilon &=2^{-23} \approx 1.1920 \times 10^{-7} \\
e_{\min} &= - 126,\quad e_{\max} = 127.
\end{aligned}$$


__IEEE Standard Double Precision (Float 64)__ $\beta =2$,$t = 53$

$$
\begin{aligned}
\epsilon_M &= 2^{-53} \approx 1.1102 \times 10^{-16}\\
\epsilon &=2^{-52} \approx 2.2204 \times 10^{-16}\\
e_{\min} &= -1022,\quad e_{\max} = 1023.
\end{aligned}$$

Izračunajmo $\epsilon$ kao najmanji pozitivni broj s plivajućim zarezom takav da je $1+\epsilon\neq 1$.
"""

# ╔═╡ 030a0e5f-f9e1-4696-af85-b890eaa129d7
begin
	b₀=1.0
	a₀=2.0
	while (b₀+a₀)!=b₀
	    a₀/=2
	    println(a₀) # Output is in terminal
	end
	a₀
end

# ╔═╡ 6065926a-21a0-4b98-b7f9-ecb96693f12d
1+a₀==1

# ╔═╡ 82a6d40f-4c81-4ff1-8e02-3f427c413f61
md"""
MATLAB naredba `eps` i Julia funkcija `eps()` daju $\epsilon = 2.2204 \times 10^{-16}$.
"""

# ╔═╡ 7c28c479-912d-4a12-bc38-d15c6a3f0501
eps()

# ╔═╡ f5a5a27d-27bc-49b5-b245-c32a1f3fe13c
# Što je ovo?
eps(64.0)

# ╔═╡ aa042b72-246d-432a-9a1c-d014da8ef957
md"""
Julia ima sustav vrsta podataka (_type system_) u kojem je tip `Float64` pod-tip tipa `AbstractFloat`, koji ima četiri pod-tipa.
To su tipovi `Float64` i `Float32`, zatim tip `Float16` koji koristi samo dva bajta računalne memorije i tip `BigFloat` koji ima mantisu od 256-bitova.
"""

# ╔═╡ dac635a5-7f41-478c-bede-5894410a0b6a
supertype(Float64)

# ╔═╡ 8c7b5693-2913-4cd3-9491-d8e57f101de8
subtypes(AbstractFloat)

# ╔═╡ 887cc5e4-ae65-4255-8a54-5425115e4618
for T in (Float16, Float32, Float64, BigFloat)
    println(eps(T))
end

# ╔═╡ b105ef44-00dc-4fe5-9732-499d14e51a51
2^(-10), 2^(-23), 2^(-52), 2^(-255)

# ╔═╡ 51356548-f58f-40a7-98b0-1b92ebeba3ed
md"""
## Osnovne operacije s plivajućim zarezom

Četiri osnovne operacije su zbrajanje ($+$), oduzimanje ($-$), množenje ($*$) i dijeljenje ($/$). Neka $\odot$ označava operaciju,

$$
\odot \in \{ + , - , *,/\}.$$

U aritmetici s plivajućim zarezom sa strojnom jednicom $\epsilon_M$, razumno je očekivati da za bilo koja dva broja s plivajućim zarezom $x$ and $y$ vrijedi (ukoliko je rezultat u rasponu brojeva s plivajućim zarezom)

$$
fl(x\;\odot\;y) = (x \; \odot\; y)\;(1 + \xi),\quad
|\xi| \leq \epsilon_M.$$

Kod dijeljenja pretpostavljamo $y \neq 0$. Svako računalo s ugrađenim IEEE standardom mora poštovati ovo pravilo. Zaokruživanje je jedno ograničenje aritmetike s plivajućim zarezom koje realna aritmetika nema. Iz ovog pravila možete lako zaključiti da će, ukoliko zbrajamo brojeve istog predznaka, množimo i dijelimo, rezultat u aritmetici plivajućeg zareza gotovo uvijek biti vrlo blizu točnom rezultatu. Poteškoće nastaju kada su $x$ i/ili $y$ već zaokruženi i imaju različite predznake te ih zbrajamo, ili imaju isti predznak te ih oduzimamo.

Drugim riječima, neka je 

$$
\tilde{x}= x(1+\delta_x), \quad \tilde{y} = y(1+\delta_y),$$

gdje su $x$ i $y$ točni rezultati nekog proračuna a $\tilde{x}$ i $\tilde{y}$ su zaokruženi rezultati u aritmetici plivajućeg zareza i vrijedi $|\delta_x|, |\delta_y| \leq \delta$ za neki mali $\delta$. Pretpostavimo i da  $x$ i $y$ imaju isti predznak. Neka je

$$
z=x-y,\quad  \tilde{z} = fl(\tilde{x} -\tilde{y}).$$

Onda je

$$
\begin{aligned}
\tilde{z} &=(\tilde{x}-\tilde{y})(1+\xi)= x(1+\delta_x)(1+\xi) -y(1+\delta_y)(1+\xi)
=x-y + \delta_z,
\end{aligned}$$

pri čemu je $|\xi| \leq \epsilon$ i

$$
\delta_z = (x-y)\xi + (x\delta_x -y\delta_y)(1+\xi).$$

Najbolja moguća ograda za $|\delta_z|$ je

$$
\begin{aligned}
|\delta_z| &\leq |x-y||\xi| + (|x||\delta_x| + |y||\delta_y|)(1+|\xi|) \\
& \leq |x-y| \epsilon_M + (|x|+|y|)\,\delta\,(1+\epsilon_M).
\end{aligned}$$

Prema tome, relativna pogreška od $z$ je

$$
\begin{aligned}
\frac{|\tilde{z}-z|}{|z|}&=\frac{|\delta_z|}{|z|}
\leq \epsilon_M + (1+\epsilon_M)\,\delta\,\frac{|x|+|y|}{|x-y|}\approx \delta \,\frac{|x|+|y|}{|x-y|}.
\end{aligned}$$

Ako je $|x-y| << |x|+|y|$, onda utjecaj zaokruživanja prilikom oduzivanja nije važan, ali pogreška iz prethodnih proračuna kod $x$ and $y$ može imati ogroman utjecaj. Taj efekat se zove __propagacija__ (_propagation_) i može dramatično promijeniti rezultat izračuna, što ćemo vidjeti kasnije na nekim primjerima. 

Zaokruživanje je prvo važno ograničenje aritmetike s plivajućim zarezom. Drugo ograničenje je raspon brojeva.
"""

# ╔═╡ 3170458d-931a-41a1-8715-b41de07aa3c6
md"""
## Rasponi brojeva

Aritmetika s plivajućim zarezom ima najveći i najmanji broj. Opišimo prvo najveći.

__Najveći broj u računalu__ $\Omega$

U bazi $2$ s mantisom od $t$ bitova najveći broj u računalu je

$$
\Omega = (1 - 2^{-t}) \cdot 2^{e_{\max+1}}$$

Brojevi apsulutno veći od $\Omega$ se spremaju kao `Inf` ($\infty$) ili `-Inf` ($-\infty$). Kažemo da se dogodio __pretek__ (_owerflow_) (➡ Ariane 5).


_IEEE Standard Single Precision_ (`Float32`)

$$
\quad \Omega = 3.4028\times 10^{38}$$

_IEEE Standard Double Precision_ (`Float64`)

$$
\Omega = 1.79777 \times 10^{308}$$

MATLAB naredba `realmax` i Julia funkcija `floatmax()` daju $\Omega$.
"""

# ╔═╡ 2ee4617f-70cf-4ecb-99b2-7167cc6b34d6
md"""
__Najmanji broj u računalu__ $\omega$

Definicija najmanjeg (pozitivnog) broja u računalu je nešto složenija.

Najmanj broj u računalu je

$$
\omega = 2^{1-t} 2^{e_{\min}}.$$

Ukoliko izračun daje broj koji je apsolutno manji od $\omega$, dogodi se ono što se zove __podtek__ (_underflow_), rezultat se postavi na $0$ ili $-0$. SAko programer to odabere, podtek može rezultirati greškom, ali u većini slučajava podteci su bezopasni. 


_IEEE Standard Single Precision_ (`Float32`):

$$
\omega = 2^{-23- 126} = 2^{-149} \approx  1.4013 \times 10^{-45}.$$

U MATLAB-u možemo koristiti naredbu `omega= eps('single')*realmin('single')`.


_IEEE Standard Double Precision_ (`Float64`):

$$
\omega= 2^{-1022-52} = 2^{-1074} \approx  4.9407 \times 10^{-324}$$

Odgovarajuća MATLAB naredba je `omega = eps*realmin`, a odgovarajuća Julia naredba je  `floatmin()*eps()`.


__Važan detalj__

Brojevi pri dnu raspona brojeva nisu normalizirani. 

MATLAB naredba `realmin` daje

$$
\omega_{koristan} \approx 2.2251 \times 10^{-308}.$$

Ovaj broj se još zove namanji _koristan_ broj s plivajućim zarezom jer je

$$
1/\omega_{koristan} \leq \Omega,$$

pri čemu je $\omega_{koristan}$ je normaliziran.

Najmanji broj s plivajućim zarezom $\omega$ je oblika

$$
0.0 \cdots 01 \times 2^{e_{\min}} \quad \cdots\quad
\mbox{postupni podtek - Gradual Underflow}$$

Prije uvođenja IEEE standarda, najmanji broj u većini računala je bio

$$
0.10 \cdots 0 \times 2^{e_{\min}} \qquad \cdots
\mbox{ normaliziran}$$

Starija računala (prije 1985) su brojeve manje od najmanjeg korisnog broja s plivajućim zarezom postavljala na nulu. Ova promjena je bila jedna on najkontroverznijih osobina IEEE standarda.

__Primjer.__ $\beta = 10$, $-5 \leq e \leq 5$

$$
\begin{aligned}
x & = 0.1957 \times 10^{-5}   \\
y & = 0.1942 \times 10^{-5}
\end{aligned}$$

Izračunajmo $fl(x - y)$. Što će se dogoditi?

$$
0.1957 \times 10^{-5}-0.1942 \times 10^{-5}  =0.0015 \times 10^{-5}$$

Prije 1985 računala bi stavila $fl(x - y)=0$.

Postupni podtek računa $fl(x - y)=0.0015 \times 10^{-5}$, odnosno garantira 
da za svaka dva broja s plivajućim zarezom $x$ i $y$ vrijedi

$$
fl(x - y) = 0 \mbox{ ako i samo ako je } x = y.$$
"""

# ╔═╡ eaf34525-9f16-468e-ab8b-54a5e0a8b40d
floatmax(27.0)

# ╔═╡ b9dd627b-5b1a-4113-bfdf-af4dc1901827
for T in (Float16, Float32, Float64, BigFloat)
    println((floatmin(T),floatmax(T)))
end

# ╔═╡ 0f7ba6b2-aede-4f09-b7c5-adca1295195b
1/floatmin(Float32),floatmax(Float32)

# ╔═╡ b8647d13-c400-4324-85f7-8994a1f2322f
for T in (Float16, Float32, Float64)
    println((floatmin(T)*eps(T)))
end

# ╔═╡ dec72eaf-d4c8-4f67-9e4e-8d71c6651796
1/floatmin()

# ╔═╡ e1b1c9aa-cff4-4155-80c2-9181be546116
1/(floatmin()*eps())

# ╔═╡ 20fa79ac-93c3-4e3f-bf31-f9c0f0e79d4a
md"""
## Binarni prikaz
"""

# ╔═╡ 1f8f4a15-0eca-4759-bee6-843306e35754
bitstring(Int32(0))

# ╔═╡ 8ca8e4e7-6450-4aaf-b208-6c0aa95ccb12
bitstring(-2)

# ╔═╡ f6f7081b-212a-4c7a-b91c-0d747f4b9fb2
bitstring(-0.0)

# ╔═╡ d7160016-c046-4ab8-b349-27f7a209d853
bitstring(0.0)

# ╔═╡ 10dcdca1-71fe-42eb-8f04-6122b3a03666
bitstring(1.0)

# ╔═╡ f85ec710-1f59-4826-a376-6672bf0552ae
bitstring(Float32(1.0))

# ╔═╡ 9b7e4278-2f1a-4fff-a6f7-2249f45aaccf
bitstring(Float32(2.0))

# ╔═╡ 40d953af-e5e8-42c0-be55-26d95e334195
2^7

# ╔═╡ 037cddff-3663-4fcb-8575-ea7a78e490b8
md"""
__Zadatak.__ Objasnite ove binarne prikaze.
"""

# ╔═╡ d0b306c4-6cbc-48dc-90ba-8ef4498f8d73
md"""
##  Posebne veličine $0$, $-0$, `Inf`,`-Inf` i `NaN`

Nula ima predznak:
"""

# ╔═╡ abe311d5-fbd1-40e2-8bce-2c341301deef
begin
	a₁=1.0
	b₁=0.0
	c₁=-b₁
	c₁,b₁==c₁
end

# ╔═╡ 88f1bf59-5eb4-4e2d-b28b-0d9b1004d5bd
a₁/b₁

# ╔═╡ be2cb230-df60-4db6-86e3-8f7d91e28a86
b₁+c₁

# ╔═╡ 36c79f63-7a8c-46c9-afa5-95bbed8fd598
begin
	d₁=a₁/b₁
	e₁=a₁/c₁
	d₁==e₁, 1/d₁==1/e₁
end

# ╔═╡ 92490cef-70af-417e-8bc8-91b1be0635cc
b₁/c₁

# ╔═╡ afc32ff2-014a-43ee-8e6a-d35a57434622
md"""
`NaN` (Not a Number) nastaje kao (vidi neodređene oblike iz Matematike 1):
"""

# ╔═╡ 03f1778e-c222-4221-a4e7-eabf1071d298
Inf+(-Inf),0*Inf, Inf/Inf, 0.0/0.0

# ╔═╡ f5f3c110-b31f-4c41-8eb5-121965a5e54d
md"""
U IEEE standardu brojevi s plivajućim zarezom i posebne veličine imaju sljedeće binarne zapise:

| Eksponent | Mantisa | Prikazuje |
| :-----    | :-----  | :-----    |
| $e=e_{\min}-1$ |  $d=0$     | $\pm 0$     |
| $e=e_{\min}-1$ |  $d\neq 0$ | $0.d\times 2^{e_\min}$ - denormalizirani brojevi |
| $e_{\min}\leq e\leq e\_\max$ |  $d$      | $1.d \times 2^e$ - standardni brojevi | 
| $e=e_{\max}+1$ |  $d=0$      |  $\pm$`Inf`     |
| $e=e_{\max}+1$ |  $d\neq 0$  |  `NaN`     |

"""

# ╔═╡ 94380193-9789-43db-8220-5ef210dc36aa
bitstring(Inf)

# ╔═╡ 949b08e6-f84b-475d-ab80-66ecc7532b8e
bitstring(-Inf)

# ╔═╡ 40702b4e-cecc-4a2c-ad88-d2613f49f71b
bitstring(0.0*Inf)

# ╔═╡ da406d8f-a48b-4a56-8f34-73f0ab837589
bitstring(0.0\0.0)

# ╔═╡ 307a27e6-8128-42be-9ba8-cbdacb498ada
md"""
IEEE aritmetika je zatvoren sustav:

 $\big\{$ floating point numbers,`Inf`,`-Inf`, `NaN`$\big\}$
$\stackrel{\odot}{\rightarrow}$
$\big\{$ floating point numbers,`Inf`,`-Inf`, `NaN` $\big\}$

bez obzira o kojoj operaciju $\odot$ se radi.

Programeri mogu korisno uopotrijebiti prethodna svojstva. Međutim, ako u zadaći iz programiranja dobijete
`NaN` ili `Inf` ili `-Inf`, vjerojatno se radi o grešci.
"""

# ╔═╡ 4e97e67b-8ef6-412d-8851-2b6635ff46e6
md"""
# Primjeri

## Korištenje razlike kvadrata

Izračunajte

$$
f(x) = \sqrt{1 + x^2} - 1, \quad \mbox{$x$ je blizu nule}.$$

Ova formula u standardnoj dvostrukoj točnosti daje $f(10^{-12}) = 0$.
"""

# ╔═╡ 88584373-9519-46fe-ae88-7199e76fb8f7
begin
	f(x)=√(1+x^2)-1
	[(x,f(x)) for x ∈ [1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12]]
end

# ╔═╡ 1215b318-8121-4adb-aa5d-139654607717
md"""
Trik s razlikom kvadrata daje

$$
\begin{aligned}
f(x) & \equiv (\sqrt{1 + x^2} - 1) \left( \frac{\sqrt{1 + x^2} + 1}{\sqrt{1 + x^2} + 1}\right) \\
& = \frac{x^2}{\sqrt{1+x^2} + 1}\equiv f_1(x),
\end{aligned}$$

odnosno,  $f_1(10^{-12}) = 0.5 \cdot 10^{-24}$. Ovaj rezultat je onoliko točan koliko možemo očekivati u standardnoj dvostrukoj točnosti.
"""

# ╔═╡ b9180215-218a-48d8-a662-2a5664c2e480
begin
	f₁(x)=x^2/(1+√(1+x^2))
    [(x,f₁(x)) for x ∈ [1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12]]
end

# ╔═╡ 209cad48-c496-41b4-8093-6e2a55cc469f
x=1e-12

# ╔═╡ 905dd1d7-5f1b-43b3-a461-651cc992bd1c
# Koristeći BigFloat
BigFloat(x)

# ╔═╡ c5d305a9-5a27-44d6-a610-50ab6d396f93
f(BigFloat(x))

# ╔═╡ 82ac21a4-384e-47e8-a964-95c8eb0889bd
md"""
## Kvadratne jednadžba

U egzaktnoj aritmetici kvadratna jednadžba

$$ax^2 + bx+c=0$$

ima korjene

$$
\begin{aligned}
x_1&=\frac{-b-\mathop{\mathrm{sign}}(b)\sqrt{b^2-4ac}}{2a} \\
x_2&\equiv\frac{-b+\mathop{\mathrm{sign}}(b)\sqrt{b^2-4ac}}{2a}= \frac{-b+\mathop{\mathrm{sign}}(b)\sqrt{b^2-4ac}}{2a}\cdot \frac{-b-\mathop{\mathrm{sign}}(b)\sqrt{b^2-4ac}}{-b-\mathop{\mathrm{sign}}(b)\sqrt{b^2-4ac}}
\\ &= \frac{2c}{-b-\mathop{\mathrm{sign}}(b)\sqrt{b^2-4ac}}\equiv x_3.
\end{aligned}$$
"""

# ╔═╡ 139debf3-10ad-48da-9c2d-a13033486fbf
begin
	a₂=2.0
	b₂=123456789.0
	c₂=4.0

	x₁(a,b,c)=(-b-sqrt(b*b-4*a*c))/(2.0*a)
	x₂(a,b,c)=(-b+sqrt(b*b-4*a*c))/(2.0*a)
	x₂ₐ(a,b,c)=(2*c)/(-b-sqrt(b*b-4*a*c))
	x₁(a₂,b₂,c₂),x₂(a₂,b₂,c₂),x₂ₐ(a₂,b₂,c₂)
end

# ╔═╡ d429eede-d77d-4709-b635-c3993fce4a47
md"""
Provjerimo koristeći `BigFloat`:
"""

# ╔═╡ ff5e3103-11cb-4438-aa19-c25af84fc3da
x₂ₐ(a₂,b₂,c₂)

# ╔═╡ ec755f39-5dba-4498-9760-07629ad00dc7
x₂(BigFloat(a₂),BigFloat(b₂),BigFloat(c₂))

# ╔═╡ c94ccf15-a0bd-4605-96f2-a27d0508fbeb
md"""
## Tangens i sinus
"""

# ╔═╡ 65f44437-4d44-4d93-a9bc-814ab9b6ee02
begin
	x₃=1e-10
	tan(x₃)-sin(x₃)
end

# ╔═╡ bf8fd692-0e01-4368-a52d-91edb4fe77c0
md"""
Međutim, trigonometrijski identiteti daju

$$
\begin{aligned}
\tan x - \sin x & = \tan x (1 - \cos x )
= \tan x (1-\cos x)\frac{1+\cos x}{1+\cos x}\\ & = \tan x \frac{1-\cos^2 x}{1+\cos x} \\
&= \tan x \sin^2 x \frac{1}{1+\cos x},
\end{aligned}$$

a Taylorova formula daje

$$
\begin{aligned}
\tan x &= x + \frac{x^3}{3} + \frac{2x^5}{15} + O(x^7) \\
\sin x &= x -\frac{x^3}{6} + \frac{x^5}{120}+O(x^7) \\
\tan x - \sin x &= \frac{x^3}{2} + \frac{7x^5}{120} +O(x^7).
\end{aligned}$$

Obe formule daju točan rezultat:
"""

# ╔═╡ d58c4727-448d-47b2-a194-2002304a4708
tan(x₃)*sin(x₃)^2/(1+cos(x₃)), x₃^3/2+7*x₃^5/120

# ╔═╡ 4035ef81-4348-4721-9448-48c569a65c9d
md"""
## Apsolutna vrijednost kompleksnog broja

Radi izbjegavanja podteka i preteka, umjesto standardne formule

$$
|z|=|x+iy|=\sqrt{x^2+y^2}$$

trebamo koristiti sljedeće formule (Objasnite!):

$$
M = \max \{ |x|,|y|\}, \quad m = \min \{ |x|,|y| \}, \quad r = \frac{m}{M}, \quad
|z| = M \sqrt{1+r^2}.$$

Ove formule su ugrađene u funkciju `abs()`.
"""

# ╔═╡ ac2caf5e-9874-44b9-ba5e-f0011940055f
z=2e+170+3e+175*im

# ╔═╡ b990e038-aff2-4d44-bb27-eb9de8c3977a
√(real(z)^2+imag(z)^2), abs(z)

# ╔═╡ c1acfea2-a60d-4433-a00b-6d3515274a18
begin
	z₁=2e-170+3e-175*im
	√(real(z₁)^2+imag(z₁)^2), abs(z₁)
end

# ╔═╡ b573a376-d60f-4f1d-b876-592dbcd47be4
md"""
__Zadatak.__ Usporedite funkciju  [hypot](https://en.wikipedia.org/wiki/Hypot) BLAS 1 funkciju `dnrm2.f`.
"""

# ╔═╡ 03679ace-1368-4be3-988c-1f0dcaa1407e
real(z)^2, imag(z)^2

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.9"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "438d35d2d95ae2c5e8780b330592b6de8494e779"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.3"

[[PlutoUI]]
deps = ["Base64", "Dates", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "Suppressor"]
git-tree-sha1 = "44e225d5837e2a2345e69a1d1e01ac2443ff9fcb"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.9"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
"""

# ╔═╡ Cell order:
# ╠═62d67e4a-ae9f-43a2-bcc8-53ee5da9e268
# ╠═0425b895-a11b-49ee-b715-228384218624
# ╟─76d37869-e20b-4211-8227-1f0616e3d8f2
# ╠═5c635357-8163-4954-949a-999dc48998f0
# ╠═9ccbb154-8618-4190-bfce-985ba66c8380
# ╟─91448849-f9f0-459b-87c0-b9fc5a386770
# ╟─1f8858f4-104d-4dd8-99de-b4bbf278e720
# ╠═030a0e5f-f9e1-4696-af85-b890eaa129d7
# ╠═6065926a-21a0-4b98-b7f9-ecb96693f12d
# ╟─82a6d40f-4c81-4ff1-8e02-3f427c413f61
# ╠═7c28c479-912d-4a12-bc38-d15c6a3f0501
# ╠═f5a5a27d-27bc-49b5-b245-c32a1f3fe13c
# ╟─aa042b72-246d-432a-9a1c-d014da8ef957
# ╠═dac635a5-7f41-478c-bede-5894410a0b6a
# ╠═8c7b5693-2913-4cd3-9491-d8e57f101de8
# ╠═887cc5e4-ae65-4255-8a54-5425115e4618
# ╠═b105ef44-00dc-4fe5-9732-499d14e51a51
# ╟─51356548-f58f-40a7-98b0-1b92ebeba3ed
# ╟─3170458d-931a-41a1-8715-b41de07aa3c6
# ╟─2ee4617f-70cf-4ecb-99b2-7167cc6b34d6
# ╠═eaf34525-9f16-468e-ab8b-54a5e0a8b40d
# ╠═b9dd627b-5b1a-4113-bfdf-af4dc1901827
# ╠═0f7ba6b2-aede-4f09-b7c5-adca1295195b
# ╠═b8647d13-c400-4324-85f7-8994a1f2322f
# ╠═dec72eaf-d4c8-4f67-9e4e-8d71c6651796
# ╠═e1b1c9aa-cff4-4155-80c2-9181be546116
# ╟─20fa79ac-93c3-4e3f-bf31-f9c0f0e79d4a
# ╠═1f8f4a15-0eca-4759-bee6-843306e35754
# ╠═8ca8e4e7-6450-4aaf-b208-6c0aa95ccb12
# ╠═f6f7081b-212a-4c7a-b91c-0d747f4b9fb2
# ╠═d7160016-c046-4ab8-b349-27f7a209d853
# ╠═10dcdca1-71fe-42eb-8f04-6122b3a03666
# ╠═f85ec710-1f59-4826-a376-6672bf0552ae
# ╠═9b7e4278-2f1a-4fff-a6f7-2249f45aaccf
# ╠═40d953af-e5e8-42c0-be55-26d95e334195
# ╟─037cddff-3663-4fcb-8575-ea7a78e490b8
# ╟─d0b306c4-6cbc-48dc-90ba-8ef4498f8d73
# ╠═abe311d5-fbd1-40e2-8bce-2c341301deef
# ╠═88f1bf59-5eb4-4e2d-b28b-0d9b1004d5bd
# ╠═be2cb230-df60-4db6-86e3-8f7d91e28a86
# ╠═36c79f63-7a8c-46c9-afa5-95bbed8fd598
# ╠═92490cef-70af-417e-8bc8-91b1be0635cc
# ╟─afc32ff2-014a-43ee-8e6a-d35a57434622
# ╠═03f1778e-c222-4221-a4e7-eabf1071d298
# ╟─f5f3c110-b31f-4c41-8eb5-121965a5e54d
# ╠═94380193-9789-43db-8220-5ef210dc36aa
# ╠═949b08e6-f84b-475d-ab80-66ecc7532b8e
# ╠═40702b4e-cecc-4a2c-ad88-d2613f49f71b
# ╠═da406d8f-a48b-4a56-8f34-73f0ab837589
# ╟─307a27e6-8128-42be-9ba8-cbdacb498ada
# ╟─4e97e67b-8ef6-412d-8851-2b6635ff46e6
# ╠═88584373-9519-46fe-ae88-7199e76fb8f7
# ╟─1215b318-8121-4adb-aa5d-139654607717
# ╠═b9180215-218a-48d8-a662-2a5664c2e480
# ╠═209cad48-c496-41b4-8093-6e2a55cc469f
# ╠═905dd1d7-5f1b-43b3-a461-651cc992bd1c
# ╠═c5d305a9-5a27-44d6-a610-50ab6d396f93
# ╟─82ac21a4-384e-47e8-a964-95c8eb0889bd
# ╠═139debf3-10ad-48da-9c2d-a13033486fbf
# ╟─d429eede-d77d-4709-b635-c3993fce4a47
# ╠═ff5e3103-11cb-4438-aa19-c25af84fc3da
# ╠═ec755f39-5dba-4498-9760-07629ad00dc7
# ╟─c94ccf15-a0bd-4605-96f2-a27d0508fbeb
# ╠═65f44437-4d44-4d93-a9bc-814ab9b6ee02
# ╟─bf8fd692-0e01-4368-a52d-91edb4fe77c0
# ╠═d58c4727-448d-47b2-a194-2002304a4708
# ╟─4035ef81-4348-4721-9448-48c569a65c9d
# ╠═ac2caf5e-9874-44b9-ba5e-f0011940055f
# ╠═b990e038-aff2-4d44-bb27-eb9de8c3977a
# ╠═c1acfea2-a60d-4433-a00b-6d3515274a18
# ╟─b573a376-d60f-4f1d-b876-592dbcd47be4
# ╠═03679ace-1368-4be3-988c-1f0dcaa1407e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
