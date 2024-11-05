### A Pluto.jl notebook ###
# v0.19.46

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

__Mantisa__ (_significand_) $d$ je oblika

$$
\begin{aligned}
	d &= d_0.d_1 \dots d_t = d_0+d_1 \beta^{-1} + d_2 \beta^{-2}
	+ \dots + d_t \beta^{-t}\\
d_i  &\in \{ 0,1\}\\
	d_0 &= 1 \qquad \mbox{ normalizirani oblik }   \\
	d_0 &= 0 \qquad \mbox{ denormalizirani oblik }   \\
\end{aligned}$$

__Decimalni dio__ (_fraction_) je $f=d_1 d_2 \dots d_t$.

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
bitstring(Int32(1))

# ╔═╡ 516e1922-4de4-4f85-a35e-313607920ff4
bitstring(0.0)

# ╔═╡ 8ca8e4e7-6450-4aaf-b208-6c0aa95ccb12
bitstring(-2)

# ╔═╡ f6f7081b-212a-4c7a-b91c-0d747f4b9fb2
bitstring(-0.0)

# ╔═╡ d7160016-c046-4ab8-b349-27f7a209d853
bitstring(0.0)

# ╔═╡ 10dcdca1-71fe-42eb-8f04-6122b3a03666
bitstring(1.0)

# ╔═╡ 290ebd38-0aca-4392-8b9c-90e4271ef19c
parse(Int, "1111111111"; base=2)

# ╔═╡ f85ec710-1f59-4826-a376-6672bf0552ae
bitstring(Float32(1.0))

# ╔═╡ ebad923b-3afb-4dd6-9e74-22eb304c0dfc
parse(Int, "1111111"; base=2)

# ╔═╡ 36d8dc90-3719-42b8-99e6-f53e05fd8695
md"
Dakle, $1.0=1\cdot 2^0$, pa binarni broj $1111111_2=127_{10}$ označava eksponent $0$.
"

# ╔═╡ 9b7e4278-2f1a-4fff-a6f7-2249f45aaccf
bitstring(Float32(2.0))

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

# ╔═╡ efb8122d-bbde-48ea-bf83-278cdc2ffe67
a₁/c₁==a₁/b₁

# ╔═╡ 2dfc2812-be24-4ee2-a619-8f2e3ad08ede
a₁/c₁

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

| Eksponent | Decimalni dio | Prikazuje |
| :-----    | :-----  | :-----    |
| $e=e_{\min}-1$ |  $f=0$     | $\pm 0$     |
| $e=e_{\min}-1$ |  $f\neq 0$ | $0.f\times 2^{e_\min}$ - denormalizirani brojevi |
| $e_{\min}\leq e\leq e\_\max$ |  $f$      | $1.f \times 2^e$ - standardni brojevi | 
| $e=e_{\max}+1$ |  $f=0$      |  $\pm$`Inf`     |
| $e=e_{\max}+1$ |  $f\neq 0$  |  `NaN`     |

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

# ╔═╡ 722ad758-0a95-4a1e-8600-86162ed8d319
b₂^2-32

# ╔═╡ f65bad48-a9ec-412e-96cc-d5aa82d9f26a
sqrt(b₂^2-32)

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
PlutoUI = "~0.7.54"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "2715a914af8a023ee857a2c9015593fe036c0e1b"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

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

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

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
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

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

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

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
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

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
version = "1.11.0"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

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
# ╠═516e1922-4de4-4f85-a35e-313607920ff4
# ╠═8ca8e4e7-6450-4aaf-b208-6c0aa95ccb12
# ╠═f6f7081b-212a-4c7a-b91c-0d747f4b9fb2
# ╠═d7160016-c046-4ab8-b349-27f7a209d853
# ╠═10dcdca1-71fe-42eb-8f04-6122b3a03666
# ╠═290ebd38-0aca-4392-8b9c-90e4271ef19c
# ╠═f85ec710-1f59-4826-a376-6672bf0552ae
# ╠═ebad923b-3afb-4dd6-9e74-22eb304c0dfc
# ╟─36d8dc90-3719-42b8-99e6-f53e05fd8695
# ╠═9b7e4278-2f1a-4fff-a6f7-2249f45aaccf
# ╟─037cddff-3663-4fcb-8575-ea7a78e490b8
# ╟─d0b306c4-6cbc-48dc-90ba-8ef4498f8d73
# ╠═abe311d5-fbd1-40e2-8bce-2c341301deef
# ╠═efb8122d-bbde-48ea-bf83-278cdc2ffe67
# ╠═2dfc2812-be24-4ee2-a619-8f2e3ad08ede
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
# ╠═722ad758-0a95-4a1e-8600-86162ed8d319
# ╠═f65bad48-a9ec-412e-96cc-d5aa82d9f26a
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
