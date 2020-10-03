### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ a481f540-027d-11eb-10b5-5bd4e827e274
md"""
# Uvod

__Numerička analiza je znanost o računanju rješenje problema koji su matematički postavljeni u polju relanih ili kompleksnih brojeva.__  

Navedimo dva primjera.
"""

# ╔═╡ c2a8c4e0-027d-11eb-26fe-1f7bf12091f4
md"""
## Računanje $\pi \approx 3.14159265358979$

Broj $\pi$ je omjer opsega kružnice i njenog promjera. 
"""

# ╔═╡ eab6ace0-027d-11eb-0f82-db588f0a253b
md"""

__Arhimedova metoda__

Upišite pravilni $n$-terokut u kružnicu radijusa $1$. Izračunajte opseg njegove gornje polovice.
To je najlakše napraviti kada je $n=2^k$. Koristeći geometriju,

\begin{equation}
p(n) = 2n \sin \: \frac{\pi}{2n}.
\end{equation}

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

\begin{align}
\sin \: \frac{\theta}{2} &= \sqrt{\frac{1-\cos \: \theta }{2}}, \\\\
\cos \: \theta &= \sqrt{1-\sin^2 \: \theta } .
\end{align}

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
steps=5

# ╔═╡ 4d59b60e-027b-11eb-02bd-5d1574308999
s=1.0

# ╔═╡ 4d66d570-027b-11eb-38a6-5d4f8eff3322
c=0.0

# ╔═╡ a47b53a0-027f-11eb-256d-f7d5ce0f97e3
for n=1:steps
	c=0.0
    s=sqrt((1-c)/2)
    c=sqrt(1-s^2)
	println("2n = ", 2^(n+1), ", približna vrijednost od π = ",2^(n+1)*s)
end

# ╔═╡ 4d68f850-027b-11eb-14a1-1741ea4ce7c2
"točan π", π

# ╔═╡ 6ac49527-05f0-4a0a-ad99-aba5d4b07c17
md"""
__Zadatak:__ Povečajte broj koraka i objasnite što se događa.
"""

# ╔═╡ 8443676d-70db-4f63-8b8a-a7e96175f813
md"""
Kasnije ćemo opisati modernije poboljšanje koje je i brže.
Uz $h=1/(2n)$ vrijedi

$$
p(n) = \frac{\sin \: \pi h}{h} = \pi -a_2 h^2 + a_4 h^4 - \cdots,
$$

gdje je $a_k = \displaystyle\frac{\pi^{k+1}}{(k+1)!}$. Rakle, ova formula teži k $\pi$ po stopi približno jednakoj
$O(h^2)$. To je _aproksimacijski problem_. Ne postoji konačan algoritam za računanje broja $\pi$, jer se radi o transcedentnom broju (iracionalan i nije korijen niti jednog polinoma s cjelobrojnim koeficijentima). 
Međutim, možemo ga aproksimirati po volji točno.
"""

# ╔═╡ bd7a0ad8-8b31-42cf-94e9-881711a3c40e
md"""
## Kvdratna jednadžba

Sljedeći problem ilustrira potpuno različite fenomene. 

Izračunajmo korijene $p(x) = ax^2 + bx +c =0$ za realne konstante $a$, $b$ i $c$.

Rješenja su

$$
x_{1,2} = \frac{ -b \pm \sqrt{b^2 - 4ac}}{2a}. \tag{1}
$$

Kako ćemo pristupiti izradi programa za računanje ovih korijena?

Trebamo riješiti pet posebnih slučajeva.

__Slučaj I.__ $a=0, b \neq 0$.

Ovo više nije kvadratna jednadžba, nego linearna. Jedino rješenje je

$$
x_1 = -c/b.
$$

__Slučaj II.__ $a=b=0$

Ako je $c \neq 0$, nema rješenja. Ako je $c=0$, bilo koji $x$ je rješenje.

Slučajevi koje smo rješavali u srednjoj školi koriste diskriminantu
$b^2 - 4ac$.

__Slučaj III.__ $b^2 -4ac < 0$. Dva kompleksna rješenje (nisu realna).

$$
x_{1,2} = -\frac{b}{2a} \pm \mathbf{i} \frac{\sqrt{4ac-b^2}}{2a}, \quad \mathbf{i}^2 = -1.
$$

__Slučaj IV.__ $b^2 - 4ac =0$. Jedan dvostruki korijen (realni)

$$
x_1 = x_2 = -\frac{b}{2a}.
$$

__Slučaj V.__ $b^2 -4ac > 0$. Dva različita realna rješenje.
Koristimo formulu (1)?
"""

# ╔═╡ b5fe78a7-5258-41ba-9544-c0cd64496ef5
a=1
b=2
c=10.0^(-17)
x₁=(-b-sqrt(b*b-4*a*c))/(2*a)
x₂=(-b+sqrt(b*b-4*a*c))/(2*a)

x₁,x₂

# ╔═╡ bca1eb35-dd3d-41af-8130-a8035551ec49
md"""
Dva realna rješenja su (prvih 17 znamenki)

$$
x_1 = -2, \quad x_2 = -5 \cdot 10^{-18}.
$$

Gornji algoritam izračuna $x_1$ točno, ali $x_2 = 0$. Standardni brojevi s plivajućim zarezom (floating-point) u dvostrukoj točnosti, `Float64`, spremaju približno $15$ decimalnih znamenki (54 binarne znamenke) i u tih 15 znamenki je $\sqrt{b^2 -4ac}-b =0$. To je zato što oduzimamo dva bliska broja, a jedan od njih je približno točan, pa je razlika isključivo "greške zaokruživanja". Jednostavno zapažanje rješava ovaj problem.

"Veliki" korijen kvadratne jednadžbe u __slučaju V__ je

$$
x_1 = \frac{ -b - \mathrm{sign}(b) \sqrt{b^2 - 4ac}}{2a},
$$

a dva korijena zadovoljavaju 

$$
x_1 x_2 = \frac{c}{a}.
$$

Primjetite da, osim unutar kvadratnog korijena, zbrajamo brojeve istog predznaka!
Nakon nekoliko transformacija imamo formulu za sitni korijen

$$
x_2 = \frac{c}{a x_1} = \frac{-2c}{ b + \mathrm{sign}(b) \sqrt{b^2 - 4ac}}.
$$

Koristeći ovu formulu oba korijena računamo blizu točnosti stroja.

U ovom primjeru imali smo egzaktnu formulu, ali u aritmetici plivajućeg zareza standarna formula daje rezultate koji se značajno razlikuju od rezultata u realnoj aritmetici.
"""

# ╔═╡ 2f084ae2-1e38-4c0b-83b6-da9f025565b1
x₁=(-b-sign(b)*sqrt(b*b-4*a*c))/(2*a)
x₂=-2*c/(b+sign(b)*sqrt(b*b-4*a*c))
x₁,x₂

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
            ximaginarni=sqrt(-Δ)/(2*a)
            xrealni=-b/(2*a)
            x₁=xrealni+im*ximaginarni
            # x₂ je kompleksno konjugirani x₁, 
            # x₂ = xrealni - im*ximaginarni
            x₂=conj(x₁)
        else
            if b==0
                # Julia lako računa s kompleksnim brojevima,
                # pa u ovom slučaju možemo koristiti formulu.
                x₁=sqrt(-c)/a
                x₂=-x₁
            else
                # Slučaj s dva različita realna korijena.
                x₁=(-b-sign(b)*sqrt(Δ))/(2*a)
                x₂=-2*c/(b+sign(b)*sqrt(Δ))
            end
        end
        poruka="korijeni su dobri"
    end
    x₁, x₂, poruka
end

# ╔═╡ cd29c4bc-9cfa-40bc-bc80-74ee5a635a94
korijeni(1,0,9)

# ╔═╡ af57717b-4dca-4b54-8272-e98fe76fd845


# ╔═╡ Cell order:
# ╟─a481f540-027d-11eb-10b5-5bd4e827e274
# ╟─c2a8c4e0-027d-11eb-26fe-1f7bf12091f4
# ╟─eab6ace0-027d-11eb-0f82-db588f0a253b
# ╟─3d725c70-027b-11eb-10cb-5b6e56311aa1
# ╠═3dea5041-85b9-4d8b-9dcf-2d39c8481314
# ╠═4d59b60e-027b-11eb-02bd-5d1574308999
# ╠═4d66d570-027b-11eb-38a6-5d4f8eff3322
# ╠═a47b53a0-027f-11eb-256d-f7d5ce0f97e3
# ╠═4d68f850-027b-11eb-14a1-1741ea4ce7c2
# ╟─6ac49527-05f0-4a0a-ad99-aba5d4b07c17
# ╟─8443676d-70db-4f63-8b8a-a7e96175f813
# ╟─bd7a0ad8-8b31-42cf-94e9-881711a3c40e
# ╠═b5fe78a7-5258-41ba-9544-c0cd64496ef5
# ╟─bca1eb35-dd3d-41af-8130-a8035551ec49
# ╠═2f084ae2-1e38-4c0b-83b6-da9f025565b1
# ╟─f8cf034f-b149-4ff9-8a7e-bc8143581d51
# ╠═0738969f-7c4f-4f8e-92de-62cb0cd6510c
# ╠═cd29c4bc-9cfa-40bc-bc80-74ee5a635a94
# ╠═af57717b-4dca-4b54-8272-e98fe76fd845
