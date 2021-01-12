### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ 7d2eceff-5dc8-4088-8504-8fa469948983
using QuadGK, FastGaussQuadrature, LinearAlgebra

# ╔═╡ 456ae14f-d92c-4561-bcfd-d1b9849b8244
using FFTW

# ╔═╡ c42b0321-1003-4935-857e-a8422fabd485
md"""

# Numerička integracija

## Newton-Cotesove formule

Funkciju $f(x):[a,b]\to\mathbb{R}$ interpoliramo polinomom stupnja $n$ kroz $n+1$ ravnomjerno raspoređenu točku te integral aproksimiramo integralom interpolacijskog polinoma. Polinom $P_n(x)$ možemo računati u Lagrangeovom obliku  
(vidi bilježnicu [NA09 Interpolacijski polinomi.jl](https://ivanslapnicar.github.io/Numericka_analiza/NA09%20Interpolacijski%20polinomi.jl)):

$$
L_k(x)=\prod_{{i=0}\atop {i\neq k}}^n \frac{x-x_i}{x_k-x_i},$$

Tada je 

$$
f(x)\approx P_n(x)=\sum_{k=0}^n f(x_k) L_k(x),$$

pa je 

$$
\int_a^b f(x)\, dx\approx \int_a^b P_n(x) \, dx=\sum_{k=0}^n f(x_k) \int_a^b L_k(x)\, dx =(b-a)\sum_{k=0}^n \omega_k f(x_k). \qquad (1)$$

Vrijedi 

$$x_i=a+\displaystyle\frac{b-a}{n}\, i.$$

Uz supstituciju $x=a+(b-a)\,t$ vrijedi

$$
\frac{x-x_i}{x_k-x_i}=\frac{nt-i}{k-i}$$

pa su __težine__ $\omega_k$ jednake

$$
\omega_k=\frac{1}{b-a}\int_a^b L_k(x)\, dx = \int_0^1 \prod_{{i=0}\atop {i\neq k}}^n \frac{nt-i}{k-i} \, dt. \qquad(2)$$
"""

# ╔═╡ 5071e7f3-d353-4429-8885-ca9431edf5a9
md"""
### Trapezna formula

Za $n=1$ formula (2) daje 

$$\omega_0=\omega_1=\frac{1}{2}$$

pa Newton-Cotesova formula (1) daje

$$\int_a^b f(x)\, dx\approx \int_a^b P_1(x) \, dx=\frac{b-a}{2}(f(a)+f(b)),$$

što je površina trapeza s vrhovima $(a,0)$, $(a,f(a))$, $b,f(b))$ i $(b,0)$.

Točno značenje $\approx$ daje nam sljedeći teorem:

__Teorem.__ Ako je funkcija $f$ neprekidna na intervalu $[a,b]$, onda postoji $c\in(a,b)$ takav da je

$$\int_a^b f(x)\, dx= \frac{b-a}{2}(f(a)+f(b))-\frac{(b-a)^3}{12}f''(c).$$

_Dokaz:_ Prema Teoremu iz bilježnice [NA10 Interpolacija funkcija.jl](https://ivanslapnicar.github.io/Numericka_analiza/NA10%20Interpolacija%20funkcija.jl.html) za $c\in(a,b)$ vrijedi 

$$
E=\int_a^b f(x)\, dx -\int_a^b P_1(x) \, dx = \int_a^b (f(x)- P_1(x)) \, dx = \int_a^b \frac{f''(c)}{2}(x-a)(x-b)\, dx.$$

Uz supstituciju $x=a+(b-a)\,t$ vrijedi

$$
E=\frac{f''(c)}{2} (b-a)^3 \int_0^1 t(t-1)\, dt= -\frac{(b-a)^3}{12}f''(c)$$

i teorem je dokazan. 

Sada interval $[a,b]$ podjelimo na $n$ jednakih podintervala,

$$[x_{i-1},x_{i}],\quad  i=1,2,\ldots,n,$$ 

i uvedimo oznake

$$
\Delta x=\frac{b-a}{n}, \quad y_i=f(x_i).$$


Na svakom podintervalu vrijedi 

$$\int_{x_{i-1}}^{x_i} f(x)\, dx= \frac{\Delta x}{2}(y_{i-1}+y_i)
-\frac{(\Delta x)^3}{12}f''(c_i),\quad c_i\in[x_{i-1},x_i].$$

Zbrajanje daje

$$
\int_a^b f(x)\, dx=I_n+E_n,$$ 

gdje je $I_n$ __trapezna formula__,

$$
I_n=\Delta x\bigg( \frac{y_0}{2} +y_1+y_2+\cdots +y_{n-1}+\frac{y_n}{2}\bigg),$$

a za __pogrešku__ $E_n$ vrijedi

$$
E_n =-\frac{(\Delta x)^3}{12}\sum_{i=1}^n f''(c_i)=
-\frac{b-a}{12}(\Delta x)^2 f''(c), \quad c\in[a,b].$$

Koristili smo činjenicu da zbog neprekidnosti funkcije $f''$ na intervalu $(a,b)$ postoji točka $c\in[a,b]$ za koju vrijedi 

$$
\frac{1}{n}\sum_{i=1}^n f''(c_i)=f''(c).$$

Također, vrijedi 

$$
|E_n|\leq \frac{b-a}{12}(\Delta x)^2 \max_{x\in(a,b)} |f''(x)|. \tag{3}$$


__Napomena__. Izvod trapezne formule i ocjene pogreške dan je u knjigama
[Numerička matematika, poglavlje 7.1](http://www.mathos.unios.hr/pim/Materijali/Num.pdf) i [Matematika 2, poglavlje 2.7.2](http://lavica.fesb.hr/mat2/predavanja/node42.html).

"""

# ╔═╡ 04f1b0f8-10bc-4aaf-bbf2-c411d740db0e
md"""
### Simpsonova formula

Za $n=2$ Newton-Cotesova formula (1) daje

$$\omega_0=\frac{1}{6},\quad \omega_1=\frac{2}{3},\quad \omega_2=\frac{1}{6}.$$

Interval $[a,b]$ podjelimo na paran broj $n$ jednakih podintervala, na svaki podinterval

$$[x_{2i-1},x_{2i+1}],\quad i=1,2,\ldots,\frac{n}{2},$$

primijenimo  Newton-Cotesovu formulu i zbrojimo, što daje __Simpsonovu formulu__:

$$
I_n=\frac{\Delta x}{3}\big(y_0 +4(y_1+y_3\cdots +y_{n-1})+2(y_2+y_4+\cdots+y_{n-2})+y_n\big).$$

Vrijedi

$$
\int_a^b f(x)\, dx =I_n+E_n,$$

pri čemu je __pogreška__ $E_n$ omeđena s

$$
|E_n|\leq \frac{b-a}{180}(\Delta x)^4 \max_{x\in(a,b)} |f^{(4)}(x)|.\tag{4}$$

Za detalje vidi knjige [Numerička matematika, poglavlje 7.3](http://www.mathos.unios.hr/pim/Materijali/Num.pdf) i [Matematika 2, poglavlje 2.7.3](http://lavica.fesb.hr/mat2/predavanja/node43.html).

"""

# ╔═╡ 7a11ccd0-68c1-46d3-9249-3cc9f867ef65
md"""
### Richardsonova ekstrapolacija

Ocjena pogreške pomoću formula (3) i (4) može biti složena. __Richardsonova ekstrapolacija__ nam omogućava da, uz određene uvjete, pogrešku procijenimo koristeći aproksimaciju integrala s $n/2$ točaka.
Ako se u ocjeni pogreške javlja član $(\Delta x)^m$, onda je pogreška približno manja od broja (vidi [Matematika 2, poglavlje 2.7.4](http://lavica.fesb.hr/mat2/predavanja/node44.html))

$$
E=\frac{\big(\frac{n}{2}\big)^m}{n^m-\big(\frac{n}{2}\big)^m}(I_n-I_{n/2}).$$

Predznak broja $E$ daje i predznak pogreške, odnosno, ako je $E>0$, tada je približno

$$
\int_a^b f(x)\, dx\in[I_n,I_n+E],$$

a ako je $E\leq 0$, tada je približno

$$
\int_a^b f(x)\, dx\in[I_n+E,I_n].$$

Kod dokazivanja ovih formula pretpostavljamo da postoji broj $\omega$ takav da je 

$$
\displaystyle I=I_n+E_n, \qquad E_n=\omega\,  (\Delta x)^m \tag{5}$$

za svaki $\Delta x$. Ova pretpostavka nije uvijek točno ispunjena, no u velikom broju slučajeva možemo smatrati da ona vrijedi. Na primjer, kod trapezne formule možemo uzeti 

$$\omega=-\frac{b-a}{12}\max_{x\in(a,b)}f''(x)$$

pa pretpostavka (5) znači da možemo uzeti (približno) isti $\omega$ za različite vrijednosti od $\Delta x$. Koristeći pretpostavku imamo

$$
\begin{aligned}
E_{n}&=I-I_n =\omega (\Delta x)^m =\omega \, \left(\frac{b-a}{n}\right)^m,\cr
E_{n/2}&=I-I_{n/2} =\omega (2\Delta x)^m =\omega\,  \left(\frac{b-a}{\frac{n}{2}}\right)^m.
\end{aligned}$$

Dakle,

$$
\displaystyle I_{n}-I_{n/2}=\omega \, (b-a)^m \left(\frac{1}{(n/2)^m}-\frac{1}{n^m}\right)$$

pa je

$$
\displaystyle \omega = \frac{n^m (n/2)^m}{(b-a)^m} \frac{I_{n}-I_{n/2}}{n^m-(n/2)^m},$$

što konačno daje

$$
\displaystyle E_{n}\approx \frac{(n/2)^m} {n^m-(n/2)^m} (I_{n}-I_{n/2})=E.$$

"""

# ╔═╡ e9054e4d-0e5a-4b26-a323-895ec36ddf73
function Trapez(f::Function,a::Number,b::Number,n::Int64)
    # n je broj intervala
    m=2
    X=range(a,stop=b,length=n+1)
    Y=map(f,X)
    Δx=(b-a)/n
    I=Δx*(Y[1]/2+sum(Y[2:end-1])+Y[end]/2)
    # Richardsonova ekstrapolacija
    I₂=2*Δx*(Y[1]/2+sum(Y[3:2:end-2])+Y[end]/2)
    E=(n/2)^m*(I-I₂)/(n^m-(n/2)^m)
    I,E
end 

# ╔═╡ 0283137b-539e-43c2-b2e4-34f1890d78b9
function Simpson(f::Function,a::Number,b::Number,n::Int64)
    # n je broj intervala, djeljiv s 4
    m=4
    X=range(a,stop=b,length=n+1)
    Y=map(f,X)
    Δx=(b-a)/n
    I=Δx/3*(Y[1]+4*sum(Y[2:2:end-1])+2*sum(Y[3:2:end-2])+Y[end])
    # Richardsonova ekstrapolacija
    I₂=2*Δx/3*(Y[1]+4*sum(Y[3:4:end-2])+2*sum(Y[5:4:end-4])+Y[end])
    E=(n/2)^m*(I-I₂)/(n^m-(n/2)^m)
    I,E
end 

# ╔═╡ 4baefc88-ba8b-43f7-8cec-48bf4d80c314
md"""
### Eliptički integral

Izračunajmo opseg elipse s polu-osima $2$ i $1$ (vidi  [Matematika 2, poglavlje 2.7.1](http://lavica.fesb.hr/mat2/predavanja/node41.html)). Parametarske jednadžbe elipse glase

$$
x=2\cos t,\quad y=\sin t,\quad t\in[0,\pi/2],$$

pa je četvrtina opsega jednaka

$$
\frac{O}{4}\int\limits_0^{\pi/2} \sqrt{(-2\sin t)^2+(\cos t)^2}\, dt
=2\int\limits_0^{\pi/2} \sqrt{1-\frac{3}{4}(\cos t)^2}\, dt.$$

Radi se o eliptičkom integralu druge vrste koji nije elementarno rješiv, ali se može naći u [tablicama](http://nvlpubs.nist.gov/nistpubs/jres/50/jresv50n1p43_A1b.pdf).

Vidimo da je 

$$O=8\int\limits_0^{\pi/2} \sqrt{1-\frac{3}{4}(\cos t)^2}\, dt \approx 8\cdot 1.21125.$$

"""

# ╔═╡ d0135ecd-9e39-4233-ab0d-74ac85ad4e2f
begin
	f₁(x)=sqrt(1-(3.0)/4*cos(x)^2)
	Trapez(f₁,0,π/2,4)
end

# ╔═╡ c0ba2cf8-7751-4bf2-8fd0-c4bb69e95504
Trapez(f₁,0,pi/2,10)

# ╔═╡ e370021d-2f28-4ec2-85c6-6bb094d242d1
Trapez(f₁,0,pi/2,24)

# ╔═╡ faab641c-5a27-4e36-a4dd-c8e82272eae9
Simpson(f₁,0,π/2,4)

# ╔═╡ c58f0660-6dd6-4045-bb33-a74351c48f73
Simpson(f₁,0,π/2,16)

# ╔═╡ a7cb847f-bc7d-4299-9b13-37a698a1f28b
Simpson(f₁,0,π/2,24)

# ╔═╡ aa07f01e-3bb9-4b47-b197-aba3c1310a2c
md"""
### Broj $\pi$

Vrijedi 

$$
\int_0^1 \frac{4}{1+x^2}\, dx=\pi.$$

Aproksimirajmo $\pi$ numeričkom integracijom i provjerimo pogrešku (vidi [Numerička matematika, poglavlje 7.3](http://www.mathos.unios.hr/pim/Materijali/Num.pdf)).

Pomoću trapezne formula možemo dobiti najviše pet točnih decimala. 
Simpsonova formula je točnija, ali je konvergencija spora.
"""

# ╔═╡ 311b6a2a-021b-4fc2-b9b0-5bafa0b1d4a6
# Za provjeru
BigFloat(π)

# ╔═╡ 1c13d010-afbf-45b4-8210-a7cfea8641a4
begin
	f₂(x)=4/(1+x^2)
	myπ=Trapez(f₂,0,1,10)
end

# ╔═╡ 94267bb3-624e-49b4-bbc8-2c24aacbd095
myπ[1]-BigFloat(π)

# ╔═╡ cb4df66a-1b00-403f-b49b-7a275dc2929f
myπ₁=Trapez(f₂,0,1,100)

# ╔═╡ 5cb0017b-1722-4f38-9f53-76ff788fce9e
myπ₁[1]-BigFloat(π)

# ╔═╡ d1b0fca4-b727-4670-ad25-3545dfb62629
myπ₂=Simpson(f₂,0,1,16)

# ╔═╡ 29e800e1-b38e-4dbf-a1f1-6084b35ef762
myπ₂[1]-BigFloat(π)

# ╔═╡ ad0d05f0-8e66-4410-956b-bf6802ada45a
myπ₃=Simpson(f₂,0,1,64)

# ╔═╡ 11c3d407-0923-4d65-85a2-b30911fd75d2
myπ₃[1]-BigFloat(π)

# ╔═╡ e7b2a060-7ae9-4418-a434-1693f73d04ac
md"""
## Gaussova kvadratura

Slično kao u formuli (1), integral aproksimiramo sumom umnožaka vrijednosti funkcije i odgovarajućih težina:

$$
\int_{a}^b \omega(x) f(x)\, dx=\sum_{k=1}^n \omega_k f(x_k),$$

gdje je $\omega(x)$ __težinska funkcija__.

Točke $x_k$ su nul-točke odgovarajućeg ortogonalnog polinoma $P_{n}(x)$ reda $n$, na primjer, __Legendreovog polinoma__ za $[a,b]=[-1,1]$ i $\omega(x)=1$ ili __Čebiševljevog polinoma__ za 
$[a,b]=[-1,1]$ i $\omega(x)=\displaystyle\frac{1}{\sqrt{1-x^2}}$ (vidi bilježnicu [NA12 Ortogonalni polinomi.ipynb](https://nbviewer.jupyter.org/github/ivanslapnicar/Numericka_analiza/blob/master/src/Jupyter/NA12%20Ortogonalni%20polinomi.ipynb)).

Težine su jednake

$$
\omega_k=\int_a^b \omega(x) \prod_{{i=1}\atop {i\neq k}}^n\frac{x-x_i}{x_k-x_i} \, dx.$$

__Pogreška__ je dana s

$$
E=\frac{f^{(2n)}(\xi)}{(2n)!}\int_a^b \omega(x) P_n^2(x)\, dx.$$

Za detalje vidi [Numerical Analysis, poglavlje 7.3](https://books.google.hr/books?id=kPDtAp3UZtIC&hl=hr&source=gbs_book_other_versions).

__Napomena__. Legendreovi i Čebiševljevi polinomi su definirani na intervalu $[-1,1]$ pa koristimo supstituciju

$$
\int_{a}^b \omega(x) f(x)\, dx = \frac{b-a}{2} \int_{-1}^1 \omega\bigg(\frac{b-a}{2}t+\frac{a+b}{2}\bigg) f\bigg(\frac{b-a}{2}t+\frac{a+b}{2}\bigg) dt.$$
"""

# ╔═╡ 656b359d-9504-4fa3-a327-d0977e10d642
mapnodes(x,a,b)=(b-a)*x/2 .+(a+b)/2

# ╔═╡ 427bfd71-5f3d-4b0b-b6f8-3753b5844a87
md"""
### Postojeće rutine

Profesionalne rutine za numeričku integraciju su složene, a većina programa ima ugrađene odgovarajuće rutine. Tako, na primjer, 

* Ṁatlab ima rutinu `quad` koja koristi adaptivnu Simpsonovu formulu, a 
* Julia u paketu [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) ima rutinu `quadgk()` i računa integral u $O(n^2)$ operacija.
* Julia također ima i paket  [`FastGaussQuadrature.jl`](https://github.com/ajt60gaibb/FastGaussQuadrature.jl) koji brzo računa točke i težine za zadani $n$ i razne težinske funkcije pa se pomoću točaka i težina lako izračuna integral u $O(n)$ operacija.
"""

# ╔═╡ cb4a8760-9ba7-4bb2-beea-9876a680829f
# ?quadgk

# ╔═╡ 4feab23e-b51b-4818-8343-a73deb67c942
# 1/8 opsega elipse
quadgk(f₁,0,π/2)

# ╔═╡ fa5a8850-c845-4328-94ca-cef66b4b52dd
# Broj π
quadgk(f₂,0,1)

# ╔═╡ ac14acbb-5e64-44ca-9922-68bd3a45b0a2
# Granice mogu biti i beskonačne
quadgk(x->exp(-x),0,Inf)

# ╔═╡ 4c7f59a4-0e11-4dd1-8d08-3c12b25223f5
varinfo(FastGaussQuadrature)

# ╔═╡ 5fbe0581-65b8-44b4-94f5-35950d579c92
# Na primjer
methods(gausschebyshev)

# ╔═╡ b9a58890-f3a2-49bb-836c-95a2971c9da7
gausschebyshev(16)

# ╔═╡ 4ae529de-6df8-4d38-917f-7b972904204d
# Sada računajmo integrale. U našem slučaju je ω(x)=1 
# pa nam treba Legendreov polinom.
ξ,ω=gausslegendre(32)

# ╔═╡ 988b6a6c-7732-47cd-915d-a26f44e1b3f0
# 1/8 opsega elipse, a=0, b=π/2
(π/2-0)/2*dot(ω,map(f₁,mapnodes(ξ,0,π/2)))

# ╔═╡ 5fa9e09b-569d-4ad4-8dd6-ffc93dc6cca9
# Broj π, a=0, b=1
(1-0)/2*dot(ω,map(f₂,mapnodes(ξ,0,1)))

# ╔═╡ 75fe9d8a-5e70-4d31-98f2-0ea82c511698
md"""
## Clenshaw-Curtisova kvadratura

Uz supstituciju $x=\cos\theta$, vrijedi

$$
I\equiv \int\limits_{-1}^1f(x)\, dx =\int\limits_0^\pi f(\cos\theta)\sin\theta \, d\theta.$$

Integral na desnoj strani se računa integriranjem Fourierovog reda parnog proširenja podintegralne funkcije:

$$
I\approx a_0+\sum_{k=1}^n \frac{2a_{2k}}{1-(2k)^2},$$

pri čemu se koeficijenti $a_k$ računaju formulom

$$
a_k=\frac{2}{\pi}\int\limits_0^\pi f(\cos\theta)\cos(k\theta)\,d\theta.$$

Koeficijenti $a_k$ se mogu računati numeričkom integracijom ili korištenjem brze Fourierove transformacije (FFT), što je puno brže.
Za detalje vidi [Clenshaw-Curtis Quadrature](https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature).

Ukoliko se integrira na intervalu $[a,b]$, prebacivanje u interval $[-1,1]$ vrši se kao u formuli (3).
"""

# ╔═╡ 3d03b422-76c1-4f9b-ad69-12c07f2d4403
function ClenshawCurtis(f::Function,a::Number,b::Number,n::Int64)
    # Implementacija koristeći numeričku integraciju
    z=Vector{Float64}(undef,n)
    g(x)=f(mapnodes(x,a,b))
    for i=1:n
        h(x)=g(cos(x))*cos(2*(i-1)*x)
        z[i]=2*quadgk(h,0,pi)[1]/pi
    end
    return (z[1]+2*sum([z[i]/(1-4*(i-1)^2) for i=2:n]))*(b-a)/2
end

# ╔═╡ cc0b62ab-821c-4b16-9c73-ec131176e296
# 1/8 opsega elipse
ClenshawCurtis(f₁,0,pi/2,8)

# ╔═╡ 6ab103b9-e138-4c57-8730-83fe6313bbb0
# Broj π
ClenshawCurtis(f₂,0,1,8)

# ╔═╡ 35dbcd72-9f1d-43a2-9134-4116eefacc97
 π

# ╔═╡ 307e81fa-5384-4aad-8fe2-7df9ab5bce0d
# "Nepravi" integral
ClenshawCurtis(x->exp(-x),0,1000,70)

# ╔═╡ 748ace02-d80f-4952-883b-de220083d58a
function ClenshawCurtisFFT(f::Function,a::Number,b::Number,n::Int64)
    # Brza implementacija pomoću fft(), 2^n je broj točaka
    g(x)=f(mapnodes(x,a,b))
    w=map(x->g(cos(x)),range(0,stop=2*pi,length=2^n))
    w[1]=(w[1]+w[end])/2
    z=real(fft(w))
    z/=2.0^(n-1)
    return (z[1]+2*sum([z[i]/(1-(i-1)^2) for i=3:2:2^(n-1)]))*(b-a)/2
end

# ╔═╡ 44406f97-1a7e-490b-9e75-3ab42e7f353f
ClenshawCurtisFFT(f₁,0,pi/2,4)

# ╔═╡ 41654b34-29c9-44d9-a0a2-d30a77d1b29c
ClenshawCurtisFFT(f₂,0,1,16),pi

# ╔═╡ 15905e38-d8f3-450a-9bcc-1d6b7a4fcae1
ClenshawCurtisFFT(x->exp(-x),0,1000,18)

# ╔═╡ a46784f5-e834-4671-9689-fd7ec3fada6b


# ╔═╡ Cell order:
# ╟─c42b0321-1003-4935-857e-a8422fabd485
# ╟─5071e7f3-d353-4429-8885-ca9431edf5a9
# ╟─04f1b0f8-10bc-4aaf-bbf2-c411d740db0e
# ╟─7a11ccd0-68c1-46d3-9249-3cc9f867ef65
# ╠═e9054e4d-0e5a-4b26-a323-895ec36ddf73
# ╠═0283137b-539e-43c2-b2e4-34f1890d78b9
# ╟─4baefc88-ba8b-43f7-8cec-48bf4d80c314
# ╠═d0135ecd-9e39-4233-ab0d-74ac85ad4e2f
# ╠═c0ba2cf8-7751-4bf2-8fd0-c4bb69e95504
# ╠═e370021d-2f28-4ec2-85c6-6bb094d242d1
# ╠═faab641c-5a27-4e36-a4dd-c8e82272eae9
# ╠═c58f0660-6dd6-4045-bb33-a74351c48f73
# ╠═a7cb847f-bc7d-4299-9b13-37a698a1f28b
# ╟─aa07f01e-3bb9-4b47-b197-aba3c1310a2c
# ╠═311b6a2a-021b-4fc2-b9b0-5bafa0b1d4a6
# ╠═1c13d010-afbf-45b4-8210-a7cfea8641a4
# ╠═94267bb3-624e-49b4-bbc8-2c24aacbd095
# ╠═cb4df66a-1b00-403f-b49b-7a275dc2929f
# ╠═5cb0017b-1722-4f38-9f53-76ff788fce9e
# ╠═d1b0fca4-b727-4670-ad25-3545dfb62629
# ╠═29e800e1-b38e-4dbf-a1f1-6084b35ef762
# ╠═ad0d05f0-8e66-4410-956b-bf6802ada45a
# ╠═11c3d407-0923-4d65-85a2-b30911fd75d2
# ╟─e7b2a060-7ae9-4418-a434-1693f73d04ac
# ╠═656b359d-9504-4fa3-a327-d0977e10d642
# ╟─427bfd71-5f3d-4b0b-b6f8-3753b5844a87
# ╠═7d2eceff-5dc8-4088-8504-8fa469948983
# ╠═cb4a8760-9ba7-4bb2-beea-9876a680829f
# ╠═4feab23e-b51b-4818-8343-a73deb67c942
# ╠═fa5a8850-c845-4328-94ca-cef66b4b52dd
# ╠═ac14acbb-5e64-44ca-9922-68bd3a45b0a2
# ╠═4c7f59a4-0e11-4dd1-8d08-3c12b25223f5
# ╠═5fbe0581-65b8-44b4-94f5-35950d579c92
# ╠═b9a58890-f3a2-49bb-836c-95a2971c9da7
# ╠═4ae529de-6df8-4d38-917f-7b972904204d
# ╠═988b6a6c-7732-47cd-915d-a26f44e1b3f0
# ╠═5fa9e09b-569d-4ad4-8dd6-ffc93dc6cca9
# ╟─75fe9d8a-5e70-4d31-98f2-0ea82c511698
# ╠═3d03b422-76c1-4f9b-ad69-12c07f2d4403
# ╠═cc0b62ab-821c-4b16-9c73-ec131176e296
# ╠═6ab103b9-e138-4c57-8730-83fe6313bbb0
# ╠═35dbcd72-9f1d-43a2-9134-4116eefacc97
# ╠═307e81fa-5384-4aad-8fe2-7df9ab5bce0d
# ╠═456ae14f-d92c-4561-bcfd-d1b9849b8244
# ╠═748ace02-d80f-4952-883b-de220083d58a
# ╠═44406f97-1a7e-490b-9e75-3ab42e7f353f
# ╠═41654b34-29c9-44d9-a0a2-d30a77d1b29c
# ╠═15905e38-d8f3-450a-9bcc-1d6b7a4fcae1
# ╠═a46784f5-e834-4671-9689-fd7ec3fada6b
