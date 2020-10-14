### A Pluto.jl notebook ###
# v0.12.3

using Markdown
using InteractiveUtils

# ╔═╡ 63842a33-b596-4f3f-bde5-7953ec4dc840
using Markdown

# ╔═╡ 53d8d8e0-0d7f-11eb-358f-3f001e5e7d78
using LinearAlgebra

# ╔═╡ 48167400-0d83-11eb-1c7d-359a2574c8b1
md"
# Gaussova eliminacija


## Općenito

Sustav $Ax=b$
se rješava u tri koraka (_bez pivotiranja_ ):

1. $A=LU\quad$ (LU rastav, $O(\frac{2}{3}n^3)$ operacija),
2. $Ly=b\quad$ (donje trokutrasti sustav, $n^2$ operacija),
3. $Ux=y\quad$ (gornje torkutasti sustav, $n^2$ operacija).

S pivotiranjem vrijedi

1. $PA=LU$
2. $Ly=P^T b$
3. $Ux=y$ 
"

# ╔═╡ 221d2474-de59-4042-918f-534305d8708f
md"""

## Primjeri

Sljedeći primjeri ukazuju na dva fenomena, jedan od kojih smo već vidjeli, dok drugi nismo. 
U ovom primjeru $\epsilon$ je vrijednost koju daje funkcija `eps()`.

Promotrimo sustav linearnih jednadžbi 

$\begin{aligned}\frac{\epsilon}{10} x_1 + x_2  = 1 \\ x_1 + x_2 = 2\end{aligned}$

Dobro približno rješenje je $x_1 = x_2 =1$. Koristimo proširenu matricu sustava:

\begin{align}
& \left(\begin{array}{cc|c} \displaystyle\displaystyle\frac{\epsilon}{10} & 1 & 1 \cr 1 & 1 & 2 \end{array}\right) \sim
\left(\begin{array}{cc|c} \displaystyle \displaystyle\frac{\epsilon}{10} & 1 & 1 \cr 0 & 1-\displaystyle\displaystyle\frac{10}{\epsilon} & 2-\displaystyle\displaystyle\frac{10}{\epsilon}\end{array}\right)
\approx \left(\begin{array}{cc|c}  \displaystyle\displaystyle\frac{\epsilon}{10} & 1 & 1 \cr 0 & -\displaystyle\displaystyle\frac{10}{\epsilon} & -\displaystyle\displaystyle\frac{10}{\epsilon}\end{array}\right).
\end{align}

Zadnja transformacija je zaokruživanje do na točnost stroja $\epsilon$.
Vrlo značajni "1" i "2" u zadnjem retku su _nestali prilikom zaokruživanja_ !

Rješavanje zadnjeg trokutastog sustava daje $x_1 = 0$, $x_2 = 1$. Uvrštavanje u originalni sustav daje $x_1+x_2 = 1$, što je "vrlo netočno".
"""

# ╔═╡ 89291e44-4ea6-4e74-99bd-a30e8ee4d895
[eps()/10 1;0 1-10/eps()]\[1; 2-10/eps()]

# ╔═╡ ed629240-6b9a-4a85-b412-97a06565a9cf
md"""
I ovdje postoji "rješenje" koje se zove _parcijalno pivotiranje_. Stavimo apsolutno najveći element koji još nije poništen u promatranom stupcu na pivotnu odnosno dijagonalnu poziciju:

\begin{align}
&\left(\begin{array}{cc|c}     1 & 1 & 2 \cr \displaystyle\frac{\epsilon}{10} & 1 & 1 \end{array}\right) \sim
\left(\begin{array}{cc|c}     1 & 1 & 2 \cr 0                   & 1-\displaystyle\frac{\epsilon}{10}&1-\displaystyle\frac{\epsilon}{10}\end{array}\right)
\approx   \left(\begin{array}{cc|c}     1 & 1 & 2 \cr 0                   & 1                    &1                    \end{array}\right)
\end{align}
"""

# ╔═╡ b8b17cd3-c828-438d-a2c0-13b1965ed778
[1 1;0 1-eps()/10]\[2; 1-eps()/10]

# ╔═╡ eb5aab00-0d73-11eb-2f45-771b2e23a5e3
md"""
Ovo je ispravno rješenje do na točnost stroja.

Ponekad promjena algoritma _ne može nikako pomoći_ ! Promotrimo sustav 

$\begin{aligned}
(1+2\epsilon)x_1 + (1+2\epsilon)x_2 &= 2 \cr
(1+\epsilon)x_1 + x_2 &=2
\end{aligned}$

čija proširena matrica glasi

$$\left(\begin{array}{cc|c}
(1+2\epsilon)&     (1+2\epsilon )&     2 \cr
(1+\epsilon)&     1   &2 \end{array}\right).$$

Pomnožimo prvi redak s 
$\alpha = (1+\epsilon)/(1+2\epsilon)= 1-\epsilon + O(\epsilon^2)$
i dodajmo drugom:

$$\left(\begin{array}{cc|c}
(1+2\epsilon)&     (1+2\epsilon )&     2 \cr
0           &     -\epsilon & 2\epsilon \end{array}\right).$$

Rješenje je $x_1 = 4$ i $x_2 =-2$, što je točno do na točnost stroja.

Međutim, mala promjena na desnoj strani daje

$\left(\begin{array}{cc|c}
(1+2\epsilon)&     (1+2\epsilon )&     2+4\epsilon \cr
(1+\epsilon)&     1   &2 +\epsilon \end{array}\right).$

Točno rješenje je $x_1=x_2=1$, ali zbog zaokruživanja dobije $x_1 =0$, $x_2 =2$. __Objasnite.__
Niti jedan trik kojeg smo vidjeli ne daje točno rješenje. 
"""

# ╔═╡ c471e93a-5c8a-433f-880a-9e26788fa601
[1+2*eps() 1+2*eps(); 1+eps() 1]\[2+4*eps(); 2+eps()]

# ╔═╡ 9983cc49-9398-4772-9d04-3f7bbf2b47a1
[BigFloat(1)+2*eps() 1+2*eps(); 1+eps() 1]\[BigFloat(2)+4*eps(); 2+eps()]

# ╔═╡ de3e2152-7ec2-4366-a6fa-dc1476d19480
md"""
__Razlog.__  IEEE aritmetika ovaj sustav zaokruži na sustav


$$\left(\begin{array}{cc|c}
(1+2\epsilon)&     (1+2\epsilon )&     2+4\epsilon \cr
(1+\epsilon)&     1   &2           \end{array}\right) $$

čija su rješenja $x_1=0$ i $x_2=2$. Ovaj problem je vrlo blizu singularnom sustavu

$$\left(\begin{array}{cc|c} 1 & 1 & 2 \cr 1 & 1 & 2\end{array}\right)$$

koji ima parametarska rješenja

$$\mathbf{x} =\begin{pmatrix} x_1 \\ x_2\end{pmatrix} = \begin{pmatrix} 1\\1\end{pmatrix}+ 
\beta \begin{pmatrix}-1 \\1\end{pmatrix}, \quad \beta \in \mathbb{R}.$$

Primijetimo da su $\begin{pmatrix} x_1 \\ x_2\end{pmatrix}= \begin{pmatrix} 1\\1\end{pmatrix}$ i $\begin{pmatrix} x_1 \\ x_2\end{pmatrix}=\begin{pmatrix} 0\\2\end{pmatrix}$ dva od tih rješenja.

__Pitanje.__ Koja je geometrijska interpretacija ovog sustava?
"""

# ╔═╡ c9bc269a-a306-4a35-acd8-aad3de58f56a
md"""
## LU rastav

$A=\begin{pmatrix}\alpha & a^T \cr b  & B \end{pmatrix}=
\begin{pmatrix} 1 & 0 \cr l & I \end{pmatrix}
\begin{pmatrix} \beta & u^T \cr 0 & C \end{pmatrix}
=LU=\begin{pmatrix} \beta & u^T \cr l\beta & lu^T+ C\end{pmatrix}$

povlači

$\beta=\alpha,\quad u=a,\quad l=b\beta^{-1},\quad C=B-lu^T=B-b\beta^{-1}a^T.$

Indukcija daje sljedeći algoritam:

"""

# ╔═╡ 66fcc372-4f27-4092-9552-8eb6e863bd4a
function mylu(A₁::Array{T}) where T # Strang, str. 100
    A=copy(A₁)
    n,m=size(A)
    # Ovo prihvaća brojeve i blok-matrice
    U=map(Float64,[zero(A[1,1]) for i=1:n, j=1:n])
    L=map(Float64,[zero(A[1,1]) for i=1:n, j=1:n])
    for k=1:n
        L[k,k]=one(A[1,1])
        for i=k+1:n
            L[i,k]=A[i,k]/A[k,k]
            for j=k+1:n
                A[i,j]=A[i,j]-L[i,k]*A[k,j]
            end
        end
        for j=k:n
            U[k,j]=A[k,j]
        end
    end
    return L,U
end

# ╔═╡ 696152e0-0d76-11eb-02d2-8b2209417bc2
import Random

# ╔═╡ b62492f0-0d7f-11eb-0e42-d377410cec70
Random.seed!(123);

# ╔═╡ 19567c60-0d82-11eb-3405-8d0312a34b5f
begin
	n=6
	A=rand(n,n)
	b=rand(n)
end

# ╔═╡ 1e51f7d0-0d82-11eb-238b-f1179d7f9a30
A

# ╔═╡ 251a1cf0-0d82-11eb-3747-cf84a824570f
L,U=mylu(A)

# ╔═╡ 2a725e60-0d82-11eb-18ce-cd017246fc46
L

# ╔═╡ 2d7fc570-0d82-11eb-06d0-b115bbbd897c
U

# ╔═╡ ddbe2df0-0d82-11eb-0a77-e9de7611ec9b
L*U-A

# ╔═╡ 946eb7cd-8a97-4aa8-880c-bfba8e6efae1
md"""
## Trokutasti sustavi
"""

# ╔═╡ 29d27f7f-d2e5-4d1f-a667-39fc266ffa17
begin
	function myU(U::Array{T},b₁::Array{T}) where T
	    b=copy(b₁)
	    n=length(b)
	    for i=n:-1:1
	       for j=n:-1:i+1
	            b[i]=b[i]-U[i,j]*b[j]
	       end
	        b[i]=b[i]/U[i,i]
	    end
	    b
	end
	
	function myL(L::Array{T},b₁::Array{T}) where T
	    b=copy(b₁)
	    n=length(b)
	    for i=1:n
	        for j=1:i-1
	            b[i]=b[i]-L[i,j]*b[j]
	        end
	        b[i]=b[i]/L[i,i]
	    end
	    b
	end
end

# ╔═╡ 6a7341a2-0d82-11eb-0831-05182f30ffe3
# Riješimo sustav koristeći ugrađenu funkciju
x=A\b

# ╔═╡ f3732010-0d82-11eb-3484-c38c2d6a1f31
# Riješimo sustav koristeći naše funkcije
y=myL(L,b)

# ╔═╡ f6da8e00-0d82-11eb-2179-b3cfdddf689a
x₁=myU(U,y)

# ╔═╡ 04064a10-0d83-11eb-2efe-7d7ca4ef35e7
# Usporedimo rješenja
x-x₁

# ╔═╡ 8a3a81c0-0df9-11eb-087a-c9298b2fd265
@which lu(A)

# ╔═╡ 07fae3c0-0dfa-11eb-02a3-b7efa4490b1c
?lu

# ╔═╡ 4225c750-b668-4331-b6b8-0509635e69c6
md"""
## Brzina

Funkcija `mylu()` je spora. Između ostalog, alocira nepotrebno tri matrice i ne računa s blok matricama.

Funkcija se može preformulirati tako da su i $L$ i $U$ spremljene u polje $A$, pri čemu se dijagonala od $L$ ne sprema jer su svi elementi jednaki 1 (vidi
[Gilbert Strang, 'Introduction to Linear Algebra, 4th Edition', Wellesley-Cambridge Press, 2009, str.100](https://books.google.hr/books?id=M19gPgAACAAJ&dq=strang%20introduction&hl=hr&source=gbs_book_other_versions))

"""

# ╔═╡ 17427400-0d83-11eb-14e2-f5d29e1650e4
function mylu₁(A₁::Array{T}) where T # Strang, str. 100
    A=copy(A₁)
    n,m=size(A)
    for k=1:n-1
        ρ=k+1:n
        A[ρ,k]=A[ρ,k]/A[k,k]
        A[ρ,ρ]=A[ρ,ρ]-A[ρ,k]*A[k,ρ]'
    end
    A
end

# ╔═╡ 0b44bb30-0d84-11eb-1d3a-dfc67cf3cf20
mylu₁(A)

# ╔═╡ 13a09fb2-0d84-11eb-15d6-e1d3b4ba658e
L

# ╔═╡ 1ccfd9c0-0d84-11eb-097b-2bfde3a9790f
U

# ╔═╡ 6f3f3257-8dd9-4d3c-b18e-cdbb37d52e2a
md"""
Usporedimo brzine ugrađene funkcije `lu()` koja koristi LAPACK-ove rutine i naše naivne funkcije `mylu()`na većoj dimenziji. 

Izvedite par puta radi točnijeg mjerenja brzine.
"""

# ╔═╡ 5da85850-0d84-11eb-091b-df4a89e4d052
A₁=rand(512,512);

# ╔═╡ 46868fc0-0d84-11eb-0bea-f9ee72af7795
lu(A₁)

# ╔═╡ 9ce0bb70-0d84-11eb-14ac-5335e9985dbf
mylu₁(A₁);

# ╔═╡ 09b70bcd-43ad-46b2-9664-2809351f9f70
md"""
### Blok varijanta

`mylu()` i $\mathsf{mylu}_1()$ su nekoliko desetaka puta sporiji od `lu()`.

Preradimo $\mathsf{mylu}_1()$ za rad s blokovima (još uvijek nemamo ugrađeno pivotiranje!):
"""

# ╔═╡ e80e7d80-0d84-11eb-2423-23867085be67
function mylu₂(A₁::Array{T}) where T # Strang, page 100
    A=copy(A₁)
    n,m=size(A)
    for k=1:n-1
        for ρ=k+1:n
            A[ρ,k]=A[ρ,k]/A[k,k]
            for l=k+1:n
                A[ρ,l]=A[ρ,l]-A[ρ,k]*A[k,l]
            end
        end
    end
    A
end

# ╔═╡ 9d740956-6b11-4c0f-bca9-58fcfa852a62
md"""
Napravimo prvo mali test, $k=2$, $l=4$:
"""

# ╔═╡ 090a3f10-0d85-11eb-0181-fdc5aa091df7
begin
	k,l=2,4   # 32,16
	Ab=[rand(k,k) for i=1:l, j=1:l]
end

# ╔═╡ 21e796b0-0dfb-11eb-2493-6fe0d3c70180
Ab[2,1]

# ╔═╡ 24c23190-0d85-11eb-382d-2536a040dfc8
A₀=mylu₂(Ab);

# ╔═╡ 3d486c20-0d85-11eb-2c81-a7d99c190907
begin
	# Provjera
	U₀=triu(A₀)
	L₀=tril(A₀)
	for i=1:maximum(size(L₀))
		L₀[i,i]=Matrix{Float64}(I,size(L₀[1,1])) # eye(L[1,1])
	end
end

# ╔═╡ 642b68b0-0d85-11eb-01d7-734195427bd9
L₀

# ╔═╡ 46ecd28e-0dfb-11eb-1e7a-2d55904f621e
L₀[1,1]

# ╔═╡ 9615e1c0-0d85-11eb-161d-39a7eff4046f
Rezidual=L₀*U₀-Ab

# ╔═╡ aa933e40-0d85-11eb-2141-d36ff3fc471b
# Pretvaranje blok matrice u običnu
unblock(A) = mapreduce(identity, hcat, [mapreduce(identity, vcat, A[:,i]) for i = 1:size(A,2)])

# ╔═╡ b307b3d0-0d85-11eb-25dc-83e2d3fbcb4f
norm(unblock(Rezidual))

# ╔═╡ 86f3ce48-73ce-4626-914f-478cf3ad1154
md"""
Probajmo veće dimenzije ($n=k\cdot l$).
"""

# ╔═╡ 83801210-b142-4d11-8eac-5eab72a181b3
md"""
Vidimo da je $\mathsf{mylu}_2()$ gotovo jednako brz kao `lu()` (na jednoj jezgri), uz napomenu da $\mathsf{mylu}_2()$ nema ugrađeno pivotiranje. 
Program još uvijek nije optimalan jer alocira previše memorije.
"""

# ╔═╡ 9a2ef752-2231-4282-aa5c-275894c21de5
md"""
## Pivotiranje

Standardne implementacije uvijek računaju Gaussovu eliminaciju koristeći __parcijalno pivotiranje__ :

u svakom koraku se retci pivotiraju (zamijene) tako da pivotni element ima najveću apsolutnu vrijednost u danom stupcu. Na taj 
način je 

$$|L_{ij}| \leq 1,$$

što u praksi dovoljno spriječava rast elemenata.
"""

# ╔═╡ 535aa570-77c2-4962-97c2-9661884a21c2
A₂=[0.00003 1;2 3]

# ╔═╡ 4e739aee-0d86-11eb-056f-8589740ddc96
L₂,U₂=mylu(A₂)

# ╔═╡ 64384480-0d86-11eb-3ac7-7f602aeaa6d8
begin
	# s pivoritranjem
	P=[0 1;1 0]
	mylu(P*A₂)
end

# ╔═╡ 416e1cd9-0fe8-4288-a1cb-b60f60139fa5
md"""
Sada probajmo ugrađenu funkciju - koristimo je precizno.
"""

# ╔═╡ bc8ce0a0-0d86-11eb-1342-413ee4eb89ef
?lu

# ╔═╡ e86a4730-0d86-11eb-266b-5b41924f61a8
F=lu(A₂)

# ╔═╡ fb2a1530-0d86-11eb-0d42-77246926cb7f
F.L*F.U == A₂[F.p, :]

# ╔═╡ 5b7d3a80-0dfe-11eb-05b2-2b3501e351fb
F.L*F.U==F.P*A₂

# ╔═╡ 68c72cf0-0dfe-11eb-37fb-d5ac2e1142e7
F.P

# ╔═╡ 4483e520-0d88-11eb-1eb4-25f4eaf1a88d
# Probajmo s prethodnom matricom A
A

# ╔═╡ 573618c0-0def-11eb-0101-e1b09c384512
F₄=lu(A)

# ╔═╡ 74b55dc0-0def-11eb-37e2-71ba59aa8295
F₄.L*F₄.U==A[F₄.p,:]

# ╔═╡ 90e4a320-0def-11eb-189e-7fc9c7b53942
# Javljaju se greške zaokruživanja
F₄.L*F₄.U-A[F₄.p,:]

# ╔═╡ 630a82aa-6998-4325-a0c9-d44f60c0df31
md"""
### Potpuno pivotiranje

Sljedeći program računa Gaussovu eliminaciju koristeći __potpuno pivotiranje__ - u svakom koraku 
se retci i stupci zamijene takoda se na pivotnu poziciju dovede element koji ima najveću 
apsolutnu vrijednost u trenutnoj podmatrici.

Rast elemenata je teoretski ograničen, ali ograda je $O(2^n)$, što nije korisno.
"""

# ╔═╡ 18ad03b0-0d87-11eb-06f9-45ac1b7e3b04
function gecp(A₁::Array{T}) where T
    # Gaussova eliminacija s potpunim pivotiranjem
    # Izlaz: Pr*L*U*Pc'=A ili Pr'*A*Pc=L*U
    A=copy(A₁)
    n,m=size(A)
    Pr=Matrix{Float64}(I,n,n)
    Pc=Matrix{Float64}(I,n,n)
    D=zeros(n)
    for i=1:n-1
        amax,indm=findmax(abs.(A[i:n,i:n]))
        imax=indm[1]+i-1
        jmax=indm[2]+i-1
        #  zamijena redaka
        if (imax != i)
            temp = Pr[:,i]
            Pr[:,i] = Pr[:,imax]
            Pr[:,imax] = temp
            temp = A[i,:]
            A[i,:] = A[imax,:]
            A[imax,:] = temp
        end
        # zamijena stupaca
        if (jmax != i)
            temp = Pc[:,i]
            Pc[:,i] = Pc[:,jmax]
            Pc[:,jmax] = temp
            temp = A[:,i]
            A[:,i] = A[:,jmax]
            A[:,jmax] = temp
        end
        # eliminacija
        D[i]=A[i,i]
        A[i+1:n,i] = A[i+1:n,i]/D[i]
        A[i+1:n,i+1:n] = A[i+1:n,i+1:n] - A[i+1:n,i]*A[i,i+1:n]'
        A[i,i+1:n]=A[i,i+1:n]/D[i]
    end
    D[n]=A[n,n]
    L=I+tril(A,-1)
    U=I+triu(A,1)
    U=Diagonal(D)*U
    return L,U,Pr,Pc
end

# ╔═╡ 1f3a42b0-0d87-11eb-1fef-8f0a35eb3cce
Lₚ,Uₚ,Pr,Pc=gecp(A)

# ╔═╡ 3d328840-0d87-11eb-3c28-5dbb05be31f8
norm(Pr*Lₚ*Uₚ*Pc'-A)

# ╔═╡ 1ea5da30-0df0-11eb-2bb6-3bcf58c68adb
yₚ=myL(Lₚ,Pr'*b)

# ╔═╡ 284e14d0-0df0-11eb-2255-7b26982e1bbf
zₚ=myU(Uₚ,yₚ)

# ╔═╡ 355977a0-0df0-11eb-0e2b-0b5161d7979e
xₚ=Pc*zₚ

# ╔═╡ 3f504770-0df0-11eb-3049-ddac0626728f
# Rezidual 
A*xₚ-b

# ╔═╡ 55d12cd0-0d87-11eb-10cc-edca8db298a1
md"""
## Točnost

Neka je zadan sustav $Ax=b$, pri čemu je matrica $A$ regularna.

Da bi primijenili koncepte iz bilježnice 
[NA04 Pogreska unatrag i stabilni algoritmi](NA04%20Pogreska%20unatrag%20i%20stabilni%20algoritmi.ipynb), potrebno je:

1. napraviti teoriju smetnje za dani problem
2. analizirati pogreške algoritma (Gaussove eliminacije)

### Teorija smetnje

Neka je 

$(A+\delta A)\hat x=(b+\delta b)$

za neki $\hat x=x+\delta x$.

Želimo ocijeniti 

$$\frac{\| \hat x - x \|}{\| x\|} \equiv \frac{\| \delta x\|}{\| x\|}.$$

Uvedimo oznake (npr. prema [Matrix Computations, poglavlje 2.6.2](https://books.google.hr/books?id=X5YfsuCWpxMC&printsec=frontcover&hl=hr#v=onepage&q&f=false))

$$\delta A=\varepsilon F, \quad \delta b=\varepsilon f, \qquad \hat x=x(\varepsilon)$$

čime smo dobili jednodimenzionalni problem 

$$(A+\varepsilon F)\,x(\varepsilon)=b+\varepsilon f$$

za neke (nepoznate) matricu $F$ i vektor $f$. 

Deriviranje po $\varepsilon$ daje

$$Fx(\varepsilon)+(A+\varepsilon F)\, \dot x(\varepsilon)=f.$$

Uvrštavanje $\varepsilon=0$ daje

$$F x+A\dot x(0)=f,$$

odnosno

$$\dot x(0)=A^{-1}(f-Fx).$$

Taylorov razvoj oko $\varepsilon=0$ glasi

$$x(\varepsilon)=x(0)+\varepsilon \dot x(0) +O(\varepsilon^2),$$

odnosno, uz zanemarivanje člana $O(\varepsilon^2)$,

$$\hat x-x=\varepsilon A^{-1}(f-Fx)=A^{-1} (\varepsilon f - \varepsilon F x) = A^{-1} (\delta b - \delta A x).$$

Svojstva norme povlače

$$\| \hat x-x\|\leq \| A^{-1} \| (\| \delta b \|  + \| \delta A \| \cdot \|  x\| ).$$

Konačno, zbog $\| b\| \leq \| A\| \| x\|$, imamo

$$\frac{\| \hat x-x\|}{\| x\|}\leq \| A\|  \cdot \| A^{-1} \| \bigg(\frac{\| \delta b \|}{\|b\|}  + \frac{\| \delta A \|}{ \|  A\|} \bigg). \tag{1}$$

Broj 

$$\kappa(A)\equiv \| A\|  \cdot \| A^{-1} \|$$ 

je __uvjetovanost__ (__kondicija__)  matrice $A$ i kazuje nam 
koliko se relativno uvećaju relativne promjene u polaznim podacima (matrici $A$ i vektoru $b$).

Pogledajmo primjer iz [Numeričke matematike, str. 42](http://www.mathos.unios.hr/pim/Materijali/Num.pdf):

"""

# ╔═╡ 06273c00-0d88-11eb-2259-230e34f04417
A₃= [0.234 0.458; 0.383 0.750]

# ╔═╡ 542b76a0-0d88-11eb-0672-c95813a3ccdc
b₃=[0.224;0.367]

# ╔═╡ 29a64d0e-0d88-11eb-0f2b-dfd116b214c4
mylu₁(A₃)

# ╔═╡ 3b73569e-0d88-11eb-271c-b983eb9cb3f5
F₃=lu(A₃)

# ╔═╡ f2d157d0-0df0-11eb-2ab9-4d55d1b5e307
x₃=A₃\b₃

# ╔═╡ fd0b7230-0df0-11eb-1443-67ea60bb2b7f
x₃[1]

# ╔═╡ 1e4c2c00-0df1-11eb-2ca2-3b4f08d92e9a
x₃[2]

# ╔═╡ 2bff8ea0-0df1-11eb-08a3-375007f3f276
begin
	δb₃=[0.00009; 0.000005]
	A₃\(b₃+δb₃)
end

# ╔═╡ 5179d372-0df1-11eb-0183-8bcc73149584
begin
	δA₃=[-0.001 0;0 0]
    x₄=(A₃+δA₃)\b₃
end

# ╔═╡ 69b8cbd0-0df1-11eb-3fcd-a9f0865efdce
cond(A₃), norm(δA₃)/norm(A₃), norm(x₄-x₃)/norm(x₃)

# ╔═╡ a47802e0-0df1-11eb-3f9f-2fe1ebc781fd
md"""
### Pogreška Gaussove eliminacije

Prema [Matrix Computations, poglavlje 3.3](https://books.google.hr/books?id=X5YfsuCWpxMC&printsec=frontcover&hl=hr#v=onepage&q&f=false), za izračunate faktore
$\hat L$ i $\hat U$ vrijedi

$$\hat L\cdot \hat U = A+\delta A$$

gdje je (nejednakost se čita po elementima matrica, $\varepsilon$ je sada točnost stroja)

$$| \delta A|\leq 3(n-1) \varepsilon (|A|+|\hat L| \cdot |\hat U|) +O(\varepsilon^2).$$

Zanemarivanje člana $O(\varepsilon^2)$ i prelazak na normu daju

$$\|\delta A \| \lesssim O(n)\varepsilon (\| A\| + \| \hat L\| \cdot \| \hat U\|),$$

pa je 

$$\frac{\|\delta A \|}{\|A\|} \lesssim O(n)\varepsilon \bigg(1+\frac{\| \hat L\| \cdot \| \hat U\|}{\|A\|}\bigg).$$

Ukoliko se Gaussova eliminacija radi s pivotiranjem, tada će najvjerojatnije zadnji kvocijent biti malen 
($\approx 1$). Također, pogreška kod rješavanja trokutastih sustava nije veća od navedene pa uvrštavanjem u (1) slijedi 
da za relativnu pogrešku izračunatog rješenja vrijedi 

$$\frac{\| \hat x-x\|}{\| x\|}\leq \kappa(A) O(n\varepsilon).$$

Zaključimo:

> _Ukoliko je kondicija matrice velika, rješenje može biti netočno._

__Primjer.__ Vandermondeova matrica

"""

# ╔═╡ f1c74560-0df1-11eb-19a7-c9ad6aed7410
begin
	nᵥ=10
	v=rand(nᵥ)
end

# ╔═╡ 026d0d00-0df2-11eb-26fd-cbc13048e56c
begin
	# Vandermondeove matrice imaju veliku kondiciju.
	V=Array{Float64}(undef,nᵥ,nᵥ)
	for i=1:nᵥ
	    V[:,i]=v.^(i-1)
	end
	V=V'
end

# ╔═╡ 2f58af40-0df2-11eb-28e4-5b1911f53b83
bᵥ=rand(nᵥ)

# ╔═╡ 3da4f680-0df2-11eb-23d4-f3fb28cdc8e7
xᵥ=V\bᵥ

# ╔═╡ 439224f0-0df2-11eb-2e57-539a3470de32
cond(V)

# ╔═╡ 50da42a0-0df2-11eb-2ded-c52a76acc155
begin
	Vbig=map(BigFloat,V)
	bbig=map(BigFloat,bᵥ)
	xbig=Vbig\bbig;
end

# ╔═╡ 56e4bd10-0df2-11eb-0152-556ef692a70e
map(Float64,norm(xbig-xᵥ)/norm(xbig))

# ╔═╡ 50ed08c1-391f-4baa-8cfa-54db04038fb1
md"""
### Umjetno loša kondicija
"""

# ╔═╡ a496c05e-0def-11eb-0ae1-83f3cdccf36e
Aᵤ=[1 1; 1 2]

# ╔═╡ aeee4670-0df2-11eb-0ef9-0bb353d8ebfe
bᵤ=[1;3]

# ╔═╡ aef0b770-0df2-11eb-3f66-8d09a5970a49
xᵤ=Aᵤ\bᵤ

# ╔═╡ aef1a1d0-0df2-11eb-0aeb-c5f670c48b32
xᵤ,cond(Aᵤ)

# ╔═╡ af144500-0df2-11eb-37e0-af9d5cd06c65
A₅=[1e-4 1e-4;1 2]

# ╔═╡ af16b600-0df2-11eb-2852-7f1e61314674
b₅=[1e-4;3]

# ╔═╡ af2b7680-0df2-11eb-2ff2-6b188cbd6b7f
x₅=A₅\b₅

# ╔═╡ af2dc072-0df2-11eb-2958-19fc8d01e81d
x₅,cond(A₅),xᵤ-x₅

# ╔═╡ bfc0ea15-556f-4109-b626-cb724ee14bfd
md"""
### Procjena kondicije

Računanje kondicije prema definiciji $\kappa(A)=\|A\| \cdot \|A^{-1}\|$ zahtijeva računanje inverzne matrice, za što je potrebno $O(n^3)$ operacija. To je isti red veličine operacija koji je potreban za rješavanje zadanog sustava. Prilikom rješavanja sustava na raspolaganju su nam trokutasti faktori $L$ i $U$, što se može iskoristiti kako bi se kondicija približno izračunala u $O(n^2)$ operacija. 
Detalji se nalaze u [Matrix Computations, poglavlje 3.5.4](https://books.google.hr/books/about/Matrix_Computations.html?id=X5YfsuCWpxMC&redir_esc=y). 
LAPACK rutina 
[dtrcon.f](http://www.netlib.org/lapack/explore-html/d9/d84/dtrcon_8f_source.html) računa inverz približne kondicije trokutaste matrice.

Izračunajmo približnu kondiciju Vandermondeove matrice iz prethodnog primjera.
 
"""

# ╔═╡ e0964d99-9170-427e-92b0-061a83362e85
md"""
Izračunajmo rezidual za Vandermondeov sustav:
"""

# ╔═╡ 79d7e431-4090-4834-8ae3-da04f5517285
md"""
Zaključujemo da je rješenje $x_v$ izračunato stabilno, odnosno s vrlo malom pogreškom unatrag u početnim podatcima. To još uvijek ne znači da je rješenje relativno vrlo točno.
"""

# ╔═╡ 8707357c-84c9-4d30-aa81-1172d7ac715e
?LAPACK.trcon!

# ╔═╡ e38b8370-0df2-11eb-0ed5-ab750f73de17
begin
	Fᵥ=lu(V)
	cond(V,1), cond(Fᵥ.L,1), cond(Fᵥ.U,1)
end

# ╔═╡ 11add0f0-0df3-11eb-2f01-5b985574b265
1 ./LAPACK.trcon!('O','L','U',Fᵥ.L),1 ./LAPACK.trcon!('O','U','N',Fᵥ.U)

# ╔═╡ 24dc3f3e-0df3-11eb-04f6-a58c63e5ba58
md"""
## Rezidual


Izračunato rješenje $\hat x$ sustava $Ax=b$ je točno rješenje nekog sličnog sustava (vidi [Afternotes on Numerical Analysis, str. 128](https://books.google.hr/books?id=w-2PWh01kWcC&printsec=frontcover&hl=hr#v=onepage&q&f=false)):


$$(A+\delta A)\,\hat x=b. \tag{1}$$

__Rezidual__ (ili __ostatak__) definiramo kao 

$$r=b-A\hat x.$$

Tada je 

$$0=b-(A+\delta A)\,\hat x=r- \delta A\,\hat x$$

pa je 

$$\| r\| = \| \delta A\,\hat x \| \leq \| \delta A\| \cdot \|\hat x \|,$$

odnosno

$$\frac{\|  \delta A\|}{\|A \|} \geq \frac{\|r\|}{\| A\| \cdot \|\hat x \|}.$$

Dakle, ako  __relativni rezidual__ 

$$\frac{r}{\| A\| \cdot \|\hat x \|}$$

ima veliku normu, tada _rješenje nije izračunato stabilno._

S druge strane, ako relativni rezidual ima malu normu, tada je _rješenje izračunato stabilno_. Naime, za

$$\delta A=\frac{r\hat x^T}{\|\hat x\|^2}$$

vrijedi (1):

$$b-(A+\delta A)\hat x=(b-A\hat x)-\delta A \hat x = r-\frac{r\hat x^T \hat x}{\|\hat x\|^2}
= r-\frac{r \|\hat x^T \hat x\|}{\|\hat x\|^2}=r-r=0.$$

Također vrijedi

$$\frac{\|  \delta A\|}{\|A \|}  \leq  \frac{\|r\|\|\hat x \|}{\| A\| \cdot \|\hat x \|^2}=
\frac{\|r\|}{\| A\| \cdot \|\hat x \|}.$$

Izračunajmo reziduale za prethodni primjer dimenzije $2$:

"""

# ╔═╡ 88f2801e-0df3-11eb-35ac-c32cb53aef8a
rᵤ=bᵤ-Aᵤ*xᵤ

# ╔═╡ e5108050-0df3-11eb-2d96-fddf9e91ef9e
norm(rᵤ)/(norm(Aᵤ)*norm(xᵤ))

# ╔═╡ f77005e0-0df3-11eb-0d2b-3b00dddef4f3
r₅=b₅-A₅*x₅

# ╔═╡ 0776ce62-0df4-11eb-1f95-3900b12d5087
norm(r₅)/(norm(A₅)*norm(x₅))

# ╔═╡ 1e40da00-0df4-11eb-3d74-03fb919b4781
md"
Izračunajmo rezidual za Vandermondeov sustav:
"

# ╔═╡ 2ec7cf00-0df4-11eb-03c6-03c37219650d
rᵥ=bᵥ-V*xᵥ

# ╔═╡ 3a45b400-0df4-11eb-3e3c-41fed6f7a499
norm(rᵥ)/(norm(V)*norm(xᵥ))

# ╔═╡ 40c1dc00-0df4-11eb-10a6-bf598057f7fb
md"
Zaključujemo da je rješenje $x_v$ izračunato stabilno, odnosno s vrlo malom pogreškom unatrag u početnim podatcima. To još uvijek ne znači da je rješenje relativno vrlo točno.
"

# ╔═╡ Cell order:
# ╟─63842a33-b596-4f3f-bde5-7953ec4dc840
# ╟─48167400-0d83-11eb-1c7d-359a2574c8b1
# ╟─221d2474-de59-4042-918f-534305d8708f
# ╠═89291e44-4ea6-4e74-99bd-a30e8ee4d895
# ╟─ed629240-6b9a-4a85-b412-97a06565a9cf
# ╠═b8b17cd3-c828-438d-a2c0-13b1965ed778
# ╟─eb5aab00-0d73-11eb-2f45-771b2e23a5e3
# ╠═c471e93a-5c8a-433f-880a-9e26788fa601
# ╠═9983cc49-9398-4772-9d04-3f7bbf2b47a1
# ╟─de3e2152-7ec2-4366-a6fa-dc1476d19480
# ╟─c9bc269a-a306-4a35-acd8-aad3de58f56a
# ╠═66fcc372-4f27-4092-9552-8eb6e863bd4a
# ╠═53d8d8e0-0d7f-11eb-358f-3f001e5e7d78
# ╠═696152e0-0d76-11eb-02d2-8b2209417bc2
# ╠═b62492f0-0d7f-11eb-0e42-d377410cec70
# ╠═19567c60-0d82-11eb-3405-8d0312a34b5f
# ╠═1e51f7d0-0d82-11eb-238b-f1179d7f9a30
# ╠═251a1cf0-0d82-11eb-3747-cf84a824570f
# ╠═2a725e60-0d82-11eb-18ce-cd017246fc46
# ╠═2d7fc570-0d82-11eb-06d0-b115bbbd897c
# ╠═ddbe2df0-0d82-11eb-0a77-e9de7611ec9b
# ╟─946eb7cd-8a97-4aa8-880c-bfba8e6efae1
# ╠═29d27f7f-d2e5-4d1f-a667-39fc266ffa17
# ╠═6a7341a2-0d82-11eb-0831-05182f30ffe3
# ╠═f3732010-0d82-11eb-3484-c38c2d6a1f31
# ╠═f6da8e00-0d82-11eb-2179-b3cfdddf689a
# ╠═04064a10-0d83-11eb-2efe-7d7ca4ef35e7
# ╠═8a3a81c0-0df9-11eb-087a-c9298b2fd265
# ╠═07fae3c0-0dfa-11eb-02a3-b7efa4490b1c
# ╟─4225c750-b668-4331-b6b8-0509635e69c6
# ╠═17427400-0d83-11eb-14e2-f5d29e1650e4
# ╠═0b44bb30-0d84-11eb-1d3a-dfc67cf3cf20
# ╠═13a09fb2-0d84-11eb-15d6-e1d3b4ba658e
# ╠═1ccfd9c0-0d84-11eb-097b-2bfde3a9790f
# ╟─6f3f3257-8dd9-4d3c-b18e-cdbb37d52e2a
# ╠═5da85850-0d84-11eb-091b-df4a89e4d052
# ╠═46868fc0-0d84-11eb-0bea-f9ee72af7795
# ╠═9ce0bb70-0d84-11eb-14ac-5335e9985dbf
# ╟─09b70bcd-43ad-46b2-9664-2809351f9f70
# ╠═e80e7d80-0d84-11eb-2423-23867085be67
# ╟─9d740956-6b11-4c0f-bca9-58fcfa852a62
# ╠═090a3f10-0d85-11eb-0181-fdc5aa091df7
# ╠═21e796b0-0dfb-11eb-2493-6fe0d3c70180
# ╠═24c23190-0d85-11eb-382d-2536a040dfc8
# ╠═3d486c20-0d85-11eb-2c81-a7d99c190907
# ╠═642b68b0-0d85-11eb-01d7-734195427bd9
# ╠═46ecd28e-0dfb-11eb-1e7a-2d55904f621e
# ╠═9615e1c0-0d85-11eb-161d-39a7eff4046f
# ╠═aa933e40-0d85-11eb-2141-d36ff3fc471b
# ╠═b307b3d0-0d85-11eb-25dc-83e2d3fbcb4f
# ╟─86f3ce48-73ce-4626-914f-478cf3ad1154
# ╟─83801210-b142-4d11-8eac-5eab72a181b3
# ╟─9a2ef752-2231-4282-aa5c-275894c21de5
# ╠═535aa570-77c2-4962-97c2-9661884a21c2
# ╠═4e739aee-0d86-11eb-056f-8589740ddc96
# ╠═64384480-0d86-11eb-3ac7-7f602aeaa6d8
# ╟─416e1cd9-0fe8-4288-a1cb-b60f60139fa5
# ╠═bc8ce0a0-0d86-11eb-1342-413ee4eb89ef
# ╠═e86a4730-0d86-11eb-266b-5b41924f61a8
# ╠═fb2a1530-0d86-11eb-0d42-77246926cb7f
# ╠═5b7d3a80-0dfe-11eb-05b2-2b3501e351fb
# ╠═68c72cf0-0dfe-11eb-37fb-d5ac2e1142e7
# ╠═4483e520-0d88-11eb-1eb4-25f4eaf1a88d
# ╠═573618c0-0def-11eb-0101-e1b09c384512
# ╠═74b55dc0-0def-11eb-37e2-71ba59aa8295
# ╠═90e4a320-0def-11eb-189e-7fc9c7b53942
# ╟─630a82aa-6998-4325-a0c9-d44f60c0df31
# ╠═18ad03b0-0d87-11eb-06f9-45ac1b7e3b04
# ╠═1f3a42b0-0d87-11eb-1fef-8f0a35eb3cce
# ╠═3d328840-0d87-11eb-3c28-5dbb05be31f8
# ╠═1ea5da30-0df0-11eb-2bb6-3bcf58c68adb
# ╠═284e14d0-0df0-11eb-2255-7b26982e1bbf
# ╠═355977a0-0df0-11eb-0e2b-0b5161d7979e
# ╠═3f504770-0df0-11eb-3049-ddac0626728f
# ╟─55d12cd0-0d87-11eb-10cc-edca8db298a1
# ╠═06273c00-0d88-11eb-2259-230e34f04417
# ╠═542b76a0-0d88-11eb-0672-c95813a3ccdc
# ╠═29a64d0e-0d88-11eb-0f2b-dfd116b214c4
# ╠═3b73569e-0d88-11eb-271c-b983eb9cb3f5
# ╠═f2d157d0-0df0-11eb-2ab9-4d55d1b5e307
# ╠═fd0b7230-0df0-11eb-1443-67ea60bb2b7f
# ╠═1e4c2c00-0df1-11eb-2ca2-3b4f08d92e9a
# ╠═2bff8ea0-0df1-11eb-08a3-375007f3f276
# ╠═5179d372-0df1-11eb-0183-8bcc73149584
# ╠═69b8cbd0-0df1-11eb-3fcd-a9f0865efdce
# ╟─a47802e0-0df1-11eb-3f9f-2fe1ebc781fd
# ╠═f1c74560-0df1-11eb-19a7-c9ad6aed7410
# ╠═026d0d00-0df2-11eb-26fd-cbc13048e56c
# ╠═2f58af40-0df2-11eb-28e4-5b1911f53b83
# ╠═3da4f680-0df2-11eb-23d4-f3fb28cdc8e7
# ╠═439224f0-0df2-11eb-2e57-539a3470de32
# ╠═50da42a0-0df2-11eb-2ded-c52a76acc155
# ╠═56e4bd10-0df2-11eb-0152-556ef692a70e
# ╟─50ed08c1-391f-4baa-8cfa-54db04038fb1
# ╠═a496c05e-0def-11eb-0ae1-83f3cdccf36e
# ╠═aeee4670-0df2-11eb-0ef9-0bb353d8ebfe
# ╠═aef0b770-0df2-11eb-3f66-8d09a5970a49
# ╠═aef1a1d0-0df2-11eb-0aeb-c5f670c48b32
# ╠═af144500-0df2-11eb-37e0-af9d5cd06c65
# ╠═af16b600-0df2-11eb-2852-7f1e61314674
# ╠═af2b7680-0df2-11eb-2ff2-6b188cbd6b7f
# ╠═af2dc072-0df2-11eb-2958-19fc8d01e81d
# ╟─bfc0ea15-556f-4109-b626-cb724ee14bfd
# ╟─e0964d99-9170-427e-92b0-061a83362e85
# ╟─79d7e431-4090-4834-8ae3-da04f5517285
# ╠═8707357c-84c9-4d30-aa81-1172d7ac715e
# ╠═e38b8370-0df2-11eb-0ed5-ab750f73de17
# ╠═11add0f0-0df3-11eb-2f01-5b985574b265
# ╟─24dc3f3e-0df3-11eb-04f6-a58c63e5ba58
# ╠═88f2801e-0df3-11eb-35ac-c32cb53aef8a
# ╠═e5108050-0df3-11eb-2d96-fddf9e91ef9e
# ╠═f77005e0-0df3-11eb-0d2b-3b00dddef4f3
# ╠═0776ce62-0df4-11eb-1f95-3900b12d5087
# ╟─1e40da00-0df4-11eb-3d74-03fb919b4781
# ╠═2ec7cf00-0df4-11eb-03c6-03c37219650d
# ╠═3a45b400-0df4-11eb-3e3c-41fed6f7a499
# ╟─40c1dc00-0df4-11eb-10a6-bf598057f7fb
