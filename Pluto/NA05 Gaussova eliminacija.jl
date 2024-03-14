### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 22dc4490-8d6b-11eb-20ce-ed84c7206167
using PlutoUI, Random, LinearAlgebra

# ╔═╡ 7b3f544c-1995-4eed-8818-37a9867375cc
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ 48167400-0d83-11eb-1c7d-359a2574c8b1
md"
# Gaussova eliminacija


# Osnovna ideja

Sustav $Ax=b$
se rješava u tri koraka (__bez pivotiranja__):

1.  $A=LU\quad$ (LU rastav, $O(\frac{2}{3}n^3)$ operacija),
2.  $Ly=b\quad$ (donje trokutrasti sustav, $n^2$ operacija),
3.  $Ux=y\quad$ (gornje torkutasti sustav, $n^2$ operacija).

S __pivotiranjem__ (redaka) imamo

1.  $PA=LU$
2.  $Ly=P b$
3.  $Ux=y$ 
"

# ╔═╡ 221d2474-de59-4042-918f-534305d8708f
md"""

## Primjeri

Sljedeći primjeri ukazuju na dva fenomena, jedan od kojih smo već vidjeli, dok drugi nismo. 
U ovom primjeru $\epsilon$ je vrijednost koju daje funkcija `eps()`.

Promotrimo sustav linearnih jednadžbi 

$\begin{aligned}\frac{\epsilon}{10} x_1 + x_2  = 1 \\ x_1 + x_2 = 2\end{aligned}$

Dobro približno rješenje je $x_1 = x_2 =1$. Koristimo proširenu matricu sustava:

$$
\begin{aligned}
& \left(\begin{array}{cc|c} \displaystyle\displaystyle\frac{\epsilon}{10} & 1 & 1 \cr 1 & 1 & 2 \end{array}\right) \sim
\left(\begin{array}{cc|c} \displaystyle \displaystyle\frac{\epsilon}{10} & 1 & 1 \cr 0 & 1-\displaystyle\displaystyle\frac{10}{\epsilon} & 2-\displaystyle\displaystyle\frac{10}{\epsilon}\end{array}\right)
\approx \left(\begin{array}{cc|c}  \displaystyle\displaystyle\frac{\epsilon}{10} & 1 & 1 \cr 0 & -\displaystyle\displaystyle\frac{10}{\epsilon} & -\displaystyle\displaystyle\frac{10}{\epsilon}\end{array}\right).
\end{aligned}$$

Zadnja transformacija je zaokruživanje do na točnost stroja $\epsilon$.
Vrlo značajni "1" i "2" u zadnjem retku su _nestali prilikom zaokruživanja_ !

Rješavanje zadnjeg trokutastog sustava daje $x_1 = 0$, $x_2 = 1$. Uvrštavanje u originalni sustav daje $x_1+x_2 = 1$, što je "vrlo netočno".
"""

# ╔═╡ 89291e44-4ea6-4e74-99bd-a30e8ee4d895
[eps()/10 1;0 1-10/eps()]\[1; 2-10/eps()]

# ╔═╡ e128967d-ee0f-4f54-bd35-8b8a55a9fbb3
cond([eps()/10 1;0 1-10/eps()])

# ╔═╡ 7b998fd0-87d4-11eb-286a-fd35bcef317c
2-10/eps()==-10/eps()

# ╔═╡ ed629240-6b9a-4a85-b412-97a06565a9cf
md"""
I ovdje postoji "rješenje" koje se zove __parcijalno pivotiranje__. Stavimo apsolutno najveći element koji još nije poništen u promatranom stupcu na pivotnu odnosno dijagonalnu poziciju:

$$
\begin{aligned}
&\left(\begin{array}{cc|c}     1 & 1 & 2 \cr \displaystyle\frac{\epsilon}{10} & 1 & 1 \end{array}\right) \sim
\left(\begin{array}{cc|c}     1 & 1 & 2 \cr 0                   & 1-\displaystyle\frac{\epsilon}{10}&1-\displaystyle\frac{2\epsilon}{10}\end{array}\right)
\approx   \left(\begin{array}{cc|c}     1 & 1 & 2 \cr 0                   & 1                    &1                    \end{array}\right)
\end{aligned}$$
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

Točno rješenje je $x_1=x_2=1$, ali zbog zaokruživanja dobije $x_1 =0$, $x_2 =2$. __Objasnite!__
Niti jedan trik kojeg smo vidjeli ne daje točno rješenje. 
"""

# ╔═╡ c471e93a-5c8a-433f-880a-9e26788fa601
[1+2*eps() 1+2*eps(); 1+eps() 1]\[2+4*eps(); 2+eps()]

# ╔═╡ e1a789e2-7402-49a5-b0f2-ed1666f1ca4c
cond([1+2*eps() 1+2*eps(); 1+eps() 1])

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
# LU rastav

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

# ╔═╡ 53d8d8e0-0d7f-11eb-358f-3f001e5e7d78
Random.seed!(127);

# ╔═╡ 19567c60-0d82-11eb-3405-8d0312a34b5f
begin
	n=6
	A=randn(n,n)
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
# Trokutasti sustavi
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
# ?lu

# ╔═╡ 4225c750-b668-4331-b6b8-0509635e69c6
md"""
# Brzina

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

# ╔═╡ cb2c4741-63a2-4969-b7ef-88ca6ecf1563
M=lu(A,Val(false))

# ╔═╡ de4e00d5-7705-4b70-9581-8b8a32059efa
M₁=lu(A)

# ╔═╡ 56cf8ac9-549e-4735-94fb-1e054b9a896e
M₁.L*M₁.U-A[M₁.p,:]

# ╔═╡ fbccb866-56d2-45cf-96ba-65f7afb4ad12
M₁.L*M₁.U-M₁.P*A

# ╔═╡ 195b7bf3-058d-4cc2-8759-5868d2bd717b
M₁.p

# ╔═╡ 6e25952c-a90d-495c-bc6e-4f67188ebe9f
M₁.P

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
mylu₁(A₁)

# ╔═╡ 09b70bcd-43ad-46b2-9664-2809351f9f70
md"""
## Blok inačica

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
	k,l=32,16 # 2,4   # 32,16
	Ab=[rand(k,k) for i=1:l, j=1:l]
end

# ╔═╡ 21e796b0-0dfb-11eb-2493-6fe0d3c70180
Ab[2,1]

# ╔═╡ 24c23190-0d85-11eb-382d-2536a040dfc8
A₀=mylu₂(Ab);

# ╔═╡ 3ad22ca1-5497-4c98-9088-517c322f06d8
A₀

# ╔═╡ 3d486c20-0d85-11eb-2c81-a7d99c190907
begin
	# Provjera
	U₀=triu(A₀)
	L₀=tril(A₀)
	for i=1:maximum(size(L₀))
		L₀[i,i]=one(L₀[1,1])
	end
end

# ╔═╡ e79269d2-532b-4929-bf19-030165807b3b
U₀

# ╔═╡ e68e0a50-87de-11eb-1348-5135365bc28f
L₀

# ╔═╡ 9615e1c0-0d85-11eb-161d-39a7eff4046f
Rezidual=L₀*U₀-Ab

# ╔═╡ d57feaf0-7000-454d-bd3d-e1ba52e36b0e
norm(Rezidual)

# ╔═╡ aa933e40-0d85-11eb-2141-d36ff3fc471b
# Pretvaranje blok matrice u običnu
unblock(A) = mapreduce(identity, hcat, [mapreduce(identity, vcat, A[:,i]) for i = 1:size(A,2)])

# ╔═╡ c64b1170-87de-11eb-1bf1-8f695cd2ac73
norm(Rezidual)/norm(Ab)

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
# Pivotiranje

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
# ?lu

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

# ╔═╡ a05104c4-48ff-4a64-a5f7-8b0557d2019b
F₄.P

# ╔═╡ 441c6f6c-187f-4d15-b9b0-1ae91acbf8f4
F₄.p

# ╔═╡ 74b55dc0-0def-11eb-37e2-71ba59aa8295
F₄.L*F₄.U==A[F₄.p,:]

# ╔═╡ 16800a1d-8ab9-48dc-ad0e-9cb11b55b00c
F₄.L*F₄.U-F₄.P*A

# ╔═╡ 90e4a320-0def-11eb-189e-7fc9c7b53942
# Javljaju se greške zaokruživanja
F₄.L*F₄.U-A[F₄.p,:]

# ╔═╡ 630a82aa-6998-4325-a0c9-d44f60c0df31
md"""
## Potpuno pivotiranje

Sljedeći program računa Gaussovu eliminaciju koristeći __potpuno pivotiranje__ - u svakom koraku 
se retci i stupci zamijene tako da se na pivotnu poziciju dovede element koji ima najveću 
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
# Točnost

Neka je zadan sustav $Ax=b$, pri čemu je matrica $A$ regularna.

Da bi primijenili koncepte iz bilježnice 
[NA04 Pogreska unatrag i stabilni algoritmi](https://ivanslapnicar.github.io/Numericka_analiza/NA04%20Pogreska%20unatrag%20i%20stabilni%20algoritmi.jl.html), potrebno je:

1. napraviti teoriju smetnje za dani problem
2. analizirati pogreške algoritma (Gaussove eliminacije)

## Teorija smetnje

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

# ╔═╡ b085eb50-87e1-11eb-32eb-59dd03302d32
norm(δb₃)/norm(b₃)

# ╔═╡ 5179d372-0df1-11eb-0183-8bcc73149584
begin
	δA₃=[-0.001 0;0 0]
    x₄=(A₃+δA₃)\b₃
end

# ╔═╡ 69b8cbd0-0df1-11eb-3fcd-a9f0865efdce
cond(A₃), norm(δA₃)/norm(A₃), norm(x₄-x₃)/norm(x₃)

# ╔═╡ a47802e0-0df1-11eb-3f9f-2fe1ebc781fd
md"""
## Pogreška Gaussove eliminacije

Prema [Matrix Computations, poglavlje 3.3](https://books.google.hr/books?id=X5YfsuCWpxMC&printsec=frontcover&hl=hr#v=onepage&q&f=false), za izračunate faktore
$\hat L$ i $\hat U$ vrijedi

$$\hat L\cdot \hat U = A+\delta A$$

gdje je (nejednakost se čita po elementima matrica, $\varepsilon$ je sada točnost stroja)

$$| \delta A|\leq 2(n-1) \varepsilon (|A|+|\hat L| \cdot |\hat U|) +O(\varepsilon^2).$$

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

"""

# ╔═╡ ecc4e74f-7fa4-4a88-9f89-eb2c0643adef
md"
### Vandermondeova matrica
"

# ╔═╡ f1c74560-0df1-11eb-19a7-c9ad6aed7410
begin
	nᵥ=10
	v=randn(nᵥ)
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

# ╔═╡ 16ce5f29-1332-48d4-ad18-56602c72a00f


# ╔═╡ 50ed08c1-391f-4baa-8cfa-54db04038fb1
md"""
## Umjetno loša kondicija
"""

# ╔═╡ a496c05e-0def-11eb-0ae1-83f3cdccf36e
Aᵤ=[1 1; 1 2]

# ╔═╡ aeee4670-0df2-11eb-0ef9-0bb353d8ebfe
bᵤ=[1;3]

# ╔═╡ aef0b770-0df2-11eb-3f66-8d09a5970a49
xᵤ=Aᵤ\bᵤ

# ╔═╡ aef1a1d0-0df2-11eb-0aeb-c5f670c48b32
xᵤ,cond(Aᵤ,1)

# ╔═╡ af144500-0df2-11eb-37e0-af9d5cd06c65
A₅=[1e-4 1e-4;1 2]
# A₅=[1 1;1e-4 2e-4]

# ╔═╡ af16b600-0df2-11eb-2852-7f1e61314674
b₅=[1e-4;3]
# b₅=[1;3e-4]

# ╔═╡ af2b7680-0df2-11eb-2ff2-6b188cbd6b7f
x₅=A₅\b₅

# ╔═╡ af2dc072-0df2-11eb-2958-19fc8d01e81d
x₅,cond(A₅),xᵤ-x₅

# ╔═╡ bfc0ea15-556f-4109-b626-cb724ee14bfd
md"""
## Procjena kondicije

Računanje kondicije prema definiciji $\kappa(A)=\|A\| \cdot \|A^{-1}\|$ zahtijeva računanje inverzne matrice, za što je potrebno $O(n^3)$ operacija. To je isti red veličine operacija koji je potreban za rješavanje zadanog sustava. Prilikom rješavanja sustava na raspolaganju su nam trokutasti faktori $L$ i $U$, što se može iskoristiti kako bi se kondicija približno izračunala u $O(n^2)$ operacija. 
Detalji se nalaze u [Matrix Computations, poglavlje 3.5.4](https://books.google.hr/books/about/Matrix_Computations.html?id=X5YfsuCWpxMC&redir_esc=y). 
LAPACK rutina 
[dtrcon.f](http://www.netlib.org/lapack/explore-html/d9/d84/dtrcon_8f_source.html) računa inverz približne kondicije trokutaste matrice.

Izračunajmo približnu kondiciju Vandermondeove matrice iz prethodnog primjera.
 
"""

# ╔═╡ 94525fd9-6847-4ff0-8662-6f6bb23f4083
# Norma inverza gornje trokutaste matrice
# GVL, str. 141, (3.5.6)
function mycondU(U::Array{T}) where T
	n=size(U,1)
	p=zeros(n)
	y=zeros(n)
	d=zeros(n)
	for k=n:-1:1
	    d[k]= p[k] < zero(T) ? one(T) :  -one(T)
		y[k]=(d[k]-p[k])/U[k,k]
		p[1:k-1]+=y[k]*U[1:k-1,k]
	end
	# Test
	println(norm(U*y-d))
	println(d)
	return norm(y,Inf)
end

# ╔═╡ 8707357c-84c9-4d30-aa81-1172d7ac715e
# ?LAPACK.trcon!

# ╔═╡ e38b8370-0df2-11eb-0ed5-ab750f73de17
begin
	Fᵥ=lu(V)
	cond(V,1), cond(Fᵥ.L,1), cond(Fᵥ.U,1)
end

# ╔═╡ 5499cc64-0f84-4ccb-adcc-f015b17d1284
mycondU(Fᵥ.U) # *norm(U,Inf)

# ╔═╡ dc8009e5-99ae-415f-97aa-f7144f822057
norm(Fᵥ.L)*norm(Fᵥ.U)/norm(V)

# ╔═╡ a1969708-ddfc-49df-aaa3-9ff994e21f72
Fᵥ.U

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

S druge strane, ako relativni rezidual ima malu normu, tada je _rješenje izračunato stabilno_. Naime, za (ovdje koristimo 2-normu)

$$\delta A=\frac{r\hat x^T}{\|\hat x\|^2_2}$$

vrijedi (1):

$$b-(A+\delta A)\hat x=(b-A\hat x)-\delta A \hat x = r-\frac{r\hat x^T \hat x}{\|\hat x\|^2_2}
= r-\frac{r \|\hat x^T \hat x\|_2}{\|\hat x\|^2_2}=r-r=0.$$

Također vrijedi

$$\frac{\|  \delta A\|_2}{\|A \|}  \leq  \frac{\|r\|_2\|\hat x \|_2}{\| A\| \cdot \|\hat x \|^2_2}=
\frac{\|r\|_2}{\| A\| \cdot \|\hat x \|_2}.$$

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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
PlutoUI = "~0.7.54"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.1"
manifest_format = "2.0"
project_hash = "519c88b955a16a6f52e4beee9c744049f942c2fe"

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

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

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
# ╠═22dc4490-8d6b-11eb-20ce-ed84c7206167
# ╠═7b3f544c-1995-4eed-8818-37a9867375cc
# ╟─48167400-0d83-11eb-1c7d-359a2574c8b1
# ╟─221d2474-de59-4042-918f-534305d8708f
# ╠═89291e44-4ea6-4e74-99bd-a30e8ee4d895
# ╠═e128967d-ee0f-4f54-bd35-8b8a55a9fbb3
# ╠═7b998fd0-87d4-11eb-286a-fd35bcef317c
# ╟─ed629240-6b9a-4a85-b412-97a06565a9cf
# ╠═b8b17cd3-c828-438d-a2c0-13b1965ed778
# ╟─eb5aab00-0d73-11eb-2f45-771b2e23a5e3
# ╠═c471e93a-5c8a-433f-880a-9e26788fa601
# ╠═e1a789e2-7402-49a5-b0f2-ed1666f1ca4c
# ╠═9983cc49-9398-4772-9d04-3f7bbf2b47a1
# ╟─de3e2152-7ec2-4366-a6fa-dc1476d19480
# ╟─c9bc269a-a306-4a35-acd8-aad3de58f56a
# ╠═66fcc372-4f27-4092-9552-8eb6e863bd4a
# ╠═53d8d8e0-0d7f-11eb-358f-3f001e5e7d78
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
# ╠═cb2c4741-63a2-4969-b7ef-88ca6ecf1563
# ╠═de4e00d5-7705-4b70-9581-8b8a32059efa
# ╠═56cf8ac9-549e-4735-94fb-1e054b9a896e
# ╠═fbccb866-56d2-45cf-96ba-65f7afb4ad12
# ╠═195b7bf3-058d-4cc2-8759-5868d2bd717b
# ╠═6e25952c-a90d-495c-bc6e-4f67188ebe9f
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
# ╠═3ad22ca1-5497-4c98-9088-517c322f06d8
# ╠═e79269d2-532b-4929-bf19-030165807b3b
# ╠═3d486c20-0d85-11eb-2c81-a7d99c190907
# ╠═e68e0a50-87de-11eb-1348-5135365bc28f
# ╠═9615e1c0-0d85-11eb-161d-39a7eff4046f
# ╠═d57feaf0-7000-454d-bd3d-e1ba52e36b0e
# ╠═aa933e40-0d85-11eb-2141-d36ff3fc471b
# ╠═c64b1170-87de-11eb-1bf1-8f695cd2ac73
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
# ╠═a05104c4-48ff-4a64-a5f7-8b0557d2019b
# ╠═441c6f6c-187f-4d15-b9b0-1ae91acbf8f4
# ╠═74b55dc0-0def-11eb-37e2-71ba59aa8295
# ╠═16800a1d-8ab9-48dc-ad0e-9cb11b55b00c
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
# ╠═b085eb50-87e1-11eb-32eb-59dd03302d32
# ╠═5179d372-0df1-11eb-0183-8bcc73149584
# ╠═69b8cbd0-0df1-11eb-3fcd-a9f0865efdce
# ╟─a47802e0-0df1-11eb-3f9f-2fe1ebc781fd
# ╟─ecc4e74f-7fa4-4a88-9f89-eb2c0643adef
# ╠═f1c74560-0df1-11eb-19a7-c9ad6aed7410
# ╠═026d0d00-0df2-11eb-26fd-cbc13048e56c
# ╠═2f58af40-0df2-11eb-28e4-5b1911f53b83
# ╠═3da4f680-0df2-11eb-23d4-f3fb28cdc8e7
# ╠═439224f0-0df2-11eb-2e57-539a3470de32
# ╠═50da42a0-0df2-11eb-2ded-c52a76acc155
# ╠═56e4bd10-0df2-11eb-0152-556ef692a70e
# ╠═16ce5f29-1332-48d4-ad18-56602c72a00f
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
# ╠═94525fd9-6847-4ff0-8662-6f6bb23f4083
# ╠═5499cc64-0f84-4ccb-adcc-f015b17d1284
# ╠═8707357c-84c9-4d30-aa81-1172d7ac715e
# ╠═e38b8370-0df2-11eb-0ed5-ab750f73de17
# ╠═dc8009e5-99ae-415f-97aa-f7144f822057
# ╠═a1969708-ddfc-49df-aaa3-9ff994e21f72
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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
