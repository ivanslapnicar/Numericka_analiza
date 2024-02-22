### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 53764ad6-7bee-4307-a599-49627451b86f
using PlutoUI, LinearAlgebra

# ╔═╡ dd774cfc-dfb2-4973-939b-0182ee63159c
TableOfContents(title="📚 Sadržaj", aside=true)

# ╔═╡ 1db27bd5-d11b-4cee-8f5e-53427cf82786
md"""
# Množenje matrica



Matrice možemo množiti na   __tri različita načina__:

$$
\begin{aligned}
\begin{bmatrix}
1& 2& 3\\
  4 &5 &6\\
  7 &8& 9
\end{bmatrix}
\begin{bmatrix}
   1&  2&  0\\
   4&  3&  2\\
   1& -1&  1
 \end{bmatrix}&=
 \begin{bmatrix}
   (1\cdot 1 + 2\cdot 4+3\cdot 1) & 1\cdot 2+2\cdot 3+3\cdot(-1))&
    (1\cdot 0 + 2\cdot 2 + 3\cdot 1)\\
   (4\cdot 1+5\cdot 4+6\cdot 1) & (4\cdot 2+5\cdot 3+6\cdot (-1)) &
   (4\cdot 0+5\cdot 2+6\cdot 1)\\
   (7\cdot 1+8\cdot 5+9\cdot 1) & (7\cdot 2+8\cdot 3+9\cdot (-1)) &
   (7\cdot 0+8\cdot 2+9\cdot 1)
 \end{bmatrix} \\ \\
\begin{bmatrix}
  1& 2& 3\\
  4 &5 &6\\
  7 &8& 9
\end{bmatrix}
\begin{bmatrix}
   1&  2&  0\\
   4&  3&  2\\
   1& -1&  1
 \end{bmatrix}& =
 \begin{bmatrix} 1\\ 4\\ 7 \end{bmatrix}
 \begin{bmatrix} 1&2&0 \end{bmatrix} +
 \begin{bmatrix} 2\\ 5\\ 8 \end{bmatrix}
 \begin{bmatrix} 4&3&2 \end{bmatrix} +
 \begin{bmatrix} 3\\ 6\\ 9 \end{bmatrix}
 \begin{bmatrix} 1&-1&1 \end{bmatrix}\\ \\
\begin{bmatrix}
  1& 2& 3\\
  4 &5 &6\\
  7 &8& 9
\end{bmatrix}
\begin{bmatrix}
   1&  2&  0\\
   4&  3&  2\\
   1& -1&  1
 \end{bmatrix}&=
 \begin{bmatrix}
   1 \begin{bmatrix} 1\\ 4\\ 7\end{bmatrix}+
   4 \begin{bmatrix} 2\\ 5\\ 8\end{bmatrix}+
   1 \begin{bmatrix} 3\\ 6\\ 9\end{bmatrix} &
   2 \begin{bmatrix} 1\\ 4\\ 7\end{bmatrix}+
   3 \begin{bmatrix} 2\\ 5\\ 8\end{bmatrix}+
   (-1) \begin{bmatrix} 3\\ 6\\ 9\end{bmatrix} &
   0 \begin{bmatrix} 1\\ 4\\ 7\end{bmatrix}+
   2 \begin{bmatrix} 2\\ 5\\ 8\end{bmatrix}+
   1 \begin{bmatrix} 3\\ 6\\ 9\end{bmatrix}
 \end{bmatrix}
\end{aligned}$$

Formule se razlikuju u načinu pristupa memoriji i, posljedično, u brzini i pogreškama zaokruživanja.
"""

# ╔═╡ d4ca4508-4d02-41f3-b343-e7d9e3d9349c
begin
	n=5
	A=ones(Int,n,n)
	B=2*ones(Int,n,n)
	A,B
end

# ╔═╡ 0f21c6e6-4aa9-4eb3-8c5e-c5366dd8f2a0
begin
	# Standardna formula
	C=zero(A)
	for i=1:n
	    for j=1:n
	        for k=1:n
	            C[i,j]=C[i,j]+A[i,k]*B[k,j]
	        end
	        @show i,j,C[i,j] # vidljivo u terminalu
	    end
	end
	C
end

# ╔═╡ c56ebccc-e0c3-4926-a7f2-6363bd0cc563
# Unutrašnja petlja: _dot
for i=1:n
    for j=1:n
        (A[i,:]⋅B[:,j])[] # vidljivo u terminalu
    end
end

# ╔═╡ 2ddf86e1-0751-4c0e-aee5-c3b1d4bed2fa
begin
	# Unutrašnja petlja: _syrk
	C₁=zero(A)
	for i=1:n
	    C₁=C₁+A[:,i]*B[i,:]'
	    println(C₁) # vidljivo u terminalu
	end
	C₁
end

# ╔═╡ 7d8c987a-3a4e-48d9-8020-2170090f9934
begin
	# Unutrašnja petlja: _axpy
	C₃=zero(A)
	for i=1:n
	    for k=1:n
	        C₃[:,i]+=A[:,k]*B[k,i]
	    end
	    @show C₃
	end
end

# ╔═╡ 8cef7e23-c539-4508-9176-37132955257d
md"""
`_axpy` je $y=\alpha x+y$
"""

# ╔═╡ 05adbae2-88a2-49ed-9ca6-129757f94266
md"""
# Basic Linear Algebra Subroutines - [BLAS](http://www.netlib.org/blas/)

Pogledajmo, na primjer, [`ddot.f`](http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f_source.html) - uočite _loop unrolling_ :
"""

# ╔═╡ 6a085347-8f1d-4dc4-8742-774c4e279cfc
md"""
$x\cdot y=\sum_i x_i\cdot y_i=x_1y_1+x_2y_2+\cdots +x_ny_n$
"""

# ╔═╡ b3f1858f-b1ec-4f60-8343-b4b12d2c60c7
md"""
```
DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     forms the dot product of two vectors.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      DDOT = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DTEMP = DTEMP + DX(IX)*DY(IY)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,5)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF (N.LT.5) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     +            DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END
```

"""

# ╔═╡ 656975d5-874e-4859-88a2-d9cd60cc871d


# ╔═╡ 85135806-d982-4ca1-b852-fc53cfd0903e
md"""
# Brzina računanja
"""

# ╔═╡ 25d321c6-26bc-4d2d-ab74-8c048560e57c
function AB(A::Array{T},B::Array{T}) where T
	# Množenje dvije kvadratne matrice
    n=size(A,1)
    C=[zero(A[1,1]) for i=1:n,j=1:n]
    for i=1:n
        for j=1:n
            for k=1:n
                C[i,j]+=A[i,k]*B[k,j]
            end
        end
    end
    C
end

# ╔═╡ 533c1833-8b02-4d45-b1fc-ad0255d0380a
begin
	import Random
	Random.seed!(123)
	n₁=1024
	A₁=rand(n₁,n₁)
	B₁=rand(n₁,n₁)
end

# ╔═╡ 1783e655-bc2f-4442-a772-2b98f5ab2e96
# Izvršite dva puta, drugo mjerenje je relevantno.
@time A₁*B₁

# ╔═╡ 67ae15b6-7c2c-40e1-8b68-c562bb9264c4
@which A₁*B₁

# ╔═╡ 698da56b-ed97-4e26-8456-403080e095e3
# ~ 5 Gfglops
operacija_u_sekundi=(2*n₁^3)/0.1

# ╔═╡ ea827322-031d-42b5-aafd-10c2d12fa469
md"""
__Zadatak.__ Izračunajte najveći $n$ za koji tri kvadratne $n\times n$ matrice `Float64` brojeva stanu u RAM vašeg računala. Koliko dugo će trajati množenje?
"""

# ╔═╡ d1dccfbf-7878-4b48-a376-f0bf9b61ad59
# Izvršite dva puta, drugo mjerenje je relevantno.
# Naš program je puno sporiji!!
@time AB(A₁,B₁);

# ╔═╡ 7dc0b94c-4955-4eef-85ee-c5c65d7ec6e2
12/0.1

# ╔═╡ 04528cf4-0b12-4164-b70b-772202176392
function ABa(A::Array{T},B::Array{T}) where T
	# Množenje dvije kvadratne matrice, unutrašnja petlja axpy
    n=size(A,1)
    C=[zero(A[1,1]) for i=1:n,j=1:n]
    for i=1:n
        for k=1:n
			for l=1:n
            	C[l,i]+=A[l,k]*B[k,i]
			end
        end
    end
end

# ╔═╡ 93d73551-f2d6-4364-9848-dcff2c847846
@time ABa(A₁,B₁);

# ╔═╡ 2015d3fc-37d2-43ee-99da-fe7029b73ee4
md"""
# Blok varijanta

Radi ubrzavanja našeg programa boljim korištenjem cache memorije trebamo računati s blok matricama (BLAS 3).

BLAS razina  |  funkcija | memorija | # operacija
:---|:-----|:---|:---
1 | $y= \alpha x+y$  | $2n$ | $2n$
2 | $y=\alpha Ax+\beta y$ | $n^2$ | $2n^2$
3 | $C=\alpha A\cdot B+\beta C$| $3n^2$ | $2n^3$

__Zadatak.__ Proučite subrutine `daxpy.f`, `dgemmv.f` i `dgemm.f`.
"""

# ╔═╡ 6fc4fe72-e7ac-4c14-aa9c-37a5459d7004
begin
	# Probajte k,l=32,16 i k,l=64,8. Izvršite dva puta.
	k,l=64,16
	Ab=[rand(k,k) for i=1:l, j=1:l]
	Bb=[randn(k,k) for i=1:l, j=1:l]
end

# ╔═╡ 669a3f33-1cfa-4ef6-930f-38771afb18c0
Ab[1,1]

# ╔═╡ 9bc2dc46-206d-46c5-b3c3-17d398d2cf6d
# This is considerably faster.
@time AB(Ab,Bb)

# ╔═╡ 4a639621-805e-4bc3-8ca1-7ca12e37f855
@time ABa(Ab,Bb)

# ╔═╡ dccef80c-4a20-49f0-a45d-073f62777432
@time Ab*Bb

# ╔═╡ 30bc11d3-78cc-422b-9792-685a388c386d
@which Ab*Bb

# ╔═╡ a7de05aa-0ed1-4e11-88de-de70d629d30b
md"""
Za još bolje ubrzanje trebali bi koristiti sve jezgre kao što to radi `*`.

Za pokretanje Julia-e u višenitnom okruženju u Linuxu upišite liniju

    export JULIA_NUM_THREADS=`nproc`

u datoteku `.bashrc` file. U Windows-ima, postavite "environment variable" JULIA_NUM_THREADS na broj procesora (run msinfo32).

"""

# ╔═╡ 209d4041-a0d9-455d-a28d-1a7cc632082f
Threads.nthreads()

# ╔═╡ 3b096734-bd85-4e21-ad81-6c1ed99e2f43
md"""
# Točnost računanja

## Osnovne operacije

Za operacije $\odot \in  \{+,*,/\}$ vrijedi (zbog jednostavnosti koristimo $\varepsilon$ umjesto $\varepsilon_M$)

$$
fl(a \odot b)=(1+\varepsilon_\odot)(a\odot b),\qquad |\varepsilon_\odot|\leq \varepsilon.$$


## Zbrajanje

Ukoliko potpuno točno zbrojimo dva broja koja imaju pogreške iz prethodnih izračuna, 
jednakost

$$a(1+\varepsilon_a)+b(1+\varepsilon_b)= (1+\varepsilon_{ab})(a+b)$$

daje

$$
\varepsilon_{ab}=\frac{a\varepsilon_a+b\varepsilon_a}{a+b},$$

odnosno

$$
|\varepsilon_{ab}|\leq \varepsilon\, \frac{|a|+  |b|}{|a+b|}.$$

Ako su  $a$ i $b$ are veliki brojevi različitih predznaka i $a-b$ je sitno, pogreška može biti ogromna  (__katastrofalno kraćenje__ ili _catastrophic cancelation_).

## Skalarni (dot) produkt

Za vektore $x$ i $y$ iz $\mathbb{R}^n$, rekurzivna primjena prethodne formule daje apsolutnu pogrešku

$$|fl(x\cdot y)-x\cdot y|\leq O(n\varepsilon) |x|\cdot |y| \tag{1}$$

i relativnu pogrešku

$$\frac{|fl(x\cdot y)-x\cdot y|}{|x\cdot y|}\leq O(n\varepsilon) \frac{|x|\cdot |y|}{|x\cdot y|}$$

Ako su vektori $x$ i $y$ gotovo okomiti, relativna pogreška može biti velika.
"""

# ╔═╡ dbefde6e-998a-4b8e-b923-884d59e172f1
md"
Dokažimo  __egzaktnu__ varijantu ocjene (1). 

__Lema 1__ [MC, Poglavlje 2.7] Ako je $n \varepsilon \leq 0.01$, onda je 

$$
|fl(x\cdot y)-x\cdot y| \leq 1.01 n\varepsilon  |x|\cdot |y|. \tag{2}$$

Za dokaz Leme 1 potreban bnam je sljedeći rezultat [ASNA, p. 68]:

__Lema 2__ Ako je $|\delta_i|\leq \varepsilon$ i $n\varepsilon \leq 0.01$, onda je 

$$
\prod_{i=1}^n (1+\delta_i) = 1+\eta_n, \qquad |\eta_n|\leq 1.01n\varepsilon.$$  

_Dokaz:_ U dokazu koristimo činjenicu da je $1+x\leq e^x$ za mali $x\geq 0$ i razvoj funkcije $e^x$ u Taylorov red.  $\square$
"

# ╔═╡ eabdfa6d-b09b-4322-bfc8-d9871bf596e8
md"
_Dokaz Leme 1:_  Označimo

$$
s_p=fl(\sum_{k=1}^p x_k\,  y_k).$$

Ako koristimo standardni algoritam, za svaki $p=2:n$ vrijedi

$$
s_p=fl(s_{p-1}+fl(x_p\, y_p))=(s_{p-1}+x_p\, y_p(1+\delta_p))(1+\epsilon_p),\quad
|\delta_p|,|\epsilon_p|\leq \varepsilon.$$

Posebno,

$$
s_1=x_1\, y_1(1+\delta_1),\qquad |\delta_1|\leq \varepsilon,$$

$$
s_2=fl(s_1+fl(x_2\, y_2))=(s_1+x_2\, y_2(1+\delta_2))(1+\epsilon_2),\quad 
|\delta_2|,|\epsilon_2|\leq \varepsilon,$$

$$
s_3=fl(s_2+fl(x_3\, y_3))=(s_2+x_3\, y_3(1+\delta_3))(1+\epsilon_3),\quad 
|\delta_3|,|\epsilon_3|\leq \varepsilon.$$
"

# ╔═╡ edc3dd30-4da8-48a3-83a2-0d6d3aa6debd
md"
Dakle,

$$
s_3=[(s_1+x_2\, y_2(1+\delta_2))(1+\varepsilon_2)+x_3\, y_3(1+\delta_3)](1+\varepsilon_3)$$

ili, uz $\epsilon_1=0$,

$$
s_3=x_1\, y_1(1+\delta_1)(1+\epsilon_1)(1+\epsilon_2)(1+\epsilon_3)+
x_2\, y_2 (1+\delta_2)(1+\epsilon_2)(1+\epsilon_3)+x_3\, y_3(1+\delta_3)(1+\epsilon_3).$$

Indukcijom slijedi

$$
s_n=fl(x\cdot y)=\sum_{k=1}^n x_k\, y_k (1+\gamma_k),$$

gdje je

$$
(1+\gamma_k)=(1+\delta_k)\prod_{j=k}^n (1+\epsilon_j).$$

Prema tome vrijedi

$$
fl(x\cdot y)=x\cdot y+\sum_{k=1}^n x_k\, y_k \gamma_k,$$

pa je

$$
|fl(x\cdot y)-x\cdot y|\leq\sum_{k=1}^n |x_k| |y_k| |\gamma_k|.$$

Prema Lemi 2 je $|\gamma_k|\leq 1.01 n\varepsilon$ za svaki $k$, pa ocjena (2) slijedi. $\square$
"

# ╔═╡ 7daface9-5571-41b1-a467-4b4936eeed08
md"""
Ocjenu iz Leme 2 možemop zapisati i na slijedeće načine:

$$
|fl(x\cdot y)-x\cdot y| \leq n\varepsilon  |x|\cdot |y| +O(\varepsilon^2)$$, 

$$
|fl(x\cdot y)-x\cdot y| \leq \phi(n)\varepsilon  |x|\cdot |y|,$$

gdje je $\phi(n)$ "umjerena" funkcija od $n$, i

$$
|fl(x\cdot y)-x\cdot y| \leq c\, n\varepsilon  |x|\cdot |y|$$

gdje je $c$ konstanta reda veličine jedan.
"""

# ╔═╡ 009af427-b54e-406b-a6f2-ba9f9589e2e5
md"
## `_axpy()`

Vrijedi

$$
fl(y+\alpha x)=y+\alpha x+z,\quad |z|\leq\varepsilon(|y|+2|\alpha x|)+O(\varepsilon^2).$$
"

# ╔═╡ 0e746d31-8024-4f04-acc8-b87d0ee9c388
md"
## Vanjski produkt

Vrijedi

$$
fl(C+uv^T)=C+uv^T+E,\quad |E|\leq \varepsilon (|C|+2|uv^T|)+O(\varepsilon^2).$$
"

# ╔═╡ 1abeedfa-325e-498a-99a5-365d8cb744c1
md"
## Spremanje matrice u memoriju

Zbog

$$
[fl(A)]_{ij}=a_{ij}(1+\epsilon_{ij}), \quad |\epsilon_{ij}|\leq \varepsilon,$$

vrijedi

$$|fl(A)-A|\leq \varepsilon|A|,$$

odnosno

$$\|fl(A)-A\|_1\leq \varepsilon\|A\|_1.$$
"

# ╔═╡ fbd0fb02-2060-437c-b1be-62cf12a530dc
md"
## Množenje matrice skalarom

$$
fl(\alpha A)=\alpha A+E,\quad |E|\leq \varepsilon |\alpha A|.$$
"

# ╔═╡ 7234f15a-b809-4cb4-a9ec-788cf8b09cb6
md"
## Zbrajanje matrica

$$fl(A+B)=(A+B)+E,\quad |E|\leq \varepsilon(|A|+|B|).$$
"

# ╔═╡ f36d6407-d921-4c80-8af3-3f8845be3f9c
md"
## Množenje matrica

Za sva tri načina množenja matrioca vrijedi

$$
fl(A\cdot B)=A\cdot B+E,\quad |E|\leq n\varepsilon |A|\cdot |B|+O(\varepsilon^2),$$

ili, u normi,

$$\|fl(A\cdot B) -A\cdot B\|_1 \leq n\varepsilon \|A\|_1 \|B\|_1 +O(\varepsilon^2),$$

ili

$$|fl(A\cdot B) -A\cdot B| \leq O(n\varepsilon) |A|\cdot |B|.$$
"

# ╔═╡ 1ff4efee-a106-46b4-855d-7b76265cb050
begin
	n₂=1_000_000
	x=rand(n₂)
	y=randn(n₂)
end

# ╔═╡ 1850009e-d407-46a7-be70-88fdb4a774be
d=x⋅y

# ╔═╡ b76c7b7c-2e9e-4985-a977-6ae52d0d4dcd
ab=abs.(x)⋅abs.(y)

# ╔═╡ 89027fca-0854-469b-947b-e68f62b3b8b7
begin
	# Check the solution using `BigFloat` numbers (70 decimal digits)
	xbig=map(BigFloat,x)
	ybig=map(BigFloat,y)
	ybig[1]
end

# ╔═╡ c1387f1e-f3ce-49fd-8daa-d270b4714987
dbig=xbig⋅ybig

# ╔═╡ 815fca39-6bf8-44fc-ab5c-5146e444eb35
# Absolute error
abserr=map(Float64,abs(d-dbig))

# ╔═╡ 24a79f79-ecc1-4dae-ab00-0e5be89d684f
# Relative error
relerr=map(Float64,abserr/d)

# ╔═╡ 23b3c19c-e8a0-4323-a501-f9b4fa53222b
begin
	n₃=256
	A₃=rand(n₃,n₃)
	B₃=randn(n₃,n₃)
end

# ╔═╡ b4831348-c018-4133-907b-1a34987a9518
begin
	Ab₃=map(BigFloat,A₃)
	Bb₃=map(BigFloat,B₃);
end

# ╔═╡ de747a60-3953-48ba-bf8c-4238beb07666
begin
	D₃=A₃*B₃
	# This lasts little longer
	Cb₃=Ab₃*Bb₃
	abserr₃=abs.(D₃-map(Float64,Cb₃))
end

# ╔═╡ d43cea1b-5574-4341-aa91-c7ea2efc5e55
abs.(A₃)*abs.(B₃)*n₃*eps()

# ╔═╡ 2b7f29bf-a440-4f07-868b-ae42f3624b05
md"""
Ovdje nije lako procijeniti relativnu pogrešku pa trebamo bolju mjeru.
U sljedećem predavanju objasnit ćemo _vektorske i matrične norme_. 
"""

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
# ╠═53764ad6-7bee-4307-a599-49627451b86f
# ╠═dd774cfc-dfb2-4973-939b-0182ee63159c
# ╟─1db27bd5-d11b-4cee-8f5e-53427cf82786
# ╠═d4ca4508-4d02-41f3-b343-e7d9e3d9349c
# ╠═0f21c6e6-4aa9-4eb3-8c5e-c5366dd8f2a0
# ╠═c56ebccc-e0c3-4926-a7f2-6363bd0cc563
# ╠═2ddf86e1-0751-4c0e-aee5-c3b1d4bed2fa
# ╠═7d8c987a-3a4e-48d9-8020-2170090f9934
# ╟─8cef7e23-c539-4508-9176-37132955257d
# ╟─05adbae2-88a2-49ed-9ca6-129757f94266
# ╟─6a085347-8f1d-4dc4-8742-774c4e279cfc
# ╟─b3f1858f-b1ec-4f60-8343-b4b12d2c60c7
# ╠═656975d5-874e-4859-88a2-d9cd60cc871d
# ╟─85135806-d982-4ca1-b852-fc53cfd0903e
# ╠═25d321c6-26bc-4d2d-ab74-8c048560e57c
# ╠═533c1833-8b02-4d45-b1fc-ad0255d0380a
# ╠═1783e655-bc2f-4442-a772-2b98f5ab2e96
# ╠═67ae15b6-7c2c-40e1-8b68-c562bb9264c4
# ╠═698da56b-ed97-4e26-8456-403080e095e3
# ╟─ea827322-031d-42b5-aafd-10c2d12fa469
# ╠═d1dccfbf-7878-4b48-a376-f0bf9b61ad59
# ╠═7dc0b94c-4955-4eef-85ee-c5c65d7ec6e2
# ╠═04528cf4-0b12-4164-b70b-772202176392
# ╠═93d73551-f2d6-4364-9848-dcff2c847846
# ╟─2015d3fc-37d2-43ee-99da-fe7029b73ee4
# ╠═6fc4fe72-e7ac-4c14-aa9c-37a5459d7004
# ╠═669a3f33-1cfa-4ef6-930f-38771afb18c0
# ╠═9bc2dc46-206d-46c5-b3c3-17d398d2cf6d
# ╠═4a639621-805e-4bc3-8ca1-7ca12e37f855
# ╠═dccef80c-4a20-49f0-a45d-073f62777432
# ╠═30bc11d3-78cc-422b-9792-685a388c386d
# ╟─a7de05aa-0ed1-4e11-88de-de70d629d30b
# ╠═209d4041-a0d9-455d-a28d-1a7cc632082f
# ╟─3b096734-bd85-4e21-ad81-6c1ed99e2f43
# ╟─dbefde6e-998a-4b8e-b923-884d59e172f1
# ╟─eabdfa6d-b09b-4322-bfc8-d9871bf596e8
# ╟─edc3dd30-4da8-48a3-83a2-0d6d3aa6debd
# ╟─7daface9-5571-41b1-a467-4b4936eeed08
# ╟─009af427-b54e-406b-a6f2-ba9f9589e2e5
# ╟─0e746d31-8024-4f04-acc8-b87d0ee9c388
# ╟─1abeedfa-325e-498a-99a5-365d8cb744c1
# ╟─fbd0fb02-2060-437c-b1be-62cf12a530dc
# ╟─7234f15a-b809-4cb4-a9ec-788cf8b09cb6
# ╟─f36d6407-d921-4c80-8af3-3f8845be3f9c
# ╠═1ff4efee-a106-46b4-855d-7b76265cb050
# ╠═1850009e-d407-46a7-be70-88fdb4a774be
# ╠═b76c7b7c-2e9e-4985-a977-6ae52d0d4dcd
# ╠═89027fca-0854-469b-947b-e68f62b3b8b7
# ╠═c1387f1e-f3ce-49fd-8daa-d270b4714987
# ╠═815fca39-6bf8-44fc-ab5c-5146e444eb35
# ╠═24a79f79-ecc1-4dae-ab00-0e5be89d684f
# ╠═23b3c19c-e8a0-4323-a501-f9b4fa53222b
# ╠═b4831348-c018-4133-907b-1a34987a9518
# ╠═de747a60-3953-48ba-bf8c-4238beb07666
# ╠═d43cea1b-5574-4341-aa91-c7ea2efc5e55
# ╟─2b7f29bf-a440-4f07-868b-ae42f3624b05
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
