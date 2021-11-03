### A Pluto.jl notebook ###
# v0.17.1

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
	C=zeros(Int,n,n)
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
	C₁=zeros(Int,n,n)
	for i=1:n
	    C₁=C₁+A[:,i]*B[i,:]'
	    println(C₁) # vidljivo u terminalu
	end
	C₁
end

# ╔═╡ da5dab82-8aa2-4bf3-9273-69bb3675508f
begin
	# Korak po korak
	C₂=zeros(n,n)
	for j=1:n
	    for k=1:n
	        for i=1:n
	            C₂[i,j]+=A[i,k]*B[k,j]
	            @show i,j,C₂[i,j]
	        end
	    end
	end
end

# ╔═╡ 7d8c987a-3a4e-48d9-8020-2170090f9934
begin
	# Unutrašnja petlja: _axpy
	C₃=zeros(n,n)
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

# ╔═╡ 85135806-d982-4ca1-b852-fc53cfd0903e
md"""
# Brzina računanja
"""

# ╔═╡ 25d321c6-26bc-4d2d-ab74-8c048560e57c
function AB(A::Array{T},B::Array{T}) where T
    (n,l)=size(A)
	(m,p)=size(B)
	if m!=l
		return("Error in dimensions")
	end
    C=zeros(n,p)
    for i=1:n
        for j=1:p
            for k=1:l
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
	n₁=512
	A₁=rand(n₁,n₁)
	B₁=rand(n₁,n₁)
end

# ╔═╡ 1783e655-bc2f-4442-a772-2b98f5ab2e96
# Izvršite dva puta, drugo mjerenje je relevantno.
@time A₁*B₁

# ╔═╡ 698da56b-ed97-4e26-8456-403080e095e3
operacija_u_sekundi=(2*n^3)/0.003

# ╔═╡ ea827322-031d-42b5-aafd-10c2d12fa469
md"""
__Zadatak.__ Izračunajte najveći $n$ za koji tri kvadratne $n\times n$ matrice `Float64` brojeva stanu u RAM vašeg računala. Koliko dugo će trajati množenje?
"""

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

# ╔═╡ be605463-3a97-4869-b9ba-426dc1822c96
function AB(A::Array{T},B::Array{T}) where T<:Array
	(n,l)=size(A)
	(m,p)=size(B)
	if m!=l
		return("Greška u dimenziji")
	end
	# Pretpostavimo blokove jednakih dimenzija
    C=[zero(A[1,1]) for i=1:n, j=1:p]
    for i=1:n
        for j=1:p
            for k=1:l
                C[i,j]=C[i,j]+A[i,k]*B[k,j]
            end
        end
    end
    C
end

# ╔═╡ d1dccfbf-7878-4b48-a376-f0bf9b61ad59
# Izvršite dva puta, drugo mjerenje je relevantno.
# Naš program je puno sporiji!!
@time AB(A₁,B₁);

# ╔═╡ 6fc4fe72-e7ac-4c14-aa9c-37a5459d7004
begin
	# Probajte k,l=32,16 i k,l=64,8. Izvršite dva puta.
	k,l=64,8
	Ab=[rand(k,k) for i=1:l, j=1:l]
	Bb=[randn(k,k) for i=1:l, j=1:l]
end

# ╔═╡ 47168e37-a8b8-4969-856e-c12486882a13
[zero(Ab[1,1]) for i=1:n, j=1:n]

# ╔═╡ 9bc2dc46-206d-46c5-b3c3-17d398d2cf6d
# This is considerably faster.
@time AB(Ab,Bb)

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

## Množenje matrica

Za matrice $A$ i $B$ reda $n$ formula (1) daje

$$|fl(A\cdot B) -A\cdot B| \leq O(n\varepsilon) |A|\cdot |B|.$$

"""

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

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

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
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

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
# ╠═53764ad6-7bee-4307-a599-49627451b86f
# ╠═dd774cfc-dfb2-4973-939b-0182ee63159c
# ╟─1db27bd5-d11b-4cee-8f5e-53427cf82786
# ╠═d4ca4508-4d02-41f3-b343-e7d9e3d9349c
# ╠═0f21c6e6-4aa9-4eb3-8c5e-c5366dd8f2a0
# ╠═c56ebccc-e0c3-4926-a7f2-6363bd0cc563
# ╠═2ddf86e1-0751-4c0e-aee5-c3b1d4bed2fa
# ╠═da5dab82-8aa2-4bf3-9273-69bb3675508f
# ╠═7d8c987a-3a4e-48d9-8020-2170090f9934
# ╟─8cef7e23-c539-4508-9176-37132955257d
# ╟─05adbae2-88a2-49ed-9ca6-129757f94266
# ╟─6a085347-8f1d-4dc4-8742-774c4e279cfc
# ╟─b3f1858f-b1ec-4f60-8343-b4b12d2c60c7
# ╟─85135806-d982-4ca1-b852-fc53cfd0903e
# ╟─25d321c6-26bc-4d2d-ab74-8c048560e57c
# ╠═533c1833-8b02-4d45-b1fc-ad0255d0380a
# ╠═1783e655-bc2f-4442-a772-2b98f5ab2e96
# ╠═698da56b-ed97-4e26-8456-403080e095e3
# ╟─ea827322-031d-42b5-aafd-10c2d12fa469
# ╠═d1dccfbf-7878-4b48-a376-f0bf9b61ad59
# ╟─2015d3fc-37d2-43ee-99da-fe7029b73ee4
# ╠═47168e37-a8b8-4969-856e-c12486882a13
# ╠═be605463-3a97-4869-b9ba-426dc1822c96
# ╠═6fc4fe72-e7ac-4c14-aa9c-37a5459d7004
# ╠═9bc2dc46-206d-46c5-b3c3-17d398d2cf6d
# ╟─a7de05aa-0ed1-4e11-88de-de70d629d30b
# ╠═209d4041-a0d9-455d-a28d-1a7cc632082f
# ╟─3b096734-bd85-4e21-ad81-6c1ed99e2f43
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
