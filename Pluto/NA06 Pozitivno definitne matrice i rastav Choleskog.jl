### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ c162d7c0-0efc-11eb-16ad-3d5ff641f6e0
using LinearAlgebra

# ╔═╡ 70942b12-0bd1-44e6-8e04-862d4948e8a0
md"""
# Pozitivno definitne matrice i rastav Choleskog


Matrica $A$ je __pozitivno definitna__ ako je simetrična, $A^T=A$, i ako su sve njene svojstvene vrijednosti pozitivne.

Pozitivno definitnu matricu možemo rastaviti (bez pivotiranja) kao

$$
A=L L^T$$

pri čemu je $L$ donje trokutasta matrica s pozitivnim dijagonalnim elementima. Taj rastav se zove
__rastav Choleskog__ (vidi [Numerička matematika, poglavlje 3.6](http://www.mathos.unios.hr/pim/Materijali/Num.pdf)).

Iz 

$$A=\begin{pmatrix}\alpha & a^T \cr a  & B \end{pmatrix}=
\begin{pmatrix} \beta & 0 \cr l & I \end{pmatrix}
\begin{pmatrix} \beta & l^T \cr 0 & C \end{pmatrix}
=\begin{pmatrix} \beta^2 & \beta l^T \cr l\beta & ll^T+ C\end{pmatrix}$$

slijedi

$$\beta=\sqrt{\alpha},\quad l=a\beta^{-1}=a\alpha^{-1/2},\quad C=B-ll^T=B-a\alpha^{-1}a^T.$$

Vidimo da treba vrijediti $\alpha>0$. Također, matrica $B$ je simetrična pa je takva i matrica $C$.

Indukcija daje sljedeći algoritam:

"""

# ╔═╡ a808c71a-97bb-4106-89fd-dcc0e902d2cb
function mychol(A₁::Matrix{T}) where T
    A=copy(A₁)
    n,m=size(A)
    for k=1:n
        A[k,k]=√A[k,k]
        for j=k+1:n
            A[k,j]=A[k,j]/A[k,k]
        end
        for j=k+1:n
            for i=k+1:j
                A[i,j]-=A[k,i]*A[k,j]
			end
        end
    end
    return triu(A)
end

# ╔═╡ ea4f0c1f-55e7-4423-a8df-b9de1889e511
begin
	A=rand(6,6)
	A=A*A'
end

# ╔═╡ 26ab5250-1377-11eb-2824-69646a3f5444
A==A'

# ╔═╡ 3c7d6960-1377-11eb-25b9-177b4c932ce6
eigvals(A)

# ╔═╡ dfea6e0a-a718-49ef-8e60-bd0ad6e204db
# Ugrađena funkcija
C=cholesky(A)

# ╔═╡ 396ba935-a1fc-4978-a64e-990bb1d78f44
# Izvadimo L iz strukture
L=C.U

# ╔═╡ e20bffc4-b17f-4a7c-98d9-182ce8eb8395
C.L

# ╔═╡ 1ff3e65c-8253-47d7-9caa-466a1359a434
# Residual 
L'*L-A

# ╔═╡ 286a74e5-6a3b-41ee-a227-cb5d88dd027e
# Naša funkcija
L₁=mychol(A)

# ╔═╡ 2e54b52f-2ba0-4de0-a7bf-73bf1c0d5b6f
L₁'*L₁-A

# ╔═╡ 7decf367-ee06-401c-b742-b910e4b647a9
L-L₁

# ╔═╡ af9b1080-8d4f-11eb-02e4-97071de47db5
md"
## Cholesky rastav s pivotiranjem

Zbog pozitivne definitnosti najveći apsolutno najveći element se nalazi na glavnoj dijagonali. Stoga se potpuno pivotiranje postiže odabirom najvećeg elementa na glavnoj dijagonali svake podmatrice.
"

# ╔═╡ bf6f0930-8d4f-11eb-28ce-7d1216521e8a
Cₚ=cholesky(A,Val(true))

# ╔═╡ d04b4a20-8d4f-11eb-030a-75b2da4964ed
Cₚ.P

# ╔═╡ efc9e18e-8d4f-11eb-0731-a7e83ca2fa59
Lₚ=Cₚ.L

# ╔═╡ e1a8f00e-8d4f-11eb-2dfa-1f02d40baf55
Lₚ*Lₚ'-Cₚ.P'*A*Cₚ.P

# ╔═╡ 13cb690e-8d50-11eb-0af4-9f2eb9733164
Lₚ*Lₚ'-A[Cₚ.p,Cₚ.p]

# ╔═╡ 8fe33ce0-8cdc-11eb-0b6a-c3632016d982
md"
## Blok Cholesky
"

# ╔═╡ 72bb3568-0528-49fa-96af-65c3e1da00cf
function mycholb(A₁::Matrix{T}) where T
    A=copy(A₁)
    n,m=size(A)
    for k=1:n
		C=cholesky(A[k,k])
        A[k,k]=C.U
        for j=k+1:n
            A[k,j]=C.L\A[k,j] # Rješavanje sustava
        end
        for j=k+1:n
            for i=k+1:j
                A[i,j]-=transpose(A[k,i])*A[k,j]
			end
        end
    end
    return triu(A)
end

# ╔═╡ ca0eb511-ba7b-42b4-89b3-8c9eb90f8fb7
# Generirajmo matricu
begin
	k,l=32,16   # 32,16
	Ab=[rand(k,k) for i=1:l, j=1:l]
	Ab=Ab'*Ab
end

# ╔═╡ 7f37a748-d31c-433f-bc66-9e7358762c1b
Ab[1,1]

# ╔═╡ dcdfcda1-74d5-42a9-bec8-0a09d87a43c5
Lb=mycholb(Ab)

# ╔═╡ 9024604d-88f7-4098-a0eb-474161bed4fe
md"
Vremena izvođenja našeg blok algoritma `mycholb()` i ugrađene funkcije `cholesky()` su slična.
"

# ╔═╡ 138188b6-1452-49cc-ba14-82309f92d0cc
Lb[1,1]

# ╔═╡ b51c6282-64f1-4b03-a2d4-009338d183d9
# Rezidual
norm(Lb'*Lb-Ab)

# ╔═╡ 73f28eab-b214-4477-af72-178016eed643
md"
Naša funkcija `mychol()` koja ne koristi blokove je 20-tak puta sporija!
"

# ╔═╡ 281cdf66-87e4-4cd1-893d-cc9779b1a813
# Prebacivanje blok matrice u standardni oblik
unblock(A) = mapreduce(identity, hcat, [mapreduce(identity, vcat, A[:,i]) for i = 1:size(A,2)])

# ╔═╡ 55f47608-829c-45ed-b6ac-eff79850c4d9
Ab₀=unblock(Ab);

# ╔═╡ 54fee016-a27f-4eda-916c-9f51413faf15
cholesky(Ab₀);

# ╔═╡ 761841cb-536f-4140-9854-d31b51669c2e
@time mychol(Ab₀)

# ╔═╡ bfd53356-1cce-47d6-b99e-4930a2a26b6d
md"
Ukoliko simetrična matrica nije pozitivno definitna, Cholesky rastav neće uspjeti  jer će u nekom trenutku trebati računati korijen negativnog broja.
"

# ╔═╡ b8584018-1371-4b03-b939-3b1f555632aa
b=randn(5,5)

# ╔═╡ cbaf7dad-5980-43a4-82c0-34bb946d9376
c=b+b'

# ╔═╡ 970c1a3b-205c-421c-a564-62b8beba625e
cholesky(c)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"
"""

# ╔═╡ Cell order:
# ╠═c162d7c0-0efc-11eb-16ad-3d5ff641f6e0
# ╟─70942b12-0bd1-44e6-8e04-862d4948e8a0
# ╠═a808c71a-97bb-4106-89fd-dcc0e902d2cb
# ╠═ea4f0c1f-55e7-4423-a8df-b9de1889e511
# ╠═26ab5250-1377-11eb-2824-69646a3f5444
# ╠═3c7d6960-1377-11eb-25b9-177b4c932ce6
# ╠═dfea6e0a-a718-49ef-8e60-bd0ad6e204db
# ╠═396ba935-a1fc-4978-a64e-990bb1d78f44
# ╠═e20bffc4-b17f-4a7c-98d9-182ce8eb8395
# ╠═1ff3e65c-8253-47d7-9caa-466a1359a434
# ╠═286a74e5-6a3b-41ee-a227-cb5d88dd027e
# ╠═2e54b52f-2ba0-4de0-a7bf-73bf1c0d5b6f
# ╠═7decf367-ee06-401c-b742-b910e4b647a9
# ╟─af9b1080-8d4f-11eb-02e4-97071de47db5
# ╠═bf6f0930-8d4f-11eb-28ce-7d1216521e8a
# ╠═d04b4a20-8d4f-11eb-030a-75b2da4964ed
# ╠═efc9e18e-8d4f-11eb-0731-a7e83ca2fa59
# ╠═e1a8f00e-8d4f-11eb-2dfa-1f02d40baf55
# ╠═13cb690e-8d50-11eb-0af4-9f2eb9733164
# ╟─8fe33ce0-8cdc-11eb-0b6a-c3632016d982
# ╠═72bb3568-0528-49fa-96af-65c3e1da00cf
# ╠═ca0eb511-ba7b-42b4-89b3-8c9eb90f8fb7
# ╠═7f37a748-d31c-433f-bc66-9e7358762c1b
# ╠═dcdfcda1-74d5-42a9-bec8-0a09d87a43c5
# ╟─9024604d-88f7-4098-a0eb-474161bed4fe
# ╠═138188b6-1452-49cc-ba14-82309f92d0cc
# ╠═b51c6282-64f1-4b03-a2d4-009338d183d9
# ╟─73f28eab-b214-4477-af72-178016eed643
# ╠═281cdf66-87e4-4cd1-893d-cc9779b1a813
# ╠═55f47608-829c-45ed-b6ac-eff79850c4d9
# ╠═54fee016-a27f-4eda-916c-9f51413faf15
# ╠═761841cb-536f-4140-9854-d31b51669c2e
# ╟─bfd53356-1cce-47d6-b99e-4930a2a26b6d
# ╠═b8584018-1371-4b03-b939-3b1f555632aa
# ╠═cbaf7dad-5980-43a4-82c0-34bb946d9376
# ╠═970c1a3b-205c-421c-a564-62b8beba625e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
