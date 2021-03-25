### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ c162d7c0-0efc-11eb-16ad-3d5ff641f6e0
using LinearAlgebra

# ╔═╡ 70942b12-0bd1-44e6-8e04-862d4948e8a0
md"""
# Pozitivno definitne matrice i rastav Choleskog


Matrica $A$ je _pozitivno definitna_ ako je simetrična, $A^T=A$, i ako su sve njene svojstvene vrijednosti pozitivne.

Pozitivno definitnu matricu možemo rastaviti (bez pivotiranja) kao

$$
A=L L^T$$

pri čemu je $L$ donje trokutasta matrica s pozitivnim dijagonalnim elementima. Taj rastav se zove
_rastav Choleskog_ (vidi [Numerička matematika, poglavlje 3.6](http://www.mathos.unios.hr/pim/Materijali/Num.pdf)).

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
## Blok matrice
"

# ╔═╡ Cell order:
# ╟─70942b12-0bd1-44e6-8e04-862d4948e8a0
# ╠═a808c71a-97bb-4106-89fd-dcc0e902d2cb
# ╠═ea4f0c1f-55e7-4423-a8df-b9de1889e511
# ╠═26ab5250-1377-11eb-2824-69646a3f5444
# ╠═3c7d6960-1377-11eb-25b9-177b4c932ce6
# ╠═c162d7c0-0efc-11eb-16ad-3d5ff641f6e0
# ╠═dfea6e0a-a718-49ef-8e60-bd0ad6e204db
# ╠═396ba935-a1fc-4978-a64e-990bb1d78f44
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
