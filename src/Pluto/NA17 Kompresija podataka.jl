### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 1dd2098b-8d4d-4fa4-a794-49badf9790c2
begin
	using Images
	using LinearAlgebra
	using Plots
	using PlutoUI
end

# ╔═╡ c94e89d0-1049-4b0d-aa26-d5a1eda436b9
md"""
# Kompresija podataka

QR rastav s pivotiranjem stupaca možemo koristiti za __kompresiju (sažimanje) podataka__.

Dijagonalni elementi matrice $R$ padaju po apsolutnoj vrijednosti pa možemo odrezati djelove matrica $Q$ i $R$ za koje smatramo da nisu značajni.

Dat ćemo primjer kompresije slike.
"""

# ╔═╡ 3d91d5e7-e134-464a-a0a4-032a28254715
img=load("../files/P8040001a.jpg")

# ╔═╡ 464d999f-f530-4cc0-b9c8-fa094d8a138c
# Opis podataka
typeof(img)

# ╔═╡ 86818c56-6d28-463a-aa87-3c0a715ef8ce
img[1,1]

# ╔═╡ ce428ce6-9bc5-4ce2-bb17-3bb8173fc937
show(img[1,1])

# ╔═╡ 1dfeac3f-51df-4fcf-919c-1675ae331189
# Razdvojimo sliku na R, G i B komponente
channels=channelview(img)

# ╔═╡ c3ae4d9d-aa93-4ad8-af1d-fad775d9c9bd
begin
	Red=channels[1,:,:]
	Green=channels[2,:,:]
	Blue=channels[3,:,:]
end

# ╔═╡ c1cbc1a9-0a1e-4ac2-90fb-e9ce3c20be62
colorview(Gray,Red)

# ╔═╡ c0cd1d7e-1fb9-4ae8-aa8d-794c36b85837
begin
	# Izračunajmo QR rastav s pivotiranjem matrice svakog kanala
	R=qr(Red,Val(true))
	G=qr(Green,Val(true))
	B=qr(Blue,Val(true));
end

# ╔═╡ 824b17f0-348b-11eb-30e2-9b12415490a3
R.R

# ╔═╡ ec0c9c42-348b-11eb-204c-8578725eba1d
R.R[570:576,570:590]

# ╔═╡ 32420d34-f066-4fdc-a488-15c4482640d8
norm(R.Q*R.R[:,invperm(R.p)]-float(Red))

# ╔═╡ c16bd7ca-a403-4f1b-adba-f0e9a028b14b
# Nacrtajmo dijagonalne elemente
scatter(1:length(diag(R.R,1)),abs.(diag(R.R)),
    title="Diagonal elements of matrix R",legend=false)

# ╔═╡ d7fd21c0-25ae-11eb-217c-51b165642c9b
@bind k Slider(10:10:200)

# ╔═╡ feebc2d6-e495-4145-8e90-3b0d7d073d04
begin
	# Izračunajmo komprimirane matrice za svaki kanal, RedC, GreenC i BlueC
	# Funkcija Matrix() je nužna radi bržeg generiranja matrice Q
	RedC=Matrix(R.Q)[:,1:k]*R.R[1:k,invperm(R.p)]
	GreenC=Matrix(G.Q)[:,1:k]*G.R[1:k,invperm(G.p)]
	BlueC=Matrix(B.Q)[:,1:k]*B.R[1:k,invperm(B.p)]
end

# ╔═╡ 955201fd-cb90-4c0f-bf50-ac1ca1850d26
# Nacrtajmo komprimiranu sliku
colorview(RGB, RedC, GreenC, BlueC)

# ╔═╡ 5bb5fe50-2466-11eb-175c-e1ff94b7840c
k

# ╔═╡ b1c18ee7-4b8f-4b44-a960-d759040a29cd
k, norm(Red-RedC)/norm(Red)

# ╔═╡ Cell order:
# ╟─c94e89d0-1049-4b0d-aa26-d5a1eda436b9
# ╠═1dd2098b-8d4d-4fa4-a794-49badf9790c2
# ╠═3d91d5e7-e134-464a-a0a4-032a28254715
# ╠═464d999f-f530-4cc0-b9c8-fa094d8a138c
# ╠═86818c56-6d28-463a-aa87-3c0a715ef8ce
# ╠═ce428ce6-9bc5-4ce2-bb17-3bb8173fc937
# ╠═1dfeac3f-51df-4fcf-919c-1675ae331189
# ╠═c3ae4d9d-aa93-4ad8-af1d-fad775d9c9bd
# ╠═c1cbc1a9-0a1e-4ac2-90fb-e9ce3c20be62
# ╠═c0cd1d7e-1fb9-4ae8-aa8d-794c36b85837
# ╠═824b17f0-348b-11eb-30e2-9b12415490a3
# ╠═ec0c9c42-348b-11eb-204c-8578725eba1d
# ╠═32420d34-f066-4fdc-a488-15c4482640d8
# ╠═c16bd7ca-a403-4f1b-adba-f0e9a028b14b
# ╠═feebc2d6-e495-4145-8e90-3b0d7d073d04
# ╠═5bb5fe50-2466-11eb-175c-e1ff94b7840c
# ╠═955201fd-cb90-4c0f-bf50-ac1ca1850d26
# ╠═d7fd21c0-25ae-11eb-217c-51b165642c9b
# ╠═b1c18ee7-4b8f-4b44-a960-d759040a29cd
