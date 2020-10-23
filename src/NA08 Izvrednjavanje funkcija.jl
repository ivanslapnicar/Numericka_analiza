### A Pluto.jl notebook ###
# v0.12.4

using Markdown
using InteractiveUtils

# ╔═╡ 08b094da-0680-4727-bdc4-83f71de41b32
using Polynomials

# ╔═╡ 5c6e1617-6682-43ec-bec5-8245f1a55d34
begin
	using Random
	Random.seed!(123)
end

# ╔═╡ 70d60c9a-ff90-4003-a2e9-b3f0b5c88fb2
md"""
# Izvrednjavanje funkcija


Računalo može izvoditi samo četiri osnovne operacije, `+`, `-`, `*` i `/` pa se sve ostale funkcije 
računaju pomoću polinoma (npr. Taylorova formula uz ocjenu ostatka ili bolje formule).

Neka je zadan polinom __stupnja__ $n$:

$$
p_n(x)=a_0+a_1 x+a_2x^2+a_3 x^3+\cdots + a_{n-1}x^{n-1}+a_n x^n,\quad a_n\neq 0.$$

## Brzina

Direktno računanje vrijednosti $p_n(x)$ treba $O(n^2)$ operacija.

Uz __pamćenje potencija__ imamo sljedeći algoritam:
"""

# ╔═╡ 5c6e0ac4-9b1e-498b-a2f5-994e2a596783
function mypolyval(p::Polynomial,x::Number)
    s=p[0]
    t=one(typeof(p[0]))
    for i=1:length(p)-1
        t*=x
        s+=p[i]*t
    end
    s
end 

# ╔═╡ d40a5115-1f38-4e4a-9bde-37408f90a5f5
p=Polynomial([1,2,3,4,5])

# ╔═╡ 58d8436a-740c-4801-aba3-92d01070a4e6
mypolyval(p,3)

# ╔═╡ 834e9612-5b50-4924-a158-1effe9046a74
mypolyval(p,π)

# ╔═╡ 2b402caa-7cbb-4402-8307-0a36af1410c5
md"""
Funkcija `mypolyval()` koristi $2n$ množenja i $n$ zbrajanja.
"""

# ╔═╡ 6da79660-153d-11eb-10f7-b9e079919759
pbig=Polynomial(rand(1000));

# ╔═╡ 923c54ff-bf4b-42df-96ae-96b22c728aa0
@time mypolyval(pbig,1.5)

# ╔═╡ e4ac8c0c-f8ec-4b95-9ca7-f8aa1b3dbdfe
@time pbig(1.5)

# ╔═╡ e3c1f285-fd03-47f5-b705-914e4cb35170
md"""
__Hornerova shema__ (Horner, 1819, Newton 1669) treba $n$ množenja i $n$ zbrajanja:
"""

# ╔═╡ b3718874-ea9f-4f5d-ad83-d22cb272b489
function myhorner(p::Polynomial,x::Number)
    s=p[end]
    for i=length(p)-2:-1:0
        # s*=x
        # s+=p[i]
        s=s*x+p[i]
    end
    s
end

# ╔═╡ 24d7653c-c02b-4239-ace6-b6b99ffea022
myhorner(p,3)

# ╔═╡ e3382c91-e26f-4636-b068-483dab80a3bc
@time myhorner(pbig,1.5)

# ╔═╡ cce27db4-32cd-4ff1-8531-880985527ff2
md"""
Hornerova shema je __optimalna__ u smislu da je općenito za izvrednjavanje polinoma $p_n(x)$ potrebno barem $n$ množenja. 

(Mogući su, naravno, posebni slučajevi, kao $x^{100}$.) 
"""

# ╔═╡ a7ab8371-c193-431d-9b1d-d3269a44459d
md"""
## Točnost

Neka je $\hat q$ vrijednost $p_n(x)$ izračunata u aritmetici s točnošću stroja $\varepsilon$. Tada vrijedi ocjena
(vidi [Accuracy and Stability of Numerical Algorithms, str. 105](https://books.google.hr/books?id=5tv3HdF-0N8C&printsec=frontcover&hl=hr#v=onepage&q&f=false)):

$$
\big|\, p_n(x)-\hat q\,\big| \leq \frac{2n\varepsilon}{1-2n\varepsilon} \sum_{i=0}^n |a_i||x|^i.$$
"""

# ╔═╡ b5ffd702-3c09-4b37-a58f-b376826fbbdb
begin
	p₁=Polynomial([1,2,3,4,5])
	myhorner(p₁,sqrt(2))
end

# ╔═╡ 4700a424-1f5e-4124-b851-c196d45dd8cc
pb₁=Polynomial(map(BigInt,[1,2,3,4,5]))

# ╔═╡ 167c82a2-014b-4ff6-b5fe-036096f00007
myhorner(pb₁,sqrt(map(BigFloat,2)))

# ╔═╡ ca410bf6-18b3-4d97-93a6-2fa6b42de453
myhorner(p₁,sqrt(200000))

# ╔═╡ fa022c78-1d6a-4490-bee9-feff0f7dbba1
myhorner(pb₁,sqrt(map(BigFloat,200000)))

# ╔═╡ 7b15e5da-8318-4b1c-9e72-07746cd8f2fe
begin
	r=[1,sqrt(2),3,4,5,6,sqrt(50)]
	p₂=fromroots(r)
end

# ╔═╡ 412c6ac6-5fc3-42a7-97a8-21d074a3b6b6
pb₂=fromroots(map(BigFloat,r))

# ╔═╡ 3f4613bb-42f4-4df8-90bc-3b134ce9910c
myhorner(p₂,sqrt(2)+0.1)

# ╔═╡ ce0209f9-8f95-42a6-b0ad-8ad4025ba74d
myhorner(pb₂,sqrt(map(BigFloat,2))+0.1)

# ╔═╡ 9bc8612e-84ce-40e0-bb76-6802b2957918
myhorner(p₂,-sqrt(10000))

# ╔═╡ 0d802709-35ac-424e-82e1-6d355dde5503
myhorner(pb₂,-sqrt(map(BigFloat,10000)))

# ╔═╡ 343bf9ea-27a3-44b8-b628-152941c5faf8


# ╔═╡ Cell order:
# ╟─70d60c9a-ff90-4003-a2e9-b3f0b5c88fb2
# ╠═08b094da-0680-4727-bdc4-83f71de41b32
# ╠═5c6e0ac4-9b1e-498b-a2f5-994e2a596783
# ╠═d40a5115-1f38-4e4a-9bde-37408f90a5f5
# ╠═58d8436a-740c-4801-aba3-92d01070a4e6
# ╠═834e9612-5b50-4924-a158-1effe9046a74
# ╟─2b402caa-7cbb-4402-8307-0a36af1410c5
# ╠═5c6e1617-6682-43ec-bec5-8245f1a55d34
# ╠═6da79660-153d-11eb-10f7-b9e079919759
# ╠═923c54ff-bf4b-42df-96ae-96b22c728aa0
# ╠═e4ac8c0c-f8ec-4b95-9ca7-f8aa1b3dbdfe
# ╟─e3c1f285-fd03-47f5-b705-914e4cb35170
# ╠═b3718874-ea9f-4f5d-ad83-d22cb272b489
# ╠═24d7653c-c02b-4239-ace6-b6b99ffea022
# ╠═e3382c91-e26f-4636-b068-483dab80a3bc
# ╟─cce27db4-32cd-4ff1-8531-880985527ff2
# ╟─a7ab8371-c193-431d-9b1d-d3269a44459d
# ╠═b5ffd702-3c09-4b37-a58f-b376826fbbdb
# ╠═4700a424-1f5e-4124-b851-c196d45dd8cc
# ╠═167c82a2-014b-4ff6-b5fe-036096f00007
# ╠═ca410bf6-18b3-4d97-93a6-2fa6b42de453
# ╠═fa022c78-1d6a-4490-bee9-feff0f7dbba1
# ╠═7b15e5da-8318-4b1c-9e72-07746cd8f2fe
# ╠═412c6ac6-5fc3-42a7-97a8-21d074a3b6b6
# ╠═3f4613bb-42f4-4df8-90bc-3b134ce9910c
# ╠═ce0209f9-8f95-42a6-b0ad-8ad4025ba74d
# ╠═9bc8612e-84ce-40e0-bb76-6802b2957918
# ╠═0d802709-35ac-424e-82e1-6d355dde5503
# ╠═343bf9ea-27a3-44b8-b628-152941c5faf8
