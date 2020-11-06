### A Pluto.jl notebook ###
# v0.12.7

using Markdown
using InteractiveUtils

# ╔═╡ a4028042-2072-11eb-0aa6-77e29a756648
begin
	using TextAnalysis
	using Arpack
	using Clustering
end

# ╔═╡ 2abb6d30-2079-11eb-2d32-47713a476605
md"
# Clustering Textbook

We will try to cluster 144 sections of textbook _Matematika 1_ - a first semester calculus with some linear algebra and analytic geometry using spectral clustering on a `DocumentTermMatrix()`.

Package `TextAnalysis.jl` has some needed functionality. 
" 

# ╔═╡ 30498550-2071-11eb-0f6b-33372a37909e
#=
for i=1:144
	download("http://www.mathematics.digital/matematika1/predavanja/node$i.html","./files/node$i.html")
end
=#

# ╔═╡ 621ed340-2072-11eb-1742-3b4d0b042fb9
# Number of documents / nodes
n=144

# ╔═╡ 022dae30-2076-11eb-225e-233a63507457
a=Array{StringDocument{String}}(undef,n)

# ╔═╡ 2f61f232-2076-11eb-0606-351e586f37f1
for i=1:n
	f=open("files/node$i.html","r")
	fr=readlines(f)
	close(f)
	frjoin=join(fr,"  ")
	sdf=StringDocument(frjoin)
	remove_corrupt_utf8!(sdf)
	prepare!(sdf, strip_html_tags | strip_articles | strip_case | strip_pronouns | strip_stopwords | strip_numbers | strip_non_letters )
	a[i]=sdf
end

# ╔═╡ 0b67ef4e-2077-11eb-2716-5fed149d4e28
c=Corpus(a)

# ╔═╡ 168f01c0-2077-11eb-20ee-e381bbb1aa61
update_lexicon!(c)

# ╔═╡ 24e5a940-2077-11eb-3071-21fb31744527
update_inverse_index!(c)

# ╔═╡ 307d09b2-2077-11eb-2f26-2d00a7bae004
c

# ╔═╡ 3d92a5b0-2077-11eb-39ba-5df49551535e
c.documents

# ╔═╡ 4296d3ae-2077-11eb-0ecb-13e1aff1268c
c.lexicon

# ╔═╡ 4996d1b0-2077-11eb-2ac0-b7296b769be5
c.h

# ╔═╡ 4c8ce030-2077-11eb-0176-9de66dab00b3
c.inverse_index

# ╔═╡ 53a7b932-2077-11eb-2bf4-1feca7dfa435
c.total_terms

# ╔═╡ 58f17bb0-2077-11eb-0939-3dff13742e14
M = DocumentTermMatrix(c)

# ╔═╡ 6ed98bc2-2077-11eb-15fa-3f753c8ed3c7
D=dtm(M)

# ╔═╡ 78d31ab0-2077-11eb-0688-7d70444e8687
T = tf_idf(D)

# ╔═╡ 88b8c6a0-2077-11eb-3f48-2f453e82c0cb
md"
We are ready for spectral clustering. We try to find 16 clusters.
"

# ╔═╡ 8cb43730-2077-11eb-08ac-451458efac8d
S,rest=svds(T,nsv=6)

# ╔═╡ a70064ae-2077-11eb-2b51-3d76c0e7e933
size(S.U)

# ╔═╡ b34b84be-2077-11eb-3a72-91c2fc5739dd
outU=kmeans(Matrix(transpose(S.U)),6)

# ╔═╡ c7a41b32-2077-11eb-3702-c728d684fcfc
outU.assignments

# ╔═╡ Cell order:
# ╟─2abb6d30-2079-11eb-2d32-47713a476605
# ╠═a4028042-2072-11eb-0aa6-77e29a756648
# ╠═30498550-2071-11eb-0f6b-33372a37909e
# ╠═621ed340-2072-11eb-1742-3b4d0b042fb9
# ╠═022dae30-2076-11eb-225e-233a63507457
# ╠═2f61f232-2076-11eb-0606-351e586f37f1
# ╠═0b67ef4e-2077-11eb-2716-5fed149d4e28
# ╠═168f01c0-2077-11eb-20ee-e381bbb1aa61
# ╠═24e5a940-2077-11eb-3071-21fb31744527
# ╠═307d09b2-2077-11eb-2f26-2d00a7bae004
# ╠═3d92a5b0-2077-11eb-39ba-5df49551535e
# ╠═4296d3ae-2077-11eb-0ecb-13e1aff1268c
# ╠═4996d1b0-2077-11eb-2ac0-b7296b769be5
# ╠═4c8ce030-2077-11eb-0176-9de66dab00b3
# ╠═53a7b932-2077-11eb-2bf4-1feca7dfa435
# ╠═58f17bb0-2077-11eb-0939-3dff13742e14
# ╠═6ed98bc2-2077-11eb-15fa-3f753c8ed3c7
# ╠═78d31ab0-2077-11eb-0688-7d70444e8687
# ╟─88b8c6a0-2077-11eb-3f48-2f453e82c0cb
# ╠═8cb43730-2077-11eb-08ac-451458efac8d
# ╠═a70064ae-2077-11eb-2b51-3d76c0e7e933
# ╠═b34b84be-2077-11eb-3a72-91c2fc5739dd
# ╠═c7a41b32-2077-11eb-3702-c728d684fcfc
