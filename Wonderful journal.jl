### A Pluto.jl notebook ###
# v0.12.3

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

# ╔═╡ c57da4c0-178f-11eb-2c26-874329cf054c
using Markdown

# ╔═╡ e0811450-178f-11eb-1669-e75839cbe2ff
using InteractiveUtils

# ╔═╡ f5bbf470-178f-11eb-3294-6f213346992e
begin
	using Pkg   
	Pkg.add.(["CSV", "DataFrames", "PlutoUI", "Shapefile", "ZipFile", "LsqFit", "Plots"])

	using CSV
	using DataFrames
	using PlutoUI
	using Shapefile
	using ZipFile
	using LsqFit
	using Plots
end

# ╔═╡ 4a5e36f0-1790-11eb-0da9-4370d41fb85d
using Dates

# ╔═╡ 5c6ad2d0-1796-11eb-08b0-3f55515186da


# ╔═╡ 74565400-1796-11eb-1333-110ec1c92e9d
md"""
## Download and load data
"""

# ╔═╡ 90703c9e-1796-11eb-3522-e3497ed8e5d7
url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv";


# ╔═╡ a23026d0-1796-11eb-17bc-5d715254365e
download(url, "covid_data.csv");

# ╔═╡ b53980ee-1796-11eb-02d3-dd4f532f5a16
md"Podaci koje koristimo su u CVS (*C*omma-*S*eparated *V*alues) fromatu.

Podatke možemo uzeti iz CSV datoteke koristeći `File` funkciju iz `CSV.jl` paketa"

# ╔═╡ 4d2c46e0-1797-11eb-0519-7ff7a6cb50e3
begin
	csv_data = CSV.File("covid_data.csv");   
	data = DataFrame(csv_data)   # it is common to use `df` as a variable name
end

# ╔═╡ 475fea10-179b-11eb-2579-4de7ca9e7f83
md"`DataFrame` je standardni način za spremanje **heterogenih podataka** u Juliji kao podaci koji se sastoje od stupaca s drugačijim tipovima vrijednosti."

# ╔═╡ f5116620-179b-11eb-2c59-9f163d0465ca
md"## Korištenje podataka"

# ╔═╡ 4da46ad0-179c-11eb-169b-5351388f02c6
md"""Promijenit ćemo nazive stupaca iz `DataFrame` u nešto kraće. To možemo napravit **in place** ili **out of place**, stvoriti novi `DataFrame`. Konvencija u Juliji za funkije koje mijenjaju argumente imaju isti naziv ali s **`!`** na kraju (zovu se "bang").

Da bi vidjeli samo prvih nekoliko redaka podataka možemo koristiti funkciju **`head`**.
"""

# ╔═╡ 0cc13de0-179c-11eb-046f-0f7ae9c527e4
begin
	data_2 = rename(data, 1 => "province", 2 => "country", 3 => "latitude", 4 => "longitude")   
	head(data_2)
end

# ╔═╡ 382c5aa0-179c-11eb-024b-a17dc15ad4be
begin
	rename!(data, 1 => "province", 2 => "country", 3 => "latitude", 4 => "longitude") 
	head(data)
end

# ╔═╡ 6fb48d20-179d-11eb-3c54-5d1e7348e30d
md"## Izvlačenje željenih podataka"

# ╔═╡ 8b82edd0-179d-11eb-3f43-8fc4fe499126
md"`DataFrame` možemo zamisliti kao matricu te koristimo slične sintakse.

Za primjer izvuć ćemo podatke iz 2. stupca (country)"

# ╔═╡ 1cf03cf0-179e-11eb-065a-8b4ebc1c2a1a
all_countries = data[:, "country"]

# ╔═╡ 28d56f40-179e-11eb-2eb8-497ab7d0f3bf
all_countries2 = data[:, :country]

# ╔═╡ 2da1c462-179e-11eb-0478-ab1b1ef57362
all_countries3 = data[:, 2]

# ╔═╡ 329577a0-179e-11eb-25c8-fd910c01f289
data[96:102, 2]

# ╔═╡ 01c811e0-179f-11eb-3b1c-230121df5dff
md"Dosta država je podijeljeno na `province` pa u `country` stupcu ima dosta ponavljanja pa da iz maknemo koristimo funkciju **`unique`**:"

# ╔═╡ 8c43140e-179e-11eb-3a4b-9f88e6f83ec2
countries = unique(all_countries)

# ╔═╡ 9563a820-179e-11eb-3ff5-33c4292ed4f8
@bind i Slider(1:length(countries), show_value=true)

# ╔═╡ ed2e7ad0-179e-11eb-3076-fbe340cf6cd2
md" $(Text(countries[i]))"

# ╔═╡ ab72f9d0-179f-11eb-08fb-39a181867aea
md"[Ovdje koristimo **string interpolation** pomoću **`$`** da stavimo u Markdown string.]"

# ╔═╡ ccaa64d0-179f-11eb-005f-3fe4b069a962
md"A i možemo koristiti **`Select`** da dobijemo padajući izbornik:"

# ╔═╡ eb0e165e-179f-11eb-0329-298e45979ce8
@bind country Select(countries)

# ╔═╡ 373a0e90-17a0-11eb-3e8a-d787679a3112
md"Da izvučemo podatke za pojedinu državu potrebno je da znamo ime države.

Možemo tražit ili možemo filtirat podatke da samo pogleda dio imena. Jedan od načina je sa **`array comprehension`**:"

# ╔═╡ e95d63ae-17a0-11eb-1264-27efa3ca4ebe
startswith("david", "d")

# ╔═╡ f5a48c22-17a0-11eb-39d9-0b70938fc47a
startswith("hello", "d")

# ╔═╡ fb1b5210-17a0-11eb-2396-3d938975bf4b
C_countries = [startswith(country, "C") for country in all_countries]

# ╔═╡ 0c83a660-17a1-11eb-23ba-c9fd2d221694
length(C_countries)

# ╔═╡ 15a1c970-17a1-11eb-3ff7-395aa3240e46
md"Kao što vidite funkcija vraća boolean. To ćemo iskorititi da indeksiramo `DataFrame`."

# ╔═╡ 552fedb0-17a1-11eb-2adb-bf36b0400aef
data[C_countries, :]

# ╔═╡ 72c21e70-17a1-11eb-00b6-1fa3a69e2ff2
md"Sad ćemo izvući podatke zasamo za  Kanadu tako da ponovno filtriramo na imenu koristeći funkciju **`filter`**.

Ovo je **higher-order function**: prvi argument je funkcija koja vraća bool (`true` `false`).
**`filter`** će vratiti redove koji ispunjavaju taj uvjet:"

# ╔═╡ 58ec6400-17a2-11eb-076d-65f8746aac4e
filter(x -> x.country == "Canada", data)

# ╔═╡ 785e2a30-17a2-11eb-0917-23d9ef1a3380
md"Ovdje koristimo **anonimnu funkciju** sa sintaksom `x -> ⋯` koja vraća što god je desno od strelice."

# ╔═╡ cac6341e-17a2-11eb-15d0-abd85466b68b
md"Da bi izvukli jedan red treba moramo imati indeks reda u `DataFrame`. **`findfirst`** funkcija traži prvi red koji ispunjava uvjet:"

# ╔═╡ 2b3ce00e-17a3-11eb-3c46-8bbfce1fe4e7
CA_row = findfirst(==("Croatia"), all_countries)

# ╔═╡ 3d398520-17a3-11eb-09f2-87635dd3f074
data[CA_row, :]

# ╔═╡ 43d3bc70-17a3-11eb-0e3c-11a5d8cef7ee
data[CA_row:CA_row, :]

# ╔═╡ 6e5897e2-17a3-11eb-0289-3dd8257f1e74
md"Izvlačimo podatke pomoću standarnog Julia **`Vector`**:"

# ╔═╡ 850caa80-17a3-11eb-0e85-951dc6272d11
CA_data = Vector(data[CA_row, 5:end])

# ╔═╡ 96626130-17a3-11eb-38a1-97015eb40c9f
scatter(CA_data, m=:o, alpha=0.5, ms=3, xlabel="day", ylabel="cumulative cases", leg=false)

# ╔═╡ ebb778f0-17a3-11eb-16ef-5f8b93ca9ae7
md"Koristimo **`scatter`** funkciju za jedan vektor tako da $x$ koordinate prime prirodne brojeve ($1$, $2$, itd.).$y$-os ispisuje samo kumulativne slučajeve do danog datuma."


# ╔═╡ 88a93cc0-17a4-11eb-3621-f97a215887f4
md"## Korištenje datuma"

# ╔═╡ 9608da10-17a4-11eb-3a78-95e12b1eaf16
md"Koristit ćemo datume umjesto brojeva dana iz snimljenih podataka."

# ╔═╡ c7d96282-17a4-11eb-26bc-37f71308f60d
column_names = names(data)

# ╔═╡ cb8e7b3e-17a4-11eb-11f6-e73984e2be10
date_strings = names(data)[5:end]

# ╔═╡ dceb84f0-17a4-11eb-2d2c-dd764dd7584e
md"Sad trebamo pretvoriti (**`parse`**) stringove u tip podataka koji može dani u `Dates.jl` paketu."

# ╔═╡ 30cce870-17a5-11eb-04eb-9f03942918ef
date_strings[1]

# ╔═╡ 3473a950-17a5-11eb-3e29-ef7d39b4a2a2
date_format = Dates.DateFormat("m/d/Y")

# ╔═╡ 3a071f00-17a5-11eb-206a-c7add444fb0d
parse(Date, date_strings[1], date_format);

# ╔═╡ 4f3d9250-17a5-11eb-3f10-b77eecb7a2ba
md"Pošto godina nije dobro napisana u originalnom obliku trebamo to ručno promijeniti."

# ╔═╡ 7d09f7a0-17a5-11eb-1853-8bee4040e75c
dates = parse.(Date, date_strings, date_format) .+ Year(2000)

# ╔═╡ a0f49bc0-17a5-11eb-0ae6-e9fa43b3d32c
begin
	plot(dates, CA_data, xrotation=45, leg=:topleft, 
	    label="Cro data", m=:o, ms=3, alpha=0.5)
	
	xlabel!("date")
	ylabel!("Kumulativni broj slučajeva u RH")
	title!("Kumulativni broj potvrđenih oboljelih slučajeva od COVID-19 virusa u RH")
end

# ╔═╡ 0048c0b0-17a6-11eb-3829-83b26cff7b06
md"## Istraživačka analiza podataka"

# ╔═╡ 0e855080-17a6-11eb-0ac6-85d689107bf8
md"
Sad ćemo pogledati broj oboljelih na dnevnoj bazi. Juslia ima **`diff`** funkciju koja računa razliku između uzastopnih argumenata vektora:
"

# ╔═╡ bf451e00-17a6-11eb-2fd8-a39032c7d106
begin
	daily_cases = diff(CA_data)
	plot(dates[2:end], daily_cases, m=:o, leg=false, xlabel="Datum", ylabel="Dnevni slučajevi u RH", alpha=0.5)   # use "o"-shaped markers
end

# ╔═╡ e6161930-17a6-11eb-3514-b974720b9c54
begin
	using Statistics
	running_mean = [mean(daily_cases[i-6:i]) for i in 7:length(daily_cases)]
end

# ╔═╡ f5697030-17a6-11eb-3db0-494e53521ff6
begin
	plot(daily_cases, label="raw daily cases")
	plot!(running_mean, m=:o, label="running weakly mean", leg=:topleft)
end

# ╔═╡ 042457c2-17a7-11eb-18bc-81744fb63e57
md"## Eksponencijalni rast"

# ╔═╡ 1c2312d0-17a7-11eb-03f4-3d30f386c083
begin
	plot(replace(daily_cases, 0 => NaN), 
		yscale=:log10, 
		leg=false, m=:o)
	
	xlabel!("Dan")
	ylabel!("Potvrđeni slučajevi u Hrvatskoj")
	title!("potvrđenih oboljelih slučajeva od COVID-19 virusa u RH")
end

# ╔═╡ 9274d680-17a7-11eb-1833-136369106530
xlims!(100, 300)

# ╔═╡ 99bb2c50-17a7-11eb-3cf9-07ac9b336700
exp_period = 100:279

# ╔═╡ a0e03e80-17a7-11eb-215e-670cedd60353
md" $(first(exp_period)) dan do $(last(exp_period)) dan :"

# ╔═╡ f14f85b0-17a7-11eb-3eee-1d37e4f1e743
dates[exp_period]

# ╔═╡ b458dfc0-17a8-11eb-2e2e-25a5a1a4fecc
md"## Uklapanje podataka"

# ╔═╡ d6872de0-17a8-11eb-28ba-33fbcd22dda8
md""" Julia `LsqFit.jl` paket ("least-squares fit") nam omogućuje da odredimo funkciju modela koja uzima vektor podataka i vektor parametara te najbolje prilagođavanje podacima."""

# ╔═╡ 4e3dffd0-17a9-11eb-172e-2d0502214e17
model(x, (c, α)) = c .* exp.(α .* x)

# ╔═╡ 53c5b5b0-17a9-11eb-2e11-45f50fb1dbbd
begin
	p0 = [0.1, 0.1]  # initial guess for parameters

	x_data = exp_period
	y_data = daily_cases[exp_period]
	
	fit = curve_fit(model, x_data, y_data, p0)
end;

# ╔═╡ 6e075be0-17a9-11eb-26d6-b1ef9775ad74
md"Koeficijenti najbolje prilagođenih podataka:"

# ╔═╡ 8b7c67b0-17a9-11eb-1a21-9f494f895381
parameters = coef(fit)

# ╔═╡ 9120581e-17a9-11eb-0bf8-fb5d8b98515a
begin
	plot(replace(daily_cases, 0 => NaN), 
		yscale=:log10, 
		leg=false, m=:o,
		xlims=(100, 300), alpha=0.5)
	
	line_range = 120:270
	plot!(line_range, model(line_range, parameters), lw=3, ls=:dash, alpha=0.7)
	
	xlabel!("Dan")
	ylabel!("Potvrđeni slučajevi u Hrvatskoj")
	title!("Potvrđenih oboljelih slučajeva od COVID-19 virusa u RH")
end

# ╔═╡ d04c8640-17a9-11eb-1485-d5099174b164
md"## Geografski podaci"

# ╔═╡ 315204b0-17aa-11eb-0d96-173a9656bfad
filter(x -> startswith(x.country, "A"), data)

# ╔═╡ 398c0950-17aa-11eb-15b3-4faac2ff1a96
province = data.province

# ╔═╡ 43326f30-17aa-11eb-32dc-dd70c304c09b
begin
	indices = ismissing.(province)
	province[indices] .= all_countries[indices]
end

# ╔═╡ 51b2bec0-17aa-11eb-0bb3-c582715fb465
begin 
	
	scatter(data.longitude, data.latitude, leg=false, alpha=0.5, ms=2)

	for i in 1:length(province)	
		annotate!(data.longitude[i], data.latitude[i], text(province[i], :center, 5, color=RGBA{Float64}(0.0,0.0,0.0,0.3)))
	end
	
	plot!(axis=false)
end

# ╔═╡ 5b43c7e0-17aa-11eb-2f79-6d19c29e9d95
data.latitude

# ╔═╡ 5fed79d0-17aa-11eb-2162-495a97f5065b
md"## Dodavanje mapa"

# ╔═╡ 7607e480-17aa-11eb-21f2-79ca3329fa42
md"Iscrtavanje zemalja koristimo pomoću  `Shapefile.jl` paketa."

# ╔═╡ 6c3a6e9e-17aa-11eb-1375-7bde1b0dba9e
begin
	zipfile = download("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/cultural/ne_110m_admin_0_countries.zip")

	r = ZipFile.Reader(zipfile);
	for f in r.files
	    println("Filename: $(f.name)")
		open(f.name, "w") do io
	    	write(io, read(f))
		end
    end
	close(r)
end

# ╔═╡ a5ed17b0-17aa-11eb-20a9-777a81339b88
shp_countries = Shapefile.shapes(Shapefile.Table("./ne_110m_admin_0_countries.shp"))

# ╔═╡ aba04a60-17aa-11eb-2a19-d10ad03403ef
plot!(shp_countries, alpha=0.2)

# ╔═╡ b55d939e-17aa-11eb-0889-19eb04c09e08
daily = max.(1, diff(Array(data[:, 5:end]), dims=2));

# ╔═╡ bce86c80-17aa-11eb-3d4d-2f5f6ec024c2
@bind day2 Slider(1:size(daily, 2), show_value=true)

# ╔═╡ 3ab74d70-17ab-11eb-0c66-1b275bbe181a
dates[day2]

# ╔═╡ 7af93150-17ab-11eb-3fef-776a3b761024
@bind day Clock(0.5)

# ╔═╡ 9c256070-17a5-11eb-3ca1-5ddd4b7a3a4c
dates[day]

# ╔═╡ 36d3d1b2-17ab-11eb-2ab9-59d56f858268
log10(maximum(daily[:, day]))

# ╔═╡ 3fb6e790-17ab-11eb-2e46-4f19c89f6ae8
world_plot = begin 
	plot(shp_countries, alpha=0.2)
	scatter!(data.longitude, data.latitude, leg=false, ms=2*log10.(daily[:, day]), alpha=0.7)
	xlabel!("zem. širina")
	ylabel!("zem. dužina")
	title!("Dnevni slučajevi po državi")
end

# ╔═╡ daec5560-17ab-11eb-03eb-07d72435158b
world_plot

# ╔═╡ Cell order:
# ╠═c57da4c0-178f-11eb-2c26-874329cf054c
# ╠═e0811450-178f-11eb-1669-e75839cbe2ff
# ╠═f5bbf470-178f-11eb-3294-6f213346992e
# ╠═4a5e36f0-1790-11eb-0da9-4370d41fb85d
# ╠═5c6ad2d0-1796-11eb-08b0-3f55515186da
# ╠═daec5560-17ab-11eb-03eb-07d72435158b
# ╟─74565400-1796-11eb-1333-110ec1c92e9d
# ╠═90703c9e-1796-11eb-3522-e3497ed8e5d7
# ╠═a23026d0-1796-11eb-17bc-5d715254365e
# ╟─b53980ee-1796-11eb-02d3-dd4f532f5a16
# ╠═4d2c46e0-1797-11eb-0519-7ff7a6cb50e3
# ╟─475fea10-179b-11eb-2579-4de7ca9e7f83
# ╟─f5116620-179b-11eb-2c59-9f163d0465ca
# ╟─4da46ad0-179c-11eb-169b-5351388f02c6
# ╠═0cc13de0-179c-11eb-046f-0f7ae9c527e4
# ╠═382c5aa0-179c-11eb-024b-a17dc15ad4be
# ╟─6fb48d20-179d-11eb-3c54-5d1e7348e30d
# ╟─8b82edd0-179d-11eb-3f43-8fc4fe499126
# ╠═1cf03cf0-179e-11eb-065a-8b4ebc1c2a1a
# ╠═28d56f40-179e-11eb-2eb8-497ab7d0f3bf
# ╠═2da1c462-179e-11eb-0478-ab1b1ef57362
# ╠═329577a0-179e-11eb-25c8-fd910c01f289
# ╟─01c811e0-179f-11eb-3b1c-230121df5dff
# ╠═8c43140e-179e-11eb-3a4b-9f88e6f83ec2
# ╠═9563a820-179e-11eb-3ff5-33c4292ed4f8
# ╠═ed2e7ad0-179e-11eb-3076-fbe340cf6cd2
# ╟─ab72f9d0-179f-11eb-08fb-39a181867aea
# ╟─ccaa64d0-179f-11eb-005f-3fe4b069a962
# ╠═eb0e165e-179f-11eb-0329-298e45979ce8
# ╟─373a0e90-17a0-11eb-3e8a-d787679a3112
# ╠═e95d63ae-17a0-11eb-1264-27efa3ca4ebe
# ╠═f5a48c22-17a0-11eb-39d9-0b70938fc47a
# ╠═fb1b5210-17a0-11eb-2396-3d938975bf4b
# ╠═0c83a660-17a1-11eb-23ba-c9fd2d221694
# ╟─15a1c970-17a1-11eb-3ff7-395aa3240e46
# ╠═552fedb0-17a1-11eb-2adb-bf36b0400aef
# ╟─72c21e70-17a1-11eb-00b6-1fa3a69e2ff2
# ╠═58ec6400-17a2-11eb-076d-65f8746aac4e
# ╟─785e2a30-17a2-11eb-0917-23d9ef1a3380
# ╟─cac6341e-17a2-11eb-15d0-abd85466b68b
# ╠═2b3ce00e-17a3-11eb-3c46-8bbfce1fe4e7
# ╠═3d398520-17a3-11eb-09f2-87635dd3f074
# ╠═43d3bc70-17a3-11eb-0e3c-11a5d8cef7ee
# ╟─6e5897e2-17a3-11eb-0289-3dd8257f1e74
# ╠═850caa80-17a3-11eb-0e85-951dc6272d11
# ╠═96626130-17a3-11eb-38a1-97015eb40c9f
# ╟─ebb778f0-17a3-11eb-16ef-5f8b93ca9ae7
# ╟─88a93cc0-17a4-11eb-3621-f97a215887f4
# ╟─9608da10-17a4-11eb-3a78-95e12b1eaf16
# ╠═c7d96282-17a4-11eb-26bc-37f71308f60d
# ╠═cb8e7b3e-17a4-11eb-11f6-e73984e2be10
# ╟─dceb84f0-17a4-11eb-2d2c-dd764dd7584e
# ╠═30cce870-17a5-11eb-04eb-9f03942918ef
# ╠═3473a950-17a5-11eb-3e29-ef7d39b4a2a2
# ╠═3a071f00-17a5-11eb-206a-c7add444fb0d
# ╟─4f3d9250-17a5-11eb-3f10-b77eecb7a2ba
# ╠═7d09f7a0-17a5-11eb-1853-8bee4040e75c
# ╠═9c256070-17a5-11eb-3ca1-5ddd4b7a3a4c
# ╠═a0f49bc0-17a5-11eb-0ae6-e9fa43b3d32c
# ╟─0048c0b0-17a6-11eb-3829-83b26cff7b06
# ╟─0e855080-17a6-11eb-0ac6-85d689107bf8
# ╠═bf451e00-17a6-11eb-2fd8-a39032c7d106
# ╠═e6161930-17a6-11eb-3514-b974720b9c54
# ╠═f5697030-17a6-11eb-3db0-494e53521ff6
# ╟─042457c2-17a7-11eb-18bc-81744fb63e57
# ╠═1c2312d0-17a7-11eb-03f4-3d30f386c083
# ╠═9274d680-17a7-11eb-1833-136369106530
# ╠═99bb2c50-17a7-11eb-3cf9-07ac9b336700
# ╟─a0e03e80-17a7-11eb-215e-670cedd60353
# ╠═f14f85b0-17a7-11eb-3eee-1d37e4f1e743
# ╟─b458dfc0-17a8-11eb-2e2e-25a5a1a4fecc
# ╟─d6872de0-17a8-11eb-28ba-33fbcd22dda8
# ╠═4e3dffd0-17a9-11eb-172e-2d0502214e17
# ╠═53c5b5b0-17a9-11eb-2e11-45f50fb1dbbd
# ╟─6e075be0-17a9-11eb-26d6-b1ef9775ad74
# ╠═8b7c67b0-17a9-11eb-1a21-9f494f895381
# ╠═9120581e-17a9-11eb-0bf8-fb5d8b98515a
# ╟─d04c8640-17a9-11eb-1485-d5099174b164
# ╠═315204b0-17aa-11eb-0d96-173a9656bfad
# ╠═398c0950-17aa-11eb-15b3-4faac2ff1a96
# ╠═43326f30-17aa-11eb-32dc-dd70c304c09b
# ╠═51b2bec0-17aa-11eb-0bb3-c582715fb465
# ╠═5b43c7e0-17aa-11eb-2f79-6d19c29e9d95
# ╠═5fed79d0-17aa-11eb-2162-495a97f5065b
# ╟─7607e480-17aa-11eb-21f2-79ca3329fa42
# ╠═6c3a6e9e-17aa-11eb-1375-7bde1b0dba9e
# ╠═a5ed17b0-17aa-11eb-20a9-777a81339b88
# ╠═aba04a60-17aa-11eb-2a19-d10ad03403ef
# ╠═b55d939e-17aa-11eb-0889-19eb04c09e08
# ╠═bce86c80-17aa-11eb-3d4d-2f5f6ec024c2
# ╠═36d3d1b2-17ab-11eb-2ab9-59d56f858268
# ╠═3ab74d70-17ab-11eb-0c66-1b275bbe181a
# ╠═3fb6e790-17ab-11eb-2e46-4f19c89f6ae8
# ╠═7af93150-17ab-11eb-3fef-776a3b761024
