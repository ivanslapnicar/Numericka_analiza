# Numerička analiza

Bilježnice za predmet __Numerička analiza__ ili __Numerička matematika__ kao prvi (uvodni) jednosemestralni predmet kakav standardno slušaju studenti STEM studija.

Bilježnice se koriste na predmetu
  __[Numerička analiza](https://nastava.fesb.unist.hr/nastava/predmeti/8183)__ koji se predaje na diplomskom studiju Računarstva (250) na [FESBu](https://www.fesb.unist.hr/).

Bilježnice su pisane u programskom jeziku [Julia](https://julialang.org) koristeći paket [Pluto.jl](https://github.com/fonsp/Pluto.jl). Prilagođene su diretnoj i on-line nastavi.

# Pregledavanje bilježnica

Unutar svojeg preglednika, bilježnice možete pregledati na poveznici
[https://ivanslapnicar.github.io/Numericka_analiza](https://ivanslapnicar.github.io/Numericka_analiza/)

#  Izvršavanje bilježnica

## Binder

1. Idite na adresu  https://ivanslapnicar.github.io/NumericalMathematics/ i odaberite želejnu bilježnicu.
2. Pritisnite `Edit or run this notebook` i odaberite `binder`. Učitat će se svi poptrebni paketi i pokrenuti bilježnica (kroz nekoliko minuta).

## Računalo

1. Klonirajte cijeli repozitorij koristeći `git` naredbu:
```
git clone https://github.com/ivanslapnicar/Numericka_analiza.git
```
Ako niste upoznati s `git` alatom možete pogledati GitHubove [stranice za pomoć](https://help.github.com/articles/set-up-git/) ili direktno preuzeti bilježnice (repozitorij) kao zip datoteku. Možete koristiti i GitHub Dekstop.

2. Instalirajte [Julia-u](https://julialang.org/downloads/). U Julia terminalu izvedite naredbe
```
> using Pkg
> Pkg.add("Pluto")
```
Prethodne naredbe je potrebno izvršiti samo jednom.

3. Server Pluto bilježnica se pokreće pomoću naredbi
```
> using Pluto
> Pluto.run()
```

Sada možete izvršavati bilježnice koje se nalaze u direktoriju `Numericka_analiza/Pluto/`.
