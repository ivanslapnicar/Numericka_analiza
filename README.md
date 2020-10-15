# Numerička analiza

Bilježnice za predmet _Numerička analiza_ ili _Numerička matematika_ kao prvi (uvodni) jednosemestralni predmet kakav standardno slušaju studenti STEM studija.

Bilježnice se koriste na predmetu
  _[Numerička analiza](https://nastava.fesb.unist.hr/nastava/predmeti/8183)_ koji se predaje na diplomskom studiju Računarstva (250) na [FESBu](https://www.fesb.unist.hr/).

Bilježnice su posebno prilagođene on-line nastavi. Bilježnice su pisane u programskom jeziku [Julia](https://julialang.org).

## Korištenje

Materijali su pisani kao [Jupyter](http://jupyter.org/) bilježnice (engl. _notebooks_) i/ili kao [Pluto](https://github.com/fonsp/Pluto.jl) bilježnice. Bilježnice možete koristiti na sljedeće načine:

### Korištenjem preglednika
Unutar svojeg preglednika, bilježnice možete pregledati na sljedećim poveznicama:
* [Jupyter notebook viewer](http://nbviewer.ipython.org/url/github.com/ivanslapnicar/Numericka_analiza/tree/master/src/)
* [Pluto](https://ivanslapnicar.github.io/Numericka_analiza/)

###  Lokalno preuzimanje i izvršavanje na vlastitom računalu
* Preuzmite bilježnice (repozitorij) korištenjem `git` naredbe:
```
git clone https://github.com/ivanslapnicar/Numericka_analiza.git
```
Ako niste upoznati s `git` alatom možete pogledati GitHubove [stranice za pomoć](https://help.github.com/articles/set-up-git/) ili direktno preuzeti bilježnice (repozitorij) kao zip datoteku.
* Instalirajte [Julia-u](https://julialang.org/downloads/). U Julia terminalu izvedite naredbe
```
> using Pkg
> Pkg.add("IJulia")
> Pkg.add("Pluto")
```
Prethodne naredbe je potrebno izvršiti samo jednom.
* Server Jupyter bilježnica se pokreće pomoću naredbi
```
> using IJulia
> notebook(detached=true)
```
a server Pluto bilježnica se pokreće pomoću naredbi
```
> using Pluto
> Pluto.run()
```

Sada možete izvršavati bilježnice koje se nalaze u direktoriju `Numericka_analiza/src`
