
using Polynomials
using SpecialMatrices
using Winston

n=6
a=0
b=pi
x=collect(linspace(a,b,n))
y=sin(x)

A=Vandermonde(x)

c=full(A)\y

p=Poly(c)

plot(x,y,"r*")

xx=linspace(a,b,100)
pS=polyval(p,xx)
sinus=sin(xx)
plot(x,y,"r*",xx,pp,xx,sinus,"b")

# maksimalne apsolutna i relativna pogreška
norm(pp[2:end-1]-sinus[2:end-1],Inf), norm((pp[2:end-1]-sinus[2:end-1])./sinus[2:end-1],Inf)

T(n,x)=cos(n*acos(x))

x1=linspace(-1,1,100)

y1=T(10,x1)

plot(x1,y1)

xn=map(Float64,[cos((2*k-1)*pi/(2*10)) for k=1:10])

yn=T(10,xn)

plot(x1,y1,xn,yn,"r*")

# Odaberimo za interpolaciju sinusa nultočke polinoma T(n,x)

xc=(a+b)/2+(b-a)/2*map(Float64,[cos((2*k-1)*pi/(2*n)) for k=1:n])

yc=sin(xc)

Ac=Vandermonde(xc)
cc=full(Ac)\yc
pc=Poly(cc)

xx=linspace(a,b,100)
pC=polyval(pc,xx)
sinus=sin(xx)
plot(x,y,"r*",xx,pp,xx,sinus,"b")

# maksimalne apsolutna i relativna pogreška
norm(pC[2:end-1]-sinus[2:end-1],Inf), norm((pC[2:end-1]-sinus[2:end-1])./sinus[2:end-1],Inf)

p1=plot(pp-sinus) # ravnomjerno rasporedjene točke
p2=plot(pC-sinus) # Čebiševljeve točke
display(p1),display(p2)


