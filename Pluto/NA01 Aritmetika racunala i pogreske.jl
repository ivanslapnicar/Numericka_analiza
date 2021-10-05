### A Pluto.jl notebook ###
# v0.16.0

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

# ╔═╡ 62d67e4a-ae9f-43a2-bcc8-53ee5da9e268
using PlutoUI

# ╔═╡ 0425b895-a11b-49ee-b715-228384218624
TableOfContents(title="📚 Table of Contents", aside=true)

# ╔═╡ 76d37869-e20b-4211-8227-1f0616e3d8f2
md"""
# Floating Point Arithmetic and Errors

## Absolute and relative error

Let $\alpha$ approximate $a$. Then

$$err=|a-\alpha| \\  relerr=\frac{err}{|a|}=\frac{|a-\alpha|}{|a|}.$$
"""

# ╔═╡ 5c635357-8163-4954-949a-999dc48998f0
md"
α = $(@bind α Slider(0:0.1:6,show_value=true))
"

# ╔═╡ 9ccbb154-8618-4190-bfce-985ba66c8380
begin
	# Try α=a:0.01:2a
	a=5.0
	err=abs(a-α)
	relerr=err/abs(a)
	α, err, relerr
end

# ╔═╡ 91448849-f9f0-459b-87c0-b9fc5a386770
md"""
## Floating Point Arithmetic

Useful book on IEEE Floating Point standard:

M. Overton, Numerical Computing with IEEE Floating Point Arithmetic, SIAM Publications, Philadephia, 2001.

Useful article: 

[David Goldberg, What Every Computer Scientist Should Know About Floating-Point Arithmetic](https://docs.oracle.com/cd/E19957-01/806-3568/ncg_goldberg.html).

## Floating Point Number System

 $x$ is a __floating point number__ if it has the form

$$x = \pm d \cdot \beta^e \quad \beta \in \{ 2,10 \}$$

__Base__ 2 is for general purpose computers, __base__ 10 is for pocket calculators.

 $e$ is the __exponent__ and satisfies

$$e_{\min} \leq e \leq e_{\max},\quad e_{\min} < 0 < e_{\max}$$

We will assume that arithmetic is in base 2, but will usually give examples in base 10.

__Mantissa__ $d$ has the form

$$
\begin{aligned}
	d &= 0.d_1 \dots d_t = d_1 \beta^{-1} + d_2 \beta^{-2}
	+ \dots + d_t \beta^{-t}\\
d  &\in \{ 0,1\}\\
	d_1 &= 1 \qquad \mbox{ normalized }   \\
	d_1 &= 0 \qquad \mbox{ unnormalized }   \\
\end{aligned}$$

Standard form for floating point numbers is normalized except at the bottom of the exponent range.

During input and output numbers are converted from binary to decimal and back.

Computer arithmetic is standardized, by the IEEE 754 standard for binary arithmetic.  All but a few modern computers follow this standard.
"""

# ╔═╡ 1f8858f4-104d-4dd8-99de-b4bbf278e720
md"""
## Machine unit and machine precision

The set

$$
\{x \colon \lfloor \log_2 \: |x| \rfloor \in [e_{min},e_{max}] \}$$

is the set of real numbers that are in the normalized range of floating point numbers. $fl(x)$ is the floating point round of $x$.

__Machine unit__ is the maximum relative distance
between a real number in the floating point range and the nearest floating point number,

$$
\epsilon = \max_{\lfloor \log_2 
\:|x|\rfloor \in
[e_{\min},e_{\max}]} \frac{|x - fl(x)|}{|x|}  = 2^{-t}$$

__Machine precision__ is the relative distance between two neighbouring floating point numbers. For $\beta=2$ obviously $\epsilon=2\epsilon_M$.

Important examples include

__IEEE Standard Single Precision (Float32)__  $\beta = 2$, $t = 24$

$$
\begin{aligned}
\epsilon_M  &= 2^{-24} \approx	5.9605 \times 10^{-8}\\
\epsilon &=2^{-23} \approx 1.1920 \times 10^{-7} \\
e_{\min} &= - 126,\quad e_{\max} = 127.
\end{aligned}$$


__IEEE Standard Double Precision (Float 64)__ $\beta =2$,$t = 53$

$$
\begin{aligned}
\epsilon &= 2^{-53} \approx 1.1102 \times 10^{-16}\\
\epsilon &=2^{-52} \approx 2.2204 \times 10^{-16}\\
e_{\min} &= -1022,\quad e_{\max} = 1023.
\end{aligned}$$

Let us compute $\epsilon$ as the smallest positive floating point number such that $1+\epsilon\neq 1$.
"""

# ╔═╡ 030a0e5f-f9e1-4696-af85-b890eaa129d7
begin
	b₀=1.0
	a₀=2.0
	while (b₀+a₀)!=b₀
	    a₀/=2
	    println(a₀) # Output is in terminal
	end
	a₀
end

# ╔═╡ 82a6d40f-4c81-4ff1-8e02-3f427c413f61
md"""
The MATLAB command `eps` and the Julia function `eps()` return $\epsilon = 2.2204 \times 10^{-16}$.
"""

# ╔═╡ 7c28c479-912d-4a12-bc38-d15c6a3f0501
eps()

# ╔═╡ f5a5a27d-27bc-49b5-b245-c32a1f3fe13c
# What is this?
eps(64.0)

# ╔═╡ aa042b72-246d-432a-9a1c-d014da8ef957
md"""
Julia, in particular, has a type system where `Float64` type is a sub-type of `AbstractFloat`, which has four sub-types. 
In addition to types `Float64` and `Float32`, there is a type `Float16` which uses only two bytes of computer memory and type `BigFloat` which has a 256-bit mantissa.  
"""

# ╔═╡ dac635a5-7f41-478c-bede-5894410a0b6a
supertype(Float64)

# ╔═╡ 8c7b5693-2913-4cd3-9491-d8e57f101de8
subtypes(AbstractFloat)

# ╔═╡ 887cc5e4-ae65-4255-8a54-5425115e4618
for T in (Float16, Float32, Float64, BigFloat)
    println(eps(T))
end

# ╔═╡ b105ef44-00dc-4fe5-9732-499d14e51a51
2^(-10), 2^(-23), 2^(-52), 2^(-255)

# ╔═╡ 51356548-f58f-40a7-98b0-1b92ebeba3ed
md"""
## Basic Floating Point Operations

We begin with the four basic arithmetic operations, addition ($+$),subtraction ($-$),multiplication ($*$), and division ($/$). Suppose that $\odot$ is an operation such that

$$
\odot \in \{ + , - , *,/\}.$$

Then, in floating point arithmetic with machine unit $\epsilon_M$, it is reasonable to expect that for any two floating point numbers $x$ and $y$, we have

$$
fl(x\;op\;y) = (x \; op\; y)\;(1 + \xi),\quad
|\xi| \leq \epsilon_M.$$

For division, we assume $y \neq 0$. Any IEEE standard computer must follow this rule.  Rounding is one of two limitations that floating point arithmetic has that real arithmetic does not have. You can quickly conclude from the above rule that as long as all that we do is add numbers of the same sign, multiply, and divide, floating point results will almost always come very close to the corresponding real arithmetic results. The difficulty occurs if we either of $x$ or $y$ is rounded, they have different signs and we add or have the same signs and we subtract. 

That is, suppose we have

$$
\tilde{x}= x(1+\delta_x), \quad \tilde{y} = y(1+\delta_y),$$

where $x$ and $y$ are the exact real results of some computation and $\tilde{x}$ and $\tilde{y}$ are rounded floating point results with $|\delta_x| |\delta_y| \leq \delta$ for some small delta.  Suppose also that $x$ and $y$ have the same sign. Let

$$
z=x-y,\quad  \tilde{z} = fl(\tilde{x} -\tilde{y}).$$

Then, 

$$
\begin{aligned}
\tilde{z} &=(\tilde{x}-\tilde{y})(1+\xi)= x(1+\delta_x)(1+\xi) -y(1+\delta_y)(1+\xi) 
=x-y + \delta_z,
\end{aligned}$$

where $|\xi| \leq \epsilon$ and

$$
\delta_z = (x-y)\xi + (x\delta_x -y\delta_y)(1+\xi).$$

The best available bound on $|\delta_z|$ is

$$
\begin{aligned}
|\delta_z| &\leq |x-y||\xi| + (|x||\delta_x| + |y||\delta_y|)(1+|\xi|) \\
& \leq |x-y| \epsilon_M + (|x|+|y|)\,\delta\,(1+\epsilon_M).
\end{aligned}$$

Thus, the relative error in $z$ is 

$$
\begin{aligned}
\frac{|\tilde{z}-z|}{|z|}&=\frac{|\delta_z|}{|z|} 
\leq \epsilon_M + (1+\epsilon_M)\,\delta\,\frac{|x|+|y|}{|x-y|}\approx \delta \,\frac{|x|+|y|}{|x-y|}.
\end{aligned}$$

If $|x-y| << |x|+|y|$, the effect of the round in the subtraction is not important, but the error from previous computations on $x$ and $y$ can have a huge effect. The effect is called __propagation__. It can dramatically change the result of a compuation! We will see this issue with some examples later in this lecture.

Rounding is the first important limitation of floating point arithmetic.  A second limitation is the number range.
"""

# ╔═╡ 3170458d-931a-41a1-8715-b41de07aa3c6
md"""
## Number Ranges

Floating point arithmetic has a largest and smallest computer number. First, the largest one.

__Largest Computer Number__ $\Omega$

In base $2$, with a $t$ bit mantissa, the largest computer number is

$$
\Omega = (1 - 2^{-t}) \cdot 2^{e_{\max+1}}$$

When numbers exceed $\Omega$, they are stored as `Inf` ($\infty$) or `-Inf` ($-\infty$). We say than an _owerflow_ occured.


_IEEE Standard Single Precision_ (`Float32`)

$$
\quad \Omega = 3.4028\times 10^{38}$$

_IEEE Standard Double Precision_ (`Float64`)

$$
\Omega = 1.79777 \times 10^{308}$$

The MATLAB command `realmax` and the Julia function `floatmax()` show this number.
"""

# ╔═╡ 2ee4617f-70cf-4ecb-99b2-7167cc6b34d6
md"""
__Smallest Computer Number__ $\omega$

The definition of the smallest computer number is somewhat more complex.

The smallest computer number is given by

$$
\omega = 2^{1-t} 2^{e_{\min}}.$$

If a computation produces a number smaller in magnitude than $\omega$, it produces what is called an __underflow__, it is set to $0$ or $-0$.  If the programmer chooses, an underflow can result in an error, but in most computations, underflows are not harmful.


_IEEE Standard Single Precision_ (`Float32`):

$$
\omega = 2^{-23- 126} = 2^{-149} \approx  1.4013 \times 10^{-45}.$$

In MATLAB, this comes from the command `omega= eps('single')*realmin('single')`.


_IEEE Standard Double Precision_ (`Float64`):

$$
\omega= 2^{-1022-52} = 2^{-1074} \approx  4.9407 \times 10^{-324}$$

The appropriate MATLAB command to get this value is `omega = eps*realmin` and the equivalent Julia command is `floatmin()*eps()`.


__Important and Subtle Point__ 

Numbers at the bottom of the exponent
range are not normalized.

MATLAB function `realmin` yields

$$
\omega_{useful} \approx 2.2251 \times 10^{-308}.$$

Some people call this the smallest USEFUL floating point number since

$$
1/\omega_{useful} \leq \Omega$$

and $\omega_{useful}$ is normalized.

Smallest floating point number, $\omega$, has the form

$$
0.0 \cdots 01 \times 2^{e_{\min}} \quad \cdots\quad
\mbox{Gradual Underflow}$$

Before the IEEE standard most computers had the smallest floating point number as

$$
0.10 \cdots 0 \times 2^{e_{\min}} \qquad \cdots
\mbox{ normalized}$$

Earlier computers, (pre-1985) set numbers below this smallest 'useful' floating point number to zero. This change was one of the more controversial features of the IEEE standard.

__Example.__ $\beta = 10$, $-5 \leq e \leq 5$

$$
\begin{aligned}
x & = 0.1957 \times 10^{-5}   \\
y & = 0.1942 \times 10^{-5}
\end{aligned}$$

Compute  $fl(x - y)$. Whar happens?

$$
0.1957 \times 10^{-5}-0.1942 \times 10^{-5}  =0.0015 \times 10^{-5}$$

Pre-1985 philosophy was to set $fl(x - y)=0$.

Gradual Underflow stores $fl(x - y)=0.0015 \times 10^{-5}$, that is, Gradual Underflow guarantees that for any two floating point numbers $x$ and $y$

$$
fl(x - y) = 0 \mbox{ if and only if } x = y.$$
"""

# ╔═╡ b9dd627b-5b1a-4113-bfdf-af4dc1901827
for T in (Float16, Float32, Float64, BigFloat)
    println((floatmin(T),floatmax(T)))
end

# ╔═╡ 0f7ba6b2-aede-4f09-b7c5-adca1295195b
1/floatmin(Float32),floatmax(Float32)

# ╔═╡ b8647d13-c400-4324-85f7-8994a1f2322f
for T in (Float16, Float32, Float64)
    println((floatmin(T)*eps(T)))
end

# ╔═╡ d0b306c4-6cbc-48dc-90ba-8ef4498f8d73
md"""
##  Special Quantities  $0$, $-0$, `Inf`,`-Inf` i `NaN`

Zero has a sign:
"""

# ╔═╡ abe311d5-fbd1-40e2-8bce-2c341301deef
begin
	a₁=1.0
	b₁=0.0
	c₁=-b₁
	c₁,b₁==c₁
end

# ╔═╡ 88f1bf59-5eb4-4e2d-b28b-0d9b1004d5bd
a₁/b₁

# ╔═╡ 36c79f63-7a8c-46c9-afa5-95bbed8fd598
begin
	d₁=a₁/b₁
	e₁=a₁/c₁
	d₁==e₁, 1/d₁==1/e₁
end

# ╔═╡ 92490cef-70af-417e-8bc8-91b1be0635cc
b₁/c₁

# ╔═╡ afc32ff2-014a-43ee-8e6a-d35a57434622
md"""
`NaN` (Not a Number) can be generated by (c.f. Calculus 1):
"""

# ╔═╡ 03f1778e-c222-4221-a4e7-eabf1071d298
Inf+(-Inf),0*Inf, Inf/Inf, 0.0/0.0

# ╔═╡ 307a27e6-8128-42be-9ba8-cbdacb498ada
md"""
IEEE Arithmetic is a closed system:

 $\big\{$ floating point numbers,`Inf`,`-Inf`, `NaN`$\big\}$ 
$\stackrel{\odot}{\rightarrow}$ 
$\big\{$ floating point numbers,`Inf`,`-Inf`, `NaN` $\big\}$

no matter what the operation $\odot$ is.

Clever programmers take advantage of these features. However, in the coding assignments in this course, if you get
`NaN` or `Inf` or `-Inf`, you have probably made an error.
"""

# ╔═╡ 20fa79ac-93c3-4e3f-bf31-f9c0f0e79d4a
md"""
## Binary Representation
"""

# ╔═╡ 1f8f4a15-0eca-4759-bee6-843306e35754
bitstring(0)

# ╔═╡ 8ca8e4e7-6450-4aaf-b208-6c0aa95ccb12
bitstring(1)

# ╔═╡ f6f7081b-212a-4c7a-b91c-0d747f4b9fb2
bitstring(0.0)

# ╔═╡ d7160016-c046-4ab8-b349-27f7a209d853
bitstring(-0.0)

# ╔═╡ 10dcdca1-71fe-42eb-8f04-6122b3a03666
bitstring(1.0)

# ╔═╡ f85ec710-1f59-4826-a376-6672bf0552ae
bitstring(Float16(1.0))

# ╔═╡ 9b7e4278-2f1a-4fff-a6f7-2249f45aaccf
bitstring(2.0)

# ╔═╡ 037cddff-3663-4fcb-8575-ea7a78e490b8
md"""
__Problem.__ Explain the above binary representations. 
"""

# ╔═╡ 4e97e67b-8ef6-412d-8851-2b6635ff46e6
md"""
## Examples

### Using difference of squares

Compute

$$
f(x) = \sqrt{1 + x^2} - 1, \quad \mbox{$x$ is near zero}.$$

This formula in standard double precision yields $f(10^{-12}) = 0$.
"""

# ╔═╡ 88584373-9519-46fe-ae88-7199e76fb8f7
begin
	f(x)=sqrt(1+x^2)-1
	[(x,f(x)) for x in [1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12]]
end

# ╔═╡ 1215b318-8121-4adb-aa5d-139654607717
md"""
The difference-of-squares trick yields

$$
\begin{aligned}
f(x) & \equiv (\sqrt{1 + x^2} - 1) \left( \frac{\sqrt{1 + x^2} + 1}{\sqrt{1 + x^2} + 1}\right) \\
& = \frac{x^2}{\sqrt{1+x^2} + 1}\equiv f_1(x),
\end{aligned}$$

that is,  $f_1(10^{-12}) = 0.5 \cdot 10^{-24}$. This answer is as accurate as we can expect in standard double precision.
"""

# ╔═╡ b9180215-218a-48d8-a662-2a5664c2e480
begin
	f₁(x)=x^2/(1+sqrt(1+x^2))
    [(x,f₁(x)) for x in [1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12]]
end

# ╔═╡ 209cad48-c496-41b4-8093-6e2a55cc469f
x=1e-12

# ╔═╡ 905dd1d7-5f1b-43b3-a461-651cc992bd1c
# Using BigFloat
BigFloat(x)

# ╔═╡ c5d305a9-5a27-44d6-a610-50ab6d396f93
f(BigFloat(x))

# ╔═╡ 82ac21a4-384e-47e8-a964-95c8eb0889bd
md"""
### Quadratic equation

In exact arithmetic, the quadratic equation

$$ax^2 + bx+c=0$$

has roots

$$
\begin{aligned}
x_1&=\frac{-b-\mathop{\mathrm{sign}}(b)\sqrt{b^2-4ac}}{2a} \\
x_2&\equiv\frac{-b+\mathop{\mathrm{sign}}(b)\sqrt{b^2-4ac}}{2a}= \frac{-b+\mathop{\mathrm{sign}}(b)\sqrt{b^2-4ac}}{2a}\cdot \frac{-b-\mathop{\mathrm{sign}}(b)\sqrt{b^2-4ac}}{-b-\mathop{\mathrm{sign}}(b)\sqrt{b^2-4ac}}
\\ &= \frac{2c}{-b-\mathop{\mathrm{sign}}(b)\sqrt{b^2-4ac}}\equiv x_3.
\end{aligned}$$
"""

# ╔═╡ 139debf3-10ad-48da-9c2d-a13033486fbf
begin
	a₂=2.0
	b₂=123456789.0
	c₂=4.0
	
	x₁(a,b,c)=(-b-sqrt(b*b-4*a*c))/(2.0*a)
	x₂(a,b,c)=(-b+sqrt(b*b-4*a*c))/(2.0*a)
	x₂ₐ(a,b,c)=(2*c)/(-b-sqrt(b*b-4*a*c))
	x₁(a₂,b₂,c₂),x₂(a₂,b₂,c₂),x₂ₐ(a₂,b₂,c₂)
end

# ╔═╡ d429eede-d77d-4709-b635-c3993fce4a47
md"""
Let us check using `BigFloat`:
"""

# ╔═╡ ff5e3103-11cb-4438-aa19-c25af84fc3da
x₂ₐ(a₂,b₂,c₂)

# ╔═╡ ec755f39-5dba-4498-9760-07629ad00dc7
x₂(BigFloat(a₂),BigFloat(b₂),BigFloat(c₂))

# ╔═╡ c94ccf15-a0bd-4605-96f2-a27d0508fbeb
md"""
### Tangent and sine
"""

# ╔═╡ 65f44437-4d44-4d93-a9bc-814ab9b6ee02
begin
	x₃=1e-10
	tan(x₃)-sin(x₃)
end

# ╔═╡ bf8fd692-0e01-4368-a52d-91edb4fe77c0
md"""
However, the trigonometric identities give

$$
\begin{aligned}
\tan x - \sin x & = \tan x (1 - \cos x ) 
= \tan x (1-\cos x)\frac{1+\cos x}{1+\cos x}\\ & = \tan x \frac{1-\cos^2 x}{1+\cos x} \\
&= \tan x \sin^2 x \frac{1}{1+\cos x},
\end{aligned}$$

and Taylor formula gives

$$
\begin{aligned}
\tan x &= x + \frac{x^3}{3} + \frac{2x^5}{15} + O(x^7) \\
\sin x &= x -\frac{x^3}{6} + \frac{x^5}{120}+O(x^7) \\
\tan x - \sin x &= \frac{x^3}{2} + \frac{7x^5}{120} +O(x^7).
\end{aligned}$$

Both formulas give accurate resut:
"""

# ╔═╡ d58c4727-448d-47b2-a194-2002304a4708
tan(x₃)*sin(x₃)^2/(1+cos(x₃)), x₃^3/2+7*x₃^5/120

# ╔═╡ 4035ef81-4348-4721-9448-48c569a65c9d
md"""
### Absolute value of a complex number

To avoid underflow or overflow, instead of using the standard formula 

$$
|z|=|x+iy|=\sqrt{x^2+y^2}$$

we must use the following formulas (Explain!):

$$
M = \max \{ |x|,|y|\}, \quad m = \min \{ |x|,|y| \}, \quad r = \frac{m}{M}, \quad 
|z| = M \sqrt{1+r^2}.$$

These formulas are encoded in the function `abs()`.
"""

# ╔═╡ ac2caf5e-9874-44b9-ba5e-f0011940055f
z=2e+170+3e+175*im

# ╔═╡ b990e038-aff2-4d44-bb27-eb9de8c3977a
√(real(z)^2+imag(z)^2), abs(z)

# ╔═╡ c1acfea2-a60d-4433-a00b-6d3515274a18
begin
	z₁=2e-170+3e-175*im
	√(real(z₁)^2+imag(z₁)^2), abs(z₁)
end

# ╔═╡ b573a376-d60f-4f1d-b876-592dbcd47be4
md"""
__Problem.__ Compare the function [hypot](https://en.wikipedia.org/wiki/Hypot) and the  BLAS 1 function `dnrm2.f`.
"""

# ╔═╡ 03679ace-1368-4be3-988c-1f0dcaa1407e
real(z)^2, imag(z)^2

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.14"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[HypertextLiteral]]
git-tree-sha1 = "72053798e1be56026b81d4e2682dbe58922e5ec9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.0"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "a8709b968a1ea6abc2dc1967cb1db6ac9a00dfb6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.0.5"

[[PlutoUI]]
deps = ["Base64", "Dates", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "d1fb76655a95bf6ea4348d7197b22e889a4375f4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.14"

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

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
"""

# ╔═╡ Cell order:
# ╠═62d67e4a-ae9f-43a2-bcc8-53ee5da9e268
# ╠═0425b895-a11b-49ee-b715-228384218624
# ╟─76d37869-e20b-4211-8227-1f0616e3d8f2
# ╟─5c635357-8163-4954-949a-999dc48998f0
# ╠═9ccbb154-8618-4190-bfce-985ba66c8380
# ╟─91448849-f9f0-459b-87c0-b9fc5a386770
# ╟─1f8858f4-104d-4dd8-99de-b4bbf278e720
# ╠═030a0e5f-f9e1-4696-af85-b890eaa129d7
# ╟─82a6d40f-4c81-4ff1-8e02-3f427c413f61
# ╠═7c28c479-912d-4a12-bc38-d15c6a3f0501
# ╠═f5a5a27d-27bc-49b5-b245-c32a1f3fe13c
# ╟─aa042b72-246d-432a-9a1c-d014da8ef957
# ╠═dac635a5-7f41-478c-bede-5894410a0b6a
# ╠═8c7b5693-2913-4cd3-9491-d8e57f101de8
# ╠═887cc5e4-ae65-4255-8a54-5425115e4618
# ╠═b105ef44-00dc-4fe5-9732-499d14e51a51
# ╟─51356548-f58f-40a7-98b0-1b92ebeba3ed
# ╟─3170458d-931a-41a1-8715-b41de07aa3c6
# ╟─2ee4617f-70cf-4ecb-99b2-7167cc6b34d6
# ╠═b9dd627b-5b1a-4113-bfdf-af4dc1901827
# ╠═0f7ba6b2-aede-4f09-b7c5-adca1295195b
# ╠═b8647d13-c400-4324-85f7-8994a1f2322f
# ╟─d0b306c4-6cbc-48dc-90ba-8ef4498f8d73
# ╠═abe311d5-fbd1-40e2-8bce-2c341301deef
# ╠═88f1bf59-5eb4-4e2d-b28b-0d9b1004d5bd
# ╠═36c79f63-7a8c-46c9-afa5-95bbed8fd598
# ╠═92490cef-70af-417e-8bc8-91b1be0635cc
# ╟─afc32ff2-014a-43ee-8e6a-d35a57434622
# ╠═03f1778e-c222-4221-a4e7-eabf1071d298
# ╟─307a27e6-8128-42be-9ba8-cbdacb498ada
# ╟─20fa79ac-93c3-4e3f-bf31-f9c0f0e79d4a
# ╠═1f8f4a15-0eca-4759-bee6-843306e35754
# ╠═8ca8e4e7-6450-4aaf-b208-6c0aa95ccb12
# ╠═f6f7081b-212a-4c7a-b91c-0d747f4b9fb2
# ╠═d7160016-c046-4ab8-b349-27f7a209d853
# ╠═10dcdca1-71fe-42eb-8f04-6122b3a03666
# ╠═f85ec710-1f59-4826-a376-6672bf0552ae
# ╠═9b7e4278-2f1a-4fff-a6f7-2249f45aaccf
# ╟─037cddff-3663-4fcb-8575-ea7a78e490b8
# ╟─4e97e67b-8ef6-412d-8851-2b6635ff46e6
# ╠═88584373-9519-46fe-ae88-7199e76fb8f7
# ╟─1215b318-8121-4adb-aa5d-139654607717
# ╠═b9180215-218a-48d8-a662-2a5664c2e480
# ╠═209cad48-c496-41b4-8093-6e2a55cc469f
# ╠═905dd1d7-5f1b-43b3-a461-651cc992bd1c
# ╠═c5d305a9-5a27-44d6-a610-50ab6d396f93
# ╟─82ac21a4-384e-47e8-a964-95c8eb0889bd
# ╠═139debf3-10ad-48da-9c2d-a13033486fbf
# ╟─d429eede-d77d-4709-b635-c3993fce4a47
# ╠═ff5e3103-11cb-4438-aa19-c25af84fc3da
# ╠═ec755f39-5dba-4498-9760-07629ad00dc7
# ╟─c94ccf15-a0bd-4605-96f2-a27d0508fbeb
# ╠═65f44437-4d44-4d93-a9bc-814ab9b6ee02
# ╟─bf8fd692-0e01-4368-a52d-91edb4fe77c0
# ╠═d58c4727-448d-47b2-a194-2002304a4708
# ╟─4035ef81-4348-4721-9448-48c569a65c9d
# ╠═ac2caf5e-9874-44b9-ba5e-f0011940055f
# ╠═b990e038-aff2-4d44-bb27-eb9de8c3977a
# ╠═c1acfea2-a60d-4433-a00b-6d3515274a18
# ╟─b573a376-d60f-4f1d-b876-592dbcd47be4
# ╠═03679ace-1368-4be3-988c-1f0dcaa1407e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
