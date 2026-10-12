---
layout: post
title: "Demystifying Spherical Harmonics: An Optimization Approach"
date: 2026-10-19
---

Ever since I first heard of spherical harmonics as an approximation in computer graphics, I think when I was starting my career in 2011 or 2012, I have felt that the way they are explained is too complicated. Most explanations start with fairly abstract mathematical constructs, without giving much intuition or justification. In addition to making it harder to learn, I think the typical approach to SH can also lead to misapplication. So I've been wanting to write this post for a while, with a different way to look at SH that I hope is more straightforward and intuitive for some people.

Let's start with the engineering problem where you might want to use a spherical harmonic approximation.

We have a whole bunch of values that are oriented in different directions. In graphics, a common example of this is the incoming light (radiance) at a particular point in space. Here is a little interactive visualization where you can either see the lighting value along each surface normal (spherical projection), or as a reflection.

\<interactive diagram of a simple reflective sphere\>

For this one specific point in space, we can call these values a function f(d\_x, d\_y, d\_z), where the d values are the x,y,z components of our direction d (a unit vector). Remember f() is the \*actual\* directional values. This is expensive. Alternatively you can pre-render a high resolution version of the spherical function, possibly storing in a cubemap or octahedral map. But this requires a lot of memory for every point in space.

If we don't need very high resolution, we may want to store a fast low resolution version, we'll call f^. A 6-pixel (one for each face) cubemap or 8-value octahedral map works in many cases, but there are downsides which I won't go into here.

So instead, let's come up with just about the simplest function approximation we can think of, a linear model:

f^(d\_x, d\_y, d\_z) \= c0 \+ cx \* d \+ cy \* dy \+ cz \* dz

Where we will store our four coefficients: c0, cx, cy, cz. Note that f^ is called "F hat" and the hat is common notation for approximations in statistics.

It might be obvious that you can think of the c0 term as an overall intensity, and the cx, cy, cz as a directional intensity, for a particular direction.

\<interactive diagram showing simple spherical approximation\>

So how do we get these four coefficients? A quick idea might be to minimize the (sum/average/integrated) square error against the original function or a high-resolution version. In other words, a good old fashioned linear regression. In this code we'll sample some random directions and regress on those.

\[python code for linear regression\]

This gives us a pretty simple spherical approximation.

\[interactive diagram of sphere with spherical projection of original image, vs. linear model, with coefficients that can be manipulated or "reset"\]

Well, we pretty much just reinvented spherical harmonic approximations\! Yup, really. This is most likely something countless people have done without even realizing that what they did was equivalent to using SH.

Maybe it would be more accurate to say that this is a valid, though inefficient, implementation of a specific SH approximation. We would call this an l \= 1, or L1, spherical harmonics approximation, where "l" is the degree of the polynomial. You can also have a l \= 0 SH approximation, with a single constant like f^() \= c\_0, or higher l values that store coefficients for higher degrees of polynomials, like l \= 2 includes the terms x^2, xy, xz, y^2, etc. each with their own coefficient. Higher degree polynomials have a bit of extra complexity we'll talk about later, but it's still just a basic linear model, with a least squares fit\!

So anyway, spherical harmonic approximations are pretty much the blandest approximation model you could think of to represent a spherical function. If this is the first time you're learning about SH, you might be thinking... why do we even use a special name for a basic linear model?

You may notice that I've been saying SH approximations rather than just SH. A big part of the answer is the difference between theory (spherical harmonics) and application (SH approximations), as I alluded to at the beginning. Also part of the answer is due to the history of the application of SH to graphics, and due to the relationships with other things in math.

But SH approximations are also special compared to other linear models, because of two interesting properties:  
1\. Very fast and simple calculation of the coefficients (encoding)  
2\. Very fast convolutions of the spherical harmonic approximation

In graphics, performance is really important, so being able to do these things fast is a big part of why SH approximations are used at all.

As I'll show shortly, the first property is because the domain we are optimizing least squares over is the uniform distribution over the unit sphere, rather than any arbitrary distribution of x,y,z values.

ENCODING THE LINEAR MODEL COEFFICIENTS

The calculation of the coefficients is often called "encoding". In a machine learning context this would often be called regression or fitting or learning or maybe compressing or mayyybe training. When it comes to spherical harmonics, the fast \*method\* of encoding is specifically called "projection" for reasons that we'll discuss later.

Then evaluating the approximation function f^ I'll usually call reconstruction, but it could also be called decoding. In a ML context it might be called inference or estimation (or prediction, but that's not really appropriate here, since we're approximating something that's known, rather than predicting something unknown).

We're going to derive the fast encoding from scratch, for both the c0 term and the cx term. The cy and cz terms are basically the same as cx. The derivation involves a small amount of calculus, but if you want, you can skip the DERIVATION section, and go to ENCODING RULES.

DERIVATION

We'll start with c0. We want to find the c0 that minimizes least squares errors. This is written as:

argmin c0 integral (f^(d) \- f(d))^2

But rather than using spherical integrals, we're going to switch to spherical \*averages\* (means), because minimizing the integral of squared error is the same as minimizing the average squared error. Formally, the average is the same as the integral divided by the surface area of the sphere:

\< f(d) \> \= integral( f(d) ) / 4\*pi

Doing this means that we can ignore the surface area of the sphere, which makes the math and programming easier. The notation for a spherical average uses the angle brackets \< \>. I am not the first person to use this to simplify SH, and I think Graham Hazel did the best job fleshing this out, but I'll link to that article in "Bridging the Gap". Anyway, so we are trying to solve:

argmin c0 \< (f^(d) \- f(d))^2 \>

To find a minimum, we need to find the value of c0 such that the first derivative of this is zero. Typically to minimize you'd also worry about local minima or finding maxima, but as we'll see, there is only one solution here. This function is concave, it has no maxima and only one minima.

d( \< (f^(d) \- f(d))^2 \> ) / dc\_0 \= 0

The derivative of an average is the average of the derivative:

\<  d( (f^(d) \- f(d))^2 ) / dc\_0 \> \= 0

Typical power and chain rule:

\< 2\*(f^(d) \- f(d))\*d(f^(d) \- f(d)) / dc\_0 \> \= 0

Now let's focus on d(f^(d) \- f(d)) / dc\_0

d(f(d)) / dc\_0 is zero, because c\_0 is unrelated to the original term.

d(c0 \+ cx \* dx \+ cy \* dy \+ cz \* dz)/dc0 is simply 1\.

Plug this back a few lines up:

\< 2\*(f^(d) \- f(d))\*1 \> \= 0

Divide by 2, and separate all the average terms, like this:

\< c0 \> \+ \< cx\*dx \> \+ \< cy \* dy \> \+ \< cz \* dz\> \= \< f(d) \>

Now the average \<c0\> is of course c0, it's a constant. \< cx\*dx \> is 0, because it's symetrically distributed over a sphere (for every dx value, there is also a \-dx value), and same for the other three terms. So we get this wonderfully intuitive and easy to calculate value.

c0 \= \< f(d) \>

The constant is just the average value of our function over the whole sphere\!

Now let's do the same for cx...

argmin c\_x \< (f^(d) \- f(d))^2 \>

skipping a few steps, we have

\< (f^(d) \- f(d))\*d(f^(d) \- f(d)) / dc\_x \> \= \< f(d) \>

Focusing on the derivative term:

d(f(d)) / dc\_x is still zero, because c\_x is unrelated to the original term.

d(c0 \+ cx \* d\_x \+ cy \* dy \+ cz \* dz)/dc\_x is simply d\_x.

So now we have:

\< (f^(d) \- f(d)\*d\_x \> \= 0

or

\< f^(d)\*d\_x \> \= \< f(d)\*d\_x \>

Expanding f^(d) and separating the average terms, we have

\< c0\*d\_x \> \+ \< cx \* dx^2 \> \+ \< cy \* dy \* dx \> \+ \< cz \* dz \* dx \> \= \< f(d)\*d\_x 

Because of spherical symmetry, every average term is zero except the dx^2 one. So we have:

\< cx \* dx^2 \> \= \< f(d)\*d\_x \>

Because d is a unit vector, we know that:

dx^2 \+ dy^2 \+ dz^2 \= 1

So we can also say

\<dx^2\> \+ \<dy^2\> \+ \<dz^2\> \= 1

Since a sphere is symetrical for all the three axes, the averages are all the same: \<d\_x^2\> \= \<d\_y^2\> \= \<d\_z^2\>, and we can substitute that in, so \<d\_x^2\> \= 1/3. That gives us:

c\_x \= 3\*\<f(d)\*d\_x\>

Now you start to get a sense for why this is called "projection", because we're projecting the function onto the x directional component.

ENCODING RULES

Now we have our super fast and easy encoding of coefficients for our linear model:

c0 \= \< f(d) \>  
c\_x \= 3\*\<f(d)\*d\_x\>  
c\_y \= 3\*\<f(d)\*d\_y\>  
c\_z \= 3\*\<f(d)\*d\_z\>

And a little Python code can be written to estimate this numerically, like:

For a long time I found spherical harmonics a bit scary. I don't know about you, but derived this way, I find it hard to imagine a simpler spherical approximation\! It's a single line of code to evaluate, and only a few lines of code to encode, and the derivation is not really hard. And on top of that, it's doing something that many engineers are already familiar with, just a simple linear model with least squares, plus a trick to encode the coefficients in a faster way.

QUESTIONING OUR ASSUMPTIONS

Now that we've derived the fast encoding, let's reconsider some important assumptions we had to make. Typical applications of spherical harmonics in graphics violate at least one of these assumptions, sometimes two, often all three. Even when the assumptions are violated, spherical harmonics are still usually pretty good \*in practice\*, but sometimes the issues can be significant.

\- Our reconstruction function is f^(d\_x, d\_y, d\_z) \= c0 \+ cx \* d \+ cy \* dy \+ cz \* dz  
\- We want to minimize least squares error, not some other error.  
\- Our function or optimization domain is uniformly distributed around the unit sphere, in all directions.

The first assumption might seem comically obvious when encoding is derived in this way, but virtually all SH implementations violate it\! For one thing, there are improvements to spherical harmonics that use different forms of reconstruction, ie completely different models. Like the ZH3 hallucination, or the Geomerics reconstruction (this also includes a convolution).

But even the most vanilla implementations of SH almost always add a clamp to prevent negative lighting values:

f^(d\_x, d\_y, d\_z) \= max( c0 \+ cx \* d \+ cy \* dy \+ cz \* dz, 0 )

And in this case, "projection" is no longer the optimal least squares solution, and overall energy becomes biased \- the average lighting is too bright. That's not always necessarily bad, and I've even seen attempts to address it introduce worse problems. But it is possible to do better in specific cases. For example, SH approximations are often used to represent irradiance. In this case, the direction is a surface normal, and the value is the reflected light (equally in all directions, over the hemisphere) for a lambertian diffuse BRDF. 

For a single point light source, we can \*perfectly\* represent the irradiance using the clamped function, with zero error. But naive SH projection is overly soft, because it's minimizing error for the wrong reconstruction function.

\<interactive diagram showing naive projection vs. clamped, where the light direction can be dragged around\>

It's also generally the case that least squares error is rarely ideal for graphics applications. A relative error is usually better, and there are even better perceptual error metrics depending on what you know about your scene.

Finally, many people make the mistake of using standard or naive SH projection when your incoming light distribution is actually a hemisphere, for example when storing directional lightmaps. This can actually lead to meaningful bias, ending up too bright even for a completely flat surface where we should be matching non-directional lightmaps.

\<diagram showing bias of hemispherical lighting\>

In practice, these may or may not be significant problems in your application, but I have seen improvements to visual quality by addressing them.

CONVOLUTION

Let's say we have an L1 spherical harmonic approximation, encoded with our 4 coefficients c0, c\_x, c\_y, c\_z. Then we have a spherical convolution g(x), which could do something like "blur" or "sharpen" the approximation. 

I'll spoil the surprise before doing the derivation this time. It turns out that for a given convolution g(x), we can simply precompute 4 coefficients, and the convolution is as simple as multiplying each of the 4 coefficients. For example, to "convolve" from incoming radiance for every direction to outgoing irradiance for every surface normal, you just multiply the c\_x, c\_y, c\_z components by 2/3. That's it\!

Now we'll derive this. Let's say that a convolution takes in a spherical function or data, a direction, and outputs a value that is a weighted integral of the original spherical spherical data.

g( f, v ) \= integral h( v dot d )\*f(d) dd

NOTE: need to find the right definition of a convolution that is "constant" over the sphere. Something about the weighting function?

Now we can also express our encoding as an integral, for example

c0( f ) \= integral f (d) / 4 pi dd

We want to get the encoding of our convolution:

c0( g( f ) ) \= integral ( integral h( v dot d )\*f(d) dd ) / 4pi dv

We can flip the nesting order of the integrals (think of them just like sums):

c0( g( f ) ) \= integral ( integral h( v dot d )\*f(d) dv ) / 4pi dd

And pull f(d) out of the inner integral because it's unrelated to dv

c0( g( f ) ) \= integral ( integral h( v dot d ) dv \* f(d) / 4pi ) dd

Now integral h( v dot d ) dv is a constant value, we'll call it k, so we have

c0( g( f ) ) \= integral ( k \* f(d) ) / 4pi dd

Pull the k out of the integral, it's a constant:

c0( g( f ) ) \= k \* integral( f(d) ) / 4pi dd

Which is simply

c0( g( f ) ) \= k \* c0( f )

This might seem a bit mathy, but I think in plain English it's intuitive. c0 is just the average intensity of our function over a sphere. A convolution multiplies the average intensity up or down, but that's all it can do, and the multiplication amount is constant (it's the integral of the weight of the convolution kernel over the sphere).

We can do basically the same rearrangements for the other coefficients, and in the end, we get very simple scale factors. 

In other words, given any convolution g( f ), we can precalculate four k values, let’s say k\_0, k\_x, k\_y, k\_z. Then the convolved approximation g(f^)(d) and the approximation of the convolution g(f)^(d) are both just c\_0 \* k\_0 \+ c\_x \* k\_x \* d\_x \+ c\_y \* k\_y \* d\_y \+ c\_z \* k\_z \* d\_z.

Again, all of this is under the assumption that our integrals are uniformly oriented around a sphere. But if we're working on, say, hemispheres or arbitrary surfaces, things can get more complicated.

I'll add one other neat detail, which is that any "kernel transform" of a sphere (this could also include rotations, squishing, etc.) can be baked into a n x n matrix for n spherical harmonic coefficients, and transforming an encoding is just a matrix multiplication. Convolution is a special case where it’s a diagonal matrix.

HIGHER DEGREE POLYNOMIALS

Now let's say we want to use second-degree polynomials. We can add terms to our linear model to get this reconstruction function, which I'll call "Fast Reconstruction Form":

f^( d\_x, d\_y, d\_z ) \= c\_0 \+ c\_x\*d\_x \+ c\_y\*d\_y \+ c\_z\*d\_x \+ c\_xy\*d\_x\*\_dy \+ c\_xz\*d\_x\*d\_z \+ c\_yz\*d\_y\*d\_z \+ c\_x^2\*d\_x^2 \+ c\_y^2\*d\_y^2

You might notice that c\_z^2\*d\_z^2 is missing. This is because we know that this is a unit vector, so x^2 \+ y^2 \+ z^2 \= 1, and so z^2 \= 1 \- x^2 \- y^2, and those aspects are redundant with the coefficients c\_0, c\_x^2, and c\_y^2.

But if we want to find the encoding rules like we did for the L1 model, things get a little more complicated. We have two options for how to approach things:

Option 1: Keep reconstruction as cheap as possible (Fast Reconstruction Form), use the simple linear model given above.  
Option 2: Keep convolution as cheap as possible in many cases, and keep reconstruction more straightforward.

Let's start with Option 1 because it's more informative as to why it’s problematic. In practice I would typically use Option 2 in the same, real-time graphics. But Option 1 I’d use in any situation where we’re actually doing gradient descent anyway (say, as part of some other ML model).

Option 1: Fast Reconstruction Form

To understand why the encoding gets more complicated, we'll do some math using a simpler model of only one axis:

f^( d\_x ) \= c\_0 \+ c\_x \* d\_x \+ c\_x^2 \* d\_x \* d\_x

Let's now rederive the least-squares solution for the c\_0 term.

\[TODO: fill in here\]

As we can see, c\_0 is no longer a simple average and some of the terms seem to cross over and affect each other. I'm not going to derive everything, but if we use Fast Reconstruction Form

f^( d\_x, d\_y, d\_z ) \= c\_0 \+ c\_x\*d\_x \+ c\_y\*d\_y \+ c\_z\*d\_x \+ c\_xy\*d\_x\*\_dy \+ c\_xz\*d\_x\*d\_z \+ c\_yz\*d\_y\*d\_z \+ c\_x^2\*d\_x^2 \+ c\_y^2\*d\_y^2

Then we end up with these encoding rules:

\[TODO: fill in here\]

A convolution is still a matrix multiply of the coefficients, but isotropic convolutions are no longer diagonal matrices. So maybe the next form would be more accurate to call "Fast Isotropic Convolution Form", but I'm just going to call it Fast Convolution Form.

Option 2: Fast Convolution Form

Instead, we're going to change our reconstruction function to keep all the terms independent:

It may not be obvious why we've done, until we derive the reconstruction rules:

Now we can see that there is a comforting similarity between the encoding and reconstruction, where we are still just "projecting" the original function onto each term of the reconstruction approximation, and multipying by a constant.

BRIDGING THE GAP

I want to highlight how this approach differs from the typical approach to spherical harmonics, and start to introduce the real mathematical definitions.

To recap here's our reconstruction math

And here's our encoding math

By contrast, a rote implementation of SH does

With encoding as

If you plug the encoding values into the reconstruction, the terms cancel, and you end up with the same result. 

Now I want to specifically highlight these sub-terms, which we'll define like this:

f\_()

Coming back to the spherical harmonics vs. SH approximations, each one of these little functions is a spherical harmonic, and this is what spherical harmonics actually are: a set of (basis) functions on a sphere. An SH approximation is the linear combination of these basis functions, with encoded coefficients.

The usefulness of this basis function representation is that it's easier to generalize to higher degree polynomials.

[https\://gpfault.net/posts/sph.html](https://gpfault.net/posts/sph.html)

After reading that other article, I also recommend Graham Hazel's "Alternative Definition of Spherical Harmonics for Lighting". The form that Hazel uses is identical to my "fast for convolution" approach in the higher degree polynomial section, and I think it is probably the best general approach to take if you need higher than L1 spherical harmonics.

[https\://grahamhazel.com/blog/2017/12/18/alternative-definition-of-spherical-harmonics-for-lighting/](https://grahamhazel.com/blog/2017/12/18/alternative-definition-of-spherical-harmonics-for-lighting/)

AI disclaimer: all of the writing here is my own, not AI's, including the derivations. However the diagrams and formatting were done by AI, and I asked AI to check for mistakes.

\* It makes me sad when I see a definition of SH that starts off by saying "they are the solution to the Legendre polynomials over the sphere". It almost feels like a statement meant to gatekeep.