---
layout: post
title: "A Direct Construction of SZ Sequence Matrices"
subtitle: "How to generate the best sample sequences for 4D, 16D, 256D..."
date: 2026-11-12
latex: true
---

Most of this post is going to be extremely niche - there might be only a dozen people who are interested in the details. But I'm going to start off with what SZ Sequences are, why I think the 4D SZ sequence should be the *default* sample sequence when implementing something like an offline path-tracer (supplanting Sobol and PMJ02), and how to easily generate them. 

After that will be the fully detailed construction and proofs.

## (0,2)-Sequences

For two-dimensional sequences in a renderer, the state-of-the-art (and default choice for many renderers) are shuffled + scrambled base-2 (0,2)-sequences. This is either called the PMJ02 sequence or Shuffle-Scrambled Sobol' (0,2)-sequence. In fact, with shuffling and scrambling, these are *identical*, they are truly the same sequence. I mean that the probability density of getting any specific sequence is exactly the same.

I won't go into the full details of shuffling and scrambling, or the generation of these sequences. I will say I tried to cover the prerequisite knowledge well in the beginning sections of ["Stochastic Generation of (t,s) Sample Sequences"](https://www.seanet.com/~myandper/egsr21.pdf), and I also strongly recommend ["Practical Hash-based Owen Scrambling"](https://jcgt.org/published/0009/04/01/paper.pdf) because that is a more suitable implementation for a lot of renderers. But I will briefly mention the stratification properties:

For a base-2 (0,2)-sequence, this means that up to the $$2^k$$ points, they are stratified on every power-of-two 2D grid. For example, for the first 8 points, they are stratified on the 8x1 grid, the 2x4 grid, the 4x2 grid, and the 1x8 grid. This interactive diagram shows the stratification. 

<div id="sz-diagram-02-stratification" class="diagram"></div>
<script defer src="{{ site.baseurl }}/post_assets/2026-11-12-sz-matrix-construction/sz-diagrams.js"></script>

This combination of stratification, but randomness within each grid cell, yields pretty much the lowest error you can get for rendering (at least at powers-of-two).

## Higher Dimensionality

However, an open question has always been, what to do for higher dimensionalities? When you're generating a path-tracer, you need random numbers for each bounce, so your sampling really has 100s of dimensions. From what I know, production path-tracers use one of the two approaches:
1. Use the higher dimensional scrambled/stochastic Sobol' sequence, up to a limited number of dimensions (4 or 5 is typical, I think).
1. "Shuffle" two scrambled/stochastic (0,2)-sequences together to make a sequence with some weaker 4D stratification.

This gives 4D or 5D sequences, and then sequences "juxtaposed" with shuffling. Juxtaposing means that to make an 8D sequence, you simply combine two 4D sequences.

"Scrambling" is a technique to randomize the points within each grid cell, while "shuffling" is a technique to randomize the order of the points (the indices), and both can be done efficiently to preserve the stratification of the sequences. Shuffling is absolutely critical with juxtaposition, otherwise two 4D sequences will be highly correlated and may never converge, even for infinite samples.

## The Canonical SZ Sequence

SZ Sequences, defined in [SZ Sequences: Binary-Based (0, 2q)-Sequences](https://dl.acm.org/doi/10.1145/3763272) by Ahmed et al., are the best possible way to do Method #2.

Ahmed et al. talk about "SZ Sequences", but in this post, I'm going to describe "The Canonical SZ Sequence" or just "The SZ Sequence". It's my belief, though I haven't proven it, that there is only one shuffled-scrambled SZ Sequence.

So the SZ Sequence is actually an **infinite-dimensional** sequence, but for rendering you'll most likely want to use only the first four dimensions. The properties are like this:

Dimensions 0 and 1 form a base-2 (0,2)-sequence.
Dimensions 2 and 3 form a base-2 (0,2)-sequence.
Dimensions 0,1,2,3 form a base-4 (0,4)-sequence.

Not only are the first and second pairs of dimensions the ideal stratification, but all 4 dimensions together also have the ideal stratification (which includes all lower dimensional projections). Not only that, but they have a high-dimensional stratification property where $$16^k$$ points are always stratified on $$2^k \times 2^k \times 2^k \times 2^k$$.

Consider the first 16 points, which have the following stratifications:

* First two dimensions: $$[16 \times 1 \times 1 \times 1]$$, $$[8 \times 2 \times 1 \times 1]$$, $$[4 \times 4 \times 1 \times 1]$$, $$[2 \times 8 \times 1 \times 1]$$, $$[1 \times 16 \times 1 \times 1]$$
* Latter two dimensions: $$[1 \times 1 \times 16 \times 1]$$, $$[1 \times 1 \times 2 \times 8]$$, $$[1 \times 1 \times 4 \times 4]$$, $$[1 \times 1 \times 8 \times 2]$$, $$[1 \times 1 \times 1 \times 16]$$
* Other 2D projections: $$[4 \times 1 \times 4 \times 1]$$, $$[4 \times 1 \times 1 \times 4]$$, $$[1 \times 4 \times 4 \times 1]$$, $$[1 \times 4 \times 1 \times 4]$$
* High-dimensional: $$[2 \times 2 \times 2 \times 2]$$

It's not really possible to do better than this in stratification. And because it's a base-2 sequence, the generation, scrambling, and shuffling can all be done with fast binary operations.

I mentioned that the SZ Sequence is infinite-dimensional. Well it actually has this expansive nesting property. The way it works is this. *All* odd-even pairs (disjoint groups of two) of dimensions are a base-2 (0,2)-sequence, that is, any two dimensions $$kd, kd+1$$ for any integer $$k \geq 0$$. *All* disjoint groups of 4, i.e. any $$kd, kd+1, kd+2, kd+3$$ are a base-4 (0,4)-sequence. All disjoint groups of 16 are a base-16 (0,16)-sequence. For example, the first 256 points are 2D stratified on the 16x16 grid for every pair of dimensions. Not just for the first 16 dimensions, but the 16 after that, and 16 after that, and so on. All disjoint groups of 256 dimensions are a base-256 (0,256)-sequence. And so on, squaring the dimensionality each time, for an infinite number of dimensions.

The higher-dimensional forms of the SZ sequence probably wouldn't be useful in a path-tracer. But they can still be very useful for offline purposes, such as precalculating some integrals. I suspect there may be other fields they are useful for (like in finance maybe, or physics).

The main contribution of this article over Ahmed et al's paper is presenting a significantly simpler and faster algorithm to construct the "generator matrices" for the SZ Sequences. This algorithm means you could pick an arbitrary dimension, say the 2^50th dimension, and easily generate the points for that dimension.

## Using SZ Sequences in your renderer

The reason that SZ Sequences are so phenomenal for a renderer is that we are often combining two inherently 2D sampling spaces anyway. For example, a 2D lens position and a 2D pixel position. Or a 2D light sample, and a 2D BRDF sample. So combining two optimal 2D sequences into an optimal 4D sequence is pretty awesome for rendering! In fact, I happened to see that OpenQMC - which is maintained primarily by Framestore - [replaced the Sobol sequence with the SZ Sequence](https://github.com/AcademySoftwareFoundation/openqmc/commit/b6a204e1c542310017541985b2a17ac6e89535bf).

Now if you want to use randomized SZ Sequences in your renderer, there are two options:
1. Use the Stochastic Generation technique to precalculate sample sequences. I've added SZ Sequences to the Github repo (TODO: add link), so you can just use that!
1. Use Burley's Hash-Based Scrambling + Shuffling. I've made a new Github Repo specifically for this purpose.

In addition these, I've added randomized SZ Sequences to Samplin Safari.

## Backstory

(if you think I'm making this up, [I actually have a Colab notebook from 2021 with the exact matrix search](https://colab.research.google.com/drive/1HvpVKyB73nc41IKj5AEm-4lMaTBL-kwX?usp=sharing))

## Prerequisite Knowledge

## Matrix Construction Algorithm

With an understanding of the *digital method* for generating QMC sequences, let's walk through how to construct a generator matrix for the SZ sequence. First, we'll look at how to generate the single matrix for any arbitrary dimension $$d$$. Then we'll look at generating all $$s$$ dimensions at once, which is even more efficient since most of the work is shared across dimensions.

To generate the $$N \times N$$ binary generator matrix $$M$$ for dimension $$d$$:
1. Initialize $$M = I_N$$ (the $$N \times N$$ identity matrix).
1. If $$d = 0$$, return $$M$$.
1. Let $$j = \lfloor \log_2(d) \rfloor$$ be the index of the highest set bit of $$d$$.
1. Generate the basis block matrices $$E_0, E_1, \dots, E_j$$ (described below).
1. For $$i = 0 \dots j$$:
    1. If bit $$i$$ of $$d$$ is set (i.e. `(d >> i) & 1`), multiply $$M$$ on the right by the $$N \times N$$ block-Pascal matrix $$\mathcal{P}(E_i)$$:
       $$M \leftarrow M \cdot \mathcal{P}(E_i)$$
1. Return $$M$$.

To generate the $$j + 1$$ basis block matrices $$E_0, E_1, \dots, E_j$$:
1. Initialize $$E_0 = \begin{bmatrix} 1 \end{bmatrix}$$ and $$n = 1$$.
1. While $$n \le j$$:
    1. Set $$\beta = E_{n-1}$$ (which is an $$n \times n$$ matrix).
    1. For $$i = 0 \dots n - 1$$ (while $$n + i \le j$$):
        1. Let $$\tilde{E}_i$$ be $$E_i$$ repeated along the main diagonal to size $$n \times n$$:
           $$\tilde{E}_i = \begin{bmatrix} E_i & & 0 \\ & \ddots & \\ 0 & & E_i \end{bmatrix}_{n \times n}$$
        1. Construct the $$2n \times 2n$$ basis matrix $$E_{n+i}$$ as the $$2 \times 2$$ block matrix:
           $$E_{n+i} = \begin{bmatrix} 0_n & \tilde{E}_i \\ \beta \tilde{E}_i & \tilde{E}_i \end{bmatrix}$$
    1. Set $$n \leftarrow 2n$$.

In function generate_sz_matrix(d, size=32), I give a Python implementation of this.

Now, there are a bunch of different ways to make this more efficient. I won't go into them in detail except for two:
1. For all practical purposes, what we really care about are the *j* Pascal matrices. These can be precalculated to something like say, 32, which make it easy to generate any of the first 2^32 dimensions.
1. There is a mathematical property that the $$\mathcal{P}(A)\mathcal{P}(B) = \mathcal{P}(A) + \mathcal{P}(B)$$ if $$A$$ and $$B$$ are the same size (TODO: is this only true for finite fields?). Many of $$E_i$$ matrices are the same size, for example $$E_2$$ and $$E_3$$ are the same size, and $$E_4$$ through $$E_7$$ are the same size. In this case, rather than multiply the matrices, you add every element modulo 2. Because the matrices are stored as 32-bit integers, this is just doing one xor per row (or is it column?), so it's very fast.

With the Pascal matrices precalculated, it is super-duper fast to generate an arbitrary one of these matrices. So I've written a simple C++ method to do it. This can generate any matrix up to dimension $$2^32-1$$, and it averages about xx nanoseconds.