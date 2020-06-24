---
layout: post
title:  "Select on RRR bitvector"
date:   2020-03-16 14:00:00 +0100
public: false
use_math: true
categories: informatics
---

Today, I will explain how to perform a select request on a bitvector compacted inside of a RRR datastructure.
Currently, I am working with Rayan Chikhi and two students on a library for bioinformatics that uses a bitvector to store information.
As our bitvector can be composed of $4^64$ bits, we have to represent in a compacted or compressed way.
This leads me to read the articles about the RRR datastrcture and [the blogpost from Alex Bow](https://alexbowe.com/rrr/) that clarify the way to perform a rank on it.
For our bioinformatics work, we will not use the RRR representation as the sd arrays are more compact but I dug to subject by curiosity to understand the technics behind the constant time select.
As the RRR original articles are hard to read, I hope that this post will help for understanding.

# Bitvector structure

A bitvector is a vector containing only 0's and 1's.
On a bitvector $b$, rank and select operations are defined as follow:

* $rank_{1}(i) = \sum_{j=0}^{i}{b_i}$ : Count the number of 1 from the begining of the vector to the ith position.
* $rank_{0}(i) = i - rank_{1}(i)$ : Count the number of 0 from the begining of the vector to the ith position.
* $select_{u}(i) = min(x \| rank_{u}(x) == i)$ : Return the position of the ith bit set to $u$.

For this post, I will assume that there is a majority of 0 in the bitvector and that we want to perform select operations on 1's.
Of course, everything is symetrical, so it's easy to centrer everything on 0's instead of 1's.

A bitvector has two main parameters: its size and the number of bit set to 1.
I will call $m$ the size of the bitvector and $n$ the number of 1.
By definition, a fully explicit bitvector occupy $m$ bits into the memory.
The RRR bitvector is spliced into small blocks of equal size $u$ (For optimality in the article $u = { {1}\over{2} } lg\,m$).
Each block represent a chunck of the universe that we call $U_i$ and there are $p$ blocks ($p = m / u$).

Basic blocs bit vector

Each of the $U_i$ are represented by two integers: the number of bit set to 1 $c_i$ and an offset identifyer ($o_i$).
The number of bit set is also called *class* of the chunck.
The are exactly ${u\choose c_i}$ different possible blocks containing $c_i$ bit set to 1.
Imagine that all these possible blocks are sorted lexicographically, then $o_i$ is the position of the block in that order.

Figure with enumeration and order of offsets for one class

Figure succint bitvector

For each block, storing $c_i$ take exactly $\lceil lg\,u \rceil$ bits regardless the block.
Storing $o_i$ take $\lceil lg{b\choose c_i} \rceil$ bits.
So, the amount of bit needed for a basic block depend on its class.
The values are computed in two different arrays: A for classes and B for offsets.
A prefix sum array is also computed on A where $psA[i] = \sum_{j=0}^{i-1}{A[j]}$.
Note that $psA[i] = rank(i \times u)$.

# Construction of the select datastructures

For rapid acess to select a few of them are precomputed and used as anchors for the search.
Instead of having select precomputed on regular indexes all along the vector, the RRR datastructure compute regular select values.
Regular select values means an array C where $C[0] = 0$ and $C[i] = select(i \times c) / u$, where c is a constant ($(lg\,p)^2$ in the article).
So, C contains the index of the blocks where the precomputed selects remain.

Then, for each segment of blocks between C[i] included and C[i+1] excpluded, we call this segment sparse if it contains at least $(lg\,p)^4$ bits.
In sparse segments, the block index of the select are exeplicitly stored.
To save space, the block indexes are relative to the begining of the segment.

Figure for sparse block

For dense segments, a tree datastructure is created with a branching factor of $\sqrt{lg\,p}$.
So, because there is at most $(lg\,p)^4$ bits, the height of the tree is bounded by a constant.
The leaves of the tree are the blocks of the segment.
The nodes are arrays containing the numbers of bit set to 1 in the subtrees.

Figure for dense block


# (Almost ?) O(1) select requests

Figure globale select

Let's consider that we want to perform a request select(x), where $0 <= x <= m$

1. First of all, we have to find the block segment where the select belongs to.
By computing $i= \lfloor x / \lfloor (lg\,p)^2 \rfloor \rfloor$ we have the index in C of the first block of the segment of interest.
C[i+1] is the index of the first block of the next segment.\\
\\
Figure of the segment with relative/real select i, i+1 positions.

2. Determine if the select value is present in the beginning of the block C[i+1].
Because the precomputed select are not neceseraly lying at the first position of a block, we have to verify the begining of the C[i+1] block.
If $rank(C[i+1] \times u) > (i+1) \times (lg\,p)^2$ (ie. $psA[C[i+1]]$), C[i+1] is the block where we have to look for the result.

3. If the select index is in the segment starting with the block C[i], we have to determine if it's sparse or dense regarding the length.
* If it's sparse, go to the begining of the sparse segment datastructure $\sigma$. The relative indexes of the blocks containing the ones are explicitely stored. So, first get the relative select value by doing $r = x - rank(C[i])$. Then get the relative index for the block of interest with $j_{rel} = \sigma[r]$. Finally the absolute index of the bloc is obtained by doing $j_{abs} = C[i] + j_{rel}$. The selected bit is inside of the block of index $j_{abs}$.
* If it's dense, get the tree corresponding to the segment $\sigma$.
From the root of the tree, perform a search into the array corresponding to the bit content of each subtree. Repeat the operation until the leaf that correspond to the relative index of the block containing the selected block. The same relative to absolute index search is performed thant the one in space segments.\\
\\
**Note**: for dense blocks, without further information, the search in an array of a node is not done in constant time. Using a partial sum array allow a dichotomic search for the right subtree. But the authors of the datastructure claim "This can be done in constant time using a lookup table" and then say that the whole structure can be strored in $O(\sqrt{lg\,p}\,lg\,lg\,p)$ but without further explanations. I am not sure how to do this step in constant time. So, if you, the reader, are able to explain me, please contact me !

4. Get the select value inside of the extracted block. Here there is no maggic, everything is precomputed in a huge lookup table. With the class and the offset of the block, we can directly jump to the precomputed row and then to the column corresponding to the selected bit inside of the block.


# Remarks

This blogpost is only due to my curisity on the theory but I think that this structure is not practical. Many improvements have be done since the first publication of this work. The purpose of the blogpost is only the clarification of the technics as the scientific papers about RRR are very complicated to read.

To understand and write this post, I was helped by Rayan Chikhi. I thank him a lot for the amount of time that he used for it.

