---
layout: post
title:  "Select on RRR bitvector"
date:   2020-03-16 14:00:00 +0100
public: false
use_math: true
categories: informatics
---

Why this blogpost ?
Explain here bitvector + rank/select
Interest of bitvectors

The complete datastructure is composed of 2 main components.
The succint bitvector itself and complementary datastructures to support rank and select operations in constant time.

A few years ago, Alex Bow did a good work explaining [how to perform rank in constant time on RRR bitvectors](https://alexbowe.com/rrr/).
This article explain well the datastructure behind the constant time rank operation.
Here I will try to explain the select datastructure as clear as Alex Bow did the rank.

First I will go back to the bitvector succintness basis and then focus on the datastructure needed for the select.

# Bitvector structure

A bitvector is a vector containing only 0's and 1's.
For this post, I will assume that there is a majority of 0 in the bitvector and that we want to rank/select the 1 values.
Of course, everything is symetrical, so it's easy to centrer everything on 0's instead of 1's.

A bitvector has two main parameters: its size and the number of bit set to 1.
I will call *m* the size of the bitvector and *n* the number of 1.
By definition, a fully explicit bitvector occupy m bits into the memory.
Let see the RRR proposition to represent the bitvector.

First, the bitvector can be splited into equal size chucks of bits.
Let call these chucks *basic blocks* and $b$ their size.
For simplicity, let's assume that $b$ is a divisor of $m$ (if not, the last block is smaller than others and can be filled with 0s at the end).
So the vector is divided into $m/b$ basic blocks.

Basic blocs bit vector

The principle to occupy less bits than the complete vector is to replace each block by two numbers.
The first one is the number of bits set to 1 into the block (also called class $c_i$ of the block $i$) and the second one is an identifyier for the block (also called an offset $o_i$).
To explicit $o_i$, you have to know $c_i$.
When you have $c_i$ bits sets into the block i, then you can have ${b\choose c_i}$ distinct possible blocks.
With a total order on the complet set of possible blocks for the $c_i$ class, you can number them.
This number is the offset $o_i$ used to represent the block.

Figure with enumeration and order of offsets for one class

Figure succint bitvector

For each block, storing $c_i$ take exactly $lg\,b$ bits regardless the block.
Storing $o_i$ take $lg{b\choose c_i}$ bits.
So, the amount of bit needed for a basic block depend on its class.

For more information on how to access quickly the data and how to perform the rank operation in constant time, refer the nice [Alex Bow blogpost](https://alexbowe.com/rrr/).


# Select using constant time

