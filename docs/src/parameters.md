# Design Parameters

## Overview

Objects of type `Parameters` are used to define objective functions and their 
respective parameters for finding optimal two-stage designs.
Every `Parameters` object also contains a [`SampleSpace`](@ref) which restricts 
the feasible region of the sample space and implements a [`score`](@ref) method 
which is to be optimized.

```@docs
Parameters

label(par::Parameters)

null(par::Parameters)

mtoer(par::Parameters)

mcrv(par::Parameters)

samplespace(par::Parameters)
```

## MESS

```@docs
MESS
```

## LiuScore

```@docs
LiuScore
```

## EB

```@docs
EB
```

## MBESS

```@docs
MBESS
```