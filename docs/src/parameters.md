# Design Parameters

## Overview

Objects of type `Parameters` are used to store parameters required for finding
optimal two-stage designs.
Every `Parameters` object contains a sample space which restricts the feasible
region of the sample space and implements a [`score`](@ref) method which is to
be optimized.

```@docs
Parameters

PointAlternative
```

## SimpleMinimalExpectedSampleSize

Currently, only a single parameter type (and thus objective function)
is implemented. A `SimpleMinimalExpectedSampleSize` objects holds information
about the null hypothesis and a point alternative on which the power is
calculated as well as allowing to specify a response rate under which the
expected sample size is to be minimized.

```@docs
SimpleMinimalExpectedSampleSize
```
