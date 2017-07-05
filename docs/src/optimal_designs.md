# Optimal Two-Stage Designs

The principle technical background in available at [arxiv.org](https://arxiv.org/abs/1605.00249).

```@docs
getoptimaldesign{T<:Integer, TS<:MathProgBase.AbstractMathProgSolver}(
    n1::T,
    parameters::Parameters,
    solver::TS;
    VERBOSE::Integer = 0
)

getoptimaldesign{T<:Integer, TS<:MathProgBase.AbstractMathProgSolver}(
    parameters::Parameters,
    solver::TS;
    VERBOSE::Integer = 1
)
```
