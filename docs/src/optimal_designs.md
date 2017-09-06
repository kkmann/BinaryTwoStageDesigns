# Optimal Two-Stage Designs

The technical background in available at [arxiv.org](https://arxiv.org/abs/1605.00249).

```@docs
optimaldesign(n1::Integer, parameters::Parameters, solver::MathProgBase.AbstractMathProgSolver)

optimaldesign(
  parameters::Parameters,
  solver::MathProgBase.AbstractMathProgSolver;
  VERBOSE::Integer = 1,
  EARLYTERMINATION::Bool = false
)
```
