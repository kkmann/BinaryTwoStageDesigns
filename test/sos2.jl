using JuMP

m = Model()
pivots = collect(linspace(0, 5, 6))
@variable(m, 0 <= lambda[piv in pivots] <= 1)
addSOS2(m, lambda)
@variable(m, f)
@variable(m, x)
@constraint(m, sum(lambda[piv]*(piv - 1)^2 for piv in pivots) == f)
@constraint(m, sum(lambda[piv]*piv for piv in pivots) == x)
@constraint(m, sum(lambda[piv] for piv in pivots) == 1)

@objective(m, Min, f)
solve(m)

getvalue(f)
getvalue(x)
