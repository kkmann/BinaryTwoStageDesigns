"""
Any planning assumption must at least contain the fields n1range as the
allowable range of n1 values, nmax as maximal final n and score(design) as
function giving the score of a design (smaller = better).
"""
abstract Parameters # only guarantees p0, nmax, n1range, alpha

abstract PointAlternative <: Parameters # must also have p1 + beta

abstract UncertainAlternative <: Parameters # must implement p0 + Beta prior on p1 but no beta
