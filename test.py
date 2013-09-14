from sage.all import *
from dirichlet_conrey import *

def test_conrey(nranges=[(20,100),(200,1000),(3000,100000)]):
  for n,r in nranges:
      for i in xrange(n):
          m = randint(1,r)
          G = DirichletGroup_conrey(m)
          inv = G.invariants()
          cinv = [ Mod(g,m).multiplicative_order() for g in G.gens() ]
          try:
              assert  prod(inv) == G.order()
              assert  tuple(cinv) == inv
          except:
              return m, inv,cinv, prod(inv)
