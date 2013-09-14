from sage.all import *
from dirichlet_conrey import *

def test_conrey(nranges=[(20,100),(200,1000),(500,500000)]):
  for n,r in nranges:
      for i in xrange(n):
          q = randint(1,r)
          G = DirichletGroup_conrey(q)
          inv = G.invariants()
          cinv = tuple([ Mod(g,q).multiplicative_order() for g in G.gens() ])
          try:
              assert  prod(inv) == G.order()
              assert  q <= 2 or G.zeta_order() == inv[0]
              assert  cinv == inv
          except:
              print 'group error'
              return q, inv, G
          if q > 2:
              m = 0
              while gcd(m,q) != 1:
                  m = randint(1,q)
              chi = G[m]
              n = randint(1,q)
              try:
                  assert chi.multiplicative_order() == Mod(m,q).multiplicative_order()
                  if gcd(n,q) != 1:
                      try:
                          assert chi.logvalue(n) == -1
                      except:
                          print 'non unit value error'
                          return chi, n, r
                  elif q < 10^6:
                      try:
                          ref = chi.sage_character()(n).n().real()
                          new = N(chi(n).real)
                          assert abs(ref - new) < 1e-5
                      except:
                          print 'unit value error'
                          return chi, n
              except:
                  print 'char error'
                  return chi, n, r
