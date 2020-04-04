# spatial-epi
A collection of bits and pieces for thinking through how to do spatial epidemic spread modelling.

In particular a series of Netlogo models, as follows. Two (very) toy models exploring self isolation 'bubbles'
+ [bubbles.nlogo](bubbles.nlogo)
+ [nested-bubbles.nlogo](nested-bubbles.nlogo)

and more substantively, a series of connected localised SEIR models
+ [distributed-seir.nlogo](distributed-seir.nlogo) has random connections among a set of equal-sized locales
+ [distributed-seir-02.nlogo](distributed-seir-02.nlogo) has more spatially coherent connections among the same
+ [distributed-seir-03.nlogo](distributed-seir-03.nlogo) is an attempt to enable the model to read in real health management zones

All three of the above have excessive mortality, which has been corrected in the later models:
+ [distributed-seir-04.nlogo](distributed-seir-04.nlogo) reverts back to the main sequence and allows for locales to vary in size and variance, under paramterisable control
+ [distributed-seir-05.nlogo](distributed-seir-05.nlogo) currently the same as 04, but will add lockdown levels and testing
