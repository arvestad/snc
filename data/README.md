* `3a.dd`: 	Three somewhat uncorrelated sequences. I meant for s0 and s1 to be correlated, but they actually are not. 
  Excel says the correlations should be:
  ```s0	s1	-0,159483023
  s0	s2	-0,233120125
  s1	s2	-0,248545745
  ```
* `3b.dd`: Two perfectly correlated sequences and one uncorrelated. We have s0 and s1 that should get correlation 1.0 
  and s2 should be uncorrelated with the others.
* `3c.dd`: Smaller data set created to expose a bug: whenever two highly correlated sequences were not adjacent in
  the similarity matrix due to being far apart in the input, then they were not found as being in the same grouping.
  
  
  
Durand_GoldStandard.txt  Identifiers taken from Song _et al_ (2008), with mapping to gene family.
