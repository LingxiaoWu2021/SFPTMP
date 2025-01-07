## FOLDER STRUCTURE:
- Dev[a]-[b]P[c]S：
  1. a: deviation parameter \Delta=a/100;
  2. b: number of stages;
  3. c: number of scenarios.

## PARAMETERS:
- Variables:
  1. CNN:  index of the instance in a group  
  2. P:    number of stages
  3. PT:   number of periods in a stage
  4. T:    number of periods
  5. I1:   number of suppliers
  6. I2:   number of customers
  7. I:    number of sites
  8. BN:   number of bids
  9. LBN:  number of bids on each arc
  10.SPN:  number of shipments
  11.SN:   number of scenarios in each stage
  12.MBSN: (maximum) number of shipments in a bid

- Vectors and Matrices:
  1.  PRO[P][SN]:         probability of each scenario at each phase 
  2.  PTN[P]:             number of periods in each stage
  3.  PTset[P][PT]:       set of periods in each phase
  4.  sup[I][P][SN][PT]:  supply (positive for supply and negative for demand) under different scenarios
  5.  iniv[I]:            initial inventory at each site
  6.  ubiv[I]:            the upper bounds of inventory at each site
  7.  len[I1][I2]:        the traveling time between two sites
  8.  c1[I]:              unit inventory holding cost at each site
  9.  c2[I]:              unit backlogging cost at each site
  10. c3[BN]:             unit variable transportation cost for each bid
  11. c4[I1][I2]:         unit transportation  cost of spot market between two sites
  12. RLBN[I1][I2]:       number of bids on each arc 
  13. LBset[I1][I2][LBN]: set of bids on each arc  
  14. frt[BN]:            freight rate of each bid 
  15. lbcap[BN]:          capacity lower bound of each bid 
  16. ubcap[BN]:          capacity upper bound of each bid 
  17. SHN[BN]:            number of shipments in each bid
  18. SHset[BN][MBSN]:    set of shipments in each bid
  19. SHsts[SPN]:         set of start times of each shipment 
  20. SHets[SPN]:         set of end times of each shipment
