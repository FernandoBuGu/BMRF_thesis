"bash script for  grill format of the trivariate with double hierarchical generalized linear model. Run in ASReml (ABG thesis)" 
!RENAME !ARG 1 2
!DOPART $1
univs model
  id 44919			!P
  tnbM	
  nbaM
  AgeCull
  ST15
  mean_parities
  hym_1stisn 1631		!A
  no_records
  ST12
  
animal.ped !MAKE
DATaF.txt !Extra 4 !maxit 1000 !Skip 1 !MVREMOVE

!PART 1
tnbM !WT no_records AgeCull ~ Trait ,
Tr.hym_1stisn at(Tr,1).mean_parities ,
!r Tr.pig_id

1 2 1
0 0 0 !S2==1
Trait 0 US !GP
9.06109 0 0.778522     

Tr.pig_id 2
Trait 0 US !GP
1.35731 0 0.0695157 
pig_id 0 AINV

VPREDICT !DEFINE
F phenvar 1:3 + 4:6
F addvar 4:6 * 1
H heritA 10 7
H heritB 12 9
R phencorr 7 8 9
R gencorr 4:6

!PART 2
nbaM !WT no_records AgeCull ~ Trait ,
Tr.hym_1stisn at(Tr,1).mean_parities ,
!r Tr.pig_id

1 2 1
0 0 0 !S2==1
Trait 0 US !GP
9.39451 0 0.778522     

Tr.pig_id 2
Trait 0 US !GP
0.756541 0 0.0695157 
pig_id 0 AINV

VPREDICT !DEFINE
F phenvar 1:3 + 4:6
F addvar 4:6 * 1
H heritA 10 7
H heritB 12 9
R phencorr 7 8 9
R gencorr 4:6


