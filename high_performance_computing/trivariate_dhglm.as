"bash script for the trivariate with double hierarchical generalized linear model. Run in ASReml (ABG thesis)" 
!DEBUG !RENAME !ARG 1 2 3 4
!DOPART $1

Here Gval is for tnb
  id 56844 				!P		
  par_cycle	32			!A
  re_isn 2					!A
  itl 	19				!I 
  ltr_pre_cbd 2				!A 
  hym_f 5046				!A 
  ltr_sire 5924				!I 
  pee 56844 				!I 
  hys_first_isn 1768			!A
  lgp_cce 2				!A
  tnb
  nba
  AgeCull 
  ST12 
  ST15 
  HR 
  Gval !=tnb !-14.06477 !*V17
  Ywt !=1. Gwt !=1. hrW !=1.

All.ped !MAKE
Data.txt !maxit 1000 !Skip 1 !MVINCLUDE 

!Part 1
!ASUV !EXTRA 100 !SLOW

# in odd iterations, we use the predicted weights for the primary response
!IF ODD !CALC W1=EXP(R2-Y2) #redefine weights for Y1

!IF EVEN !CALC S1=1./W1; H0=MIN(H1/S1, .9999); Z2=MAX(R1*R1,.0001)/(1-H0)
!IF EVEN !CALC Y2=LOG(S1)+(Z2-S1)/S1 #redefine Y2
!IF EVEN !CALC W2=(1-H0)/2 #redefine weights for Y2

!ASSIGN USIV 1.71788 0.005 0.03 0.005 0.005 0.01979 !GPPPPPP
!ASSIGN sire 0.02396 0.005 0.03 0.0 0.0 0.000001 !GPPPPPP
!ASSIGN pe 0.66795 0.005 0.03 0.0 0.0 0.000001 !GPPPPPP

tnb Gval AgeCull !Weight Ywt !WT Gwt !WT hrW ~ Trait ,
at(Tr,1).par_cycle ,
at(Tr,1).re_isn ,
at(Tr,1).itl ,
at(Tr,1).ltr_pre_cbd ,
at(Tr,1).hym_f ,
at(Tr,2).hys_first_isn at(Tr,2).lgp_cce ,
!r us(Trait,$USIV).pig_id ,
us(Trait,$sire).ltr_sire ,
us(Trait,$pe).pee !f mv
residual units.id(Trait)     !VARIANCESCALE

!Part 2
!ASUV !EXTRA 100 !SLOW

# in odd iterations, we use the predicted weights for the primary response
!IF ODD !CALC W1=EXP(R2-Y2) #redefine weights for Y1

!IF EVEN !CALC S1=1./W1; H0=MIN(H1/S1, .9999); Z2=MAX(R1*R1,.0001)/(1-H0)
!IF EVEN !CALC Y2=LOG(S1)+(Z2-S1)/S1 #redefine Y2
!IF EVEN !CALC W2=(1-H0)/2 #redefine weights for Y2

!ASSIGN USIV 1.71788 0.005 0.03 0.005 0.005 0.70351 !GPPPPPP
!ASSIGN sire 0.02396 0.005 0.03 0.0 0.0 0.000001 !GPPPPPP
!ASSIGN pe 0.66795 0.005 0.03 0.0 0.0 0.000001 !GPPPPPP

tnb Gval HR !Weight Ywt !WT Gwt !WT hrW ~ Trait ,
at(Tr,1).par_cycle ,
at(Tr,1).re_isn ,
at(Tr,1).itl ,
at(Tr,1).ltr_pre_cbd ,
at(Tr,1).hym_f ,
at(Tr,2).hys_first_isn at(Tr,2).lgp_cce ,
!r us(Trait,$USIV).pig_id ,
us(Trait,$sire).ltr_sire ,
us(Trait,$pe).pee !f mv
residual units.id(Trait)     !VARIANCESCALE

!Part 3
!ASUV !EXTRA 100 !SLOW

# in odd iterations, we use the predicted weights for the primary response
!IF ODD !CALC W1=EXP(R2-Y2) #redefine weights for Y1

!IF EVEN !CALC S1=1./W1; H0=MIN(H1/S1, .9999); Z2=MAX(R1*R1,.0001)/(1-H0)
!IF EVEN !CALC Y2=LOG(S1)+(Z2-S1)/S1 #redefine Y2
!IF EVEN !CALC W2=(1-H0)/2 #redefine weights for Y2

!ASSIGN USIV 1.71788 0.005 0.03 0.005 0.005 0.07906 !GPPPPPP
!ASSIGN sire 0.02396 0.005 0.03 0.0 0.0 0.000001 !GPPPPPP
!ASSIGN pe 0.66795 0.005 0.03 0.0 0.0 0.000001 !GPPPPPP

tnb Gval HR !Weight Ywt !WT Gwt !WT hrW ~ Trait ,
at(Tr,1).par_cycle ,
at(Tr,1).re_isn ,
at(Tr,1).itl ,
at(Tr,1).ltr_pre_cbd ,
at(Tr,1).hym_f ,
at(Tr,2).hys_first_isn at(Tr,2).lgp_cce ,
!r us(Trait,$USIV).pig_id ,
us(Trait,$sire).ltr_sire ,
us(Trait,$pe).pee !f mv
residual units.id(Trait)     !VARIANCESCALE

!Part 4
!ASUV !EXTRA 100 !SLOW

# in odd iterations, we use the predicted weights for the primary response
!IF ODD !CALC W1=EXP(R2-Y2) #redefine weights for Y1

!IF EVEN !CALC S1=1./W1; H0=MIN(H1/S1, .9999); Z2=MAX(R1*R1,.0001)/(1-H0)
!IF EVEN !CALC Y2=LOG(S1)+(Z2-S1)/S1 #redefine Y2
!IF EVEN !CALC W2=(1-H0)/2 #redefine weights for Y2

!ASSIGN USIV 1.71788 0.005 0.03 0.005 0.005 0.00102 !GPPPPPP
!ASSIGN sire 0.02396 0.005 0.03 0.0 0.0 0.000001 !GPPPPPP
!ASSIGN pe 0.66795 0.005 0.03 0.0 0.0 0.000001 !GPPPPPP

tnb Gval HR !Weight Ywt !WT Gwt !WT hrW ~ Trait ,
at(Tr,1).par_cycle ,
at(Tr,1).re_isn ,
at(Tr,1).itl ,
at(Tr,1).ltr_pre_cbd ,
at(Tr,1).hym_f ,
at(Tr,2).hys_first_isn at(Tr,2).lgp_cce ,
!r us(Trait,$USIV).pig_id ,
us(Trait,$sire).ltr_sire ,
us(Trait,$pe).pee !f mv
residual units.id(Trait)     !VARIANCESCALE

