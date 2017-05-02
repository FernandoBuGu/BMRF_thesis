# Fernando_BMRF_thesis
Protein function prediction

In this repository the code to calculate AUC from code in https://github.com/jwbargsten/bmrf 
Folders "R" and "reference_implementation_yiannis_kourmpetis" are from https://github.com/jwbargsten/bmrf.

For validation, the test is a portion of GO terms for which all proteins except a minimum are masked. 

Results.txt shows results in 100 experiments. The parameters that I have changed are:
minGOsize (also a parameter in code https://github.com/jwbargsten/bmrf)
Whether only Go terms with experimental evidence scores are considered ("1" if so, "0" otherwise)
Whether only Go terms related to biological processes are considered ("1" if so, "0" otherwise)
The portion of GO terms to be masked to minGOsize
The size of the network (#conexions):
  S1 (10973), S2 (26879), S3 (26774), S4 (64519)
  Int: "integrated network" (242,504)
  complet:all subsets and integrated together (111,390 after removing duplicates)
  
 Results that I considered (columns in Results.txt):
 mean AUC among all GO terms for which prediction was possible within those in the train set;
 GOs with AUC: # GO terms for which prediction was possible within those in the train set;
 sd_AUCs: sdeviation;
 no.go.above.AUC80. # GO terms for which AUC>80%;
 no.unnanotated.prots. #labels that were not included in the train set;
 total.prots in the network;
 total.predicted.go. #GO terms for which prediction was made using bmrf.R. Only for those in test set, AUC was computed.
 
