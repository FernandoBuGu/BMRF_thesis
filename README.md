# Fernando_BMRF_thesis
Protein function prediction

Code to calculate AUC from code in https://github.com/jwbargsten/bmrf 
Folders "R" and "reference_implementation_yiannis_kourmpetis" are from https://github.com/jwbargsten/bmrf.

For validation, the test is a portion of GO terms for which all proteins except a minimum are masked. 

results.txt shows results in 108 experiments. For each, I average over 10 replicates. The parameters that I have changed are:
minGOsize (also a parameter in code https://github.com/jwbargsten/bmrf);
Whether only GO terms related to biological processes are considered ("1" if so, "0" otherwise);
Whether only GO terms with experimental evidence scores are considered ("1" if so, "0" otherwise);
The portion of GO terms to be masked to minGOsize;
The size of the network (#conexions):
  S1 (10,973), S2 (26,879), S3 (26,774), S4 (64,519);
  Int: "integrated network" (242,504);
  complet:all subsets and integrated together (111,390 after removing duplicates)
  ppi data from old code (401,820)
  
 Results that I considered (columns in Results.txt):
 mean AUC among all GO terms for which prediction was possible within those in the train set;
 GOs with AUC: # GO terms for which prediction was possible within those in the train set;
 sd_AUCs: sdeviation;
 no.go.above.AUC80. # GO terms for which AUC>80%;
 no.unnanotated.prots. #labels that were not included in the train set;
 total.prots in the network;
 total.predicted.go. #GO terms for which prediction was made using bmrf.R. Only for those in test set, AUC was computed.
 
