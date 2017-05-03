#!/usr/bin/env python
#Request to extract objects from webpages
#argv[1] : Codes_labels_inData
import sys
import re
import urllib2,sys
#from sys import argv
import pandas as pd
from tempfile import TemporaryFile
#from rpy2.robjects import pandas2ri
#pandas2ri.activate()

#from rpy2.robjects.packages import importr


with open("Codes_labels_inData_short") as f:
    lines = f.read().splitlines()

allCodesR = lines

D={}
for i in allCodesR:
    url = 'http://www.uniprot.org/uniprot/%s'%(i)
    u = []
    try:
        html = urllib2.urlopen(url).read()
        d=re.findall(r"IPR\d{6}",html)
        u = list(set(d))
    except:
        pass
    D[i] = u


df = pd.DataFrame.from_dict(D, orient='index')


#base = importr('base')
## call an R function on a Pandas DataFrame
#base.summary(df)

#df.to_csv("domains_no_argv")



