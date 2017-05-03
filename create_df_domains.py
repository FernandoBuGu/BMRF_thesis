#!/usr/bin/env python
#Request to extract objects from webpages
#argv[1] : Codes_labels_inData
import sys
import re
import urllib2,sys
#from sys import argv
import pandas as pd
from tempfile import TemporaryFile

with open("Codes_labels_inData") as f:
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
#print(D)


#read_dictionary = np.load('my_file.npy').item()
#print(read_dictionary)

df = pd.DataFrame.from_dict(D, orient='index')
#df.to_csv("domains_no_argv")



