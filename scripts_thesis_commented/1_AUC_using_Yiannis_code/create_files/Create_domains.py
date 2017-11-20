#!/usr/bin/env python

#Script to extract the domain information from the uniprot web page. 

    #MSc thesis bioinfomratics WUR. Protein function prediction for poorly annotated species.
    #Author: Fernando Bueno Gutierrez
    #email1: fernando.buenogutierrez@wur.nl
    #email2: fernando.bueno.gutie@gmail.com

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

df = pd.DataFrame.from_dict(D, orient='index')

