#!/usr/bin/env python

from __future__ import division
from sys import argv
import numpy
import os
import subprocess


# Author : Souvik Das
# Registration number : 920107172010

# calculate M50
# parsefrom small_
# reun to get uptpu.
# what is uncovered?



def parsing_ref(ref):
    # in:fasta file containing one seq(ref)
    # out: dicc of type {header1:seq1}
    dicc = {}
    values = []
    for line in ref:
        if line.startswith(">"):
            header = line
        else:
            values.append(line)
    v = ''.join(values).replace("\n", "")
    dicc[header] = v
    d = str(dicc.values())
    return d

def parsing_comp(m_comp_file):
    # in:fasta file containing some seqs(comp)
    # out: dicc of type {header1:seq1, header2:seq2} and list of lengs (1 value per sequence)
    huge_list = m_comp_file.read()
    H = huge_list.replace(">", "@@@>")
    splited = H.split("@@@")
    headers = []
    all_values = []
    lengs = []
    for record in splited:
        values = []
        record_splited_in_lines = record.split("\n")
        for lin in record_splited_in_lines:
            if lin.startswith(">"):
                headers.append(lin)
            else:
                values.append(lin)
        v = ''.join(values).replace("\n", "")
        lengs.append(len(v))
        all_values.append(v)
    dicc = dict(zip(headers, lengs[1:]))
    return dicc



def c_median(dicc):
    M50v = numpy.median(numpy.array(dicc.values()))
    return M50v


def find_M50index(dicc, M50v):
    seqs_sorted = sorted(dicc.values())
    above_median = []
    for idx, seq in enumerate(seqs_sorted):
        if seq > M50v:
            above_median.append(idx)
    return above_median[0]


def running_v(ref, m_comp_file, out_file="out.lastz"):
    if not os.path.exists(out_file):
        cmd = 'lastz %s %s --format=general > %s'%(ref, m_comp_file, out_file)
        subprocess.check_output(cmd, shell=True)
	return 1

def parse_out(outp_lastz):
#In: output file from lastz
#OUT: list of coord start (REF), and list of end (REF)
	starts = []
	ends = []
	

	a = outp_lastz.read()
	print a

	for line in outp_lastz:
		line_sp = line.split()
		print line_sp
		starts.append(line_sp[4])
		ends.append(line_sp[5])
	print starts
	return starts[1:], ends[1:]

def list_miss(dicc_ref,starts,ends):
#In: dicc of type {header1:seq1}
#In: list of coord start (REF), and list of end (REF)
	r = dicc_ref.replace("\n","")
	leng = len(r)
	L = []
	for i in range(len(starts)):
		st = starts[i]
		nd = ends[i]
		bitt = r[int(st):int(nd)]
		L.append(bitt)
	leng_all = len("".join(L))
	print leng
	return L, leng, leng_all 

def printing(seqs,M50v,M50i,leng,leng_all,starts,ends):
	print("%s : TOTAL=miss ;N50 SIZE=NA; N50 INDEX=NA"%(argv[1]))
	print("%s : TOTAL= miss ;N50 SIZE=%s; N50 INDEX=%s"%(argv[2], M50v, M50i))
	print
	print("Uncovered regions:")
	#for i in range(len(starts)):
	#	print("%s : %s seqs"%(starts[i],ends[i]))
	print
	print("Number of uncovered regions: %s"%(len(starts)))
	print("Number of uncovered bases: %s"%(len(leng_all)))
	
	#len(L)
	


if __name__ == '__main__':
    ref = open(argv[1])
    m_comp_file = open(argv[2])
    dicc_ref = parsing_ref(ref)
    dicc = parsing_comp(m_comp_file)
    M50v = c_median(dicc)
    M50i = find_M50index(dicc, M50v)
    #running_v(argv[1], argv[2])
    out_parsed_output =parse_out(open(argv[3]))
    starts = out_parsed_output[0]
    ends = out_parsed_output[1]
    seqs_leng_leng_all = list_miss(dicc_ref,starts,ends)
    seqs = seqs_leng_leng_all[0]
    leng = seqs_leng_leng_all[1]
    leng_all = seqs_leng_leng_all[2]
    #R = printing(seqs,M50v,M50i,leng,leng_all,starts,ends)
	

	
