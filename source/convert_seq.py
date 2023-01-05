import numpy as np
import sys
import scipy.io as sio
import time
from multiprocessing import Pool
import os
import pdb
import Levenshtein as L
import all_functions as func
from functools import partial
import multiprocessing as mp

reads = []
name = []
level = []
qlevels = []
kmermap = []
rc = []

def qlevels_dict(qlevel):
    dict_qlevels={}
    for i in range(len(qlevel)):
        dict_qlevels[len(qlevel[i])-1] = qlevel[i]
    return dict_qlevels

def quantize_seq(i):
    global reads,name,level,qlevels,kmermap,rc
    kmer_k = 6
    seq_name = name[i]
    seq = reads[i]
    current_seq = []
    reads_q = []
    loop = 1
    if rc==True:
        rc_current_seq = []
        rc_reads_q = []
        loop = 2
    for k in range(0,loop):
        if(k==1):
            seq = func.revcom(seq)
        current_mean = np.zeros(len(seq)-kmer_k+1)
        for j in range(0,len(seq)-kmer_k+1):
            try:
                current_mean[j] = kmermap[seq[j:j+6].upper()][0]
            except:
                current_mean[j] = -1
        if(k==0):
            current_seq = current_mean
        else:
            rc_current_seq = current_mean
    reads_q = func.get_quantized_seq(current_seq,qlevels[level])
    if rc==True:
        rc_reads_q = func.get_quantized_seq(rc_current_seq,qlevels[level])
        return [seq_name,reads_q,rc_reads_q]
    else:
        return [seq_name,reads_q]

if (__name__ == "__main__"):
    filename = str(sys.argv[1])
    level = int(sys.argv[2])
    rc = int(sys.argv[3])
    if rc == 0:
        rc = False
    else:
        rc = True
    MAX_PROCESS = int(sys.argv[4])
    out_dir = str(sys.argv[5])
    reads_ref_indicator = int(sys.argv[6])
    if reads_ref_indicator == 0:
        reads,num_reads,name = func.get_reads_from_fasta(filename)
    else:
        reads,num_reads,name = func.get_genome_from_fasta(filename)
    import pickle
    with open('qlevels.txt','rb') as f:
        qlevel = pickle.load(f)
    qlevels = qlevels_dict(qlevel)
    kmermap = func.get_kmer_map('qmer_map/r9.4_6mer_nucleotide_template_model.txt')
    print('\tConverting nucleotide seqs to quantized seqs\n\tNumber of seqs = {}\t\tlevel = {}\trc = {}'.format(num_reads,level,rc))
    start = time.time()
    output = []
    # f = partial(quantize_seq,reads_name=name,reads=reads,level=2,qlevels=qlevels_dict(qlevels),kmermap=kmermap)
    p = Pool(MAX_PROCESS)
    output += p.map(quantize_seq,range(len(reads)))
    p.close()
    p.join()
    print('\tTime taken to convert = {}'.format(time.time()-start))
    print('\tPerforming sanity check...')
    print_name = []
    print_reads = []
    if rc == True:
        print_rc_reads = []
    for i in range(len(output)):
        print_name.append(output[i][0])
        print_reads.append(output[i][1])
        if rc == True:
            print_rc_reads.append(output[i][2])
        if len(reads[i])-len(output[i][1]) != 5:
            exit(1)
    print('\tWriting files...')
    if os.path.isdir(out_dir) == False:
        os.system('mkdir {}'.format(out_dir))
    if reads_ref_indicator == 0:
        outfilename = sys.argv[7]#out_dir+'reads_q{}.fasta'.format(level)
        if rc == True:
            outfilename_rc = sys.argv[8]#out_dir+'rc_reads_q{}.fasta'.format(level)
    else:
        outfilename = sys.argv[7]#out_dir+'ref_q{}.fasta'.format(level)
        if rc == True:
            outfilename_rc = sys.argv[8]#out_dir+'rc_ref_q{}.fasta'.format(level)
    func.print_reads_to_fasta(print_reads,outfilename,print_name)
    if rc == True:
        func.print_reads_to_fasta(print_rc_reads,outfilename_rc,print_name)
    print('\tDone!')
