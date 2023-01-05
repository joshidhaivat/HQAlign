import sys
import os
import all_functions as func
import numpy as np
import Levenshtein as L
import multiprocessing as mp
from multiprocessing import Pool
import time
from functools import partial
import pdb
from cigar import Cigar
import re
import math

ref = {}
rc_ref = {}
reads = {}
rc_reads = {}
acgt_alignments = {}
kmer = []
new_dir = []
outfile = []
ref_dict = {}
reads_dict = {}
minimap2_threads = 30
batch_size = 50000

# Analyse the alignment function
def analyse_alignment(filename):
    fid = open(filename,'r')
    out = []
    for aline in fid:
        cols = aline.split('\t')
        read_length = int(cols[1])
        read_start = int(cols[2])
        read_end = int(cols[3])
        ref_info = cols[5]
        ref_info = ref_info.split('_')
        chr_num = int(ref_info[1])
        chr_length = int(ref_info[2][7:])
        chr_start = int(ref_info[3][6:])
        chr_end = int(ref_info[4][4:])
        chr_strand = int(ref_info[5][7])
        if chr_strand==0:
            cols[4] = '-'
            cols[2] = str(read_length-read_end)
            cols[3] = str(read_length-read_start)
        cols[5] = 'ref_'+str(chr_num)
        cols[6] = str(chr_length)
        cols[7] = str(int(cols[7])+chr_start)
        cols[8] = str(int(cols[8])+chr_start)
        out.append('\t'.join(cols))
    fid.close()
    return out

def reverse_cigar(cigar_str):
    cigar_str = list(Cigar(cigar_str).items())
    out = ''
    cigar_str = cigar_str[::-1]
    for acig in cigar_str:
        out += str(acig[0])
        out += str(acig[1])
    return out

def reverse_cigar_SA(cigar_str):
    cigar_str = list(Cigar(cigar_str).items())
    out = ''
    rev_cigar_str = cigar_str[::-1]
    for i in range(len(cigar_str)):
        if (i == 0) or (i == len(cigar_str)-1):
            acig = rev_cigar_str[i]
        else:
            acig = cigar_str[i]
        out += str(acig[0])
        out += str(acig[1])
    return out

def cigar_info(cigar_str):
    cigar_str = list(Cigar(cigar_str).items())
    M = 0
    D = 0
    N = 0
    for acig in cigar_str:
        if acig[1] == 'M':
            M += int(acig[0])
        elif acig[1] == 'D':
            D += int(acig[0])
        elif acig[1] == 'N':
            N += int(acig[0])
    return M+D+N-2

def reverse_MD(md_str):
    pattern = '(\d+)|(\^[A-Za-z]+)|([A-Za-z])'
    md_str = re.split(pattern,md_str)
    md_str = md_str[::-1]
    out = ''
    for achar in md_str:
        if achar != None:
            out += str(achar)
    return out

def analyse_alignment_sam(filename1,outfilename,filename2=None):
    if filename2 == None:
        filtering = False
    else:
        filtering = True
    fid1 = open(filename1,'r')
    out = []
    paf_alignments = []
    if filtering:
        fid2 = open(filename2,'r')
        for i,aline in enumerate(fid2):
            paf_alignments.append(aline)
        fid2.close()
    header_lines = 0
    for i,aline in enumerate(fid1):
        if filtering==True:
            cols = aline.split('\t')
            if cols[0][0] == '@':
                header_lines += 1
                continue
            # cols1 = paf_alignments[i-header_lines].split('\t')
            flag = int(cols[1])
            if flag == 4:
                continue
            read_info = cols[0]
            ref_info = cols[2]
            ref_info = ref_info.split('&')
            chr_num = ref_info[1]
            chr_length = int(ref_info[2][7:])
            chr_start = int(ref_info[3][6:])
            chr_end = int(ref_info[4][4:])
            chr_strand = int(ref_info[5][7])
            cols[0] = read_info
            if chr_strand == 0:
                cols[1] = str(16+flag)
                align_length = cigar_info(cols[5])
                cols[3] = str(chr_length-(chr_start+int(cols[3])+align_length))
                cols[5] = str(reverse_cigar(cols[5]))
            else:
                cols[3] = str(int(cols[3])+chr_start)
            cols[2] = str(chr_num)
            for ii in range(11,len(cols)):
                if cols[ii][:5] == 'MD:Z:':
                    if chr_strand == 0:
                        mdstring = cols[ii][5:]
                        cols[ii] = 'MD:Z:'+reverse_MD(mdstring)
                    else:
                        continue
                if cols[ii][:5] == 'SA:Z:':
                    supp_list = cols[ii][5:].split(';')
                    for jj in range(len(supp_list)-1):
                        supp_info = supp_list[jj].split(',')
                        ref_info = supp_info[0].split('&')
                        chr_num = ref_info[1]
                        chr_length = int(ref_info[2][7:])
                        chr_start = int(ref_info[3][6:])
                        chr_end = int(ref_info[4][4:])
                        chrom_strand = int(ref_info[5][7])
                        supp_info[0] = chr_num
                        if chrom_strand==1:
                            supp_info[1] = str(int(supp_info[1])+chr_start)
                            supp_info[2] = '+'
                        else:
                            align_length = cigar_info(supp_info[3])
                            supp_info[1] = str(chr_length-(int(supp_info[1])+chr_start+align_length))
                            supp_info[2] = '-'
                            supp_info[3] = reverse_cigar(supp_info[3])
                        supp_list[jj] = ','.join(supp_info)
                    cols[ii] = 'SA:Z:'+';'.join(supp_list)
                    continue
            out.append('\t'.join(cols))
        else:
            cols = aline.split('\t')
            if cols[0][0] == '@':
                header_lines += 1
                continue
            flag = int(cols[1])
            if flag == 4:
                continue
            bin_flag = bin(flag)
            read_info = cols[0]
            ref_info = cols[2]
            ref_info = ref_info.split('&')
            chr_num = ref_info[1]
            chr_length = int(ref_info[2][7:])
            chr_start = int(ref_info[3][6:])
            chr_end = int(ref_info[4][4:])
            chr_strand = int(ref_info[5][7])
            cols[0] = read_info
            if chr_strand == 0:
                if flag >= 16 and bin_flag[-5] == '1':
                    cols[1] = str(flag-16)
                else:
                    cols[1] = str(flag+16)
                cols[3] = str(chr_length-(chr_start+int(cols[3])+cigar_info(cols[5])))
                cols[5] = str(reverse_cigar(cols[5]))
            else:
                cols[3] = str(int(cols[3])+chr_start)
            cols[2] = str(chr_num)
            for ii in range(11,len(cols)):
                if cols[ii][:5] == 'MD:Z:':
                    if chr_strand == 0:
                        mdstring = cols[ii][5:]
                        cols[ii] = 'MD:Z:'+reverse_MD(mdstring)
                    else:
                        continue
                if cols[ii][:5] == 'SA:Z:':
                    supp_list = cols[ii][5:].split(';')
                    for jj in range(len(supp_list)-1):
                        supp_info = supp_list[jj].split(',')
                        if supp_info[2] == '+':
                            local_strand = 1
                        else:
                            local_strand = 0
                        ref_info = supp_info[0].split('&')
                        chr_num = ref_info[1]
                        chr_length = int(ref_info[2][7:])
                        chr_start = int(ref_info[3][6:])
                        chr_end = int(ref_info[4][4:])
                        chrom_strand = int(ref_info[5][7])
                        supp_info[0] = chr_num
                        if  not chrom_strand^local_strand:
                            #supp_info[1] = str(int(supp_info[1])+chr_start)
                            supp_info[2] = '+'
                        else:
                            #supp_info[1] = str(chr_length-(int(supp_info[1])+chr_start+cigar_info(supp_info[3])))
                            supp_info[2] = '-'
                            #supp_info[3] = reverse_cigar(supp_info[3])
                        if chrom_strand == 1:
                            supp_info[1] = str(int(supp_info[1])+chr_start)
                        else:
                            supp_info[1] = str(chr_length-(int(supp_info[1])+chr_start+cigar_info(supp_info[3])))
                            supp_info[3] = reverse_cigar_SA(supp_info[3])
                        supp_list[jj] = ','.join(supp_info)
                    cols[ii] = 'SA:Z:'+';'.join(supp_list)
                    continue
            out.append('\t'.join(cols))
    fid1.close()
    fid = open(outfilename,'a')
    for aline in out:
        fid.write(aline)
    fid.close()
    return

def filter_chrom(primary_list,new_member):
    if primary_list == []:
        primary_list.append(new_member)
    else:
        write = 0
        ovp_index = []
        for ii,amember in enumerate(primary_list):
            if amember[0] != new_member[0] or amember[4] != new_member[4]:
                continue
            else:
                ovp = min(amember[3],new_member[3])-max(amember[2],new_member[2])
                if ovp < 0:
                    continue
                else:
                    primary_list[ii][2] = min(amember[2],new_member[2])
                    primary_list[ii][3] = max(amember[3],new_member[3])
                    write += 1
                    ovp_index.append(ii)
        if write == 0:
            primary_list.append(new_member)
        if write > 1:
            start_index = []
            end_index = []
            # print('new mem: {}'.format(new_member))
            for i,ii in enumerate(ovp_index):
                # print('primary list {}-{}={}: {}'.format(ii,i,ii-i,primary_list[ii-i]))
                start_index.append(primary_list[ii-i][2])
                end_index.append(primary_list[ii-i][3])
                primary_list.remove(primary_list[ii-i])
            new_member[2] = min(start_index)
            new_member[3] = max(end_index)
            # print(new_member)
            # pdb.set_trace()
            primary_list.append(new_member)
    return primary_list

def hybrid_alignment(i):
    global ref,rc_ref,reads,rc_reads,acgt_alignments,kmer,new_dir,ref_dict,reads_dict,outfile,minimap2_threads,batch_size
    print_read_name = []
    print_ref_name = []
    print_read = []
    print_rc_read = []
    print_ref = []
    for_list = []
    pid = str(os.getpid())
    subset_keys = [*acgt_alignments][batch_size*(i):min(len([*acgt_alignments]),batch_size*(i+1))]
    for akey in subset_keys:
        alignments_list = acgt_alignments[akey]
        for analign in alignments_list:
            r_num = analign[0]
            r_length = int(analign[1])
            r_start = int(analign[2])
            r_end = int(analign[3])
            if analign[4] == '+':
                strand = 1
            else:
                strand = 0
            ref_chr = analign[5]
            chr_length = int(analign[6])
            ref_start = int(analign[7])
            ref_end = int(analign[8])
            matched_bases = int(analign[9])
            total_bases = int(analign[10])
            mapping_quality = int(analign[11])
            percent_unaligned = (r_start+r_length-r_end)/r_length + 0.25
            bases_appended = int(percent_unaligned*r_length)
            name_string = '>'+str(r_num)+'\n'
            if name_string not in print_read_name:
                print_read_name.append(name_string)
                print_read.append(reads[r_num])
                print_rc_read.append(rc_reads[r_num])
            if bases_appended>ref_start and ref_end+bases_appended>chr_length:
                new_mem = [ref_chr,chr_length,0,chr_length,strand]
            elif bases_appended>ref_start and ref_end+bases_appended<=chr_length:
                new_mem = [ref_chr,chr_length,0,ref_end+bases_appended,strand]
            elif bases_appended<=ref_start and ref_end+bases_appended>chr_length:
                new_mem = [ref_chr,chr_length,ref_start-bases_appended,chr_length,strand]
            else:
                new_mem = [ref_chr,chr_length,ref_start-bases_appended,ref_end+bases_appended,strand]
            for_list = filter_chrom(for_list,new_mem)
    for amem in for_list:
        if amem[4] == 1:
            print_ref_name.append('>ref&'+str(amem[0])+'&length='+str(amem[1])+'&start='+str(amem[2])+'&end='+str(amem[3])+'&strand='+str(amem[4])+'\n')
            print_ref.append(ref[amem[0]][amem[2]:amem[3]])
        else:
            print_ref_name.append('>ref&'+str(amem[0])+'&length='+str(amem[1])+'&start='+str(amem[1]-amem[3])+'&end='+str(amem[1]-amem[2])+'&strand='+str(amem[4])+'\n')
            print_ref.append(rc_ref[amem[0]][amem[1]-amem[3]:amem[1]-amem[2]])
    # Print reads and ref section to files
    func.print_reads_to_fasta(print_read,new_dir+'test_read_'+pid+'.fasta',print_read_name)
    func.print_reads_to_fasta(print_rc_read,new_dir+'test_rc_read_'+pid+'.fasta',print_read_name)
    func.print_ref_to_fasta(print_ref,new_dir+'test_ref_'+pid+'.fasta',print_ref_name)
    # Perform alignment using minimap2
    os.system('softwares/minimap2-2.24_hqalign/minimap2 -ax map-ont -t '+str(minimap2_threads)+' --MD -k '+str(kmer)+' -q '+new_dir+'test_rc_read_'+pid+'.fasta '+new_dir+'test_ref_'+pid+'.fasta '+new_dir+'test_read_'+pid+'.fasta > '+new_dir+'test_align_'+pid+'.sam')
    analyse_alignment_sam(new_dir+'test_align_'+pid+'.sam',outfile)
    return

def load_required_reads(filename,required_reads_list=[]):
    genome = []
    filename = str(filename)
    fid = open(filename,'r')
    seq=''
    name=[]
    for aline in fid:
        if(aline[0]=='>'):
            name.append(aline)
            if not seq:
                continue
            else:
                genome.append(seq)
            seq=''
        else:
            seq+=aline[:-1]
    genome.append(seq)
    fid.close()
    r = {}
    for i,aname in enumerate(name):
        akey = aname[1:-1].split(' ')[0]
        if (akey in required_reads_list) or (required_reads_list == []):
            r[akey] = genome[i]
    return r

if (__name__ == "__main__"):
    t1 = time.time()
    foldername = str(sys.argv[1])
    ref_filename = str(sys.argv[2])
    rc_ref_filename = str(sys.argv[3])
    initial_filename = str(sys.argv[4])
    forward_reads_dir = str(sys.argv[5])
    initial_kmer = int(sys.argv[6])
    kmer = int(sys.argv[7])
    PROCESS = int(sys.argv[8])
    sam_filename = str(sys.argv[9])
    outfile = str(sys.argv[10])#foldername+'align_hq2_'+str(initial_kmer)+'_'+str(kmer)+'_list-3.sam'
    reverse_reads_dir = str(sys.argv[11])
    ref = load_required_reads(ref_filename)
    rc_ref = load_required_reads(rc_ref_filename)
    fid = open(initial_filename,'r')
    for aline in fid:
        cols = aline.split('\t')
        try:
            acgt_alignments[cols[0]].append(cols[:12])
        except:
            acgt_alignments[cols[0]] = []
            acgt_alignments[cols[0]].append(cols[:12])
    fid.close()
    for afilename in os.listdir(forward_reads_dir):
        readfilename = '{}{}'.format(forward_reads_dir,afilename)
        reads.update(load_required_reads(readfilename,acgt_alignments.keys()))
    for afilename in os.listdir(reverse_reads_dir):
        readfilename = '{}{}'.format(reverse_reads_dir,afilename)
        rc_reads.update(load_required_reads(readfilename,acgt_alignments.keys()))
    # acgt_alignments = func.extract_from_paf_v1(initial_filename,reads_dict,ref_dict)
    threads1 = 1 #math.ceil(PROCESS/30)
    minimap2_threads = math.floor(PROCESS/threads1)
    alignments = []
    new_dir = 'test_dbg_alignments/'
    if os.path.isdir(new_dir):
        os.system('rm -rf {}'.format(new_dir))
    os.system('mkdir '+new_dir)
    wfid = open(outfile,'w')
    append_lines = []
    fid = open(sam_filename,'r')
    for aline in fid:
        if aline[0]=='@':
            wfid.write(aline)
    fid.close()
    wfid.close()
    for i in range(math.ceil(len([*acgt_alignments])/batch_size)):
        hybrid_alignment(i)
    os.system('rm -rf '+new_dir)
    print('Hybrid alignment finished in {:.2f} seconds'.format(time.time()-t1))