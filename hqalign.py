import os,sys,getopt,pdb,argparse

def parsearguments():
    argParser = argparse.ArgumentParser()
    argParser.add_argument("-r", "--ref", help="reference genome filename in fasta format", required=True)
    argParser.add_argument("-i", "--reads", help="location of directory of read files in fasta format", required=True)
    argParser.add_argument("-o", "--output", help="location of directory of output files", required=True)
    argParser.add_argument("-t", "--threads", type=int, help="maximum number of parallel threads (default=4)")
    argParser.add_argument("-k", "--kmer", type=int, help="minimizer length for hybrid step (default=18)")
    args = argParser.parse_args()
    if args.ref == None:
        print("Please include a ref filename")
        sys.exit(1)
    elif args.reads == None:
        print("Please include a directory location for reads")
        sys.exit(2)
    elif args.output == None:
        print("Please include an output directory location for alignment files")
        sys.exit(3)
    if args.threads == None:
        args.threads = 4
    if args.kmer == None:
        args.kmer = 18
    return args

def main():
    args = parsearguments()
    if not os.path.isdir(args.output):
        os.system("mkdir {}".format(args.output))
    else:
        print('Output directory already exists. Please create a new directory for output...')
        sys.exit(4)
    ######### Initial alignments using minimap2 ##########
    print('\tPerforming initial alignment using minimap2...')
    outdir = '{}{}'.format(args.output,'minimap2_alignments/')
    os.system('mkdir {}'.format(outdir))
    for afile in os.listdir(args.reads):
        if afile.endswith('.fasta'):
            read_file = '{}{}'.format(args.reads,afile)
            samfile = '{}{}.sam'.format(outdir,afile.split('.fasta')[0])
            paffile = '{}{}.paf'.format(outdir,afile.split('.fasta')[0])
            command_line = 'minimap2 -ax map-ont --MD -t {} {} {} > {}'.format(args.threads,args.ref,read_file,samfile)
            os.system(command_line)
            command_line = 'source/softwares/k8-0.2.4/k8-Linux source/softwares/minimap2-2.24/paftools.js sam2paf {} > {}'.format(samfile,paffile)
            os.system(command_line)
    ######### Convert nucleotide deqs to Q3 seqs ##########
    print('\tConverting nucleotide ref to Q3 ref')
    outdir = '{}{}'.format(args.output,'Q3_ref/')
    os.system('mkdir {}'.format(outdir))
    outfile = '{}ref_q3.fasta'.format(outdir)
    outfile_rc = '{}rc_ref_q3.fasta'.format(outdir)
    command_line = 'python source/convert_seq.py {} {} {} {} {} {} {} {}'.format(args.ref,3,1,args.threads,outdir,1,outfile,outfile_rc)
    os.system(command_line)

    print('\tConverting nucleotide reads to Q3 reads')
    outdir = '{}{}'.format(args.output,'Q3_reads/')
    os.system('mkdir {}'.format(outdir))
    os.system('mkdir {}/forward_reads/'.format(outdir))
    os.system('mkdir {}/reverse_reads/'.format(outdir))
    for afile in os.listdir(args.reads):
        if afile.endswith('.fasta'):
            read_file = '{}{}'.format(args.reads,afile)
            outfile = '{}forward_reads/{}'.format(outdir,afile)
            outfile_rc = '{}reverse_reads/rc_{}'.format(outdir,afile)
            command_line = 'python source/convert_seq.py {} {} {} {} {} {} {} {}'.format(read_file,3,1,args.threads,outdir,0,outfile,outfile_rc)
            os.system(command_line)
    ######### Hybrid step using modified minimap2 pipeline ##########
    print('\tPerforming the hybrid step using modified minimap2 pipeline')
    workdir = '{}minimap2_alignments/'.format(args.output)
    outdir = '{}hq3_alignments/'.format(args.output)
    os.system('mkdir {}'.format(outdir))
    ref_file = '{}Q3_ref/ref_q3.fasta'.format(args.output)
    rc_ref_file = '{}Q3_ref/rc_ref_q3.fasta'.format(args.output)
    reads_file_dir = '{}Q3_reads/forward_reads/'.format(args.output)
    rc_reads_file_dir = '{}Q3_reads/reverse_reads/'.format(args.output)
    init_kmer = 15
    kmer = args.kmer
    threads = args.threads
    sam_file = ''
    for afile in os.listdir('{}minimap2_alignments/'.format(args.output)):
        if afile.endswith('.sam'):
            sam_file = '{}minimap2_alignments/{}'.format(args.output,afile)
            break
    for afilename in os.listdir(workdir):
        if afilename.endswith('.paf'):
            init_file = '{}{}'.format(workdir,afilename)
            outfile = '{}{}.sam'.format(outdir,afilename.split('.paf')[0])
            command_line = 'python source/hybrid_aligner_hq3_list_WGS_v2.py {} {} {} {} {} {} {} {} {} {} {}'.format(outdir,ref_file,rc_ref_file,
            init_file,reads_file_dir,init_kmer,kmer,threads,sam_file,outfile,rc_reads_file_dir)
            os.system(command_line)
    print('\tFinished HQAlign pipeline!')
    return

if (__name__ == "__main__"):
    main()
