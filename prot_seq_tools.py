"""
Israel Desta
Oct 28, 2021
Purpose of this module is to manipulate protein sequences in different ways
"""
import os, sys
import pdb
from math import ceil
import random
from Bio import SeqIO

def save_fasta_from_pdb(pdbfile,fastafile,chains=None):
    pdbid = os.path.splitext(os.path.basename(pdbfile))[0][:4]
    ff = open(fastafile,'w')
    with open(pdbfile) as pf:
        for record in SeqIO.parse(pf, 'pdb-atom'):
            if chains == None:
                print('>%s:%s'%(pdbid,record.id[-1]), file=ff)
                print (record.seq, file=ff)
            elif record.id[-1] in chains:
                print('>%s:%s'%(pdbid,record.id[-1]), file=ff)
                print (record.seq, file=ff)

    ff.close()

def wrap_save_fasta_from_pdb(pdbtxtfile, fastadir=None):
    with open(pdbtxtfile) as ptf:
        ptflines = ptf.read().splitlines()

    for line in ptflines:
        try:
            pdbfile,chains = line.strip().split()
        except ValueError:
            pdbfile = line.strip().split()
            chains=None
        if chains == None:
            suff = '.fasta'
        else:
            suff = '%s.fasta'%chains
        fastaname = os.path.splitext(os.path.basename(pdbfile))[0][:4]
        if fastadir==None:
            fastafile = os.path.join(os.path.dirname(pdbfile),'%s_%s'%(fastaname,suff))
        else:
            fastafile = os.path.join(fastadir,'%s_%s'%(fastaname,suff))
        save_fasta_from_pdb(pdbfile,fastafile,chains=chains)
            

def wrap_write_fasta_file(cases_list_file, output_dest, sep=' '):
    """
    output_dest can be text file with dest for each case
    OR can be a directory where all the fasta files are saved
    """
    with open(cases_list_file) as clf:
        clflines = clf.read().splitlines()
    
    if os.path.isdir(output_dest):
        destdir = output_dest
    else:
        with open(output_dest) as dr:
            drlines = dr.read().splitlines()
    
    for cind, line in enumerate(clflines):
        info = line.strip().split(sep)
        pdbid, chains, seq = info[:3]
        if not os.path.isdir(output_dest):
            destdir = drlines[cind]
        fastafile = os.path.join(destdir,'%s_%s.fasta'%(pdbid,chains))
        if os.path.exists(fastafile): continue
        write_fasta_file([chains], [seq], fastafile, pdbid=pdbid)

def write_fasta_file(chains_list, seq_list, output, pdbid='trgt'):
    """
    writes a fasta file from inputs of 
    @pdbid string of PDB ID
    @chains list of chain names eg. [AB,C,D]
    @seq list of sequences - must be same length as chains
    @output path to output file
    """
    with open(output, 'w') as fp:
        for chain, seq in zip(chains_list, seq_list):
            fp.write(">%s_%s\n"%(pdbid,chain))
            fp.write("%s\n"%seq)


def randomize_prot_seq(prot_seq, seqid=50):
    """
    Given a protein sequence it creates a randomly perturbed
    prot seq with seqID of X% [default is 50%] to given seq
    """
    #random.seed()
    len_seq = len(prot_seq)
    ident_seq_ind = random.sample(range(len_seq), ceil((seqid/100)*len_seq))
    ident_seq_res = [prot_seq[i] for i in ident_seq_ind]
    residue_opts = ['A', 'G', 'I', 'L', 'P', 'V', 'F', 'W', \
                    'Y', 'D', 'E', 'R', 'H', 'K', 'S', 'T', \
                    'C', 'M', 'N', 'Q']
    #change_inds = [i for i in range(len_seq) if i not in ident_seq_ind]
    new_seq = ''
    dbg_seq = ''
    #pdb.set_trace()
    #print (len_seq)
    #print (ident_seq_ind)
    #print (len(ident_seq_ind))
    for ind, res in enumerate(prot_seq):
        if ind in ident_seq_ind:
            new_seq += res
            dbg_seq += res
        else:
            res_opts = [i for i in residue_opts if i!=res]
            new_res = random.choice(res_opts)
            new_seq += new_res
            dbg_seq += '|'
    #print (dbg_seq)
    #print (prot_seq)
    return new_seq

def wrap_randomize_prot_seq(seq_listfile, dest_dir, seqid=50, sep=' '):
    with open(seq_listfile) as sql:
        sqlines = sql.read().splitlines()
    output_file = open(os.path.join(dest_dir, 'randomized_seqs'), 'w')
    for line in sqlines:
        if line.startswith("#"): continue
        #print (line.strip().split(sep))
        if len(line.strip().split(sep)) == 3:
            pdbid,chains,seq = line.strip().split(sep)
        elif len(line.strip().split(sep)) == 4:
            pdbid,chains,seq,seqid = line.strip().split(sep)
        case_dest = os.path.join(dest_dir, '%s_%s'%(pdbid,chains))
        if not os.path.exists(case_dest):
            os.mkdir(case_dest)
        fasta_opt = os.path.join(case_dest, '%s_%s.fasta'%(pdbid,chains))
        newseq = randomize_prot_seq(seq.strip(), seqid=float(seqid))
        write_fasta_file([chains], [newseq], fasta_opt, pdbid=pdbid)
        print (sep.join([pdbid,chains,newseq,str(seqid)]), file = output_file)

    output_file.close() 
