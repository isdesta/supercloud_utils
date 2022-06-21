"""
Israel Desta
Oct 21, 2021
Purpose of this module is to serve as a way to manipulate pdb data
"""
import os, sys
from Bio import SeqIO
from Bio.PDB import PDBParser, Select, MMCIFParser, PDBList
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import PDBIO
import subprocess
import pymol
from pymol import cmd
from prot_seq_tools import write_fasta_file
import pdb
import string
import glob
import requests
import statistics, json

pp = PDBParser()
mp = MMCIFParser()

def wrap_avg_bfac(pdbdir,optfile,keyname):
    pdbfiles = glob.glob(os.path.join(pdbdir,'*.pdb'))
    pdbfiles.sort()
    opt_dict = dict()
    opt_dict[keyname] = dict()
    for pdbf in pdbfiles:
        avg_bfac = average_bfac(pdbf)
        mod = os.path.basename(pdbf)
        opt_dict[keyname][mod] = avg_bfac

    with open(optfile, 'w') as of:
        json.dump(opt_dict, of)


def average_bfac(pdbfile):
    with open(pdbfile) as pf:
        pflines = pf.read().splitlines()
    
    bfac_vals, prev_resnum = [], ''
    for line in pflines:
        if line.startswith('ATOM'):
            resnum, bfac = line[22:26].strip(), line[60:66].strip()
            if resnum != prev_resnum:
                bfac_vals.append(float(bfac))
                prev_resnum = resnum
    return statistics.mean(bfac_vals)

def wrap_categorize_oligstate(txtfile, dest_dir):
    with open(txtfile) as tf:
        tf_lines = tf.read().splitlines()
    
    stoich_divs = ['monomers', 'monomers_wnonpoly', 'homomers', 'homomers_wnonpoly', 'heteromers', 'heteromers_wnonpoly']
    dest_files = []
    for ind,st in enumerate(stoich_divs):
        dt = os.path.join(dest_dir,st)
        dest_files.append(open(dt,'w'))

    for line in tf_lines:
        pdbid = line[:4]
        ostate = categorize_oligstate(pdbid)

        if ostate == 'monomer':
            print (line,file=dest_files[0])
        elif ostate == 'monomers_wnonpoly':
            print (line,file=dest_files[1])
        elif ostate == 'homomer':
            print (line,file=dest_files[2])
        elif ostate == 'homomers_wnonpoly':
            print (line,file=dest_files[3])
        elif ostate == 'heteromer':
            print (line,file=dest_files[4])
        elif ostate == 'heteromers_wnonpoly':
            print (line,file=dest_files[5])
    
    for fs in dest_files:
        fs.close()

def categorize_oligstate(pdbid):
    url = 'https://data.rcsb.org/rest/v1/core/entry/%s'%pdbid
    resp = requests.get(url)
    caseinfo = resp.json()
    #print (caseinfo)
    #print (caseinfo.keys())

    if caseinfo['rcsb_entry_info']['deposited_nonpolymer_entity_instance_count']==0:
        if caseinfo['rcsb_entry_info']['deposited_polymer_entity_instance_count'] == 1:
            return 'monomer'
        elif caseinfo['rcsb_entry_info']['deposited_polymer_entity_instance_count'] > 1 and \
        len(caseinfo['rcsb_entry_container_identifiers']['polymer_entity_ids']) == 1:
            return 'homomer'
        else:
            return 'heteromer'
    else:
        if caseinfo['rcsb_entry_info']['deposited_polymer_entity_instance_count'] == 1:
            return 'monomers_wnonpoly'
        elif caseinfo['rcsb_entry_info']['deposited_polymer_entity_instance_count'] > 1 and \
        len(caseinfo['rcsb_entry_container_identifiers']['polymer_entity_ids']) == 1:
            return 'homomers_wnonpoly'
        else:
            return 'heteromers_wnonpoly'



def clear_insertions(pdbfile, newpdbfile):
    npf = open(newpdbfile, 'w')
    with open(pdbfile) as pf:
        pflines = pf.read().splitlines()
    for line in pflines:
        if not line.startswith('ATOM'): continue
        if line[26] == ' ':
            print (line, file=npf)
    npf.close()

def wrap_clear_insertions(pdbdir, cleared_dest):
    orig_pdbs = glob.glob(os.path.join(pdbdir,'*.pdb'))

    for op in orig_pdbs:
        opname = os.path.basename(op)
        newpdb = os.path.join(cleared_dest, opname)
        clear_insertions(op, newpdb)


def wrap_strip_side_chains(pdbmap_txtfile, keep_cb=False, dest_dir=None):
    with open(pdbmap_txtfile) as ptf:
        ptflines = ptf.read().splitlines()

    for line in ptflines:
        info_line = line.strip().split()
        if len(info_line) == 3:
            pdbf, newpdbf, chains = info_line
            chainlist = list(chains)
        elif len(info_line) == 2:
            pdbf, chains = info_line
            chainlist = list(chains)
            pdbname = os.path.splitext(os.path.basename(pdbf))[0]
            newpdbf = 'bb_%s.pdb'%pdbname
        elif len(info_line) == 1:
            pdbf = line.strip()
            pdbname = os.path.splitext(os.path.basename(pdbf))[0]
            newpdbf = 'bb_%s.pdb'%pdbname
            chainlist = []
        else:
            continue
        if dest_dir==None:
            dest_dir = os.path.abspath(os.path.dirname(pdbf))
        newpdbfile = os.path.join(dest_dir,newpdbf)
        strip_side_chains(pdbf, newpdbfile, chains=chainlist, inc_cb=keep_cb)
            

def strip_side_chains(pdbfile, newpdbfile, chains=[], inc_cb=False):
    npf = open(newpdbfile, 'w')
    with open(pdbfile) as pf:
        pflines = pf.read().splitlines()

    if inc_cb:
        bb_atoms = ['N', 'CA', 'O', 'C', 'CB']
    else:
        bb_atoms = ['N', 'CA', 'O', 'C']
    for line in pflines:
        if not line.startswith('ATOM'): continue
        atomname = line[13:17].strip()
        chainid = line[21]
        #pdb.set_trace()
        if atomname in bb_atoms and (len(chains)==0 or chainid in chains):
            print (line, file=npf)
        if len(chains)!= 0 and chainid not in chains:
            print (line, file=npf)
    npf.close()
    return npf

def save_pdb_to_local(pdbid, destdir, pdbpath=None, cifpath=None, chains=None):
    if pdbpath==None and cifpath==None:
        cifdir = os.getcwd()
        PDBList().retrieve_pdb_file(pdb_code=pdbid, pdir=cifdir)
        cifpath = os.path.join(cifdir,'%s.cif'%pdbid.lower())
        #print (cifpath)
    #pdb.set_trace()
    if cifpath != None:
        struc = mp.get_structure(pdbid, cifpath)
        #destdir = os.path.dirname(cifpath)
    elif pdbpath != None:
        struc = pp.get_structure(pdbid,pdbpath)
        #destdir = os.path.dirname(pdbpath)
    io = PDBIO()
    io.set_structure(struc)
    if chains==None:
        suffix = '.pdb'
    else:
        suffix = '_%s.pdb'%''.join(chains)
    io.save(os.path.join(destdir,'%s%s'%(pdbid,suffix)), select=SelectChains(chains))

def wrap_save_pdb_to_local(pdblist_textfile, destdir):
    with open(pdblist_textfile) as ptf:
        ptflines = ptf.read().splitlines()
    for line in ptflines:
        chains, pdbpath, cifpath = None, None, None
        try:
            pdbid, chains, stpath = line.strip().split()
            if stpath.endswith('cif'):
                cifpath = stpath
            elif stpath.endswith('pdb'):
                pdbpath = stpath
        except ValueError:
            pdbid, chains = line.strip().split()
        except ValueError:
            pdbid = line.strip().split()
        except:
            continue

        save_pdb_to_local(pdbid,destdir,pdbpath=pdbpath,cifpath=cifpath,chains=chains)

def reorder_chains(pdbfile,newpdbfile,chainorder):
    new = open(newpdbfile, 'w')
    pdbid = os.path.splitext(os.path.basename(pdbfile))[0][:4]
    pdbdir = os.path.dirname(pdbfile)
    struc = pp.get_structure(pdbid, pdbfile)
    io = PDBIO()
    io.set_structure(struc)
    for chain in chainorder:
        tempfile = os.path.join(pdbdir,'%s_%s'%(pdbid,chain))
        io.save(tempfile, select=SelectChains(list(chain)))
        with open(tempfile) as tf:
            tflines = tf.readlines()
            if tflines[-1].startswith('END'):
                new.write(''.join(tflines[:-1]))
            else:
                new.write(''.join(tflines))
        os.remove(tempfile)
    new.write('END')
    new.close()

def wrap_reorder_chains(pdblist_txtfile):
    with open(pdblist_txtfile) as ptf:
        ptflines = ptf.read().splitlines()

    for line in ptflines:
        try:
            pdbfile,newpdbfile,chainorder = line.strip().split()
        except ValueError:
            newpdbfile,chainorder = line.strip().split()
            pdbfile = os.path.join(os.getcwd(),os.path.basename(newpdbfile))
        except:
            continue
        reorder_chains(pdbfile, newpdbfile, chainorder)


def pdb_to_cif(pdbfile, cifopt):
    name = os.path.splitext(os.path.basename(pdbfile))[0]
    struc = p.get_structure(name, pdbfile)
    io = MMCIFIO()
    io.set_structure(struc)
    io.save(cifopt)

def wrap_pdb_to_cif(path, dest=None):
    if os.path.isfile(path):
        pdbpaths, destlist = [], []
        with open(path) as txt:
            lines = txt.read().splitlines()
        for ln in lines:
            inf = ln.strip().split()
            if len(inf) > 1:
                pdbpaths.append(inf[0])
                dest = inf[1]
            elif dest == None and len(inf) == 1:
                dest = os.path.abspath(os.path.dirname(inf[0]))
            destlist.append(dest)
    else:
        pdbpaths = glob.glob(os.path.join(path, '*.pdb'))
        if dest==None:
            dest = os.path.abspath(path)
        destlist = [dest]*len(pdbpaths)
    for pdbf,casedest in zip(pdbpaths,destlist):
        name = os.path.splitext(os.path.basename(pdbf))[0]
        cifopt = os.path.join(casedest,'%s.cif'%name)
        pdb_to_cif(pdbf, cifopt)

def rename_pdb(pdbfile, optdir, trnslt_dictionary= None, \
        chain_opt=list(string.ascii_uppercase)):
    pdbname = os.path.basename(pdbfile)[:-4]
    cmd.load(pdbfile, 'trgt')
    if trnslt_dictionary == None:
        curr_chains = cmd.get_chains('trgt')
        for old_chain, new_chain in zip(curr_chains,chain_opt[:len(curr_chains)]):
            cmd.alter("trgt and chain %s"%old_chain, "chain=\'%s\'"%new_chain)
    else:
        for old_chain, new_chain in trnslt_dictionary.items():
            cmd.alter("trgt and chain %s"%old_chain, "chain=\'%s\'"%new_chain)
    #if optdir == None:
    #    optdir = os.path.dirname(pdbfile)
    cmd.save(os.path.join(optdir,'%s_renamed.pdb'%pdbname))
    cmd.delete('*')

def chain_translate(oldchains, newchains):
    if len(oldchains) != len(newchains):
        assert len(oldchains) != len(newchains), "The oldchains and \
                new chains lengths do not match."
    else:
        trnslt_dict = dict()
        for ochain,nchain in zip(oldchains,newchains):
            trnslt_dict[ochain] = nchain

    return trnslt_dict

def wrap_rename_pdb_fromtext(infotext, optpath=None):
    with open(infotext) as itf:
        info_lines = itf.read().splitlines()
    
    for line in info_lines:
        infr = line.strip().split()
        # pdbpath, [oldchains], [corresponding_newchains], [optpath]
        if optpath==None and len(infr) != 4:
            assert optpath==None, "Output path or directory is not given"
        if len(infr) >= 3:
            translate_dict = chain_translate(infr[1], infr[2])
            try:
                pdbpath, oldchains, new_chains, optpath = infr
            except:
                pdbpath, oldchains, new_chains = infr
            rename_pdb(pdbpath, optpath, trnslt_dictionary=translate_dict)
        else:
            try:
                pdbpath, oldchains = infr
                translate_dict = chain_translate(oldchains,\
                        list(string.ascii_uppercase)[:len(oldchains)])
                rename_pdb(pdbpath, optpath, trnslt_dictionary=translate_dict)
            except:
                pdbpath = infr[0]
                rename_pdb(pdbpath, optpath)
        
        
def wrap_rename_pdb(pdbpath, optpath, chain_options=None):
    if os.path.isdir(pdbpath):
        pdbfiles = glob.glob(os.path.join(pdbpath,'*.pdb'))
    else:
        assert "The given path is not a directory."
    if os.path.isdir(optpath):
        optdirs = [optpath]
    else:
        with open(optpath) as of:
            optdirs = of.read().splitlines()

    if len(pdbfiles) > 1 and len(optdirs) == 1:
        optdirs = [optpath] * len(pdbfiles)
    
    if chain_options != None:
        with open(chain_options) as cof:
            coflines = cof.read().splitlines()

        if len(pdbfiles) == len(optdirs) == len(coflines):
            for pdbf, optd, copts in zip(pdbfiles,optdirs,coflines):
                rename_pdb(pdbf, optd, copts)
    else:
        if len(pdbfiles) == len(optdirs):
            for pdbf, optd in zip(pdbfiles,optdirs):
                rename_pdb(pdbf, optd)

def get_seqs_and_chains_from_fastastr(fastastr, chainlist):
    chains, seqs = [], []
    for ch in fastastr.split('>'):
        if ch == '': continue
        name_and_seq = ch.strip().split('\n')
        chain_name = name_and_seq[0][-1:] # just chain name
        seq = "".join(name_and_seq[1:])
        
        if chainlist=='all' and len(seq) > 0:
            chains.append(chain_name)
            seqs.append(seq)
        elif chain_name in chainlist and len(seq) > 0:
            chains.append(chain_name)
            seqs.append(seq)

    return chains, seqs

def get_seqs_and_chains_from_file(pdb_file):
    cmd.load(pdbfile, 'trgt')
    fasta_seqs = cmd.get_fastastr('trgt')
    chains, seqs = get_seqs_and_chains_from_fastastr(fasta_seqs)
    
    cmd.delete('trgt')
    return chains, seqs

def save_pdb_from_id(pdbid_listfile, dest_dirfile, sep=' '):
    """
    @param pdbid_listfile is a text file with PDBIDs
           if chains are specified then only those chains
           will be saved
           eg. 1FW3_AG
    @param dest_dirfile can be a text file with dirs for
           each of the pdbs in the first @param, OR
           it can be a directory in which all the pdbs 
           will be saved
    [@param] sep is the char with which the ID and the chain
             names are separated. Default is space
    """
    with open(pdbid_listfile) as plf:
        plflines = plf.read().splitlines()
    
    if os.path.isdir(dest_dirfile):
        destdir = dest_dirfile
    else:
        with open(dest_dirfile) as dr:
            drlines = dr.read().splitlines()
    cmd.delete('*') 
    for pind, case in enumerate(plflines):
        if not os.path.isdir(dest_dirfile):
            destdir = drlines[pind]
        inf = case.strip().split(sep)
        pdbid = inf[0]
        cmd.fetch(pdbid, 'tarProtein')
        if len(inf) == 1:
            chainList='all'
            output_file = os.path.join(destdir, '%s.pdb'%pdbid)
        else:
            chainList = list(inf[1])
            output_file = os.path.join(destdir, '%s_%s.pdb'%(pdbid,inf[1]))
            all_chains = cmd.get_chains('tarProtein')
            del_chains = [i for i in all_chains if i not in chainList]
            print ("Desired chains are %s"%chainList)
            if len(del_chains) > 0:
                print ("Deleting Chains %s"%','.join(del_chains))
                cmd.select('chain %s'%','.join(del_chains))
                cmd.remove('sele')
        cmd.save(output_file)
        cmd.delete('*')

class SelectChains(Select):
    """
    accept model, chain(s), residue(s) and/or atom(s) specified
    """
    def __init__(self,chainlist=None):
        self.chain_list = chainlist
    
    def accept_chain(self, chain):
        #print (self.chain_list)
        #print (chain.get_id())
        return (chain.get_id() in self.chain_list)

def get_seq_from_pdbid(pdbID, chain_list='all', opt_file=None):
    cmd.fetch(pdbID, 'tar')
    fasta_seqs = cmd.get_fastastr('tar')
    chains, seqs = get_seqs_and_chains_from_fastastr(fasta_seqs, chain_list)

    if opt_file != None:
        write_fasta_file(chains, seqs, opt_file, pdbid=pdbID)
    cmd.delete('tar')
    return chains, seqs

def wrap_get_seq_from_id(pdbid_listfile, dest_dir, sep=' '):
    with open(pdbid_listfile, 'r') as plf:
        plflines = plf.read().splitlines()
    for case in plflines:
        inf = case.strip().split(sep)
        pdbid = inf[0]
        if len(inf) == 1:
            chainList='all'
            output_file = os.path.join(dest_dir, '%s.fasta'%pdbid)
        else:
            chainList = list(inf[1])
            output_file = os.path.join(dest_dir, '%s_%s.fasta'%(pdbid,inf[1]))
        get_seq_from_pdbid(pdbid, chain_list=chainList, opt_file=output_file)
