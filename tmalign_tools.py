"""
Israel Desta
Oct 21, 2021
Puspose of this module to serve as an extension of the TMalign executable.
It offers multiple functions to batch run TMalign.
"""
import os, sys
import subprocess
import json
import pdb

class TmAlign:
    def __init__(self, p='/home/gridsan/idesta/bin/utils/tmalign/TMalign'):
        # initialize with just the path to the TMalign executable
        # if user has the executable in a different path, they can assign
        self.tapath = p

    def default(self, chain1, chain2, opt=None, \
                length=None, avg=False, \
                shrt=False, lng=False):
        # default usage of tmalign (requires two chains)
        # also gives users option to normalize over different lengths
        # AND/OR to output the aligned structure into a rasmol file
        cmd = [self.tapath]
        cmd += [chain1, chain2]
        if length != None: # user determined length for tm-align normalising
            cmd += ["-L", length]
        if avg == True: # take avg length of the two chains
            cmd += ["-a"]
        if shrt == True: # take the shorter of the 2 chains
            cmd += ["-b"]
        if lng == True: # take the longer of the 2 chains
            cmd += ["-c"]
        if opt != None: # output the result into a rasmol file
            cmd += ["-o", opt]
        
        print (" ".join(cmd), file = sys.stderr)
        result = subprocess.run(cmd, stdout=subprocess.PIPE, universal_newlines=True, check=True)
        #print (result.stdout)
        return result.stdout 
    
    @staticmethod
    def parse_res(tmalign_opt):
        # static method to parse tmalign output into a dictionary
        tmlines = tmalign_opt.splitlines()
        parsed_tmalign = dict()

        for line in tmlines:
            if line.startswith("Name of Chain_1"):
                name = line.strip().split(':')[1]
                parsed_tmalign['Chain1'] = os.path.abspath(name.strip())
            if line.startswith("Name of Chain_2"):
                name = line.strip().split(':')[1]
                parsed_tmalign['Chain2'] = os.path.abspath(name.strip())
            if line.startswith("Length of Chain_1"):
                lngth = line.strip().split(':')[1]
                parsed_tmalign['length_C1'] = int(lngth.strip().split()[0])
            if line.startswith("Length of Chain_2"):
                lngth = line.strip().split(':')[1]
                parsed_tmalign['length_C2'] = int(lngth.strip().split()[0])
            if line.startswith("Aligned"):
                inf = line.strip().split(',')
                al_lngth = int(inf[0].split('=')[1].strip())
                rmsd = float(inf[1].split('=')[1].strip()) 
                seq_id = float(inf[2].split('=')[2].strip())
                parsed_tmalign['aligned_length'] = al_lngth
                parsed_tmalign['rmsd'] = rmsd
                parsed_tmalign['seq_identity'] = seq_id
            if line.startswith("TM") and line.endswith("Chain_1)"):
                tm_score1 = line.strip().split()[1]
                parsed_tmalign['TM_score1'] = float(tm_score1)
            if line.startswith("TM") and line.endswith("Chain_2)"):
                tm_score2 = line.strip().split()[1]
                parsed_tmalign['TM_score2'] = float(tm_score2)
            if line.startswith("(\":\""):
                al_ind = tmlines.index(line)
        # alignment (3 separate lines) is saved as a list
        parsed_tmalign['alignment'] = tmlines[al_ind+1:al_ind+4]

        return parsed_tmalign
    
    @staticmethod
    def extract_info_from_parsed_data(json_file, info_keys=['rmsd'], opt=sys.stdout):
        """
        This function extracts one or more information from the parsed
        tmalign results which it exptects to be saved in json format,
        i.e a list of dictionaries. 
        by default, it will assume that chain1 is the reference file. 
        outputs a text file with the following fmt
        chain 2, info_key1, info_key2...
        """
        with open(json_file) as jf:
            data = json.load(jf)
        for case in data:
            line = [case['Chain2']]
            for key in info_keys:
                line.append(str(case[key]))
            print (','.join(line), file=opt)

    @staticmethod
    def def_options(norm_lngth):
        lngth_opts = ['short', 'long', 'average']
        lngth_val = None
        shrt_val, lng_val, avg_val = False, False, False
        if norm_lngth==None:
            return lngth_val, shrt_val, lng_val, avg_val
        if norm_lngth=='short':
            shrt_val = True
        elif norm_lngth=='long':
            lng_val = True
        elif norm_lngth=='average':
            avg_val = True
        elif type(norm_lngth)==int:
            lngth_val = norm_lngth
        elif type(norm_lngth)!=int or norm_lngth not in lngth_opts:
            return "Please enter an integer or \'short\', \'long\' or \'average\'"
        
        return lngth_val, shrt_val, lng_val, avg_val

    def all_to_all(self, direc, files_list, norm_lngth=None):
        """
        given a directory and the list of pdb files inside
        this function performs an all-to-all tmalignment
        and outputs the results for each
        NOTE:list should have full name of file including pdb ext
        NOTE: norm_lngth can be an integer, or 'short', 'long', 'average'
        """
        lngth_val, shrt_val, lng_val, avg_val = self.def_options(norm_lngth)
        all_to_all_list = []
        with open(files_list, 'r') as fl:
            flines = fl.read().splitlines()
        for c1 in flines:
            c1_file = os.path.join(direc, c1)
            for c2 in flines:
                c2_file = os.path.join(direc, c2)
                res = self.default(c1_file, c2_file, length=lngth_val, \
                shrt=shrt_val, lng = lng_val, avg = avg_val)
                all_to_all_list.append(self.parse_res(res))
        
        return all_to_all_list

    def one_to_all(self, chain1, direc2, file2_list, norm_lngth=None):
        """
        given a directory of pdb files, their list and a reference pdb file
        this function performs a one-to-all tmalignment
        and outputs the results for each
        NOTE:list should have full name of file including pdb ext
        NOTE: norm_lngth can be an integer, or 'short', 'long', 'average'
        """
        lngth_val, shrt_val, lng_val, avg_val = self.def_options(norm_lngth)
        one_to_all_list = []
        with open(file2_list, 'r') as fl:
            flines = fl.read().splitlines()
        for c2 in flines:
            c2_file = os.path.join(direc2, c2)
            res = self.default(chain1, c2_file, length=lngth_val, \
            shrt=shrt_val, lng = lng_val, avg = avg_val)
            one_to_all_list.append(self.parse_res(res))
        
        return one_to_all_list

    def all_to_one(self, direc1, file1_list, chain2, norm_lngth=None):
        """
        given a directory of pdb files, their list and a reference pdb file
        this function performs a all-to-one tmalignment
        and outputs the results for each
        NOTE:list should have full name of file including pdb ext
        NOTE: norm_lngth can be an integer, or 'short', 'long', 'average'
        """
        lngth_val, shrt_val, lng_val, avg_val = self.def_options(norm_lngth)
        all_to_one_list = []
        with open(file1_list, 'r') as fl:
            flines = fl.read().splitlines()
        for c1 in flines:
            c1_file = os.path.join(direc1, c1)
            #pdb.set_trace()
            res = self.default(chain2, c1_file, length=lngth_val, \
            shrt=shrt_val, lng = lng_val, avg = avg_val)
            all_to_one_list.append(self.parse_res(res))
        
        return all_to_one_list

