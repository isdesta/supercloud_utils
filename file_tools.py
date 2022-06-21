"""
Israel Oct 29, 2021
base module for manipulation of files
"""
import os, sys
import json
import glob
import pdb

def create_pdblist_file_in_subdirs(par_dir, dest_dir=None):
    """
    Given a parent directory, this function writes names of pdbfiles 
    in each of its subdirs onto a text file
    @param par_dir is the path to the parent directory with multiple
           subdirs which have pdb files within them
    """
    files_list = os.listdir(par_dir)
    subdirs = [os.path.join(par_dir,i) for i in files_list if os.path.isdir(i)]
    for sub in subdirs:
        #print (sub)
        sub_files = os.listdir(sub)
        pdbfiles = [i for i in sub_files if i.endswith('.pdb')]
        subname = os.path.basename(sub)
        if dest_dir==None:
            opt_file = os.path.join(sub,'%s_pdblist'%subname)
        else:
            opt_file = os.path.join(dest_dir, '%s_pdblist'%subname)
        #print (opt_file)
        opt = open(opt_file, 'w')
        for case in pdbfiles:
            print (os.path.basename(case), file=opt)
        opt.close()

def parse_to_json(text_file, output=None, sep=','):
    with open(text_file) as tf:
        tflines = tf.read().splitlines()
    cases_list = []
    keys = tflines[0].strip().split(sep)
    for line in tflines[1:]:
        line_dets = line.strip().split(sep)
        if len(keys) != len(line_dets):continue
        case_dict = dict()
        for kind,k in enumerate(keys):
            if k == '': continue
            if k == 'pdbid' or k=='chain': # for cases with chain_ids that are digits
                case_dict[k.strip()] = line_dets[kind].strip()
            else:
                try:
                    case_dict[k.strip()] = float(line_dets[kind])
                except:
                    case_dict[k.strip()] = line_dets[kind].strip()
        cases_list.append(case_dict)
    
    if output != None:
        with open(output, 'w') as opt:
            json.dump(cases_list,opt)
    else:
        return cases_list

def merge_jsons(outjson,jsonfilelist):
    """
    merges multiple jsons by combining dictionaries with the same keys
    assumes, values are dictionaries
    """
    print (jsonfilelist)
    jsonlist = []
    for jfl in jsonfilelist:
        with open(jfl) as ipt_json:
            data = json.load(ipt_json)
        jsonlist.append(data)
    comb_dict = jsonlist[0]
    #pdb.set_trace()
    for jl in jsonlist[1:]:
        for key,val in comb_dict.items():
            try:
                newdict = {**jl[key], **comb_dict[key]}
                comb_dict[key] = newdict
            except:
                continue
    
    with open(outjson,'w') as oj:
        json.dump(comb_dict,oj)



def merge_json_to_multjsons(idjson, multijson_dir, dest_json, suffix):
    """
    merges one json that has id of the json files in multijson_dir
    with each of the jsons inside the multijson_dir
    """
    multijson_list = glob.glob(os.path.join(multijson_dir, '*%s.json'%suffix))
    with open(idjson) as ij:
        ij_list = json.load(ij)
    dj_list = []
    #pdb.set_trace()
    ct = 0
    for case in ij_list:
        try:
            ident = case['ids']
        except:
            ident = case['CaseName']
        #ident = case['ids']
        case_json = os.path.join(multijson_dir, '%s%s.json'%(ident,suffix))
        if not os.path.exists(case_json):
            #pdb.set_trace()
            ct += 1
            continue
        with open(case_json) as cj:
            cj_list = json.load(cj)
        for cj_case in cj_list:
            newdict = {**cj_case, **case}
            dj_list.append(newdict)
    print ("cases with no TM results: %d"%ct)
    #pdb.set_trace()
    with open(dest_json, 'w') as dj:
        json.dump(dj_list, dj)

def listjson_to_dictjson(listjson, identifier, optjson):
    """
    convert json file in list format to dict format with an identifier
    """
    with open(listjson) as lj:
        lj_list = json.load(lj)
    lj_dict = dict()
    for case in lj_list:
        #if case['Chain1'] == case['Chain2']: continue
        #pdb.set_trace()
        #print (case[identifier].strip())
        if isinstance(case[identifier], str):
            case_key = os.path.basename(case[identifier].strip())
        else:
            case_key = str(case[identifier])
        lj_dict[case_key] = case

    with open(optjson, 'w') as oj:
        json.dump(lj_dict, oj)

def getjson_from_multiplejsons(jsonlist, optjson, suffix='_ranked', \
        keys= []):
    """
    from two dict jsons, get a new json that combines 
    some keys from them. The ref json should be placed first
    """
    jsninfo = []
    for jsnfile in jsonlist:
        with open(jsnfile) as jf:
            jfdict = json.load(jf)
            jsninfo.append(jfdict)
    optdict = dict()
    ct = 0
    for key,val in jsninfo[0].items():
        #pdb.set_trace()
        key_dets = key.strip().split('_')
        if len(key_dets) < 3:continue
        ct += 1
        key_start = key[:key.find(suffix)]
        #key_start = '_'.join(key_dets[:3])
        #print (key_start)
        
        #pdb.set_trace()
        if key_start not in jsninfo[1] and key not in jsninfo[1]:
            #pdb.set_trace()
            continue
        try:
            sec_val = jsninfo[1][key_start]
        except:
            sec_val = jsninfo[1][key]
        #pdb.set_trace()
        optdict[key] = dict()
        optdict[key]['model_name'] = key
        str_rank = key[key.find(suffix)+8]
        act_rank = int(str_rank) +1
        optdict[key]['rank'] = act_rank
        for kd in keys:
            new_kd = kd
            if kd == 'pLDDT' and kd not in val \
                    and kd not in sec_val:
                kd = 'pLDDT_%d'%act_rank
            if kd == 'pTMipTM' and kd not in val \
                    and kd not in sec_val:
                kd = 'ptmiptm_%d'%act_rank
            if kd in val:
                optdict[key][new_kd] = val[kd]
            elif kd in sec_val:
                optdict[key][new_kd] = sec_val[kd]
        #print (optdict)
        #pdb.set_trace()
    print (ct)
    print (len(jsninfo[0]))
    print (len(jsninfo[1]))
    print (len(optdict))
    with open(optjson, 'w') as oj:
        json.dump(optdict,oj)


if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--nargs', nargs='+')
    nargs = []
    for _, value in parser.parse_args()._get_kwargs():
        if value is not None:
            nargs += value
    '''
    parser.add_argument("--pdblistfile", action='store_true', \
    help = "writes file of pdb names inside each subdir of given\
    parent directory. takes two arguments")
    parser.add_argument("--parsetext", action='store_true', \
    help="parses txtfile into a json file with 1st line as keys")
    parser.add_argument("--merge_json", action='store_true',\
    help="merges 'ids' containing json with multiple json files that \
    have other info of the same key.")
    parser.add_argument("--list_to_dict", action='store_true', \
    help="converts listjson to dictjson using given identifier")
    parser.add_argument("--getjson_from_dictjson", action='store_true',\
    help="makes new json with given keys from two given json in form of\
    dictionaries")
    args = parser.parse_args()
    pdb.set_trace()
    '''
    
    if nargs[0] == 'pdblistfile':
        pardir = nargs[1]
        destdir = None
        if len(nargs) > 2:
            destdir = nargs[2]
        create_pdblist_file_in_subdirs(pardir, dest_dir=destdir)
    elif nargs[0] =='parsetext':
        txtfile = nargs[1]
        opt, sep_str = None, ' '
        if len(nargs) == 3:
            opt, sep_str = nargs[2], ' '
        elif len(nargs) == 4:
            opt, sep_str = nargs[2], nargs[3]
        parse_to_json(txtfile, output=opt, sep=sep_str)
    elif nargs[0] == 'merge_json_to_mj':
        idjson, multijson_dir = nargs[1], nargs[2]
        dest_json, suffix = nargs[3], nargs[4]
        merge_json_to_multjsons(idjson, multijson_dir, dest_json, suffix)
    elif nargs[0] =='list_to_dict':
        lstjson, ident, optjson = nargs[1], nargs[2], nargs[3]
        listjson_to_dictjson(lstjson, ident, optjson)
    elif nargs[0] == 'getjson_from_dictjson':
        refjson,secjson,optjson = nargs[1], nargs[2], nargs[3]
        inp_keys = []
        if len(nargs) > 4:
            for na in nargs[3:]:
                inp_keys.append(na)
        getjson_from_multiplejsons([refjson,secjson],optjson,keys=inp_keys)
    elif nargs[0] == 'merge_jsons':
        outjson = nargs[1]
        json_files = []
        if len(nargs) > 2:
            for na in nargs[2:]:
                json_files.append(na)
        merge_jsons(outjson, json_files)
    else:
        print ("Please choose between the functions you want to use")
        print ("<--pdblistfile> OR <--parsetext> OR <--merge_json>")
        print ("<--list_to_dict> OR <--getjson_from_dictjson>")
