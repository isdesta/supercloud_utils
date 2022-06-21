"""
Israel Desta
Nov 8, 2021
Run tmalign program with different functionality
"""
from usalign_tools import *
import pdb

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("tatype", help="enter 'def' or 'all_to_all'\
                        or 'one_to_all' or 'all_to_one'")
    parser.add_argument("--norm_len", default=None, \
    help="Please enter an integer or \'average\'")
    parser.add_argument("--dir1dir2")
    parser.add_argument("--dir1")
    parser.add_argument("--dir2")
    parser.add_argument("--chain1")
    parser.add_argument("--chain2")
    parser.add_argument("--files_list")
    parser.add_argument("--mm_type", type=int, default=0)
    parser.add_argument("--ter_type", type=int, default=2)
    parser.add_argument("--sup", action='store_true')
    parser.add_argument("--wrap", action='store_true')
    parser.add_argument("--wrap_list", default=None)
    parser.add_argument("--wrap_dest", default=None)
    args = parser.parse_args()
    
    choice = args.tatype
    len_choice = args.norm_len
    if choice == 'def':
        lngth_val, avg_val = USAlign().def_options(len_choice)
        fin = USAlign().default(args.chain1, args.chain2, \
              sup=args.sup, length=lngth_val, avg=avg_val)
        algmt = USAlign().parse_res(fin)
        print (fin, "\n")
        print (algmt)
    elif choice == 'all_to_all':
        if args.wrap and args.wrap_list != None:
            with open(args.wrap_list) as wl:
                wlines = wl.read().splitlines()
            for line in wlines:
                files_list = os.path.join(line,"%s_pdblist"%line)
                fin = USAlign().all_to_all(line, files_list, norm_lngth=len_choice)
                res_file = os.path.join(args.wrap_dest,'%s_statsfile.json'%line)
                with open(res_file, 'w') as rf:
                    json.dump(fin, rf)
        elif args.wrap and args.wrap_list == None:
            print ("Please enter a list of subdirs to analyze")
        else:
            fin = USAlign().all_to_all(args.dir1dir2, args.files_list, norm_lngth=len_choice)
            print (fin)
    elif choice == 'one_to_all':
        if args.wrap and args.wrap_list != None:
            #pdb.set_trace()
            with open(args.wrap_list) as wl:
                wlines = wl.read().splitlines()
            #print (len(wlines))
            for ind,line in enumerate(wlines):
                if ind%100 == 0: print (ind, line)
                files_list = os.path.join(line,"%s_pdblist"%line)
                chain1 = os.path.join(line,"%s_renamed.pdb"%line.upper())
                try:
                    fin = USAlign().one_to_all(chain1, line, files_list,\
                            norm_lngth=len_choice, mm_type=args.mm_type, \
                            ter_type=args.ter_type)
                except:
                    continue
                res_file = os.path.join(args.wrap_dest,'%s_statsfile.json'%line)
                with open(res_file, 'w') as rf:
                    json.dump(fin, rf)
        elif args.wrap and args.wrap_list == None:
            print ("Please enter a list of subdirs to analyze")
        else:
            fin = USAlign().one_to_all(args.chain1, args.dir2, args.files_list, norm_lngth=len_choice)
            print (fin)
    elif choice == 'all_to_one':
        print ("All to one")
        if args.wrap and args.wrap_list != None:
            with open(args.wrap_list) as wl:
                wlines = wl.read().splitlines()
            for line in wlines:
                files_list = os.path.join(line,"%s_pdblist"%line)
                chain2 = os.path.join(line,"%s.pdb"%line.upper())
                try:
                    fin = USAlign().all_to_one(line, files_list, chain2, norm_lngth=len_choice)
                except:
                    print ("% failed to be scores by TMalign"%line)
                    continue
                res_file = os.path.join(args.wrap_dest,'%s_statsfile.json'%line)
                with open(res_file, 'w') as rf:
                    json.dump(fin, rf)
        elif args.wrap and args.wrap_list == None:
            print ("Please enter a list of subdirs to analyze")
        else:
            fin = USAlign().all_to_one(args.dir1, args.files_list, args.chain2, norm_lngth=len_choice)
            print (fin)
        
    else:
        print ("Please enter only one of the following options:")
        print ("'def' or 'all_to_all' or 'one_to_all' or 'all_to_one'")

