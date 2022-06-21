import pdb
import json
import glob
import os, sys

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser=ArgumentParser(description="parses dTERMen output file \
            to extract original seq, dTERMen seq , recovery & energy")
    parser.add_argument("dTERMen_res_path", type=str, nargs='+', \
            help="text file or path(s) with *.out dtermen res files inside")
    parser.add_argument("output_dir", type=str, \
            help="saves json and csv file to the dir specified")
    parser.add_argument("--output_name", type=str, \
            default="parsed_dTERMen")
    args = parser.parse_args()

    id_idx = 2
    # 3AV9_AB_Y_24_1_dTERMen.18425506
    # for uniq_ids like this, extract just the "3AV9_AB_Y" comp using id_idx=3

    drp_inp = args.dTERMen_res_path
    res_files = []
    if len(drp_inp) == 1 and os.path.isfile(drp_inp[0]):
        with open(drp_inp[0]) as dif:
            res_files = dif.read().splitlines()
        # each line contains full paths to dTERMen .out files
    elif len(drp_inp) > 1 and os.path.isfile(drp_inp[0]):
        print ("Please enter a text file with the paths to all res files\
                or one/multiple dirs with .out files inside them")
    else:
        for given_path in drp_inp:
            cur_resf = glob.glob(os.path.join(given_path,'*.out'))
            if len(cur_resf) == 0:
                cur_resf = glob.glob(os.path.join(given_path,'*.txt'))
            res_files.extend(cur_resf)
    #pdb.set_trace()
    opt_csv = open("%s.csv"%os.path.join(args.output_dir,\
            args.output_name), 'w')
    print ("case_ID,dtermen_seq,real_seq,dtermen_rec,avg_engy", \
            file=opt_csv)
    res_json, list_json = dict(), []
    for res in res_files:
        job_fin = False
        dirname = os.path.basename(os.path.dirname(res))
        resname = os.path.splitext(os.path.basename(res))[0]
        #uniq_id = "%s_%s"%(dirname,resname)
        
        with open(res) as rsf:
            rsf_lines = rsf.read().splitlines()
        if len(rsf_lines) < 1:
            print (res)
            continue
        caseid = rsf_lines[-1]
        for line in rsf_lines:
            if line.endswith('lowest-energy sequence'):
                dtermen_seq = line.strip().split()[0]
            elif line.endswith('original sequence'):
                real_seq = line.strip().split()[0]
            elif line.endswith('recovery'):
                dtermen_rec = float(line.strip().split()[0][:-1])
            elif line.startswith('mean energy'):
                avg_engy = float(line.strip().split()[-1])
            elif line.endswith('elapsed'):
                job_fin = True
        uniq_id = "%s_%s"%(caseid,resname)
        if not job_fin: continue
        cased = dict()
        cased['dtermen_seq'] = dtermen_seq
        cased['real_seq'] = real_seq
        cased['dtermen_recovery'] =  dtermen_rec
        cased['avg_energy'] = avg_engy
        res_json[uniq_id] = cased
        cased['ids'] = '_'.join(uniq_id.split('_')[:id_idx])
        list_json.append(cased)
        info_line = [uniq_id, dtermen_seq, real_seq, dtermen_rec, \
                avg_engy]
        print (",".join([str(i) for i in info_line]), file=opt_csv)

    opt_csv.close()
    with open("%sdict.json"%os.path.join(args.output_dir, \
            args.output_name), 'w') as oj:
        json.dump(res_json,oj)
    with open("%slist.json"%os.path.join(args.output_dir, \
            args.output_name), 'w') as oj:
        json.dump(list_json,oj)
