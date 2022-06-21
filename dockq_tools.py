import pdb
import glob, json
import os, sys
import subprocess
DQPATH = "/home/gridsan/idesta/bin/"

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser=ArgumentParser(description="Wrapper script for DockQ \
            - Quality measure for protein-protein docking models")
    parser.add_argument("model_paths", nargs=1, type=str, help=\
            "text file with list of paths to model files or dirs with \
            model files inside")
    parser.add_argument("native_paths", nargs=1, type=str, help=\
            "text file with list of paths to native files that correspond\
            directly to the model file or dir with models files inside \
            on the same line of the <model_paths>")
    parser.add_argument('res_output',type=str,default=None,help=\
            "output file for the DockQ scores as a csv. it will also \
            save the results into a diff json file in the same dir")
    parser.add_argument('-short',default=False,action='store_true',\
            help='short output')
    parser.add_argument('-fix_numbering',default=False,action='store_true',\
            help='fix number of model based on native using the perl script')
    parser.add_argument('-quiet',default=True,action='store_true',\
            help='keep quiet!')
    parser.add_argument('-useCA',default=False,action='store_true',\
            help='use CA instead of backbone')
    parser.add_argument('-capripeptide',default=False,action='store_true',\
            help='use CAPRI peptide metrics')
    parser.add_argument('-perm1',default=False,action='store_true',\
            help='use all chain1 permutations to find maximum DockQ \
            (number of comparisons is n! = 24, if combined with \
            -perm2 there will be n!*m! combinations')
    parser.add_argument('-perm2',default=False,action='store_true',\
            help='use all chain2 permutations to find maximum DockQ \
            (number of comparisons is n! = 24, if combined with \
            -perm1 there will be n!*m! combinations')
    parser.add_argument('-start_line',type=int,default=0,help=\
            "starting line number of the two given text files")
    parser.add_argument('-end_line',type=int,default=None,help=\
            "ending line number of the two given text files. NOTE:\
            if default (-1) the last line IS INCLUDED")

    args = parser.parse_args()
    DOCKQ_PATH = os.path.join(DQPATH,'./DockQ/DockQ.py')
    """
    NOTE: the native_paths and model_paths text files should contain
    the chain IDs of each case listed in them
    model_chain1 -> type=str,nargs='+', help=pdb chain order to group together partner 1
    model_chain2 -> type=str,nargs='+', help=pdb chain order to group together partner 2 (complement to partner 1 if undef)
    native_chain1 -> type=str,nargs='+', help=pdb chain order to group together from native partner 1
    native_chain2 -> type=str,nargs='+', help=pdb chain order to group together from native partner 2 (complement to partner 1 if undef)
    """
    #print (args.model_paths)
    #print (args.native_paths)
    #pdb.set_trace()
    optf = open(args.res_output, 'w')
    with open(args.model_paths[0]) as mpf, open(args.native_paths[0]) as npf:
        modlines = mpf.read().splitlines()
        natlines = npf.read().splitlines()
    mod_considered = modlines[args.start_line:args.end_line]
    nat_considered = natlines[args.start_line:args.end_line]
    result_json = dict()
    print ("Modelfile, Nativefile, DockQ, Fnat, iRMS, LRMS, Fnonnat", file=optf)
    for mod_dets, nat_dets in zip(mod_considered, nat_considered):
        modpath, model_chain1, model_chain2 = mod_dets.strip().split()
        natpath, native_chain1, native_chain2 = nat_dets.strip().split()

        if os.path.isdir(modpath):
            modfiles = glob.glob(os.path.join(modpath,'*_ranked_*.pdb'))
        else:
            modfiles = [modpath]

        for modfile in modfiles:
            cmd = [DOCKQ_PATH]
            #print (modfile)
            if args.fix_numbering:
                FIX_NUM = os.path.join(DQPATH, './DockQ/scripts/fix_numbering.pl')
                cmdfix = [FIX_NUM]
                cmdfix += [modfile, natpath]
                print (" ".join(cmdfix))
                process = subprocess.Popen(
                        cmdfix, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = process.communicate()
                retcode = process.wait()

                if retcode:
                    raise RuntimeError(
                            'Fix numbering failed:\nstdout:\n%s\n\nstderr:\n%s\n' % (
                                stdout.decode('utf-8'), stderr[:100_000].decode('utf-8')))
                destdir = os.path.dirname(modfile)
                modname = os.path.basename(modfile)
                modfile = os.path.join(destdir,'%s.fixed'%modname)
            cmd += [modfile, natpath]
            cmd += ['-native_chain1'] + list(native_chain1)
            cmd += ['-model_chain1'] + list(model_chain1)
            cmd += ['-native_chain2'] + list(native_chain2)
            cmd += ['-model_chain2'] + list(model_chain2)
            
            if args.quiet:
                cmd += ['-quiet']
            if args.short:
                cmd += ['-short']
            if args.useCA:
                cmd += ['-useCA']
            if args.perm1:
                cmd += ['-perm1']
            if args.perm2:
                cmd += ['-perm2']
            if args.capripeptide:
                cmd += ['-capri_peptide']

            print (" ".join(cmd))
            process = subprocess.Popen(cmd,\
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()
            if stdout.decode('utf-8') == '':
                if args.perm1:
                    cmd.remove('-perm1')
                if args.perm2:
                    cmd.remove('-perm2')
                try:
                    #print ("permutation of chains is not allowed")
                    process = subprocess.Popen(cmd,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    stdout, stderr = process.communicate()
                except:
                    print ("DockQ evaluation of Model %s DID NOT WORK"%modfile)

            print (stdout.decode('utf-8'))
            #pdb.set_trace()
            res = stdout.decode('utf-8').strip().split()
            dockq, fnat, irms, lrms, fnonnat = res[1], res[3], res[5], res[7], res[9]
            modelpath, nativepath = res[10], res[11]
            resdata = [modelpath, nativepath, dockq, fnat, irms, lrms, fnonnat]
            print (','.join([str(i) for i in resdata]), file=optf)
            resdict = dict()
            resdict['native'] = nativepath
            resdict['dockq'] = float(dockq)
            resdict['fnat'] = float(fnat)
            resdict['irms'] = float(irms)
            resdict['lrms'] = float(lrms)
            resdict['fnonnat'] = float(fnonnat)
            resdict['model'] = modelpath
            result_json[os.path.basename(modelpath)] = resdict
    optf.close()
    json_file = os.path.join(os.path.dirname(args.res_output), '%s.json'%args.res_output)
    with open(json_file, 'w') as jf:
        json.dump(result_json,jf)

        

