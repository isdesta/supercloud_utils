==============================================================================
   US-align: universal structure alignment of monomeric and complex proteins
   and nucleic acids

   This program was written by Chengxin Zhang at Yang Zhang lab,
   Department of Computational Medicine and Bioinformatics,
   University of Michigan, 100 Washtenaw Ave, Ann Arbor, MI 48109-2218.
   Please report issues to yangzhanglab@umich.edu

   References to cite:
   S Gong, C Zhang, Y Zhang. Bioinformatics, btz282 (2019)
   Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005)

   DISCLAIMER:
     Permission to use, copy, modify, and distribute this program for 
     any purpose, with or without fee, is hereby granted, provided that
     the notices on the head, the reference information, and this
     copyright notice appear in all copies or substantial portions of 
     the Software. It is provided "as is" without express or implied 
     warranty.

   *************** updating history ********************************
   2012/01/24: A C/C++ code of TM-align was constructed by J Yang
   2016/05/21: Several updates of this program were made by J Wu, including
              (1) fixed several compiling bugs
              (2) made I/O of C/C++ version consistent with the Fortran version
              (3) added outputs including full-atom and ligand structures
              (4) added options of '-i', '-I' and '-m'
   2016/05/25: fixed a bug on PDB file reading
   2018/06/04: Several updates were made by C Zhang, including
              (1) Fixed bug in reading PDB files with negative residue index,
                  at the expense of the '-o' option now only being able to
                  output superposed structure instead of full rasmol script.
              (2) Implemented the fTM-align algorithm (by the '-fast' option)
                  as described in R Dong, S Pan, Z Peng, Y Zhang, J Yang
                  (2018) Nucleic acids research. gky430.
              (3) Included option to perform TM-align against a whole 
                  folder of PDB files. A full list of options not available
                  in the Fortran version can be explored by TMalign -h
   2018/07/27: Added the -byresi option for TM-score superposition without
               re-alignment as in TMscore and TMscore -c
   2018/08/07: Added the -dir option
   2018/08/14: Added the -split option
   2018/08/16: Added the -infmt1, -infmt2 options.
               TMalign can now read .gz and .bz2 compressed files.
   2018/10/20: C Zhang and S Gong updated the RNA alignment part of
               the program. Changes include:
              (1) new d0 calculation for RNA.
              (2) secondary structure assignment for RNA.
              (3) automatic detection of molecule type (protein vs RNA).
   2019/01/07: C Zhang added support for PDBx/mmCIF format.
   2019/02/09: Fixed asymmetric alignment bug.
   2019/03/17: Added the -cp option for circular permutation
   2019/03/27: Added the -mirror option for mirror structure alignment
   2019/04/25: The RNA-align algorithm was published by Bioinformatics
   2019/07/24: Fixed bug in displaying matching residues.
               Added GDT and MaxSub to TMscore program.
   2019/08/18: Prevent excessive circular permutation alignment by -cp.
   2020/05/19: Add back rasmol output
   2020/12/12: Fixed bug in double precision coordinate mmcif alignment
   2021/01/07: Fixed bug in TMscore -c
   2021/05/29: Remove unnecessary depedency on malloc.h, which prevent
               compilation on Mac OS
   2021/08/17: Complete implementation of MMalign
   2021/10/03: Support Windows
===============================================================================

=========================
 How to install US-align
=========================
To compile the program in your Linux computer, simply enter

    make

or

    g++ -static -O3 -ffast-math -lm -o USalign USalign.cpp

The '-static' flag should be removed on Mac OS, which does not support
building static executables.

USalign compiled on Linux, Mac OS and Linux Subsystem for Windows (WSL2) on
Windows 10 onwards can read both uncompressed files and gz compressed
files, provided that the "zcat" command is available. On the other hand, due
to the lack of POSIX support on Windows, US-align natively compiled on Windows
without WSL2 cannot parse gz compressed files.

=====================
 How to use US-align
=====================
You can run the program without arguments to obtain a brief instruction

    ./USalign structure1.pdb structure2.pdb

A full list of available options can be explored by:

    ./USalign -h
