#!/usr/bin/env python

import os
import sys

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "usage: automate_guidemaker_cfg_files.py <cladeslist> <examplefile>"
        sys.exit(0)
    
    listfile = open(sys.argv[1],"r")
    examplefile = open(sys.argv[2],"r")
    exampledata = [l.strip() for l in examplefile]

    outdir = "guidefiles/"
    try:
        os.mkdir(outdir)
    except OSError:
        pass
    
    for clades in [c.strip() for c in listfile]:
        print clades
        outfile = open(outdir + "_".join([name.strip() for name in clades.split(",")])+".guide","w")
        for line in exampledata:
            outfile.write(line.replace("{}",clades)+"\n")
        outfile.close()
