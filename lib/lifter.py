import subprocess
import os
from lib.event_reader import *

class Pair:
    def __init__(self, id1, id2, m):
        self.id1 = id1
        self.id2 = id2
        self.m = m


def coord2bed(dic, outfile):
    logging.info("Converting coordinates to bed...")
    oh = open(outfile, 'w')
    for k,v in dic.items():
        for i in v:
            for c1,c2 in i.coord.items():
                if type(c2) is list:
                    for elem in c2:
                        oh.write("%s\t%s\t%s\t%s\n" % ("chr" + i.chr, elem, int(elem) + 1, i.id + "_" + c1))
                else:
                    oh.write("%s\t%s\t%s\t%s\n" %("chr"+i.chr, c2, int(c2)+1, i.id+"_"+c1))


def compare_coords(d1, d2, d3, w):
    match = {}
    for k,v in d2.items():
        if (float(d2[k]) - w) <= float(d3[k]) <= (float(d2[k]) + w):
            match[k] = (d1[k],d3[k])
    return(match)


def compare_coords_exon(d1, d2, d3, w):
    match = {}
    for k,v in d2.items():
        if k == "e1" | k == "s2":
            match[k] = ""
            for i in len(d2[k]):
                for j in len(d3[k]):
                    if (float(d2[k][i]) - w) <= float(d3[k][j]) <= (float(d2[k][i]) + w):
                        match[k] = (d1[k][i], d3[k][j])
        else:
            if (float(d2[k]) - w) <= float(d3[k]) <= (float(d2[k]) + w):
                match[k] = ([d1[k]],d3[k])
    return(match)


def assign_coordinates(dic, file):
    logging.info("Reading lifted coodinates...\n")
    with open(file, "r") as f:
        for l in f:
            fields = l.rstrip("\n").split("\t")
            gene = fields[3].split(";")[0]
            key = fields[3].split("_")[1]
            id = fields[3].split("_")[0]
            lchr = fields[0].lstrip("chr")
            lcoord = fields[1]

            for i in dic[gene]:
                if i.id == id:
                    i.get_lifted(lchr, key, lcoord)


def get_pairs(dic1, dic2, w):
    a = []
    for k1,v1 in dic1.items():
        for i1 in v1:
            values = list(i1.liftcoord.values())
            if values.count(float('nan')) < 2:
                #Improve performance by providing a list of 1 to 1 orthologues (?)
                for k2,v2 in dic2.items():
                    for i2 in v2:
                        if (i2.chr == i1.liftchr) & (i2.type == "E"):
                            m = compare_coords(i1.coord,i1.liftcoord, i2.coord, w)
                            if "s1" in m.keys() and "e2" in m.keys():
                                a.append(Pair(i1.id, i2.id, m))
                        elif (i2.chr == i1.liftchr) & (i2.type == i1.type):
                            m = compare_coords(i1.coord,i1.liftcoord, i2.coord, w)
                            if len(m) == len(i1.liftcoord):
                                a.append(Pair(i1.id, i2.id, m))

    return(a)


def print_pairs(list, outfile):
    fh = open(outfile, "w")
    for item in list:
        printable=""
        for k,v in item.m.items():
            print(k,v)
            printable=printable+k+":"+"-".join(v)+";"
        fh.write("%s\t%s\t%s\n" %(item.id1,item.id2,printable))

def lift(dic1, dic2, liftpath, chainz, outdir, w):
    bedfile = outdir + "_event_coord.bed"
    coord2bed(dic1, bedfile)
    logging.info("Lifting coordinates...")
    lifted = outdir + "_lifted.bed"
    unlifted = outdir + "_unlifted.bed"
    if not os.path.isfile(lifted):
        subprocess.call([liftpath, bedfile, chainz, lifted, unlifted])
    assign_coordinates(dic1, lifted)
    pair_list = get_pairs(dic1, dic2, w)
    return(pair_list)






