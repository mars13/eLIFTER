import os
import re
import subprocess
import glob
from lib import event_reader
from lib import lifter


#Class definitions
class Microexon:
    def __init__(self, line, w):
        fields = line.rstrip("\n").split("\t")
        self.id = fields[0]
        self.chr = fields[1]
        self.donors = fields[0].split("-")[0].split(":")[-1].split(",")
        self.acceptors = fields[0].split("-")[2].split(":")[0].split(",")
        self.start = int(fields[2])
        self.end = int(fields[3])
        self.window = w
        self.strand = fields[4]
        self.gene = fields[0].split(";")[0]
        self.seq = fields[5]
        self.lchr = ""
        self.lstart = ""
        self.lend = ""
        self.lacceptors = []
        self.ldonors = []

    def print_ME(self):
        printable= ("%s\t%s\t%d\t%d\t%s\n" %(self.gene, self.chr, self.start, self.end, self.seq))
        return(printable)

    def print_LIFT(self):
        printable= ("%s\t%s\t%s\n" %(self.lchr, self.lstart, self.lend))
        return(printable)


    def print_bed_coord(self):
        printables = []
        for d in self.donors:
            printables.append(self.chr + "\t" + str(d) + "\t" + str((int(d)+1)) + "\t" +"COORD;"+ self.chr + ":" + str(d) + "-" + str((int(d)+1)))
        printables.append(self.chr + "\t" + str(self.start) + "\t" + str(self.end) + "\t" + "COORD;" + self.chr + ":" + str(self.start) + "-" + str(sel
f.end))
        for a in self.acceptors:
            printables.append(self.chr + "\t" + str(a) + "\t" + str((int(a) + 1)) + "\t" + "COORD;" + self.chr + ":" + str(a) + "-" + str((int(a)+1)))
        return(printables)


    def search_lifted(self, dict):
        for d in self.donors:
            id = self.chr + ":" + str(d) + "-" + str((int(d)+1))
            try:
                self.ldonors.append(dict[id][1])
            except:
                self.ldonors.append("NA")
        for a in self.acceptors:

            id = self.chr + ":" + str(a) + "-" + str((int(a)+1))
            try:
                self.ldonors.append(dict[id][1])
            except:
                self.ldonors.append("NA")
        id = self.chr + ":" + str(self.start) + "-" + str(self.end)
        try:
            self.lchr = dict[id][0]
            self.lstart = dict[id][1]
            self.lend = dict[id][2]
        except:
            self.lchr = "NA"
            self.lstart = "NA"
            self.lstart = "NA"


class Synt:
    def __init__(self, header, seq1, seq2):
        fields = header.rstrip("\n").split(" ")
        self.start = int(fields[2])
        self.end = int(fields [3])
        self.chr = fields [4]
        self.mstart = int(fields[5])
        self.mend = int(fields[6])
        self.mstrand = fields[7]
        self.aln_len = int(fields[8])
        self.seq1= seq2.rstrip("\n")
        self.seq2= seq1.rstrip("\n")

    def print_synt(self, s,e):
        printable = ("%d\t%d\t%s\t%d\t%d\t%d\n%s\n%s\n" %(self.start, self.end, self.chr, self.mstart, self.mend, self.aln_len, self.seq1[s:e], self.se
q2[s:e]))
        return(printable)

class Pairs:
    def __init__(self, human, mouse):
        # Human and mouse are lists with: [gene, chr, start, end, strand, seq] (for each species)
        self.humanCoord = "chr" + str(human[1]) + ":" + str(human[2]) + "-" + str(human[3])
        self.mouseCoord = "chr" + str(mouse[1]) + ":" + str(mouse[2]) + "-" + str(mouse[3])
        self.humanGene = human[0]
        self.humanChr = human[1]
        self.humanStart = human[2]
        self.humanEnd = human[3]
        self.humanStrand = human[4]
        self.humanSeq = human[5].upper()
        self.mouseGene = mouse[0]
        self.mouseChr = mouse[1]
        self.mouseStart = mouse[2]
        self.mouseEnd = mouse[3]
        self.mouseStrand = mouse[4]
        self.mouseSeq = mouse[5].upper()

    def compare(self, otherPair):
        if self.humanCoord == otherPair.humanCoord and self.mouseCoord == otherPair.mouseCoord:
            return True
        else:
            return False

    def print_pair(self):
        printable = "%s\t%s\t%s\t%s\t%s\t%s" %(self.humanGene, self.mouseGene, self.humanCoord, self.mouseCoord, self.humanSeq, self.mouseSeq)
        return(printable)

###################
# Other functions #
###################


def get_ME(file, window):
    print("Reading microexon information...\n")
    all = {}
    with open(file, "r") as f:
        for l in f:
            fields = l.rstrip("\n").split("\t")
            gene = fields[0].split(";")[0]
            if gene in all:
                all[gene].append(Microexon(l, window))
            else:
                all[gene] = [Microexon(l, window)]
    return all


def buffer_synt(dir):
    print("Buffering synteny info...\n")
    dict = {}
    for file in glob.glob(dir + "/*"):
        chr = file.split(".")[0].split("/")[-1]
        print(chr)
        dict.setdefault(chr, {})
        with open(file, "r") as f:
            for l in f:
                if l.startswith("#"):
                    next(f)
                elif l[0].isdigit():
                    header = l.rstrip("\n")
                    fields = header.split(" ")
                    seq1 = next(f)
                    seq2 = next(f)
                    dict[chr][fields[0]] = Synt(header, seq1, seq2)
    return(dict)


def get_pairs(mex, synt):
    print("Retrieving orthology pairs information\n")
    pairs_list = []
    c=0
    for key,me in mex.items():
        for i in me:
            l = abs(i.end - i.start) + 1
            for n,aln in synt[i.chr].items():
                if i.strand == "+" and aln.mstrand == "+":
                    if i.start >= aln.start and i.end <= aln.end:
                        w = abs(i.start - aln.start) + 1
                        b1 = 0
                        b2 = 0
                        a1 = 0
                        a2 = 0
                        for x in range(0, len(aln.seq1)):
                            if x < w:
                                if aln.seq1[x] == "-":
                                    b1 += 1
                                    w += 1
                                elif aln.seq2[x] == "-":
                                    b2 += 1
                            elif x >= w and x < l:
                                if aln.seq1[x] == "-":
                                    a1 += 1
                                    l += 1
                                elif aln.seq2[x] == "-":
                                    a2 += 1

                        print(w)
                        print(i.strand)
                        print(aln.mstrand)
                        print(i.print_ME())
                        print(aln.print_synt(0,len(aln.seq1) + 1))
                        c += 1
    print(c)



def parse_table(file):
    pairs_list = []
    with open(file, "r") as f:
        next(f)
        for line in f:
            l = line.rstrip("\n").split("\t")
            human = [l[2].split("-")[0], l[4].split(":")[0].lstrip("chr"), l[4].split(":")[1].split("-")[0], l[4].split(":")[1].split("-")[1], "", ""]
            mouse = [l[1].split("-")[0], l[7].split(":")[0].lstrip("chr"), l[7].split(":")[1].split("-")[0], l[7].split(":")[1].split("-")[1], "", ""]
            pairs_list.append(Pairs(human, mouse))
    return pairs_list


def lift(dic):
    f = "tmp_coord.bed"
    if not os.path.isfile(f):
        fh = open(f, "w")
        for gene, me in dic.items():
            for i in me:
                for n in i.print_bed_coord():
                    fh.write("%s\n" %n)
        fh.close()
    if not os.path.isfile("tmp_lifted.bed"):
        subprocess.call(["/home/mreixachs/liftOver", f, "/home/mreixachs/brain_gonzalez/annotation/hg19ToMm10.over.chain.gz", "tmp_lifted.bed","tmp_unl
ifted.bed"])

    sdic = {}
    with(open("tmp_lifted.bed", "r")) as f:
        for l in f:
            line = l.strip("\n").split("\t")
            id = line[-1].lstrip("COORD;")
            sdic[id] = [line[0], line[1], line[2]]
    print(sdic)

    for gene, me in dic.items():
        for i in me:
            i.search_lifted(sdic)
            print(i.print_ME())
            print(i.print_LIFT())



def main():
    #####ME######
    human = get_events(human_ME,window=20)
    lift(human)

    #synt_dict = buffer_synt(synt_dir)
    #synt_pairs = get_pairs(human, synt_dict)




if __name__ == "__main__":
    sp1 = "./files/"
    sp2="./files/"
    # Table of known pairs of orthologous ME
    synt_dir = "/home/mreixachs/brain_gonzalez/annotation/hg19_mm10_synteny"
    output = "human_mouse_ort_blast_ME.txt"


    human_SE = "../general_files/microexons/human_SE_full.tab"
    mouse_SE = "../general_files/microexons/mouse_SE_full.tab"

    # Ensembl Fasta annotation for all genes
    human_fa_all = "human_genes_SE.fa"
    mouse_fa_all = "mouse_genes_SE.fa"
    output_SE = "human_mouse_ort_blast_SE.txt"
    main()

