import logging
import re

class Event:
    def __init__(self, id, gene, type, chr, coord, strand):
        self.id = id
        self.gene = gene
        self.type = type
        self.chr = chr
        self.strand = strand
        self.coord = coord
        self.liftchr = ""
        self.liftcoord= {key: float('nan') for key in self.coord}
        self.paired = {key: float('nan') for key in self.coord}

    def get_lifted(self, chr, key, coord):
        self.liftchr = chr
        if self.type == "E" and (key == "s1" or key == "e2"):
            if self.liftcoord[key] == float('nan'):
                self.liftcoord[key] = [coord]
            else:
                self.liftcoord[key].append(coord)
        else:
            self.liftcoord[key] = coord


def parse_event(line):
    id = line.rstrip("\n")
    f = re.split("\W+", id)
    gene = f[0]
    type = f[1]
    chr = f[2]
    coord = {}

    if type == "SE":
        coord["e1"] = f[3]
        coord["s2"] = f[4]
        coord["e2"] = f[5]
        coord["s3"] = f[6]
        strand = f[7]

    elif type == "MX":
        coord["e1"] = f[3]
        coord["s2"] = f[4]
        coord["e2"] = f[5]
        coord["s4"] = f[6]
        coord["e1"] = f[7]
        coord["s3"] = f[8]
        coord["e3"] = f[9]
        coord["s4"] = f[10]
        strand = f[11]

    elif type == "A5":
        coord["e2"] = f[3]
        coord["s3"] = f[4]
        coord["e1"] = f[5]
        coord["s3"] = f[6]
        strand = f[7]

    elif type == "A3":
        coord["e1"] = f[3]
        coord["s2"] = f[4]
        coord["e1"] = f[5]
        coord["s3"] = f[6]
        strand = f[7]

    elif type == "RI":
        coord["s1"] = f[3]
        coord["e1"] = f[4]
        coord["s2"] = f[5]
        coord["e2"] = f[6]
        strand = f[7]

    elif type == "AF":
        coord["s1"] = f[3]
        coord["e1"] = f[4]
        coord["s3"] = f[5]
        coord["s2"] = f[6]
        coord["e2"] = f[7]
        coord["s3"] = f[8]
        strand = f[9]

    elif type == "AL":
        coord["e1"] = f[3]
        coord["s2"] = f[4]
        coord["e2"] = f[5]
        coord["e1"] = f[6]
        coord["s3"] = f[7]
        coord["e3"] = f[8]
        strand = f[9]


    # EXONS - SE like cases with multiple e1/s3 allowed in order to sum up exon information
    elif type == "E":
        coord["s1"] = id.split(":")[2].split("-")[0].split(",")
        coord["e1"] = id.split(":")[2].split("-")[0].split(":")[0]
        coord["s2"] = id.split(":")[2].split("-")[0].split(":")[1]
        coord["e2"] = id.split(":")[2].split("-")[3].split(",")
        strand = id.split(":")[3]

    #TO DO
    #Junctions - two coordinates
    #elif type == "J":

    return Event(id, gene, type, chr, coord, strand)


def get_events(file):
    logging.info("Reading microexon information...\n")
    all = {}
    with open(file, "r") as f:
        for l in f:
            fields = l.rstrip("\n").split("\t")
            gene = fields[0].split(";")[0]
            event = parse_event(l)
            if gene in all:
                all[gene].append(event)
            else:
                all[gene] = [event]
    return all
