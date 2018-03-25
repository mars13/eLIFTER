import logging
import re


class Event:
    def __init__(self, id, gene, type, chr, coord, strand, tag):
        self.id = id
        self.gene = gene
        self.type = type
        self.chr = chr
        self.strand = strand
        self.coord = coord
        self.liftchr = ""
        self.liftcoord= {key: float(0) for key in self.coord}
        self.paired = {key: float(0) for key in self.coord}

    def get_lifted(self, chr, key, coord, tag):
        self.liftchr = chr
        if tag == "E" and (key == "e1" or key == "s2"):
            if self.liftcoord[key] == float(0):
                self.liftcoord[key] = [coord]
            else:
                self.liftcoord[key].append(coord)
        else:
            self.liftcoord[key] = coord


def parse_event(line, tag):
    id = line.rstrip("\n")
    f = re.split(";|:", id)
    gene = f[0]
    type = f[1]
    chr = f[2]
    coord = {}

    if type == "SE":
        coord["s1"] = f[3].split("-")[1]
        coord["e2"] = f[4].split("-")[0]

        if tag == "E":
            coord["e1"] = f[3].split("-")[0].split(",")
            coord["s2"] = f[4].split("-")[1].split(",")
        else:
            coord["e1"] = f[3].split("-")[0]
            coord["s2"] = f[4].split("-")[1]
        strand = f[5]

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
    #TO DO
    #Junctions - two coordinates
    #elif type == "J":

    return Event(id, gene, type, chr, coord, strand, tag)


def get_events(file, tag):
    logging.info("Reading microexon information...\n")
    all = {}
    with open(file, "r") as f:
        for l in f:
            fields = l.rstrip("\n").split("\t")
            gene = fields[0].split(";")[0]
            event = parse_event(l, tag)
            if gene in all:
                all[gene].append(event)
            else:
                all[gene] = [event]
    return all
