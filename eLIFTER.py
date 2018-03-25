import argparse
from lib.lifter import *

parser = argparse.ArgumentParser(description='Obtains pairs orthologous events, exons or junctions using liftOver')
parser.add_argument("-1","--sp1", type=str, required=True,
                    help='events id for species 1')
parser.add_argument("-2","--sp2", type=str, required=True,
                    help='events id for species 2')
parser.add_argument("-l", "--liftover", default="./lib/liftOver",
                    help='liftover path')
parser.add_argument("-c", "--chainz", default="./files/hg19ToMm10.over.chain.gz",
                    help='chainz file path')
parser.add_argument("-t", "--tag", default="",
                    help='tag for exon or junction analysis if necessary')
parser.add_argument("-o", "--outpref", required=True,
                    help='outfile path + prefix')

args = parser.parse_args()


def main():
    sp1 = get_events(sp1_events, tag)
    sp2 = get_events(sp2_events, tag)
    pair_list = lift(sp1, sp2, tag, liftpath, chainz, outpx, 0)
    print_pairs(pair_list, outpx+"_lifted_events.txt")
    # for k,v in sp1.items():
    #     for i in v:
    #         print(i.chr)
    #         print(i.coord)
    #         print(i.liftchr)
    #         print(i.liftcoord)

if __name__ == "__main__":
    #sp1_events = "/home/marina/lab_projects/eLIFTER/Homo_sapiens.GRCh37.85.CDS.NO_PSEUDOGENE.CLEAN.events_SE_strict.merged.txt"
    #sp2_events ="/home/marina/lab_projects/eLIFTER/Mus_musculus.GRCm38.85.CDS.NO_PSEUDOGENE.CLEAN.events_SE_strict.merged.txt"
    #chainz = "/home/marina/lab_projects/eLIFTER/files/hg19ToMm10.over.chain.gz"
    #liftpath = "/home/marina/lab_projects/eLIFTER/lib/liftOver"
    #outpx = "/home/marina/lab_projects/eLIFTER/testing_me"
    #tag="E"

    sp1=args.sp1
    sp2=args.sp2
    chainz=args.chainz
    liftpath=args.liftover
    tag=args.tag
    outpx=args.outpref

    main()
