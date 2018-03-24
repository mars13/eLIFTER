import argparse
from lib.lifter import *

parser = argparse.ArgumentParser(description='Obtains pairs orthologous events, exons or junctions using liftOver')
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='an integer for the accumulator')


def main():
    sp1 = get_events(sp1_events)
    sp2 = get_events(sp2_events)
    pair_list = lift(sp1, sp2, liftpath, chainz, outpx, 0)
    print_pairs(pair_list, outpx+"_lifted_events.txt")
    # for k,v in sp1.items():
    #     for i in v:
    #         print(i.chr)
    #         print(i.coord)
    #         print(i.liftchr)
    #         print(i.liftcoord)

if __name__ == "__main__":
    sp1_events = "/home/marina/lab_projects/eLIFTER/files/human_events.txt"
    sp2_events ="/home/marina/lab_projects/eLIFTER/files/mouse_events.txt"
    chainz = "/home/marina/lab_projects/eLIFTER/files/hg19ToMm10.over.chain.gz"
    liftpath = "/home/marina/lab_projects/eLIFTER/lib/liftOver"
    outpx = "/home/marina/lab_projects/eLIFTER/testing"

    main()
