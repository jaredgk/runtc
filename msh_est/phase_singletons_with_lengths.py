import sys
import argparse
import gzip

def createParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcfname")
    parser.add_argument("left_lengths")
    parser.add_argument("right_lengths")
    parser.add_argument("--mutrate",dest="mutrate",type=float,default=1e-8)
    parser.add_argument("--statfile",dest="statfile",type=str)
    return parser

def findSingleton(la):
    for i in range(9,len(la)):
        if '1' in la[i][0:3]:
            return i
    raise Exception("No singleton found")


def parseChi(f1,f2):
    if ':' in f1:
        fa1 = f1.split(':')
        fa2 = f2.split(':')
        phys1 = (int(fa1[0]) if '*' not in fa1[0] else int(fa1[0][:-1]))
        phys2 = (int(fa2[0]) if '*' not in fa2[0] else int(fa2[0][:-1]))
        gen1 = (float(fa1[1]) if '*' not in fa1[1] else float(fa1[1][:-1]))
        gen2 = (float(fa2[1]) if '*' not in fa2[1] else float(fa2[1][:-1]))
        phys = phys1+phys2
        gen = gen1 + gen2
        return phys, gen
        #chi = float(phys1)*mutrate+float(phys2)*mutrate+gen1+gen2
        #return chi
    phys1 = (int(f1) if '*' not in f1 else int(f1[:-1]))
    phys2 = (int(f2) if '*' not in f2 else int(f2[:-1]))
    return phys1+phys2, None

def computeChi(f1,f2,mutrate):
    phys, gen = parseChi(f1,f2)
    if gen is None:
        return float(phys)*mutrate
    return float(phys)*mutrate + gen

def getStatString(current_left_la,current_right_la,geno):
    if geno[0] != '0':
        physl,genl = parseChi(current_left_la[-2],current_right_la[-2])
        physr,genr = parseChi(current_left_la[-1],current_right_la[-1])
    else:
        physr,genr = parseChi(current_left_la[-2],current_right_la[-2])
        physl,genl = parseChi(current_left_la[-1],current_right_la[-1])
    if genl is None:
        return str(physl)+'\t'+str(physr)
    return str(physl)+'\t'+str(physr)+'\t'+str(genl)+'\t'+str(genr)

def phaseSingleton(orig_geno,current_left_la,current_right_la,mutrate):
    chi_1 = computeChi(current_left_la[-2],current_right_la[-2],mutrate)
    chi_2 = computeChi(current_left_la[-1],current_right_la[-1],mutrate)
    new_geno = ('0|1' if chi_1 > chi_2 else '1|0')
    return new_geno, (new_geno[0] == orig_geno[0])


def phase_with_lengths(sysargs):
    parser = createParser()
    args = parser.parse_args(sysargs)
    if args.vcfname[-3:] == '.gz':
        vcf = gzip.open(args.vcfname,'r')
    else:
        vcf = open(args.vcfname,'r')
    left_decode = False
    right_decode = False
    if args.left_lengths[-3:] == '.gz':
        left_f = gzip.open(args.left_lengths,'r')
        left_decode = True
    else:
        left_f = open(args.left_lengths,'r')
    if args.right_lengths[-3:] == '.gz':
        right_f = gzip.open(args.right_lengths,'r')
        right_decode = True
    else:
        right_f = open(args.right_lengths,'r')

    if args.statfile is not None:
        statf = open(args.statfile,'w')

    if left_decode:
        left_line = left_f.readline().decode()
    else:
        left_line = left_f.readline()
    if right_decode:
        right_line = right_f.readline().decode()
    else:
        right_line = right_f.readline()

    current_left_la = left_line.strip().split()
    current_right_la = right_line.strip().split()
    current_pos = int(current_left_la[0])
    for line in vcf:
        if line[0] == '#':
            sys.stdout.write(line)
            continue
        la = line.strip().split()
        pos = int(la[1])
        if pos != current_pos:
            sys.stdout.write(line)
            continue
        singleton_idx = findSingleton(la)
        for i in range(len(la)):
            if i != singleton_idx:
                sys.stdout.write(la[i])
            else:
                new_geno, changed = phaseSingleton(la[i],current_left_la,current_right_la,args.mutrate)
                sys.stdout.write(new_geno)
                if args.statfile is not None:
                    statf.write(la[1]+'\t'+str(singleton_idx)+'\t'+getStatString(current_left_la,current_right_la,la[i])+'\n')
            if i != len(la)-1:
                sys.stdout.write('\t')
            else:
                sys.stdout.write('\n')
        if left_decode:
            left_line = left_f.readline().decode()
        else:
            left_line = left_f.readline()
        if right_decode:
            right_line = right_f.readline().decode()
        else:
            right_line = right_f.readline()
        current_left_la = left_line.strip().split()
        current_right_la = right_line.strip().split()

        try:
            current_pos = int(current_left_la[0])
        except IndexError:
            current_pos = 0
            continue


if __name__ == '__main__':
    phase_with_lengths(sys.argv[1:])
        
    
