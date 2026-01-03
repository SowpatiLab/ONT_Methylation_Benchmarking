from Bio import SeqIO
import re, argparse

def get_args():
    parser = argparse.ArgumentParser(description=(f"""This script can be used to aggregate rockfish readwise output
    to generate the site-wise methylation output
    """
),formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input-readwise', type=str, help='f5c readwise file', required=True)
    parser.add_argument('-r', '--reference', type=str, help='reference fasta file', required=True)
    parser.add_argument('-o', '--output', type=str, help='output file to write to', required=True)
    
    return parser.parse_args()

if __name__=="__main__":
    args = get_args()
    input_file, reference, out_file = args.input_readwise, args.reference, args.output

    fasta = SeqIO.index(reference,'fasta')
    print('reference loaded!')
    INCREMENTS = { '-': 1, '+': 0, }
    def getRefMotif(fsta, chrom, p1, strand):
        incr1, incr2 = INCREMENTS[strand], INCREMENTS[strand]+1
        if strand=="+": return fsta[chrom].seq[p1+incr1-8:p1+incr2+8]
        else: return fsta[chrom].seq[p1+incr1-(8):p1+incr2+8].reverse_complement()

    with open(input_file) as f:
        with open(out_file, 'w') as fw:
            header = next(f)
            fw.write(header)
            
            for line in f:
                cols = line.split('\t')
                p1, p2 = INCREMENTS[cols[1]], 1+INCREMENTS[cols[1]]
                if int(cols[9])==1: 
                    fw.write(f"{cols[0]}\t{cols[1]}\t{int(cols[2])+p1}\t{int(cols[2])+p2}" + '\t' + '\t'.join(cols[4:]))
                else:
                    splits = [ int(cols[2])-8+i.start() for i in re.finditer('CG', cols[10])]
                    for split in splits:
                        fw.write(f"{cols[0]}\t{cols[1]}\t{split+p1}\t{split+p2}" + '\t' + '\t'.join(cols[4:-2]) + '\t' + f'1\t{getRefMotif(fasta, cols[0], split, cols[1])}' + '\n')
