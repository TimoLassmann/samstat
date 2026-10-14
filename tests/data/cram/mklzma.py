import random
random.seed(24)
u = 'ACGT'*10 + 'TTGCAAGGCC'
ref = u*100
M5 = 'd1f100f64f7fccbffef69e01113684eb'
out = ['@HD\tVN:1.6\tSO:coordinate', '@SQ\tSN:chr1\tLN:5000\tM5:' + M5]
recs = []
for i in range(300):
    pos = random.randint(0, 5000-100)
    seq = list(ref[pos:pos+100])
    md = ''
    run = 0
    nm = 0
    for j in range(100):
        if random.random() < 0.03:
            alt = random.choice([c for c in 'ACGT' if c != seq[j]])
            md += str(run) + seq[j]
            run = 0
            seq[j] = alt
            nm += 1
        else:
            run += 1
    md += str(run)
    q = ''.join(chr(33 + random.randint(2, 40)) for _ in range(100))
    name = 'r%05d_%s' % (i, ''.join(random.choice('abcdefghijklmnop') for _ in range(24)))
    mapq = random.choice([0, 5, 15, 30, 60])
    recs.append((pos + 1, '\t'.join([name, '0', 'chr1', str(pos+1), str(mapq), '100M', '*', '0', '0', ''.join(seq), q, 'NM:i:%d' % nm, 'MD:Z:' + md])))
recs.sort(key=lambda r: r[0])
open('lz.sam', 'w').write('\n'.join(out + [r[1] for r in recs]) + '\n')
print('reads with MAPQ >= 30:', sum(1 for r in recs if int(r[1].split('\t')[4]) >= 30))

# a long, repeated tag: this is what makes the LZMA codec the smallest choice for a block
s = ''.join(random.choice('ACGTN') for _ in range(3000))
lines = open('lz.sam').read().splitlines()
open('lz2.sam', 'w').write('\n'.join(l if l.startswith('@') else l + '\tXX:Z:' + s for l in lines) + '\n')
