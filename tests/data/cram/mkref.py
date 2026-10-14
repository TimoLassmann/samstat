import hashlib
u = 'ACGT'*10 + 'TTGCAAGGCC'
seq = u*100
open('ref.fa','w').write('>chr1\n' + '\n'.join(seq[i:i+60] for i in range(0,len(seq),60)) + '\n')
m5 = hashlib.md5(seq.encode()).hexdigest()
open('ref/'+m5,'w').write(seq)
print(m5)
