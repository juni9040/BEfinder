import numpy as np
import pandas as pd
ref_seq = 'CCGCAGTCTACGCTTTGTGTTCCAG'
result = None

for i in range(9,18):
    f = open(f'result/eVLP/S{i}/S{i}.extendedFrags.fastq')
    data = f.readlines()
    seq_list = list()
    num = 0
    AG_conv = 0
    for line in data:
        line = line.strip()
        line = line[4:]
        if line.startswith('CCAGCCCCATCTGTCAAACT') and len(line)==200:
            target=line[79:104]
            target_list = list(target)
            seq_list.append(target_list)
            num += 1


    #BEfinder ver 2.0 : result with all base 
    seq_vector = np.array(seq_list)
    ref_vector = np.array(list(ref_seq))
    print(seq_vector.shape, ref_vector.shape)
    base_array = np.array([[['A']], [['T']], [['G']], [['C']]])
    seq_vector = seq_vector==base_array
    ref_vector = ref_vector==base_array
    mut_vector = (seq_vector!=ref_vector) & seq_vector
    print(seq_vector.shape, ref_vector.shape, mut_vector.shape)
    out = [np.array(list(ref_seq))]
    out = np.concatenate((out, seq_vector.sum(axis=1)),axis=0)
    
    df = pd.DataFrame(out)
    df.index=['sample'+str(i),'A','T','G','C']
    df.to_csv('result/eVLP/'+str(i)+'result.csv', header=True, index=True)


''' ver 1.0
    mut_vector = seq_vector == np.array(list(ref_seq))

    
    mut_list = seq_vector[np.where(mut_vector==True)]
    print(seq_vector, mut_vector, mut_list)
    for base in mut_list:
        if base=='C':AG_conv+=1
    print(i, num, AG_conv, np.sum(mut_vector, axis=0)/num*100)
    output = np.sum(mut_vector, axis=0)/num*100
    if result is None:
        result = output
    else:
        result = np.concatenate((result, output))
print(result.reshape(8,25))
np.savetxt('result1127.csv', result.reshape(8,25), delimiter=',')
'''