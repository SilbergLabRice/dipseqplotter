#######################
## Dependendencies ####
#######################
import re, sys
import UTIL.UTIL as util
import Bio.SeqIO
import numpy as np
import pandas as pd
from matplotlib import pyplot as pp
from matplotlib import rc

#######################
### Set Plot Style ####
#######################

font = 'Helvetica'
rc('font',**{'family':'sans-serif','sans-serif':[font]})
rc('ytick',**{'major.size': 5})
rc('ytick',**{'minor.size': 3.5})
rc('xtick',**{'major.size': 5})
rc('ytick',**{'labelsize': 10})
rc('xtick',**{'labelsize': 10})
rc('axes',**{'labelsize': 12})
rc('svg', **{'fonttype': 'none'})
rc('mathtext', **{'default' : 'regular'}) #change latex font to regular font

#######################
## Data Input/Output ##
#######################
files = ['T17_dipseq_Tile1_R1.csv', 'T19_dipseq_Tile1_R2.csv','T20_dipseq_Tile2_R1.csv','T21_dipseq_Tile2_R2.csv','T22_dipseq_Tile3_R1.csv','T23_dipseq_Tile3_R2.csv', 'T24_dipseq_Tile4_R1.csv', 'T25_dipseq_Tile4_R2.csv']

counts = pd.DataFrame(index = range(1,171), columns = ['counts'])
counts['counts'] = [0]*170
print(counts)

for i in files:
    insertions_file = i
    insertions = pd.read_csv(insertions_file)
    insertions['residue_aa_idx'] = (insertions['insertion_site']/3)-1

    print('\n')
    print('There are %s total reads with 5p or 3p fixed_ends'  % len(insertions))

    tag_5p = insertions.loc[insertions['fixed_seq_end'] == '5p']
    tag_3p = insertions.loc[insertions['fixed_seq_end'] == '3p']
    print('There are %s total reads with 5p fixed_ends'  % len(tag_5p))
    print('There are %s total reads with 3p fixed_ends'  % len(tag_3p))

    print('\n')
    print('There are %s reads with 5p and 3p fixed_ends'  % len(tag_5p.loc[tag_5p['read_id'].isin(tag_3p['read_id'])]))
    print('There are %s reads with 5p but no 3p fixed_ends'  % len(tag_5p.loc[~tag_5p['read_id'].isin(tag_3p['read_id'])]))
    print('There are %s reads with 3p but no 5p fixed_ends'  % len(tag_3p.loc[~tag_3p['read_id'].isin(tag_5p['read_id'])]))

    print('\n')
    matched_pairs = tag_5p.merge(tag_3p, how = 'inner', on = ['read_id', 'insertion_idx'])
    nonmatched_pairs = tag_5p.merge(tag_3p, how = 'outer', on = ['read_id', 'insertion_idx'])
    print('There are %s reads with matching insertion_idx calls between the 5p and 3p fixed_ends'  % len(matched_pairs))
    print('There are %s reads with non-matching insertion_idx calls between the 5p and 3p fixed_ends'  % len(nonmatched_pairs))

    print('\n')

    ##Filter reads
    inframe_forward = matched_pairs.loc[(matched_pairs['forward_insertion_x']== True) & (matched_pairs['in_frame_insertion_x'] == True)]
    print('There are %s matched pairs with inframe insertions'  % len(inframe_forward))


    #insert is based on the dipseq insert file and length is a parameter in __init__.py
    insert_length = inframe_forward['insert_end_idx_x'] - inframe_forward['insert_start_idx_x'] 
    inframe_forward_inserted = inframe_forward.loc[insert_length == 11]
    print('There are %s matched pairs with correct length insertions'  % len(inframe_forward_inserted))


    #Filter for sequences without a linker 
    linker_length = inframe_forward['linker_end_idx_x'] - inframe_forward['linker_start_idx_x'] 
    inframe_forward_inserted_nolink = inframe_forward_inserted.loc[linker_length == 0]
    print('There are %s matched pairs with correct length insertions with no linker'  % len(inframe_forward_inserted_nolink))




    filtered = inframe_forward_inserted_nolink
    filtered.to_csv('filtered.csv')
    # print(filtered.loc[filtered['residue_aa_idx_x'] > 44].to_string())
    #count_residues

    count = filtered['residue_aa_idx_x'].value_counts().sort_index()

    print(count)
    # counts.add(count, fill_value=0)
    counts = pd.concat([count, counts], axis = 1).sum(axis = 1)
    # print(count.values.tolist())
    # counts['counts'] = counts['counts'] + count.values.tolist()

print(counts)

T1_res_start = 1
T1_res_end = 44
T2_res_start = 45
T2_res_end = 86
T3_res_start = 87
T3_res_end = 128 
T4_res_start = 129
T4_res_end = 170 


T1_Fwd_Primer_Res = 1
T1_Rev_Primer_Res = 88
T2_Fwd_Primer_Res = 26
T2_Rev_Primer_Res = 129

print((T2_Rev_Primer_Res - T2_res_start)*3)
#######################
###### Plotting #######
#######################
# define plot window size 
figWidth = 12
figHeight = 3
fig = pp.figure(figsize=(figWidth, figHeight))
ax = pp.subplot(1,2,1)


# ax.axvline(T4_res_start, ls = ':', color = (0,0,1,0.8))
# ax.axvline(T4_res_end, ls = ':', color = (0,0,1,0.8))
# ax.axvline(T2_Fwd_Primer_Res, color = (1,0,0,1))
# ax.axvline(T2_Rev_Primer_Res, color = (1,0,0,1))

ax.bar(counts.index, counts.values, width = 1, color = [0,0,0,1])

ax.set_title('Filtered reads (n= %i)' % (np.sum(counts.values)), color = [0,0,0,1])

ax.set_ylabel('Read Count')
ax.set_xlabel('Residue')

ax.set_xlim(0,170)

ax.set_yscale('log')
ax.set_ylim(10**0,10**4)

pp.savefig('insertions_tmp.png', bbox_inches="tight")

