import pysam

old_sam_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/m6A/release_v1/output/mAFiA.reads.sam'
# new_sam_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/mAFiA_v1.1/mAFiA.reads.sam'
new_sam_file = '/home/adrian/Data/TRR319_RMaP/Project_BaseCalling/Adrian/results/mAFiA_v1.1_fixed/mAFiA.reads.sam'


# with pysam.Samfile(old_sam_file, 'r') as old_sam:
#     with pysam.Samfile(new_sam_file, 'r') as new_sam:
#         for old_read in old_sam.fetch():
#             if old_read.query_name!='e729fd60-a0ae-407b-be0b-724b4da49758':
#                 continue
#             old_mod_bases = old_read.modified_bases
#             if (old_mod_bases is not None) and len(old_mod_bases):
#                 for new_read in new_sam.fetch():
#                     if new_read.query_name==old_read.query_name:
#                         new_mod_bases = new_read.modified_bases
#                         print(old_read.query_name, old_mod_bases, new_mod_bases)
#                         break

read_len = {}
new_collection = {}
with pysam.Samfile(new_sam_file, 'r') as new_sam:
    for new_read in new_sam.fetch():
        new_mod_bases = new_read.modified_bases_forward
        if new_mod_bases:
            new_collection[new_read.query_name] = new_mod_bases
            read_len[new_read.query_name] = len(new_read.seq)

old_collection = {}
with pysam.Samfile(old_sam_file, 'r') as old_sam:
    for old_read in old_sam.fetch():
        old_mod_bases = old_read.modified_bases_forward
        if old_mod_bases:
            old_collection[old_read.query_name] = old_mod_bases

for read_id, old_mod_bases in old_collection.items():
    if read_id in new_collection.keys():
        # print(set(list(old_mod_bases.values())[0])==set(list(new_collection[read_id].values())[0]))
        if set(list(old_mod_bases.values())[0])!=set(list(new_collection[read_id].values())[0]):
            print(list(old_mod_bases.values())[0], list(new_collection[read_id].values())[0])

            # old_pos = [tup[0] for tup in list(old_mod_bases.values())[0]]
            # new_pos = [tup[0] for tup in list(new_collection[read_id].values())[0]]
            # corr_pos = [read_len[read_id]-this_pos-1 for this_pos in new_pos]
            # print(set(old_pos).issubset(set(new_pos)))