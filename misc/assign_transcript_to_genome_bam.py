import pandas as pd
import numpy as np
from Bio import SeqIO
import pysam
from tqdm import tqdm
import os

workspace = '/home/darthvader/Data/TRR319_RMaP/Project_B01/Adrian/JK_HEK293_DMSO_1_2_RTA'
in_mod_file = os.path.join(workspace, '_mAFiA/mAFiA.sites.bed')
genome_bam_file = os.path.join(workspace, '_mAFiA/mAFiA.reads.bam')
transcriptome_ref_file = '/home/darthvader/Data/GRCh38_102/GRCh38.cdna.all.fa'
transcriptome_bam_file = os.path.join(workspace, 'transcriptome_mapped.bam')
out_mod_file = os.path.join(workspace, '_mAFiA/mAFiA.sites.bed.annotated.new')

df_mod = pd.read_csv(in_mod_file, sep='\t')

ref = {}
for record in SeqIO.parse(transcriptome_ref_file, 'fasta'):
    ref[record.id] = str(record.seq)

genome_bam = pysam.AlignmentFile(genome_bam_file, 'rb')
genome_index = pysam.IndexedReads(genome_bam)
genome_index.build()

transcriptome_bam = pysam.AlignmentFile(transcriptome_bam_file, 'rb')
transcriptome_index = pysam.IndexedReads(transcriptome_bam)
transcriptome_index.build()


def get_transcript_modRatios(modProbs, thresh=127):
    numReads_modRatios = {}
    for xscript, list_modProbs in modProbs.items():
        numReads_modRatios[xscript] = (len(list_modProbs), int(np.mean(np.array(list_modProbs)>thresh)*100.0))
    return numReads_modRatios


min_coverage = 10
min_splitting = 20


sites_collected = 0
transcripts_collected = 0
for _, row in tqdm(df_mod.iterrows()):
    chrom = row['chrom']
    chromStart = row['chromStart']
    strand = row['strand']
    ref5mer = row['ref5mer']

    transcript_modProbs = {}
    for pileupcolumn in genome_bam.pileup(chrom, chromStart, chromStart+1, truncate=True):
        if pileupcolumn.reference_pos == chromStart:
            for pileupread in pileupcolumn.pileups:
                read = pileupread.alignment
                if read.flag not in [0, 16]:
                    continue
                read_id = read.query_name
                mod_bases = list(read.modified_bases_forward.values())
                if len(mod_bases)==0:
                    continue
                mod_pos, mod_prob = mod_bases[0][0]
                pred_motif = read.get_forward_sequence()[mod_pos-2:mod_pos+3]

                for this_iter in transcriptome_index.find(read_id):
                    if this_iter.flag!=0:
                        continue
                    ref_seq = ref[this_iter.reference_name]
                    this_read_pos = mod_pos
                    ref_shift = 0
                    this_cigar = this_iter.cigar.copy()
                    while this_read_pos>0:
                        next_cigar_tuple = this_cigar.pop(0)
                        remaining = min(next_cigar_tuple[1], this_read_pos)
                        if next_cigar_tuple[0]==0:   # match
                            this_read_pos -= remaining
                            ref_shift += remaining
                        elif next_cigar_tuple[0]==1:   # insertion
                            this_read_pos -= remaining
                        elif next_cigar_tuple[0]==2:   # deletion
                            ref_shift += next_cigar_tuple[1]
                        elif next_cigar_tuple[0]==4:   # soft-clip
                            this_read_pos -= remaining
                    this_ref_pos = this_iter.reference_start + ref_shift

                    # print(this_iter.flag, ref5mer, ref_seq[this_ref_pos-3:this_ref_pos+4])

                    if (this_ref_pos < this_iter.reference_end) and (ref5mer in ref_seq[this_ref_pos - 3:this_ref_pos + 4]):
                        transcript = this_iter.reference_name
                        if transcript not in transcript_modProbs.keys():
                            transcript_modProbs[transcript] = [mod_prob]
                        else:
                            transcript_modProbs[transcript].append(mod_prob)
                        break

    transcript_modProbs = {k: v for (k, v) in transcript_modProbs.items() if len(v)>=min_coverage}
    if (len(transcript_modProbs)<2):
        continue

    transcript_modRatios = get_transcript_modRatios(transcript_modProbs)
    all_modRatios = [val[1] for val in transcript_modRatios.values()]
    if (np.max(all_modRatios) - np.min(all_modRatios)) < min_splitting:
        continue

    print(row)
    print(transcript_modRatios)

    df_out = pd.DataFrame()
    transcriptName = []
    transcriptCoverage = []
    transcriptModRatio = []
    for this_transcript, (this_coverage, this_modRatio) in transcript_modRatios.items():
        transcriptName.append(this_transcript)
        transcriptCoverage.append(this_coverage)
        transcriptModRatio.append(this_modRatio)
        transcripts_collected += 1

    indMax = np.argmax(transcriptModRatio)
    indMin = np.argmin(transcriptModRatio)

    deltaModRatio = transcriptModRatio[indMax] - transcriptModRatio[indMin]
    deltaCoverage = int(np.abs((transcriptCoverage[indMax] - transcriptCoverage[indMin])) / (transcriptCoverage[indMax] + transcriptCoverage[indMin]) * 100.0)

    new_row = row.copy()
    new_row['transcriptName'] = ','.join(transcriptName)
    new_row['transcriptCoverage'] = ','.join([str(cov) for cov in transcriptCoverage])
    new_row['transcriptModRatio'] = ','.join([str(mr) for mr in transcriptModRatio])
    new_row['deltaModRatio'] = deltaModRatio
    new_row['deltaCoverage'] = deltaCoverage
    df_out = pd.concat([df_out, pd.DataFrame(new_row).T])

    if os.path.exists(out_mod_file):
        df_out.to_csv(out_mod_file, sep='\t', index=False, header=False, mode='a')
    else:
        df_out.to_csv(out_mod_file, sep='\t', index=False, mode='w')
    sites_collected += 1
print(f'Collected {sites_collected} sites, {transcripts_collected} transcripts')

genome_bam.close()
transcriptome_bam.close()