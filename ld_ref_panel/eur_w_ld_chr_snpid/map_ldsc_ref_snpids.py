#!/usr/bin/env python3

import numpy as np
import pandas as pd
from pyfaidx import Fasta

fa = Fasta('/hpfs/userws/chiouj02/resources/reference_genomes/GRCh37.fa')
snplist = pd.read_table('w_hm3.snplist', sep='\t', header=0, index_col=0)
#ukb = pd.read_table('/hpfs/projects/ukbiobank/genotypes/UKB_500k_imputed/ukb_imp.snps', sep='\t', header=0)
#ukb = ukb.drop_duplicates(subset='rsid')
#ukb.index = ukb['rsid']
ldscs = pd.DataFrame()
#for i in range(1,23):
#    ldsc = pd.read_table(f'{i}.l2.ldscore.gz', sep='\t', header=0)
#    ldscs = pd.concat([ldscs, ldsc], ignore_index=True)
#    snps = ldsc['SNP'].map(ukb['GRCh37'])
#    ldsc.loc[snps.notnull(), 'SNP'] = snps[snps.notnull()]
#    snplist_chr = snplist.loc[snplist.index.isin(ldsc['SNP'])]
#    snplist_chr['CHR'] = snplist_chr.index.map(ldsc.set_index('SNP')['CHR'])
#    snplist_chr['BP'] = snplist_chr.index.map(ldsc.set_index('SNP')['BP'])
#    snplist_chr['CHR'] = snplist_chr.index.map(ldsc.set_index('SNP')['CHR'])
#    snplist_chr['REF'] = snplist_chr.apply(lambda x: (fa[f'chr{i}'][x['BP']-1]).seq.upper(), axis=1)
#    snplist_chr['ALT'] = np.where(snplist_chr['A1']==snplist_chr['REF'], snplist_chr['A2'], snplist_chr['A1'])
#    snplist_chr['SNPID'] = snplist_chr[['CHR','BP','REF','ALT']].astype(str).apply(lambda x: f'{x["CHR"]}:{x["BP"]}:{x["REF"]}:{x["ALT"]}', axis=1)
#    ldsc.loc[ldsc['SNP'].isin(snplist_chr.index), 'SNP'] = ldsc.loc[ldsc['SNP'].isin(snplist_chr.index), 'SNP'].map(snplist_chr['SNPID'])
#    ldsc.to_csv(f'new_{i}.l2.ldscore.gz', sep='\t', header=True, index=False, compression='gzip')

ldscs = pd.DataFrame()
for i in range(1,23):
    ldsc = pd.read_table(f'../eur_w_ld_chr_rsid/{i}.l2.ldscore.gz', sep='\t', header=0)
    ldscs = pd.concat([ldscs, ldsc], ignore_index=True)
ldscs.index = ldscs['SNP']
snplist = snplist.loc[snplist.index.isin(ldscs.index)]
snplist['CHR'] = snplist.index.map(ldscs['CHR'])
snplist['BP'] = snplist.index.map(ldscs['BP'])
snplist['REF'] = snplist.apply(lambda x: (fa[f'chr{x["CHR"]}'][x['BP']-1]).seq.upper(), axis=1)
snplist['ALT'] = np.where(snplist['A1']==snplist['REF'], snplist['A2'], snplist['A1'])
snplist['RSID'] = snplist.index
snplist.index = snplist['CHR'].astype(str) + ':' + snplist['BP'].astype(str) + ':' + snplist['REF'] + ':' + snplist['ALT']
snplist.to_csv('w_hm3.snplist', sep='\t', header=True, index=True)
