## Design a capture-based DNA-seq panel for hematologic cancers

NOTE: *We are abandoning this panel design in favor of using probes in our freezer from [this group purchase](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10424130/table/tbl3).*

GOAL: Highly sensitive detection of somatic mutations in myeloid and lymphoid specimens, for screening and monitoring of patients with hematologic cancers. Compared to the [pan-cancer panel](https://github.com/ucladx/pancan-panel-design), this panel will have 10x deeper sequencing with a higher rate of duplicate reads per fragment, and UMI-based collapsing of duplicates into higher quality reads for increased specificity at low variant allele fractions (VAFs).

### Overview of required targets:
* Coding exons: Of all major transcripts of ~200 genes
* Non-coding: Cancer-associated variants from ClinVar
* Fusions: Common breakpoints for important gene fusions
* Germline SNPs: Commonly heterozygous ~1600 SNPs for CNVs

### Coding exons

Shortlisted 218 genes for a new panel in [`data/cancer_genes_review.txt`](data/cancer_genes_review.txt). Mutations per Mbp was calculated using TCGA+TARGET MuTect2 MAFs from NCI GDC, and gene sizes from Gencode v35. Mean gain/loss was calculated using TCGA+TARGET Gistic2 gene-level absolute CN from NCI GDC.

Extract the gene names and their Ensembl ENSG IDs:
```bash
mkdir -p data targets
cut -f1,2 data/cancer_genes_review.txt | grep -v Gene | sort > data/exon_targets_gene_list.txt
```

Manually added two more genes per request from colleagues: DKC1 and ERG

Create a BED file for these genes' coding regions with 2bp flanks, using their Gencode basic transcripts except level 3 (not verified nor curated):
```bash
curl -L ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.basic.annotation.gff3.gz -o data/gencode.v40.basic.annotation.gff3.gz
gzip -dc data/gencode.v40.basic.annotation.gff3.gz | grep -w "$(cut -f2 data/exon_targets_gene_list.txt)" | perl -a -F'\t' -ne '%t=map{split("=")} split(";",$F[8]); if(($t{gene_type} eq "protein_coding" and $F[2] eq "CDS" and $t{level} ne "3" and $t{ID}!~m/PAR/) or ($t{gene_type}=~/lncRNA|miRNA|pseudogene/ and $F[2] eq "exon")){$F[3]-=3; $F[4]+=2; print join("\t",@F[0,3,4],$t{gene_name}.":".$F[2],@F[5,6])."\n"}' | sort -s -k1,1V -k2,2n -k3,3n | bedtools merge -i - -c 4 -o distinct > targets/exon_targets_hg38.bed
```

Try different combinations of transcripts/exons/annotations in steps above and measure total target space:
```bash
awk -F"\t" '{sum+=$3-$2} END {print sum}' targets/exon_targets_hg38.bed
```

650581 bps in total; 651439 if we include level 3 transcripts; 693463 if we used 8bp flanks

### Non-coding

For each gene in the panel, target the ends of 5'UTRs of MANE transcripts where mutations could cause loss of function:
```bash
gzip -dc data/gencode.v40.basic.annotation.gff3.gz | grep -w "$(cut -f2 data/exon_targets_gene_list.txt)" | perl -a -F'\t' -ne '%t=map{split("=")} split(";",$F[8]); if($t{gene_type} eq "protein_coding" and $F[2] eq "five_prime_UTR" and $t{tag}=~/MANE/ and $t{ID}!~m/PAR/){print join("\t",$F[0],$F[3]-2,($F[3]+118<$F[4]?$F[3]+118:$F[4]),$t{gene_name}.":5pUTR")."\n".join("\t",$F[0],($F[4]-118>$F[3]?$F[4]-118:$F[3]),$F[4]+2,$t{gene_name}.":5pUTR")."\n"}' | sort -s -k1,1V -k2,2n -k3,3n | bedtools merge -i - -c 4 -o distinct > targets/non_coding_targets_hg38.bed
```

Fetch loci of Pathogenic and Likely Pathogenic (P/LP) mutations with decent evidence from ClinVar related to cancer:
```bash
curl -sL 'https://api.genome.ucsc.edu/getData/track?genome=hg38;track=clinvarMain' | jq -r '.clinvarMain[] | [.chrom,.chromStart,.chromEnd,._variantId,._clinSignCode,.reviewStatus,.phenotypeList] | @tsv' | perl -a -F'\t' -ne 'print if($F[4]=~m/^P|LP$/ && $F[5]=~m/practice guideline|expert panel|multiple submitters/i) && $F[6]=~m/cancer|lynch|neoplas|tumor|adenoma|carcinoma|li-fraumeni|polyposis|hippel-lindau/i' | cut -f1-4 | sort -su -k1,1V -k2,2n -k3,3n > data/clinvar_plp_muts_hg38.bed
```

Target the P/LP ClinVar variants that are <1200bp to the targeted exons, but not overlapping:
```bash
bedtools closest -d -a data/clinvar_plp_muts_hg38.bed -b targets/exon_targets_hg38.bed -g /srv/ref/hg38.fa.fai | awk -v OFS='\t' '{sub(/:CDS$/,":ClinVar",$8); if($9>0 && $9<1200) print $1, $2, $3, $8}' >> targets/non_coding_targets_hg38.bed
sort -su -k1,1V -k2,2n -k3,3n targets/non_coding_targets_hg38.bed -o targets/non_coding_targets_hg38.bed
```

### Fusions

Target reasonably small introns of `ABL1, KMT2A, MYC, NTRK1, RARA` that are likely to contain breakpoints for important gene fusions

Download hg19 splice junction loci from FusionGDB2 for the genes we care about, and convert them to an hg38 BED file:
```bash
curl -L https://compbio.uth.edu/FusionGDB2/tables/feature_importance_scores_top1pct_regions_only_assembled.txt | tr '_' '\t' > data/fusion_gdb2_splice_jns_hg19.txt
(cut -f1,3,4 data/fusion_gdb2_splice_jns_hg19.txt; cut -f2,5,6 data/fusion_gdb2_splice_jns_hg19.txt) | sort -u | grep -wE "^(ABL1|KMT2A|MYC|NTRK1|RARA)" | awk -v OFS='\t' '{print $2, $3-1, $3, $1}' | sort -s -k1,1V -k2,2n -k3,3n > data/fusion_splice_jns_hg19.bed
liftOver data/fusion_splice_jns_hg19.bed ~/src/vcf2maf/data/hg19_to_hg38.chain data/fusion_splice_jns_hg38.bed unmapped.bed
rm -f data/fusion_splice_jns_hg19.bed unmapped.bed
```

For each splice junction, target the adjacent intron unless it is 10Kbp or longer:
```bash
gzip -dc data/gencode.v40.basic.annotation.gff3.gz | perl -a -F'\t' -ne '%t=map{split("=")} split(";",$F[8]); if($t{tag}=~/MANE/ and $F[2] eq "CDS" and $t{ID}!~m/PAR/){print join("\t",@F[0,3,4])."\n"}' | sort -su -k1,1V -k2,2n -k3,3n | bedtools complement -L -i - -g /srv/ref/hg38.fa.fai | bedtools window -w 2 -a data/fusion_splice_jns_hg38.bed -b - | awk -v OFS='\t' '{if($7-$6<10000) print $5, $6, $7, $4":FusionSite"}' | sort -su -k1,1V -k2,2n -k3,3n | bedtools merge -i - -c 4 -o distinct > targets/fusion_targets_hg38.bed
```

### Germline SNPs

Based on testing done for the pan-cancer panel, fetch the list of well performing SNPs for CNV calling:
```bash
curl -L https://raw.githubusercontent.com/ucladx/pancan-panel-design/dbf4ef9/data/snp_candidates_good_grch38.bed -o data/snp_candidates_good_hg38.bed
```

Target 447 SNPs within exon targets or <240bp near them, which gives 135 genes at least 1 SNP likely to help detect LOH (don't choose >1 SNP <200bp apart):
```bash
bedtools window -w 240 -a targets/exon_targets_hg38.bed -b data/snp_candidates_good_hg38.bed | cut -f5-8 | sort -su -k1,1V -k2,2n -k3,3n | bedtools spacing -i - | awk -F'\t' '{if($5>=200) print}' | cut -f1-4 | sed 's/$/:SNP_LOH/' > targets/snp_targets_hg38.bed
```

Add 1153 more SNPs (1600 total, sufficient for large CNVs) that are most distant from their nearest SNPs:
```bash
grep SNP_LOH$ targets/snp_targets_hg38.bed | bedtools slop -b 200 -g /srv/ref/hg38.fa.fai -i | bedtools subtract -a data/snp_candidates_good_hg38.bed -b - | bedtools spacing -i - | sort -k7,7rn | head -n1153 | cut -f1-4 | sed 's/$/:SNP_CNV/' >> targets/snp_targets_hg38.bed
sort -s -k1,1V -k2,2n -k3,3n targets/snp_targets_hg38.bed -o targets/snp_targets_hg38.bed
```

### Probe design

Estimate how many probes will be needed for 1x tiling (targets <=120bp get one 120bp probe each, others get total bps ÷ 120):
```bash
cat targets/*_targets_hg38.bed | sort -s -k1,1V -k2,2n -k3,3n | bedtools merge -i - -c 4 -o distinct | awk -F"\t" '{len=$3-$2; sum+=(len<120?120:len)} END {print sum/120}'
```
