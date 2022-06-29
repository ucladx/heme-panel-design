## Design a hybridization capture-based cancer DNA sequencing panel for hematologic cancers

### Overview of required targets:
* Coding exons: Of all major transcripts of ~200 genes
* Non-coding: Cancer-associated variants from ClinVar
* Fusions: Common breakpoints for important gene fusions
* Germline SNPs: Commonly heterozygous ~1600 SNPs for CNVs

Compared to the [pan-cancer panel](https://github.com/ucladx/pancan-panel-design), this panel will have 10x deeper sequencing with a higher rate of duplicate reads per fragment, and UMI-based collapsing for increased variant calling precision at low allele fractions.

### Coding exons

Shortlisted 218 genes for a new panel in [`data/cancer_genes_review.txt`](data/cancer_genes_review.txt). Mutations per Mbp was calculated using TCGA+TARGET MuTect2 MAFs from NCI GDC, and gene sizes from Gencode v35. Mean gain/loss was calculated using TCGA+TARGET Gistic2 gene-level absolute CN from NCI GDC.

Extract the gene names and their Ensembl ENSG IDs:
```bash
mkdir -p data targets
cut -f1,2 data/cancer_genes_review.txt | grep -v Gene | sort > data/exon_targets_gene_list.txt
```

Create a BED file for these genes' coding regions with 2bp flanks, using their Gencode basic transcripts except level 3 (not verified nor curated):
```bash
curl -L ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/gencode.v40.basic.annotation.gff3.gz -o data/gencode.v40.basic.annotation.gff3.gz
gzip -dc data/gencode.v40.basic.annotation.gff3.gz | grep -w "$(cut -f2 data/exon_targets_gene_list.txt)" | perl -a -F'\t' -ne '%t=map{split("=")} split(";",$F[8]); if(($t{gene_type} eq "protein_coding" and $F[2] eq "CDS" and $t{level} ne "3" and $t{ID}!~m/PAR/) or ($t{gene_type}=~/lncRNA|miRNA|pseudogene/ and $F[2] eq "exon")){$F[3]-=3; $F[4]+=2; print join("\t",@F[0,3,4],$t{gene_name}.":".$F[2],@F[5,6])."\n"}' | sort -s -k1,1V -k2,2n -k3,3n | bedtools merge -i - -c 4 -o distinct > targets/exon_targets_hg38.bed
```

Try different combinations of transcripts/exons/annotations in steps above and measure total target space:
```bash
awk -F"\t" '{sum+=$3-$2} END {print sum}' targets/exon_targets_hg38.bed
```

647414 bps in total; 648272 if we include level 3 transcripts; 689066 if we used 8bp flanks

### Non-coding

For each gene in the panel, target the ends of 5'UTRs of MANE transcripts where mutations could cause loss of function:
```bash
gzip -dc data/gencode.v40.basic.annotation.gff3.gz | grep -w "$(cut -f2 data/exon_targets_gene_list.txt)" | perl -a -F'\t' -ne '%t=map{split("=")} split(";",$F[8]); if($t{gene_type} eq "protein_coding" and $F[2] eq "five_prime_UTR" and $t{tag}=~/MANE/ and $t{ID}!~m/PAR/){print join("\t",$F[0],$F[3]-2,($F[3]+118<$F[4]?$F[3]+118:$F[4]),$t{gene_name}.":5pUTR")."\n".join("\t",$F[0],($F[4]-118>$F[3]?$F[4]-118:$F[3]),$F[4]+2,$t{gene_name}.":5pUTR")."\n"}' | sort -s -k1,1V -k2,2n -k3,3n | bedtools merge -i - -c 4 -o distinct >> targets/non_coding_targets_hg38.bed
sort -su -k1,1V -k2,2n -k3,3n targets/non_coding_targets_hg38.bed -o targets/non_coding_targets_hg38.bed
```

Fetch loci of Pathogenic and Likely Pathogenic (P/LP) mutations with decent evidence from ClinVar related to cancer:
```bash
curl -sL 'https://api.genome.ucsc.edu/getData/track?genome=hg38;track=clinvarMain' | jq -r '.clinvarMain[] | [.chrom,.chromStart,.chromEnd,._variantId,._clinSignCode,.reviewStatus,.phenotypeList] | @tsv' | perl -a -F'\t' -ne 'print if($F[4]=~m/^P|LP$/ && $F[5]=~m/practice guideline|expert panel|multiple submitters/i) && $F[6]=~m/cancer|lynch|neoplas|tumor|adenoma|carcinoma|li-fraumeni|polyposis|hippel-lindau/i' | cut -f1-4 | sort -su -k1,1V -k2,2n -k3,3n > data/clinvar_plp_muts_hg38.bed
```

Target the subset of P/LP ClinVar mutations that do not overlap targeted exons:
```bash
bedtools subtract -a data/clinvar_plp_muts_hg38.bed -b targets/exon_targets_hg38.bed | sed 's/$/:ClinVar/' | sort -su -k1,1V -k2,2n -k3,3n >> targets/non_coding_targets_hg38.bed
sort -s -k1,1V -k2,2n -k3,3n targets/non_coding_targets_hg38.bed -o targets/non_coding_targets_hg38.bed
```

### Fusions

* Target the most common fusion breakpoints for `ABL1, ALK, BCR, CBFB, ETV6, KMT2A, MYC, MYH11, NTRK1, NTRK2, NTRK3, NUP98, PBX1, PML, RARA, RUNX1, TCF3`

Download hg19 breakpoints from FusionGDB2 for the genes we care about, and convert them to an hg38 BED file:
```bash
curl -L https://compbio.uth.edu/FusionGDB2/tables/feature_importance_scores_top1pct_regions_only_assembled.txt | tr '_' '\t' > data/fusion_gdb2_breakpoints_hg19.txt
(cut -f1,3,4 data/fusion_gdb2_breakpoints_hg19.txt; cut -f2,5,6 data/fusion_gdb2_breakpoints_hg19.txt) | sort -u | grep -wE "^(ABL1|ALK|BCR|CBFB|ETV6|KMT2A|MYC|MYH11|NTRK1|NTRK2|NTRK3|NUP98|PBX1|PML|RARA|RUNX1|TCF3)" | awk -v OFS='\t' '{print $2, $3-1, $3, $1":FusionSite"}' | sort -s -k1,1V -k2,2n -k3,3n > targets/fusion_site_targets_hg19.bed
liftOver targets/fusion_site_targets_hg19.bed ~/src/vcf2maf/data/hg19_to_hg38.chain targets/fusion_site_targets_hg38.bed unmapped.bed
rm -f targets/fusion_site_targets_hg19.bed unmapped.bed
```

### Germline SNPs

Fetch the list of SNPs shortlisted for the pan-cancer panel:

Target XX SNPs within exon targets or <240bp near them, which gives XX genes at least 1 SNP likely to help detect LOH (don't choose >1 SNP <200bp apart):

Add XX more SNPs (1600 total, sufficient for LOH) that are most distant from their nearest SNPs:

### Probe design

Combined all targets into a single merged BED file:
```bash
cat targets/*_targets_hg38.bed | sort -s -k1,1V -k2,2n -k3,3n > data/ucla_mdl_targets_hg38.bed
```

Estimated how many probes will be needed for 1x tiling (targets <=120bp get one 120bp probe each, others get total bps ÷ 120):
```bash
bedtools merge -i data/ucla_mdl_targets_hg38.bed -c 4 -o distinct | awk -F"\t" '{len=$3-$2; sum+=(len<120?120:len)} END {print sum/120}'
```
