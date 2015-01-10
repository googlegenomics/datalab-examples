
# Genome-wide association study (GWAS)

This notebook demonstrates conducting a genome-wide association study using the
public 1000 Genomes dataset stored in BigQuery.

### In this notebook you will
* Explore the 1000 Genomes dataset available in BigQuery
* Use the `%%bq_sql` statement to write and execute SQL statements within the
notebook
* Extract data from BigQuery and create a local dataset that can be manipulated
in Python for visualization and further analysis

### Attribution
This notebook is based on the [Google
Genomics](https://cloud.google.com/genomics/) BigQuery examples here:
* [GWAS Chi-squared Test](https://github.com/googlegenomics/bigquery-
examples/blob/master/1000genomes/sql/gwas-pattern-chi-squared-test.sql)
* [GWAS z-test](https://github.com/googlegenomics/bigquery-
examples/blob/master/1000genomes/sql/gwas-pattern-two-proportion-z-test.sql)
----

## Experiment design and dataset construction

In this experiment, we'll be identifying variant positions within chromosome 12
that differ significantly between the case and control groups. The case group
for the purposes of this notebook will be individuals from the "EAS" (East
Asian) super population.

To work with genomics data stored within BigQuery, we'll import the gcp.bigquery
library, which provides a set of high-level wrappers for manipulating SQL within
the notebook, issuing queries to BigQuery and translating result sets into
Pandas dataframes.


    import gcp.bigquery as bq

### Table selection and summary

Variant data from the 1000 genomes dataset is publicly accessible within
BigQuery. We can view the schema for the variants table via the gcp.bigquery
library.


    variants_table = bq.table('genomics-public-data:1000_genomes.variants')
    variants_table.schema()




<table><tr><th>name</th><th>data_type</th><th>mode</th><th>description</th></tr><tr><td>reference_name</td><td>STRING</td><td>NULLABLE</td><td>An identifier from the reference genome or an angle-bracketed ID String pointing to a contig in the assembly file.</td></tr><tr><td>start</td><td>INTEGER</td><td>NULLABLE</td><td>The reference position, with the first base having position 0.</td></tr><tr><td>end</td><td>INTEGER</td><td>NULLABLE</td><td>INFO=&lt;ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record"&gt;</td></tr><tr><td>reference_bases</td><td>STRING</td><td>NULLABLE</td><td>Each base must be one of A,C,G,T,N (case insensitive). Multiple bases are permitted. The value in the POS field refers to the position of the first base in the String.</td></tr><tr><td>alternate_bases</td><td>STRING</td><td>REPEATED</td><td>List of alternate non-reference alleles called on at least one of the samples. ("at least one" not true for this dataset)</td></tr><tr><td>quality</td><td>FLOAT</td><td>NULLABLE</td><td>phred-scaled quality score for the assertion made in ALT.</td></tr><tr><td>filter</td><td>STRING</td><td>REPEATED</td><td>PASS if this position has passed all filters, i.e. a call is made at this position. Otherwise, if the site has not passed all filters, a list of codes for filters that fail.</td></tr><tr><td>names</td><td>STRING</td><td>REPEATED</td><td>List of unique identifiers for the variant where available.</td></tr><tr><td>call</td><td>RECORD</td><td>REPEATED</td><td>Per-sample measurements.</td></tr><tr><td>call.call_set_id</td><td>STRING</td><td>NULLABLE</td><td>The id of the callset from which this data was exported from the Google Genomics Variants API.</td></tr><tr><td>call.call_set_name</td><td>STRING</td><td>NULLABLE</td><td>Sample identifier.</td></tr><tr><td>call.genotype</td><td>INTEGER</td><td>REPEATED</td><td>List of genotypes.</td></tr><tr><td>call.phaseset</td><td>STRING</td><td>NULLABLE</td><td>If this value is null, the data is unphased.  Otherwise it is phased.</td></tr><tr><td>call.genotype_likelihood</td><td>FLOAT</td><td>REPEATED</td><td>List of genotype likelihoods.</td></tr><tr><td>call.DP</td><td>INTEGER</td><td>NULLABLE</td><td>FORMAT=&lt;ID=DP,Number=1,Type=Integer,Description="# high-quality bases"&gt;</td></tr><tr><td>call.DS</td><td>FLOAT</td><td>NULLABLE</td><td>FORMAT=&lt;ID=DS,Number=1,Type=Float,Description="Genotype dosage from MaCH/Thunder"&gt;</td></tr><tr><td>call.FT</td><td>STRING</td><td>NULLABLE</td><td></td></tr><tr><td>call.GQ</td><td>STRING</td><td>NULLABLE</td><td>FORMAT=&lt;ID=GQ,Number=1,Type=Float,Description="Genotype quality"&gt;</td></tr><tr><td>call.PL</td><td>INTEGER</td><td>REPEATED</td><td>FORMAT=&lt;ID=PL,Number=.,Type=Integer,Description="List of Phred-scaled genotype likelihoods, number of values is (#ALT+1)*(#ALT+2)/2"&gt;</td></tr><tr><td>call.SP</td><td>INTEGER</td><td>NULLABLE</td><td>FORMAT=&lt;ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value"&gt;</td></tr><tr><td>AA</td><td>STRING</td><td>NULLABLE</td><td>INFO=&lt;ID=AA,Number=1,Type=String,Description="Ancestral Allele, ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/technical/reference/ancestral_alignments/README"&gt;</td></tr><tr><td>AC</td><td>INTEGER</td><td>REPEATED</td><td>INFO=&lt;ID=AC,Number=.,Type=Integer,Description="Alternate Allele Count"&gt;</td></tr><tr><td>AC1</td><td>INTEGER</td><td>NULLABLE</td><td>INFO=&lt;ID=AC1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele count (no HWE assumption)"&gt;</td></tr><tr><td>AF</td><td>FLOAT</td><td>NULLABLE</td><td>INFO=&lt;ID=AF,Number=1,Type=Float,Description="Global Allele Frequency based on AC/AN"&gt;</td></tr><tr><td>AF1</td><td>FLOAT</td><td>NULLABLE</td><td>INFO=&lt;ID=AF1,Number=1,Type=Float,Description="Max-likelihood estimate of the first ALT allele frequency (assuming HWE)"&gt;</td></tr><tr><td>AFR_AF</td><td>FLOAT</td><td>NULLABLE</td><td>INFO=&lt;ID=AFR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from AFR based on AC/AN"&gt;</td></tr><tr><td>AMR_AF</td><td>FLOAT</td><td>NULLABLE</td><td>INFO=&lt;ID=AMR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from AMR based on AC/AN"&gt;</td></tr><tr><td>AN</td><td>INTEGER</td><td>NULLABLE</td><td>INFO=&lt;ID=AN,Number=1,Type=Integer,Description="Total Allele Count"&gt;</td></tr><tr><td>ASN_AF</td><td>FLOAT</td><td>NULLABLE</td><td>INFO=&lt;ID=ASN_AF,Number=1,Type=Float,Description="Allele Frequency for samples from ASN based on AC/AN"&gt;</td></tr><tr><td>AVGPOST</td><td>FLOAT</td><td>NULLABLE</td><td>INFO=&lt;ID=AVGPOST,Number=1,Type=Float,Description="Average posterior probability from MaCH/Thunder"&gt;</td></tr><tr><td>CIEND</td><td>INTEGER</td><td>REPEATED</td><td>INFO=&lt;ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants"&gt;</td></tr><tr><td>CIPOS</td><td>INTEGER</td><td>REPEATED</td><td>INFO=&lt;ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants"&gt;</td></tr><tr><td>DP</td><td>INTEGER</td><td>NULLABLE</td><td>INFO=&lt;ID=DP,Number=1,Type=Integer,Description="Raw read depth"&gt;</td></tr><tr><td>DP4</td><td>INTEGER</td><td>REPEATED</td><td>INFO=&lt;ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases"&gt;</td></tr><tr><td>ERATE</td><td>FLOAT</td><td>NULLABLE</td><td>INFO=&lt;ID=ERATE,Number=1,Type=Float,Description="Per-marker Mutation rate from MaCH/Thunder"&gt;</td></tr><tr><td>EUR_AF</td><td>FLOAT</td><td>NULLABLE</td><td>INFO=&lt;ID=EUR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from EUR based on AC/AN"&gt;</td></tr><tr><td>FQ</td><td>FLOAT</td><td>NULLABLE</td><td>INFO=&lt;ID=FQ,Number=1,Type=Float,Description="Phred probability of all samples being the same"&gt;</td></tr><tr><td>G3</td><td>FLOAT</td><td>REPEATED</td><td>INFO=&lt;ID=G3,Number=3,Type=Float,Description="ML estimate of genotype frequencies"&gt;</td></tr><tr><td>HOMLEN</td><td>INTEGER</td><td>NULLABLE</td><td>INFO=&lt;ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints"&gt;</td></tr><tr><td>HOMSEQ</td><td>STRING</td><td>NULLABLE</td><td>INFO=&lt;ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints"&gt;</td></tr><tr><td>HWE</td><td>FLOAT</td><td>NULLABLE</td><td>INFO=&lt;ID=HWE,Number=1,Type=Float,Description="Chi^2 based HWE test P-value based on G3"&gt;</td></tr><tr><td>LDAF</td><td>FLOAT</td><td>NULLABLE</td><td>INFO=&lt;ID=LDAF,Number=1,Type=Float,Description="MLE Allele Frequency Accounting for LD"&gt;</td></tr><tr><td>MQ</td><td>INTEGER</td><td>NULLABLE</td><td>INFO=&lt;ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads"&gt;</td></tr><tr><td>PV4</td><td>FLOAT</td><td>REPEATED</td><td>INFO=&lt;ID=PV4,Number=4,Type=Float,Description="P-values for strand bias, baseQ bias, mapQ bias and tail distance bias"&gt;</td></tr><tr><td>RSQ</td><td>FLOAT</td><td>NULLABLE</td><td>INFO=&lt;ID=RSQ,Number=1,Type=Float,Description="Genotype imputation quality from MaCH/Thunder"&gt;</td></tr><tr><td>SNPSOURCE</td><td>STRING</td><td>REPEATED</td><td>INFO=&lt;ID=SNPSOURCE,Number=.,Type=String,Description="indicates if a snp was called when analysing the low coverage or exome alignment data"&gt;</td></tr><tr><td>SOURCE</td><td>STRING</td><td>REPEATED</td><td></td></tr><tr><td>SVLEN</td><td>INTEGER</td><td>NULLABLE</td><td>INFO=&lt;ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles"&gt;</td></tr><tr><td>SVTYPE</td><td>STRING</td><td>NULLABLE</td><td>INFO=&lt;ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant"&gt;</td></tr><tr><td>THETA</td><td>FLOAT</td><td>NULLABLE</td><td>INFO=&lt;ID=THETA,Number=1,Type=Float,Description="Per-marker Transition rate from MaCH/Thunder"&gt;</td></tr><tr><td>VT</td><td>STRING</td><td>NULLABLE</td><td>INFO=&lt;ID=VT,Number=1,Type=String,Description="indicates what type of variant the line represents"&gt;</td></tr></table>



### Classifying per-call variant positions into variant/non-variant groups

We can tally the reference/alternate allele accounts for *individual* variant
positions within chromosome 12. The field `call.genotype` is an integer ranging
from `[0, num_alternate_bases]`. A value of zero indicates that the genotype for
the call is the same as the reference (i.e., non-variant). A value of 1 would
indicate that the genotype for the call is the 1st value in the list of
alternate bases (likewise for values >1).


    %%bq_sql allele_counts
    
    SELECT 
        reference_name,
        start,
        reference_bases,
        alternate_bases,
        end,
        vt,    
        call.call_set_name AS call_set_name,
        (0 = first_allele) + (0 = second_allele) AS ref_count,
        (1 = first_allele) + (1 = second_allele) AS alt_count,    
    FROM
      FLATTEN((
        SELECT
          reference_name,
          start,
          reference_bases,
          alternate_bases,
          end,
          vt,
          call.call_set_name,
          NTH(1, call.genotype) WITHIN call AS first_allele,
          NTH(2, call.genotype) WITHIN call AS second_allele,
        FROM
          $variants_table
        WHERE
          reference_name = '12' -- i.e., chromosome 12
        HAVING
          -- Exclude calls _where one _or both alleles were NOT called (-1)
          0 <= first_allele
          AND 0 <= second_allele
        ),
        call)

Let's verify that our allele counts match our expectations before moving on. For
any given row, the alternate + reference counts should sum to 2 for this
experiment.


    allele_counts.sample()




<div>Number of rows: 5</div><div>Query job ID  : job_Vzesy7ezHD4LagW372EEklNp4Pc</div>
<table><tr><th>call_set_name</th><th>end</th><th>reference_bases</th><th>alt_count</th><th>reference_name</th><th>alternate_bases</th><th>ref_count</th><th>start</th><th>vt</th></tr><tr><td>HG00261</td><td>26584346</td><td>G</td><td>2</td><td>12</td><td>A</td><td>0</td><td>26584345</td><td>SNP</td></tr><tr><td>HG00593</td><td>26584346</td><td>G</td><td>2</td><td>12</td><td>A</td><td>0</td><td>26584345</td><td>SNP</td></tr><tr><td>NA12749</td><td>26584346</td><td>G</td><td>2</td><td>12</td><td>A</td><td>0</td><td>26584345</td><td>SNP</td></tr><tr><td>HG00150</td><td>26584346</td><td>G</td><td>2</td><td>12</td><td>A</td><td>0</td><td>26584345</td><td>SNP</td></tr><tr><td>NA19675</td><td>26584346</td><td>G</td><td>2</td><td>12</td><td>A</td><td>0</td><td>26584345</td><td>SNP</td></tr></table>



### Assigning case and control groups

Now we can join our allele counts with metadata available in the sample info
table. We'll use this sample metadata to split the set of genomes into case and
control groups based upon the super population group.


    sample_info_table = bq.table('genomics-public-data:1000_genomes.sample_info')


    %%bq_sql exp_groups
    
    SELECT
      super_population
      , IF('EAS' = super_population, TRUE, FALSE) AS is_case
      , call_set_name
      , reference_name
      , start
      , reference_bases
      , alternate_bases
      , end
      , vt
      , ref_count
      , alt_count
    FROM $allele_counts AS allele_counts
    JOIN $sample_info_table AS samples
      ON allele_counts.call_set_name = samples.sample


    exp_groups.sample()




<div>Number of rows: 5</div><div>Query job ID  : job_6A1Gh67SXh03LEqLYNfE1K4tQnc</div>
<table><tr><th>call_set_name</th><th>end</th><th>reference_bases</th><th>alt_count</th><th>super_population</th><th>alternate_bases</th><th>is_case</th><th>ref_count</th><th>start</th><th>vt</th><th>reference_name</th></tr><tr><td>HG00261</td><td>113541904</td><td>C</td><td>1</td><td>EUR</td><td>A</td><td>False</td><td>1</td><td>113541903</td><td>SNP</td><td>12</td></tr><tr><td>HG00593</td><td>113541904</td><td>C</td><td>0</td><td>EAS</td><td>A</td><td>True</td><td>2</td><td>113541903</td><td>SNP</td><td>12</td></tr><tr><td>NA12749</td><td>113541904</td><td>C</td><td>1</td><td>EUR</td><td>A</td><td>False</td><td>1</td><td>113541903</td><td>SNP</td><td>12</td></tr><tr><td>HG00150</td><td>113541904</td><td>C</td><td>0</td><td>EUR</td><td>A</td><td>False</td><td>2</td><td>113541903</td><td>SNP</td><td>12</td></tr><tr><td>NA19675</td><td>113541904</td><td>C</td><td>0</td><td>AMR</td><td>A</td><td>False</td><td>2</td><td>113541903</td><td>SNP</td><td>12</td></tr></table>



The variants table contains a few different types of variant: structural
variants ("SV"), indels ("INDEL") and SNPs ("SNP").


    %%bq_sql
    SELECT 
      vt,
      COUNT(*)
    FROM $exp_groups
    GROUP BY vt




<div>Number of rows: 3</div><div>Query job ID  : job_EfYz-DpQC6BUSQCZ5naWC7BjZfw</div>
<table><tr><th>f0_</th><th>vt</th></tr><tr><td>1921583664</td><td>SNP</td></tr><tr><td>73853052</td><td>INDEL</td></tr><tr><td>745836</td><td>SV</td></tr></table>



For the purposes of this experiment, let's limit the variants to only SNPs.


    %%bq_sql snps
    SELECT * 
    FROM $exp_groups
    WHERE vt='SNP'


    snps.sample()




<div>Number of rows: 5</div><div>Query job ID  : job_xskHq0jxzjZ_qcG5eePx-smYODs</div>
<table><tr><th>call_set_name</th><th>end</th><th>reference_bases</th><th>alt_count</th><th>super_population</th><th>alternate_bases</th><th>is_case</th><th>ref_count</th><th>start</th><th>vt</th><th>reference_name</th></tr><tr><td>HG00261</td><td>86787245</td><td>A</td><td>0</td><td>EUR</td><td>T</td><td>False</td><td>2</td><td>86787244</td><td>SNP</td><td>12</td></tr><tr><td>HG00593</td><td>86787245</td><td>A</td><td>0</td><td>EAS</td><td>T</td><td>True</td><td>2</td><td>86787244</td><td>SNP</td><td>12</td></tr><tr><td>NA12749</td><td>86787245</td><td>A</td><td>0</td><td>EUR</td><td>T</td><td>False</td><td>2</td><td>86787244</td><td>SNP</td><td>12</td></tr><tr><td>HG00150</td><td>86787245</td><td>A</td><td>0</td><td>EUR</td><td>T</td><td>False</td><td>2</td><td>86787244</td><td>SNP</td><td>12</td></tr><tr><td>NA19675</td><td>86787245</td><td>A</td><td>0</td><td>AMR</td><td>T</td><td>False</td><td>2</td><td>86787244</td><td>SNP</td><td>12</td></tr></table>



### Tallying reference/alternate allele counts for case/control groups

Now that we've assigned each call set to either the case or the control group,
we can tally up the counts of reference and alternate alleles within each of our
assigned case/control groups, for each variant position, like so:


    %%bq_sql grouped_counts
    
    SELECT
        reference_name,
        start,
        end,
        reference_bases,
        alternate_bases,
        vt,
    
        SUM(ref_count + alt_count) AS allele_count,
        SUM(ref_count) AS ref_count,
        SUM(alt_count) AS alt_count,
        SUM(IF(TRUE = is_case,
            INTEGER(ref_count + alt_count),
            0)) AS case_count,
        SUM(IF(FALSE = is_case,
            INTEGER(ref_count + alt_count),
            0)) AS control_count,
        SUM(IF(TRUE = is_case,
            ref_count,
            0)) AS case_ref_count,
        SUM(IF(TRUE = is_case,
            alt_count,
            0)) AS case_alt_count,
        SUM(IF(FALSE = is_case,
            ref_count,
            0)) AS control_ref_count,
        SUM(IF(FALSE = is_case,
            alt_count,
            0)) AS control_alt_count
    FROM $snps
    GROUP BY
        reference_name,
        start,
        end,
        reference_bases,
        alternate_bases,
        vt


    print grouped_counts.sql

    
    SELECT
        reference_name,
        start,
        end,
        reference_bases,
        alternate_bases,
        vt,
    
        SUM(ref_count + alt_count) AS allele_count,
        SUM(ref_count) AS ref_count,
        SUM(alt_count) AS alt_count,
        SUM(IF(TRUE = is_case,
            INTEGER(ref_count + alt_count),
            0)) AS case_count,
        SUM(IF(FALSE = is_case,
            INTEGER(ref_count + alt_count),
            0)) AS control_count,
        SUM(IF(TRUE = is_case,
            ref_count,
            0)) AS case_ref_count,
        SUM(IF(TRUE = is_case,
            alt_count,
            0)) AS case_alt_count,
        SUM(IF(FALSE = is_case,
            ref_count,
            0)) AS control_ref_count,
        SUM(IF(FALSE = is_case,
            alt_count,
            0)) AS control_alt_count
    FROM (SELECT * 
    FROM (
    SELECT
      super_population
      , IF('EAS' = super_population, TRUE, FALSE) AS is_case
      , call_set_name
      , reference_name
      , start
      , reference_bases
      , alternate_bases
      , end
      , vt
      , ref_count
      , alt_count
    FROM (
    SELECT 
        reference_name,
        start,
        reference_bases,
        alternate_bases,
        end,
        vt,    
        call.call_set_name AS call_set_name,
        (0 = first_allele) + (0 = second_allele) AS ref_count,
        (1 = first_allele) + (1 = second_allele) AS alt_count,    
    FROM
      FLATTEN((
        SELECT
          reference_name,
          start,
          reference_bases,
          alternate_bases,
          end,
          vt,
          call.call_set_name,
          NTH(1, call.genotype) WITHIN call AS first_allele,
          NTH(2, call.genotype) WITHIN call AS second_allele,
        FROM
          [genomics-public-data:1000_genomes.variants]
        WHERE
          reference_name = '12' -- i.e., chromosome 12
        HAVING
          -- Exclude calls _where one _or both alleles were NOT called (-1)
          0 <= first_allele
          AND 0 <= second_allele
        ),
        call)) AS allele_counts
    JOIN [genomics-public-data:1000_genomes.sample_info] AS samples
      ON allele_counts.call_set_name = samples.sample)
    WHERE vt='SNP')
    GROUP BY
        reference_name,
        start,
        end,
        reference_bases,
        alternate_bases,
        vt


Again, validate that the results are sensical for the group level counts (still
per variant position).


    grouped_counts.sample()




<div>Number of rows: 5</div><div>Query job ID  : job_QL_sawtyV0C23udPCRQbZ7Rbc5Y</div>
<table><tr><th>control_alt_count</th><th>case_alt_count</th><th>end</th><th>reference_bases</th><th>alt_count</th><th>reference_name</th><th>alternate_bases</th><th>allele_count</th><th>ref_count</th><th>start</th><th>control_count</th><th>vt</th><th>case_ref_count</th><th>control_ref_count</th><th>case_count</th></tr><tr><td>252</td><td>45</td><td>99812732</td><td>T</td><td>297</td><td>12</td><td>C</td><td>2184</td><td>1887</td><td>99812731</td><td>1612</td><td>SNP</td><td>527</td><td>1360</td><td>572</td></tr><tr><td>5</td><td>0</td><td>99603035</td><td>T</td><td>5</td><td>12</td><td>C</td><td>2184</td><td>2179</td><td>99603034</td><td>1612</td><td>SNP</td><td>572</td><td>1607</td><td>572</td></tr><tr><td>0</td><td>2</td><td>99119418</td><td>C</td><td>2</td><td>12</td><td>G</td><td>2184</td><td>2182</td><td>99119417</td><td>1612</td><td>SNP</td><td>570</td><td>1612</td><td>572</td></tr><tr><td>2</td><td>1</td><td>99429980</td><td>G</td><td>3</td><td>12</td><td>T</td><td>2184</td><td>2181</td><td>99429979</td><td>1612</td><td>SNP</td><td>571</td><td>1610</td><td>572</td></tr><tr><td>2</td><td>0</td><td>99875616</td><td>T</td><td>2</td><td>12</td><td>C</td><td>2184</td><td>2182</td><td>99875615</td><td>1612</td><td>SNP</td><td>572</td><td>1610</td><td>572</td></tr></table>



## Quantify the statistical significance at each variant positions

We can quantify the statistical significance of each variant position using the
Chi-squared test. Furthermore, we can restrict our result set to *only*
statistically significant variant positions for this experiment by ranking each
position by its statistical signficance (decreasing) and thresholding the
results for significance at `p <= 5e-8` (chi-squared score >= 29.7).


    %%bq_sql results
    
    SELECT
      reference_name,
      start,
      end,
      reference_bases,
      alternate_bases,
      vt,
      case_count,
      control_count,
      allele_count,
      ref_count,
      alt_count,
      case_ref_count,
      case_alt_count,
      control_ref_count,
      control_alt_count,
      ROUND(
        POW(case_ref_count - (ref_count/allele_count)*case_count,
          2)/((ref_count/allele_count)*case_count) +
        POW(control_ref_count - (ref_count/allele_count)*control_count,
          2)/((ref_count/allele_count)*control_count) +
        POW(case_alt_count - (alt_count/allele_count)*case_count,
          2)/((alt_count/allele_count)*case_count) +
        POW(control_alt_count - (alt_count/allele_count)*control_count,
          2)/((alt_count/allele_count)*control_count),
        3)
      AS chi_squared_score
    FROM $grouped_counts
    WHERE
      # For chi-squared, expected counts must be at least 5 for each group
      (ref_count/allele_count)*case_count >= 5.0
      AND (ref_count/allele_count)*control_count >= 5.0
      AND (alt_count/allele_count)*case_count >= 5.0
      AND (alt_count/allele_count)*control_count >= 5.0
    HAVING
      # Chi-squared critical value for df=1, p-value=5*10^-8 is 29.71679
      chi_squared_score >= 29.71679
    ORDER BY
      chi_squared_score DESC,
      allele_count DESC

A quick sanity check of the results reveals that the positions deemed
significant do in fact have significantly different case/control counts for the
alternate/reference bases.


    results.sample()




<div>Number of rows: 5</div><div>Query job ID  : job_SmjoYfx0dFM_NYxPo2w0uoWCdrc</div>
<table><tr><th>control_alt_count</th><th>case_alt_count</th><th>end</th><th>reference_bases</th><th>alt_count</th><th>reference_name</th><th>alternate_bases</th><th>allele_count</th><th>ref_count</th><th>start</th><th>control_count</th><th>vt</th><th>case_ref_count</th><th>chi_squared_score</th><th>control_ref_count</th><th>case_count</th></tr><tr><td>19</td><td>352</td><td>110571374</td><td>A</td><td>371</td><td>12</td><td>G</td><td>2184</td><td>1813</td><td>110571373</td><td>1612</td><td>SNP</td><td>220</td><td>1090.781</td><td>1593</td><td>572</td></tr><tr><td>139</td><td>456</td><td>22509095</td><td>T</td><td>595</td><td>12</td><td>C</td><td>2184</td><td>1589</td><td>22509094</td><td>1612</td><td>SNP</td><td>116</td><td>1076.666</td><td>1473</td><td>572</td></tr><tr><td>142</td><td>457</td><td>22525353</td><td>A</td><td>599</td><td>12</td><td>C</td><td>2184</td><td>1585</td><td>22525352</td><td>1612</td><td>SNP</td><td>115</td><td>1071.835</td><td>1470</td><td>572</td></tr><tr><td>143</td><td>457</td><td>22524241</td><td>T</td><td>600</td><td>12</td><td>C</td><td>2184</td><td>1584</td><td>22524240</td><td>1612</td><td>SNP</td><td>115</td><td>1068.856</td><td>1469</td><td>572</td></tr><tr><td>93</td><td>413</td><td>33240568</td><td>C</td><td>506</td><td>12</td><td>G</td><td>2184</td><td>1678</td><td>33240567</td><td>1612</td><td>SNP</td><td>159</td><td>1046.758</td><td>1519</td><td>572</td></tr></table>



### Computing Chi-squared statistics in BigQuery vs Python vs R

Let's compare these BigQuery-computed Chi-squared scores to ones calculated via
Python's statistical packages


    import numpy as np
    from scipy.stats import chi2_contingency
    
    chi2, p, dof, expected = chi2_contingency(np.array([ 
        [220, 352], # case 
        [1593, 19]  # control
    ]))
    
    print 'Python Chi-sq score = %.3f' % chi2
    print 'BigQuery was 1090.781'

    Python Chi-sq score = 1086.505
    BigQuery was 1090.781


Comparing Python statistics to those computed in R
```
> chisq.test(rbind(c(220, 352), c(1593, 19)))

        Pearson's Chi-squared test with Yates' continuity correction

data:  rbind(c(220, 352), c(1593, 19))
X-squared = 1086.505, df = 1, p-value < 2.2e-16
```


    %%bq_sql top1k_snps
    SELECT
      chi_squared_score
      , control_ref_count 
      , control_alt_count 
      , case_ref_count
      , case_alt_count
    FROM $results
    LIMIT 1000


    chisq_df = top1k_snps.results().to_dataframe()
    chisq_df.shape
    chisq_df.ix[1]




    case_alt_count        456.000
    case_ref_count        116.000
    chi_squared_score    1076.666
    control_alt_count     139.000
    control_ref_count    1473.000
    Name: 1, dtype: float64




    for i in xrange(len(chisq_df)):
        # construct a contingency table for the reference/alternate allele counts
        ctable = np.array([
            [chisq_df.ix[i]['case_ref_count'],    chisq_df.ix[i]['case_alt_count']],
            [chisq_df.ix[i]['control_ref_count'], chisq_df.ix[i]['control_alt_count']]])
        
        py_chi2, pvalue, dof, expected = chi2_contingency(ctable)
        bq_chi2 = chisq_df.ix[i]['chi_squared_score']
        print 'BQ: %.3f, SciPy: %.3f: delta=%.3f' % (bq_chi2, py_chi2, abs(bq_chi2 - py_chi2))


    BQ: 1090.781, SciPy: 1086.505: delta=4.276
    BQ: 1076.666, SciPy: 1073.082: delta=3.584
    BQ: 1071.835, SciPy: 1068.266: delta=3.569
    BQ: 1068.856, SciPy: 1065.294: delta=3.562
    BQ: 1046.758, SciPy: 1043.029: delta=3.729
    BQ: 1041.396, SciPy: 1037.664: delta=3.732
    BQ: 1035.231, SciPy: 1031.721: delta=3.510
    BQ: 1033.299, SciPy: 1029.576: delta=3.723
    BQ: 1025.179, SciPy: 1021.494: delta=3.685
    BQ: 1022.474, SciPy: 1019.082: delta=3.392
    BQ: 1017.172, SciPy: 1013.468: delta=3.704
    BQ: 1011.200, SciPy: 1007.545: delta=3.655
    BQ: 1011.200, SciPy: 1007.545: delta=3.655
    BQ: 1011.200, SciPy: 1007.545: delta=3.655
    BQ: 1009.720, SciPy: 1006.322: delta=3.398
    BQ: 1007.764, SciPy: 1004.401: delta=3.363
    BQ: 1005.137, SciPy: 1001.780: delta=3.357
    BQ: 1003.169, SciPy: 999.524: delta=3.645
    BQ: 1000.886, SciPy: 997.425: delta=3.461
    BQ: 999.162, SciPy: 995.522: delta=3.640
    BQ: 999.162, SciPy: 995.522: delta=3.640
    BQ: 999.162, SciPy: 995.522: delta=3.640
    BQ: 992.621, SciPy: 988.998: delta=3.623
    BQ: 989.933, SciPy: 986.489: delta=3.444
    BQ: 987.079, SciPy: 983.642: delta=3.437
    BQ: 983.945, SciPy: 980.579: delta=3.366
    BQ: 975.998, SciPy: 972.658: delta=3.340
    BQ: 973.262, SciPy: 969.827: delta=3.435
    BQ: 971.840, SciPy: 968.489: delta=3.351
    BQ: 971.335, SciPy: 967.627: delta=3.708
    BQ: 970.550, SciPy: 967.144: delta=3.406
    BQ: 969.843, SciPy: 966.426: delta=3.417
    BQ: 967.848, SciPy: 964.511: delta=3.337
    BQ: 967.399, SciPy: 964.182: delta=3.217
    BQ: 967.399, SciPy: 964.182: delta=3.217
    BQ: 966.391, SciPy: 963.170: delta=3.221
    BQ: 965.033, SciPy: 961.821: delta=3.212
    BQ: 963.831, SciPy: 960.499: delta=3.332
    BQ: 961.104, SciPy: 957.838: delta=3.266
    BQ: 955.470, SciPy: 951.882: delta=3.588
    BQ: 955.470, SciPy: 951.882: delta=3.588
    BQ: 944.493, SciPy: 941.217: delta=3.276
    BQ: 944.493, SciPy: 941.217: delta=3.276
    BQ: 940.991, SciPy: 937.402: delta=3.589
    BQ: 935.395, SciPy: 932.137: delta=3.258
    BQ: 934.483, SciPy: 930.912: delta=3.571
    BQ: 932.851, SciPy: 929.598: delta=3.253
    BQ: 932.626, SciPy: 929.072: delta=3.554
    BQ: 932.626, SciPy: 929.072: delta=3.554
    BQ: 931.249, SciPy: 927.687: delta=3.562
    BQ: 930.112, SciPy: 926.941: delta=3.171
    BQ: 928.029, SciPy: 924.476: delta=3.553
    BQ: 927.348, SciPy: 923.791: delta=3.557
    BQ: 924.825, SciPy: 921.412: delta=3.413
    BQ: 918.873, SciPy: 915.714: delta=3.159
    BQ: 918.363, SciPy: 914.600: delta=3.763
    BQ: 916.726, SciPy: 913.515: delta=3.211
    BQ: 911.137, SciPy: 907.596: delta=3.541
    BQ: 907.489, SciPy: 904.262: delta=3.227
    BQ: 907.319, SciPy: 903.936: delta=3.383
    BQ: 905.927, SciPy: 902.091: delta=3.836
    BQ: 905.755, SciPy: 902.488: delta=3.267
    BQ: 904.733, SciPy: 901.180: delta=3.553
    BQ: 903.432, SciPy: 900.055: delta=3.377
    BQ: 902.672, SciPy: 898.913: delta=3.759
    BQ: 902.672, SciPy: 898.913: delta=3.759
    BQ: 902.672, SciPy: 898.913: delta=3.759
    BQ: 902.672, SciPy: 898.913: delta=3.759
    BQ: 902.672, SciPy: 898.913: delta=3.759
    BQ: 899.937, SciPy: 896.646: delta=3.291
    BQ: 899.737, SciPy: 896.320: delta=3.417
    BQ: 899.737, SciPy: 896.320: delta=3.417
    BQ: 899.737, SciPy: 896.320: delta=3.417
    BQ: 899.737, SciPy: 896.320: delta=3.417
    BQ: 898.285, SciPy: 894.457: delta=3.828
    BQ: 898.285, SciPy: 894.457: delta=3.828
    BQ: 897.239, SciPy: 893.924: delta=3.315
    BQ: 895.325, SciPy: 891.589: delta=3.736
    BQ: 895.325, SciPy: 891.589: delta=3.736
    BQ: 895.017, SciPy: 891.267: delta=3.750
    BQ: 892.007, SciPy: 888.600: delta=3.407
    BQ: 892.007, SciPy: 888.600: delta=3.407
    BQ: 891.865, SciPy: 888.836: delta=3.029
    BQ: 890.401, SciPy: 887.278: delta=3.123
    BQ: 889.210, SciPy: 884.950: delta=4.260
    BQ: 887.948, SciPy: 884.561: delta=3.387
    BQ: 887.881, SciPy: 884.510: delta=3.371
    BQ: 886.601, SciPy: 882.330: delta=4.271
    BQ: 883.554, SciPy: 880.228: delta=3.326
    BQ: 881.242, SciPy: 877.898: delta=3.344
    BQ: 880.329, SciPy: 876.827: delta=3.502
    BQ: 880.160, SciPy: 876.875: delta=3.285
    BQ: 877.334, SciPy: 873.892: delta=3.442
    BQ: 876.436, SciPy: 872.729: delta=3.707
    BQ: 876.436, SciPy: 872.729: delta=3.707
    BQ: 876.436, SciPy: 872.729: delta=3.707
    BQ: 875.283, SciPy: 871.874: delta=3.409
    BQ: 874.403, SciPy: 870.813: delta=3.590
    BQ: 873.597, SciPy: 870.324: delta=3.273
    BQ: 872.311, SciPy: 868.595: delta=3.716
    BQ: 872.048, SciPy: 868.318: delta=3.730
    BQ: 871.580, SciPy: 868.242: delta=3.338
    BQ: 869.004, SciPy: 865.590: delta=3.414
    BQ: 868.927, SciPy: 865.491: delta=3.436
    BQ: 868.677, SciPy: 865.332: delta=3.345
    BQ: 868.471, SciPy: 865.075: delta=3.396
    BQ: 867.728, SciPy: 864.380: delta=3.348
    BQ: 867.058, SciPy: 863.798: delta=3.260
    BQ: 865.922, SciPy: 862.661: delta=3.261
    BQ: 865.804, SciPy: 862.295: delta=3.509
    BQ: 865.706, SciPy: 862.228: delta=3.478
    BQ: 865.388, SciPy: 862.410: delta=2.978
    BQ: 862.038, SciPy: 858.691: delta=3.347
    BQ: 861.624, SciPy: 858.540: delta=3.084
    BQ: 860.816, SciPy: 857.106: delta=3.710
    BQ: 860.816, SciPy: 857.106: delta=3.710
    BQ: 859.677, SciPy: 856.050: delta=3.627
    BQ: 858.701, SciPy: 855.471: delta=3.230
    BQ: 858.701, SciPy: 855.471: delta=3.230
    BQ: 858.701, SciPy: 855.471: delta=3.230
    BQ: 858.701, SciPy: 855.471: delta=3.230
    BQ: 858.701, SciPy: 855.471: delta=3.230
    BQ: 858.701, SciPy: 855.471: delta=3.230
    BQ: 858.228, SciPy: 854.886: delta=3.342
    BQ: 857.528, SciPy: 854.297: delta=3.231
    BQ: 857.528, SciPy: 854.297: delta=3.231
    BQ: 857.528, SciPy: 854.297: delta=3.231
    BQ: 857.528, SciPy: 854.297: delta=3.231
    BQ: 857.040, SciPy: 853.335: delta=3.705
    BQ: 857.040, SciPy: 853.335: delta=3.705
    BQ: 857.040, SciPy: 853.335: delta=3.705
    BQ: 857.040, SciPy: 853.335: delta=3.705
    BQ: 856.605, SciPy: 852.863: delta=3.742
    BQ: 856.594, SciPy: 852.813: delta=3.781
    BQ: 854.883, SciPy: 851.658: delta=3.225
    BQ: 854.883, SciPy: 851.658: delta=3.225
    BQ: 854.883, SciPy: 851.658: delta=3.225
    BQ: 853.983, SciPy: 850.128: delta=3.855
    BQ: 853.405, SciPy: 849.711: delta=3.694
    BQ: 853.270, SciPy: 849.569: delta=3.701
    BQ: 853.270, SciPy: 849.569: delta=3.701
    BQ: 853.270, SciPy: 849.569: delta=3.701
    BQ: 853.270, SciPy: 849.569: delta=3.701
    BQ: 851.789, SciPy: 848.830: delta=2.959
    BQ: 851.535, SciPy: 847.906: delta=3.629
    BQ: 851.072, SciPy: 847.853: delta=3.219
    BQ: 851.072, SciPy: 847.853: delta=3.219
    BQ: 851.072, SciPy: 847.853: delta=3.219
    BQ: 849.903, SciPy: 846.683: delta=3.220
    BQ: 849.504, SciPy: 845.808: delta=3.696
    BQ: 847.267, SciPy: 844.053: delta=3.214
    BQ: 847.267, SciPy: 844.053: delta=3.214
    BQ: 847.267, SciPy: 844.053: delta=3.214
    BQ: 847.267, SciPy: 844.053: delta=3.214
    BQ: 847.267, SciPy: 844.053: delta=3.214
    BQ: 847.267, SciPy: 844.053: delta=3.214
    BQ: 847.267, SciPy: 844.053: delta=3.214
    BQ: 847.267, SciPy: 844.053: delta=3.214
    BQ: 846.026, SciPy: 842.348: delta=3.678
    BQ: 845.744, SciPy: 842.052: delta=3.692
    BQ: 844.640, SciPy: 841.433: delta=3.207
    BQ: 844.063, SciPy: 840.774: delta=3.289
    BQ: 843.240, SciPy: 839.601: delta=3.639
    BQ: 843.194, SciPy: 840.134: delta=3.060
    BQ: 842.022, SciPy: 838.822: delta=3.200
    BQ: 840.819, SciPy: 837.216: delta=3.603
    BQ: 840.291, SciPy: 837.117: delta=3.174
    BQ: 839.483, SciPy: 835.849: delta=3.634
    BQ: 839.483, SciPy: 835.849: delta=3.634
    BQ: 839.413, SciPy: 836.220: delta=3.193
    BQ: 839.245, SciPy: 835.397: delta=3.848
    BQ: 838.008, SciPy: 834.822: delta=3.186
    BQ: 837.950, SciPy: 834.915: delta=3.035
    BQ: 837.790, SciPy: 834.063: delta=3.727
    BQ: 837.777, SciPy: 834.034: delta=3.743
    BQ: 836.282, SciPy: 832.317: delta=3.965
    BQ: 835.531, SciPy: 832.613: delta=2.918
    BQ: 835.284, SciPy: 831.450: delta=3.834
    BQ: 835.074, SciPy: 831.249: delta=3.825
    BQ: 834.998, SciPy: 830.765: delta=4.233
    BQ: 834.931, SciPy: 831.274: delta=3.657
    BQ: 832.507, SciPy: 829.468: delta=3.039
    BQ: 832.042, SciPy: 828.192: delta=3.850
    BQ: 831.900, SciPy: 828.910: delta=2.990
    BQ: 830.302, SciPy: 826.568: delta=3.734
    BQ: 829.874, SciPy: 826.291: delta=3.583
    BQ: 829.874, SciPy: 826.291: delta=3.583
    BQ: 828.480, SciPy: 824.865: delta=3.615
    BQ: 826.903, SciPy: 823.725: delta=3.178
    BQ: 826.659, SciPy: 822.961: delta=3.698
    BQ: 826.572, SciPy: 822.843: delta=3.729
    BQ: 826.572, SciPy: 822.843: delta=3.729
    BQ: 826.572, SciPy: 822.843: delta=3.729
    BQ: 826.572, SciPy: 822.843: delta=3.729
    BQ: 826.572, SciPy: 822.843: delta=3.729
    BQ: 824.814, SciPy: 821.395: delta=3.419
    BQ: 824.172, SciPy: 820.749: delta=3.423
    BQ: 824.172, SciPy: 820.749: delta=3.423
    BQ: 824.172, SciPy: 820.749: delta=3.423
    BQ: 823.427, SciPy: 820.175: delta=3.252
    BQ: 822.436, SciPy: 819.181: delta=3.255
    BQ: 820.020, SciPy: 816.474: delta=3.546
    BQ: 819.133, SciPy: 815.428: delta=3.705
    BQ: 819.103, SciPy: 815.061: delta=4.042
    BQ: 818.083, SciPy: 814.220: delta=3.863
    BQ: 816.963, SciPy: 812.898: delta=4.065
    BQ: 815.477, SciPy: 811.743: delta=3.734
    BQ: 812.937, SciPy: 809.407: delta=3.530
    BQ: 812.907, SciPy: 809.292: delta=3.615
    BQ: 812.527, SciPy: 809.144: delta=3.383
    BQ: 811.866, SciPy: 808.318: delta=3.548
    BQ: 811.724, SciPy: 808.004: delta=3.720
    BQ: 811.700, SciPy: 807.988: delta=3.712
    BQ: 810.627, SciPy: 807.055: delta=3.572
    BQ: 809.600, SciPy: 806.080: delta=3.520
    BQ: 809.585, SciPy: 806.426: delta=3.159
    BQ: 809.489, SciPy: 806.116: delta=3.373
    BQ: 808.748, SciPy: 805.531: delta=3.217
    BQ: 808.023, SciPy: 804.339: delta=3.684
    BQ: 807.994, SciPy: 804.286: delta=3.708
    BQ: 807.994, SciPy: 804.286: delta=3.708
    BQ: 807.994, SciPy: 804.286: delta=3.708
    BQ: 807.994, SciPy: 804.286: delta=3.708
    BQ: 807.994, SciPy: 804.286: delta=3.708
    BQ: 807.994, SciPy: 804.286: delta=3.708
    BQ: 807.994, SciPy: 804.286: delta=3.708
    BQ: 807.994, SciPy: 804.286: delta=3.708
    BQ: 807.994, SciPy: 804.286: delta=3.708
    BQ: 807.994, SciPy: 804.286: delta=3.708
    BQ: 807.994, SciPy: 804.286: delta=3.708
    BQ: 807.994, SciPy: 804.286: delta=3.708
    BQ: 807.994, SciPy: 804.286: delta=3.708
    BQ: 807.994, SciPy: 804.286: delta=3.708
    BQ: 807.860, SciPy: 804.965: delta=2.895
    BQ: 807.172, SciPy: 803.811: delta=3.361
    BQ: 806.638, SciPy: 803.065: delta=3.573
    BQ: 805.134, SciPy: 801.608: delta=3.526
    BQ: 804.968, SciPy: 801.341: delta=3.627
    BQ: 804.437, SciPy: 800.899: delta=3.538
    BQ: 804.437, SciPy: 800.899: delta=3.538
    BQ: 804.320, SciPy: 800.640: delta=3.680
    BQ: 804.293, SciPy: 800.606: delta=3.687
    BQ: 804.293, SciPy: 800.606: delta=3.687
    BQ: 804.293, SciPy: 800.606: delta=3.687
    BQ: 804.293, SciPy: 800.606: delta=3.687
    BQ: 804.292, SciPy: 800.589: delta=3.703
    BQ: 804.292, SciPy: 800.589: delta=3.703
    BQ: 804.292, SciPy: 800.589: delta=3.703
    BQ: 804.284, SciPy: 800.589: delta=3.695
    BQ: 804.284, SciPy: 800.589: delta=3.695
    BQ: 804.108, SciPy: 800.564: delta=3.544
    BQ: 803.491, SciPy: 799.935: delta=3.556
    BQ: 803.491, SciPy: 799.935: delta=3.556
    BQ: 803.491, SciPy: 799.935: delta=3.556
    BQ: 803.491, SciPy: 799.935: delta=3.556
    BQ: 801.752, SciPy: 798.505: delta=3.247
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.426, SciPy: 797.905: delta=3.521
    BQ: 801.071, SciPy: 797.545: delta=3.526
    BQ: 801.071, SciPy: 797.545: delta=3.526
    BQ: 801.066, SciPy: 797.102: delta=3.964
    BQ: 800.730, SciPy: 797.197: delta=3.533
    BQ: 800.622, SciPy: 796.914: delta=3.708
    BQ: 800.622, SciPy: 796.914: delta=3.708
    BQ: 800.595, SciPy: 796.912: delta=3.683
    BQ: 800.595, SciPy: 796.912: delta=3.683
    BQ: 800.595, SciPy: 796.912: delta=3.683
    BQ: 800.595, SciPy: 796.912: delta=3.683
    BQ: 800.595, SciPy: 796.912: delta=3.683
    BQ: 800.595, SciPy: 796.912: delta=3.683
    BQ: 800.595, SciPy: 796.912: delta=3.683
    BQ: 800.300, SciPy: 797.179: delta=3.121
    BQ: 800.088, SciPy: 796.543: delta=3.545
    BQ: 798.090, SciPy: 794.580: delta=3.510
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.723, SciPy: 794.207: delta=3.516
    BQ: 797.636, SciPy: 794.539: delta=3.097
    BQ: 797.401, SciPy: 793.441: delta=3.960
    BQ: 797.369, SciPy: 793.847: delta=3.522
    BQ: 797.369, SciPy: 793.847: delta=3.522
    BQ: 796.927, SciPy: 793.257: delta=3.670
    BQ: 796.701, SciPy: 793.167: delta=3.534
    BQ: 794.771, SciPy: 791.271: delta=3.500
    BQ: 794.771, SciPy: 791.271: delta=3.500
    BQ: 794.771, SciPy: 791.271: delta=3.500
    BQ: 794.771, SciPy: 791.271: delta=3.500
    BQ: 794.392, SciPy: 790.886: delta=3.506
    BQ: 794.392, SciPy: 790.886: delta=3.506
    BQ: 794.392, SciPy: 790.886: delta=3.506
    BQ: 794.392, SciPy: 790.886: delta=3.506
    BQ: 794.392, SciPy: 790.886: delta=3.506
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 794.025, SciPy: 790.514: delta=3.511
    BQ: 793.922, SciPy: 790.830: delta=3.092
    BQ: 793.922, SciPy: 790.830: delta=3.092
    BQ: 793.922, SciPy: 790.830: delta=3.092
    BQ: 793.922, SciPy: 790.830: delta=3.092
    BQ: 793.672, SciPy: 790.155: delta=3.517
    BQ: 793.613, SciPy: 790.150: delta=3.463
    BQ: 793.332, SciPy: 789.809: delta=3.523
    BQ: 793.216, SciPy: 789.526: delta=3.690
    BQ: 791.872, SciPy: 788.388: delta=3.484
    BQ: 791.618, SciPy: 788.443: delta=3.175
    BQ: 791.447, SciPy: 788.362: delta=3.085
    BQ: 791.311, SciPy: 788.483: delta=2.828
    BQ: 791.077, SciPy: 787.582: delta=3.495
    BQ: 791.077, SciPy: 787.582: delta=3.495
    BQ: 791.077, SciPy: 787.582: delta=3.495
    BQ: 791.077, SciPy: 787.582: delta=3.495
    BQ: 791.077, SciPy: 787.582: delta=3.495
    BQ: 791.077, SciPy: 787.582: delta=3.495
    BQ: 791.077, SciPy: 787.582: delta=3.495
    BQ: 791.077, SciPy: 787.582: delta=3.495
    BQ: 790.698, SciPy: 787.198: delta=3.500
    BQ: 790.698, SciPy: 787.198: delta=3.500
    BQ: 790.698, SciPy: 787.198: delta=3.500
    BQ: 790.698, SciPy: 787.198: delta=3.500
    BQ: 790.698, SciPy: 787.198: delta=3.500
    BQ: 790.698, SciPy: 787.198: delta=3.500
    BQ: 790.333, SciPy: 786.826: delta=3.507
    BQ: 790.333, SciPy: 786.826: delta=3.507
    BQ: 790.333, SciPy: 786.826: delta=3.507
    BQ: 790.083, SciPy: 786.131: delta=3.952
    BQ: 789.438, SciPy: 786.614: delta=2.824
    BQ: 787.891, SciPy: 784.336: delta=3.555
    BQ: 787.779, SciPy: 784.294: delta=3.485
    BQ: 787.388, SciPy: 783.898: delta=3.490
    BQ: 787.388, SciPy: 783.898: delta=3.490
    BQ: 787.388, SciPy: 783.898: delta=3.490
    BQ: 787.388, SciPy: 783.898: delta=3.490
    BQ: 787.388, SciPy: 783.898: delta=3.490
    BQ: 787.388, SciPy: 783.898: delta=3.490
    BQ: 787.388, SciPy: 783.898: delta=3.490
    BQ: 787.388, SciPy: 783.898: delta=3.490
    BQ: 786.646, SciPy: 783.144: delta=3.502
    BQ: 786.646, SciPy: 783.144: delta=3.502
    BQ: 786.646, SciPy: 783.144: delta=3.502
    BQ: 786.646, SciPy: 783.144: delta=3.502
    BQ: 786.646, SciPy: 783.144: delta=3.502
    BQ: 786.646, SciPy: 783.144: delta=3.502
    BQ: 786.646, SciPy: 783.144: delta=3.502
    BQ: 786.503, SciPy: 782.899: delta=3.604
    BQ: 785.930, SciPy: 782.232: delta=3.698
    BQ: 785.705, SciPy: 782.889: delta=2.816
    BQ: 785.705, SciPy: 782.889: delta=2.816
    BQ: 784.699, SciPy: 781.663: delta=3.036
    BQ: 784.497, SciPy: 781.023: delta=3.474
    BQ: 784.095, SciPy: 780.615: delta=3.480
    BQ: 782.964, SciPy: 779.467: delta=3.497
    BQ: 782.964, SciPy: 779.467: delta=3.497
    BQ: 782.964, SciPy: 779.467: delta=3.497
    BQ: 782.566, SciPy: 778.952: delta=3.614
    BQ: 782.182, SciPy: 778.506: delta=3.676
    BQ: 781.659, SciPy: 778.201: delta=3.458
    BQ: 780.818, SciPy: 777.349: delta=3.469
    BQ: 780.818, SciPy: 777.349: delta=3.469
    BQ: 780.124, SciPy: 776.952: delta=3.172
    BQ: 779.735, SciPy: 776.884: delta=2.851
    BQ: 778.540, SciPy: 775.200: delta=3.340
    BQ: 778.295, SciPy: 775.495: delta=2.800
    BQ: 778.275, SciPy: 774.765: delta=3.510
    BQ: 777.910, SciPy: 775.130: delta=2.780
    BQ: 777.557, SciPy: 774.098: delta=3.459
    BQ: 777.143, SciPy: 773.679: delta=3.464
    BQ: 776.354, SciPy: 772.879: delta=3.475
    BQ: 776.319, SciPy: 773.520: delta=2.799
    BQ: 776.319, SciPy: 773.520: delta=2.799
    BQ: 776.319, SciPy: 773.520: delta=2.799
    BQ: 775.978, SciPy: 772.497: delta=3.481
    BQ: 775.615, SciPy: 772.128: delta=3.487
    BQ: 775.483, SciPy: 771.894: delta=3.589
    BQ: 775.483, SciPy: 771.894: delta=3.589
    BQ: 775.203, SciPy: 772.494: delta=2.709
    BQ: 775.005, SciPy: 771.607: delta=3.398
    BQ: 773.490, SciPy: 770.527: delta=2.963
    BQ: 773.138, SciPy: 769.178: delta=3.960
    BQ: 772.951, SciPy: 769.534: delta=3.417
    BQ: 772.551, SciPy: 768.788: delta=3.763
    BQ: 772.138, SciPy: 768.567: delta=3.571
    BQ: 771.879, SciPy: 768.491: delta=3.388
    BQ: 771.820, SciPy: 768.235: delta=3.585
    BQ: 771.820, SciPy: 768.235: delta=3.585
    BQ: 771.791, SciPy: 768.676: delta=3.115
    BQ: 771.600, SciPy: 768.112: delta=3.488
    BQ: 771.444, SciPy: 768.339: delta=3.105
    BQ: 771.084, SciPy: 767.647: delta=3.437
    BQ: 770.648, SciPy: 767.204: delta=3.444
    BQ: 769.559, SciPy: 766.531: delta=3.028
    BQ: 769.482, SciPy: 766.286: delta=3.196
    BQ: 769.226, SciPy: 766.117: delta=3.109
    BQ: 768.161, SciPy: 764.581: delta=3.580
    BQ: 768.026, SciPy: 764.439: delta=3.587
    BQ: 767.819, SciPy: 764.958: delta=2.861
    BQ: 767.316, SciPy: 764.540: delta=2.776
    BQ: 766.670, SciPy: 763.569: delta=3.101
    BQ: 766.670, SciPy: 763.569: delta=3.101
    BQ: 766.670, SciPy: 763.569: delta=3.101
    BQ: 766.670, SciPy: 763.569: delta=3.101
    BQ: 766.599, SciPy: 763.295: delta=3.304
    BQ: 766.129, SciPy: 762.608: delta=3.521
    BQ: 765.147, SciPy: 761.735: delta=3.412
    BQ: 765.081, SciPy: 761.793: delta=3.288
    BQ: 765.041, SciPy: 762.246: delta=2.795
    BQ: 764.787, SciPy: 761.771: delta=3.016
    BQ: 764.768, SciPy: 761.922: delta=2.846
    BQ: 764.284, SciPy: 760.806: delta=3.478
    BQ: 763.943, SciPy: 760.942: delta=3.001
    BQ: 763.843, SciPy: 759.696: delta=4.147
    BQ: 763.769, SciPy: 760.341: delta=3.428
    BQ: 763.366, SciPy: 760.343: delta=3.023
    BQ: 763.219, SciPy: 760.050: delta=3.169
    BQ: 762.977, SciPy: 759.585: delta=3.392
    BQ: 761.326, SciPy: 758.538: delta=2.788
    BQ: 761.131, SciPy: 758.122: delta=3.009
    BQ: 761.131, SciPy: 758.122: delta=3.009
    BQ: 760.725, SciPy: 757.147: delta=3.578
    BQ: 760.566, SciPy: 757.149: delta=3.417
    BQ: 760.474, SciPy: 757.385: delta=3.089
    BQ: 760.048, SciPy: 757.045: delta=3.003
    BQ: 760.006, SciPy: 757.273: delta=2.733
    BQ: 759.733, SciPy: 756.891: delta=2.842
    BQ: 759.457, SciPy: 756.673: delta=2.784
    BQ: 759.457, SciPy: 756.673: delta=2.784
    BQ: 758.852, SciPy: 755.413: delta=3.439
    BQ: 758.830, SciPy: 755.319: delta=3.511
    BQ: 758.830, SciPy: 755.319: delta=3.511
    BQ: 758.830, SciPy: 755.319: delta=3.511
    BQ: 758.830, SciPy: 755.319: delta=3.511
    BQ: 758.830, SciPy: 755.319: delta=3.511
    BQ: 758.830, SciPy: 755.319: delta=3.511
    BQ: 758.830, SciPy: 755.319: delta=3.511
    BQ: 758.830, SciPy: 755.319: delta=3.511
    BQ: 758.830, SciPy: 755.319: delta=3.511
    BQ: 758.251, SciPy: 755.268: delta=2.983
    BQ: 758.161, SciPy: 755.403: delta=2.758
    BQ: 758.045, SciPy: 755.285: delta=2.760
    BQ: 758.045, SciPy: 755.285: delta=2.760
    BQ: 757.336, SciPy: 753.874: delta=3.462
    BQ: 755.733, SciPy: 752.957: delta=2.776
    BQ: 755.733, SciPy: 752.957: delta=2.776
    BQ: 755.454, SciPy: 751.955: delta=3.499
    BQ: 754.698, SciPy: 751.179: delta=3.519
    BQ: 754.698, SciPy: 751.179: delta=3.519
    BQ: 754.698, SciPy: 751.179: delta=3.519
    BQ: 754.698, SciPy: 751.179: delta=3.519
    BQ: 754.698, SciPy: 751.179: delta=3.519
    BQ: 753.726, SciPy: 750.172: delta=3.554
    BQ: 753.049, SciPy: 749.392: delta=3.657
    BQ: 752.341, SciPy: 749.326: delta=3.015
    BQ: 752.155, SciPy: 748.365: delta=3.790
    BQ: 752.096, SciPy: 748.608: delta=3.488
    BQ: 752.096, SciPy: 748.608: delta=3.488
    BQ: 752.096, SciPy: 748.608: delta=3.488
    BQ: 752.028, SciPy: 749.261: delta=2.767
    BQ: 752.022, SciPy: 748.145: delta=3.877
    BQ: 751.418, SciPy: 747.480: delta=3.938
    BQ: 751.345, SciPy: 748.190: delta=3.155
    BQ: 751.299, SciPy: 747.792: delta=3.507
    BQ: 751.299, SciPy: 747.792: delta=3.507
    BQ: 750.838, SciPy: 747.318: delta=3.520
    BQ: 749.369, SciPy: 745.724: delta=3.645
    BQ: 748.756, SciPy: 745.280: delta=3.476
    BQ: 748.342, SciPy: 745.583: delta=2.759
    BQ: 748.193, SciPy: 745.393: delta=2.800
    BQ: 747.919, SciPy: 744.423: delta=3.496
    BQ: 747.919, SciPy: 744.423: delta=3.496
    BQ: 747.919, SciPy: 744.423: delta=3.496
    BQ: 747.887, SciPy: 744.025: delta=3.862
    BQ: 747.812, SciPy: 743.878: delta=3.934
    BQ: 747.021, SciPy: 744.125: delta=2.896
    BQ: 746.562, SciPy: 743.811: delta=2.751
    BQ: 746.562, SciPy: 743.811: delta=2.751
    BQ: 746.047, SciPy: 742.232: delta=3.815
    BQ: 746.047, SciPy: 742.232: delta=3.815
    BQ: 745.710, SciPy: 742.079: delta=3.631
    BQ: 745.128, SciPy: 741.657: delta=3.471
    BQ: 744.897, SciPy: 740.953: delta=3.944
    BQ: 744.836, SciPy: 741.358: delta=3.478
    BQ: 744.722, SciPy: 741.309: delta=3.413
    BQ: 744.557, SciPy: 741.073: delta=3.484
    BQ: 743.841, SciPy: 740.504: delta=3.337
    BQ: 743.284, SciPy: 739.943: delta=3.341
    BQ: 743.254, SciPy: 739.536: delta=3.718
    BQ: 743.254, SciPy: 739.536: delta=3.718
    BQ: 743.254, SciPy: 739.536: delta=3.718
    BQ: 743.254, SciPy: 739.536: delta=3.718
    BQ: 742.807, SciPy: 739.365: delta=3.442
    BQ: 742.449, SciPy: 739.741: delta=2.708
    BQ: 742.074, SciPy: 738.456: delta=3.618
    BQ: 742.074, SciPy: 738.456: delta=3.618
    BQ: 742.074, SciPy: 738.456: delta=3.618
    BQ: 741.693, SciPy: 737.626: delta=4.067
    BQ: 741.313, SciPy: 738.486: delta=2.827
    BQ: 741.100, SciPy: 737.692: delta=3.408
    BQ: 740.909, SciPy: 738.087: delta=2.822
    BQ: 740.787, SciPy: 737.491: delta=3.296
    BQ: 740.672, SciPy: 737.629: delta=3.043
    BQ: 739.561, SciPy: 736.041: delta=3.520
    BQ: 738.843, SciPy: 735.400: delta=3.443
    BQ: 738.514, SciPy: 734.925: delta=3.589
    BQ: 738.457, SciPy: 734.844: delta=3.613
    BQ: 738.425, SciPy: 735.379: delta=3.046
    BQ: 737.903, SciPy: 734.505: delta=3.398
    BQ: 736.359, SciPy: 732.384: delta=3.975
    BQ: 735.930, SciPy: 732.504: delta=3.426
    BQ: 735.666, SciPy: 732.850: delta=2.816
    BQ: 735.666, SciPy: 732.850: delta=2.816
    BQ: 734.976, SciPy: 731.635: delta=3.341
    BQ: 734.268, SciPy: 731.569: delta=2.699
    BQ: 733.592, SciPy: 730.889: delta=2.703
    BQ: 733.162, SciPy: 729.426: delta=3.736
    BQ: 733.162, SciPy: 729.426: delta=3.736
    BQ: 732.895, SciPy: 729.169: delta=3.726
    BQ: 732.770, SciPy: 728.864: delta=3.906
    BQ: 732.425, SciPy: 728.719: delta=3.706
    BQ: 732.138, SciPy: 728.245: delta=3.893
    BQ: 732.138, SciPy: 728.245: delta=3.893
    BQ: 730.560, SciPy: 727.878: delta=2.682
    BQ: 730.560, SciPy: 727.878: delta=2.682
    BQ: 730.378, SciPy: 726.920: delta=3.458
    BQ: 730.131, SciPy: 727.436: delta=2.695
    BQ: 730.115, SciPy: 727.079: delta=3.036
    BQ: 729.295, SciPy: 725.573: delta=3.722
    BQ: 729.295, SciPy: 725.573: delta=3.722
    BQ: 729.295, SciPy: 725.573: delta=3.722
    BQ: 729.049, SciPy: 725.338: delta=3.711
    BQ: 728.924, SciPy: 725.425: delta=3.499
    BQ: 728.878, SciPy: 726.199: delta=2.679
    BQ: 728.824, SciPy: 725.123: delta=3.701
    BQ: 728.279, SciPy: 725.512: delta=2.767
    BQ: 727.635, SciPy: 724.036: delta=3.599
    BQ: 726.898, SciPy: 724.193: delta=2.705
    BQ: 726.897, SciPy: 723.137: delta=3.760
    BQ: 726.384, SciPy: 723.661: delta=2.723
    BQ: 725.987, SciPy: 722.514: delta=3.473
    BQ: 725.228, SciPy: 721.531: delta=3.697
    BQ: 725.023, SciPy: 721.337: delta=3.686
    BQ: 724.868, SciPy: 721.066: delta=3.802
    BQ: 724.706, SciPy: 721.386: delta=3.320
    BQ: 723.701, SciPy: 720.954: delta=2.747
    BQ: 723.372, SciPy: 720.674: delta=2.698
    BQ: 723.260, SciPy: 720.579: delta=2.681
    BQ: 723.260, SciPy: 720.579: delta=2.681
    BQ: 723.029, SciPy: 720.348: delta=2.681
    BQ: 722.663, SciPy: 718.929: delta=3.734
    BQ: 721.733, SciPy: 717.923: delta=3.810
    BQ: 721.636, SciPy: 717.944: delta=3.692
    BQ: 721.636, SciPy: 717.944: delta=3.692
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.615, SciPy: 718.922: delta=2.693
    BQ: 721.492, SciPy: 718.797: delta=2.695
    BQ: 721.317, SciPy: 718.639: delta=2.678
    BQ: 721.317, SciPy: 718.639: delta=2.678
    BQ: 720.882, SciPy: 718.201: delta=2.681
    BQ: 720.882, SciPy: 718.201: delta=2.681
    BQ: 720.784, SciPy: 716.916: delta=3.868
    BQ: 720.768, SciPy: 717.407: delta=3.361
    BQ: 720.546, SciPy: 716.930: delta=3.616
    BQ: 720.460, SciPy: 716.862: delta=3.598
    BQ: 720.181, SciPy: 716.751: delta=3.430
    BQ: 719.863, SciPy: 717.173: delta=2.690
    BQ: 719.863, SciPy: 717.173: delta=2.690
    BQ: 719.863, SciPy: 717.173: delta=2.690
    BQ: 719.863, SciPy: 717.173: delta=2.690
    BQ: 719.863, SciPy: 717.173: delta=2.690
    BQ: 719.863, SciPy: 717.173: delta=2.690
    BQ: 719.863, SciPy: 717.173: delta=2.690
    BQ: 719.848, SciPy: 717.175: delta=2.673
    BQ: 719.720, SciPy: 715.969: delta=3.751
    BQ: 719.608, SciPy: 716.934: delta=2.674
    BQ: 719.608, SciPy: 716.934: delta=2.674
    BQ: 719.608, SciPy: 716.934: delta=2.674
    BQ: 719.049, SciPy: 715.592: delta=3.457
    BQ: 718.115, SciPy: 715.429: delta=2.686
    BQ: 718.115, SciPy: 715.429: delta=2.686
    BQ: 718.115, SciPy: 715.429: delta=2.686
    BQ: 718.049, SciPy: 714.361: delta=3.688
    BQ: 717.982, SciPy: 715.295: delta=2.687
    BQ: 717.982, SciPy: 715.295: delta=2.687
    BQ: 717.982, SciPy: 715.295: delta=2.687
    BQ: 717.982, SciPy: 715.295: delta=2.687
    BQ: 717.948, SciPy: 715.042: delta=2.906
    BQ: 717.904, SciPy: 715.234: delta=2.670
    BQ: 717.904, SciPy: 715.234: delta=2.670
    BQ: 717.844, SciPy: 714.166: delta=3.678
    BQ: 717.671, SciPy: 715.000: delta=2.671
    BQ: 717.671, SciPy: 715.000: delta=2.671
    BQ: 717.649, SciPy: 714.863: delta=2.786
    BQ: 717.649, SciPy: 714.863: delta=2.786
    BQ: 717.649, SciPy: 714.863: delta=2.786
    BQ: 717.582, SciPy: 714.175: delta=3.407
    BQ: 717.450, SciPy: 714.777: delta=2.673
    BQ: 717.450, SciPy: 714.777: delta=2.673
    BQ: 716.760, SciPy: 713.479: delta=3.281
    BQ: 716.387, SciPy: 713.403: delta=2.984
    BQ: 716.371, SciPy: 713.690: delta=2.681
    BQ: 716.371, SciPy: 713.690: delta=2.681
    BQ: 716.204, SciPy: 713.537: delta=2.667
    BQ: 716.204, SciPy: 713.537: delta=2.667
    BQ: 716.204, SciPy: 713.537: delta=2.667
    BQ: 716.107, SciPy: 713.422: delta=2.685
    BQ: 715.970, SciPy: 713.049: delta=2.921
    BQ: 715.966, SciPy: 713.298: delta=2.668
    BQ: 714.980, SciPy: 711.515: delta=3.465
    BQ: 714.508, SciPy: 711.845: delta=2.663
    BQ: 714.508, SciPy: 711.845: delta=2.663
    BQ: 714.508, SciPy: 711.845: delta=2.663
    BQ: 714.508, SciPy: 711.845: delta=2.663
    BQ: 714.508, SciPy: 711.845: delta=2.663
    BQ: 714.508, SciPy: 711.845: delta=2.663
    BQ: 714.508, SciPy: 711.845: delta=2.663
    BQ: 714.508, SciPy: 711.845: delta=2.663
    BQ: 714.508, SciPy: 711.845: delta=2.663
    BQ: 714.358, SciPy: 711.677: delta=2.681
    BQ: 714.358, SciPy: 711.677: delta=2.681
    BQ: 714.265, SciPy: 711.601: delta=2.664
    BQ: 714.261, SciPy: 710.587: delta=3.674
    BQ: 714.261, SciPy: 710.587: delta=3.674
    BQ: 714.132, SciPy: 711.447: delta=2.685
    BQ: 714.076, SciPy: 710.413: delta=3.663
    BQ: 713.374, SciPy: 709.767: delta=3.607
    BQ: 713.374, SciPy: 709.767: delta=3.607
    BQ: 713.374, SciPy: 709.767: delta=3.607
    BQ: 713.323, SciPy: 709.724: delta=3.599
    BQ: 713.273, SciPy: 709.701: delta=3.572
    BQ: 713.024, SciPy: 709.043: delta=3.981
    BQ: 712.916, SciPy: 709.162: delta=3.754
    BQ: 712.612, SciPy: 709.935: delta=2.677
    BQ: 712.333, SciPy: 709.671: delta=2.662
    BQ: 712.319, SciPy: 708.951: delta=3.368
    BQ: 712.037, SciPy: 709.329: delta=2.708
    BQ: 711.932, SciPy: 709.003: delta=2.929
    BQ: 711.774, SciPy: 708.988: delta=2.786
    BQ: 711.577, SciPy: 708.624: delta=2.953
    BQ: 711.127, SciPy: 708.472: delta=2.655
    BQ: 710.887, SciPy: 707.208: delta=3.679
    BQ: 710.887, SciPy: 707.208: delta=3.679
    BQ: 710.876, SciPy: 708.219: delta=2.657
    BQ: 710.682, SciPy: 707.013: delta=3.669
    BQ: 710.682, SciPy: 707.013: delta=3.669
    BQ: 710.682, SciPy: 707.013: delta=3.669
    BQ: 710.623, SciPy: 707.983: delta=2.640
    BQ: 710.581, SciPy: 706.645: delta=3.936
    BQ: 709.711, SciPy: 707.060: delta=2.651
    BQ: 709.524, SciPy: 706.309: delta=3.215
    BQ: 709.107, SciPy: 705.200: delta=3.907
    BQ: 709.107, SciPy: 705.200: delta=3.907
    BQ: 708.987, SciPy: 705.248: delta=3.739
    BQ: 708.986, SciPy: 706.350: delta=2.636
    BQ: 708.078, SciPy: 704.880: delta=3.198
    BQ: 707.851, SciPy: 704.598: delta=3.253
    BQ: 707.809, SciPy: 703.862: delta=3.947
    BQ: 707.502, SciPy: 704.853: delta=2.649
    BQ: 707.371, SciPy: 704.041: delta=3.330
    BQ: 707.333, SciPy: 704.417: delta=2.916
    BQ: 707.313, SciPy: 703.639: delta=3.674
    BQ: 707.295, SciPy: 704.172: delta=3.123
    BQ: 706.983, SciPy: 703.202: delta=3.781
    BQ: 706.923, SciPy: 703.269: delta=3.654
    BQ: 706.875, SciPy: 703.909: delta=2.966
    BQ: 706.759, SciPy: 703.114: delta=3.645
    BQ: 706.759, SciPy: 703.114: delta=3.645
    BQ: 706.372, SciPy: 703.758: delta=2.614
    BQ: 706.170, SciPy: 702.581: delta=3.589
    BQ: 705.831, SciPy: 703.171: delta=2.660
    BQ: 705.550, SciPy: 701.647: delta=3.903
    BQ: 705.550, SciPy: 701.647: delta=3.903
    BQ: 705.273, SciPy: 701.850: delta=3.423
    BQ: 705.046, SciPy: 702.135: delta=2.911
    BQ: 705.046, SciPy: 702.135: delta=2.911
    BQ: 705.046, SciPy: 702.135: delta=2.911
    BQ: 705.046, SciPy: 702.135: delta=2.911
    BQ: 704.734, SciPy: 702.048: delta=2.686
    BQ: 703.690, SciPy: 700.670: delta=3.020
    BQ: 703.424, SciPy: 700.759: delta=2.665
    BQ: 703.322, SciPy: 700.655: delta=2.667
    BQ: 703.159, SciPy: 699.982: delta=3.177
    BQ: 703.044, SciPy: 699.413: delta=3.631
    BQ: 702.966, SciPy: 700.318: delta=2.648
    BQ: 702.652, SciPy: 699.059: delta=3.593
    BQ: 702.622, SciPy: 699.393: delta=3.229
    BQ: 702.515, SciPy: 699.734: delta=2.781
    BQ: 702.294, SciPy: 699.324: delta=2.970
    BQ: 702.029, SciPy: 698.682: delta=3.347
    BQ: 701.997, SciPy: 698.097: delta=3.900
    BQ: 701.997, SciPy: 698.097: delta=3.900
    BQ: 701.997, SciPy: 698.097: delta=3.900
    BQ: 701.996, SciPy: 698.585: delta=3.411
    BQ: 701.983, SciPy: 698.749: delta=3.234
    BQ: 701.716, SciPy: 699.075: delta=2.641
    BQ: 701.544, SciPy: 698.685: delta=2.859
    BQ: 700.122, SciPy: 697.521: delta=2.601
    BQ: 699.821, SciPy: 695.726: delta=4.095
    BQ: 699.789, SciPy: 696.143: delta=3.646
    BQ: 699.785, SciPy: 696.461: delta=3.324
    BQ: 699.722, SciPy: 696.340: delta=3.382
    BQ: 699.624, SciPy: 695.988: delta=3.636
    BQ: 699.370, SciPy: 696.696: delta=2.674
    BQ: 699.070, SciPy: 695.846: delta=3.224
    BQ: 699.070, SciPy: 695.846: delta=3.224
    BQ: 698.935, SciPy: 695.680: delta=3.255
    BQ: 698.697, SciPy: 694.567: delta=4.130
    BQ: 697.958, SciPy: 694.244: delta=3.714
    BQ: 697.407, SciPy: 694.477: delta=2.930
    BQ: 697.090, SciPy: 693.408: delta=3.682
    BQ: 697.089, SciPy: 693.221: delta=3.868
    BQ: 696.860, SciPy: 693.801: delta=3.059
    BQ: 696.755, SciPy: 694.054: delta=2.701
    BQ: 696.228, SciPy: 692.587: delta=3.641
    BQ: 696.228, SciPy: 692.587: delta=3.641
    BQ: 696.211, SciPy: 693.381: delta=2.830
    BQ: 695.918, SciPy: 692.297: delta=3.621
    BQ: 695.470, SciPy: 692.878: delta=2.592
    BQ: 695.470, SciPy: 692.878: delta=2.592
    BQ: 695.387, SciPy: 692.137: delta=3.250
    BQ: 695.182, SciPy: 691.055: delta=4.127
    BQ: 695.182, SciPy: 691.055: delta=4.127
    BQ: 694.329, SciPy: 691.677: delta=2.652
    BQ: 694.104, SciPy: 691.144: delta=2.960
    BQ: 693.926, SciPy: 691.337: delta=2.589
    BQ: 693.802, SciPy: 690.113: delta=3.689
    BQ: 693.168, SciPy: 689.988: delta=3.180
    BQ: 692.432, SciPy: 689.193: delta=3.239
    BQ: 692.017, SciPy: 689.415: delta=2.602
    BQ: 691.920, SciPy: 688.180: delta=3.740
    BQ: 691.844, SciPy: 688.600: delta=3.244
    BQ: 691.671, SciPy: 687.547: delta=4.124
    BQ: 691.362, SciPy: 687.473: delta=3.889
    BQ: 691.063, SciPy: 687.896: delta=3.167
    BQ: 690.910, SciPy: 686.966: delta=3.944
    BQ: 690.723, SciPy: 688.099: delta=2.624
    BQ: 690.363, SciPy: 687.445: delta=2.918
    BQ: 690.257, SciPy: 686.834: delta=3.423
    BQ: 690.001, SciPy: 686.141: delta=3.860
    BQ: 689.305, SciPy: 685.664: delta=3.641
    BQ: 689.305, SciPy: 685.664: delta=3.641
    BQ: 689.221, SciPy: 686.284: delta=2.937
    BQ: 689.221, SciPy: 686.284: delta=2.937
    BQ: 689.221, SciPy: 686.284: delta=2.937
    BQ: 689.121, SciPy: 685.489: delta=3.632
    BQ: 688.577, SciPy: 685.611: delta=2.966
    BQ: 688.530, SciPy: 685.258: delta=3.272
    BQ: 688.421, SciPy: 684.846: delta=3.575
    BQ: 688.233, SciPy: 685.629: delta=2.604
    BQ: 688.205, SciPy: 685.626: delta=2.579
    BQ: 687.825, SciPy: 683.940: delta=3.885
    BQ: 687.783, SciPy: 685.207: delta=2.576
    BQ: 687.415, SciPy: 684.084: delta=3.331
    BQ: 687.179, SciPy: 684.596: delta=2.583
    BQ: 687.033, SciPy: 684.173: delta=2.860
    BQ: 687.033, SciPy: 684.173: delta=2.860
    BQ: 687.017, SciPy: 683.729: delta=3.288
    BQ: 686.623, SciPy: 683.281: delta=3.342
    BQ: 686.588, SciPy: 682.808: delta=3.780
    BQ: 686.463, SciPy: 682.607: delta=3.856
    BQ: 685.960, SciPy: 683.358: delta=2.602
    BQ: 685.471, SciPy: 682.322: delta=3.149
    BQ: 685.471, SciPy: 682.322: delta=3.149
    BQ: 685.471, SciPy: 682.322: delta=3.149
    BQ: 685.471, SciPy: 682.322: delta=3.149
    BQ: 685.471, SciPy: 682.322: delta=3.149
    BQ: 685.471, SciPy: 682.322: delta=3.149
    BQ: 685.418, SciPy: 682.752: delta=2.666
    BQ: 685.153, SciPy: 682.501: delta=2.652
    BQ: 684.052, SciPy: 681.480: delta=2.572
    BQ: 684.051, SciPy: 681.354: delta=2.697
    BQ: 684.051, SciPy: 681.354: delta=2.697
    BQ: 684.009, SciPy: 680.853: delta=3.156
    BQ: 683.596, SciPy: 679.729: delta=3.867
    BQ: 683.528, SciPy: 680.956: delta=2.572
    BQ: 683.484, SciPy: 680.202: delta=3.282
    BQ: 682.710, SciPy: 679.367: delta=3.343
    BQ: 682.645, SciPy: 679.929: delta=2.716
    BQ: 682.517, SciPy: 679.896: delta=2.621
    BQ: 682.511, SciPy: 679.939: delta=2.572
    BQ: 682.146, SciPy: 678.697: delta=3.449
    BQ: 681.951, SciPy: 678.807: delta=3.144
    BQ: 681.951, SciPy: 678.807: delta=3.144
    BQ: 681.951, SciPy: 678.807: delta=3.144
    BQ: 681.951, SciPy: 678.807: delta=3.144
    BQ: 681.562, SciPy: 678.607: delta=2.955
    BQ: 681.562, SciPy: 678.607: delta=2.955
    BQ: 681.562, SciPy: 678.607: delta=2.955
    BQ: 681.232, SciPy: 677.694: delta=3.538
    BQ: 681.181, SciPy: 677.232: delta=3.949
    BQ: 680.962, SciPy: 678.393: delta=2.569
    BQ: 680.953, SciPy: 677.686: delta=3.267
    BQ: 680.798, SciPy: 678.104: delta=2.694
    BQ: 680.491, SciPy: 677.340: delta=3.151
    BQ: 680.421, SciPy: 677.033: delta=3.388
    BQ: 680.144, SciPy: 677.558: delta=2.586
    BQ: 680.144, SciPy: 677.558: delta=2.586
    BQ: 679.993, SciPy: 677.431: delta=2.562
    BQ: 678.899, SciPy: 676.317: delta=2.582
    BQ: 678.899, SciPy: 676.317: delta=2.582
    BQ: 678.666, SciPy: 675.755: delta=2.911
    BQ: 678.494, SciPy: 674.876: delta=3.618
    BQ: 677.920, SciPy: 675.360: delta=2.560
    BQ: 677.697, SciPy: 674.173: delta=3.524
    BQ: 677.480, SciPy: 674.913: delta=2.567
    BQ: 677.430, SciPy: 674.169: delta=3.261
    BQ: 677.430, SciPy: 674.169: delta=3.261
    BQ: 677.179, SciPy: 673.803: delta=3.376
    BQ: 677.037, SciPy: 674.028: delta=3.009
    BQ: 676.891, SciPy: 673.509: delta=3.382
    BQ: 676.557, SciPy: 673.878: delta=2.679
    BQ: 676.557, SciPy: 673.878: delta=2.679
    BQ: 676.421, SciPy: 673.295: delta=3.126
    BQ: 676.401, SciPy: 673.759: delta=2.642
    BQ: 676.325, SciPy: 673.151: delta=3.174
    BQ: 675.705, SciPy: 673.131: delta=2.574
    BQ: 675.671, SciPy: 672.541: delta=3.130
    BQ: 675.055, SciPy: 671.312: delta=3.743
    BQ: 674.980, SciPy: 671.992: delta=2.988
    BQ: 673.675, SciPy: 670.558: delta=3.117
    BQ: 673.654, SciPy: 670.284: delta=3.370
    BQ: 672.917, SciPy: 669.797: delta=3.120
    BQ: 672.917, SciPy: 669.797: delta=3.120
    BQ: 672.917, SciPy: 669.797: delta=3.120
    BQ: 672.837, SciPy: 670.212: delta=2.625
    BQ: 672.149, SciPy: 669.581: delta=2.568
    BQ: 671.715, SciPy: 668.611: delta=3.104
    BQ: 671.439, SciPy: 668.198: delta=3.241
    BQ: 671.439, SciPy: 668.198: delta=3.241
    BQ: 670.942, SciPy: 667.834: delta=3.108
    BQ: 670.694, SciPy: 667.196: delta=3.498
    BQ: 670.561, SciPy: 667.996: delta=2.565
    BQ: 670.436, SciPy: 667.077: delta=3.359
    BQ: 670.198, SciPy: 666.333: delta=3.865
    BQ: 670.198, SciPy: 666.333: delta=3.865
    BQ: 670.191, SciPy: 667.625: delta=2.566
    BQ: 670.176, SciPy: 667.065: delta=3.111
    BQ: 669.419, SciPy: 666.305: delta=3.114
    BQ: 668.691, SciPy: 666.172: delta=2.519
    BQ: 668.602, SciPy: 666.040: delta=2.562
    BQ: 668.594, SciPy: 665.565: delta=3.029
    BQ: 668.221, SciPy: 665.122: delta=3.099
    BQ: 667.612, SciPy: 664.680: delta=2.932
    BQ: 667.213, SciPy: 663.667: delta=3.546
    BQ: 667.098, SciPy: 664.014: delta=3.084
    BQ: 666.966, SciPy: 664.237: delta=2.729
    BQ: 666.393, SciPy: 663.143: delta=3.250
    BQ: 666.393, SciPy: 663.143: delta=3.250
    BQ: 665.928, SciPy: 662.818: delta=3.110
    BQ: 665.877, SciPy: 662.571: delta=3.306
    BQ: 665.437, SciPy: 662.882: delta=2.555
    BQ: 664.642, SciPy: 662.092: delta=2.550
    BQ: 664.263, SciPy: 661.423: delta=2.840
    BQ: 663.859, SciPy: 661.308: delta=2.551
    BQ: 663.655, SciPy: 660.168: delta=3.487
    BQ: 663.600, SciPy: 660.095: delta=3.505
    BQ: 663.442, SciPy: 660.635: delta=2.807
    BQ: 663.413, SciPy: 660.064: delta=3.349
    BQ: 663.024, SciPy: 660.317: delta=2.707
    BQ: 662.442, SciPy: 659.338: delta=3.104
    BQ: 661.513, SciPy: 658.794: delta=2.719
    BQ: 661.393, SciPy: 658.206: delta=3.187
    BQ: 661.125, SciPy: 658.398: delta=2.727
    BQ: 661.052, SciPy: 658.521: delta=2.531


## Analyzing the GWAS results

First, how many statistically significant variant positions did we find?


    %%bq_sql 
    SELECT COUNT(*) AS num_significant_snps
    FROM $results




<div>Number of rows: 1</div><div>Query job ID  : job_DLEZFiOvJpAr7CfdRadEada1r30</div>
<table><tr><th>num_significant_snps</th></tr><tr><td>230064</td></tr></table>



We now have a dataset that is sufficiently small to fit into memory on our
instance, so let's pull the top 1000 SNP positions locally. Since we only need a
subset of the columns, we can project our data first to remove unneeded columns.


    %%bq_sql sig_snps_dataset
    SELECT * FROM (
        SELECT
          reference_name
          , start
          , reference_bases
          , alternate_bases
          , chi_squared_score
        FROM $results
        LIMIT 1000
    )
    ORDER BY start asc


    sig_snps = sig_snps_dataset.results().to_dataframe()
    sig_snps[:10]




<div style="max-height:1000px;max-width:1500px;overflow:auto;">
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>alternate_bases</th>
      <th>chi_squared_score</th>
      <th>reference_bases</th>
      <th>reference_name</th>
      <th>start</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td> G</td>
      <td> 747.021</td>
      <td> A</td>
      <td> 12</td>
      <td> 1595431</td>
    </tr>
    <tr>
      <th>1</th>
      <td> T</td>
      <td> 831.900</td>
      <td> A</td>
      <td> 12</td>
      <td> 1600977</td>
    </tr>
    <tr>
      <th>2</th>
      <td> T</td>
      <td> 714.980</td>
      <td> C</td>
      <td> 12</td>
      <td> 2129565</td>
    </tr>
    <tr>
      <th>3</th>
      <td> C</td>
      <td> 757.336</td>
      <td> T</td>
      <td> 12</td>
      <td> 2132494</td>
    </tr>
    <tr>
      <th>4</th>
      <td> G</td>
      <td> 826.659</td>
      <td> A</td>
      <td> 12</td>
      <td> 2134197</td>
    </tr>
    <tr>
      <th>5</th>
      <td> A</td>
      <td> 709.524</td>
      <td> C</td>
      <td> 12</td>
      <td> 2138687</td>
    </tr>
    <tr>
      <th>6</th>
      <td> C</td>
      <td> 772.951</td>
      <td> A</td>
      <td> 12</td>
      <td> 2139264</td>
    </tr>
    <tr>
      <th>7</th>
      <td> C</td>
      <td> 708.078</td>
      <td> T</td>
      <td> 12</td>
      <td> 2140624</td>
    </tr>
    <tr>
      <th>8</th>
      <td> T</td>
      <td> 793.613</td>
      <td> C</td>
      <td> 12</td>
      <td> 2141908</td>
    </tr>
    <tr>
      <th>9</th>
      <td> G</td>
      <td> 772.551</td>
      <td> T</td>
      <td> 12</td>
      <td> 6505528</td>
    </tr>
  </tbody>
</table>
</div>



Let's visualize the distribution of significant SNPs along the length of the
chromosome. The y-value of the charts indicates the Chi-squared score: larger
values are more significant.


    import matplotlib.pyplot as plt
    import seaborn as sns
    
    #g = sns.distplot(sig_snps['start'], rug=False, hist=False, kde_kws=dict(bw=0.1))
    fig, ax = plt.subplots()
    ax.scatter(sig_snps['start'], sig_snps['chi_squared_score'], alpha=0.3, c='red')
    ax.set_ylabel('Chi-squared score')
    ax.set_xlabel('SNP position (bp)')




    <matplotlib.text.Text at 0x7f260bf4a810>




![png](output_48_1.png)


Let's zoom in on one region that contains a large number of very significant
SNPs:


    fig, ax = plt.subplots()
    ax.scatter(sig_snps['start'], sig_snps['chi_squared_score'], alpha=0.5, c='red')
    ax.set_xlim([3.3e7, 3.5e7])
    ax.set_ylabel('Chi-squared score')
    ax.set_xlabel('SNP position (bp)')




    <matplotlib.text.Text at 0x7f260bed0dd0>




![png](output_50_1.png)


## Further exploration of statistically significant variant positions

We can take our analysis further by mapping selected variant positions back to
the chromosome and visualizing call sets and reads. Let's retrieve the top SNP
identified when ranked by the Chi-squared score:


    %%bq_sql top_snp
    SELECT start
    FROM $sig_snps_dataset
    ORDER BY chi_squared_score desc
    LIMIT 1


    top_snp.results()




<div>Number of rows: 1</div><div>Query job ID  : job_E5lc6akBdK-DYYAPL_eDsl06O0g</div>
<table><tr><th>start</th></tr><tr><td>110571373</td></tr></table>



Grab an arbitrary set of 10 callset IDs for rendering in the genome browser.


    %%bq_sql callset_ids
    SELECT * FROM (
      SELECT call.call_set_id AS callset_id
      FROM $variants_table
      GROUP BY callset_id)
    LIMIT 10


    callsets_df = callset_ids.results().to_dataframe()
    callsets = list(callsets_df['callset_id'])


    # TODO: move this to some library function
    from IPython.display import HTML
    def gabrowse(dataset, reference_name, start_position, callset_ids):
        callsets_query_params = ''.join('&callsetId=%s&cBackend=GOOGLE' % callset_id for callset_id in callset_ids)
        url = ('https://gabrowse.appspot.com/#=&readsetId=CMvnhpKTFhDd-bD9uZqjgw0&backend=GOOGLE&location=12%3A'
             + str(start_position)
             + callsets_query_params)
        return HTML('<iframe src="%s" width=1024 height=800></iframe>' % url)

Now we can render the call sets and reads for the selected SNP position by
embedding the GABrowse application directly in our notebook.


    gabrowse('1000genomes', '12', 110571373, callsets)




<iframe src="https://gabrowse.appspot.com/#=&readsetId=CMvnhpKTFhDd-bD9uZqjgw0&backend=GOOGLE&location=12%3A110571373&callsetId=10473108253681171589-0&cBackend=GOOGLE&callsetId=10473108253681171589-1&cBackend=GOOGLE&callsetId=10473108253681171589-2&cBackend=GOOGLE&callsetId=10473108253681171589-3&cBackend=GOOGLE&callsetId=10473108253681171589-4&cBackend=GOOGLE&callsetId=10473108253681171589-5&cBackend=GOOGLE&callsetId=10473108253681171589-6&cBackend=GOOGLE&callsetId=10473108253681171589-7&cBackend=GOOGLE&callsetId=10473108253681171589-8&cBackend=GOOGLE&callsetId=10473108253681171589-9&cBackend=GOOGLE" width=1024 height=800></iframe>



# Summary

This notebook illustrated how to conduct a GWAS experiment using variant data
stored within the Google Genomics BigQuery tables, retrieve a local copy of the
top results and visualize the data with Python libraries.


    
