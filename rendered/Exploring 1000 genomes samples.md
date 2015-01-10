
# Exploring the 1000 Genomes dataset

This notebook demonstrates exploring the public 1000 Genomes dataset stored in
BigQuery.

### In this notebook you will
* Explore the 1000 Genomes dataset available in BigQuery
* Use the `%%bq_sql` statement to write and execute SQL statements within the
notebook
* Extract data from BigQuery and create a local dataset that can be manipulated
in Python
* Refine and pivot your local dataset via the Pandas Python library
* Visualize different aspects of your dataset via the Matplotlib and Seaborn
Python libraries

### Attribution
This notebook is based on the [Google
Genomics](https://cloud.google.com/genomics/) BigQuery examples here:
* [1000 Genomes](https://github.com/googlegenomics/bigquery-
examples/tree/master/1000genomes/data-stories/exploring-the-phenotypic-data)

----

### Installing ad-hoc dependencies


    !pip install seaborn

## Working with data in BigQuery

First let's take a look at the [1000 Genomes](http://www.1000genomes.org/)
sample dataset table. We can use the `gcp.bigquery` Python package to fetch the
table schema.


    import gcp.bigquery as bq
    samples_table = bq.table('genomics-public-data:1000_genomes.sample_info')
    samples_table.schema()




<table><tr><th>name</th><th>data_type</th><th>mode</th><th>description</th></tr><tr><td>Sample</td><td>STRING</td><td>NULLABLE</td><td>Sample ID</td></tr><tr><td>Family_ID</td><td>STRING</td><td>NULLABLE</td><td>Family ID</td></tr><tr><td>Population</td><td>STRING</td><td>NULLABLE</td><td>3 letter population code</td></tr><tr><td>Population_Description</td><td>STRING</td><td>NULLABLE</td><td>Description of Population</td></tr><tr><td>Gender</td><td>STRING</td><td>NULLABLE</td><td>Gender</td></tr><tr><td>Relationship</td><td>STRING</td><td>NULLABLE</td><td>Relationship to other members of the family</td></tr><tr><td>Unexpected_Parent_Child</td><td>STRING</td><td>NULLABLE</td><td>sample id for unexpected parent child relationships</td></tr><tr><td>Non_Paternity</td><td>STRING</td><td>NULLABLE</td><td>sample ids for annotated non paternal relationships</td></tr><tr><td>Siblings</td><td>STRING</td><td>NULLABLE</td><td>sample ids for any siblings</td></tr><tr><td>Grandparents</td><td>STRING</td><td>NULLABLE</td><td>sample ids for any grand parents</td></tr><tr><td>Avuncular</td><td>STRING</td><td>NULLABLE</td><td>sample ids for any avuncular relationships</td></tr><tr><td>Half_Siblings</td><td>STRING</td><td>NULLABLE</td><td>sample ids for any half siblings</td></tr><tr><td>Unknown_Second_Order</td><td>STRING</td><td>NULLABLE</td><td>sample ids for any unknown second order relations</td></tr><tr><td>Third_Order</td><td>STRING</td><td>NULLABLE</td><td>sample ids for any third order cryptic relations. As mentioned above, this analysis was not as widely run as the other relatedness analyses and as such there may still be unannotated third order relations in the set</td></tr><tr><td>In_Low_Coverage_Pilot</td><td>BOOLEAN</td><td>NULLABLE</td><td>The sample is in the low coverage pilot experiment</td></tr><tr><td>LC_Pilot_Platforms</td><td>STRING</td><td>NULLABLE</td><td>low coverage pilot sequencing platforms </td></tr><tr><td>LC_Pilot_Centers</td><td>STRING</td><td>NULLABLE</td><td>low coverage pilot sequencing centers</td></tr><tr><td>In_High_Coverage_Pilot</td><td>BOOLEAN</td><td>NULLABLE</td><td>The sample is in the high coverage pilot</td></tr><tr><td>HC_Pilot_Platforms</td><td>STRING</td><td>NULLABLE</td><td>high coverage sequencing platforms</td></tr><tr><td>HC_Pilot_Centers</td><td>STRING</td><td>NULLABLE</td><td>high coverage sequencing centers</td></tr><tr><td>In_Exon_Targetted_Pilot</td><td>BOOLEAN</td><td>NULLABLE</td><td>The Sample is in the exon targetted pilot experiment</td></tr><tr><td>ET_Pilot_Platforms</td><td>STRING</td><td>NULLABLE</td><td>exon targetted sequencing platforms,</td></tr><tr><td>ET_Pilot_Centers</td><td>STRING</td><td>NULLABLE</td><td>exon targetted sequencing centers,</td></tr><tr><td>Has_Sequence_in_Phase1</td><td>BOOLEAN</td><td>NULLABLE</td><td>Has sequence low coverage sequence in the 20101123.sequence.index file or exome sequence in the 20110522 sequence index file</td></tr><tr><td>Phase1_LC_Platform</td><td>STRING</td><td>NULLABLE</td><td>phase1 low coverage sequencing platforms</td></tr><tr><td>Phase1_LC_Centers</td><td>STRING</td><td>NULLABLE</td><td>phase1 low coverage sequencing centers</td></tr><tr><td>Phase1_E_Platform</td><td>STRING</td><td>NULLABLE</td><td>phase1 exome sequencing platforms</td></tr><tr><td>Phase1_E_Centers</td><td>STRING</td><td>NULLABLE</td><td>phase1 exome sequencing centers</td></tr><tr><td>In_Phase1_Integrated_Variant_Set</td><td>BOOLEAN</td><td>NULLABLE</td><td>The sample is genotyped in the phase1 integrated call set on autosomes and chrX</td></tr><tr><td>Has_Phase1_chrY_SNPS</td><td>BOOLEAN</td><td>NULLABLE</td><td>The sample is genotyped in the chrY phase1 snp set</td></tr><tr><td>Has_phase1_chrY_Deletions</td><td>BOOLEAN</td><td>NULLABLE</td><td>The sample is genotyepd in the chrY phase1 deletions</td></tr><tr><td>Has_phase1_chrMT_SNPs</td><td>BOOLEAN</td><td>NULLABLE</td><td>The sample is genotyped in the phase1 chrMT snps</td></tr><tr><td>Main_project_LC_Centers</td><td>STRING</td><td>NULLABLE</td><td>low coverage sequencing centers for final sequencing round</td></tr><tr><td>Main_project_LC_platform</td><td>STRING</td><td>NULLABLE</td><td>low coverage sequencing platform for final sequencing round</td></tr><tr><td>Total_LC_Sequence</td><td>FLOAT</td><td>NULLABLE</td><td>The total amount of low coverage sequence available</td></tr><tr><td>LC_Non_Duplicated_Aligned_Coverage</td><td>FLOAT</td><td>NULLABLE</td><td>The non duplicated aligned coverage for the low coverage sequence data.  This was calculated using the ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/alignment_indices/20130502.low_coverage.alignment.index.bas.gz file, the (mapped bases - duplicated bases) was summed for each sample and divided by 2.75GB and rounded to 2dp</td></tr><tr><td>Main_Project_E_Centers</td><td>STRING</td><td>NULLABLE</td><td>Exome sequencing centers for the final sequencing round</td></tr><tr><td>Main_Project_E_Platform</td><td>STRING</td><td>NULLABLE</td><td>Exome sequencing platform for the final sequencing round</td></tr><tr><td>Total_Exome_Sequence</td><td>FLOAT</td><td>NULLABLE</td><td>The total amount of exome sequence available</td></tr><tr><td>X_Targets_Covered_to_20x_or_greater</td><td>FLOAT</td><td>NULLABLE</td><td>The percentage of targets covered to 20x or greater as calculated by the picard function CalculateHsMetrics using these targets ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/exome_pull_down_targets_phases1_and_2/20120518.consensus.annotation.bed</td></tr><tr><td>VerifyBam_E_Omni_Free</td><td>FLOAT</td><td>NULLABLE</td><td>Value from UMich's VerifyBamID BAM QC program http://genome.sph.umich.edu/wiki/VerifyBamID.  The Free measures use a statistical model based on the haplotypes discovered by the chip. The Chip measure considers the genotypes available for that individual from that chip. We use greater than 3% as a cut off for our our low coverage samples and greater than 3.5% for our exome samples.</td></tr><tr><td>VerifyBam_E_Affy_Free</td><td>FLOAT</td><td>NULLABLE</td><td>Value from UMich's VerifyBamID BAM QC program http://genome.sph.umich.edu/wiki/VerifyBamID.  The Free measures use a statistical model based on the haplotypes discovered by the chip. The Chip measure considers the genotypes available for that individual from that chip. We use greater than 3% as a cut off for our our low coverage samples and greater than 3.5% for our exome samples.</td></tr><tr><td>VerifyBam_E_Omni_Chip</td><td>FLOAT</td><td>NULLABLE</td><td>Value from UMich's VerifyBamID BAM QC program http://genome.sph.umich.edu/wiki/VerifyBamID.  The Free measures use a statistical model based on the haplotypes discovered by the chip. The Chip measure considers the genotypes available for that individual from that chip. We use greater than 3% as a cut off for our our low coverage samples and greater than 3.5% for our exome samples.</td></tr><tr><td>VerifyBam_E_Affy_Chip</td><td>FLOAT</td><td>NULLABLE</td><td>Value from UMich's VerifyBamID BAM QC program http://genome.sph.umich.edu/wiki/VerifyBamID.  The Free measures use a statistical model based on the haplotypes discovered by the chip. The Chip measure considers the genotypes available for that individual from that chip. We use greater than 3% as a cut off for our our low coverage samples and greater than 3.5% for our exome samples.</td></tr><tr><td>VerifyBam_LC_Omni_Free</td><td>FLOAT</td><td>NULLABLE</td><td>Value from UMich's VerifyBamID BAM QC program http://genome.sph.umich.edu/wiki/VerifyBamID.  The Free measures use a statistical model based on the haplotypes discovered by the chip. The Chip measure considers the genotypes available for that individual from that chip. We use greater than 3% as a cut off for our our low coverage samples and greater than 3.5% for our exome samples.</td></tr><tr><td>VerifyBam_LC_Affy_Free</td><td>FLOAT</td><td>NULLABLE</td><td>Value from UMich's VerifyBamID BAM QC program http://genome.sph.umich.edu/wiki/VerifyBamID.  The Free measures use a statistical model based on the haplotypes discovered by the chip. The Chip measure considers the genotypes available for that individual from that chip. We use greater than 3% as a cut off for our our low coverage samples and greater than 3.5% for our exome samples.</td></tr><tr><td>VerifyBam_LC_Omni_Chip</td><td>FLOAT</td><td>NULLABLE</td><td>Value from UMich's VerifyBamID BAM QC program http://genome.sph.umich.edu/wiki/VerifyBamID.  The Free measures use a statistical model based on the haplotypes discovered by the chip. The Chip measure considers the genotypes available for that individual from that chip. We use greater than 3% as a cut off for our our low coverage samples and greater than 3.5% for our exome samples.</td></tr><tr><td>VerifyBam_LC_Affy_Chip</td><td>FLOAT</td><td>NULLABLE</td><td>Value from UMich's VerifyBamID BAM QC program http://genome.sph.umich.edu/wiki/VerifyBamID.  The Free measures use a statistical model based on the haplotypes discovered by the chip. The Chip measure considers the genotypes available for that individual from that chip. We use greater than 3% as a cut off for our our low coverage samples and greater than 3.5% for our exome samples.</td></tr><tr><td>LC_Indel_Ratio</td><td>FLOAT</td><td>NULLABLE</td><td>Both Indel ratios are the ratio of insertions to deletions found in that sample using a quick test (based on samtools). If the ratio is higher than 5 the sample is withdrawn.</td></tr><tr><td>E_Indel_Ratio</td><td>FLOAT</td><td>NULLABLE</td><td>Both Indel ratios are the ratio of insertions to deletions found in that sample using a quick test (based on samtools). If the ratio is higher than 5 the sample is withdrawn.</td></tr><tr><td>LC_Passed_QC</td><td>BOOLEAN</td><td>NULLABLE</td><td>These are binary flags showing if the sample passed QC, All samples which have passed QC have bam files. Only samples which have both exome and low coverage data are found under ftp/data and listed in the standard alignment index. The small number of other samples are found in ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/phase3_EX_or_LC_only_alignment/</td></tr><tr><td>E_Passed_QC</td><td>BOOLEAN</td><td>NULLABLE</td><td>These are binary flags showing if the sample passed QC, All samples which have passed QC have bam files. Only samples which have both exome and low coverage data are found under ftp/data and listed in the standard alignment index. The small number of other samples are found in ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/phase3_EX_or_LC_only_alignment/</td></tr><tr><td>In_Final_Phase_Variant_Calling</td><td>BOOLEAN</td><td>NULLABLE</td><td>Any sample which has both LC and E QC passed bams is in the final analysis set</td></tr><tr><td>Has_Omni_Genotypes</td><td>BOOLEAN</td><td>NULLABLE</td><td>Omni Genotypes in ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20120131_omni_genotypes_and_intensities/Omni25_genotypes_2141_samples.b37.vcf.gz   </td></tr><tr><td>Has_Axiom_Genotypes</td><td>BOOLEAN</td><td>NULLABLE</td><td>Axiom Genotypes in ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20110210_Affymetrix_Axiom/Affymetrix_Axiom_DB_2010_v4_b37.vcf.gz         </td></tr><tr><td>Has_Affy_6_0_Genotypes</td><td>BOOLEAN</td><td>NULLABLE</td><td>Affy 6.0 Genotypes in  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20121128_corriel_p3_sample_genotypes/</td></tr><tr><td>Has_Exome_LOF_Genotypes</td><td>BOOLEAN</td><td>NULLABLE</td><td>ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20121009_broad_exome_chip/ALL.wgs.broad_exome_lof_indel_v2.20121009.snps_and_indels.snpchip.genotypes.vcf.gz</td></tr><tr><td>EBV_Coverage</td><td>FLOAT</td><td>NULLABLE</td><td>This was calculated by looking at the alignment of the data to NC_007605 in the low coverage bam files and using that to calculate coverage   </td></tr><tr><td>DNA_Source_from_Coriell</td><td>STRING</td><td>NULLABLE</td><td>This was the annotated DNA Source from Coriell     </td></tr><tr><td>Has_Sequence_from_Blood_in_Index</td><td>BOOLEAN</td><td>NULLABLE</td><td>In the later stages of the project some populations has multiple study ids, one to indicate sequencing from blood. This data for each sample has not been treated independently in the alignment process but when there is both LCL and Blood sourced data they are both together in single bams</td></tr><tr><td>Super_Population</td><td>STRING</td><td>NULLABLE</td><td>From ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/20131219.superpopulations.tsv</td></tr><tr><td>Super_Population_Description</td><td>STRING</td><td>NULLABLE</td><td>From ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/20131219.superpopulations.tsv</td></tr></table>



We can see in the table schema that a number of different annotations exist for
each genomic sample. To get a feel for the dataset, let's see how the samples
are distributed across populations and super populations.

The `%%bq_sql` statement allows us to write SQL within our notebook and execute
it within BigQuery. In the query below, we'll refer to our previously defined
`samples_table` Python variable as `$samples_table`.


    %%bq_sql pops
    SELECT
      population,
      population_description,
      super_population,
      super_population_description,
      COUNT(population) AS population_count,
    FROM
      $samples_table
    GROUP BY
      population,
      population_description,
      super_population,
      super_population_description

We've just defined the query so far, nothing has been executed yet. We can refer
to this defined query by the name after the `%%bq_sql` declaration ("pops" in
this case). We can fire off the query and print the results like so:


    pops.results()




<div>Number of rows: 26</div><div>Query job ID  : job_tslQM-tp9voSmTgyOv0Pu6L3-hc</div>
<table><tr><th>population_description</th><th>super_population</th><th>population_count</th><th>super_population_description</th><th>population</th></tr><tr><td>British in England and Scotland</td><td>EUR</td><td>107</td><td>European</td><td>GBR</td></tr><tr><td>Finnish in Finland</td><td>EUR</td><td>105</td><td>European</td><td>FIN</td></tr><tr><td>Southern Han Chinese, China</td><td>EAS</td><td>171</td><td>East Asian</td><td>CHS</td></tr><tr><td>Puerto Rican in Puerto Rico</td><td>AMR</td><td>150</td><td>American</td><td>PUR</td></tr><tr><td>Chinese Dai in Xishuangbanna, China</td><td>EAS</td><td>109</td><td>East Asian</td><td>CDX</td></tr><tr><td>Colombian in Medellin, Colombia</td><td>AMR</td><td>148</td><td>American</td><td>CLM</td></tr><tr><td>Iberian populations in Spain</td><td>EUR</td><td>162</td><td>European</td><td>IBS</td></tr><tr><td>Peruvian in Lima, Peru</td><td>AMR</td><td>130</td><td>American</td><td>PEL</td></tr><tr><td>Punjabi in Lahore,Pakistan</td><td>SAS</td><td>158</td><td>South Asian</td><td>PJL</td></tr><tr><td>Kinh in Ho Chi Minh City, Vietnam</td><td>EAS</td><td>124</td><td>East Asian</td><td>KHV</td></tr><tr><td>African Caribbean in Barbados</td><td>AFR</td><td>123</td><td>African</td><td>ACB</td></tr><tr><td>Gambian in Western Division, The Gambia</td><td>AFR</td><td>180</td><td>African</td><td>GWD</td></tr><tr><td>Esan in Nigeria</td><td>AFR</td><td>173</td><td>African</td><td>ESN</td></tr><tr><td>Bengali in Bangladesh</td><td>SAS</td><td>144</td><td>South Asian</td><td>BEB</td></tr><tr><td>Mende in Sierra Leone</td><td>AFR</td><td>128</td><td>African</td><td>MSL</td></tr><tr><td>Sri Lankan Tamil in the UK</td><td>SAS</td><td>128</td><td>South Asian</td><td>STU</td></tr><tr><td>Indian Telugu in the UK</td><td>SAS</td><td>118</td><td>South Asian</td><td>ITU</td></tr><tr><td>Utah residents with Northern and Western European ancestry</td><td>EUR</td><td>183</td><td>European</td><td>CEU</td></tr><tr><td>Yoruba in Ibadan, Nigeria</td><td>AFR</td><td>186</td><td>African</td><td>YRI</td></tr><tr><td>Han Chinese in Bejing, China</td><td>EAS</td><td>108</td><td>East Asian</td><td>CHB</td></tr><tr><td>Japanese in Tokyo, Japan</td><td>EAS</td><td>105</td><td>East Asian</td><td>JPT</td></tr><tr><td>Luhya in Webuye, Kenya</td><td>AFR</td><td>116</td><td>African</td><td>LWK</td></tr><tr><td>African Ancestry in Southwest US</td><td>AFR</td><td>112</td><td>African</td><td>ASW</td></tr><tr><td>Mexican Ancestry in Los Angeles, California</td><td>AMR</td><td>107</td><td>American</td><td>MXL</td></tr><tr><td>Toscani in Italy</td><td>EUR</td><td>112</td><td>European</td><td>TSI</td></tr><tr><td>Gujarati Indian in Houston,TX</td><td>SAS</td><td>113</td><td>South Asian</td><td>GIH</td></tr></table>



In order to further analyze our query results locally within the IPython
notebook, let's convert the result set into a [Pandas
dataframe](http://pandas.pydata.org/pandas-
docs/dev/generated/pandas.DataFrame.html):


    pops_df = pops.results().to_dataframe()
    pops_df[:5]




<div style="max-height:1000px;max-width:1500px;overflow:auto;">
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>population</th>
      <th>population_count</th>
      <th>population_description</th>
      <th>super_population</th>
      <th>super_population_description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td> GBR</td>
      <td> 107</td>
      <td>     British in England and Scotland</td>
      <td> EUR</td>
      <td>   European</td>
    </tr>
    <tr>
      <th>1</th>
      <td> FIN</td>
      <td> 105</td>
      <td>                  Finnish in Finland</td>
      <td> EUR</td>
      <td>   European</td>
    </tr>
    <tr>
      <th>2</th>
      <td> CHS</td>
      <td> 171</td>
      <td>         Southern Han Chinese, China</td>
      <td> EAS</td>
      <td> East Asian</td>
    </tr>
    <tr>
      <th>3</th>
      <td> PUR</td>
      <td> 150</td>
      <td>         Puerto Rican in Puerto Rico</td>
      <td> AMR</td>
      <td>   American</td>
    </tr>
    <tr>
      <th>4</th>
      <td> CDX</td>
      <td> 109</td>
      <td> Chinese Dai in Xishuangbanna, China</td>
      <td> EAS</td>
      <td> East Asian</td>
    </tr>
  </tbody>
</table>
</div>



Pandas dataframes have a [`dataframe.plot()`](http://pandas.pydata.org/pandas-do
cs/dev/generated/pandas.DataFrame.plot.html?highlight=dataframe.plot#pandas.Data
Frame.plot) method that allows us to quickly render a number of standard
visualizations. Let's draw a bar chart of the per-population sample counts:


    import seaborn as sns
    pops_sum_df.plot(kind='bar', y='population_count')




    <matplotlib.axes._subplots.AxesSubplot at 0x7fdc743717d0>




![png](output_14_1.png)


Now that we have a dataframe object, we can compute arbitrary rollups and
aggregations locally. For example, let's aggregate the count of samples from the
population level to the super population level. We'll do this via the
[`dataframe.groupby`](http://pandas.pydata.org/pandas-docs/dev/generated/pandas.
DataFrame.groupby.html?highlight=groupby#pandas.DataFrame.groupby) operation,
which is generally of the form:

```dataframe.groupby("column name").aggregation_func()```


    superpop_df = pops_df.groupby('super_population').sum()
    superpop_df.plot(kind='bar')




    <matplotlib.axes._subplots.AxesSubplot at 0x7fdc70e9c050>




![png](output_16_1.png)


We see that the distribution of genome samples across population and super
population is relatively uniform, the exception being the AFR (African) super
population.

## Genome sample metrics

Let's explore a few quantitative attributes of the genomic samples. We'll keep
our super population sample classification, but also add some metrics around the
extent to which given samples were sequenced, for exomic and low coverage
regions. We'll further annotate our samples with a tag for the laboratory/center
that produced performed the sequencing.


    %%bq_sql metrics
    select 
        super_population
        , total_lc_sequence -- lc=low coverage
        , total_exome_sequence
        , Main_Project_E_Centers
    from $samples_table
    where 
        total_exome_sequence is not null
        and total_lc_sequence is not null
        and Main_Project_E_Centers is not null
        and Main_Project_E_Centers != 'BCM,BGI' -- remove this single outlier

Again, we can convert these results to a Pandas dataframe for local analysis and
visualization


    df = metrics.results().to_dataframe()
    df[:10]




<div style="max-height:1000px;max-width:1500px;overflow:auto;">
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Main_Project_E_Centers</th>
      <th>super_population</th>
      <th>total_exome_sequence</th>
      <th>total_lc_sequence</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td> WUGSC</td>
      <td> EUR</td>
      <td>  9580861000</td>
      <td> 14714705600</td>
    </tr>
    <tr>
      <th>1</th>
      <td>   BCM</td>
      <td> EUR</td>
      <td>  5856213714</td>
      <td> 30285067857</td>
    </tr>
    <tr>
      <th>2</th>
      <td>   BCM</td>
      <td> EUR</td>
      <td>  8125011458</td>
      <td> 25134554475</td>
    </tr>
    <tr>
      <th>3</th>
      <td>    BI</td>
      <td> EUR</td>
      <td> 15144863140</td>
      <td> 42906753668</td>
    </tr>
    <tr>
      <th>4</th>
      <td> WUGSC</td>
      <td> EUR</td>
      <td> 11754108800</td>
      <td> 22115533074</td>
    </tr>
    <tr>
      <th>5</th>
      <td> WUGSC</td>
      <td> EUR</td>
      <td>  9003406100</td>
      <td> 21621189887</td>
    </tr>
    <tr>
      <th>6</th>
      <td> WUGSC</td>
      <td> EUR</td>
      <td>  9145497300</td>
      <td> 16331033600</td>
    </tr>
    <tr>
      <th>7</th>
      <td>   BCM</td>
      <td> EUR</td>
      <td>  5851507215</td>
      <td> 22262720677</td>
    </tr>
    <tr>
      <th>8</th>
      <td>   BCM</td>
      <td> EUR</td>
      <td>  6640380441</td>
      <td> 24284331308</td>
    </tr>
    <tr>
      <th>9</th>
      <td> WUGSC</td>
      <td> EUR</td>
      <td>  9595593000</td>
      <td> 23778576248</td>
    </tr>
  </tbody>
</table>
</div>



To get a feel for the quantitative sample attributes, let's see how their values
are distributed among the dataset samples by plotting histograms. Note that a
histogram of any dataframe attribute can be easily rendered with the pattern
`dataframe.attribute_name.hist()`


    df.hist(alpha=0.5, bins=30)




    array([[<matplotlib.axes._subplots.AxesSubplot object at 0x7fdc71e434d0>,
            <matplotlib.axes._subplots.AxesSubplot object at 0x7fdc72f87410>]], dtype=object)




![png](output_24_1.png)


What does the joint distribution of these two traits look like?


    g = sns.jointplot('total_exome_sequence', 'total_lc_sequence', df)


![png](output_26_0.png)


Only a very slight positive correlation. Let's further annotate our scatter
chart by rendering each mark with a color according to its super population
assignment.


    g = sns.lmplot('total_exome_sequence', 'total_lc_sequence', hue='super_population', data=df, fit_reg=False)


![png](output_28_0.png)


Now, let's take the same plot as above, by facet our results based upon the
genomic sequencing center that produced it to look for inter-center variability
in the dataset.


    g = sns.lmplot('total_exome_sequence', 'total_lc_sequence', hue='super_population', col='Main_Project_E_Centers', col_wrap=2, data=df, fit_reg=False)


![png](output_30_0.png)


The WUSGC (the Genome Institute at Washington University) shows a small outlier
cluster that is distinct relative to the other centers. The BCM (Baylor College
of Medicine) facet appears the least variable within the exome sequencing
dimension.

Are there any super population trends here? We can facet our data a second time,
this time by super population to dig deeper.


    g = sns.lmplot('total_exome_sequence', 'total_lc_sequence', hue='super_population', col='super_population', row='Main_Project_E_Centers', data=df, fit_reg=False)


![png](output_32_0.png)


With the data now faceted by super population as well as genomic center, the
variability of the Broad institute data becomes more apparent. Furthermore, the
outlier cluster noted earlier in the WUGSC facet appears to be primarily
constituted by the EAS and AMR super populations, with no representation from
the SAS super population.

## Summary

You should have the tools necessary for interacting with genomic data stored in
BigQuery, be able to create curated local datasets and finally, visualize and
explore your result sets.


    
