# Example of steps for filtering pbsv output and calling genotypes

function_file <- 'pb_SV_analysis_functions.r'
source(function_file)

# directory containing the output from pbsv
data_dir <- '/DIRCTORY/WITH/PBSV/OUTPUTFILES/'

# name of VCF file that is output from pbsv
combo_vcf_short <- 'NAME_OF_COMBINED_FILE.vcf'
combo_vcf_file <- paste(data_dir, combo_vcf_short, sep = '')

# input vcf file and process into format for following steps
combo1_vcf <- make_combo_indel_df(combo_vcf_file)

# Filter pbsv output vcf
#   see below for overview of function - documentation is also in function_file
combo1_precall_df <- precalling_SV_filtering(sv_geno_df = combo1_vcf,
  mer_length = 8, per_mn_cutoff = 0.7, per_bn_pure_cutoff = 0.5,
  per_bn_multi_cutoff = 0.6, dist_cut = 1000, sd_cut = 3, use_sd = T)

# Call genotypes using filtered output
#   see below for overview of function - documentation is also in function_file
combo1_genotypes_1 <- generate_filtered_genotypes(combo_df = combo1_precall_df,
  min_sv_coverage = 10, het_ratio_cut = 0.25, min_minor_allele_count = 2,
  max_hom_ratio = 0.05, est_error = 0.032, est_penetrance = 0.35,
  good_score = 0.9, great_score = 0.98, same_geno_bonus = 0.67,
  ss_min_cov = 20)

# use filtered genotypes for remaining analysis

quite(save = 'no')

# Explanation of 'precalling_SV_filtering'
  # Run the pre-genotype calling filtering steps
  #  These include: 1) removed SVs with too many repeate k-mers; 2) removing
  #    SVs at duplicated postions; 3) Remove SVs that are too close together;
  #    4) remove SVs with missing default genotypes; 5) remove SV's with
  #    excess sequencing coverage
  # INPUTS #
  # sv_geno_df = data.frame from VCF that includes the columns <ALT> which
  #                contains the sequence for tha ALTernate allele; <REF> which
  #                contains the sequence for the REFerence allele; <type>
  #                which indicates whether the SV is INSertion or DELetion; and
  #                <sv_length> which contains the length of the SV.
  #                The SV sequence for 'INS' is in the 'ALT' column, sequence
  #                for 'DEL' is in the 'REF' column
  #                Also includes the columns <full_name> 
  #                which is character string of Chromosome and positon pasted
  #                together
  # mer_length = the total length of the nucleotide repeat sequences that
  #                are interested in
  # per_mn_cutoff = the percentage cutoff for a mononucleotide repeat. If the 
  #                  % of an SV's sequence is the repeat sequence, then it will
  #                  be removed
  # per_bn_pure_cutoff = the percentage cutoff for a single binucleotide 
  #                     repeat. If the % of an SV's sequence is 
  #                     a single type of  binucleotide repeat sequence, 
  #                     then it will be removed
  # per_bn_multi_cutoff = the percentage cutoff for two type of binucleotide
  #                      repeat. If two binucleotide repeate sequences make up
  #                      a % of an SV sequence above this cutoff, then
  #                      the SV will be removed
  # dist_cut = distance to look within for other SV's that are nearby; shortest
  #              SV withing this distance is removed
  # use_sd = use a standard deviation-based cutoff for identifying SVs with
  #            excessive sequencing coverage
  # sd_cut = if <use_sd == T>, then the number of SD's away from the mean used
  #            to identify SV's with excess coverage
  # num_coverage_cut = if <use_sd == F>, then then sequencing coverage used as
  #                      the cutoff for id'ing SV's with excess coverage
  # OUTPUT #
  # data.frame in same format at sv_geno_df but with SV filtered out from
  #  encoded filtering steps
  ########

# Explanation of 'generate_filtered_genotypes'
  # Generate genotypes straight from SV dataframe and filter based on
  #   probability-based genotype scores and minimum coverage for at least one
  #   sample
  #   Steps: 1) Generate info_list; 2) Generate preliminary genotypes;
  #     3) Assign probability-based genotype scores; 4) Adjust genotype score
  #     by genotypes in other samples; 5) Filter based on minimum coverage in
  #     at least one sample
  # INPUTS #
  # combo_df = data.frame including info about indels from the combined VCF;
  #              should be first generated using make_combo_indel_df(), 
  #              then filtered using the remove_missing_combo_SVs() or 
  #              the precalling_SV_filtering() function
  # min_sv_coverage = the required minimum coverage/total read depth for an SV;
  #                     if coverage is below <min_sv_coverage>, then call
  #                     genotype as NA
  # het_ratio_cut = the ratio of A/R or R/A above which the genotype is called
  #                   as heterozygous (or "1")
  #                   ex: if <het_ratio_cut> = 0.15, any SV's with a R/A ratio
  #                   between 0.15 and 0.85 are called as heterozygous
  # min_minor_allele_count = the required minimum number of reads for the
  #                            allele with fewer reads in order for a genotype
  #                            to be called heterozygous. If read ratio is
  #                            above <het_ratio_cut> but below 
  #                            <min_minor_allele_count>, then genotype is NA
  # max_hom_ratio = the maximum read ratio for a homozygous genotype; if
  #                   the read ratio is between <max_hom_ratio> and 
  #                   <het_ratio_cut>, then genotype is NA 
  # OUTPUTS #
  # matrix containing filtered, numerical genotypes for each sample;
  #   columns = samples, rows = SVs  
  ############

