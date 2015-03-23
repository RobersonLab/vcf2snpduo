# vcf2snpduo
Simple python script for converting a valid, multi-sample VCF file into a SNPduo compatible format.

## Description
The [SNPduo](http://github.com/RobersonLab/snpduoweb) web-based tool is capable of generating plots of identity-by-state as well as tabulating IBS state counts. However, it was designed to analyze SNP array data. The included script represents an easy way to convert a multi-sample VCF into a format readable by SNPduo. Data output goes to stdout, log messages go to stderr.

## Basic usage

### Run example
Filter by AD
```bash
python vcf2snpduo.py --minDepth 20 --minSampleCount 4 myInput.vcf 1>custom_duo_format.csv 2>duo.log
```

Filter by DP
```bash
python vcf2snpduo.py --minDepth 20 --minSampleCount 4 --filterDP myInput.vcf 1>custom_duo_format.csv 2>duo.log
```

No depth filter
```bash
python vcf2snpduo.py --noDepthFilter myInput.vcf 1>custom_duo_format.csv 2>duo.log
```

**Important note** The depth filtering for this tool uses the *individual* level data, i.e. the DP field in INFO for all samples at a site is not used. Not all VCF files have both allele depth (AD) and depth (DP) per sample, so choose based on your file format. If neither are present or depth filtering is not requested, use the `--noDepthFilter` option.

### Running in SNPduo
Follow instructions for the SNPduo repository if you aren't running the tool via the [Pevsner lab](http://pevsnerlab.kennedykrieger.org/SNPduo/) web-interface. Choose CSV delimiting character and Custom Data Format.

## Important notes
* The higher the density of genotype calls, the better the resolution. Not optimal for low-coverage exomes.
* When exporting genotype calls using tools such as the [Genome Analysis ToolKit](https://www.broadinstitute.org/gatk/), --excludeNonVariants should **not** be invoked.
* Performance may be improved by including only variants concordant with known polymorphic (MAF >1%) sites
