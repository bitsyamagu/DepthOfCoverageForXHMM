# DepthOfCoverageForXHMM

## Summary

This program is specifically designed to generate input files for XHMM, offering a fast and efficient solution for calculating coverage from genetic sequencing data. 

This efficient tool is designed to rapidly generate XHMM-compatible input files with minimal computational load, distinguishing itself from GATK by its lightweight and fast performance. It's specifically engineered to operate with lower resource consumption, allowing multiple instances to run concurrently without significantly impacting system performance. This feature makes it exceptionally suitable for large-scale genetic studies, where it can process vast sample cohorts simultaneously.

## Requiirement

- htslib
  - We tested with htslib-1.19.1
- C++
  - We tested with g++ (GCC) 8.5.0 20210514

## Usage

- Input BAM file
  - BAM files should contain only one sample and SM tag designating its name
- Target region file(SureSelect, Twist etc)
  - BED format or interval_list format

```
mtdoc --bam <bamFile> --bed <bedFile> --out <outputFile> [--threads <threads>]
```
  --bed supports both BED and target region format('.interva_list')

  --threads default: 8

