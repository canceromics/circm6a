<img src="icon.png" align="right" />

# Circm6A

Circm6A is a powerful tool for detection  m6A modification of circular RNA(circRNA).



## Table of Contents
* [Requirements](#Requirements)
* [Installation](#Installation)
* [Quick Start Guide](#QuickStart)
* [Usage ](#Usage)
* [Example](#Example)
* [Output Headers](#OutputHeaders)
* [License](#License)

## Requirements

* JDK 8

## Installation

* This tool can be installed by instructions as follows:

```
git clone https://github.com/canceromics/circm6a.git
cd circm6a/huntcircRNA
javac -d ./ -classpath ./lib/htsjdk-2.10.1.jar ./src/main/*.java
jar -cvmf META-INF/MANIFEST.MF circm6a.jar *
```
The tool is generated as circm6a.jar in this directory.

## QuickStart

* Start from bam files of IP sample and Input sample

```
java -Xmx16g -jar circm6A.jar -ip IP.bam -input Input.bam -o ./example
```

Running this instruction will result in getting a file named example_circ.bed. `./example` means output_dir/file_prefix

## Usage

* More details of this tool can be found with -h parameter

```
java -Xmx16g -jar circm6A.jar -h
For usage:
	-ip	IP sam/bam file searching back junctions
	-input	INPUT sam/bam file searching back junctions, calling IP peaks with IP file)
	-trim	trimed sam/bam file calling peak instead of INPUT file
	-g	a reference genome file which is the same as the mapping one
	-r	a reference gencode file in gtf format
	-o	prefix of out files
	-rrna	enable rRNA remove (can specify a bed file)
	-circ	provide circ bed file to call peak
	-adjust	choose method for p-value adjusting [bon|bh](bon means Bonferroni, bh means Benjamini and Hochberg)
	-sup	min support number of a back junction (default is 1)
	-back	background size in peak calling (0 for whole genome, default is 0)
	-window	window size in peak calling (default is 25)
	-peak	min peak length in peak calling (default is 100)
	-cl	cutoff of distance bewteen junction points (default is 100)
	-clmax	upper bound cutoff of distance bewteen junction points (default is 200000)
	-dev	max deviation permitted in searching back junctions (default is 5)
	-pt	threshold for p-value in peak calling (default is 0.05)
	-detail	output circ detail file
	-pair	enable pair end examination of back junction
	-uniq	enable using only uniq mapping reads
	-h	show help text
```

* Tips: -g enables GT-AG filter of circRNA detecting. -r enables exon boundary filter of circRNA detecting.

## Example

  ```
java -Xmx16g -jar circm6A.jar -ip ../test_data/HeLa_eluate_rep_1.chr22.bam -input ../test_data/HeLa_input_rep_1.chr22.bam -r ../test_data/gencode_chr22.gtf -g ../test_data/hg38_chr22.fa -o ../test_data/test_Hela
  ```

## OutputHeaders

* Here are definitions of headers in output file named `(output_dir/file_prefix)_circ.bed`

| Field       | Description                           |
| :---------- | :------------------------------------ |
| chrom       | Chromosome                            |
| chromStart       | Start of circular RNA                 |
| chromEnd         | End of circular RNA                   |
| name        | Circular RNA/Junction reads           |
| score       | Flag of fusion junction realignment   |
| strand      | + or - for strand                     |
| thickStart  | No meaning                            |
| thickEnd    | No meaning                            |
| itemRgb     | 0,0,0                                 |
| blockCount   | Number of exons                       |
| blockSizes   | Exon sizes                            |
| blockStarts | Exon offsets                          |



## License
Licensed GPLv3 for open source use or contact zuoLab (zuozhx@sysucc.org.cn) for commercial use.
