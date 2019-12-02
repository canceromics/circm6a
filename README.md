# Circm6A

Circm6A is a powerful tool for detection  m6A modification of circular RNA(circRNA).



## Table of Contents
* [Requirements](#Requirements)
* [Installation](#Installation)
* [Quick Start Guide](#QuickStart)
* [Usage ](#Usage)
* [Authors](#Authors)
* [Citation](#Citation)
* [License](#License)

## Requirements

* JDK 8

## Installation
```
git clone https://github.com/canceromics/circm6a.git
cd circm6a/huntcircRNA
javac -d ./ -classpath ./lib/htsjdk-2.10.1.jar ./src/main/*.java
jar -cvmf META-INF/MANIFEST.MF circm6a.jar *
```

## QuickStart
* Start from bam files of IP sample and Input sample

```
java -Xmx16g -jar circm6A.jar -ip IP.bam -input Input.bam -r Gencode.gtf -g hg38.fa -o ${out_dir}/${file_prefix}
```

## Usage

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

### Example

  ```
java -Xmx16g -jar circm6A.jar -ip IP.bam -input Input.bam -r Gencode.gtf -g hg38.fa -o ${out_dir}/${file_prefix} -rrna -pair  -dev 5 -sup 5
  ```

### Output
* ${out_dir}/${file_prefix}_circ.bed

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




## Citation

XXXXXXX


## License
