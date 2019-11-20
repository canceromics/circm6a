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
* Gradle 3.0+
* A working terminal
* At least some knowledge

## Installation
```
git clone https://github.com/XXXX
cd Circm6A
```

## QuickStart
* input 

```
java -Xmx16g -jar Circm6A.jar -ip IP.bam -input Input.bam -r Gencode.gtf -g hg38.fa -o ${out_dir}/${file_prefix}
```

## Usage

* The standard way to run GATK4 tools is via the **`gatk`** wrapper script located in the root directory of a clone of this repository.
    * Requires Python 2.6 or greater (this includes Python 3.x)
    * You need to have built the GATK as described in the [Building GATK4](#building) section above before running this script.
    * There are several ways `gatk` can be run:
        * Directly from the root of your git clone after building
        * By extracting the zip archive produced by `./gradlew bundle` to a directory, and running `gatk` from there
        * Manually putting the `gatk` script within the same directory as fully-packaged GATK jars produced by `./gradlew localJar` and/or `./gradlew sparkJar`
        * Defining the environment variables `GATK_LOCAL_JAR` and `GATK_SPARK_JAR`, and setting them to the paths to the GATK jars produced by `./gradlew localJar` and/or `./gradlew sparkJar` 
    * `gatk` can run non-Spark tools as well as Spark tools, and can run Spark tools locally, on a Spark cluster, or on Google Cloud Dataproc.
    * ***Note:*** running with `java -jar` directly and bypassing `gatk` causes several important system properties to not get set, including htsjdk compression level!
    
* For help on using `gatk` itself, run **`./gatk --help`**

* To print a list of available tools, run **`./gatk --list`**.
    * Spark-based tools will have a name ending in `Spark` (eg., `BaseRecalibratorSpark`). Most other tools are non-Spark-based.

* To print help for a particular tool, run **`./gatk ToolName --help`**.

* To run a non-Spark tool, or to run a Spark tool locally, the syntax is: **`./gatk ToolName toolArguments`**.

* Tool arguments that allow multiple values, such as -I, can be supplied on the command line using a file with the extension ".args". Each line of the file should contain a
  single value for the argument.

### Example

  ```
java -Xmx16g -jar Circm6A.jar -ip IP.bam -input Input.bam -r Gencode.gtf -g hg38.fa -o ${out_dir}/${file_prefix}
  ```

####  Passing JVM options to gatk

* To pass JVM arguments to GATK, run `gatk` with the `--java-options` argument: 

    ```
    ./gatk --java-options "-Xmx4G -XX:+PrintGCDetails" <rest of command>
    ```
### Output
* output_dir/quant/quant.txt

| Field       | Description                           |
| :---------- | :------------------------------------ |
| chrom       | Chromosome                            |
| start       | Start of circular RNA                 |
| end         | End of circular RNA                   |
| name        | Circular RNA/Junction reads           |
| score       | Flag of fusion junction realignment   |
| strand      | + or - for strand                     |
| thickStart  | No meaning                            |
| thickEnd    | No meaning                            |
| itemRgb     | 0,0,0                                 |
| exonCount   | Number of exons                       |
| exonSizes   | Exon sizes                            |
| exonOffsets | Exon offsets                          |
| readNumber  | Number of junction reads              |
| circType    | Type of circular RNA                  |
| geneName    | Name of gene                          |
| isoformName | Name of isoform                       |
| index       | Index of exon or intron               |
| flankIntron | Left intron/Right intron              |



## Citation

XXXXXXX


## License
