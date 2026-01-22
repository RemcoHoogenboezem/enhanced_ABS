# enhanced_ABS\
\
Annotate annovar or VEP file with statistics from one or more bam files using a sophisticated realignment algorithm for indel loci.\
\
Several system libraries are required to compile this project:\
On debian and ubuntu: libhts-dev, libparasail-dev, zlib1g-dev, libbz2-dev, liblzma-dev, libcurl4-gnutls-dev or libcurl4-openssl-dev, libssl-dev, libdeflate-dev, libomp-dev\
\
To compile use cmake:\
cmake -DCMAKE_BUILD_TYPE=Release -S ./\
make\
\
Run options:\
enhanced_ABS v2.0 [options] > output.txt (output always to std::out)\
\
|Short option|Long option|Type|Description|
|---|---|---|---|
|-f|--fasta-file            |text |Single fasta file(required)                                                             |
|-V|--vcf-file              |text |Single vcf file (optional)                                                              |
|-a|--annovar-file          |text |Single annovar file (required either annovar-file or vep-file)                          |
|-e|--vep-file              |text |Single VEP file (required either annovar-file or vep-file)                              |
|-b|--bam-files             |text |One or more bam files (required)                                                        |
|-m|--min-match-length      |int  |Minimum match length (optional default=15)                                              |
|-T|--pileup-tolerance      |int  |Pileup tolerance number of extra bases around the variant position (optional default=5) |
|-s|--min-hq-base-score     |int  |Minimum base score for filtered statistics (optional,default=30)                        |
|-S|--min-hq-alignment-score|int  |Minimum alignment score for filtered statistics (optional,default=40)                   |
|-r|--min-alignment-rate    |float|Minimum Smith-Waterman alignment rate (optional,default=0.95)                           |
|-t|--threads               |int  |Number of threads to use (optional default=1)                                           |
|-d|--count-duplicates      |void |If specified duplicates fragments are used in the statistics (optional default=false)   |
|-u|--count-secondary       |void |If specified secondary fragments are used in the statistics (optional default=false)    |
|-v|--verbose               |void |If specified be verbose (optional default=false)                                        |
|-h|--help                  |void |This help                                                                               |
