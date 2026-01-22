cmake_minimum_required(VERSION 3.5)

include_directories(htslib/include)
link_directories(htslib/lib)

add_compile_options(-Wall -mavx2)

project(enhanced_ABS LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(enhanced_ABS main.cpp annotate_bam_statistics.cpp annotate_bam_statistics.h annovar_file.cpp annovar_file.h vep_file.cpp vep_file.h variant_entry.h fasta_file.cpp fasta_file.h bam_file.cpp bam_file.h vcf_file.cpp vcf_file.h ssw.cpp ssw.h progress_bar.h)

find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
target_link_libraries (enhanced_ABS PRIVATE libparasail.a libhts.a z m bz2 lzma curl crypto deflate Threads::Threads OpenMP::OpenMP_CXX)
