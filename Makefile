BAM_LIB=-L/home/bdownie/src/bamtools/lib
BAM_INC=-I/home/bdownie/src/bamtools/include/
GCC=g++ -O3 -Wall

remove_3prime_5prime_variants: remove_3prime_5prime_variants.cpp
	$(GCC) $(BAM_LIB) $(BAM_INC) -lz -lbamtools -o remove_3prime_5prime_variants remove_3prime_5prime_variants.cpp
