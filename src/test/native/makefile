.PHONY: test

SHELL=bash

test: bin/gatk-bwa
	diff ../resources/ERR000097_part.sam <($^ ../resources/human_g1k_v37.chr17_1Mb.fasta ../resources/ERR000097_part_1.fastq.gz ../resources/ERR000097_part_2.fastq.gz) -q && echo 'gatkbwa:TEST Pass!'

-include ../../main/native/makefile

bin/gatk-bwa: src/gatk-bwa.cpp ../../main/native/obj/org_broadinstitute_hellbender_tools_spark_bwa_NativeBwaSparkEngine.o ../../main/native/obj/preprocess.o ../../main/native/$(BWA_PATH)/libbwa.a
	@mkdir -p $(@D)
	cd ../../main/native;$(CXX) $(CXXFLAGS) $(BWA_FLAGS) $(HTSLIB_FLAGS) $(HTSLIB_LDFLAGS) $(ZLIB_LDFLAGS) $(HDFS_LDFLAGS) $(JNI_LDFLAGS) ../../test/native/$^ -o ../../test/native/$@


