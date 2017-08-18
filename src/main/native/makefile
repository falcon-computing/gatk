CXX = g++
CXXFLAGS = -std=c++11 -O3 -fPIC -pthread -Wall

JNI_FLAGS = -I${JAVA_HOME}/include -I${JAVA_HOME}/include/linux
JNI_LDFLAGS = -L${JAVA_HOME}/jre/lib/amd64/server -ljvm

BWA_FLAGS = -I$(BWA_PATH)
BWA_FLAGS += -I$(BWA_FLOW_PATH) -I$(BWA_FLOW_PATH)/src
BWA_FLAGS += -DUSE_HTSLIB 	# Current Java code requires BAM output from bwa,
							# which requires the USE_HTSLIB macro.

HTSLIB_FLAGS = -I$(HTSLIB_PATH)
HTSLIB_LDFLAGS = -L$(HTSLIB_PATH) -lhts

ZLIB_LDFLAGS = -lz

HDFS_FLAGS = -I$(HADOOP_PATH)/include
HDFS_LDFLAGS = -L$(HADOOP_PATH)/lib/native -lhdfs

BWA_JNI_ENTRY_CLASS = org.broadinstitute.hellbender.tools.spark.bwa.NativeBwaSparkEngine
BWA_JNI_ENTRY = $(subst .,_,$(BWA_JNI_ENTRY_CLASS))
BWA_LIB_NAME = gatkbwa
BWA_FLOW_PATH = bwa-flow
BWA_PATH = $(BWA_FLOW_PATH)/bwa

HTSLIB_PATH ?= htslib
CXXOPTS_PATH ?= cxxopts
HADOOP_PATH ?= hadoop

SPLIT_JNI_ENTRY_CLASS = org.broadinstitute.hellbender.tools.spark.bwa.NativeSplitFastqSparkEngine
SPLIT_JNI_ENTRY = $(subst .,_,$(SPLIT_JNI_ENTRY_CLASS))
SPLIT_LIB_NAME = splitfastq

.PHONY: all, clean

all: lib/lib$(BWA_LIB_NAME).so lib/lib$(SPLIT_LIB_NAME).so bin/split-fastq

clean:
	$(MAKE) -C $(BWA_PATH) $@
	rm -rf obj bin lib "src/$(BWA_JNI_ENTRY).h"

# srcs
src/$(BWA_JNI_ENTRY).h:
	javah -cp ../../../build/install/gatk/lib/\* -d $(@D) $(BWA_JNI_ENTRY_CLASS)
src/$(SPLIT_JNI_ENTRY).h:
	javah -cp ../../../build/install/gatk/lib/\* -d $(@D) $(SPLIT_JNI_ENTRY_CLASS)

# objs
obj/split-fastq.o: src/split-fastq.cpp
	@mkdir -p $(@D)
	$(CXX) -c $(CXXFLAGS) $(HDFS_FLAGS) $< -o $@

obj/split-fastq-main.o: src/split-fastq-main.cpp
	@mkdir -p $(@D)
	$(CXX) -c $(CXXFLAGS) -I$(CXXOPTS_PATH)/include $< -o $@

obj/$(BWA_JNI_ENTRY).o: src/$(BWA_JNI_ENTRY).cpp src/$(BWA_JNI_ENTRY).h $(BWA_FLOW_PATH)/src/bwa_wrapper.h
	@mkdir -p $(@D)
	$(CXX) -c $(CXXFLAGS) $(JNI_FLAGS) $(BWA_FLAGS) $(HTSLIB_FLAGS) $(HDFS_FLAGS) $< -o $@
	
obj/preprocess.o: $(BWA_FLOW_PATH)/src/preprocess.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(BWA_FLAGS) $(HTSLIB_FLAGS) -MM -MP -MT $@ -MF $(@:.o=.d) $<
	$(CXX) -c $(CXXFLAGS) $(BWA_FLAGS) $(HTSLIB_FLAGS) $< -o $@

-include obj/preprocess.d

obj/$(SPLIT_JNI_ENTRY).o: src/$(SPLIT_JNI_ENTRY).cpp src/$(SPLIT_JNI_ENTRY).h
	@mkdir -p $(@D)
	$(CXX) -c $(CXXFLAGS) $(JNI_FLAGS) $(HDFS_FLAGS) $< -o $@

# libs
lib/lib$(BWA_LIB_NAME).so: obj/$(BWA_JNI_ENTRY).o obj/preprocess.o $(BWA_PATH)/libbwa.a
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(HTSLIB_LDFLAGS) $(HDFS_LDFLAGS) $(JNI_LDFLAGS) -shared $^ -o $@

$(BWA_PATH)/libbwa.a:
	$(MAKE) -C $(BWA_PATH) libbwa.a

lib/lib$(SPLIT_LIB_NAME).so: obj/$(SPLIT_JNI_ENTRY).o obj/split-fastq.o
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(HDFS_LDFLAGS) $(JNI_LDFLAGS) $(ZLIB_LDFLAGS) -shared $^ -o $@

# bins
bin/split-fastq: obj/split-fastq.o obj/split-fastq-main.o
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(ZLIB_LDFLAGS) $(HDFS_LDFLAGS) $(JNI_LDFLAGS) $^ -o $@

