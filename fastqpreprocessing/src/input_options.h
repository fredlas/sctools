#ifndef __SCTOOLS_FASTQPREPROCESSING_INPUT_OPTIONS_H_
#define __SCTOOLS_FASTQPREPROCESSING_INPUT_OPTIONS_H_
/**
 *  @file   input_options.h
 *  @brief  Utility functions for input options processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/

#include "utilities.h"

#include <string>
#include <vector>

const unsigned int kDefaultNumAlignsPerThread = 1000000;

struct InputOptionsFastqReadStructure
{
  // I1, R1 and R2 files name
  std::vector<std::string> I1s, R1s, R2s;

  // Bead Barcode list
  std::string white_list_file;

  std::string output_format;

  // Bam file size to split by (in GB)
  double bam_size = 1.0;

  std::string read_structure;

  std::string sample_id;
};


// Structure to hold input options for fastqprocess
struct InputOptionsFastqProcess
{
  // I1, R1 and R2 files name
  std::vector<std::string> I1s, R1s, R2s;

  // Barcode white list file
  std::string white_list_file;

  std::string output_format;

  // chemistry dependent (V2/V3) barcode and UMI length
  int barcode_length = -1;
  int umi_length = -1;

  // Bam file size to split by (in GB)
  double bam_size = 1.0;

  std::string sample_id;
};


// Structure to hold input options for tagsort
struct InputOptionsTagsort
{
  // Initialize some of the values
  InputOptionsTagsort()
  {
    bam_input = "";
    gtf_file = "";
    temp_folder =  std::string("/tmp/");
    nthreads = 1;
    metric_type = "";
  }
  std::string metric_type;
  bool output_sorted_info = false;
  bool compute_metric = false;
  // name of the bam file
  std::string bam_input;
  // name of the gtf file
  std::string gtf_file;
  // temp folder for disk sorting
  std::string temp_folder = "/tmp/";

  std::string metric_output_file;
  // sorted tsv output file
  std::string sorted_output_file;

  unsigned int alignments_per_thread = kDefaultNumAlignsPerThread;
  unsigned int nthreads;
  std::string barcode_tag;
  std::string umi_tag;
  std::string gene_tag;

  // order of the tags to sort by
  std::unordered_map<std::string, unsigned int> tag_order;
};

InputOptionsFastqProcess readOptionsFastqProcess(int argc, char** argv);

InputOptionsTagsort readOptionsTagsort(int argc, char** argv);

InputOptionsFastqReadStructure readOptionsFastqSlideseq(int argc, char** argv);

InputOptionsFastqReadStructure readOptionsFastqMetrics(int argc, char** argv);

int64_t getNumBlocks(InputOptionsFastqProcess const& options);
int64_t getNumBlocks(InputOptionsFastqReadStructure const& options);

#endif // __SCTOOLS_FASTQPREPROCESSING_INPUT_OPTIONS_H_
