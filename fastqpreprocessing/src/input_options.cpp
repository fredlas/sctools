/**
 *  @file   input_options.cpp
 *  @brief  functions for optons and input checking
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/
#include "input_options.h"

#include <experimental/filesystem>
#include <regex>
#include <getopt.h>
#include <cassert>
#include <iostream>
#include <cmath>

namespace fs = std::experimental::filesystem;
using std::string;

int64_t filesize(string const& filename)
{
  FILE* f = fopen(filename.c_str(), "rb");

  int64_t size = 0;
  if (fseek(f, 0, SEEK_END) == 0)
    size = ftell(f);
  fclose(f);
  return size;
}

void printFileInfo(std::vector<string> const& fastqs,
                   string const& type)
{
  if (fastqs.size())
  {
    std::cout << "INFO " << type << " files:" << std::endl;
    for (unsigned int i= 0; i < fastqs.size(); i++)
    {
      if (fs::exists(fastqs[i].c_str()))
      {
        std::cout << "\t " << fastqs[i]  <<  " exists, file size "
                  <<  filesize(fastqs[i])  <<  std::endl;
      }
      else
      {
        std::cout << "ERROR " << fastqs[i] << " is missing!\n";
        std::cerr << "ERROR " << fastqs[i] << " is missing!\n";
        exit(1);
      }
    }
  }
}

int64_t get_num_blocks(std::vector<string> const& I1s,
                     std::vector<string> const& R1s,
                     std::vector<string> const& R2s, double bam_size)
{
  assert(R1s.size() == R2s.size());
  double tot_size = 0;
  for (unsigned int i = 0; i < R1s.size(); i++)
  {
    assert(I1s.empty() || I1s.size() == R1s.size());
    if (!I1s.empty())
      tot_size += filesize(I1s[i]);

    std::cout << "file " << R1s[i] << " : " << filesize(R1s[i]) << " bytes" << std::endl;
    tot_size += filesize(R1s[i]);
    tot_size += filesize(R2s[i]);
  }

  const int GiB = 1024*1024*1024;
  return std::ceil((tot_size / GiB) / bam_size);
}

int64_t get_num_blocks(InputOptionsFastqProcess const& options)
{
  return get_num_blocks(options.I1s, options.R1s, options.R2s, options.bam_size);
}

int64_t get_num_blocks(INPUT_OPTIONS_FASTQ_READ_STRUCTURE const& options)
{
  return get_num_blocks(options.I1s, options.R1s, options.R2s, options.bam_size);
}

/** @copydoc readOptionsTagsort */
INPUT_OPTIONS_TAGSORT readOptionsTagsort(int argc, char** argv)
{
  INPUT_OPTIONS_TAGSORT options;
  int c;
  int i;

  static struct option long_options[] =
  {
    /* These options set a flag. */
    {"compute-metric",             no_argument,       0, 'm'},
    {"output-sorted-info",         no_argument,       0, 'n'},
    /* These options don’t set a flag.
       We distinguish them by their indices. */
    {"bam-input",                  required_argument, 0, 'b'},
    {"gtf-file",                   required_argument, 0, 'a'},
    {"temp-folder",                required_argument, 0, 't'},
    {"sorted-output",              required_argument, 0, 'o'},
    {"metric-output",              required_argument, 0, 'M'},
    {"alignments-per-thread",      required_argument, 0, 'p'},
    {"nthreads",                   required_argument, 0, 'T'},
    {"barcode-tag",                required_argument, 0, 'C'},
    {"umi-tag",                    required_argument, 0, 'U'},
    {"gene-tag",                   required_argument, 0, 'G'},
    {"metric-type",                required_argument, 0, 'K'},
    {"mitochondrial-gene-names-filename", required_argument, 0, 'g'},
    {0, 0, 0, 0}
  };

  // help messages when the user types -h
  const char* help_messages[] =
  {
    "compute metric, metrics are computed if this option is provided [optional]",
    "sorted output file is produced if this option is provided [optional]",
    "input bam file [required]",
    "gtf file (unzipped) required then metric type is cell [required with metric cell]",
    "temp folder for disk sorting [options: default /tmp]",
    "sorted output file [optional]",
    "metric file, the metrics are output in this file  [optional]",
    "number of alignments per thread [optional: default 1000000], if this number is increased then more RAM is required but reduces the number of file splits",
    "number of threads [optional: default 1]",
    "barcode-tag the call barcode tag [required]",
    "umi-tag the umi tag [required]: the tsv file output is sorted according the tags in the options barcode-tag, umi-tag or gene-tag",
    "gene-tag the gene tag [required]",
    "metric type, either \"cell\" or \"gene\" [required]",
    "file listing gene names, one per line, that the program should care about. [required, may omit if you want mouse or human]"
  };

  string metric_type_str;

  /* getopt_long stores the option index here. */
  int option_index = 0;
  int curr_size = 0;
  while ((c = getopt_long(argc, argv,
                          "b:a:t:no:mM:p:T:C:U:G:K:",
                          long_options,
                          &option_index)) !=- 1)
  {
    // process the option or arguments
    switch (c)
    {
    case 'm':
      options.compute_metric = true;
      break;
    case 'n':
      options.output_sorted_info = true;
      break;
    case 0:
      /* If this option set a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
        break;
      printf("option %s", long_options[option_index].name);
      if (optarg)
        printf(" with arg %s", optarg);
      printf("\n");
      break;
    case 'b':
      options.bam_input = string(optarg);
      break;
    case 'a':
      options.gtf_file = string(optarg);
      break;
    case 't':
      options.temp_folder = string(optarg);
      break;
    case 'o':
      options.sorted_output_file = string(optarg);
      break;
    case 'M':
      options.metric_output_file = string(optarg);
      break;
    case 'p':
      options.alignments_per_batch = atoi(optarg);
      break;
    case 'T':
      options.nthreads = atoi(optarg);
      break;
    case 'C':
      options.barcode_tag = string(optarg);
      curr_size = options.tag_order.size();
      options.tag_order[string(optarg)] = curr_size;
      break;
    case 'U':
      options.umi_tag = string(optarg);
      curr_size = options.tag_order.size();
      options.tag_order[string(optarg)] = curr_size;
      break;
    case 'G':
      options.gene_tag = string(optarg);
      curr_size = options.tag_order.size();
      options.tag_order[string(optarg)] = curr_size;
      break;
    case 'K':
      metric_type_str = string(optarg);
      break;
    case 'g':
      options.mitochondrial_gene_names_filename = string(optarg);
      break;
    case '?':
    case 'h':
      i = 0;
      printf("Usage: %s [options] \n", argv[0]);
      while (long_options[i].name != 0)
      {
        printf("\t--%-20s  %-25s  %-35s\n", long_options[i].name,
               long_options[i].has_arg == no_argument?
               "no argument" : "required_argument",
               help_messages[i]);
        i = i + 1;
      }
      /* getopt_long already printed an error message. */
      exit(0);
      break;
    default:
      abort();
    }
  }

  // Check the options
  // either metric computation or the sorted tsv file must be produced
  if (!options.output_sorted_info && !options.compute_metric)
    crash("ERROR: The choice of either the sorted alignment info or metric computation must be specified");

  if (options.compute_metric && options.metric_output_file.empty())
    crash("ERROR: Must specify --metric-output when specifying --compute-metric");

  if (options.output_sorted_info && options.sorted_output_file.empty())
    crash("ERROR: Must specify --sorted-output when specifying --output-sorted-info");

  if (metric_type_str == "cell")
    options.metric_type = MetricType::Cell;
  else if (metric_type_str == "gene")
    options.metric_type = MetricType::Gene;
  else if (metric_type_str == "umi")
    options.metric_type = MetricType::Umi;
  else
    crash("ERROR: --metric-type must be \"cell\", \"gene\", or \"umi\"");

  // TODO TODO also needed for umi?
  // if metric type is cell then the gtf file must be provided
  if (options.metric_type == MetricType::Cell && options.gtf_file.empty())
    crash("ERROR: The gtf file name must be provided with metric_type \"cell\"");

  // the gtf file should not be gzipped
  std::regex reg1(".gz$", std::regex_constants::icase);
  if (std::regex_search(options.gtf_file, reg1))
    crash("ERROR: The gtf file must not be gzipped");

  // bam input file must be there
  if (options.bam_input.empty())
    crash("ERROR: Must specify a input file name");

  // check for input file
  if (!fs::exists(options.bam_input.c_str()))
    crash("ERROR: bam_input " + options.bam_input + " is missing!");

  // check for the temp folder
  if (!fs::exists(options.temp_folder.c_str()))
    crash("ERROR: temp folder " + options.temp_folder + " is missing!");

  // check for three distinct tags, barcode, umi and gene_id tags
  if (options.tag_order.size() != 3)
    crash("ERROR:  Must have three distinct tags");
  bool seen_tag_index[3] = { false, false, false };
  for (auto [tag, index] : options.tag_order)
  {
    if (index < 0 || index > 2)
      crash("Invalid tag index " + std::to_string(index) + "; must be 0 1 or 2");
    else
      seen_tag_index[index] = true;
  }
  if (!(seen_tag_index[0] && seen_tag_index[1] && seen_tag_index[2]))
    crash("Need tag indices 0 1 and 2");

  // The size of a set of aligments for in-memory sorting must be positive
  if (options.alignments_per_batch < 1000)
    crash("ERROR: The number of alignments per thread must be at least 1000");

  // The number of threads must be between 1 and kMaxTagsortThreads
  if (options.nthreads > kMaxTagsortThreads || options.nthreads < 1)
    crash("ERROR: The number of threads must be between 1 and " + std::to_string(kMaxTagsortThreads));

  return options;
}


/** @copydoc readOptionsFastqProcess */
InputOptionsFastqProcess readOptionsFastqProcess(int argc, char** argv)
{
  InputOptionsFastqProcess options;
  int c;
  int i;
  bool verbose_flag = false;

  static struct option long_options[] =
  {
    /* These options set a flag. */
    {"verbose",           no_argument,       0, 'v'},
    /* These options don’t set a flag.
       We distinguish them by their indices. */
    {"barcode-length",    required_argument, 0, 'b'},
    {"umi-length",        required_argument, 0, 'u'},
    {"bam-size",          required_argument, 0, 'B'},
    {"sample-id",         required_argument, 0, 's'},
    {"I1",                required_argument, 0, 'I'},
    {"R1",                required_argument, 0, 'R'},
    {"R2",                required_argument, 0, 'r'},
    {"white-list",        required_argument, 0, 'w'},
    {"output-format",     required_argument, 0, 'F'},
    {0, 0, 0, 0}
  };

  // help messages when the user types -h
  const char* help_messages[] =
  {
    "verbose messages  ",
    "barcode length [required]",
    "UMI length [required]",
    "output BAM file in GB [optional: default 1 GB]",
    "sample id [required]",
    "I1 [optional]",
    "R1 [required]",
    "R2 [required]",
    "whitelist (from cellranger) of barcodes [required]",
    "output-format : either FASTQ or BAM [required]",
  };


  /* getopt_long stores the option index here. */
  int option_index = 0;
  while ((c = getopt_long(argc, argv,
                          "b:u:B:s:I:R:r:w:F:v",
                          long_options,
                          &option_index)) !=- 1
        )
  {
    // process the option or arguments
    switch (c)
    {
    case 'v':
      verbose_flag = true;
      break;
    case 0:
      /* If this option set a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
        break;
      printf("option %s", long_options[option_index].name);
      if (optarg)
        printf(" with arg %s", optarg);
      printf("\n");
      break;
    case 'b':
      options.barcode_length = atoi(optarg);
      break;
    case 'u':
      options.umi_length = atoi(optarg);
      break;
    case 'B':
      options.bam_size = atof(optarg);
      break;
    case 's':
      options.sample_id = string(optarg);
      break;
    case 'I':
      options.I1s.push_back(string(optarg));
      break;
    case 'R':
      options.R1s.push_back(string(optarg));
      break;
    case 'r':
      options.R2s.push_back(string(optarg));
      break;
    case 'w':
      options.white_list_file = string(optarg);
      break;
    case 'F':
      options.output_format = string(optarg);
      break;
    case '?':
    case 'h':
      i = 0;
      printf("Usage: %s [options] \n", argv[0]);
      while (long_options[i].name != 0)
      {
        printf("\t--%-20s  %-25s  %-35s\n", long_options[i].name,
               long_options[i].has_arg == no_argument?
               "no argument" : "required_argument",
               help_messages[i]);
        i = i + 1;
      }
      /* getopt_long already printed an error message. */
      return options;
    default:
      abort();
    }
  }

  if ((options.R1s.size() != options.R2s.size()))
  {
    crash("ERROR: Unequal number of R1 and R2 fastq files in input: R1: " +
          std::to_string(options.R1s.size()) + ", R2: " + std::to_string(options.R2s.size()));
  }

  if (options.R1s.empty())
    crash("ERROR: No R1 file provided");

  if (options.I1s.size() != options.R1s.size() && !options.I1s.empty())
    crash("ERROR: Must provide as many I1 input files as R1 input files, or else no I1 input files at all.");

  if (options.bam_size <= 0)
    crash("ERROR: Size of a bam file (in GB) cannot be negative or 0");

  if (options.sample_id.empty())
    crash("ERROR: Must provide a sample id or name");

  if (options.output_format!="FASTQ" && options.output_format!="BAM")
    crash("ERROR: output-format must be either FASTQ or BAM");

  if (options.barcode_length <= 0)
    crash("ERROR: Barcode length must be a positive integer");

  if (options.umi_length <= 0)
    crash("ERROR: UMI length must be a positive integer");

  if (verbose_flag)
  {
    if (!options.I1s.empty())
      printFileInfo(options.I1s, string("I1"));
    if (!options.R1s.empty())
      printFileInfo(options.R1s, string("R1"));
    if (!options.R2s.empty())
      printFileInfo(options.R2s, string("R2"));
  }

  return options;
}

INPUT_OPTIONS_FASTQ_READ_STRUCTURE readOptionsFastqSlideseq(int argc, char** argv)
{
  INPUT_OPTIONS_FASTQ_READ_STRUCTURE options;
  int c;
  int i;
  bool verbose_flag = false;

  static struct option long_options[] =
  {
    /* These options set a flag. */
    {"verbose",           no_argument,       0, 'v'},
    /* These options don’t set a flag.
       We distinguish them by their indices. */
    {"bam-size",          required_argument, 0, 'B'},
    {"read-structure",    required_argument, 0, 'S'},
    {"sample-id",         required_argument, 0, 's'},
    {"I1",                required_argument, 0, 'I'},
    {"R1",                required_argument, 0, 'R'},
    {"R2",                required_argument, 0, 'r'},
    {"white-list",        required_argument, 0, 'w'},
    {"output-format",     required_argument, 0, 'F'},
    {0, 0, 0, 0}
  };

  // help messages when the user types -h
  const char* help_messages[] =
  {
    "verbose messages  ",
    "output BAM file in GB [optional: default 1 GB]",
    "read structure [required]",
    "sample id [required]",
    "I1 [optional]",
    "R1 [required]",
    "R2 [required]",
    "whitelist (from cellranger) of barcodes [required]",
    "output-format : either FASTQ or BAM [required]",
  };


  /* getopt_long stores the option index here. */
  int option_index = 0;
  while ((c = getopt_long(argc, argv,
                          "B:S:s:I:R:r:w:F:v",
                          long_options,
                          &option_index)) !=- 1)
  {
    // process the option or arguments
    switch (c)
    {
    case 'v':
      verbose_flag = true;
      break;
    case 0:
      /* If this option set a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
        break;
      printf("option %s", long_options[option_index].name);
      if (optarg)
        printf(" with arg %s", optarg);
      printf("\n");
      break;
    case 'B':
      options.bam_size = atof(optarg);
      break;
    case 'S':
      options.read_structure = string(optarg);
      break;
    case 's':
      options.sample_id = string(optarg);
      break;
    case 'I':
      options.I1s.push_back(string(optarg));
      break;
    case 'R':
      options.R1s.push_back(string(optarg));
      break;
    case 'r':
      options.R2s.push_back(string(optarg));
      break;
    case 'w':
      options.white_list_file = string(optarg);
      break;
    case 'F':
      options.output_format = string(optarg);
      break;
    case '?':
    case 'h':
      i = 0;
      printf("Usage: %s [options] \n", argv[0]);
      while (long_options[i].name != 0)
      {
        printf("\t--%-20s  %-25s  %-35s\n", long_options[i].name,
               long_options[i].has_arg == no_argument?
               "no argument" : "required_argument",
               help_messages[i]);
        i = i + 1;
      }
      /* getopt_long already printed an error message. */
      return options;
    default:
      abort();
    }
  }

  if ((options.R1s.size() != options.R2s.size()))
  {
    crash("ERROR: Unequal number of R1 and R2 fastq files in input: R1: " +
          std::to_string(options.R1s.size()) + ", R2: " + std::to_string(options.R2s.size()));
  }

  if (options.R1s.empty())
    crash("ERROR: No R1 file provided");

  if (options.I1s.size() != options.R1s.size() && !options.I1s.empty())
    crash("ERROR: Must provide as many I1 input files as R1 input files, or else no I1 input files at all.");

  if (options.bam_size <= 0)
    crash("ERROR: Size of a bam file (in GB) cannot be negative or 0");

  if (options.sample_id.empty())
    crash("ERROR: Must provide a sample id or name");

  if (options.output_format!="FASTQ" && options.output_format!="BAM")
    crash("ERROR: output-format must be either FASTQ or BAM");

  if (options.read_structure.empty())
    crash("ERROR: Must provide read structures");

  if (verbose_flag)
  {
    if (!options.R1s.empty())
      printFileInfo(options.R1s, string("R1"));
    if (!options.R2s.empty())
      printFileInfo(options.R2s, string("R2"));
  }

  return options;
}


INPUT_OPTIONS_FASTQ_READ_STRUCTURE readOptionsFastqMetrics(int argc, char** argv)
{
  INPUT_OPTIONS_FASTQ_READ_STRUCTURE options;
  int c;
  int i;
  bool verbose_flag = false;

  static struct option long_options[] =
  {
    /* These options set a flag. */
    {"verbose",           no_argument,       0, 'v'},
    /* These options don’t set a flag.
       We distinguish them by their indices. */
    {"read-structure",    required_argument, 0, 'S'},
    {"sample-id",         required_argument, 0, 's'},
    {"R1",                required_argument, 0, 'R'},
    {"white-list",        required_argument, 0, 'w'},
    {0, 0, 0, 0}
  };

  // help messages when the user types -h
  const char* help_messages[] =
  {
    "verbose messages  ",
    "read structure [required]",
    "sample id [required]",
    "R1 [required]",
    "whitelist of cell/bead barcodes [required]",
  };


  /* getopt_long stores the option index here. */
  int option_index = 0;
  while ((c = getopt_long(argc, argv,
                          "S:s:R:w:v",
                          long_options,
                          &option_index)) !=- 1
        )
  {
    // process the option or arguments
    switch (c)
    {
    case 'v':
      verbose_flag = 1;
      break;
    case 0:
      /* If this option set a flag, do nothing else now. */
      if (long_options[option_index].flag != 0)
        break;
      printf("option %s", long_options[option_index].name);
      if (optarg)
        printf(" with arg %s", optarg);
      printf("\n");
      break;
    case 'S':
      options.read_structure = string(optarg);
      break;
    case 's':
      options.sample_id = string(optarg);
      break;
    case 'R':
      options.R1s.push_back(string(optarg));
      break;
    case 'w':
      options.white_list_file = string(optarg);
      break;
    case '?':
    case 'h':
      i = 0;
      printf("Usage: %s [options] \n", argv[0]);
      while (long_options[i].name != 0)
      {
        printf("\t--%-20s  %-25s  %-35s\n", long_options[i].name,
               long_options[i].has_arg == no_argument?
               "no argument" : "required_argument",
               help_messages[i]);
        i = i + 1;
      }
      /* getopt_long already printed an error message. */
      return options;
    default:
      abort();
    }
  }

  if (options.R1s.empty())
    crash("ERROR: No R1 file provided");

  if (options.read_structure.empty())
    crash("ERROR: Must provide read structures");

  if (options.sample_id.empty())
    crash("ERROR: Must provide a sample id or name");

  if (verbose_flag && !options.R1s.empty())
    printFileInfo(options.R1s, string("R1"));

  return options;
}
