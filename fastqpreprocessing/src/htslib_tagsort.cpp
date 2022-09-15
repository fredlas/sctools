/**
 *  @file   htslib_tagsort.cpp
 *  @brief  functions for file processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/

constexpr int kThreshold = 30; // qual score threshold

#include "htslib_tagsort.h"

#include <algorithm>
#include <tuple>
#include <cstdint>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <thread>
#include <mutex>

#include <cassert>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <unordered_map>
#include <memory>
#include <set>

extern "C" {
  bam_hdr_t* sam_hdr_read(samFile*);   //read header
  htsFile* hts_open(const char* fn, const char* mode);
}


/*
  @brief get the int tag or -1

*/
inline int get_itag_or_default(bam1_t* aln, const char* tagname, int default_value)
{
  uint8_t* p;
  int tag_value = -1;
  if ((p = bam_aux_get(aln, tagname)) == nullptr)
    tag_value = default_value;
  else
    tag_value = bam_aux2i(p);

  return  tag_value;
}

/*
  @brief get the string tag or the default

*/
inline char* get_Ztag_or_default(bam1_t* aln, const char* tagname, char* default_value)
{
  uint8_t* p;
  char* tag_value = nullptr;
  if ((p = bam_aux_get(aln, tagname)) == nullptr)
    tag_value = default_value;
  else
  {
    tag_value = bam_aux2Z(p);
    if (strcmp(tag_value, "-") == 0)
      tag_value = default_value;
  }
  return  tag_value;
}

using TAGTUPLE = std::tuple<
    TRIPLET /*  barcode umi and gene_id, not necessarily in that order */,
    std::string /* reference */, TODO also stringbank this?
    std::string /* biotype */, TODO also stringbank this?
    int /* pos */,
    int /*rev strand   1 for yes, 0 otherwise*/,
    float /*avg barcode qual score */,
    float /* frac of barcode qual score >30 */,
    float /*avg qual seq */,
    float /*fract of >30 score qual seq*/,
    int /*NH*/,
    int /*perfect molecule barcode, 1 is yes, 0 otherwise*/,
    int /*spliced reads 1 yes, 0 otherwise*/,
    int /*is duplicate */,
    int /*perfect cell barcode 1 is yes, 0 otherwise*/,
    float /* fraction of umi qual score > 30 */
    >;

void parseOneAlignment(std::vector<TAGTUPLE>* tuple_records, bam1_t* aln,
                       INPUT_OPTIONS_TAGSORT& options, const bam_hdr_t* bam_hdr,
                       TagTripletManager& triplet_manager)
{
  // "consts" that the library doesn't allow to be const.
  char empty[] = "";
  char none[] = "None";
  char nochr[] = "*";

  // extract the barcodes corrected and  corrected
  char* barcode = get_Ztag_or_default(aln, options.barcode_tag.c_str(), none);
  char* barcode_raw = get_Ztag_or_default(aln, "CR", empty);

  // to be called perfect, the corrected and raw barcodes should match
  int perfect_cell_barcode = (strcmp(barcode, barcode_raw) == 0) ? 1 : 0;

  // barcode quality score
  char* barcode_qual = get_Ztag_or_default(aln, "CY", empty);

  //average barcode across the query and the fraction of barcodes above threshold
  float sum_barcode_qual = 0;
  float num_bp_above_threshold = 0;
  size_t len = strlen(barcode_qual);
  for (unsigned int k = 0; k < len; k++)
  {
    // barcodes qual strings are in ASCII symbols subtracting 33 gives the phred qual score
    uint8_t qual_score = (((uint8_t)barcode_qual[k]) - 33);
    sum_barcode_qual += qual_score;
    if (qual_score > kThreshold)
      num_bp_above_threshold += 1;
  }
  int avg_cell_barcode_qual = sum_barcode_qual / (float)len; // TODO truncation intended?
  float cell_barcode_qual_above_threshold = (float)num_bp_above_threshold / (float)len;

  // corrected molecule barcodes (UMIs)
  char* umi = get_Ztag_or_default(aln, options.umi_tag.c_str(), none);
  // raw molecule barcodes
  char* umi_raw = get_Ztag_or_default(aln, "UR", empty);

  // to be called perfect, the corrected and raw molecular barcodes should match
  int perfect_molecule_barcode = (strcmp(umi, umi_raw) == 0) ? 1 : 0;

  // qual score for molecular barcodes
  char* umi_qual = get_Ztag_or_default(aln, "UY", empty);

  float sum_umi_qual = 0;
  float num_umi_above_threshold = 0;
  len = strlen(umi_qual);
  for (unsigned int k = 0; k < len; k++)
  {
    // molecular barcodes qual strings are in ASCII symbols subtracting 33 gives the phred qual score
    sum_umi_qual += ((uint8_t)umi_qual[k] -33);
    if (((uint8_t)umi_qual[k] - 33) > kThreshold)
      num_umi_above_threshold += 1;
  }
  float frac_umi_qual_above_threshold = (float)num_umi_above_threshold / (float)len;

  char* gene_id = get_Ztag_or_default(aln, options.gene_tag.c_str(), none);
  char* location_tag = get_Ztag_or_default(aln, "XF", empty);

  int nh_num = get_itag_or_default(aln, "NH", -1);

  const char* chr = (aln->core.tid == -1) ? nochr : bam_hdr->target_name[aln->core.tid];

  uint32_t pos = aln->core.pos; // position.
  uint32_t isrev = bam_is_rev(aln) ? 1 : 0;   // is reverse stand
  uint32_t is_duplicate = ((aln->core.flag & BAM_FDUP) != 0) ? 1 : 0;

  // sequence quality score
  float avg_sequence_qual = 0, sum_qual = 0;
  float qual_above_threshold = 0;
  uint8_t* qual_seq = bam_get_qual(aln);  // pointer to the qual data
  len = aln->core.l_qseq; //length of qual seq.
  for (unsigned int k = 0; k < len; k++)
  {
    // the qual string are already in phred scores
    sum_qual += qual_seq[k];
    if (qual_seq[k] > kThreshold)
      qual_above_threshold += 1;
  }
  avg_sequence_qual = sum_qual / (float)len;
  qual_above_threshold = qual_above_threshold / (float)len;

  uint32_t* cigar = bam_get_cigar(aln);
  // see if it is spliced, i.e., N appears in the CIGAR string
  uint32_t spliced_read = 0;
  for (unsigned int k = 0; k < aln->core.n_cigar; k++)
  {
    uint32_t op = cigar[k] & BAM_CIGAR_MASK;
    if (op == 3 && (cigar[k] >> BAM_CIGAR_SHIFT) != 0)
    {
      spliced_read = 1;
      break;
    }
  }

  tuple_records->emplace_back(
      triplet_manager.makeTriplet(barcode, umi, gene_id), /* triplet of tags */
      std::string(chr),  /* record[0] */
      std::string(location_tag), /* record[1] */
      pos,   /* record [2] */
      isrev, /* record[3] */
      avg_cell_barcode_qual,  /* record[4] */
      cell_barcode_qual_above_threshold, /* record[5] */
      avg_sequence_qual, /* record[6] */
      qual_above_threshold, /* record[7]  */
      nh_num, /* record[8] */
      perfect_molecule_barcode, /* record[9] */
      spliced_read, /* record[10] */
      is_duplicate, /* record[11] */
      perfect_cell_barcode,  /* record[12] */
      frac_umi_qual_above_threshold /* record[13] */);
}

// Generates a random alphanumeric string (AZaz09) of a fixed length.
constexpr int kStringLen = 40;
std::string randomString()
{
  auto randchar = []() -> char
  {
    const char charset[] =
    "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz";
    const size_t max_index = (sizeof(charset) - 1);
    return charset[ rand() % max_index ];
  };
  std::string str(kStringLen, 0);
  std::generate_n(str.begin(), kStringLen, randchar);
  return str;
}

/**
 * @brief This function takes a vector of tuples of the tags, sorts them
 * in the dictionary order of the tags and then writes these in the same
 * order to a txt file
 *
 * @details
 * The function take the vector of tags tuples, writes the sorted tuples into
 * a file. The filename is generated randomly (enought to avoid collision with other files)
 * in the temp folder specified.
 *
 * @param tuple_records: vector<TAGTUPLE> &, reference to a vector of TAGTUPLES
 * @return a string for the random file name
*/
std::string sortAndWriteToPartialTxtFile(std::vector<TAGTUPLE> const& tuple_records,
                                         std::string const& tmp_folder,
                                         TagTripletManager& m)
{
  using std::get;

  std::string tempfile_filename = tmp_folder + "/" + randomString() + ".txt";
  std::ofstream outfile(tempfile_filename);

  // Sort by triplet, maintaining each triplet's pre-sorted index...
  std::vector<std::pair<TRIPLET, int>> index_pairs;
  for (size_t i = 0; i < tuple_records.size(); i++)
    index_pairs.emplace_back(get<0>(tuple_records[i]), i);

  std::sort(index_pairs.begin(), index_pairs.end(),
  [&m](std::pair<TRIPLET, int> const& a, std::pair<TRIPLET, int> const& b)
  {
    std::string const& a0 = m.getString0(a.first);
    std::string const& b0 = m.getString0(b.first);
    int comp0 = a0.compare(b0)
    if (comp0 != 0)
      return comp0 < 0;

    std::string const& a1 = m.getString1(a.first);
    std::string const& b1 = m.getString1(b.first);
    int comp1 = a1.compare(b1)
    if (comp1 != 0)
      return comp1 < 0;

    std::string const& a2 = m.getString2(a.first);
    std::string const& b2 = m.getString2(b.first);
    int comp2 = a2.compare(b2)
    return comp2 < 0;
  });

  // ...then write the triplets in sorted order, linked up with the tuple_record
  // located at its pre-sorted index, i.e. the record for the triplet.
  for (auto& [triplet, record_index] : index_pairs)
  {
    // TODO?
    // what if you ran out of disk space ???? NEED TO add logic
    outfile << m.getString0(triplet) /*  first tag */ << "\t"
            << m.getString1(triplet) /*  second tag */ << "\t"
            << m.getString2(triplet) /*  third tag */ << "\t"
            << get<1>(tuple_records[record_index]) /* record[0] */  << "\t"
            << get<2>(tuple_records[record_index]) /* record[1] */  << "\t"
            << get<3>(tuple_records[record_index]) /* record[2] */  << "\t"
            << get<4>(tuple_records[record_index]) /* record[3] */  << "\t"
            << get<5>(tuple_records[record_index]) /* record[4] */  << "\t"
            << get<6>(tuple_records[record_index]) /* record[5] */  << "\t"
            << get<7>(tuple_records[record_index]) /* record[6] */  << "\t"
            << get<8>(tuple_records[record_index]) /* record[7] */  << "\t"
            << get<9>(tuple_records[record_index]) /* record[8] */  << "\t"
            << get<10>(tuple_records[record_index]) /* record[9] */ << "\t"
            << get<11>(tuple_records[record_index]) /* record[10] */ << "\t"
            << get<12>(tuple_records[record_index]) /* record[11] */ << "\t"
            << get<13>(tuple_records[record_index]) /* record[12] */ << "\t"
            << get<14>(tuple_records[record_index]) /* record[13] */ << "\n";
  }

  return tempfile_filename;
}

// Manages worker threads' access to reading the input file. Any worker can take
// any line from the file, but only one thread can be reading at a time. This
// class lets them take turns: they call readAlignments() whenever they want
// more data, and it blocks until they can have it.
class AlignmentReader
{
public:
  explicit AlignmentReader(INPUT_OPTIONS_TAGSORT options) : options_(options)
  {
    if (options_.alignments_per_batch > 4294967295)
      crash("alignments_per_batch too large; max allowed is 4294967295");

    if ((sam_file_ptr_ = hts_open(options.bam_input.c_str(),"r")) == nullptr)
      crash(options.bam_input + ": cannot open file.");

    bam_hdr_ = sam_hdr_read(sam_file_ptr_); //read header

    allocateAlignmentBuffers();
  }

  ~AlignmentReader()
  {
    for (unsigned int i = 0; i < options_.nthreads; i++)
    {
      for (unsigned int k = 0; k < options_.alignments_per_batch; k++)
        bam_destroy1(aln_arr_[i][k]);
      free(aln_arr_[i]);
    }
    free(aln_arr_);
    sam_hdr_destroy(bam_hdr_);
    hts_close(sam_file_ptr_);
  }

  // Blocks until it's this thread's turn to read, then reads a batch of alignments.
  // Returns a pointer to the alignment ptr array, and number of alignment ptrs in the array.
  std::pair<bam1_t**, unsigned int> readAlignments(int thread_index)
  {
    const std::lock_guard<std::mutex> lock(mutex_);
    unsigned int cur_num_read = 0;
    while (cur_num_read < options_.alignments_per_batch)
    {
      if (sam_read1(sam_file_ptr_, bam_hdr_, aln_arr_[thread_index][cur_num_read]) == 0)
        cur_num_read++;
      else
        break;
    }
    total_aligns_read_ += cur_num_read;
    batches_read_++;
    std::cout << "Finished reading batch number: " << batches_read_ << std::endl;
    return std::make_pair(aln_arr_[thread_index], cur_num_read);
  }

  void addToPartialFilenames(std::vector<std::string> names)
  {
    const std::lock_guard<std::mutex> lock(mutex_);
    for (std::string name : names)
      partial_filenames_.push_back(name);
  }

  std::vector<std::string> partial_filenames() const { return partial_filenames_; }
  bam_hdr_t* bam_hdr() const { return bam_hdr_; }
  uint64_t total_aligns_read() const { return total_aligns_read_; }

private:
  void allocateAlignmentBuffers()
  {
    assert(options_.nthreads <= kMaxTagsortThreads);
    std::string msg = "Now allocating alignment buffers. If the program crashes "
                      "here, it probably ran out of memory...";
    std::cout << msg << std::endl;
    std::cerr << msg << std::endl;

    aln_arr_ = (bam1_t***)malloc(sizeof(bam1_t**) * options_.nthreads);
    for (unsigned int i = 0; i < options_.nthreads; i++)
    {
      aln_arr_[i] = (bam1_t**)malloc(sizeof(bam1_t*) * options_.alignments_per_batch);
      for (unsigned int k = 0; k < options_.alignments_per_batch; k++)
        aln_arr_[i][k] = bam_init1(); //initialize an alignment
    }
    std::string done_msg = "Successfully allocated alignment buffers.";
    std::cout << done_msg << std::endl;
    std::cerr << done_msg << std::endl;
  }

  std::mutex mutex_;
  uint64_t total_aligns_read_ = 0;
  uint64_t batches_read_ = 0;
  INPUT_OPTIONS_TAGSORT options_;
  samFile* sam_file_ptr_ = nullptr;
  bam_hdr_t* bam_hdr_ = nullptr;
  bam1_t*** aln_arr_ = nullptr;
  std::vector<std::string> partial_filenames_;
};

void partialSortWorkerThread(int my_thread_index, AlignmentReader* alignment_reader,
                             TagOrder tag_order, INPUT_OPTIONS_TAGSORT options)
{
  std::vector<std::string> my_partial_filenames;
  bam_hdr_t* bam_hdr = alignment_reader->bam_hdr();
  while (true)
  {
    TagTripletManager triplet_manager(options); TODO rebuild per chunk, or keep across chunks? (will strings be repeated?)
    std::vector<TAGTUPLE> tuple_records;

    auto [aln_ptr_array, alns_length] = alignment_reader->readAlignments(my_thread_index);
    if (alns_length == 0)
      break;

    for (unsigned int i = 0; i < alns_length; i++)
      parseOneAlignment(&tuple_records, aln_ptr_array[i], options, bam_hdr, triplet_manager);

    my_partial_filenames.push_back(
        sortAndWriteToPartialTxtFile(tuple_records, options.temp_folder, triplet_manager));
  }
  alignment_reader->addToPartialFilenames(my_partial_filenames);
}

/**
 * @brief From the input bam create a list of txt files with the records (lines)
 * sorted according to the * tags
 *
 * @details
 * The input bam file is read chunk by chunk, sorted by the tags and the written
 * out as a text file in the sorted manner.
 *
 * @param options: INPUT_OPTIONS_TAGSORT the inputs to the program
 * @return a vector containing the file paths of the partial files
*/
std::vector<std::string> create_sorted_file_splits_htslib(INPUT_OPTIONS_TAGSORT options)
{
  std::cout << "Running htslib" << std::endl;
  AlignmentReader alignment_reader(options);

  TagOrder tag_order = getTagOrder(options);

  std::vector<std::thread> worker_threads;
  for (int i = 0; i < options.nthreads; i++)
  {
    worker_threads.emplace_back(partialSortWorkerThread,
                                i, &alignment_reader, tag_order, options);
  }
  for (auto& worker_thread : worker_threads)
    worker_thread.join();

  std::cout << "Read " << alignment_reader.total_aligns_read() << " records in batches of "
            << options.alignments_per_batch << std::endl;

  return alignment_reader.partial_filenames();
}  // function

