#include <algorithm>
#include <cassert>
#include <cctype>
#include <chrono>
#include <fstream>
#include <queue>
#include <regex>
#include <string>
#include <unordered_set>

#include "htslib_tagsort.h"
#include "metricgatherer.h"

constexpr int kLinesToReadInAChunk = 1000;

// Hands you lines one at a time from a file, reading in nice big efficient
// chunks whenever necessary.
class PartialFile
{
public:
  PartialFile(std::string const& filename) : file_(filename)
  {
    if (!file_)
      crash("ERROR failed to open the file " + filename);
    fillBuffer();
  }

  bool stillHasData()
  {
    if (sorted_.empty())
      fillBuffer();
    return !sorted_.empty();
  }

  // It is expected that you call stillHasData() immediately before every call
  // of takeNext(), to keep PartialFile's internal buffer populated.
  std::string takeNext()
  {
    std::string ret = sorted_.front();
    sorted_.pop();
    return ret;
  }

private:
  void fillBuffer()
  {
    int alignments_filled = 0;
    std::string line;
    while (std::getline(file_, line))
    {
      sorted_.push(line);
      if (++alignments_filled >= kLinesToReadInAChunk)
        break;
    }
  }

  std::ifstream file_;
  std::queue<std::string> sorted_;
};

// Helper class for merging several sorted lists into one, using
// std::priority_queue to keep track of which item to merge next.
class Merger
{
public:
  Merger() : heap_(greater_than_) {}

  // returns the appropriate next item, and an index into your PartialFile
  // array, from which you should giveInput() the next item into this Merger.
  // (If that PartialFile is now out of data, then you don't have to call
  //  giveInput() this time).
  std::pair<std::string, int> receiveNextOutput()
  {
    auto ret = heap_.top();
    heap_.pop();
    return ret;
  }

  void giveInput(std::string input, int src_file_ind)
  {
    heap_.emplace(input, src_file_ind);
  }

  bool empty() const { return heap_.empty(); }

private:
  std::function<bool(std::pair<std::string, int> const&,
                     std::pair<std::string, int> const&)>
      greater_than_ =
          [](std::pair<std::string, int> const& a,
             std::pair<std::string, int> const& b)
          {
            return a.first > b.first;
          };

  std::priority_queue<std::pair<std::string, int>,
                      std::vector<std::pair<std::string, int>>,
                      decltype(greater_than_) > heap_;
};

// Optional wrapper around MetricGatherer.
class MaybeMetrics
{
public:
  MaybeMetrics(INPUT_OPTIONS_TAGSORT const& options)
    : do_we_care_(options.compute_metric)
  {
    if (!do_we_care_)
      return;

    if (options.metric_type == MetricType::Cell)
    {
      mg_ = std::make_unique<CellMetricGatherer>(
          options.gtf_file, options.mitochondrial_gene_names_filename);
    }
    else if (options.metric_type == MetricType::Gene)
      mg_ = std::make_unique<GeneMetricGatherer>();
    else if (options.metric_type == MetricType::Umi)
      mg_ = std::make_unique<UmiMetricGatherer>();
    else
      crash("new MetricType enum value is not yet handled by MetricGatherer!");

    file_.open(options.metric_output_file);
    file_ << mg_->getHeader() << "\n";
  }

  void ingestLine(std::string line)
  {
    if (do_we_care_)
      mg_->ingestLine(line, file_);
  }

  void writeFinalMetrics()
  {
    if (do_we_care_)
    {
      mg_->output_metrics(file_);
      mg_->output_metrics_extra(file_);
    }
  }

private:
  const bool do_we_care_;
  std::unique_ptr<MetricGatherer> mg_{};
  std::ofstream file_;
};


// returns number of alignments processed
int mergePartialFiles(INPUT_OPTIONS_TAGSORT const& options,
                            std::vector<std::string> const& partial_files)
{
  int num_alignments = 0;

  MaybeMetrics maybe_metrics(options);
  std::ofstream fout;
  if (options.output_sorted_info)
    fout.open(options.sorted_output_file);

  // We are going to merge the contents of each of these into a single sorted file.
  std::vector<PartialFile> sorted_partial_files;
  for (std::string const& fname : partial_files)
    sorted_partial_files.emplace_back(fname);

  // We use a heap to track which item to put into the final sorted file next:
  // we start it with one item from each input file.
  Merger merger;
  for (int i = 0; i < sorted_partial_files.size(); i++)
    if (sorted_partial_files[i].stillHasData())
      merger.giveInput(sorted_partial_files[i].takeNext(), i);

  while (!merger.empty())
  {
    // 'line' is the next item to go to the sorted file.
    auto [line, i] = merger.receiveNextOutput();

    if (options.output_sorted_info)
      fout << line << "\n";
    maybe_metrics.ingestLine(line);
    num_alignments++;

    // Now that 'line' has been removed from the heap, it needs to be replaced
    // with another item from the same file (i) that it came from.
    // (If file i is out of data, then we don't need to replace).
    if (sorted_partial_files[i].stillHasData())
      merger.giveInput(sorted_partial_files[i].takeNext(), i);
  }

  maybe_metrics.writeFinalMetrics();

  if (options.output_sorted_info)
  {
    std::cout << "Wrote "<< num_alignments << " alignments in sorted order to "
              << options.sorted_output_file << std::endl;
  }
  return num_alignments;
}

void warnIfNo_mitochondrial_gene_names_filename(INPUT_OPTIONS_TAGSORT const& options)
{
  if (options.mitochondrial_gene_names_filename.empty())
  {
    std::string msg =
"*** WARNING! You did not specify --mitochondrial_gene_names_filename.\n"
"Therefore, we fell back to selecting only genes beginning with 'mt-' (case\n"
"insensitive). Please write a list of all gene names you're interested in into\n"
"a file, and pass the filename with --mitochondrial_gene_names_filename.";
    std::cout << msg << std::endl;
    std::cerr << msg << std::endl;
  }
}

/* Flag set by ‘--verbose’. */
int main(int argc, char** argv)
{
  INPUT_OPTIONS_TAGSORT options = readOptionsTagsort(argc, argv);
  warnIfNo_mitochondrial_gene_names_filename(options);

  std::cout << "bam input " << options.bam_input << std::endl;
  std::cout << "temp folder " << options.temp_folder << std::endl;
  std::cout << "sorted output file " <<  options.sorted_output_file << std::endl;
  std::cout << "metric output file " <<  options.metric_output_file << std::endl;
  std::cout << "temp folder " << options.alignments_per_batch << std::endl;
  std::cout << "tags:" << std::endl;

  for (auto const& [tag, tag_order_num] : options.tag_order)
    std::cout << "\t" << tag << "\t" << tag_order_num << std::endl;

  /* first create a list of sorted, and simplified sorted files */
  std::vector<std::string> partial_files = create_sorted_file_splits_htslib(options);

  /* now merge the sorted files to create one giant sorted file by using
    a head to compare the values based on the tags used  */
  std::cout << "Merging " <<  partial_files.size() << " sorted files!"<< std::endl;

  int num_alignments = mergePartialFiles(options, partial_files);
  std::cout << "Processed " << num_alignments << " alignments." << std::endl;

  // we no longer need the partial files
  for (unsigned int i=0; i < partial_files.size(); i++)
    if (remove(partial_files[i].c_str()) != 0)
      std::cerr << "Warning: error deleting file " << partial_files[i] << std::endl;

  warnIfNo_mitochondrial_gene_names_filename(options);
  return 0;
}
