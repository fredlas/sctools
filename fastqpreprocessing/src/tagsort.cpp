/**
 *  @file   tagsort.cpp
 *  @brief  functions for file processing
 *  @author Kishori M. Konwar
 *  @date   2021-08-11
 ***********************************************/

#include "htslib_tagsort.h"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <chrono>
#include <fstream>
#include <priority_queue>
#include <regex>
#include <string>

constexpr int kDataBufferSize = 1000;

struct Context
{
  std::vector<std::vector<std::string>> data;

  std::vector<long int> file_offset;
  std::vector<int> data_size;
  std::vector<int> ptrs;
  std::vector<bool> isempty;
  int i = -1;
  int num_active_files = 0;
  const unsigned int num_parts_;

  Context(unsigned int num_parts) : num_parts_(num_parts)
  {
    // set the file offsets to 0
    for (auto i=0; i < this->num_parts_; i++)
      this->file_offset.push_back(0);

    // set the isempty for each file to false
    for (auto i=0; i < this->num_parts_; i++)
      this->isempty.push_back(false);

    // set a vector of vectors of data for each file
    for (auto i=0; i < this->num_parts_; i++)
      this->data.push_back(std::vector<std::string>());

    // set the data_size of the buffer for each file to 0
    for (auto i=0; i < this->num_parts_; i++)
      this->data_size.push_back(0);

    // set the pointer to f each buffer to kDataBufferSize
    for (auto i=0; i < this->num_parts_; i++)
      this->ptrs.push_back(kDataBufferSize);
  }

  void print_status()
  {
    std::cout << "Contx status " << std::endl;
    for (auto i=0; i < this->num_parts_; i++)
    {
      this->i = i;
      std::cout << "\t" << this->i << "\t" << this->data[this->i].size() << "\t"
                << this->data_size[this->i] << "\t"  << this->ptrs[this->i] << std::endl;
    }
  }

  void clear()
  {
    data_size.clear();
    ptrs.clear();
    isempty.clear();
  }
};

using QUEUETUPLE = std::tuple<std::string, int, int>;

extern std::vector<std::string> partial_files;
int filling_counter = 0;

inline std::string ltrim(std::string& s)
{
  auto it = find_if_not(s.begin(), s.end(), [](int c) { return isspace(c); });
  s.erase(s.begin(), it);
  return s;
}

inline std::string rtrim(std::string& s)
{
  auto it = find_if_not(s.rbegin(), s.rend(), [](int c) { return isspace(c); }).base();
  s.erase(it, s.end());
  return s;
}

unsigned int split_buffer_to_fields(const std::string& str, char* line, char** fields, char delim)
{
  // copy the string to a buffer to split by tab
  str.copy(line, str.size(), 0);
  line[str.size()]='\0';
  char* c = line;

  unsigned int k = 0;
  fields[k] = c;
  while (*c!='\0')
  {
    if (*c == delim)
    {
      *c='\0';
      fields[++k] = c + 1;
    }
    c++;
  }
  return k + 1;
}


std::set<std::string> get_mitochondrial_gene_names(const std::string& gtf_filename)
{
  char field_buffer[1000];
  char* fields[20];

  char attrib_buffer[1000];
  char* attribs[20];

  char keyval_buffer[1000];
  char* keyvals[20];

  std::set<std::string> mitochondrial_gene_ids;
  std::ifstream input_file(gtf_filename);
  if (!input_file)
    crash("ERROR failed to open the GTF file " + gtf_filename);

  for (std::string line; std::getline(input_file, line);)
  {
    // skip comment lines
    if (std::regex_search(line, std::regex("^#")))
      continue;

    int num_fields = split_buffer_to_fields(line, field_buffer, fields, '\t');
    // must have at least 8 fields
    assert(num_fields >= 8);

    // skip the line unless it is a gene
    if (std::string(fields[2]).compare(std::string("gene"))!=0)
      continue;

    // split the attributes field separated by ";"
    int num_attribs = split_buffer_to_fields(std::string(fields[8]), attrib_buffer, attribs, ';');

    std::string gene_name("");
    std::string gene_id("");
    // now examine each of the attribute name-value pairs
    for (int k=0; k < num_attribs; k++)
    {
      // splied each attribute name-value pair by the space
      if (std::string(attribs[k]).size() == 0)
        continue;
      std::string attrib = std::move(std::string(attribs[k]));
      attrib = ltrim(attrib);

      int n = split_buffer_to_fields(attrib, keyval_buffer, keyvals, ' ');
      if (n!=2)
        throw std::runtime_error("Expect 2 field but found " + std::to_string(n) + " fields");

      // the second element in the pair is the value string
      std::string value = std::string(keyvals[1]);
      // remove the " (quotes) from the beginning and end of the value string
      value.erase(std::remove_if(value.begin(), value.end(), [](unsigned char c)
      {
        return c=='\"';
      }), value.end());

      if (std::string(keyvals[0]).compare("gene_id")==0)
        gene_id = value;
      if (std::string(keyvals[0]).compare("gene_name")==0)
        gene_name = value;
    }
    if (gene_name.compare("")==0)
      throw ("Malformed GTF file detected. Record is of type gene but does not have a gene_name in line" + line);

    if (std::regex_search(gene_name, std::regex("^mt-", std::regex_constants::icase)))
      mitochondrial_gene_ids.insert(gene_id);
  }
  std::cout << "Number of mitochondrial genes found " << mitochondrial_gene_ids.size() << std::endl;
  return mitochondrial_gene_ids;
}


/*
 * @brief fills the buffer for the files
 *
 * @param contx is the context of the file
*/
void fill_buffer(Context& contx)
{
  contx.data[contx.i].clear();
  int k = 0;

  std::ifstream input_file(partial_files[contx.i]);
  if (!input_file)
    crash("ERROR failed to open the file " + partial_files[contx.i]);

  input_file->seekg(contx.file_offset[contx.i]);

  // the order of the loop condition is iportant first make sure if you can accomodate then try to read,
  // otherwise it might create a read but never processed
  for (std::string line; k < kDataBufferSize && std::getline(input_file, line); k++)
  {
    contx.data[contx.i].push_back(line);
    filling_counter += 1;
  }
  assert(contx.data[contx.i].size() <= kDataBufferSize);

  contx.file_offset[contx.i] = input_file.tellg();

  contx.data_size[contx.i] = contx.data[contx.i].size();

  if (contx.data_size[contx.i] != 0)
  {
    contx.ptrs[contx.i] = 0;
    contx.isempty[contx.i] = false;
  }
  else
  {
    contx.ptrs[contx.i] = kDataBufferSize;
    contx.isempty[contx.i] = true;
  }

#ifdef DEBUG
  std::cout << "-->" << std::endl;
  for (int m = 0; m < contx.num_parts_; m++)
    std::cout << "\t" << m << " : " << contx.data_size[m] << " : " << contx.ptrs[m] << std::endl;
#endif

}

void mergeSortedPartialFiles(const InputOptionsTagsort& options)
{
  const std::string& sorted_output_file = options.sorted_output_file;
  const std::string& metric_type  = options.metric_type;
  const std::string& metric_output_file = options.metric_output_file;

  std::set<std::string> mitochondrial_genes;
  if (!options.gtf_file.empty())
    mitochondrial_genes = get_mitochondrial_gene_names(options.gtf_file);

  // input the buffer size and partial files
  Context contx(partial_files.size());
  auto cmp = [](const QUEUETUPLE &a, const  QUEUETUPLE &b)
  {
    return std::get<0>(a) > std::get<0>(b);
  };
  std::priority_queue<QUEUETUPLE, std::vector<QUEUETUPLE>,  decltype(cmp) > heap(cmp);

  for (auto i=0; i < contx.num_parts_; i++)
  {
    contx.i = i;
    fill_buffer(contx);
  }

  std::regex rgx("\t");
  std::sregex_token_iterator end;

  // create the heap from the first batch loaded data
  contx.num_active_files = 0;
  for (auto i=0; i< contx.num_parts_; i++)
  {
    contx.i = i;
    if (contx.ptrs[i] != kDataBufferSize)
    {
      std::sregex_token_iterator iter(contx.data[i][contx.ptrs[i]].begin(),
                                      contx.data[i][contx.ptrs[i]].end(), rgx, -1);
      std::stringstream comp_tag;
      for (auto k = 0; k < 3 && iter != end; ++iter, k++)
      {
        if (k>0)
          comp_tag << "\t";
        comp_tag << *iter;
      }
      heap.push(QUEUETUPLE(comp_tag.str(), i, contx.ptrs[i]));
      contx.ptrs[i]++;
      contx.num_active_files += 1;
    }
  }

  //  now merge by pop an push
  ofstream fout;
  if (options.compute_metric) // TODO i think this is a mistake, and should actually be options.output_sorted_info
    fout.open(sorted_output_file.c_str());

  // pop and push from the heap
  int num_alignments = 0;
  int i, j;

  Metrics* metric_gatherer;
  MetricType metric_type_enum;
  if (metric_type.compare("cell")==0)
  {
    metric_gatherer = new CellMetrics;
    metric_type_enum = MetricType::Cell;
  }
  else if (metric_type.compare("gene")==0)
  {
    metric_gatherer = new GeneMetrics;
    metric_type_enum = MetricType::Gene;
  }
  else
    crash("Expected metric_type 'cell' or 'gene', got: " + metric_type);

  metric_gatherer->clear();

  ofstream fmetric_out;
  if (options.compute_metric)
  {
    fmetric_out.open(metric_output_file.c_str());
    fmetric_out << metric_gatherer->getHeader() << std::endl;
  }

  // TODO just write directly to fout
  std::stringstream str(stringstream::out | stringstream::binary);
  std::string prev_comp_tag = "";
  while (!heap.empty())
  {
    // read the top
    QUEUETUPLE qtuple = heap.top();
    std::string curr_comp_tag = get<0>(qtuple);
    assert(prev_comp_tag.compare(curr_comp_tag) <= 0);

#ifdef DEBUG
    contx.print_status();
    if (prev_comp_tag.compare(curr_comp_tag) <= 0)
      std::cout << "Expected " << prev_comp_tag << "\n\t\t" << curr_comp_tag << std::endl;
    else
      crash("Anomaly " + prev_comp_tag + "\n\t\t" + curr_comp_tag);
#endif
    i = get<1>(qtuple);  //buffer no
    j = get<2>(qtuple);  //the pointer into the ith buffer array

    heap.pop();

    // start writing in chunks from the stream buffer
    if (num_alignments%kDataBufferSize==0)
    {
      if (options.output_sorted_info)
      {
        fout.write(str.str().c_str(), str.str().length());
        str.clear();
        str.str("");
      }
    }

    // load into stream buffer
    std::string field  = contx.data[i][j];
    if (options.output_sorted_info)
      str << field << std::endl;

    if (options.compute_metric)
      metric_gatherer->parse_line(field, fmetric_out, mitochondrial_genes, metric_type_enum);
    num_alignments += 1;

    // if ismpty is true means the file has been fully read
    if (!contx.isempty[i] && contx.ptrs[i] == contx.data_size[i])
    {
      contx.i = i;
      fill_buffer(contx);
    }

    // make sure it is not empty
    if (contx.data_size[i] > 0)
    {
      std::sregex_token_iterator iter(contx.data[i][contx.ptrs[i]].begin(),
                                      contx.data[i][contx.ptrs[i]].end(), rgx, -1);
      std::stringstream comp_tag;
      for (auto k = 0; k < 3 && iter!=end; ++iter, k++)
      {
        if (k > 0)
          comp_tag << "\t";
        comp_tag << *iter;
      }
      heap.push(QUEUETUPLE(comp_tag.str(), i, contx.ptrs[i]));
      contx.ptrs[i]++;
    }
    else     // one more file is fully read
      contx.num_active_files -= 1;

    if (num_alignments % 1000000 == 0)
      std::cout << "num alns read " << num_alignments << std::endl;

    prev_comp_tag = curr_comp_tag;
  }

  // process the final line
  metric_gatherer->finalize(mitochondrial_genes);
  metric_gatherer->output_metrics(fmetric_out);
  metric_gatherer->output_metrics_extra(fmetric_out);
  delete metric_gatherer;

  // close the metric file
  if (options.compute_metric)
    fmetric_out.close();

  // write out the remaining data
  if (options.output_sorted_info)
  {
    fout.write(str.str().c_str(), str.str().length());
    str.str("");
    str.clear();
  }

  // close output files as there is no more to write
  if (options.output_sorted_info)
    fout.close();

  std::cout << "Written "<< num_alignments << " alignments in total" << std::endl;
  contx.clear();
}

/* Flag set by ‘--verbose’. */
int main(int argc, char** argv)
{
  InputOptionsTagsort options = readOptionsTagsort(argc, argv);

  std::cout << "bam input " << options.bam_input << std::endl;
  std::cout << "temp folder " << options.temp_folder << std::endl;
  std::cout << "sorted output file " <<  options.sorted_output_file << std::endl;
  std::cout << "metric output file " <<  options.metric_output_file << std::endl;
  std::cout << "temp folder " << options.alignments_per_thread << std::endl;
  std::cout << "tags:" << std::endl;

  for (auto it = options.tag_order.begin(); it != options.tag_order.end(); it++)
    std::cout << "\t" << it->first << "\t" << it->second << std::endl;

  /* first create a list of sorted, and simplified sorted files */
  create_sorted_file_splits_htslib(options);

  /* now merge the sorted files to create one giant sorted file by using
    a head to compare the values based on the tags used  */
  std::cout << "Merging " <<  partial_files.size() << " sorted files!"<< std::endl;

  mergeSortedPartialFiles(options);

  // we no longer need the partial files
  for (unsigned int i=0; i < partial_files.size(); i++)
    if (remove(partial_files[i].c_str()) != 0)
      std::cerr << string("Error deleting file") << partial_files[i] << std::endl;

  partial_files.clear();
  std::cout << "Aligments " <<  filling_counter << " loaded to buffer " << std::endl;

  return 0;
}

