#ifndef __OPTIMUS_UTILITES__
#define __OPTIMUS_UTILITES__

/**
 *  @file   utilities.h
 *  @brief  Utility functions for file processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/

#include <string>
#include <vector>
#include <unordered_map>

// structure for correcting the barcodes
struct WhiteListData
{
  // an unordered map from whitelist barcodes and 1-mutations
  // to the index of the correct barcode
  std::unordered_map<std::string, int64_t> mutations;
  // vector of whitelist barcodes
  std::vector<std::string> barcodes;
};

/**
 * @brief Build barcode correction map white list barcodes & mutations
 *
 * @details
 * A barcode is computed by checking if it is either in the white
 * list or 1-mutation away from any white listed barcode. To check
 * whether a barcode is correct or to correct it, if 1-mutation away from
 * a barcode in the white list, we build a
 * a map is created with the barcodes and the 1-mutation. The keys are
 * barcodes or mutation and the values are index of the crrect barcode
 *
 * @param whilte_list_file  white list file from 10x genomics' cellranger
 * @return a stricture containing the barcode/1-mutation barcode to index
 *         of the correct barcode
*/
WhiteListData readWhiteList(std::string const& white_list_file);

/**
 * @brief Print system error and exit
 *
 * @param msg  error string to print
*/
void crashWithPerror(std::string msg);

void crash(std::string msg);

// The barcodes/umis/gene_ids have the same string occuring many times, and
// we expect to be sorting hundreds of millions total. So, we want to avoid
// duplicates to save memory. This class encapsulates that.
// Usage: whenever you read a string, feed it to getToken(), and hold onto the
//        returned token in place of the string. Whenever you want to access
//        actual string data, feel that token into getString().
// (Not thread safe; each thread should have its own StringBank).
// (Can only hold 2^32 different strings).
class StringBank
{
public:
  using Token = uint32_t;
  Token getToken(std::string s);
  std::string const& getString(Token t) const;

private:
  std::vector<std::string> data_;
  std::unordered_map<std::string, Token> indices_;
};

using TRIPLET = std::tuple<StringBank::Token, StringBank::Token, StringBank::Token>;

class TagTripletManager
{
public:
  explicit TagTripletManager(std::unordered_map<std::string, unsigned int> options_tag_order);

  TRIPLET makeTriplet(std::string barcode, std::string umi, std::string gene_id);

  std::string const& getString0(TRIPLET triplet);
  std::string const& getString1(TRIPLET triplet);
  std::string const& getString2(TRIPLET triplet);

private:
  enum class TagOrder { BUG, BGU, UBG, UGB, GUB, GBU };
  TagOrder tag_order_;
  StringBank sb_b_;
  StringBank sb_u_;
  StringBank sb_g_;
  // filled by sb_b_, sb_u_, and sb_g_ in the order specified by tag_order_.
  StringBank* sb0_;
  StringBank* sb1_;
  StringBank* sb2_;
};

#endif
