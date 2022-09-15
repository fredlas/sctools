/**
 *  @file   utilities.cpp
 *  @brief  Utility functions for file processing
 *  @author Kishori Konwar
 *  @date   2021-08-11
 ***********************************************/

#include "utilities.h"

#include <fstream>
#include <iostream>

/** @copydoc readWhiteList */
WhiteListData readWhiteList(std::string const& white_list_file)
{
  const char ATCG[] = {'A', 'C', 'G', 'T', 'N'};

  std::ifstream file(white_list_file);
  if (!file.is_open())
    crash("Couldn't open whitelist file " + white_list_file);

  WhiteListData white_list_data;
  int k = 0;
  // read data from file object and put it into string.
  for (std::string tp; getline(file, tp); )
  {
    white_list_data.barcodes.push_back(tp);

    for (unsigned int i=0; i < tp.size(); i++)
    {
      for (int j=0; j < 5; j++)
      {
        char c = tp[i];
        tp[i] = ATCG[j];
        // If the mutation we're writing is already present, we just overwrite
        // what was there with the current.
        // This is done to have the same values for corrected barcodes
        // as in the python implementation.
        white_list_data.mutations[tp] = k;
        tp[i] = c;
      }
    }

    // -1 suggests it is already a whitelisted barcode
    // This is used, instead of the actual index, because when
    // the barcode is seen with -1 then no correction is necessary.
    // Avoids lots of map lookups, as most barcodes are not erroneous.
    white_list_data.mutations[tp] = -1;
    k++;
  }

  return white_list_data;
}

/** @copydoc crashWithPerror */
void crashWithPerror(std::string msg)
{
  perror(msg.c_str());
  exit(1);
}

void crash(std::string msg)
{
  std::cout << msg << std::endl;
  std::cerr << msg << std::endl;
  exit(1);
}


Token StringBank::getToken(std::string s)
{
  if (auto it = indices_.find(s) ; it != indices_.end())
    return data_[it->second];
  data_.push_back(s);
  indices_[s] = data_.size() - 1;
  return data_.size() - 1;
}
std::string const& StringBank::getString(Token t) const
{
  return data_.find(t)->second;
}

TagTripletManager::TagTripletManager(std::unordered_map<std::string, unsigned int> options_tag_order)
{
  assert(options.tag_order.size() == 3);
  // the order of the three tags are define by the order of the supplied input arguments
  // tag.order [tag_name] -> order map
  if (options.tag_order[options.barcode_tag] == 0 &&
      options.tag_order[options.gene_tag] == 1 &&
      options.tag_order[options.umi_tag] == 2)
  {
    tag_order_ = TagOrder::BGU;
    sb0_ = &sb_b_; sb1_ = &sb_g_; sb2_ = &sb_u_;
  }
  else if (options.tag_order[options.umi_tag] == 0 &&
            options.tag_order[options.barcode_tag] == 1 &&
            options.tag_order[options.gene_tag] == 2)
  {
    tag_order_ = TagOrder::UBG;
    sb0_ = &sb_u_; sb1_ = &sb_b_; sb2_ = &sb_g_;
  }
  else if (options.tag_order[options.umi_tag] == 0 &&
            options.tag_order[options.gene_tag] == 1 &&
            options.tag_order[options.barcode_tag] == 2)
  {
    tag_order_ = TagOrder::UGB;
    sb0_ = &sb_u_; sb1_ = &sb_g_; sb2_ = &sb_b_;
  }
  else if (options.tag_order[options.gene_tag] == 0 &&
            options.tag_order[options.umi_tag] == 1 &&
            options.tag_order[options.barcode_tag] == 2)
  {
    tag_order_ = TagOrder::GUB;
    sb0_ = &sb_g_; sb1_ = &sb_u_; sb2_ = &sb_b_;
  }
  else if (options.tag_order[options.gene_tag] == 0 &&
            options.tag_order[options.barcode_tag] == 1 &&
            options.tag_order[options.umi_tag] == 2)
  {
    tag_order_ = TagOrder::GBU;
    sb0_ = &sb_g_; sb1_ = &sb_b_; sb2_ = &sb_u_;
  }
  else
  {
    tag_order_ = TagOrder::BUG;
    sb0_ = &sb_b_; sb1_ = &sb_u_; sb2_ = &sb_g_;
  }
}

TRIPLET TagTripletManager::makeTriplet(std::string barcode, std::string umi, std::string gene_id)
{
  switch (tag_order_)
  {
    case TagOrder::BUG: return TRIPLET(sb_b_.getToken(barcode), sb_u_.getToken(umi), sb_g_.getToken(gene_id));
    case TagOrder::BGU: return TRIPLET(sb_b_.getToken(barcode), sb_g_.getToken(gene_id), sb_u_.getToken(umi));
    case TagOrder::UBG: return TRIPLET(sb_u_.getToken(umi), sb_b_.getToken(barcode), sb_g_.getToken(gene_id));
    case TagOrder::UGB: return TRIPLET(sb_u_.getToken(umi), sb_g_.getToken(gene_id), sb_b_.getToken(barcode));
    case TagOrder::GUB: return TRIPLET(sb_g_.getToken(gene_id), sb_u_.getToken(umi), sb_b_.getToken(barcode));
    case TagOrder::GBU: return TRIPLET(sb_g_.getToken(gene_id), sb_b_.getToken(barcode), sb_u_.getToken(umi));
    default: crash("no such TagOrder"); return TRIPLET(1,2,3);
  }
}

std::string const& TagTripletManager::getString0(TRIPLET triplet)
{
  return sb0_->getString(std::get<0>(triplet));
}
std::string const& TagTripletManager::getString1(TRIPLET triplet)
{
  return sb1_->getString(std::get<1>(triplet));
}
std::string const& TagTripletManager::getString2(TRIPLET triplet)
{
  return sb2_->getString(std::get<2>(triplet));
}
