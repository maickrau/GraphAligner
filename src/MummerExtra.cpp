//For some reason, the version of mummer on bioconda doesn't have these two methods??
//Manually compile them to prevent undefined references in linking.

// https://github.com/mummer4/mummer/blob/master/src/essaMEM/sparseSA.cpp

// ----------------------

#include <mummer/sparseSA.hpp>
#include <mummer/fasta.hpp>

namespace mummer::mummer {

bool sparseSA::save(const std::string &prefix) const {
  //print auxiliary information
  if(!sparseSA_aux::save(prefix + ".aux"))
    return false;
  //print sa
  if(!SA.save(prefix + ".sa"))
    return false;
  if(!LCP.save(prefix + ".lcp"))
    return false;
  if(hasSufLink && !ISA.save(prefix + ".isa")) //print ISA if nec
    return false;
  if(hasChild){ //print child if nec
    const std::string child = prefix + ".child";
    std::ofstream child_s (child.c_str(), std::ios::binary);
    unsigned int sizeCHILD = CHILD.size();
    child_s.write((const char*)&sizeCHILD,sizeof(sizeCHILD));
    child_s.write((const char*)&CHILD[0],sizeCHILD*sizeof(int));
    if(!child_s.good()) return false;
  }
  if(hasKmer){ //print kmer if nec
    const std::string kmer = prefix + ".kmer";
    std::ofstream kmer_s (kmer.c_str(), std::ios::binary);
    unsigned int sizeKMR = KMR.size();
    kmer_s.write((const char*)&sizeKMR,sizeof(sizeKMR));
    kmer_s.write((const char*)&KMR[0],sizeKMR*sizeof(saTuple_t));
    if(!kmer_s.good()) return false;
  }
  return true;
}

bool sparseSA::load(const std::string &prefix) {
  // Load auxiliary infomation
  if(!sparseSA_aux::load(prefix + ".aux"))
    return false;
  //read sa
  if(!SA.load(prefix + ".sa"))
    return false;
  //read LCP
  LCP.sa = &SA;
  if(!LCP.load(prefix + ".lcp"))
    return false;
  if(hasSufLink && !ISA.load(prefix + ".isa")) //read ISA if nec
    return false;
  if(hasChild){ //read child if nec
    const std::string child = prefix + ".child";
    std::ifstream     child_s (child.c_str(), std::ios::binary);
    unsigned int      sizeCHILD;
    child_s.read((char*)&sizeCHILD,sizeof(sizeCHILD));
    CHILD.resize(sizeCHILD);
    child_s.read((char*)&CHILD[0],sizeCHILD*sizeof(int));
    if(!child_s.good()) return false;
  }
  if(hasKmer){ //read kmer table if nec
    const std::string kmer = prefix + ".kmer";
    std::ifstream     kmer_s (kmer.c_str(), std::ios::binary);
    unsigned int      sizeKMR;
    kmer_s.read((char*)&sizeKMR,sizeof(sizeKMR));
    KMR.resize(sizeKMR);
    kMerTableSize = sizeKMR;
    kmer_s.read((char*)&KMR[0],sizeKMR*sizeof(saTuple_t));
    if(!kmer_s.good()) return false;
  }

  return true;
}

}