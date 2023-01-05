#ifndef MY_UTILS_H
#define MY_UTILS_H
#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <unordered_set>
#include <fstream>
#include <cmath>
#include <string.h>
#include <vector>
#include <algorithm>
#include <initializer_list>
#include <cstdlib>

const unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

template <typename T=__uint128_t> std::string decimal_string(T n) {
  std::string result;

  do {
    result.push_back("0123456789"[n % 10]);
    n /= 10;
  } while (n);

  std::reverse(result.begin(), result.end());
  return result;
}

template <typename T=__uint128_t> T kmer_to_int(char* kmer, int k){
  T res = 0;
  for(size_t i = 0; i < k; i++) {
    const size_t c = seq_nt4_table[(size_t)*kmer++];
    if(c & 4)
      return 0;
    res = (res << 2) | c;
  }
  return res;
}

template <typename T=__uint128_t> static inline T hash64(T key, T mask)
{
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

template <typename T=__uint128_t> std::vector<bool> getSetFromStringFile(std::string filename, int k){
  std::ifstream infile(filename);
  std::string line;
  std::vector<bool> kmerset(pow(4, k),false);
  while(std::getline(infile, line)){
     T mer = kmer_to_int(&line[0],k);
     kmerset[mer]=true;
  }
  infile.close();
  return kmerset;
}

std::vector<bool> getSetFromIntFile(std::string filename, int k){
  std::ifstream infile(filename);
  std::vector<bool> kmerset(pow(4, k),false);
  uint64_t mer;
  while(infile >> mer){
     kmerset[mer]=true;
  }
  infile.close();
  return kmerset;
}

char* generateRandomSeq(int seqlen, const char* outFile,unsigned seed) {
    srand(seed+time(0));
    char* seq = new char[seqlen];
    static const char alphabet[] =
        "A"
        "C"
        "G"
        "T";

    for (int i = 0; i < seqlen; ++i) {
        seq[i] = alphabet[rand() % (sizeof(alphabet) - 1)];
    }
    std::ofstream seqFile;
    seqFile.open(outFile);
    seqFile << seq << std::endl;
    seqFile.close();

    return seq;
}

std::vector<char> getFileSeqVec(const char* inFile){
	std::ifstream seqFile(inFile);
	if (!seqFile.good()) {
		std::cout << "Error opening input file: " << inFile << std::endl;
		exit(1);
  }
	std::vector<char> seqv;
	char base;
	//seqFile.get(base);
	while (seqFile >> std::noskipws >> base) {
		if(base == 'A' || base == 'C' || base == 'G' || base == 'T') {
			seqv.push_back(base);
		}
		//seqFile.get(base);
	}
	return seqv;
}

unsigned readSeqLen(const char* inFile){
	std::vector<char> seqv = getFileSeqVec(inFile);
	return seqv.size();
}

char* readSeq(const char* inFile){
	std::vector<char> seqv = getFileSeqVec(inFile);
  std::cout << seqv.size() << std::endl;
  uint64_t len = seqv.size();
  char* seq = new char[len];
  for (int i = 0; i < len; i++) seq[i] = seqv[i];
  return seq;
}

char* sampleSeq(const char* inFile, unsigned seqlen, unsigned readLen){
	std::random_device rand_dev;
	std::mt19937 generator(rand_dev());
	std::uniform_int_distribution<unsigned>  distr(0, readLen-seqlen+1);
	unsigned rand_start = distr(generator);

	std::vector<char> seqv = getFileSeqVec(inFile);
	char* seq = new char[seqlen];
  for (int i = 0; i < seqlen; i++) seq[i] = seqv[rand_start+i];
  return seq;
}

double computeAverage(std::vector<double>& values){
	double avg = 0.0;
	for(double v : values){
		avg += v;
	}
	return avg/double(values.size());
}

double computeStd(std::vector<double> values){
 double avg = computeAverage(values);
 double stddev = 0.0;
 for(double v: values){
	 stddev += (v-avg)*(v-avg);
 }
return sqrt(stddev/double(values.size()));
}


#endif
