#ifndef MINIMIZERS_H
#define MINIMIZERS_H
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
#include <random>
#include <cassert>

#include "minimizer_utils.h"


class Minimizers {
  protected:
    unsigned k,L;
    __uint128_t mask;
    __uint128_t rand_hash;
  private:
    char * seq;
    uint64_t seqlen;
    std::vector<uint64_t> selected_inds;
    std::vector<__uint128_t> hashed_kmer_buf;
    unsigned buf_ind;
    unsigned buf_start;

    void init_min(){

        if(k > 63){
          throw std::runtime_error("Minimizers limited to K <= 63 by 128 bit primitive types. Use BigMinimizers instead");
        }

        selected_inds.reserve(2*seqlen/L);

        hashed_kmer_buf.resize(L);
        buf_ind = 0;
        buf_start = 0;

        mask = (__uint128_t{1} << 2*k)-1;

        std::random_device rand_dev;
        std::mt19937 generator(rand_dev());
        std::uniform_int_distribution<uint64_t>  top_distr(0, UINT64_MAX);//(1ULL << (2*k)%64)-1);
        std::uniform_int_distribution<uint64_t>  bot_distr(0, UINT64_MAX);//(1ULL << min(2*k,64)-1);
        rand_hash = top_distr(generator);
        rand_hash = (rand_hash << 64 | bot_distr(generator))&mask;
    }

/**
* NOTE: For now we will assume that all letters are nucleotides (and in caps!)
* TODO: Turn get first window into while that restarts if run into non-nt letter
*	whenever run into non-nt letter get new window starting from the next letter.
**/

    uint64_t getFirstWindow(uint64_t &curr_ind, __uint128_t& kmer, __uint128_t & minval){
        uint64_t minind = curr_ind;
        kmer = kmer_to_int(seq+curr_ind,k);
        __uint128_t hashed_mer = myMinHash(kmer);//kmer^rand_hash;
        minval = hashed_mer;
        uint64_t nextchar_ind = curr_ind+k;
        curr_ind++;

        for(int i = 0; i < L-k; i++, curr_ind++, nextchar_ind++){
            kmer = (kmer << 2 | seq_nt4_table[seq[nextchar_ind]])&mask;
            hashed_mer = myMinHash(kmer);//^rand_hash;

            hashed_kmer_buf[buf_ind] = hashed_mer;
            ++ buf_ind;

            if(hashed_mer < minval){
                minval = hashed_mer;
                minind = curr_ind;
                buf_ind = 0;
            }
        }
        return minind;
    }

    virtual __uint128_t myMinHash(__uint128_t kmer){
      return kmer ^ rand_hash;//hash64(kmer,mask);//kmer ^ rand_hash;
    }

    uint64_t getNewWindowFromBuf(__uint128_t & minval){
      //TODO: DOES THIS REALLY SPEED THINGS UP?
      // *OR* SHOULD WE JUST RECOMPUTE THE HASH VALUES FROM THE LAST MIN
       assert(buf_ind - buf_start == L-k+1);
       minval = hashed_kmer_buf[buf_start];
       uint64_t minind = buf_start;
       for(int i = buf_start+1; i < buf_ind; i++){
          if(hashed_kmer_buf[i] < minval){
              minval = hashed_kmer_buf[i];
              minind = i;
          }
       }
       // move the buffer
       for (int i = 0, j=minind+1; j < buf_ind; i++,j++){
         hashed_kmer_buf[i] = hashed_kmer_buf[j];
       }
       buf_start = 0;
       buf_ind = buf_ind - minind-1;
       return minind;
   }


  public:
    Minimizers(char * seq, uint64_t seqlen, unsigned k,unsigned L) : seq(seq), seqlen(seqlen), k(k), L(L){

        init_min();
    }


    void getMinimizers(){
        uint64_t currind = 0;
        uint64_t minind = 0;
        __uint128_t minval;

        // get the minimum of the first window
        __uint128_t kmer;// = kmer_to_int(seq+currind,k);
        uint64_t last_minind = getFirstWindow(currind,kmer,minval);
        __uint128_t hashed_mer;
        selected_inds.push_back(last_minind);

        //  now get for all subsequent windows
        uint64_t nextchar_ind = currind+k;
        for(; currind < seqlen-k+1; currind++, nextchar_ind++){
            kmer = (kmer << 2 | seq_nt4_table[seq[nextchar_ind]])&mask;
            hashed_mer = myMinHash(kmer);//^rand_hash;

            hashed_kmer_buf[buf_ind] = hashed_mer;
            ++ buf_ind;

            // if previous min dropped out of the window, find new min
            if(currind>last_minind+L-k){
                last_minind = last_minind+1+getNewWindowFromBuf(minval);
                selected_inds.push_back(last_minind);
            }

            else if(hashed_mer < minval){
                minval = hashed_mer;
                minind = currind;
                last_minind = currind;
                selected_inds.push_back(minind);
                buf_ind = 0;
                buf_start = 0;
            }
       }
    }

    std::vector<uint64_t> returnMinimizers(){
       return selected_inds;
    }

    double getDensity(){
      return selected_inds.size()/(double(seqlen-k+1));
    }

    std::unordered_map<std::string, uint64_t> getBinLoads(){

        std::unordered_map<std::string, uint64_t> bin_counts;
        bin_counts.reserve(seqlen);
        std::unordered_set<std::string> seen_windows;
        seen_windows.reserve(seqlen);


        uint64_t last_processed_window = 0;
        uint64_t next_window_ind; //start of the first window that overlaps the next minimum
        __uint128_t curr_min = myMinHash(kmer_to_int(seq+selected_inds[0],k));
        std::string currmin_str(seq+selected_inds[0],k);
        __uint128_t next_min;
        for(int i = 0; i<selected_inds.size()-1;i++){
          next_min = myMinHash(kmer_to_int(seq+selected_inds[i+1],k));
          if(curr_min > next_min){
            next_window_ind = selected_inds[i+1]+k-L;
          }
          else{
            next_window_ind = selected_inds[i]+1;
          }


          for(int w = last_processed_window;w<next_window_ind;w++){
  //std::cout << "w: " <<w << std::endl;
              std::string currwindow(seq+w,L);
              if(seen_windows.insert(currwindow).second){ //is unique
                  ++bin_counts[currmin_str];
              }
          }
          curr_min = next_min;
          currmin_str = std::string(seq+selected_inds[i+1],k);
          last_processed_window = next_window_ind;
        }
        // last minimizer
        for (int w = last_processed_window; w<seqlen-L;w++){
            std::string currwindow(seq+w,L);
            if(seen_windows.insert(currwindow).second){ //is unique
                ++bin_counts[currmin_str];
            }
        }
        // for(auto c: bin_counts){
        //   bin_counts[c.first]=c.second;
        // }
        return bin_counts;
    }

    double getMaxBinLoad(){
      std::unordered_map<std::string, uint64_t> binLoads = getBinLoads();
      uint64_t max_load = 0;
      for(auto& it : binLoads){
          if(it.second>max_load) {max_load = it.second;}
      }
      return max_load;
    }

};


class DecyclingMinimizers: public Minimizers {

  protected:

    __uint128_t nondec_mask;
    __uint128_t twobit_mask;
    double u;
    std::vector<std::vector<double> > weights;
    std::vector<char> kmerArray;

    double computeShiftSum(__uint128_t kmer){
      double shiftsum = 0.0;
      kmer = kmer>>2;
      for(int i = 0; i < k-1; i++) {
        shiftsum += weights[kmer & twobit_mask][i];
        kmer = kmer >> 2;
      }
      return shiftsum;
    }

    void getKmerArray(__uint128_t kmer){
      for(int i = k-1;i>=0;i--){
        kmerArray[i] = kmer & twobit_mask;
        kmer = kmer >> 2;
      }
    }

    bool inDec(__uint128_t kmer){
      //TODO: Is this actually faster, or is it better to just convert to an array and iterate over the array (see original code)
      double sum = 0.0;
      __uint128_t shiftmer = kmer;
      for(int i = 0; i < k - 1; i++){
        sum += weights[shiftmer & twobit_mask][i];
        shiftmer = shiftmer >> 2;
      }
      if(sum < -0.00000001){ return false;} // R node

      else if(sum > 0.00000001){ // L node - check if it's the first
        double shiftsum = computeShiftSum(kmer);
        if(shiftsum > 0.00000001) {return false;} // not the first L node
        else {return true;} // the first L node
      }
      else{ // I node - check if it's an I-cycle and if yes, check if this is the lexicographically smallest rotation
        double shiftsum = computeShiftSum(kmer);
        getKmerArray(kmer);
        if(shiftsum > 0.00000001 || shiftsum < -0.00000001) {return false;} // not all I
        // an I-cycle, check if it's the smallest shift
        int i = 0, j = 1;
        for(; j < k; j++){
          if(kmerArray[j]<kmerArray[i]) {return false;}
          if(kmerArray[j]>kmerArray[i]) {i = 0;}
          else {i++;}
          if(i==0 || i==k) {return true;}
        }
        for(j=0;j<k;j++){
          if(kmerArray[j]<kmerArray[i]) {return false;}
          if(kmerArray[j]>kmerArray[i]) {i = 0;}
          else {i++;}
          if(i==0 || i==k) {return true;}
        }
      }
    }


    __uint128_t myMinHash(__uint128_t kmer){
        __uint128_t hashed_mer = kmer ^ rand_hash;//hash64(kmer,mask);//
        if(!inDec(kmer)){
          hashed_mer = hashed_mer | nondec_mask;
        }
          return hashed_mer;
    }

    void getWeights(){
      for(int j = 0; j < 4; j++){
        std::vector<double> jweights;
        for(int i = k-1; i >= 0; i--){
          jweights.push_back(sin(i*u)*j); //NOTE: weights is in opposite order
        }
        weights.push_back(jweights);
      }
    }

  public:
    DecyclingMinimizers(char * seq, uint64_t seqlen, unsigned k,unsigned L)
                                            : Minimizers(seq, seqlen, k, L) {
      u = 3.1415926535898*2.0 / double(k);
      nondec_mask = __uint128_t{1} << 2*k;
      twobit_mask = __uint128_t{3};
      kmerArray.resize(k);
      getWeights();
    }

};

class DoubleDecyclingMinimizers : public DecyclingMinimizers{
public:
  DoubleDecyclingMinimizers(char * seq, uint64_t seqlen, unsigned k,unsigned L)
                                          : DecyclingMinimizers(seq, seqlen, k, L) {
    nondec_mask = __uint128_t{1} << 2*k+1;
    symdec_mask = __uint128_t{1} << 2*k;
  }
protected:
    __uint128_t symdec_mask;

    __uint128_t myMinHash(__uint128_t kmer){
        __uint128_t hashed_mer = kmer ^ rand_hash;//hash64(kmer,mask);//
        unsigned setMembership = inDoubleDec(kmer);
        if(setMembership == 0){
          hashed_mer = hashed_mer | nondec_mask;
        }
        else if(setMembership == 1){ // in the symmetric set
          hashed_mer = hashed_mer | symdec_mask;
        }
          return hashed_mer;
    }


    unsigned inDoubleDec(__uint128_t kmer){
      //TODO: Is this actually faster, or is it better to just convert to an array and iterate over the array (see original code)
      double sum = 0.0;
      __uint128_t shiftmer = kmer;
      for(int i = 0; i < k - 1; i++){
        sum += weights[shiftmer & twobit_mask][i];
        shiftmer = shiftmer >> 2;
      }
      double shiftsum = computeShiftSum(kmer);
      if(sum > 0.00000001 && shiftsum < 0.00000001) {return 2;} //the first L node
      if(sum < -0.00000001 && shiftsum > -0.00000001){return 1;} //the first R node
      if( (sum < 0.00000001 && sum > -0.00000001) &&
          (shiftsum < 0.00000001 && shiftsum > -0.00000001) ){ //an I-cycle, check if it's the lexicographically smallest rotation

        getKmerArray(kmer);

        int i = 0, j = 1;
        for(; j < k; j++){
          if(kmerArray[j]<kmerArray[i]) {return 0;}
          if(kmerArray[j]>kmerArray[i]) {i = 0;}
          else {i++;}
          if(i==0 || i==k) {return 2;}
        }
        for(j=0;j<k;j++){
          if(kmerArray[j]<kmerArray[i]) {return 0;}
          if(kmerArray[j]>kmerArray[i]) {i = 0;}
          else {i++;}
          if(i==0 || i==k) {return 2;}
        }
      }
      return 0;
    }
};



class SetBasedMinimizers : public Minimizers{
  private:
    std::vector<bool> kmer_set;
    std::string setfile;
    __uint128_t nonset_mask;

    __uint128_t myMinHash(__uint128_t kmer){
        __uint128_t hashed_mer = kmer ^ rand_hash;//hash64(kmer,mask);//
        if(!kmer_set[kmer]){
          hashed_mer = hashed_mer | nonset_mask;
        }
          return hashed_mer;
    }

  public:
    SetBasedMinimizers(char * seq, __uint128_t seqlen, unsigned k,unsigned L, std::string setfile) : Minimizers(seq, seqlen, k, L), setfile(setfile){
        nonset_mask = __uint128_t{1} << 2*k;
        kmer_set = getSetFromIntFile(setfile,k);
    }
};



#endif
