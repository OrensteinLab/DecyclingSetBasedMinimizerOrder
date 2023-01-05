#ifndef BIGMINIMIZERS_H
#define BIGMINIMIZERS_H
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
#include "gmpxx.h"

#include "minimizer_utils.h"


class BigMinimizers {
  protected:
    unsigned k,L;
    mpz_class mask;
    mpz_class rand_hash;
  private:
    char * seq;
    uint64_t seqlen;
    std::vector<uint64_t> selected_inds;
    std::vector<mpz_class> hashed_kmer_buf;
    unsigned buf_ind;
    unsigned buf_start;

    void init_min(){

        // if(k > 63){
        //   throw std::runtime_error("Minimizers limited to K <= 63 by 128 bit primitive types. Use BigMinimizers instead");
        // }
        if(k<64){
          std::cout << "NOTE: This class uses GMP to represent long kmers." <<std::endl;
          std::cout << "For k < 64, regular Minimizers class is 10x faster." << std::endl;
        }

        selected_inds.reserve(2*seqlen/L);

        hashed_kmer_buf.resize(L);
        buf_ind = 0;
        buf_start = 0;

        mask = (mpz_class{1} << 2*k)-1;

        std::random_device rand_dev;
        std::mt19937 generator(rand_dev());
        std::uniform_int_distribution<uint64_t>  distr(0, UINT64_MAX);
        rand_hash = distr(generator);
        for(int i=64;i<2*k;i+=64){
            rand_hash = (rand_hash << 64 | distr(generator))&mask;
        }
    }

/**
* NOTE: For now we will assume that all letters are nucleotides (and in caps!)
* TODO: Turn get first window into while that restarts if run into non-nt letter
*	whenever run into non-nt letter get new window starting from the next letter.
**/

    uint64_t getFirstWindow(uint64_t &curr_ind, mpz_class& kmer, mpz_class& hashed_mer, mpz_class & minval){
        uint64_t minind = curr_ind;
        kmer = kmer_to_int<mpz_class>(seq+curr_ind,k);
        myMinHash(kmer,hashed_mer);//kmer^rand_hash;
        minval = hashed_mer;
        uint64_t nextchar_ind = curr_ind+k;
        curr_ind++;

        for(int i = 0; i < L-k; i++, curr_ind++, nextchar_ind++){
            kmer = (kmer << 2 | seq_nt4_table[seq[nextchar_ind]])&mask;
            myMinHash(kmer,hashed_mer);//^rand_hash;

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

    virtual void myMinHash(mpz_class& kmer, mpz_class& hashed_mer){
      hashed_mer = kmer ^ rand_hash;//hash64(kmer,mask);//kmer ^ rand_hash;
    }

    uint64_t getNewWindowFromBuf(mpz_class & minval){
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
    BigMinimizers(char * seq, uint64_t seqlen, unsigned k,unsigned L) : seq(seq), seqlen(seqlen), k(k), L(L){

        init_min();
    }


    void getMinimizers(){
        uint64_t currind = 0;
        uint64_t minind = 0;
        mpz_class minval;

        // get the minimum of the first window
        mpz_class kmer;// = kmer_to_int(seq+currind,k);
        mpz_class hashed_mer;
        uint64_t last_minind = getFirstWindow(currind,kmer,hashed_mer,minval);
        selected_inds.push_back(last_minind);

        //  now get for all subsequent windows
        uint64_t nextchar_ind = currind+k;
        for(; currind < seqlen-k+1; currind++, nextchar_ind++){
            kmer = (kmer << 2 | seq_nt4_table[seq[nextchar_ind]])&mask;
            myMinHash(kmer,hashed_mer);//^rand_hash;

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
        mpz_class curr_min;
        mpz_class kmer = kmer_to_int<mpz_class>(seq+selected_inds[0],k);
        myMinHash(kmer,curr_min);
        std::string currmin_str(seq+selected_inds[0],k);
        mpz_class next_min;
        for(int i = 0; i<selected_inds.size()-1;i++){
          kmer = kmer_to_int<mpz_class>(seq+selected_inds[i+1],k);
          myMinHash(kmer,next_min);
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


class BigDecyclingMinimizers: public BigMinimizers {

  protected:

    mpz_class nondec_mask;
    mpz_class twobit_mask;
    double u;
    std::vector<std::vector<double> > weights;
    std::vector<char> kmerArray;

    double computeShiftSum(mpz_class kmer){
      double shiftsum = 0.0;
      kmer = kmer>>2;
      for(int i = 0; i < k-1; i++) {
        shiftsum += weights[(kmer.get_ui() & twobit_mask.get_ui())][i];
        kmer = kmer >> 2;
      }
      return shiftsum;
    }

    void getKmerArray(mpz_class kmer){
      for(int i = k-1;i>=0;i--){
        kmerArray[i] = (kmer.get_ui() & twobit_mask.get_ui());
        kmer = kmer >> 2;
      }
    }

    bool inDec(mpz_class kmer){
      //TODO: Is this actually faster, or is it better to just convert to an array and iterate over the array (see original code)
      double sum = 0.0;
      mpz_class shiftmer = kmer;
      for(int i = 0; i < k - 1; i++){
        sum += weights[(shiftmer.get_ui() & twobit_mask.get_ui())][i];
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


    void myMinHash(mpz_class& kmer,mpz_class& hashed_mer){
        hashed_mer = kmer ^ rand_hash;//hash64(kmer,mask);//
        if(!inDec(kmer)){
          hashed_mer = hashed_mer | nondec_mask;
        }
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
    BigDecyclingMinimizers(char * seq, uint64_t seqlen, unsigned k,unsigned L)
                                            : BigMinimizers(seq, seqlen, k, L) {
      u = 3.1415926535898*2.0 / double(k);
      nondec_mask = mpz_class{1} << 2*k;
      twobit_mask = mpz_class{3};
      kmerArray.resize(k);
      getWeights();
    }

};

class BigDoubleDecyclingMinimizers : public BigDecyclingMinimizers{
public:
  BigDoubleDecyclingMinimizers(char * seq, uint64_t seqlen, unsigned k,unsigned L)
                                          : BigDecyclingMinimizers(seq, seqlen, k, L) {
    nondec_mask = mpz_class{1} << 2*k+1;
    symdec_mask = mpz_class{1} << 2*k;
  }
protected:
    mpz_class symdec_mask;

    void myMinHash(mpz_class& kmer, mpz_class& hashed_mer){
        hashed_mer = kmer ^ rand_hash;//hash64(kmer,mask);//
        unsigned setMembership = inDoubleDec(kmer);
        if(setMembership == 0){
          hashed_mer = hashed_mer | nondec_mask;
        }
        else if(setMembership == 1){ // in the symmetric set
          hashed_mer = hashed_mer | symdec_mask;
        }
    }


    unsigned inDoubleDec(mpz_class kmer){
      //TODO: Is this actually faster, or is it better to just convert to an array and iterate over the array (see original code)
      double sum = 0.0;
      mpz_class shiftmer = kmer;
      for(int i = 0; i < k - 1; i++){
        sum += weights[(shiftmer.get_ui() & twobit_mask.get_ui())][i];
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





#endif
