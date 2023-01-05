#include "minimizers.h"
#include "big_minimizers.h"
#include <chrono>
#include <fstream>

/*
 * params: k, L, method, seqFile (file to write or read from), len (sequence length), repeats, setFile
 */

 void printUsage(std::string prog){
   std::cout << "Usage: " << prog << " -k <k> -L <L> -r/--reps -l/--seqlen -m/--method -seq/--seqfile -o/--outfile [-set/--setfile]" << std::endl;
   std::cout << "\treps:    # of randomized repeats" <<std::endl;
   std::cout << "\tseqlen:  Length of random sequence to generate." <<std::endl;
   std::cout << "\t         If input sequence is read from file and length is non-zero, this length will be sampled from the input." <<std::endl;
   std::cout << "\tmethod:  Order to use. Options: 'random', 'docks', 'pasha', 'uhs', 'decycling', 'double'" <<std::endl;
   std::cout << "\tseqfile: The sequence is read from this file, or if it doesn't exist, the random sequence will be saved to this file" <<std::endl;
   std::cout << "\toutfile: File to append output statistics to." << std::endl;
   std::cout << "\tsetfile: The file with the set to use. Required when method is 'set' " <<std::endl;
 }



int main (int argc, char* argv[]) {

  unsigned k, L, seqlen,reps;
  std::string method, seqfile, setfile, outfile;

  for (int i = 1; i < argc; ++i){
        std::string param = argv[i];
        if(param == "-k"){k = atoi(argv[++i]);}
        else if(param == "-L"){L = atoi(argv[++i]);}
        else if(param == "-r" || param == "--reps"){reps = atoi(argv[++i]);}
        else if(param == "-l" || param == "--seqlen"){seqlen = atoi(argv[++i]);}
        else if(param == "-m" || param == "--method"){method = argv[++i];}
        else if(param == "-seq" || param == "--seqfile"){seqfile = argv[++i];}
        else if(param == "-set" || param == "--setfile"){setfile = argv[++i];}
        else if(param == "-o" || param == "--outfile"){outfile = argv[++i];}
        else{ printUsage(argv[0]);
              exit(0);
            }
  }



  std::vector<double> densities;
  std::vector<double> runtimes;
  std::vector<double> loads;

  //get the sequence
  for(unsigned i = 0; i < reps; i++){

      char *seq;

      std::string rep_file(seqfile);
      rep_file += "_" + std::to_string(i) + ".txt";
      //if the sequence file exists, read it
      if (FILE *file = fopen(seqfile.c_str(), "r")) {
          fclose(file);
          //fixed sequence
          if(seqlen == 0){
            seq = readSeq(seqfile.c_str());
            seqlen = readSeqLen(seqfile.c_str());
          }
          //sample from fixed sequence
          else{
            unsigned seqfilelen = readSeqLen(seqfile.c_str());
            if(seqfilelen < seqlen){
              throw std::runtime_error("seqlen longer than saved sequence");
            }
            if(seqfilelen == seqlen){
              seq = readSeq(seqfile.c_str());
            }
            else{
              seq = sampleSeq(seqfile.c_str(),seqlen,seqfilelen);
            }
          }
      }
      // random sequence
      else{
        seq = generateRandomSeq(seqlen,rep_file.c_str(),(i+42)*42);
      }

      // Instantiate the desired scheme
      BigMinimizers* bigmin;
      Minimizers* min;
      if(method == "random"){
        if(k>63)
          bigmin = new BigMinimizers(seq,seqlen,k,L);
        else
          min = new Minimizers(seq,seqlen,k,L);
      }
      else if(method == "set" || method == "uhs" || method == "docks" || method == "pasha" || method == "decycling_set" || !setfile.empty()){
        // if(k>63)
        //   bigmin = new SetBasedMinimizers(seq,seqlen,k,L,setfile);
        // else
          min = new SetBasedMinimizers(seq,seqlen,k,L,setfile);
      }
      else if(method =="decycling"){
        if(k>63)
          bigmin = new BigDecyclingMinimizers(seq,seqlen,k,L);
        else
          min = new DecyclingMinimizers(seq,seqlen,k,L);
      }
      else if(method=="double"){
        if(k>63)
          bigmin = new BigDoubleDecyclingMinimizers(seq,seqlen,k,L);
        else
          min = new DoubleDecyclingMinimizers(seq,seqlen,k,L);
      }

      auto start = std::chrono::high_resolution_clock::now();

      std::cout << "Getting " << method << " minimizers. Iteration " << i << std::endl;

      if(k>63)
        bigmin->getMinimizers();
      else
        min->getMinimizers();

      auto end = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
      runtimes.push_back(duration.count()/1000.0);

      if(k>63){
        densities.push_back(bigmin->getDensity());
        loads.push_back(double(bigmin->getMaxBinLoad()));
      }
      else{
        densities.push_back(min->getDensity());
        loads.push_back(double(min->getMaxBinLoad()));
      }
      if(k>63)
        delete bigmin;
      else
        delete min;
      delete seq;
  }

  double avg_density = computeAverage(densities);
  double density_std = computeStd(densities);
  double avg_runtime = computeAverage(runtimes);
  double runtime_std = computeStd(runtimes);
  double avg_load = computeAverage(loads);
  double load_std = computeStd(loads);

  std::ofstream out;

  out.open(outfile.c_str(), std::ios_base::app);

  out << k << "\t" << L << "\t" << method << "\t" << avg_density << "\t" << density_std << "\t"
            << avg_runtime << "\t" << runtime_std << "\t" << avg_load << "\t" << load_std << std::endl;

  return 0;

}
