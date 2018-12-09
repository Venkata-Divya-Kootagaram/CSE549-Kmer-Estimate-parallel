//#include <google/sparse_hash_map>
//#include <google/dense_hash_map>
//#include "MurmurHash3.cpp"
#include <iostream>
#include <ctime>
#include <climits>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <time.h>
#include "metrohash64.cpp"
#include <stdint.h>
#include <unordered_map>
#include <iomanip>
#include <cmath>
#include <stdlib.h>
#include <cassert>
#include <string.h>
#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <math.h>
#include <sys/time.h>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <list>
#include <stack>
#include <limits.h>
#include <map>
#include <bitset>
#include <ctime>
#include <queue>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <cstring>
#include <iostream>
#include <random>
#include <cinttypes>
#include "ntHashIterator.hpp"

#include <pthread.h>
#include <semaphore.h>
#include <unistd.h>

#define SPP_MIX_HASH 1
#include "sparsepp/spp.h"

using spp::sparse_hash_map;

#define SIZE 10
#define PRODUCER_LOOPS 2
#define CONSUMER_LOOPS 2

using namespace std;
//KSEQ_INIT(gzFile, gzread)
KSEQ_INIT(int, read)
std::map<char, char> mapp = {{'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'}, {'N', 'N'}};
    
kseq_t *seq, *seq2;
int n;
int NUMB_THREADS = 4;

//ntHashIterator itr1;
ntHashIterator *buffer;


void *producer2(void *thread_n) {
    int l;
    int thread_numb = *(int *)thread_n;
    l =  kseq_read(seq);
    //cout<< "thread end "<<thread_numb << seq->seq.s << "\n";
    //cout << "thread : " << thread_numb << ", seq : " <<seq->seq.s<<endl;    
    ntHashIterator itr2(seq->seq.s, 1, n);
    //cout <<"size of bffer : "<< sizeof(buffer);
    buffer[thread_numb] = itr2;
    pthread_exit(0);
}

static const int MultiplyDeBruijnBitPosition[32] =
{
  0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
};
unsigned trailing_zeros(unsigned n) {
    return n ? __builtin_ctz(n) : -1;
}

static const char basemap[255] =
    {   
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*   0 -   9 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  10 -  19 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  20 -  29 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  30 -  39 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  40 -  49 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /*  50 -  59 */
        '\0', '\0', '\0', '\0', '\0',  'T', '\0',  'G', '\0', '\0', /*  60 -  69 */
        '\0',  'C', '\0', '\0', '\0', '\0', '\0', '\0',  'N', '\0', /*  70 -  79 */
        '\0', '\0', '\0', '\0',  'A',  'A', '\0', '\0', '\0', '\0', /*  80 -  89 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0',  't', '\0',  'g', /*  90 -  99 */
        '\0', '\0', '\0',  'c', '\0', '\0', '\0', '\0', '\0', '\0', /* 100 - 109 */
        '\0', '\0', '\0', '\0', '\0', '\0',  'a',  'a', '\0', '\0', /* 110 - 119 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 120 - 129 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 130 - 139 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 140 - 149 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 150 - 159 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 160 - 169 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 170 - 179 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 180 - 189 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 190 - 199 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 200 - 209 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 210 - 219 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 220 - 229 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 230 - 239 */
        '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', /* 240 - 249 */
        '\0', '\0', '\0', '\0', '\0'                                /* 250 - 254 */
    };


unsigned trailing_zeros(uint64_t n) {
    return n ? __builtin_ctzll(n) : -1;
}

void printHelp()
{
    
    cout << "KmerEst [options] -f <fasta/fastq> -k <k-mer length>  -s <sample size> -t <no of threads> -o <output file>"    << endl
    << "  -h               help"                                   << endl
    << "  -f <file>       Input sequence file "                << endl
    << "  -k <k-mer size >        kmer size (default 31) "        << endl
    << "  -s <sample size>        sample size (default 25m)"        << endl
    << "  -t <no of threads>        no of threads (default 4)"        << endl
     << "  -c coverage>       coverage (default 64)"        << endl
    << "  -o         	  Prefix of the Output file " << endl;
    
    exit(0);
}
 
int main(int argc, char** argv)
{
    
    if(argc == 1){
		cout << argv[0] << " -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -c <coverage> -t <no. of threads> -o <out.txt>" << endl;
		exit(0);
    }
    n = 31;
    int s = 25000000;
    int cov = 64;
    string f = "", outf = "";
    for (int c = 1; c < argc; c++){
        if (!strcmp(argv[c], "-h"))       { printHelp(); }
        else if (!strcmp(argv[c], "-k"))     { n = atoi(argv[c+1]); c++; }
        else if (!strcmp(argv[c], "-f"))    { f = argv[c+1]; c++; }
        else if (!strcmp(argv[c], "-s"))    { s = atoi(argv[c+1]); c++; }
        else if (!strcmp(argv[c], "-c"))    { cov = atoi(argv[c+1]); c++; }
        else if (!strcmp(argv[c], "-t"))    { NUMB_THREADS = atoi(argv[c+1]); c++; }
        else if (!strcmp(argv[c], "-o")) { outf = argv[c+1]; c++; }
    }
        
   if (f.empty()  || outf.empty()){
    	printHelp();
    }
    //int n = atoi(argv[2]);
    FILE *fp;
    int l;
    //fp = fopen(argv[1], "r");
    fp = fopen(f.c_str(), "r");
    if( fp == Z_NULL){
		cout<<"File: "<< f  << " does not exist" <<endl;
		exit(1);
    }
    seq = kseq_init(fileno(fp));
    seq2 = seq;
    int k = s;
    int th = 0;
    uint64_t total = 0, no_kmers = 0;
    int count = 0;
    uint64_t hash=0, fhVal=0, rhVal=0; 
    
    //clock_t begin = clock(); 

    while ((l = kseq_read(seq)) >= 0) {
        ++total;
      }
    fclose(fp);
    //ntHashIterator buff[NUMB_THREADS];
    //buffer = buff;

    typedef sparse_hash_map<uint64_t, uint32_t> SMap;
    vector<SMap> MAP(64);
    cout << "read the Sequences .. " << endl;
    cout<<"Number of threads: "<< NUMB_THREADS<<endl;

    fp = fopen(f.c_str(), "r");
    if( fp == Z_NULL){
    	cout<<"File: "<< f  << " does not exist" <<endl;
    	exit(1);
    }
    seq = kseq_init(fileno(fp));
    seq2 = seq;
    k = s;


    ntHashIterator itr;

    int p = 0;
    int NUMB_THREADS2 = NUMB_THREADS;
    for(p = 0; p <= total/NUMB_THREADS2; p++){
    	if(p == total/NUMB_THREADS2) { //for last iteration of 'p'
        	if(total%NUMB_THREADS2 == 0) //for remainder number of sequences
          		break;
        	else {
          		NUMB_THREADS = total%NUMB_THREADS2;  		
        	}
      	}

        ntHashIterator buff[NUMB_THREADS];
        buffer = buff;

		pthread_t thread[NUMB_THREADS];
		int thread_numb[NUMB_THREADS];

        clock_t begin = clock();

		for (int i = 0; i < NUMB_THREADS; i++) {
			thread_numb[i] = i;
			pthread_create(thread + i, // pthread_t *t
							NULL, // const pthread_attr_t *attr
							producer2, // void *(*start_routine) (void *)
							thread_numb+ i);  // void *arg
      	}
		for (int i = 0; i < NUMB_THREADS; i++)
			pthread_join(thread[i], NULL);

		/* Timing threads
        clock_t end = clock();
        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
        cout<<"Time taken:"<<elapsed_secs<<endl; */
        
		for (int i = 0; i < NUMB_THREADS; i++){
		    while (buffer[i] != buffer[i].end()) {
		        hash = (*buffer[i])[0];
		        ++no_kmers;
		        uint8_t tz = trailing_zeros(hash);
		        if(tz >= th){
		            if(MAP[tz].find(hash) != MAP[tz].end())
		            	MAP[tz][hash] += 1; 
		            else{ // insert if not there 
						MAP[tz].insert(make_pair(hash, 1)); 
						++count;  // insert if not there 
						//cout << "\r" << "count: " << count << flush;// << endl;
						if(count == k){
							int cnt = MAP[th].size();
							count = count - cnt;
							SMap().swap(MAP[th]);
							++th;
						}
		            }
		        }
		  		++buffer[i];
			}
		}
        buffer = NULL;
	} 

    cout << "th: " << th << endl;
    cout << "No. of sequences: " << total << endl;
    FILE *fo = fopen(outf.c_str(), "w");
    uint32_t csize = 0; //MAP.size();
    for(int i=th; i<64; i++) csize += MAP[i].size();
    unsigned long F0 = csize * pow(2, (th));
    cout << "F0: " << F0 << endl;
    fprintf(fo, "F1\t%lu\n", no_kmers);
    fprintf(fo, "F0\t%lu\n", F0);
    cout << endl;
    cout << "total: " << total << endl;
    cout << "no_kmer: " << no_kmers << endl;
    unsigned long *freq = new unsigned long[cov];
    for(int i=1; i<=cov; i++) freq[i] = 0;
    unsigned long tot = 0;
    int xx = 0;
    for(int i=th; i<64; i++){
		for(auto& p: MAP[i]){
			if(p.second <= cov) freq[p.second]++;
		}
    }
    cout << "th: " << th << endl;
    for(int i=1; i<=cov; i++){
    	unsigned long fff = (freq[i]*pow(2, th));
    	fprintf(fo, "f%d\t%lu\n", i, fff); 
    }
    fclose(fo);

    /* Timing results
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    cout<<"Time taken:"<<elapsed_secs<<endl; */
    return 0;
        
}