Compile and Run
------------------------------
Compile:


		g++ -o kmerEst kmerCountEstimateParallel.cpp -std=c++11 -O3 -march=native

Run:

		./KmerEst -f <seq.fa> -k  <kmerLen> -s <minHeap_Size> -c <coverage> -t <number of threads> -o <out.txt>
