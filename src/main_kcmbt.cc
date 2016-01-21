

#include <iostream>
#include <chrono>
#include <zlib.h>
#include <stdio.h>
#include <ctime>
#include <fstream>
#include <sys/resource.h>
#include <cassert>
#include <vector>
#include <algorithm>
#include <iterator>

#include "burst_sort_kmer.h"
#include "kseq.h"

using namespace std;



using Clock = std::chrono::high_resolution_clock;
using TimePoint = Clock::time_point;

const int kMaxBuff = 1 << 10;


KSEQ_INIT(gzFile, gzread)

uint64_t total_kmer_count = 0, total_uniq_count = 0, tkc1 = 0, tkc3 = 0, max_uniq = 0, t_u = 0;
TrieNode* root[kTotalLayer * total_tree];

string in_file_name, out_file_name ("out");


void HowToUse() {
	cout << "\n\n\n";

	string h("KCMBT (k-mer Counter based on Multiple Burst Trees) 1.0.0");
	cout << h << endl;
	vector<char> c_arr(h.size(), '=');
	copy(c_arr.begin(), c_arr.end(), ostream_iterator<char>(cout));
	cout << endl;

	cout << "Usage: \n\t" << "kcmbt -k <k-mer length> -i <@file_listing_fastq_files or fastq_file> -o <output_file_name>" << endl;
	cout << "Example: \n\t" << "kcmbt -k 28 -i srr.fastq -o srr" << endl;
	cout << "Parameters: \n";
	cout << "\t-k\t\tk-mer length (10 <= k <= 32, default 28)\n";
	cout << "\t-i\t\tinput file in fastq format (start with @ if the file contains a list of fastq files)\n";
	cout << "\t-o\t\toutput file (a binary file; please run hr_kcmbt to generate human readable output)\n";

	cout << "\n\n\n";
}

void SortKmer(TrieNode* tree_root, uint64_t tree_ind, FILE* out_file) {
	uint64_t kmer_uniq = 0;
	TraverseFinal(tree_root, 0, kmer_uniq, tree_ind);

#ifdef DEBUG
	//bucket_fill += bucket_fill_arr[root->ind];
	//bucket_size += bucket_size_arr[root->ind];
	//std::cout << root->ind << ": " << bucket_fill_arr[root->ind] << "\t" <<  bucket_size_arr[root->ind] << "\t" << (bucket_fill_arr[root->ind] * 1.0 /  bucket_size_arr[root->ind]) << "\n";
#endif

	total_uniq_count += kmer_uniq;
	//cout << tree_ind << ":\t" << kmer_uniq << "\t" << total_uniq_count << endl;

	if (kmer_uniq > max_uniq)
		max_uniq = kmer_uniq;

	pre_arr[tree_ind + 1] = kmer_uniq;
	fwrite(kmer_arr, sizeof(uint64_t), kmer_uniq, out_file);
}

uint64_t ff1 = 0;

void SortInsertKmer(int tree_type, int ind) {
	uint64_t kmer_uniq = 0;
	Traverse(root[tree_type * total_tree + ind], 0, kmer_uniq);
	if (kmer_uniq > max_uniq)
		max_uniq = kmer_uniq;
	t_u += kmer_uniq;
	uint64_t k_len =  2 * (kmer_len - tree_pow);
	uint64_t mask_kmer = ((1ULL << k_len) - 1) << (64 - k_len);
	uint64_t tree_mask = total_tree - 1;

	uint64_t* buffer __attribute__((aligned(64))) = new uint64_t[kMaxBuff * total_tree];
	uint64_t pos_arr[total_tree];
	std::memset(pos_arr, 0, total_tree * sizeof(uint64_t));

	for (int i = 0; i < kmer_uniq; ++i) {
		uint64_t k_c = kmer_arr[i] & k_count_mask[tree_type];
		uint64_t k_p = kmer_arr[i] & k_kmer_mask[tree_type];
		for (int j = 0; j <= tree_type; ++j) {
			uint64_t kmer = ((k_p << 2 * j) & mask_kmer) | k_c;
			int tree_ind = ((ind << 2 * j) | (k_p >> (64 - 2 * j))) & tree_mask;
			buffer[tree_ind * kMaxBuff + (pos_arr[tree_ind]++)] = kmer;
			if (pos_arr[tree_ind] >= kMaxBuff) {
				InsertBatch(root[0 * total_tree + tree_ind], &buffer[tree_ind * kMaxBuff], pos_arr[tree_ind]);
				pos_arr[tree_ind] = 0;
			}
#ifdef TEST1
			if ((tree_type == 3) && (i > 100000) && (ind >= 1) && k_c > 2)
				std::cout << i << ":" << j << "\t" <<  kmer_arr[i] << " " << kmer << " " << tree_ind << " " << k_c << " " << k_p << std::endl;
#endif

		}
#ifdef TEST1
		if ((tree_type == 3) && (i > 100000) && (ind >= 1) && k_c > 2)
			exit(1);
#endif
	}

	for (uint64_t i = 0; i < total_tree; ++i) {
		InsertBatch(root[0 * total_tree + i], &buffer[i * kMaxBuff], pos_arr[i]);
	}

	delete[] buffer;
}

bool Compare(const BigKmer& k1, const BigKmer& k2) {
	return k1.kmer < k2.kmer;
}

void CountBigKmer() {
	cout << "CountBigKmer " << over_c << endl;
	int n = bkmer_arr.size(), j = 0;
	uint64_t* bkmerw_arr = new uint64_t[n];
	sort(bkmer_arr.begin(), bkmer_arr.end(), Compare);
	for (int i = 0; i < n; ++i) {
		bkmerw_arr[j++] = bkmer_arr[i].kmer;
		bkmerw_arr[j] = bkmer_arr[i].count + k_count_mask[0];
		while (i + 1 < n && bkmer_arr[i].kmer == bkmer_arr[i + 1].kmer) {
			bkmerw_arr[j] += bkmer_arr[i + 1].count;
			++i;
		}
		++j;
	}
	string big ("big_" + out_file_name);
	FILE* big_file = fopen(big.c_str(), "wb");
	fwrite(bkmerw_arr, sizeof(uint64_t), j, big_file);
	fclose(big_file);

#ifdef DEBUG
	cout << n << " big kmer count " << j << endl;
#endif
	total_uniq_count += (j >> 1);

	delete[] bkmerw_arr;
}


int main(int argc, char** argv) {
	if (argc % 2 == 0) {
		HowToUse();
		return EXIT_FAILURE;
	}

	for (int i = 1; i < argc; i += 2) {
		if (strcmp(argv[i], "-k") == 0)
			kmer_len = std::atoi(argv[i + 1]); // command line kmer size. default = 28
		else if (strcmp(argv[i], "-i") == 0)
			in_file_name = argv[i + 1];
		else if (strcmp(argv[i], "-o") == 0)
			out_file_name = argv[i + 1];
	}

	if (in_file_name.empty()) {
		cout << "Please specify input file";
		HowToUse();
		return EXIT_FAILURE;
	}


	kTotalMaskBitArr = new uint64_t[kTotalLayer]; // how many bits are left not for kmer
	k_count_mask = new uint64_t[kTotalLayer]; // mask for counter bits right aligned
	k_kmer_mask = new uint64_t[kTotalLayer]; // mask for kmer+x bits left aligned
	for (int i = 0; i < kTotalLayer; ++i) {
		kTotalMaskBitArr[i] = (32 - kmer_len + tree_pow - i) << 1;
		k_count_mask[i] = (1ULL << kTotalMaskBitArr[i]) - 1;
		k_kmer_mask[i] = ~k_count_mask[i];
	}

	mask = total_tree - 1; // tree index mask
	k_mask = (~0ULL) >> ((32 - kmer_len) << 1); // kmer+0 mask right aligned

#ifdef DEBUG
	std::cout << kTotalMaskBitArr[0] << " " << kTotalMaskBitArr[1] << " " << kTotalMaskBitArr[2] << " " << kTotalMaskBitArr[3] << " " <<  k_count_mask[0] << " " <<   k_kmer_mask[0] << " " <<      mask << " " <<  k_mask << std::endl;
#endif

	// mapping letters to int
	char codes[256];
	for(int i = 0; i < 256; ++i)
		codes[i] = -1;
	codes['A'] = codes['a'] = 0;
	codes['C'] = codes['c'] = 1;
	codes['G'] = codes['g'] = 2;
	codes['T'] = codes['t'] = 3;
	// ===

	// time count starts
	TimePoint tp1 = Clock::now();
	std::clock_t cput1 = std::clock();

	// extract input file names
	std::vector<std::string> in_file_arr;
	if (in_file_name[0] == '@') {
		std::ifstream in_file(in_file_name.substr(1, in_file_name.length() - 1));
		std::string file_str;
		while (in_file >> file_str) {
			in_file_arr.push_back(file_str);
			std::cout << file_str << "\n";
		}
		in_file.close();
	}
	else
		in_file_arr.push_back(in_file_name);

	// ===

	temp_mem = new uint64_t[10 * kMaxBucketSize]; // temporary memory for sorting

	uint64_t total_entry = kTotalLayer * total_tree;
	uint64_t* buffer = new uint64_t[total_entry * kMaxBuff];// __attribute__((aligned(64))); // kmer buffers, while full, insert into trees
	int pos_arr[total_entry]; // pos in buffer
	int kmer_count_arr[total_entry]; // count kmers in each tree


	for (uint64_t i = 0; i < kTotalLayer; ++i) {
		for (uint64_t j = 0; j < total_tree; ++j) {
			uint64_t ind = i * total_tree + j;
			root[ind] = new TrieNode(i);
			root[ind]->InitRoot();
			pos_arr[ind] = 0;
			kmer_count_arr[ind] = 0;
		}
	}

#ifdef DEBUG
	std::cout << k_mask << "\n";
#endif

	uint64_t k_mask_arr[kAlphabetSize], r_shift_arr[kAlphabetSize], l_shift_arr[kAlphabetSize];
	for (int i = 0; i < kAlphabetSize; ++i) {
		k_mask_arr[i] = ~0ULL >> ((32 - kmer_len - i) << 1); // mask for kmer+x left aligned
		r_shift_arr[i] = (kmer_len + i - tree_pow) << 1; // right shift for tree kmer+x
		l_shift_arr[i] = (32 - kmer_len - i + tree_pow) << 1; // left shift for tree kmer+x
	}

	uint64_t kmer, kmer_rev, kmer_bit;
	uint64_t len, kmer_count = 0, tree_ind, c_layer, r_pos;
	uint64_t prev_k_less;

	for (int f = 0; f < in_file_arr.size(); ++f) { // read each file and generate kmers
		gzFile fp;
		kseq_t *seq;
		fp = gzopen(in_file_arr[f].c_str(), "r");
		seq = kseq_init(fp);

		while (kseq_read(seq) >= 0) { // read read sequences
			prev_k_less = kPrevInit; // 2 means unassigned
			r_pos = 0;
			c_layer = 0;
			len = 0;
			kmer = 0ULL;
			kmer_rev = 0ULL;

			char *s = seq->seq.s;
			for (int i=0; i<seq->seq.l; ++i, ++s) {
				int c = codes[*s];
				if (c < 0) { //'N'
					if (c_layer > 0) {
						c_layer--;
						bool b = prev_k_less == 1? true : false; // true means kmer is minimal, false means rev kmer is small
						tree_ind = b * (kmer >> r_shift_arr[c_layer]) + !b * (kmer_rev >> r_shift_arr[c_layer]);
						tree_ind &= mask;
						kmer_bit = b * ((kmer << l_shift_arr[c_layer]) | kInitCount) + !b * ((kmer_rev << l_shift_arr[c_layer]) | kInitCount);
						buffer[(c_layer * total_tree + tree_ind) * kMaxBuff + (pos_arr[c_layer * total_tree + tree_ind]++)] = kmer_bit;
						if (pos_arr[c_layer * total_tree + tree_ind] >= kMaxBuff) { // does not affect that much
							InsertBatch(root[(c_layer * total_tree + tree_ind)], &buffer[(c_layer * total_tree + tree_ind) * kMaxBuff], pos_arr[c_layer * total_tree + tree_ind]);
							kmer_count_arr[c_layer * total_tree + tree_ind] += pos_arr[c_layer * total_tree + tree_ind];
							pos_arr[c_layer * total_tree + tree_ind] = 0;
						}
					}

					prev_k_less = kPrevInit;
					r_pos = 0;
					c_layer = 0;
					len = 0;
					kmer = 0ULL;
					kmer_rev = 0ULL;
					continue;
				}
				++len;
				kmer = (kmer << 2) | c;
				kmer_rev = kmer_rev | ((3ULL - c) << r_pos);
				r_pos += 2;
				if (len >= kmer_len) {
					uint64_t kmer0 = kmer & k_mask;
					uint64_t kmer_rev0 = (kmer_rev >> 2 * c_layer) & k_mask;
					bool b = kmer0 < kmer_rev0;
					bool cond1 = (b && prev_k_less == 1) || (!b && prev_k_less == 0); // same kmer as before
					bool cond2 = prev_k_less > 1 || (c_layer < kTotalLayerMinusOne && cond1 ); // check whether continue with this extended kmer
					c_layer += cond2;


					if (!cond2) {
						c_layer += cond1;
						c_layer--;

						tree_ind = cond1 * (prev_k_less * (kmer >> r_shift_arr[c_layer])   + !prev_k_less * (kmer_rev >> r_shift_arr[c_layer]))
								+ !cond1 * (prev_k_less * ((kmer >> r_shift_arr[c_layer] + 2))  + !prev_k_less * (kmer_rev >> r_shift_arr[c_layer]));
						tree_ind &= mask;
						kmer_bit = cond1 * (prev_k_less * ((kmer << l_shift_arr[c_layer]) | kInitCount) + !prev_k_less * ((kmer_rev << l_shift_arr[c_layer]) | kInitCount))
								+ !cond1 * (prev_k_less * (((kmer << (l_shift_arr[c_layer] - 2)) & k_kmer_mask[c_layer]) | kInitCount) + !prev_k_less * ((kmer_rev << l_shift_arr[c_layer]) | kInitCount));
						buffer[(c_layer * total_tree + tree_ind) * kMaxBuff + (pos_arr[c_layer * total_tree + tree_ind]++)] = kmer_bit;
						if (pos_arr[c_layer * total_tree + tree_ind] >= kMaxBuff) { // does not affect that much
							InsertBatch(root[(c_layer * total_tree + tree_ind)], &buffer[(c_layer * total_tree + tree_ind) * kMaxBuff], pos_arr[c_layer * total_tree + tree_ind]);
							kmer_count_arr[c_layer * total_tree + tree_ind] += pos_arr[c_layer * total_tree + tree_ind];
							pos_arr[c_layer * total_tree + tree_ind] = 0;
						}

						kmer = kmer & k_mask;
						kmer_rev = (kmer_rev >> (c_layer + 1) * 2) & k_mask;
						bool cond3 = c_layer == kTotalLayerMinusOne;
						c_layer = 1 - cond3;
						r_pos = (kmer_len - cond3) << 1;
						prev_k_less = b + cond3 * 2;
					}
					else
						prev_k_less = b;

					++kmer_count;
				}
			}

			if (c_layer > 0) {
				c_layer--;
				bool b = prev_k_less == 1? true : false;
				tree_ind = b * (kmer >> r_shift_arr[c_layer]) + !b * (kmer_rev >> r_shift_arr[c_layer]);
				tree_ind &= mask;
				kmer_bit = b * ((kmer << l_shift_arr[c_layer]) | kInitCount) + !b * ((kmer_rev << l_shift_arr[c_layer]) | kInitCount);
				buffer[(c_layer * total_tree + tree_ind) * kMaxBuff + (pos_arr[c_layer * total_tree + tree_ind]++)] = kmer_bit;
				if (pos_arr[c_layer * total_tree + tree_ind] >= kMaxBuff) { // does not affect that much
					InsertBatch(root[(c_layer * total_tree + tree_ind)], &buffer[(c_layer * total_tree + tree_ind) * kMaxBuff], pos_arr[c_layer * total_tree + tree_ind]);
					kmer_count_arr[c_layer * total_tree + tree_ind] += pos_arr[c_layer * total_tree + tree_ind];
					pos_arr[c_layer * total_tree + tree_ind] = 0;
				}
			}
		}

		kseq_destroy(seq);
		gzclose(fp);
		putchar('#');
		fflush(stdout);
	}

	int max_len = -1;
	for (uint64_t i = 0; i < kTotalLayer; ++i) {
		for (uint64_t j = 0; j < total_tree; ++j) {
			InsertBatch(root[i * total_tree + j], &buffer[(i * total_tree + j) * kMaxBuff], pos_arr[i * total_tree + j]);
			kmer_count_arr[i * total_tree + j] += pos_arr[i * total_tree + j];
			if (kmer_count_arr[i * total_tree + j] > max_len)
				max_len = kmer_count_arr[i * total_tree + j];
		}
	}


	uint64_t t_sum = 0;
	for (uint64_t i = 0; i < kTotalLayer; ++i) {
		uint64_t t_s = 0;
		for (uint64_t j = 0; j < total_tree; ++j) {
				t_s += kmer_count_arr[i * total_tree + j];
		}
		std::cout << i << ": " << t_s << std::endl;
		t_sum += (i + 1) * t_s;
	}

	max_len <<= 1;
	if (max_len > 200 * 1024 * 1024)
		max_len = 200 * 1024 * 1024;

	kmer_arr = new uint64_t[max_len];

	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
#ifdef DEBUG
	std::cout << "mem used: " << usage.ru_maxrss << std::endl;
	std::cout << " kmer: " << kmer_count << "\n";
#endif

	delete[] buffer;

	// time for tree traversal starts
	TimePoint tp2 = Clock::now();
	std::clock_t cput2 = std::clock();

	for (uint64_t i = 1; i < kTotalLayer; ++i) {
		t_u=0;
		for (uint64_t j = 0; j < total_tree; ++j) {
			SortInsertKmer(i, j);
		}
#ifdef DEBUG
		putchar(i + 'a');
		std::cout << t_u << std::endl;
		fflush(stdout);
#endif
	}

	getrusage(RUSAGE_SELF, &usage);
#ifdef DEBUG
	std::cout << "next mem used: " << usage.ru_maxrss << std::endl;
#endif

	TimePoint tp3 = Clock::now();
	std::clock_t cput3 = std::clock();

	pre_arr = new uint64_t[10000];
	pre_arr[0] = (kmer_len << 32) | tree_pow;

	string suf ("suf_" + out_file_name);
	FILE* out_file = fopen(suf.c_str(), "wb");

	for (uint64_t i = 0; i < total_tree; ++i)
		SortKmer(root[i], i, out_file);

	fclose(out_file);

	string pre ("pre_" + out_file_name);
	FILE* pre_file = fopen(pre.c_str(), "wb");
	fwrite(pre_arr, sizeof(uint64_t), total_tree + 1, pre_file);
	fclose(pre_file);


	delete[] pre_arr;
	delete[] kmer_arr;
	delete[] temp_mem;

	CountBigKmer();


	TimePoint tp4 = Clock::now();
	std::clock_t cput4 = std::clock();



#ifdef DEBUG
	std::cout << bucket_fill << " out " << bucket_size << "\t ratio: " << (1.0 * bucket_fill / bucket_size) << "\n";
	std::cout << tc1 << " node: " << num_node << "\t" << num_bucket << std::endl;
	std::cout << "sort: " << time_radm << "\tin_sort: " << tmm1 << "\ttr_sort: " << tmm2 << "\n";
#endif
#ifdef DEBUG
	std::cout << " over " << over_c << "\t" << max_uniq << " :: " << tkc5 << " " <<  tkc4 << " " <<  tkc3 << ":" << tkc2 << ":" << tkc1 << " uniq: " << total_uniq_count << "\tkmer: " << kmer_count << "\n";
#endif
	std::cout << "unique kmer: " << total_uniq_count << "\ttotal kmer: " << kmer_count << "\n";
	std::cout << "insertion time: wall " << double(std::chrono::duration_cast<std::chrono::milliseconds>(tp2-tp1).count()) / 1000 << "\tcpu " << (cput2 - cput1) / (double)CLOCKS_PER_SEC << "\n";
	std::cout << "2nd phase time: wall " << double(std::chrono::duration_cast<std::chrono::milliseconds>(tp3-tp2).count()) / 1000 << "\tcpu " << (cput3 - cput2) / (double)CLOCKS_PER_SEC << "\n";
	std::cout << "3rd phase time: wall " << double(std::chrono::duration_cast<std::chrono::milliseconds>(tp4-tp3).count()) / 1000 << "\tcpu " << (cput4 - cput3) / (double)CLOCKS_PER_SEC << "\n";
	std::cout << "Total time: wall " << double(std::chrono::duration_cast<std::chrono::milliseconds>(tp4-tp1).count()) / 1000 << "\tcpu " << (cput4 - cput1) / (double)CLOCKS_PER_SEC << "\n";

#ifdef DEBUG
	for (int i = 0; i < 30; ++i)
		if (depth_c[i])
			std::cout << i << ":\t" << depth_c[i] << std::endl;
#endif
	std::time_t curr_time = std::time(nullptr);
	std::cout << "date_time: " << std::asctime(std::localtime(&curr_time)) << std::endl;
}
