#ifndef BURST_SORT_KMER_H_
#define BURST_SORT_KMER_H_


#include "msd_sort.h"
#include "common.h"



#ifdef DEBUG
uint64_t bucket_fill = 0;
uint64_t bucket_size = 0;
uint64_t node_stat_arr[total_tree] = {0};
uint64_t bucket_fill_arr[total_tree] = {0};
uint64_t bucket_size_arr[total_tree] = {0};
uint64_t tc1 = 0;
uint64_t num_bucket = 0, num_node = 0;
double time_node = 0.0, tmm1 = 0.0, tmm2 = 0.0;
std::chrono::time_point<std::chrono::high_resolution_clock> tpn1, tpn2, tpp1, tpp2;
uint64_t depth_c[30] = {0};
int ct = 0;
#endif

const int shr_arr[32] =
{
		62, 60, 58, 56, 54, 52, 50, 48, 46, 44, 42, 40, 38, 36, 34, 32, 30, 28, 26, 24, 22, 20, 18, 16, 14, 12, 10, 8
};


struct TrieNode;

struct TrieNode {
	uint64_t* bucket_arr[kAlphabetSize] __attribute__((aligned(64))) = {nullptr, nullptr, nullptr, nullptr};
	TrieNode* trie_arr[kAlphabetSize] = {nullptr, nullptr, nullptr, nullptr};
	uint64_t size_arr[kAlphabetSize] = {0, 0, 0, 0};
	uint64_t pos_arr[kAlphabetSize] = {0, 0, 0, 0};
	uint64_t sorted_pos_arr[kAlphabetSize] = {0, 0, 0, 0};
	bool trie_flag[kAlphabetSize] = {false, false, false, false};
	int type;
	TrieNode(int t) : type(t) {}
	void InitRoot() {
		for (uint64_t i = 0; i < kAlphabetSize; ++i) {
#ifdef MALLOC
			bucket_arr[i]  = (uint64_t *) malloc(kRootBucketSize * sizeof(uint64_t));//new uint64_t[kRootBucketSize];
#endif
#ifdef NEW
			bucket_arr[i]  = new uint64_t[kRootBucketSize];
#endif
			size_arr[i] = kRootBucketSize;
		}
	}

	void SplitBucket(uint64_t* bucket, uint64_t buc_size, uint64_t depth) {
		uint64_t num_shift = (31 - depth) << 1;
		int ch = (bucket[0] >> num_shift) & 3;
		int s = 0;
		uint64_t i = 1;
		for (; i < buc_size; ++i) {
			if (ch != ((bucket[i] >> num_shift) & 3)) {
				pos_arr[ch] = sorted_pos_arr[ch] = i - s;
				size_arr[ch] = pos_arr[ch] << 1;
				if (size_arr[ch] < kMinBucketSize)
					size_arr[ch] = kMinBucketSize;
#ifdef MALLOC
				bucket_arr[ch] = (uint64_t *) malloc(size_arr[ch] * sizeof(uint64_t));//new uint64_t[size_arr[ch]];
#endif
#ifdef NEW
				bucket_arr[ch] = new uint64_t[size_arr[ch]];
#endif
				std::move(&bucket[s], &bucket[i], &bucket_arr[ch][0]);
				s = i;
				ch = (bucket[i] >> num_shift) & 3;
			}
		}
		pos_arr[ch] = sorted_pos_arr[ch] = i - s;
		size_arr[ch] = pos_arr[ch] << 1;
		if (size_arr[ch] < kMinBucketSize)
			size_arr[ch] = kMinBucketSize;
#ifdef MALLOC
		bucket_arr[ch] = (uint64_t *) malloc(size_arr[ch] * sizeof(uint64_t));//new uint64_t[size_arr[ch]];
#endif
#ifdef NEW
		bucket_arr[ch] = new uint64_t[size_arr[ch]];
#endif

		std::move(&bucket[s], &bucket[i], &bucket_arr[ch][0]);
		s = i;
		ch = (bucket[i] >> num_shift) & 3;

		for (uint64_t i = 0; i < kAlphabetSize; ++i) {
			if (size_arr[i] <= 0) {
				size_arr[i] = kMinBucketSize;
#ifdef MALLOC
				bucket_arr[i] = (uint64_t *) malloc(size_arr[i] * sizeof(uint64_t));//new uint64_t[size_arr[i]];
#endif
#ifdef NEW
				bucket_arr[i] = new uint64_t[size_arr[i]];
#endif

			}
		}
#ifdef MALLOC
		free(bucket);
#endif
#ifdef NEW
		delete[] bucket;
#endif
		bucket = nullptr;
	}
};

uint64_t tkc2 = 0, tkc4 =  0, tkc5 = 0;
void InsertBatch(TrieNode* root, uint64_t* l_kmer_arr, uint64_t kmer_total) {

	for (uint64_t i = 0; i < kmer_total; ++i) {
		uint64_t kmer = l_kmer_arr[i];
		uint64_t depth = 0;
		uint64_t ch = (kmer >> shr_arr[depth]) & 3;
		TrieNode* node = root;
		while (node->trie_flag[ch]) {
			node = node->trie_arr[ch];
			depth += 1;
			ch = (kmer >> shr_arr[depth]) & 3;
		}
		uint64_t* bucket = node->bucket_arr[ch];
		bucket[node->pos_arr[ch]++] = kmer;
		if (node->pos_arr[ch] >= node->size_arr[ch]) {
			if (depth >= kmer_len) {
				std::cout << "MAX DEPTH" << std::endl;
	#ifdef DEBUG
				++tc1;
	#endif
				// resize array to array * 2
				uint64_t s = node->size_arr[ch] << 1;
#ifdef MALLOC
				node->bucket_arr[ch] = (uint64_t *) realloc(node->bucket_arr[ch], s * sizeof(uint64_t));
#endif
#ifdef NEW
				uint64_t* temp = new uint64_t[s];
				std::move(&bucket[0], &bucket[node->pos_arr[ch]], &temp[0]);
				delete[] node->bucket_arr[ch];
				node->bucket_arr[ch] = temp;
#endif
				node->size_arr[ch] = s;
			}
			else if (node->size_arr[ch] < kMaxBucketSize) {
	#ifdef DEBUG
				++tc1;
	#endif
				uint64_t kmer_in_bucket = MergeAndDivideMSD(bucket, node->sorted_pos_arr[ch], node->pos_arr[ch], depth, node->type);
				kmer_in_bucket = Merge(bucket, node->sorted_pos_arr[ch], kmer_in_bucket, node->type);
				node->pos_arr[ch] = node->sorted_pos_arr[ch] = kmer_in_bucket;
				uint64_t s = node->pos_arr[ch] << 1;
				if (s > kMaxBucketSize)
					s = kMaxBucketSize;
#ifdef MALLOC
				node->bucket_arr[ch] = (uint64_t *) realloc(node->bucket_arr[ch], s * sizeof(uint64_t));
#endif
#ifdef NEW
				uint64_t* temp = new uint64_t[s];
				std::move(&bucket[0], &bucket[node->pos_arr[ch]], &temp[0]);
				delete[] node->bucket_arr[ch];
				node->bucket_arr[ch] = temp;
#endif
				node->size_arr[ch] = s;
			}
			else {
	#ifdef TIME3
					tpn1 = std::chrono::high_resolution_clock::now();
	#endif
				uint64_t kmer_in_bucket = MergeAndDivideMSD(bucket, node->sorted_pos_arr[ch], node->pos_arr[ch], depth, node->type);
				kmer_in_bucket = Merge(bucket, node->sorted_pos_arr[ch], kmer_in_bucket, node->type);
				node->trie_arr[ch] = new TrieNode(node->type);
				node->trie_flag[ch] = true;
				node->trie_arr[ch]->SplitBucket(bucket, kmer_in_bucket, depth + 1);
	#ifdef TIME3
					tpn2 = std::chrono::high_resolution_clock::now();
					tmm1 += std::chrono::duration_cast<std::chrono::microseconds>(tpn2 - tpn1).count();
	#endif
			}
		}
	}
}

static void Traverse(TrieNode* node, uint64_t depth, uint64_t& kmer_uniq)
{
	for (uint64_t i = 0; i < kAlphabetSize; ++i) {
		if (node->trie_flag[i])
			Traverse(node->trie_arr[i], depth + 1, kmer_uniq);
		else {
			uint64_t* bucket = node->bucket_arr[i];
			if (node->pos_arr[i] > 0) {
#ifdef TIME2
				tpn1 = std::chrono::high_resolution_clock::now();
#endif

				uint64_t kmer_in_bucket = node->pos_arr[i];
				if (node->sorted_pos_arr[i] < node->pos_arr[i]) {
					kmer_in_bucket = MergeAndDivideMSD(bucket, node->sorted_pos_arr[i], node->pos_arr[i], depth, node->type);
					if (node->sorted_pos_arr[i])
						kmer_in_bucket = Merge(bucket, node->sorted_pos_arr[i], kmer_in_bucket, node->type);
				}
				std::move(&bucket[0], &bucket[kmer_in_bucket], &kmer_arr[kmer_uniq]);
				kmer_uniq += kmer_in_bucket;

#ifdef TIME2
				tpn2 = std::chrono::high_resolution_clock::now();
				tmm2 += std::chrono::duration_cast<std::chrono::microseconds>(tpn2 - tpn1).count();
#endif

#ifdef DEBUG
				depth_c[depth] ++;
				//bucket_fill_arr[node->ind] += kmer_in_bucket;
				//bucket_size_arr[node->ind] += node->size_arr[i];
#endif
			}
#ifdef DEBUG
			++num_bucket;
#endif


#ifdef MALLOC
			free(bucket);
#endif
#ifdef NEW
			delete[] bucket;
#endif
			bucket = nullptr;
		}
	}
#ifdef DEBUG
	++num_node;
#endif
	delete node;
	node = nullptr;
}


static void TraverseFinal(TrieNode* node, uint64_t depth, uint64_t& kmer_uniq, uint64_t tree_ind)
{
	for (uint64_t i = 0; i < kAlphabetSize; ++i) {
		if (node->trie_flag[i])
			TraverseFinal(node->trie_arr[i], depth + 1, kmer_uniq, tree_ind);
		else {
			uint64_t* bucket = node->bucket_arr[i];
			if (node->pos_arr[i] > 0) {
#ifdef TIME2
				tpn1 = std::chrono::high_resolution_clock::now();
#endif

				uint64_t kmer_in_bucket = node->pos_arr[i];
				if (node->sorted_pos_arr[i] < node->pos_arr[i]) {
					//kmer_in_bucket = MergeAndDivideMSD(bucket, 0, node->pos_arr[i], depth, node->type);
					kmer_in_bucket = MergeAndDivideMSD(bucket, node->sorted_pos_arr[i], node->pos_arr[i], depth, node->type);
					if (node->sorted_pos_arr[i] > 0)
						kmer_in_bucket = MergeFinal(bucket, node->sorted_pos_arr[i], kmer_in_bucket, tree_ind);
					else
						kmer_in_bucket = MergeFinalFull(bucket, kmer_in_bucket, tree_ind);

				}

#ifdef TIME2
				tpn2 = std::chrono::high_resolution_clock::now();
				tmm2 += std::chrono::duration_cast<std::chrono::microseconds>(tpn2 - tpn1).count();
#endif

				std::move(&bucket[0], &bucket[kmer_in_bucket], &kmer_arr[kmer_uniq]);
				kmer_uniq += kmer_in_bucket;

#ifdef DEBUG
				depth_c[depth] ++;
				//bucket_fill_arr[node->ind] += kmer_in_bucket;
				//bucket_size_arr[node->ind] += node->size_arr[i];
#endif
			}
#ifdef DEBUG
			++num_bucket;
#endif


#ifdef MALLOC
			free(bucket);
#endif
#ifdef NEW
			delete[] bucket;
#endif
			bucket = nullptr;
		}
	}
#ifdef DEBUG
	++num_node;
#endif
	delete node;
	node = nullptr;
}

#endif
