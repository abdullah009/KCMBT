

#ifndef MSD_SORT8_H_
#define MSD_SORT8_H_

#include <vector>
#include <cmath>
#include <chrono>
#include <cstring>
#include <algorithm>
#include "common.h"


double time_radm = 0.0;
std::chrono::time_point<std::chrono::high_resolution_clock> tprm1, tprm2;

void DoRadixSort(uint64_t* data_begin, uint64_t* data_end, uint64_t *temp_begin, const uint64_t depth, int type) {
	uint64_t *temp_end = temp_begin + (data_end - data_begin);

	bool odd = 0;
	uint64_t limit = (64 - 2 * depth);
	for (uint64_t i = kTotalMaskBitArr[type]; i < limit; i += 8) {
		uint64_t count_arr[0x100] = {};
		for (uint64_t *p = data_begin; p != data_end; ++p)
			++count_arr[(*p >> i) & 0xFF];
		uint64_t* bucket_arr[0x100];
		uint64_t* q = temp_begin;
		for (uint64_t j = 0; j < 0x100; q += count_arr[j++])
			bucket_arr[j] = q;
		for (uint64_t* p = data_begin; p != data_end; ++p)
			*bucket_arr[(*p >> i) & 0xFF]++ = *p;
		std::swap(data_begin, temp_begin);
		std::swap(data_end, temp_end);

		odd ^= 1;
	}

	if (odd) {
		std::swap(data_begin, temp_begin);
		std::swap(data_end, temp_end);
		std::move(temp_begin, temp_end, data_begin);
	}
}

void RadixSortLSD8(uint64_t *begin, uint64_t *end, uint64_t *begin1, uint64_t maxshift, int type) {
	bool odd = 0;
	uint64_t *end1 = begin1 + (end - begin);
	//std::cout << type << " " << kTotalMaskBitArr[type] << " " << maxshift << std::endl;
	for (uint64_t shift = kTotalMaskBitArr[type]; shift < maxshift; shift += 8) {
		size_t count[0x100] = {};
		for (uint64_t *p = begin; p != end; p++)
			count[(*p >> shift) & 0xFF]++;
		uint64_t *bucket[0x100], *q = begin1;
		for (int i = 0; i < 0x100; q += count[i++])
			bucket[i] = q;
		for (uint64_t *p = begin; p != end; p++)
			*bucket[(*p >> shift) & 0xFF]++ = *p;
		std::swap(begin, begin1);
		std::swap(end, end1);
		odd ^= 1;
	}

	if (odd) {
		std::swap(begin, begin1);
		std::swap(end, end1);
	}
	else
		std::move(begin, end, begin1);
}

void DoMSDSort6(uint64_t *begin, uint64_t *end, uint64_t *begin1, uint64_t depth, int type)
{
	uint64_t shift = 58 - 2 * depth;
    uint64_t *end1 = begin1 + (end - begin);
    size_t count[0x40] = {};
    for (uint64_t *p = begin; p != end; p++)
        count[(*p >> shift) & 0x3F]++;
    uint64_t *bucket[0x40], *obucket[0x40], *q = begin1;
    for (int i = 0; i < 0x40; q += count[i++])
        obucket[i] = bucket[i] = q;
    for (uint64_t *p = begin; p != end; p++)
        *bucket[(*p >> shift) & 0x3F]++ = *p;
    for (int i = 0; i < 0x40; ++i)
    	RadixSortLSD8(obucket[i], bucket[i], begin + (obucket[i] - begin1), shift, type);
}

uint64_t MergeAndDivideMSD(uint64_t* bucket, uint64_t start, uint64_t end, uint64_t depth, int type) {
#ifdef TIME1
	tprm1 = std::chrono::high_resolution_clock::now();
#endif
	//std::cout << type << " " << kTotalMaskBitArr[type] << " d: " << depth << std::endl;
	DoMSDSort6(&bucket[start], &bucket[end], temp_mem, depth, type);
	//DoRadixSort(&bucket[start], &bucket[end], temp_mem, depth, type);
	//std::sort(&bucket[start], &bucket[end]);
#ifdef TIME1
	tprm2 = std::chrono::high_resolution_clock::now();
	time_radm += std::chrono::duration_cast<std::chrono::microseconds>(tprm2 - tprm1).count();
#endif

	uint64_t* first = &bucket[start], *last = &bucket[end];
	uint64_t* result = first;
	uint64_t count = *result & k_count_mask[type];
	uint64_t kmer_part = *result & k_kmer_mask[type];
	uint64_t count_mask = k_count_mask[type], kmer_mask = k_kmer_mask[type];

	while (++first != last) {
		uint64_t first_c = (*first & count_mask);
		uint64_t first_k = (*first & kmer_mask);
		if ((kmer_part != first_k) || (count + first_c >= count_mask)) {
			*result = kmer_part | count;
			*(++result) = *first;
			count = first_c;
			kmer_part = first_k;
		}
		else
			count += first_c;
	}
	*result = kmer_part | count;
	return std::distance(&bucket[0], ++result);
}

uint64_t Merge(uint64_t* bucket, uint64_t start2, uint64_t end2, int type) {
	uint64_t* first1 = &bucket[0], *last1 = &bucket[start2];
	uint64_t* first2 = &bucket[start2], *last2 = &bucket[end2];
	uint64_t* out = new uint64_t[end2];
	uint64_t* result = out;
	uint64_t dist = 0;
	while (true) {
		if (first1 == last1) {
			dist = std::distance(&out[0], result) + std::distance(first2, last2);
			std::copy(first2, last2, result);
			break;
		}
		if (first2 == last2) {
			dist = std::distance(&out[0], result) + std::distance(first1, last1);
			std::copy(first1, last1, result);
			break;
		}
		uint64_t k1 = *first1 & k_kmer_mask[type];
		uint64_t k2 = *first2 & k_kmer_mask[type];
		if (k1 < k2)
			*result++ = *first1++;
		else if (k1 > k2)
			*result++ = *first2++;
		else {
			uint64_t c1 = *first1 & k_count_mask[type];
			uint64_t c2 = *first2 & k_count_mask[type];
			if (c1 + c2 < k_count_mask[type]) {
				*result++ = k1 | (c1 + c2);
				first1++;
				first2++;
			}
			else {
				*result++ = *first1++;
				*result++ = *first2++;
			}
		}
	}
	std::move(&out[0], &out[dist], &bucket[0]);
	delete[] out;
	out = nullptr;
	return dist;
}


uint64_t MergeFinal(uint64_t* bucket, uint64_t start2, uint64_t end2, uint64_t tree_ind) {
	uint64_t count_mask = k_count_mask[0], kmer_mask = k_kmer_mask[0];
	uint64_t* out = new uint64_t[end2];
	uint64_t* result = out;

	uint64_t* first1 = &bucket[0], *last1 = &bucket[start2];
	uint64_t* temp1 = first1;
	++first1;
	while (first1 != last1) {
		if ((*temp1 & kmer_mask) != (*first1 & kmer_mask)) {
			*result++ = *temp1;
			temp1 = first1;
			++first1;
		}
		else {
			if (((*temp1 & count_mask) + (*first1 & count_mask)) >= count_mask) { // check & remove this line
				over_c++;
				BigKmer bkmer {(tree_ind << ((kmer_len - tree_pow) << 1)) | ((*temp1 & kmer_mask) >> ((32 - kmer_len + tree_pow) << 1)), (*temp1 & count_mask) + (*first1 & count_mask) - count_mask};
				bkmer_arr.push_back(bkmer);
				*temp1 = (*temp1 & kmer_mask) | count_mask;
			}
			else
				*temp1 = (*temp1 & kmer_mask) | ((*temp1 & count_mask) + (*first1 & count_mask));
			++first1;
		}
	}
	*result++ = *temp1;

	uint64_t* s1 = out;
	uint64_t* e1 = result;
	uint64_t* s2 = result;


	uint64_t* first2 = &bucket[start2], *last2 = &bucket[end2];
	uint64_t* temp2 = first2;
	++first2;
	while (first2 != last2) {
		if ((*temp2 & kmer_mask) != (*first2 & kmer_mask)) {
			*result++ = *temp2;
			temp2 = first2;
			++first2;
		}
		else {
			if (((*temp2 & count_mask) + (*first2 & count_mask)) > count_mask) { // check & remove this line
				over_c++;
				BigKmer bkmer {(tree_ind << ((kmer_len - tree_pow) << 1)) | ((*temp2 & kmer_mask) >> ((32 - kmer_len + tree_pow) << 1)), (*temp2 & count_mask) + (*first2 & count_mask) - count_mask};
				bkmer_arr.push_back(bkmer);
				*temp2 = (*temp2 & kmer_mask) | count_mask;
			}
			else
				*temp2 = (*temp2 & kmer_mask) | ((*temp2 & count_mask) + (*first2 & count_mask));
			++first2;
		}
	}
	*result++ = *temp2;

	uint64_t* e2 = result;

	result = &bucket[0];
	uint64_t dist = 0;
	while (true) {
		if (s1 == e1) {
			dist = std::distance(&bucket[0], result) + std::distance(s2, e2);
			std::copy(s2, e2, result);
			break;
		}
		if (s2 == e2) {
			dist = std::distance(&bucket[0], result) + std::distance(s1, e1);
			std::copy(s1, e1, result);
			break;
		}
		uint64_t k1 = *s1 & kmer_mask;
		uint64_t k2 = *s2 & kmer_mask;
		if (k1 < k2) {
			//if (__builtin_expect((*s1 & count_mask) < count_mask, 1))
			if ((*s1 & count_mask) < count_mask)
				*result++ = *s1;
			++s1;
		}
		else if (k1 > k2) {
			//if (__builtin_expect((*s2 & count_mask) < count_mask, 1))
			if ((*s2 & count_mask) < count_mask)
				*result++ = *s2;
			++s2;
		}
		else {
			uint64_t c1 = *s1 & count_mask;
			uint64_t c2 = *s2 & count_mask;
			if (c1 + c2 < count_mask)
				*result++ = k1 | (c1 + c2);
			else {
				BigKmer bkmer {(tree_ind << ((kmer_len - tree_pow) << 1)) | ((k1 & kmer_mask) >> ((32 - kmer_len + tree_pow) << 1)), c1 + c2 - count_mask};
				bkmer_arr.push_back(bkmer);
				over_c++;
			}
			s1++;
			s2++;
		}
	}
	delete[] out;
	out = nullptr;
	return dist;
}


uint64_t MergeFinalFull(uint64_t* bucket, uint64_t end, uint64_t tree_ind) {
	uint64_t count_mask = k_count_mask[0], kmer_mask = k_kmer_mask[0];
	uint64_t* out = new uint64_t[end];
	uint64_t* result = out;

	uint64_t* first1 = &bucket[0], *last1 = &bucket[end];
	uint64_t* temp1 = first1;
	++first1;
	while (first1 != last1) {
		if ((*temp1 & kmer_mask) != (*first1 & kmer_mask)) {
			//if (__builtin_expect((*temp1 & count_mask) != count_mask, 1))
			if ((*temp1 & count_mask) != count_mask)
				*result++ = *temp1;
			else {

			}
			temp1 = first1;
			++first1;
		}
		else {
			if (((*temp1 & count_mask) + (*first1 & count_mask)) >= count_mask) { // check & remove this line
				over_c++;
				BigKmer bkmer {(tree_ind << ((kmer_len - tree_pow) << 1)) | ((*temp1 & kmer_mask) >> ((32 - kmer_len + tree_pow) << 1)), (*temp1 & count_mask) + (*first1 & count_mask) - count_mask};
				bkmer_arr.push_back(bkmer);
				*temp1 = (*temp1 & kmer_mask) | count_mask;
			}
			else
				*temp1 = (*temp1 & kmer_mask) | ((*temp1 & count_mask) + (*first1 & count_mask));
			++first1;
		}
	}
	//if (__builtin_expect((*temp1 & count_mask) != count_mask, 1))
	if ((*temp1 & count_mask) != count_mask)
		*result++ = *temp1;

	uint64_t dist =  std::distance(&out[0], result);
	std::move(&out[0], &out[dist], &bucket[0]);
	delete[] out;
	out = nullptr;
	return dist;
}

#endif
