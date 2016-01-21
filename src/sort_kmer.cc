
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <fstream>
#include <limits>
#include <cstdlib>

using namespace std;

struct Kmer {
	string kmer;
	uint64_t count;
};

bool Compare(const Kmer& k1, const Kmer& k2) {
	return k1.kmer < k2.kmer;
}

int main(int argc, char** argv) {
	if (argc != 3 && argc != 5) {
		cerr << "./a.out in out begin_count end_count\n";
		return EXIT_FAILURE;
	}
	string in_str = argv[1], out_str = argv[2];
	int begin_count = 0, end_count = numeric_limits<int>::max();
	if (argc > 3) {
		begin_count = atoi(argv[3]);
		end_count = atoi(argv[4]);
	}

	if (begin_count < 0 || end_count < 0 || begin_count > end_count) {
		cerr << "count information error\n";
		return EXIT_FAILURE;
	}

	ifstream in(in_str);
	vector<Kmer> kmer_arr;
	string s;
	uint64_t c;
	while (in >> s >> c) {
		if (c >= begin_count && c <= end_count)
			kmer_arr.push_back({s, c});
	}
	in.close();

	sort(kmer_arr.begin(), kmer_arr.end(), Compare);

	ofstream out(out_str);
	for (auto kmer : kmer_arr)
		out << kmer.kmer << "\t" << kmer.count << "\n";
	out.close();

	return 0;
}
