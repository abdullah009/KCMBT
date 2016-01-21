# KCMBT: A very fast _k_-mer counter

Copyright 2015 Abdullah-Al Mamun <br />
abdullah.am.cs (at) engr.uconn.edu

 	KCMBT is a free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    KCMBT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with KCMBT. If not, see <http://www.gnu.org/licenses/>.

## What is KCMBT
KCMBT (_k_-mer Counter based on Multiple Burst Trees) is a very fast single-threaded _k_-mer counting algorithm. It uses cache efficient burst tries to store _k_-mers. Experimental results show that it outperforms all well-known algorithms.

## Compilation
To compile the source code, type

```
make
``` 

It will create two executable files in bin directory:

* _kcmbt_ generates binary files containing _k_-mers with their counts, 
* _kcmbt\_dump_ produces human readable file having _k_-mers with their counts.

## Usage
To go to bin directory:

	cd bin

Please run kcmbt at first, and then run kcmbt\_dump. kcmbt\_dump uses kcmbt generated files as input, and outputs list of human readable _k_-mers with their counts.

To run kcmbt, use

```
./kcmbt -k <k-mer length> -i <@file_listing_fastq_files or fastq_file> -o <output_file_name>
Parameters:
	-k:	k-mer length (10 <= k <= 32, default 28) 
	-i:	input file in fastq format (start with @ if the file contains a list of fastq files)
	-o:	output file (a binary file; please run kcmbt_dump to generate human readable output)
```
Example: ```kcmbt -k 28 -i srr.fastq -o srr```

To run kcmbt_dump, use

```
./kcmbt_dump -s <start_count> -e <end_count> -i <binary_file_name> -o <output_file_name>
Parameters:
	-s:	start count of k-mers (default 1)
	-e:	end count of k-mers (default 2^32)
	-i:	input file (must be the output_file_name used in kcmbt)
	-o:	output file (default "kmer_list")
```
Example: ```./kcmbt_dump -s 2 -e 1000 -i srr -o kmer_list```


## Contact
For questions, suggestions, bugs, and other issues, please contact:

```
Abdullah-Al Mamun
abdullah.am.cs (at) engr.uconn.edu
```
