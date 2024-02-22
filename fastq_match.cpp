/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#include <iostream>
#include <zlib.h>
#include "fastq_match.h"

/**
 * @brief Construct a new Fastq Match:: Fastq Match object
 *
 * @param options Execution options.
 * @param bitwiseOperation Bitwise operation.
 */
FastqMatch::FastqMatch(Options *options, BitwiseOperation *bitwiseOperation)
{
	this->options = options;
	this->bitwiseOperation = bitwiseOperation;
}

/**
 * @brief Destroy the Fastq Match:: Fastq Match object
 *
 */
FastqMatch::~FastqMatch()
{
}

/**
 * @brief Read the fastq.gz file.
 *
 * @param fastqFile FASTQ file
 * @param merCounter Counter of each mer
 * @param merTotalCounter Mer total counter per file
 * @return Counter of each mer for parallel processing
 */
std::unordered_map<std::string, unsigned int> FastqMatch::read_fastqFile(
	const std::string &fastqFile,
	const std::unordered_map<std::string, unsigned int> &merCounter,
	u_int64_t &merTotalCounter) const
{
	// File mode
	const gzFile file = gzopen(fastqFile.c_str(), "rb");
	if (!file)
	{
		std::cerr << "[Error] Could not open (" << fastqFile << ")." << std::endl;
		std::exit(1);
	}

	const unsigned int kmerLen = this->options->kmer;
	const unsigned int max_buff = this->options->max_read_length + 2;
	char buff[max_buff];
	std::string aLine[4];
	unsigned int nLine = 0;
	u_int64_t readCounter = 0;
	std::vector<std::string> fastqData;
	std::unordered_map<std::string, unsigned int> merLocalCounter;

	while (gzgets(file, buff, max_buff) != Z_NULL)
	{
		aLine[nLine++] = std::string(buff);
		if (nLine == 4)
		{
			nLine = 0;
			if (aLine[0][0] != '@' || aLine[2][0] != '+')
			{
				aLine[0].pop_back();
				std::cerr << "[Error] Could not get sequence (" << aLine[0] << ")." << std::endl;
				exit(1);
			}
			else
			{
				if (aLine[1].length() > kmerLen)
				{
					// Delete line break (\n)
					aLine[1].pop_back();
					fastqData.push_back(aLine[1]);
					if (fastqData.size() > this->options->fastq_read_lines)
					{
						this->count_match(fastqFile, fastqData, merCounter,
										  merLocalCounter, merTotalCounter, readCounter);
						fastqData.clear();
					}
				}
			}
		}
	}

	this->count_match(fastqFile, fastqData, merCounter,
					  merLocalCounter, merTotalCounter, readCounter);
	fastqData.clear();
	gzclose(file);
	return merLocalCounter;
}

//============================================================================//
// Private function
//============================================================================//
/**
 * @brief Count k-mer.
 *
 * @param fastqFile FASTQ file
 * @param fastqData FASTQ data
 * @param merCounter Counter of each mer
 * @param merLocalCounter Counter of each mer for parallel processing
 * @param merTotalCounter Mer total counter per file
 * @param readCounter Read counter per file
 */
void FastqMatch::count_match(
	const std::string &fastqFile, std::vector<std::string> &fastqData,
	const std::unordered_map<std::string, unsigned int> &merCounter,
	std::unordered_map<std::string, unsigned int> &merLocalCounter,
	u_int64_t &merTotalCounter, u_int64_t &readCounter) const
{
	const unsigned int kmer = this->options->kmer;
	const unsigned int mask = this->options->max_chunk_array;
	const unsigned int chunk_length = this->options->chunk_length;
	unsigned char *dna2bit = this->bitwiseOperation->get_dna2bit();
	unsigned char *chunk = this->bitwiseOperation->get_chunk();
	unsigned int dnabit, j;
	std::string mer;

#ifdef _OPENMP
#pragma omp parallel for num_threads(this->options->inner_parallel) private(dnabit, j, mer) \
	reduction(+ : merTotalCounter)
#endif
	for (size_t i = 0; i < fastqData.size(); i++)
	{
#ifdef _OPENMP
#pragma omp critical(match)
#endif
		if (++readCounter % this->options->log_output_interval == 0)
		{
			std::cerr << fastqFile << ": parsing " << readCounter
					  << " reads (k-mer match)." << std::endl;
		}

		dnabit = dna2bit[(unsigned char)fastqData[i][0]];
		for (j = 1; j < chunk_length - 1; j++)
		{
			dnabit = (dnabit << 2) + dna2bit[(unsigned char)fastqData[i][j]];
		}

		for (j = 0; j <= fastqData[i].length() - kmer; j++)
		{
			dnabit = (dnabit << 2) + dna2bit[(unsigned char)fastqData[i][chunk_length - 1 + j]];
			dnabit = dnabit & mask;
			if (chunk[dnabit] == 1 || dnabit == mask)
			{
				mer.assign(fastqData[i], j, kmer);
				if (merCounter.find(mer) != merCounter.end())
				{
#ifdef _OPENMP
#pragma omp atomic
#endif
					merLocalCounter[mer]++;
				}
			}
			merTotalCounter++;
		}
	}
}
