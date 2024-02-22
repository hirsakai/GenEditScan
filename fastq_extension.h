/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#ifndef FASTQ_EXTENSION_
#define FASTQ_EXTENSION_

#include <string>
#include <unordered_map>
#include "bitwise_operation.h"

/**
 * @brief Input the read data for the extension analysis.
 *
 */
class FastqExtension
{
public:
	/**
	 * @brief Construct a new Fastq Extension object
	 *
	 * @param options Execution options.
	 * @param bitwiseOperation Bitwise operation.
	 */
	FastqExtension(Options *options, BitwiseOperation *bitwiseOperation);

	/**
	 * @brief Destroy the Fastq Extension object
	 *
	 */
	virtual ~FastqExtension();

	/**
	 * @brief Read the fastq.gz file.
	 *
	 * @param fastqFile FASTQ file
	 * @param merPair Mer pairs at each end
	 * @param merTotalCounter Mer total counter per file
	 * @return Mer pairs at each end for parallel processing
	 */
	std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> read_fastqFile(
		const std::string &fastqFile,
		const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &merPair,
		u_int64_t &merTotalCounter) const;

private:
	/**
	 * @brief Execution options.
	 *
	 */
	Options *options;

	/**
	 * @brief Bitwise operation.
	 *
	 */
	BitwiseOperation *bitwiseOperation;

	/**
	 * @brief Count k-mer.
	 *
	 * @param fastqFile FASTQ file
	 * @param fastqData FASTQ data
	 * @param merPair Mer pairs at each end
	 * @param merLocalPair Mer pairs at each end for parallel processing
	 * @param merTotalCounter Mer total counter per file
	 * @param readCounter Read counter per file
	 */
	void count_extension(
		const std::string &fastqFile, std::vector<std::string> &fastqData,
		const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &merPair,
		std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &merLocalPair,
		u_int64_t &merTotalCounter, u_int64_t &readCounter) const;
};
#endif /* FASTQ_EXTENSION_ */
