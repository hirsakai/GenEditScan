/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#ifndef FASTQ_MATCH_H_
#define FASTQ_MATCH_H_

#include <string>
#include <unordered_map>
#include "bitwise_operation.h"

/**
 * @brief Input the read data for the match analysis.
 */
class FastqMatch
{
public:
	/**
	 * @brief Construct a new Fastq Match object
	 *
	 * @param options Execution options.
	 * @param bitwiseOperation FFor bitwise operation. DNA expressed in 2 bits.
	 */
	FastqMatch(Options *options, BitwiseOperation *bitwiseOperation);

	/**
	 * @brief Destroy the Fastq Match object
	 *
	 */
	virtual ~FastqMatch();

	/**
	 * @brief Read the fastq.gz file.
	 *
	 * @param fastqFile FASTQ file
	 * @param merCounter Counter of each mer
	 * @param merTotalCounter Mer total counter per file
	 * @return Counter of each mer for parallel processing
	 */
	std::unordered_map<std::string, unsigned int> read_fastqFile(
		const std::string &fastqFile,
		const std::unordered_map<std::string, unsigned int> &merCounter,
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
	 * @param merCounter Counter of each mer
	 * @param merLocalCounter Counter of each mer for parallel processing
	 * @param merTotalCounter Mer total counter per file
	 * @param readCounter Read counter per file
	 */
	void count_match(
		const std::string &fastqFile, std::vector<std::string> &fastqData,
		const std::unordered_map<std::string, unsigned int> &merCounter,
		std::unordered_map<std::string, unsigned int> &merLocalCounter,
		u_int64_t &merTotalCounter, u_int64_t &readCounter) const;
};
#endif /* FASTQ_MATCH_H_ */
