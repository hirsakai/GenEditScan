/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#ifndef KMER_EXTENSION_H_
#define KMER_EXTENSION_H_

#include "bitwise_operation.h"
#include "statistics_file.h"
#include "fastq_extension.h"

/**
 * @brief Extension analysis of k-mer.
 *
 */
class KmerExtension
{
public:
	/**
	 * @brief Construct a new Kmer Extension object
	 *
	 * @param options Execution options.
	 * @param bitwiseOperation Bitwise operation.
	 * @param statisticsFile Create statistics files.
	 */
	KmerExtension(Options *options, BitwiseOperation *bitwiseOperation,
				  StatisticsFile *statisticsFile);

	/**
	 * @brief Destroy the Kmer Extension object
	 *
	 */
	virtual ~KmerExtension();

	/**
	 * @brief Execute extension analysis of k-mer.
	 *
	 */
	void execution() const;

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
	 * @brief Run statistical analysis.
	 *
	 */
	StatisticsFile *statisticsFile;

	/**
	 * @brief Input the read data for the extension analysis.
	 *
	 */
	FastqExtension *fastqExtension;

	/**
	 * @brief Set k-mer in hash table.
	 *
	 * @param fdr FDR on vector array
	 * @param merCounter Mer counter at each end
	 * @param chunk For bitwise operation. Chunk array.
	 * @return Count of target mer
	 */
	unsigned int set_merCounter(
		std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &merCounter) const;

	/**
	 * @brief Create chunk array.
	 *
	 * @param merCounter Mer counter at each end
	 * @param chunk For bitwise operation. Chunk array.
	 */
	void create_chunk(
		const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &merCounter) const;
};
#endif /* KMER_EXTENSION_H_ */
