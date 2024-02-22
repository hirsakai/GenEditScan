/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#ifndef KMER_MATCH_H_
#define KMER_MATCH_H_

#include <map>
#include "bitwise_operation.h"
#include "statistics_file.h"
#include "fastq_match.h"

/**
 * @brief Match analysis of k-mer.
 *
 */
class KmerMatch
{
public:
	/**
	 * @brief Construct a new Kmer Match object
	 *
	 * @param options Execution options.
	 * @param bitwiseOperation Bitwise operation.
	 * @param statisticsFile Create statistics files.
	 */
	KmerMatch(Options *options, BitwiseOperation *bitwiseOperation,
			  StatisticsFile *statisticsFile);

	/**
	 * @brief Destroy the Kmer Match object
	 *
	 */
	virtual ~KmerMatch();

	/**
	 * @brief Execute match analysis of k-mer.
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
	 * @brief Input the read data for the match analysis.
	 *
	 */
	FastqMatch *fastqMatch;

	/**
	 * @brief Set position frequencies and write merFreq.txt files.
	 *
	 * @param mutantMerCounter Mutant mer counter
	 * @param wildTypeMerCounter Wild type mer counter
	 */
	void control_freqFile(
		const std::unordered_map<std::string, unsigned int> &mutantMerCounter,
		const std::unordered_map<std::string, unsigned int> &wildTypeMerCounter) const;

	/**
	 * @brief Set position frequencies.
	 *
	 * @param merCounter Counter of each mer
	 * @param vectorPosPair Position and k-mer complementary pair on vector
	 * @return Position frequencies
	 */
	std::vector<unsigned int> set_posFreq(
		const std::unordered_map<std::string, unsigned int> &merCounter,
		std::map<unsigned int, std::pair<std::string, std::string>> &vectorPosPair) const;

	/**
	 * @brief Create merFreq.txt file.
	 *
	 * @param merCounter Counter of each mer
	 * @param type '.mutant' or '.wildtype'
	 */
	void create_merFreqFile(const std::unordered_map<std::string, unsigned int> &merCounter,
							const std::string type) const;
};
#endif /* KMER_MATCH_H_ */
