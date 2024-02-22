/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#ifndef STATISTICS_FILE_H_
#define STATISTICS_FILE_H_

#include <map>
#include <tuple>
#include "gtest.h"
#include "outside_data.h"

/**
 * @brief Create statistics files.
 *
 */
class StatisticsFile
{
public:
	/**
	 * @brief Construct a new StatisticsFile object
	 *
	 * @param optios Execution options.
	 */
	StatisticsFile(Options *optios);

	/**
	 * @brief Destroy the StatisticsFile object
	 *
	 */
	virtual ~StatisticsFile();

	/**
	 * @brief Set mer total count.
	 *
	 * @param mutant_mer_total Count of mutant total mer
	 * @param wildType_mer_total Count of wild type total mer
	 */
	void set_merCounter(const u_int64_t mutant_mer_total,
						const u_int64_t wildType_mer_total) const
	{
		this->gtest->set_merCounter(mutant_mer_total, wildType_mer_total);
	}

	/**
	 * @brief Create the statistics.txt file.
	 *
	 */
	void create_statisticsFile() const;

	/**
	 * @brief Create the outside.txt file.
	 *
	 * @param mutantMerPair Mutant mer pairs at each end
	 * @param wildTypeMerPair Wild type mer pairs at each end
	 */
	void create_outsideFile(
		const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &mutantMerPair,
		const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &wildTypeMerPair) const;

	// Setter / Getter

	void set_vectorArray(const std::string &vectorArray)
	{
		this->vectorArray = vectorArray;
	}

	std::string &get_vectorArray()
	{
		return this->vectorArray;
	}

	void set_vectorPosPair(
		const std::unordered_map<unsigned int, std::pair<std::string, std::string>> &vectorPosPair)
	{
		this->vectorPosPair = vectorPosPair;
	}

	std::unordered_map<unsigned int, std::pair<std::string, std::string>> get_vectorPosPair() const
	{
		return this->vectorPosPair;
	}

	void set_mutantPosFreq(const std::vector<unsigned int> &mutantPosFreq)
	{
		this->mutantPosFreq = mutantPosFreq;
	}

	void set_wildTypePosFreq(const std::vector<unsigned int> &posFreq)
	{
		this->wildTypePosFreq = posFreq;
	}

	std::unordered_map<unsigned int, double> get_fdr() const
	{
		return this->gtest->get_fdr();
	};

private:
	/**
	 * @brief Execution options.
	 *
	 */
	Options *options;

	/**
	 * @brief Run the G-test.
	 *
	 */
	Gtest *gtest;

	/**
	 * @brief Bases of the circular vector genome
	 *
	 */
	std::string vectorArray;

	/**
	 * @brief Position and k-mer complementary pair on vector
	 *
	 */
	std::unordered_map<unsigned int, std::pair<std::string, std::string>> vectorPosPair;

	/**
	 * @brief Position frequency of mutant
	 *
	 */
	std::vector<unsigned int> mutantPosFreq;

	/**
	 * @brief Position frequency of wild type
	 *
	 */
	std::vector<unsigned int> wildTypePosFreq;

	/**
	 * @brief Create outside data.
	 *
	 * @param mutantMerPair Mutant mer pairs at each end
	 * @param wildTypeMerPair Wild type mer pairs at each end
	 * @return Outside data
	 */
	std::tuple<unsigned int, std::unordered_map<unsigned int, size_t>, OutsideData> create_outsideData(
		const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &mutantMerPair,
		const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &wildTypeMerPair) const;
};
#endif /* STATISTICS_FILE_H_ */
