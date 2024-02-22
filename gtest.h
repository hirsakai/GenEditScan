/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#ifndef GTEST_H_
#define GTEST_H_

#include <tuple>
#include <unordered_map>
#include "options.h"

/**
 * @brief Run the G-test.
 *
 */
class Gtest
{
public:
	/**
	 * @brief Construct a new Gtest object
	 *
	 * @param options Execution options.
	 */
	Gtest(Options *options);

	/**
	 * @brief Destroy the Gtest object
	 */
	virtual ~Gtest();

	/**
	 * @brief Set mer total count.
	 *
	 * @param mutant_mer_total Count of mutant total mer
	 * @param wildType_mer_total Count of wild type total mer
	 */
	void set_merCounter(const u_int64_t mutant_mer_total,
						const u_int64_t wildType_mer_total)
	{
		this->mutant_mer_total = (double)mutant_mer_total;
		this->wildType_mer_total = (double)wildType_mer_total;
	}

	/**
	 * @brief Calculate G-value for k-mer match analysis.
	 *
	 * @param mutantPosFreq Position frequency of mutant
	 * @param wildTypePosFreq Position frequency of wild type
	 */
	void kmer_match(const std::vector<unsigned int> &mutantPosFreq,
					const std::vector<unsigned int> &wildTypePosFreq);

	/**
	 * @brief Calculate FDR using the Benjamini-Hochberg method.
	 *
	 * @param pval P-values
	 * @return FDR
	 */
	std::unordered_map<unsigned int, std::unordered_map<unsigned int, double>>
	fdr_extension(const std::unordered_map<unsigned int, std::unordered_map<unsigned int, double>> &pval) const;

	/**
	 * @brief Calculate G-value for k-mer extension analysis.
	 *
	 * @param mutant_count Count of mutant match mer
	 * @param wildType_count Count of wild type match mer
	 * @return G-value
	 */
	std::tuple<double, double> kmer_extension(const unsigned int mutant_count,
											  const unsigned int wildType_count) const;

	// Getter

	std::unordered_map<unsigned int, double> get_gval() const
	{
		return this->gval;
	};

	std::unordered_map<unsigned int, double> get_pval() const
	{
		return this->pval;
	};

	std::unordered_map<unsigned int, double> get_fdr() const
	{
		return this->fdr;
	};

	std::unordered_map<unsigned int, double> get_bon() const
	{
		return this->bon;
	};

private:
	/**
	 * @brief Execution options.
	 */
	Options *options;

	/**
	 * @brief Count of mutant total mer
	 */
	double mutant_mer_total;

	/**
	 * @brief Count of wild type total mer
	 */
	double wildType_mer_total;

	/**
	 * @brief G-value on vector array
	 */
	std::unordered_map<unsigned int, double> gval;

	/**
	 * @brief P-value on vector array
	 */
	std::unordered_map<unsigned int, double> pval;

	/**
	 * @brief FDR on vector array (Benjamini-Hochberg)
	 */
	std::unordered_map<unsigned int, double> fdr;

	/**
	 * @brief Bonferroni on vector array
	 */
	std::unordered_map<unsigned int, double> bon;

	/**
	 * @brief Williams's correction of G-value.
	 *
	 * @param mutant_total_log The product of log(count of mutant total mer) and (count of mutant total mer)
	 * @param wildType_total_log The product of log(count of wild type total mer) and (count of wild type total mer)
	 * @param total Count of mutant and wild type total mer
	 * @param q3 The product of log(total) and total
	 * @param qcomm Correction parameter
	 * @param mutant_match Count of mutant match mer
	 * @param wildType_match Count of wild type match mer
	 * @return Corrected G-value
	 */
	double adjusted_g(
		const double mutant_total_log, const double wildType_total_log, const double total,
		const double q3, const double qcomm,
		const double mutant_match, const double wildType_match) const;

	/**
	 * @brief Calculate FDR using the Benjamini-Hochberg method.
	 */
	void fdr_match();
};
#endif /* GTEST_H_ */
