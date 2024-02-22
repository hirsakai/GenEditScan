/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#include <map>
#include <cmath>
#include <algorithm>
#include "gtest.h"

/**
 * @brief Complemented Chi square.
 *
 * @param double Degree of freedom
 * @param double G-value (df = 1)
 *
 * @return P-value
 */
extern "C" double chdtrc(double, double);

/**
 * @brief Construct a new Gtest:: Gtest object
 *
 * @param options Execution options.
 */
Gtest::Gtest(Options *options)
{
	this->options = options;
}

/**
 * @brief Destroy the Gtest:: Gtest object
 *
 */
Gtest::~Gtest()
{
}

/**
 * @brief Calculate G-value for k-mer match analysis.
 *
 * @param mutantPosFreq Position frequency of mutant
 * @param wildTypePosFreq Position frequency of wild type
 */
void Gtest::kmer_match(const std::vector<unsigned int> &mutantPosFreq,
					   const std::vector<unsigned int> &wildTypePosFreq)
{
	// for total
	const double mutant_mer_total_log = log(this->mutant_mer_total) * this->mutant_mer_total;
	const double wildType_mer_total_log = log(this->wildType_mer_total) * this->wildType_mer_total;
	const double mer_total = this->mutant_mer_total + this->wildType_mer_total;
	const double mer_q3 = log(mer_total) * mer_total;
	const double mer_qcomm = (mer_total / this->mutant_mer_total + mer_total / this->wildType_mer_total - 1.0) / (6.0 * mer_total);

	std::map<std::pair<unsigned int, unsigned int>, double> gval_stock;
	std::map<std::pair<unsigned int, unsigned int>, double> pval_stock;
	std::map<std::pair<unsigned int, unsigned int>, double> bon_stock;
	const size_t vector_len = mutantPosFreq.size();

	for (size_t i = 0; i < vector_len; i++)
	{
		std::pair<unsigned int, unsigned int> target = std::make_pair(mutantPosFreq[i], wildTypePosFreq[i]);
		if (gval_stock.find(target) != gval_stock.end())
		{
			this->gval[i] = gval_stock[target];
			this->pval[i] = pval_stock[target];
			this->bon[i] = bon_stock[target];
		}
		else
		{
			// for total
			const double mutant_mer_match = mutantPosFreq[i];
			const double wildType_mer_match = wildTypePosFreq[i];

			if (mutant_mer_match * this->wildType_mer_total > wildType_mer_match * this->mutant_mer_total)
			{
				// G-value
				this->gval[i] = this->adjusted_g(
					mutant_mer_total_log, wildType_mer_total_log, mer_total,
					mer_q3, mer_qcomm, mutant_mer_match, wildType_mer_match);
				// P-value; To avoid igamc underflow error.
				// When Gval is greater than 170, Pval is less than 1.175494e-38 (float limit).
				this->pval[i] = this->gval[i] > 0.0 ? (this->gval[i] < 170.0 ? chdtrc(1.0, this->gval[i]) : 0.0) : 1.0;
				// Bonferroni
				this->bon[i] = std::min(this->pval[i] * vector_len, 1.0);
			}
			else
			{
				this->gval[i] = 0.0;
				this->pval[i] = 1.0;
				this->bon[i] = 1.0;
			}

			gval_stock[target] = this->gval[i];
			pval_stock[target] = this->pval[i];
			bon_stock[target] = this->bon[i];
		}
	}

	// Calculate FDR using the Benjamini-Hochberg method.
	this->fdr_match();
}

/**
 * @brief Calculate G-value for k-mer extension analysis.
 *
 * @param mutant_count Count of mutant match mer
 * @param wildType_count Count of wild type match mer
 * @return G-value
 */
std::tuple<double, double> Gtest::kmer_extension(const unsigned int mutant_count,
												 const unsigned int wildType_count) const
{
	const double mutant_mer_total_log = log(this->mutant_mer_total) * this->mutant_mer_total;
	const double wildType_mer_total_log = log(this->wildType_mer_total) * this->wildType_mer_total;
	const double mer_total = this->mutant_mer_total + this->wildType_mer_total;
	const double mer_q3 = log(mer_total) * mer_total;
	const double mer_qcomm = (mer_total / this->mutant_mer_total + mer_total / this->wildType_mer_total - 1.0) / (6.0 * mer_total);

	const double mutant_mer_match = (double)mutant_count;
	const double wildType_mer_match = (double)wildType_count;

	if (mutant_mer_match * this->wildType_mer_total >= wildType_mer_match * this->mutant_mer_total)
	{
		// G-value
		const double g = this->adjusted_g(
			mutant_mer_total_log, wildType_mer_total_log, mer_total,
			mer_q3, mer_qcomm, mutant_mer_match, wildType_mer_match);
		// P-value; To avoid igamc underflow error.
		// When Gval is greater than 170, Pval is less than 1.175494e-38 (float limit).
		const double p = g > 0.0 ? (g < 170.0 ? chdtrc(1.0, g) : 0.0) : 1.0;
		return {g, p};
	}
	else
	{
		return {0.0, 1.0};
	}
}

/**
 * @brief Calculate FDR using the Benjamini-Hochberg method.
 *
 * @param pval P-values
 * @return FDR
 */
std::unordered_map<unsigned int, std::unordered_map<unsigned int, double>>
Gtest::fdr_extension(const std::unordered_map<unsigned int, std::unordered_map<unsigned int, double>> &pval) const
{
	std::vector<std::pair<double, unsigned int>> v;
	unsigned int ic = 0;
	for (auto itr = pval.begin(); itr != pval.end(); itr++)
	{
		for (auto itr_p = itr->second.begin(); itr_p != itr->second.end(); itr_p++)
		{
			v.push_back(std::make_pair(itr_p->second, ic++));
		}
	}
	sort(v.begin(), v.end());

	std::unordered_map<unsigned int, double> f;
	const double vector_len = (double)v.size();
	double vector_pos = 1.0;
	double pval_prev = v[0].first;
	double fdr_prev = std::min(pval_prev * vector_len, 1.0);

	for (auto itr = v.begin(); itr != v.end(); itr++)
	{
		if (itr->first == pval_prev)
		{
			f[itr->second] = fdr_prev;
		}
		else
		{
			pval_prev = itr->first;
			f[itr->second] = std::min(itr->first * vector_len / vector_pos, 1.0);
			fdr_prev = f[itr->second];
		}
		vector_pos += 1.0;
	}

	std::unordered_map<unsigned int, std::unordered_map<unsigned int, double>> fdr_map;
	ic = 0;
	for (auto itr = pval.begin(); itr != pval.end(); itr++)
	{
		for (auto itr_p = itr->second.begin(); itr_p != itr->second.end(); itr_p++)
		{
			fdr_map[itr->first][itr_p->first] = f[ic++];
		}
	}
	return fdr_map;
}

//============================================================================//
// Private function
//============================================================================//
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
double Gtest::adjusted_g(
	const double mutant_total_log, const double wildType_total_log,
	const double total, const double q3, const double qcomm,
	const double mutant_match, const double wildType_match) const
{
	const double mutant_match_log = log(mutant_match) * mutant_match;
	const double wildType_match_log = log(wildType_match) * wildType_match;

	const double mutant_notmatch = this->mutant_mer_total - mutant_match;
	const double mutant_notmatch_log = log(mutant_notmatch) * mutant_notmatch;
	const double wildType_notmatch = this->wildType_mer_total - wildType_match;
	const double wildType_notmatch_log = log(wildType_notmatch) * wildType_notmatch;

	const double match = mutant_match + wildType_match;
	const double match_log = log(match) * match;
	const double notmatch = mutant_notmatch + wildType_notmatch;
	const double notmatch_log = log(notmatch) * notmatch;

	double q1;
	if (mutant_match == 0.0 && wildType_match == 0.0)
	{
		q1 = mutant_notmatch_log + wildType_notmatch_log;
	}
	else if (wildType_match == 0.0)
	{
		q1 = mutant_match_log + mutant_notmatch_log + wildType_notmatch_log;
	}
	else if (mutant_match == 0.0)
	{
		q1 = mutant_notmatch_log + wildType_match_log + wildType_notmatch_log;
	}
	else
	{
		q1 = mutant_match_log + mutant_notmatch_log + wildType_match_log + wildType_notmatch_log;
	}

	const double q2 = match == 0.0 ? mutant_total_log + wildType_total_log + notmatch_log : mutant_total_log + wildType_total_log + match_log + notmatch_log;

	// Compute G
	const double g = 2.0 * (q1 - q2 + q3);

	// Williams's correction for 2 x 2 table is ...
	const double q = match == 0.0 ? 1.0 + (total / notmatch - 1.0) * qcomm : 1.0 + (total / match + total / notmatch - 1.0) * qcomm;
	return g / q;
}

/**
 * @brief Calculate FDR using the Benjamini-Hochberg method.
 */
void Gtest::fdr_match()
{
	std::vector<std::pair<double, unsigned int>> v;
	for (auto itr = this->pval.begin(); itr != this->pval.end(); itr++)
	{
		v.push_back(std::make_pair(itr->second, itr->first));
	}
	sort(v.begin(), v.end());

	const double vector_len = (double)v.size();
	double vector_pos = 1.0;
	double pval_prev = v[0].first;
	double fdr_prev = std::min(pval_prev * vector_len, 1.0);

	for (auto itr = v.begin(); itr != v.end(); itr++)
	{
		if (itr->first == pval_prev)
		{
			this->fdr[itr->second] = fdr_prev;
		}
		else
		{
			pval_prev = itr->first;
			this->fdr[itr->second] = std::min(itr->first * vector_len / vector_pos, 1.0);
			fdr_prev = this->fdr[itr->second];
		}
		vector_pos += 1.0;
	}
}
