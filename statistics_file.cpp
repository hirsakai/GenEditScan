/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#include <fstream>
#include <algorithm>
#include "statistics_file.h"
#include "complementary.h"

/**
 * @brief Construct a new StatisticsFile:: StatisticsFile object
 *
 * @param options Execution options.
 */
StatisticsFile::StatisticsFile(Options *options)
{
	this->options = options;
	this->gtest = new Gtest(this->options);
}

/**
 * @brief Destroy the StatisticsFile:: StatisticsFile object
 *
 */
StatisticsFile::~StatisticsFile()
{
}

/**
 * @brief Create the statistics.txt file.
 *
 */
void StatisticsFile::create_statisticsFile() const
{
	const std::string statisticsTxt = this->options->out_prefix + ".statistics.txt";
	std::ofstream ofs(statisticsTxt.c_str());
	if (!ofs)
	{
		std::cerr << "[Error] Could not open (" << statisticsTxt << ")." << std::endl;
		std::exit(1);
	}
	ofs << "#K-mer\t" << this->options->kmer << std::endl;
	ofs << "#Pos\tSeq\tMutant\tWildType\tGval\tPval\tFDR\tBonferroni\n";

	// Calculate G-value for k-mer match analysis.
	this->gtest->kmer_match(this->mutantPosFreq, this->wildTypePosFreq);

	//========== Output ==========//
	for (size_t i = 0; i < this->mutantPosFreq.size(); i++)
	{
		ofs << i + 1 << "\t"
			<< this->vectorArray[i] << "\t"
			<< this->mutantPosFreq[i] << "\t"
			<< this->wildTypePosFreq[i] << "\t"
			<< (float)this->gtest->get_gval()[i] << "\t"
			<< (float)this->gtest->get_pval()[i] << "\t"
			<< (float)this->gtest->get_fdr()[i] << "\t"
			<< (float)this->gtest->get_bon()[i] << std::endl;
	}
	ofs.close();
}

/**
 * @brief Create the outside.txt file.
 *
 * @param mutantMerPair Mutant mer pairs at each end
 * @param wildTypeMerPair Wild type mer pairs at each end
 */
void StatisticsFile::create_outsideFile(
	const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &mutantMerPair,
	const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &wildTypeMerPair) const
{
	Complementary complementary;
	std::string fdr_str = std::to_string((float)this->options->threshold_fdr);
	while (fdr_str[fdr_str.length() - 1] == '0')
	{
		fdr_str.pop_back();
	}

	const std::string outsideFile = this->options->out_prefix + ".outside.txt";
	std::ofstream ofs(outsideFile.c_str());
	if (!ofs)
	{
		std::cerr << "[Error] Could not open (" << outsideFile << ")." << std::endl;
		std::exit(1);
	}

	auto [number_of_extensions, table_size, outsideData] = this->create_outsideData(mutantMerPair, wildTypeMerPair);

	// Calculate FDR using the Benjamini-Hochberg method.
	std::unordered_map<unsigned int, std::unordered_map<unsigned int, double>>
		fdr_extension = this->gtest->fdr_extension(outsideData.pval);

	ofs << "#K-mer\t"
		<< this->options->kmer
		<< "\tFDR\t"
		<< this->options->threshold_fdr
		<< "\tBases\t"
		<< this->options->bases_on_each_side
		<< std::endl;

	for (unsigned i = 0; i < this->vectorArray.length() - this->options->kmer; i++)
	{
		std::string kmer = this->vectorArray.substr(i, this->options->kmer);
		if (this->gtest->get_fdr()[i] <= this->options->threshold_fdr)
		{
			ofs << i + 1 << "\t"
				<< table_size[i] << "\t"
				<< kmer << "\t"
				<< this->mutantPosFreq[i] << "\t"
				<< this->wildTypePosFreq[i] << "\t"
				<< (float)this->gtest->get_gval()[i] << "\t"
				<< (float)this->gtest->get_pval()[i] << "\t"
				<< (float)this->gtest->get_fdr()[i] << "\t"
				<< (float)this->gtest->get_bon()[i]
				<< std::endl;

			for (size_t j = 0; j < outsideData.left_chain[i].size(); j++)
			{
				ofs << outsideData.left_chain[i][j] << "\t"
					<< outsideData.right_chain[i][j] << "\t"
					<< outsideData.mutant_count[i][j] << "\t"
					<< outsideData.wildType_count[i][j] << "\t"
					<< outsideData.left_chain[i][j]
					<< kmer
					<< outsideData.right_chain[i][j] << "\t"
					<< (float)outsideData.gval[i][j] << "\t"
					<< (float)outsideData.pval[i][j] << "\t"
					<< (float)fdr_extension[i][j] << "\t"
					<< (float)std::min(outsideData.pval[i][j] * number_of_extensions, 1.0)
					<< std::endl;
			}
		}
	}
	ofs.close();
}

//============================================================================//
// Private function
//============================================================================//
/**
 * @brief Create outside data.
 *
 * @param mutantMerPair Mutant mer pairs at each end
 * @param wildTypeMerPair Wild type mer pairs at each end
 * @return Outside data
 */
std::tuple<unsigned int, std::unordered_map<unsigned int, size_t>, OutsideData> StatisticsFile::create_outsideData(
	const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &mutantMerPair,
	const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &wildTypeMerPair) const
{
	// Obtain the complementary sequence of k-mer.
	Complementary complementary;
	// Outside data.
	OutsideData outsideData;
	// Number of outside data per k-mer
	std::unordered_map<unsigned int, size_t> table_size;

	std::map<std::pair<unsigned int, unsigned int>, double> gval_stock;
	std::map<std::pair<unsigned int, unsigned int>, double> pval_stock;
	unsigned int number_of_extensions = 0;

	for (size_t i = 0; i < this->vectorArray.length() - this->options->kmer; i++)
	{
		if (this->gtest->get_fdr()[i] <= this->options->threshold_fdr)
		{
			const std::string mer_plus = this->vectorPosPair.at(i).first;
			const std::string mer_minus = this->vectorPosPair.at(i).second;

			std::map<std::pair<std::string, std::string>, unsigned int> mutant_side_pair_count;
			std::map<std::pair<std::string, std::string>, unsigned int> wildType_side_pair_count;

			if (mutantMerPair.find(mer_plus) != mutantMerPair.end())
			{
				for (auto itr = mutantMerPair.at(mer_plus).begin(); itr != mutantMerPair.at(mer_plus).end(); ++itr)
				{
					mutant_side_pair_count[std::make_pair(itr->first, itr->second)]++;
				}
			}

			if (wildTypeMerPair.find(mer_plus) != wildTypeMerPair.end())
			{
				for (auto itr = wildTypeMerPair.at(mer_plus).begin(); itr != wildTypeMerPair.at(mer_plus).end(); ++itr)
				{
					wildType_side_pair_count[std::make_pair(itr->first, itr->second)]++;
				}
			}

			if (mer_plus != mer_minus)
			{
				if (mutantMerPair.find(mer_minus) != mutantMerPair.end())
				{
					for (auto itr = mutantMerPair.at(mer_minus).begin(); itr != mutantMerPair.at(mer_minus).end(); ++itr)
					{
						// Obtain the complementary sequence of k-mer.
						std::string revMer1 = complementary.mer(itr->second);
						std::string revMer2 = complementary.mer(itr->first);
						mutant_side_pair_count[std::make_pair(revMer1, revMer2)]++;
					}
				}

				if (wildTypeMerPair.find(mer_minus) != wildTypeMerPair.end())
				{
					for (auto itr = wildTypeMerPair.at(mer_minus).begin(); itr != wildTypeMerPair.at(mer_minus).end(); ++itr)
					{
						// Obtain the complementary sequence of k-mer.
						std::string revMer1 = complementary.mer(itr->second);
						std::string revMer2 = complementary.mer(itr->first);
						wildType_side_pair_count[std::make_pair(revMer1, revMer2)]++;
					}
				}
			}

			std::vector<std::pair<unsigned int, std::pair<std::string, std::string>>> v;
			for (auto itr = mutant_side_pair_count.begin(); itr != mutant_side_pair_count.end(); itr++)
			{
				v.push_back({itr->second, itr->first});
			}

			sort(v.rbegin(), v.rend());
			table_size[i] = v.size();
			unsigned int ic = 0;

			for (auto itr_v = v.begin(); itr_v != v.end(); ++itr_v)
			{
				outsideData.left_chain[i].push_back(itr_v->second.first);
				outsideData.right_chain[i].push_back(itr_v->second.second);

				unsigned int mutant_count =
					mutant_side_pair_count.find(itr_v->second) != mutant_side_pair_count.end()
						? mutant_side_pair_count[itr_v->second]
						: 0;
				unsigned int wildType_count =
					wildType_side_pair_count.find(itr_v->second) != wildType_side_pair_count.end()
						? wildType_side_pair_count[itr_v->second]
						: 0;

				outsideData.mutant_count[i].push_back(mutant_count);
				outsideData.wildType_count[i].push_back(wildType_count);

				std::pair<unsigned int, unsigned int> target = std::make_pair(mutant_count, wildType_count);
				if (gval_stock.find(target) != gval_stock.end())
				{
					outsideData.gval[i][ic] = gval_stock[target];
					outsideData.pval[i][ic] = pval_stock[target];
				}
				else
				{
					// G-test
					auto [g, p] = this->gtest->kmer_extension(mutant_count, wildType_count);
					outsideData.gval[i][ic] = g;
					outsideData.pval[i][ic] = p;

					gval_stock[target] = g;
					pval_stock[target] = p;
				}

				ic++;
				number_of_extensions++;
			}
		}
	}
	return {number_of_extensions, table_size, outsideData};
}
