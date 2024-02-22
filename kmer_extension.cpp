/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#include "kmer_extension.h"
#include "complementary.h"

/**
 * @brief Construct a new Kmer Extension:: Kmer Extension object
 *
 * @param options Execution options.
 * @param bitwiseOperation Bitwise operation.
 * @param statisticsFile Create statistics files.
 */
KmerExtension::KmerExtension(Options *options, BitwiseOperation *bitwiseOperation,
							 StatisticsFile *statisticsFile)
{
	this->options = options;
	this->bitwiseOperation = bitwiseOperation;
	this->statisticsFile = statisticsFile;

	// FASTQ extension
	this->fastqExtension = new FastqExtension(this->options, this->bitwiseOperation);
}

/**
 * @brief Destroy the Kmer Extension:: Kmer Extension object
 *
 */
KmerExtension::~KmerExtension()
{
}

/**
 * @brief Execute extension analysis of k-mer.
 *
 */
void KmerExtension::execution() const
{
	std::cout << "\n---------- Extension analysis of k-mer (FDR <= "
			  << this->options->threshold_fdr << ") ----------" << std::endl;

	// Mutant mer pairs at each end
	std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> mutantMerCounter;

	// Set k-mer pairs
	if (this->set_merCounter(mutantMerCounter) == 0)
	{
		std::cout << "Count of target mer    = 0" << std::endl;
		return;
	}

	// Wild type mer pairs at each end
	std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> wildTypeMerCounter(mutantMerCounter);

	// Counter
	u_int64_t mutantMerTotalCounter = 0;
	u_int64_t wildTypeMerTotalCounter = 0;

	// Number of fastq files
	const size_t nMutant = this->options->mutant_files.size();

	// Mer pair
	std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> merPair;

	// Total mer counter
	u_int64_t merTotalCounter = 0;

#ifdef _OPENMP
#if _OPENMP < 202011
	omp_set_nested(1);
#endif
	omp_set_max_active_levels(2);
	omp_set_dynamic(0);
#pragma omp parallel for num_threads(this->options->outer_parallel) private(merTotalCounter, merPair) \
	reduction(+ : mutantMerTotalCounter, wildTypeMerTotalCounter)
#endif
	for (size_t i = 0; i < this->options->number_of_samples(); i++)
	{
		merTotalCounter = 0;
		if (i < nMutant)
		{
			// Read the fastq.gz file (mutant_files)
			merPair = this->fastqExtension->read_fastqFile(this->options->mutant_files[i],
														   mutantMerCounter, merTotalCounter);
			mutantMerTotalCounter += merTotalCounter;
#ifdef _OPENMP
#pragma omp critical(mutantPair)
#endif
			for (auto itr = merPair.begin(); itr != merPair.end(); ++itr)
			{
				for (auto itr_second = itr->second.begin();
					 itr_second != itr->second.end(); ++itr_second)
				{
					mutantMerCounter[itr->first].push_back(
						std::make_pair(itr_second->first, itr_second->second));
				}
			}
		}
		else
		{
			// Read the fastq.gz file (wildType_files)
			merPair = this->fastqExtension->read_fastqFile(this->options->wildType_files[i - nMutant],
														   wildTypeMerCounter, merTotalCounter);
			wildTypeMerTotalCounter += merTotalCounter;
#ifdef _OPENMP
#pragma omp critical(wildTypePair)
#endif
			for (auto itr = merPair.begin(); itr != merPair.end(); ++itr)
			{
				for (auto itr_second = itr->second.begin();
					 itr_second != itr->second.end(); ++itr_second)
				{
					wildTypeMerCounter[itr->first].push_back(
						std::make_pair(itr_second->first, itr_second->second));
				}
			}
		}
	}

	std::cout << "Count of mutant mer    = " << mutantMerTotalCounter << std::endl;
	std::cout << "Count of wild type mer = " << wildTypeMerTotalCounter << std::endl;

	if (this->options->threshold_fdr >= 0.0)
	{
		// Set mer total count.
		this->statisticsFile->set_merCounter(mutantMerTotalCounter, wildTypeMerTotalCounter);

		// Write the output.txt file.
		this->statisticsFile->create_outsideFile(mutantMerCounter, wildTypeMerCounter);
	}
}

//============================================================================//
// Private function
//============================================================================//
/**
 * @brief Set k-mer in hash table.
 *
 * @param mutantMerCounter Mutant mer pairs at each end
 * @return Count of target mer
 */
unsigned int KmerExtension::set_merCounter(
	std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &mutantMerCounter) const
{
	Complementary complementary;
	const std::string vectorArray = this->statisticsFile->get_vectorArray();
	const std::unordered_map<unsigned int, double> fdr = this->statisticsFile->get_fdr();

	// Set k-mer pairs
	std::vector<std::pair<std::string, std::string>> listPair;
	for (auto itr = fdr.begin(); itr != fdr.end(); ++itr)
	{
		if (itr->second <= this->options->threshold_fdr)
		{
			const std::string mer = vectorArray.substr(itr->first, this->options->kmer);
			// Obtain the complementary sequence of k-mer.
			const std::string revMer = complementary.mer(mer);
			mutantMerCounter[mer] = listPair;
			mutantMerCounter[revMer] = listPair;
		}
	}
	// Create chunk array.
	this->create_chunk(mutantMerCounter);

	return mutantMerCounter.size();
}

/**
 * @brief Create chunk array.
 *
 * @param merCounter Mer counter at each end
 */
void KmerExtension::create_chunk(
	const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &merCounter) const
{
	unsigned char *dna2bit = this->bitwiseOperation->get_dna2bit();
	unsigned char *chunk = this->bitwiseOperation->get_chunk();

	for (unsigned int i = 0; i < this->options->max_chunk_array; i++)
	{
		chunk[i] = 0;
	}

	for (auto itr = merCounter.begin(); itr != merCounter.end(); ++itr)
	{
		unsigned int dnabit = dna2bit[(unsigned char)itr->first[0]];
		for (unsigned int i = 1; i < this->options->chunk_length; i++)
		{
			dnabit = (dnabit << 2) + dna2bit[(unsigned char)itr->first[i]];
		}
		if (dnabit != this->options->max_chunk_array)
		{
			chunk[dnabit] = 1;
		}
	}
}
