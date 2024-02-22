/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#include <fstream>
#include "kmer_match.h"
#include "vector_sequence.h"

/**
 * @brief Construct a new Kmer Match:: Kmer Match object
 *
 * @param options Execution options.
 * @param bitwiseOperation Bitwise operation.
 * @param statisticsFile Create statistics files.
 */
KmerMatch::KmerMatch(Options *options, BitwiseOperation *bitwiseOperation,
					 StatisticsFile *statisticsFile)
{
	this->options = options;
	this->bitwiseOperation = bitwiseOperation;
	this->statisticsFile = statisticsFile;

	// FASTQ match
	this->fastqMatch = new FastqMatch(this->options, this->bitwiseOperation);
}

/**
 * @brief Destroy the Kmer Match:: Kmer Match object
 *
 */
KmerMatch::~KmerMatch()
{
}

/**
 * @brief Execute match analysis of k-mer.
 *
 */
void KmerMatch::execution() const
{
	std::cout << "\n---------- Match analysis of k-mer (K-mer = "
			  << this->options->kmer << ") ----------" << std::endl;

	// Position and k-mer complementary pair on vector
	std::unordered_map<unsigned int, std::pair<std::string, std::string>> vectorPosPair;

	// Mutant mer counter
	std::unordered_map<std::string, unsigned int> mutantMerCounter;

	// Read the vector file
	VectorSequence *vectorSequence = new VectorSequence(this->options, this->bitwiseOperation);
	const std::string vectorArray = vectorSequence->read_vectorFile(mutantMerCounter, vectorPosPair);

	this->statisticsFile->set_vectorArray(vectorArray);
	this->statisticsFile->set_vectorPosPair(vectorPosPair);

	// Wild type mer counter
	std::unordered_map<std::string, unsigned int> wildTypeMerCounter(mutantMerCounter);

	// Counter
	u_int64_t mutantMerTotalCounter = 0;
	u_int64_t wildTypeMerTotalCounter = 0;

	// Number of fastq files
	const size_t nMutant = this->options->mutant_files.size();

	// Mer counter
	std::unordered_map<std::string, unsigned int> merCounter;

	// Total mer counter
	u_int64_t merTotalCounter = 0;

#ifdef _OPENMP
#if _OPENMP < 202011
	omp_set_nested(1);
#endif
	omp_set_max_active_levels(2);
	omp_set_dynamic(0);
#pragma omp parallel for num_threads(this->options->outer_parallel) private(merTotalCounter, merCounter) \
	reduction(+ : mutantMerTotalCounter, wildTypeMerTotalCounter)
#endif
	for (size_t i = 0; i < this->options->number_of_samples(); i++)
	{
		merTotalCounter = 0;
		if (i < nMutant)
		{
			// Read the fastq.gz file (mutant_files)
			merCounter = this->fastqMatch->read_fastqFile(this->options->mutant_files[i],
														  mutantMerCounter, merTotalCounter);
			mutantMerTotalCounter += merTotalCounter;
#ifdef _OPENMP
#pragma omp critical(mutant)
#endif
			for (auto itr = mutantMerCounter.begin(); itr != mutantMerCounter.end(); ++itr)
			{
				mutantMerCounter[itr->first] += merCounter[itr->first];
			}
		}
		else
		{
			// Read the fastq.gz file (wildType_files)
			merCounter = this->fastqMatch->read_fastqFile(this->options->wildType_files[i - nMutant],
														  wildTypeMerCounter, merTotalCounter);
			wildTypeMerTotalCounter += merTotalCounter;
#ifdef _OPENMP
#pragma omp critical(wildType)
#endif
			for (auto itr = wildTypeMerCounter.begin(); itr != wildTypeMerCounter.end(); ++itr)
			{
				wildTypeMerCounter[itr->first] += merCounter[itr->first];
			}
		}
	}

	std::cout << "Count of mutant mer    = " << mutantMerTotalCounter << std::endl;
	std::cout << "Count of wild type mer = " << wildTypeMerTotalCounter << std::endl;

	this->control_freqFile(mutantMerCounter, wildTypeMerCounter);

	// Set mer total count.
	this->statisticsFile->set_merCounter(mutantMerTotalCounter, wildTypeMerTotalCounter);

	// Write the statistics.txt file.
	this->statisticsFile->create_statisticsFile();
}

//============================================================================//
// Private function
//============================================================================//
/**
 * @brief Set position frequencies and write merFreq.txt files.
 *
 * @param mutantMerCounter Mutant mer counter
 * @param wildTypeMerCounter Wild type mer counter
 */
void KmerMatch::control_freqFile(
	const std::unordered_map<std::string, unsigned int> &mutantMerCounter,
	const std::unordered_map<std::string, unsigned int> &wildTypeMerCounter) const
{
//========== Output ==========//
#ifdef _OPENMP
#pragma omp parallel sections
	{
#pragma omp section
		{
#endif
			std::map<unsigned int, std::pair<std::string, std::string>> vectorPosPair(
				this->statisticsFile->get_vectorPosPair().begin(),
				this->statisticsFile->get_vectorPosPair().end());

			// Set position frequencies.
			const std::vector<unsigned int> mutantPosFreq = this->set_posFreq(mutantMerCounter, vectorPosPair);
			const std::vector<unsigned int> wildTypePosFreq = this->set_posFreq(wildTypeMerCounter, vectorPosPair);
			this->statisticsFile->set_mutantPosFreq(mutantPosFreq);
			this->statisticsFile->set_wildTypePosFreq(wildTypePosFreq);

#ifdef _OPENMP
		}
#pragma omp section
		{
#endif
			// Create merFreq.txt file.
			this->create_merFreqFile(mutantMerCounter, ".mutant");
			this->create_merFreqFile(wildTypeMerCounter, ".wildtype");
#ifdef _OPENMP
		}
	}
#endif
}

/**
 * @brief Set position frequencies.
 *
 * @param merCounter Counter of each mer
 * @param vectorPosPair Position and k-mer complementary pair on vector
 * @return Position frequencies
 */
std::vector<unsigned int> KmerMatch::set_posFreq(
	const std::unordered_map<std::string, unsigned int> &merCounter,
	std::map<unsigned int, std::pair<std::string, std::string>> &vectorPosPair) const
{
	std::vector<unsigned int> posBothFreq;
	for (std::map<unsigned int, std::pair<std::string, std::string>>::iterator itr = vectorPosPair.begin();
		 itr != vectorPosPair.end(); ++itr)
	{
		const unsigned int merBothCounter = merCounter.at(itr->second.first) +
										 merCounter.at(itr->second.second);
		posBothFreq.push_back(merBothCounter);
	}
	return posBothFreq;
}

/**
 * @brief Create merFreq.txt file.
 *
 * @param merCounter Counter of each mer
 * @param type '_mutant' or '_wildtype'
 */
void KmerMatch::create_merFreqFile(const std::unordered_map<std::string, unsigned int> &merCounter,
								   const std::string type) const
{
	const std::string outfile = this->options->out_prefix + type + ".merFreq.txt";
	std::ofstream ofs(outfile.c_str());
	if (!ofs)
	{
		std::cerr << "[Error] Could not open (" << outfile << ")." << std::endl;
		std::exit(1);
	}

	std::map<std::string, unsigned int> sortedCount(merCounter.begin(), merCounter.end());

	for (std::map<std::string, unsigned int>::iterator itr = sortedCount.begin();
		 itr != sortedCount.end(); ++itr)
	{
		ofs << itr->first << '\t' << itr->second << std::endl;
	}
	ofs.close();
}
