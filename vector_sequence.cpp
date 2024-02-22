/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#include <fstream>
#include <algorithm>
#include "vector_sequence.h"
#include "complementary.h"

/**
 * @brief Construct a new Vector Sequence:: Vector Sequence object
 *
 * @param options Execution options.
 * @param bitwiseOperation Bitwise operation.
 */
VectorSequence::VectorSequence(Options *options, BitwiseOperation *bitwiseOperation)
{
	this->options = options;
	this->bitwiseOperation = bitwiseOperation;
}

/**
 * @brief Destroy the Vector Sequence:: Vector Sequence object
 *
 */
VectorSequence::~VectorSequence()
{
}

/**
 * @brief Read the fasta file.
 *
 * @param merCounter Counter of each mer
 * @param posPair Position and k-mer complementary pair on vector
 * @return Vector sequence
 */
std::string VectorSequence::read_vectorFile(
	std::unordered_map<std::string, unsigned int> &merCounter,
	std::unordered_map<unsigned int, std::pair<std::string, std::string>> &posPair) const
{
	std::ifstream ifs(this->options->vector_file.c_str());
	if (!ifs)
	{
		std::cerr << "[Error] Could not open vector file("
				  << this->options->vector_file << ")." << std::endl;
		std::exit(1);
	}

	std::string sequence;
	std::string str;

	while (getline(ifs, str))
	{
		if (str[0] == '>')
		{
			if (sequence.length() > 0)
			{
				break;
			}
		}
		else
		{
			// Delete line breaks
			if (str[str.length() - 1] == '\n')
			{
				str.pop_back();
			}
			if (str[str.length() - 1] == '\r')
			{
				str.pop_back();
			}
			sequence += str;
		}
	}
	// Set k-mer in hash table.
	this->set_merCounter(sequence, merCounter, posPair);

	// Create chunk array.
	this->create_chunk(merCounter);
	return sequence;
}

//============================================================================//
// Private function
//============================================================================//
/**
 * @brief Set k-mer in hash table.
 *
 * @param sequence Vector sequence
 * @param merCounter Counter of each mer
 * @param posPair Position and k-mer complementary pair on vector
 */
void VectorSequence::set_merCounter(
	std::string &sequence,
	std::unordered_map<std::string, unsigned int> &merCounter,
	std::unordered_map<unsigned int, std::pair<std::string, std::string>> &posPair) const
{
	Complementary complementary;
	const unsigned int kmer = this->options->kmer;
	const unsigned int vector_length = sequence.length();

	if (vector_length >= kmer)
	{
		sequence += sequence.substr(0, kmer - 1);
		// Convert to upper case
		transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);
		for (unsigned int i = 0; i < vector_length; i++)
		{
			std::string mer = sequence.substr(i, kmer);
			// Obtain the complementary sequence of k-mer.
			std::string revMer = complementary.mer(mer);
			merCounter[mer] = 0;
			merCounter[revMer] = 0;
			posPair[i] = std::make_pair(mer, revMer);
		}
	}
	else
	{
		std::cerr << "[Error] Vector is shorter than k-mer.\n";
		std::exit(1);
	}
}

/**
 * @brief Create chunk array.
 *
 * @param merCounter Counter of each mer
 */
void VectorSequence::create_chunk(
	const std::unordered_map<std::string, unsigned int> &merCounter) const
{
	unsigned char *dna2bit = this->bitwiseOperation->get_dna2bit();
	unsigned char *chunk = this->bitwiseOperation->get_chunk();

	for (unsigned i = 0; i < this->options->max_chunk_array; i++)
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
