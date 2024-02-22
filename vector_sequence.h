/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#ifndef VECTOR_SEQUENCE_H_
#define VECTOR_SEQUENCE_H_

#include <string>
#include <unordered_map>
#include "bitwise_operation.h"

/**
 * @brief Input vector sequences.
 */
class VectorSequence
{
public:
	/**
	 * @brief Construct a new Vector Sequence object
	 *
	 * @param options Execution options.
	 * @param bitwiseOperation Bitwise operation.
	 */
	VectorSequence(Options *options, BitwiseOperation *bitwiseOperation);

	/**
	 * @brief Destroy the Vector Sequence object
	 *
	 */
	virtual ~VectorSequence();

	/**
	 * @brief Read the fasta file.
	 *
	 * @param merCounter Counter of each mer
	 * @param posPair Position and k-mer complementary pair on vector
	 * @return Vector sequence
	 */
	std::string read_vectorFile(
		std::unordered_map<std::string, unsigned int> &merCounter,
		std::unordered_map<unsigned int, std::pair<std::string, std::string>> &posPair) const;

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
	 * @brief Set k-mer in hash table.
	 *
	 * @param sequence Vector sequence
	 * @param merCounter Counter of each mer
	 * @param posPair Position and k-mer complementary pair on vector
	 */
	void set_merCounter(
		std::string &sequence, std::unordered_map<std::string, unsigned int> &merCounter,
		std::unordered_map<unsigned int, std::pair<std::string, std::string>> &posPair) const;

	/**
	 * @brief Create chunk array.
	 *
	 * @param merCounter Counter of each mer
	 */
	void create_chunk(const std::unordered_map<std::string, unsigned int> &merCounter) const;
};
#endif /* VECTOR_SEQUENCE_H_ */
