/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#ifndef BITWISE_OPERATION_H_
#define BITWISE_OPERATION_H_

#include "options.h"

/**
 * @brief Bitwise operation.
 *
 */
class BitwiseOperation
{
public:
	/**
	 * @brief Construct a new Bitwise Operation object
	 *
	 */
	BitwiseOperation(Options *options);

	/**
	 * @brief Destroy the Bitwise Operation object
	 *
	 */
	virtual ~BitwiseOperation();

	// Getter

	unsigned char *get_dna2bit() const
	{
		return this->dna2bit;
	}

	unsigned char *get_chunk() const
	{
		return this->chunk;
	}

private:
	/**
	 * @brief For bitwise operation. DNA expressed in 2 bits.
	 *
	 */
	unsigned char *dna2bit;

	/**
	 * @brief For bitwise operation. Chunk array.
	 *
	 */
	unsigned char *chunk;
};
#endif /* BITWISE_OPERATION_H_ */
