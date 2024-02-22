/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#include "bitwise_operation.h"

/**
 * @brief Construct a new Bitwise Operation:: Bitwise Operation object
 *
 */
BitwiseOperation::BitwiseOperation(Options *options)
{
	this->dna2bit = new unsigned char[128];
	for (unsigned int i = 0; i < 128; i++)
	{
		this->dna2bit[i] = 0;
	}
	this->dna2bit[84] = 0; // T or others
	this->dna2bit[67] = 1; // C
	this->dna2bit[65] = 2; // A
	this->dna2bit[71] = 3; // G
	this->chunk = new unsigned char[options->max_chunk_array];
}

/**
 * @brief Destroy the BitwiseOperation:: BitwiseOperation object
 *
 */
BitwiseOperation::~BitwiseOperation()
{
	delete[] this->dna2bit;
	delete[] this->chunk;
}
