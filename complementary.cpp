/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#include <algorithm>
#include "complementary.h"

/**
 * @brief Construct a new Complementary:: Complementary object
 *
 */
Complementary::Complementary()
{
}

/**
 * @brief Destroy the Complementary:: Complementary object
 *
 */
Complementary::~Complementary()
{
}

/**
 * @brief Obtain the complementary sequence of k-mer.
 *
 * @param mer K-mer sequence
 * @return Complementary sequence of k-mer
 */
std::string Complementary::mer(const std::string &mer) const
{
	std::string revMer(mer);
	std::reverse(revMer.begin(), revMer.end());
	for (size_t i = 0; i < mer.length(); ++i)
	{
		switch (revMer[i])
		{
		case 'A':
			revMer[i] = 'T';
			break;
		case 'T':
			revMer[i] = 'A';
			break;
		case 'G':
			revMer[i] = 'C';
			break;
		case 'C':
			revMer[i] = 'G';
			break;
		default:
			break;
		}
	}
	return revMer;
}
