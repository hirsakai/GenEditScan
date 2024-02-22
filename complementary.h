/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#ifndef COMPLEMENTARY_H_
#define COMPLEMENTARY_H_

#include <string>

/**
 * @brief Obtain the complementary sequence of k-mer.
 *
 */
class Complementary
{
public:
	/**
	 * @brief Construct a new Complementary object
	 *
	 */
	Complementary();

	/**
	 * @brief Destroy the Complementary object
	 *
	 */
	virtual ~Complementary();

	/**
	 * @brief Obtain the complementary sequence of k-mer.
	 *
	 * @param mer K-mer sequence
	 * @return Complementary sequence of k-mer
	 */
	std::string mer(const std::string &mer) const;
};
#endif /* COMPLEMENTARY_H_ */
