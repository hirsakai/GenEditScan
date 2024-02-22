/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#ifndef OUTSIDE_DATA_H_
#define OUTSIDE_DATA_H_

#include <unordered_map>

/**
 * @brief Outside data.
 *
 */
struct OutsideData
{
	std::unordered_map<unsigned int, std::unordered_map<unsigned int, double>> gval;
	std::unordered_map<unsigned int, std::unordered_map<unsigned int, double>> pval;
	std::unordered_map<unsigned int, std::vector<std::string>> left_chain;
	std::unordered_map<unsigned int, std::vector<std::string>> right_chain;
	std::unordered_map<unsigned int, std::vector<unsigned int>> mutant_count;
	std::unordered_map<unsigned int, std::vector<unsigned int>> wildType_count;
};
#endif /* OUTSIDE_DATA_H_ */
