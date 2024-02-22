/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#ifndef OPTIONS_H_
#define OPTIONS_H_

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * Execution options.
 */
class Options
{
public:
	/**
	 * @brief Construct a new Options object
	 */
	Options()
	{
		this->start_time = std::chrono::system_clock::now();
	}

	/**
	 * @brief Destroy the Options object
	 */
	virtual ~Options()
	{
	}

	// Calculation mode (kmer or oft)
	std::string calc_mode;

	// Vector file
	std::string vector_file;

	// Mutant files
	std::vector<std::string> mutant_files;

	// Wild type files
	std::vector<std::string> wildType_files;

	// K-mer
	unsigned int kmer = 20;

	// Threshold by FDR
	double threshold_fdr = 0.01;

	// Number of bases on each side
	unsigned int bases_on_each_side = 5;

	// Output prefix
	std::string out_prefix = "out_prefix";

	// Maximum read length
	unsigned int max_read_length = 512;

	// Number of lines of Fastq file to be read in memory
	unsigned int fastq_read_lines = 10000000;

	// Log output interval
	unsigned int log_output_interval = 1000000;

	// Number of threads
	unsigned int threads = 0;

	// OpenMP outer parallel
	unsigned int outer_parallel = 2;

	// OpenMP inner parallel
	unsigned int inner_parallel = 1;

	// int(32 bit) / (2 bit/base) = 16 bases
	const unsigned int MAX_CHUNKLENGTH = 16;

	// MIN_CHUNKLENGTH > 0
	const unsigned int MIN_CHUNKLENGTH = 8;

	// Chunk length
	unsigned int chunk_length;

	// Array length required for specified chunk length
	unsigned int max_chunk_array;

	// start time
	std::chrono::system_clock::time_point start_time;

	/**
	 * @brief Get start time.
	 *
	 * @return Start time
	 */
	std::string get_start()
	{
		const std::time_t t = std::chrono::system_clock::to_time_t(this->start_time);
		return static_cast<const std::ostringstream &>(std::ostringstream() << std::put_time(localtime(&t), "%F %X")).str();
	}

	/**
	 * @brief Get the current time.
	 *
	 * @return Current time
	 */
	std::string get_now()
	{
		auto now = std::chrono::system_clock::now();
		const std::time_t t = std::chrono::system_clock::to_time_t(now);
		return static_cast<const std::ostringstream &>(std::ostringstream() << std::put_time(localtime(&t), "%F %X")).str();
	}

	/**
	 * @brief Get end time.
	 *
	 * @return End time
	 */
	std::string get_elapsed()
	{
		auto elapsed = std::chrono::system_clock::now() - this->start_time;
		auto seconds_only = std::chrono::duration_cast<std::chrono::seconds>(elapsed);
		auto hours = std::chrono::duration_cast<std::chrono::hours>(elapsed);
		elapsed -= hours;
		auto minutes = std::chrono::duration_cast<std::chrono::minutes>(elapsed);
		elapsed -= minutes;
		auto seconds = std::chrono::duration_cast<std::chrono::seconds>(elapsed);

		std::ostringstream ostr;

		if (hours.count() > 1)
		{
			ostr << hours.count() << " hours ";
		}
		else if (hours.count() == 1)
		{
			ostr << hours.count() << " hour ";
		}

		if (minutes.count() > 1)
		{
			ostr << minutes.count() << " minutes ";
		}
		else if (minutes.count() == 1)
		{
			ostr << minutes.count() << " minute ";
		}

		if (seconds.count() > 1)
		{
			ostr << seconds.count() << " seconds ";
		}
		else if (seconds.count() == 1)
		{
			ostr << seconds.count() << " second ";
		}

		ostr << "(" << seconds_only.count() << " seconds)";
		return static_cast<const std::ostringstream &>(ostr).str();
	}

	// Number of mutant and wild type files
	unsigned int number_of_samples() const
	{
		return this->mutant_files.size() + this->wildType_files.size();
	}

	/**
	 * @brief Echo the settings.
	 *
	 * @param version Program version
	 */
	void output(const std::string &version)
	{
		std::cout << version << std::endl;
		std::cout << "Start time     : " << this->get_start() << std::endl;
		std::cout << "\n---------- K-mer analysis settings ----------" << std::endl;
		std::cout << "Vector file = " << this->vector_file << std::endl;
		std::cout << "Mutant files:\n";
		for (auto itr = mutant_files.begin(); itr != mutant_files.end(); ++itr)
		{
			std::cout << "              " << *itr << std::endl;
		}
		std::cout << "Wild type files:\n";
		for (auto itr = wildType_files.begin(); itr != wildType_files.end(); ++itr)
		{
			std::cout << "              " << *itr << std::endl;
		}
		std::cout << "K-mer                         = " << this->kmer << std::endl;
		std::cout << "Threshold by FDR              = " << this->threshold_fdr << std::endl;
		std::cout << "Number of bases on each side  = " << this->bases_on_each_side << std::endl;
		std::cout << "Output prefix                 = " << this->out_prefix << std::endl;
		std::cout << "Maximum read length           = " << this->max_read_length << std::endl;
		std::cout << "Number of lines of Fastq file" << std::endl;
		std::cout << "         to be read in memory = " << this->fastq_read_lines << std::endl;
		std::cout << "Log output interval           = " << this->log_output_interval << std::endl;
		std::cout << std::flush;

#ifdef _OPENMP
		unsigned int num_threads;
		unsigned int inner_threads;
#pragma omp parallel
		num_threads = this->threads == 0 ? omp_get_num_threads() : this->threads;
		this->outer_parallel = std::min(this->number_of_samples(), num_threads);
		inner_threads = num_threads / this->outer_parallel;
		this->inner_parallel = std::max(inner_threads, 1u);
		std::cout << "Number of threads             = " << num_threads << std::endl;
		std::cout << "OpenMP outer parallel         = " << this->outer_parallel << std::endl;
		std::cout << "OpenMP inner parallel         = " << this->inner_parallel << std::endl;
		std::cout << std::flush;
#endif
	}
};
#endif /* OPTIONS_H_ */
