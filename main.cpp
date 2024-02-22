//============================================================================//
// Name        : GenEditScan
// Author      : NARO
// Version     : 202402
// Copyright   : (c) 2018 National Agriculture and Food Research Organization (NARO)
// Description : K-mer analysis tool
//============================================================================//
#include <getopt.h>
#include <math.h>
#include <cstring>
#include <sstream>
#include "bitwise_operation.h"
#include "statistics_file.h"
#include "kmer_match.h"
#include "kmer_extension.h"

/**
 * @brief Split the delimiter separator.
 *
 * @param fastq_files Fastq files connected by delimiters
 * @param delimiter delimiter (comma)
 * @return Fastq files
 */
std::vector<std::string> split(const std::string fastq_files, const char delimiter)
{
	std::vector<std::string> fastq_vector;
	std::stringstream ss(fastq_files);
	while (ss.good())
	{
		std::string substr;
		getline(ss, substr, delimiter);
		fastq_vector.push_back(substr);
	}
	return fastq_vector;
}

/**
 * @brief Print help menu.
 *
 * @param options Execution options.
 * @param version Program version
 * @param execute Executable program name
 */
void help(const Options &options, const std::string &version, const std::string &execute)
{
	std::cerr << version << std::endl;
	std::cerr << "Usage : " << execute << " kmer [options]\n";
	std::cerr << "\n[required]\n";
	std::cerr << "-v | --vector   : Vector file\n";
	std::cerr << "-m | --mutant   : Mutant files (connect with comma)\n";
	std::cerr << "-w | --wild     : Wild type files (connect with comma)\n";
	std::cerr << "\n[optional]\n";
	std::cerr << "-k | --kmer     : K-mer (" << options.kmer << ")\n";
	std::cerr << "-f | --fdr      : Threshold by FDR (" << options.threshold_fdr << ")\n";
	std::cerr << "-b | --bases    : Number of bases on each side (" << options.bases_on_each_side << ")\n";
	std::cerr << "-o | --out      : Output prefix (" << options.out_prefix << ")\n";
	std::cerr << "-t | --threads  : Number of threads (all threads)\n";
	std::cerr << "-l | --length   : Maximum read length (" << options.max_read_length << ")\n";
	std::cerr << "-r | --read     : Number of lines of Fastq file to be read in memory (" << options.fastq_read_lines << ")\n";
	std::cerr << "-i | --interval : Log output interval (" << options.log_output_interval << ")\n";
	std::cerr << "-h | --help     : Print this menu\n";
}

/**
 * @brief Main function.
 *
 * @param argc Number of arguments
 * @param argv Arguments
 * @return Exit code
 */
int main(int argc, char *argv[])
{
	/**
	 * Execution options.
	 */
	Options options;

	const char DELIMITER = ',';

	const std::string version = "Program version: GenEditScan-1.0.0";

	const struct option long_options[] = {
		{"vector", required_argument, NULL, 'v'},
		{"mutant", required_argument, NULL, 'm'},
		{"wild", required_argument, NULL, 'w'},
		{"kmer", required_argument, NULL, 'k'},
		{"fdr", required_argument, NULL, 'f'},
		{"bases", required_argument, NULL, 'b'},
		{"out", required_argument, NULL, 'o'},
		{"threads", required_argument, NULL, 't'},
		{"read", required_argument, NULL, 'r'},
		{"length", required_argument, NULL, 'l'},
		{"interval", required_argument, NULL, 'i'},
		{"help", required_argument, NULL, 'h'},
		{0, 0, 0, 0}};

	try
	{
		int c;
		int long_index;
		unsigned int kmer;
		while ((c = getopt_long(argc, argv, "v:m:w:k:f:b:o:t:r:l:i:h::", long_options, &long_index)) != -1)
		{
			switch (c)
			{
			case 'v':
				options.vector_file = optarg;
				break;
			case 'm':
				options.mutant_files = split(optarg, DELIMITER);
				break;
			case 'w':
				options.wildType_files = split(optarg, DELIMITER);
				break;
			case 'k':
				kmer = std::stoi(optarg);
				if (kmer >= options.MIN_CHUNKLENGTH)
				{
					options.kmer = kmer;
				}
				else
				{
					std::cerr << "[Error] K-mer (" << kmer << ") must be >= " << options.MIN_CHUNKLENGTH << "." << std::endl;
					return EXIT_FAILURE;
				}
				break;
			case 'f':
				options.threshold_fdr = std::stod(optarg);
				break;
			case 'b':
				options.bases_on_each_side = std::stoi(optarg);
				break;
			case 'o':
				options.out_prefix = optarg;
				break;
			case 't':
				options.threads = std::atoi(optarg);
				break;
			case 'r':
				options.fastq_read_lines = std::stoi(optarg);
				break;
			case 'l':
				options.max_read_length = std::stoi(optarg);
				break;
			case 'i':
				options.log_output_interval = std::stoi(optarg);
				break;
			case 'h':
				help(options, version, argv[0]);
				return EXIT_FAILURE;
			default:
				std::cerr << "[Error] Unrecognized option " << c << " " << optopt << std::endl;
				return EXIT_FAILURE;
			}
		}

		if (optind == argc || strcmp(argv[optind], "kmer") != 0 || options.vector_file.length() == 0
		 || options.mutant_files.size() == 0 || options.wildType_files.size() == 0)
		{
			help(options, version, argv[0]);
			return EXIT_FAILURE;
		}

		options.calc_mode = argv[optind];
		options.chunk_length = std::min(options.kmer, options.MAX_CHUNKLENGTH);
		options.max_chunk_array = (unsigned int)(pow(2.0, 2.0 * options.chunk_length) - 1.0);
		options.output(version);
	}
	catch (const std::exception &e)
	{
		std::cerr << "[Error] " << argv[0] << ": " << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	/**
	 * Execute the k-mer analysis
	 */
	try
	{
		/**
		 * Bitwise operation.
		 */
		BitwiseOperation *bitwiseOperation = new BitwiseOperation(&options);

		/**
		 * Create statistics files.
		 */
		StatisticsFile *statisticsFile = new StatisticsFile(&options);

		/**
		 * K-mer match analysis
		 */
		KmerMatch *kmerMatch = new KmerMatch(&options, bitwiseOperation, statisticsFile);
		kmerMatch->execution();
		delete kmerMatch;

		/**
		 * K-mer extension analysis
		 */
		KmerExtension *kmerExtension = new KmerExtension(&options, bitwiseOperation,
														 statisticsFile);
		kmerExtension->execution();
		delete kmerExtension;

		delete statisticsFile;
		delete bitwiseOperation;

		std::cout << "\nEnd time    : " << options.get_now() << std::endl;
		std::cout << "Elapsed time: " << options.get_elapsed() << std::endl;
		return EXIT_SUCCESS;
	}
	catch (const std::bad_alloc &e)
	{
		std::cerr << "[Error] " << argv[0] << ": out of memory." << std::endl;
		return EXIT_FAILURE;
	}
	catch (const std::exception &e)
	{
		const char *s = e.what();
		if (*s)
		{
			std::cerr << "[Error] " << argv[0] << ": " << s << std::endl;
		}
		return EXIT_FAILURE;
	}
}
