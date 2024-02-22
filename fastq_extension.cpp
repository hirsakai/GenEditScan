/*
 * GenEditScan
 * Copyright 2018 National Agriculture and Food Research Organization (NARO)
 */
#include <iostream>
#include <zlib.h>
#include "fastq_extension.h"

/**
 * @brief Construct a new Fastq Extension:: Fastq Extension object
 *
 * @param options Execution options.
 * @param bitwiseOperation Bitwise operation.
 */
FastqExtension::FastqExtension(Options *options, BitwiseOperation *bitwiseOperation)
{
    this->options = options;
    this->bitwiseOperation = bitwiseOperation;
}

/**
 * @brief Destroy the Fastq Extension:: Fastq Extension object
 *
 */
FastqExtension::~FastqExtension()
{
}

/**
 * @brief Read the fastq.gz file.
 *
 * @param fastqFile FASTQ file
 * @param merCounter Mer counter at each end
 * @param merTotalCounter Mer total counter per file
 * @return Mer pairs at each end for parallel processing
 */
std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> FastqExtension::read_fastqFile(
    const std::string &fastqFile,
    const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &merCounter,
    u_int64_t &merTotalCounter) const
{
    // File mode
    const gzFile file = gzopen(fastqFile.c_str(), "rb");
    if (!file)
    {
        std::cerr << "[Error] Could not open (" << fastqFile << ")." << std::endl;
        std::exit(1);
    }

#ifdef _OPENMP
#pragma omp single nowait
#endif
    std::cout << "Count of target mer    = " << merCounter.size() << std::endl;

    const unsigned int kmer = this->options->kmer;
    const unsigned int nbase = this->options->bases_on_each_side;
    const unsigned int max_buff = this->options->max_read_length + 2;
    char buff[max_buff];
    std::string aLine[4];
    unsigned int nLine = 0;
    u_int64_t readCounter = 0;
    std::vector<std::string> fastqData;
    std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> merLocalPair;

    while (gzgets(file, buff, max_buff) != Z_NULL)
    {
        aLine[nLine++] = std::string(buff);
        if (nLine == 4)
        {
            nLine = 0;
            if (aLine[1].length() > kmer + nbase * 2)
            {
                // Delete line break (\n)
                aLine[1].pop_back();
                fastqData.push_back(aLine[1]);
                if (fastqData.size() > this->options->fastq_read_lines)
                {
                    this->count_extension(fastqFile, fastqData, merCounter,
                                          merLocalPair, merTotalCounter, readCounter);
                    fastqData.clear();
                }
            }
        }
    }

    this->count_extension(fastqFile, fastqData, merCounter, merLocalPair,
                          merTotalCounter, readCounter);
    fastqData.clear();
    gzclose(file);
    return merLocalPair;
}

//============================================================================//
// Private function
//============================================================================//
/**
 * @brief Count k-mer.
 *
 * @param fastqFile FASTQ file
 * @param fastqData FASTQ data
 * @param merCounter Mer counter at each end
 * @param merLocalPair Mer pairs at each end for parallel processing
 * @param merTotalCounter Mer total counter per file
 * @param readCounter Read counter per file
 */
void FastqExtension::count_extension(
    const std::string &fastqFile, std::vector<std::string> &fastqData,
    const std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &merCounter,
    std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> &merLocalPair,
    u_int64_t &merTotalCounter, u_int64_t &readCounter) const
{
    const unsigned int kmer = this->options->kmer;
    const unsigned int nbase = this->options->bases_on_each_side;
    const unsigned int mask = this->options->max_chunk_array;
    const unsigned int chunk_length = this->options->chunk_length;
    unsigned char *dna2bit = this->bitwiseOperation->get_dna2bit();
    unsigned char *chunk = this->bitwiseOperation->get_chunk();
    unsigned int dnabit, j;
    std::string mer, p5, p3;

#ifdef _OPENMP
#pragma omp parallel for num_threads(this->options->inner_parallel) private(dnabit, j, mer, p5, p3) \
    reduction(+ : merTotalCounter)
#endif
    for (size_t i = 0; i < fastqData.size(); i++)
    {
#ifdef _OPENMP
#pragma omp critical(extension)
#endif
        if (++readCounter % this->options->log_output_interval == 0)
        {
            std::cerr << fastqFile << ": parsing " << readCounter << " reads (k-mer extension)." << std::endl;
        }

        dnabit = dna2bit[(unsigned char)fastqData[i][0]];
        for (j = 1; j < chunk_length - 1; j++)
        {
            dnabit = (dnabit << 2) + dna2bit[(unsigned char)fastqData[i][j]];
        }

        for (j = 0; j < nbase; j++)
        {
            dnabit = (dnabit << 2) + dna2bit[(unsigned char)fastqData[i][chunk_length - 1 + j]];
        }

        for (j = nbase; j <= fastqData[i].length() - kmer - nbase; j++)
        {
            dnabit = (dnabit << 2) + dna2bit[(unsigned char)fastqData[i][chunk_length - 1 + j]];
            dnabit = dnabit & mask;
            if (chunk[dnabit] == 1 || dnabit == mask)
            {
                mer.assign(fastqData[i], j, kmer);
                if (merCounter.find(mer) != merCounter.end())
                {
                    p5 = fastqData[i].substr(j - nbase, nbase);
                    p3 = fastqData[i].substr(j + kmer, nbase);
#ifdef _OPENMP
#pragma omp critical(push)
#endif
                    merLocalPair[mer].push_back(std::make_pair(p5, p3));
                }
            }
            merTotalCounter++;
        }
    }
}
