#include <fstream>
#include <getopt.h>
#include <iostream>
#include <map>
#include <omp.h>
#include <string>
#include <unordered_map>
#include <vector>

#include <btllib/seq.hpp>
#include <btllib/seq_reader.hpp>
#include <btllib/mi_bloom_filter.hpp>
#include <Sequence/Translate.hpp>

#include "AAHash.hpp"

size_t calc_optimal_size(size_t entries, unsigned hash_num, double occupancy)
{
    size_t non64ApproxVal =
      size_t(-double(entries) * double(hash_num) / log(1.0 - occupancy));
    return non64ApproxVal + (64 - non64ApproxVal % 64);
}


std::vector<std::string> sixframe_translate(const std::string& dna) {
    std::vector<std::string> protein;
    std::string rev_dna = btllib::get_reverse_complement(dna);
    protein.push_back(Sequence::Translate(dna.begin(), dna.end()));
    protein.push_back(Sequence::Translate(rev_dna.begin(), rev_dna.end()));
    protein.push_back(Sequence::Translate(dna.begin() + 1, dna.end()));
    protein.push_back(Sequence::Translate(rev_dna.begin() + 1, rev_dna.end()));
    protein.push_back(Sequence::Translate(dna.begin() + 2, dna.end()));
    protein.push_back(Sequence::Translate(rev_dna.begin() + 2, rev_dna.end()));
    return protein;
}

// main function that takes in command line arguments using getopt
int main(int argc, char** argv) {
    // declare variables
    int opt;
    int option_index = 0;
    int verbose_flag = 0;
    int help_flag = 0;
    size_t threads = 1;
    std::string input_file = "";
    std::string reference_path = "";
    std::string output_prefix = "_";
    uint8_t hash_num = 1;
    uint8_t kmer_size = 10;
    uint64_t genome_size = 0;
    static struct option long_options[] = {
        {"help", no_argument, &help_flag, 1},
        {"verbose", no_argument, &verbose_flag, 1},
        {"threads", required_argument, 0, 't'},
        {"input", required_argument, 0, 'i'},
        {"reference", required_argument, 0, 'r'},
        {"output", required_argument, 0, 'o'},
        {"hash", required_argument, 0, 'h'},
        {"kmer", required_argument, 0, 'k'},
        {"genome", required_argument, 0, 'g'},
        {0, 0, 0, 0}
    };


    // loop through command line arguments
    while ((opt = getopt_long(argc, argv, "g:h:i:t:o:r:k:", long_options, &option_index)) != -1) {
        switch (opt) {
            case 0:
                if (long_options[option_index].flag != 0) {
                    break;
                }
                std::cout << "option " << long_options[option_index].name;
                if (optarg) {
                    std::cout << " with arg " << optarg;
                }
                std::cout << std::endl;
                break;
            case 't':
                threads = std::stoul(optarg);
                break;
            case 'i':
                input_file = optarg;
                break;
            case 'r':
                reference_path = optarg;
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 'h':
                hash_num = std::stoi(optarg);
                break;
            case 'k':
                kmer_size = std::stoi(optarg);
                break;
            case 'g':
                genome_size = std::stoul(optarg);
                break;
            case '?':
                break;
            default:
                std::cout << "Unknown option: " << opt << std::endl;
                break;
        }

    }

    // print help message with required arguments
    if (help_flag) {
        std::cerr << "Usage: " << argv[0] << " [options]" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << "  -h, --help\t\t\tPrint this help message" << std::endl;
        std::cerr << "  -i, --input\t\t\tInput file name" << std::endl;
        std::cerr << "  -o, --output\t\t\tOutput prefix" << std::endl;
        std::cerr << "  -r, --reference\t\tReference path" << std::endl;
        std::cerr << "  -t, --threads\t\t\tNumber of threads to use (default: 1)" << std::endl;
        std::cerr << "  -v, --verbose\t\t\tVerbose output" << std::endl;
        exit(0);
    }

    // print error message if input file is not provided
    if (input_file.empty()) {
        std::cerr << "Input file is required. Use -h or --help for more information." << std::endl;
        exit(1);
    }

    // print error message if reference path is not provided
    if (reference_path.empty()) {
        std::cerr << "Reference path is required. Use -h or --help for more information." << std::endl;
        exit(1);
    }

    // print error message if threads is not provided
    if (threads == 0) {
        std::cerr << "Threads is required. Use -h or --help for more information." << std::endl;
        exit(1);
    }

    // print error message if genome size is not provided
    if (genome_size == 0) {
        std::cerr << "Genome size is required. Use -h or --help for more information." << std::endl;
        exit(1);
    }

    // print log of current parameters if verbose flag is set
    if (verbose_flag) {
        std::cerr << "Input file: " << input_file
                    << "Output prefix: " << output_prefix
                    << "Reference path: " << reference_path
                    << "Threads: " << threads
                    << "Hash number: " << hash_num
                    << "Kmer size: " << kmer_size
                    << "Genome size: " << genome_size << std::endl;
    }




    omp_set_num_threads(threads);
    std::ofstream output_file(output_prefix + ".results.tsv");

    if (verbose_flag) {
        std::cerr << "Reading reference file: " << reference_path << std::endl;
    }
    btllib::MIBloomFilter<uint32_t> mi_bf(calc_optimal_size(genome_size / 3 * 6, hash_num, 0.1), hash_num);

    std::unordered_map<std::string, uint32_t> seq_ID_to_miBf_ID;
    {
        uint32_t miBf_ID = 1;
        btllib::SeqReader reader(reference_path, btllib::SeqReader::Flag::LONG_MODE);
        for (const auto record : reader) {            
            // insert record.id into seq_ID_to_miBf_ID with value miBf_ID
            seq_ID_to_miBf_ID[record.id] = miBf_ID;
        }
    }
    
    for (int stage = 0; stage < 3; stage++) {
        btllib::SeqReader reader(reference_path, btllib::SeqReader::Flag::LONG_MODE);
#pragma omp parallel
        for (const auto record : reader) {
            std::vector<std::string> protein = sixframe_translate(record.seq);
            auto& miBf_ID = seq_ID_to_miBf_ID[record.id];
            for (uint8_t i = 0; i < protein.size(); i++) {
                AAHash itr(protein[i], hash_num, kmer_size);
                while (itr != AAHash::end()) {
                    if (stage == 0) {
                        mi_bf.insert_bv(*itr);
                    } else if (stage == 1) {
                        mi_bf.insert_id(*itr, miBf_ID);
                    } else {
                        mi_bf.insert_saturation(*itr, miBf_ID);
                    }
                    ++itr;
                }
            }
        }

    }


    output_file << "name\thits\texpected_hits" << std::endl;

    btllib::SeqReader reader(input_file, btllib::SeqReader::Flag::LONG_MODE);
    if (verbose_flag) {
        std::cerr << "Reading input file: " << input_file << std::endl;
    }
#pragma omp parallel
  for (const auto record : reader) {
    std::vector<std::string> protein = sixframe_translate(record.seq);
    std::map<uint32_t, size_t> id_to_hits;
    size_t expected_hits = protein[0].size() - kmer_size + 1;
    for (uint8_t i = 0; i < protein.size(); i++) {
                AAHash itr(protein[i], hash_num, kmer_size);
                while (itr != AAHash::end()) {
                    auto temp_ID_hits =  mi_bf.get_id(*itr); // change this to avoid reallocating memory
                    for (auto& ID_hits : temp_ID_hits) {
                        if (id_to_hits.find(ID_hits) == id_to_hits.end()) {
                            id_to_hits[ID_hits] = 1;
                        } else {
                            id_to_hits[ID_hits]++;
                        }
                    }
                }
    
    }
    auto max_hits = std::max_element(id_to_hits.begin(), id_to_hits.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
#pragma omp critical 
    output_file << record.id << "\t" << max_hits->second << "\t" << expected_hits << std::endl;
  }





}

