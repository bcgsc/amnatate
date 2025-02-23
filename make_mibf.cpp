#include <algorithm>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <iomanip>
#include <map>
#include <omp.h>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <btllib/aahash.hpp>
#include <btllib/seq.hpp>
#include <btllib/seq_reader.hpp>
#include <btllib/mi_bloom_filter.hpp>
#include <Sequence/Translate.hpp>

size_t calc_optimal_size(size_t entries, unsigned hash_num, double occupancy)
{
    size_t non64ApproxVal =
        size_t(-double(entries) * double(hash_num) / log(1.0 - occupancy));
    return non64ApproxVal + (64 - non64ApproxVal % 64);
}

btllib::MIBloomFilter<uint64_t> make_small_mibf(const std::string& seq, const size_t hash_num, const size_t kmer_size) 
{
    size_t filter_size = seq.size() * 3;
    if (filter_size < 100000) {
        filter_size = 100000;
    }
    btllib::MIBloomFilter<uint64_t> mi_bf(calc_optimal_size(filter_size, 1, 0.1), 1);
    for (int stage = 0; stage < 3; stage++)
    {

        btllib::AAHash itr(seq, hash_num, kmer_size, 1);
        btllib::AAHash itr2(seq, hash_num, kmer_size, 2);
        btllib::AAHash itr3(seq, hash_num, kmer_size, 3);
        size_t miBf_ID = 1;
        

        while (itr.roll() && itr2.roll() && itr3.roll())
        {
            if (stage == 0)
            {
                mi_bf.insert_bv(itr.hashes());
                mi_bf.insert_bv(itr2.hashes());
                mi_bf.insert_bv(itr3.hashes());
            }
            else if (stage == 1)
            {
                uint64_t new_ID = (uint64_t)miBf_ID << 32 | itr.get_pos();
                mi_bf.insert_id(itr.hashes(), new_ID);
                mi_bf.insert_id(itr2.hashes(), new_ID);
                mi_bf.insert_id(itr3.hashes(), new_ID);
            }
            else
            {
                uint64_t new_ID = (uint64_t)miBf_ID << 32 | itr.get_pos();
                mi_bf.insert_saturation(itr.hashes(), new_ID);
                mi_bf.insert_saturation(itr2.hashes(), new_ID);
                mi_bf.insert_saturation(itr3.hashes(), new_ID);
            }
        }
        
        if (stage == 0)
        {
            mi_bf.complete_bv_insertion();
        }
    }
    return mi_bf;

}

// main function that takes in command line arguments using getopt
int main(int argc, char **argv)
{
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
        {0, 0, 0, 0}};

    // loop through command line arguments
    while ((opt = getopt_long(argc, argv, "g:h:i:t:o:r:k:", long_options, &option_index)) != -1)
    {
        switch (opt)
        {
        case 0:
            if (long_options[option_index].flag != 0)
            {
                break;
            }
            std::cout << "option " << long_options[option_index].name;
            if (optarg)
            {
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
    if (help_flag)
    {
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
    if (input_file.empty())
    {
        std::cerr << "Input file is required. Use -h or --help for more information." << std::endl;
        exit(1);
    }

    // print error message if reference path is not provided
    if (reference_path.empty())
    {
        std::cerr << "Reference path is required. Use -h or --help for more information." << std::endl;
        exit(1);
    }

    // print error message if threads is not provided
    if (threads == 0)
    {
        std::cerr << "Threads is required. Use -h or --help for more information." << std::endl;
        exit(1);
    }

    // print log of current parameters if verbose flag is set
    if (verbose_flag)
    {
        std::cerr << "Input file: " << input_file << "\n"
                  << "Output prefix: " << output_prefix << "\n"
                  << "Reference path: " << reference_path << "\n"
                  << "Threads: " << threads << "\n"
                  << "Hash number: " << (uint64_t)hash_num << "\n"
                  << "Kmer size: " << (uint64_t)kmer_size << "\n"
                  << "Genome size: " << genome_size << std::endl;
    }

    omp_set_num_threads(threads);

    std::unordered_map<std::string, uint32_t> seq_ID_to_miBf_ID;


    if (verbose_flag)
    {
        std::cerr << "Reading reference file: " << reference_path << std::endl;
    }

    // read through reference file which is a fasta file and count the number of characters in the sequences and assign it genome_size
    {
    btllib::SeqReader ref_reader(reference_path, btllib::SeqReader::Flag::LONG_MODE);
    for (const auto record : ref_reader)
    {
        genome_size += record.seq.size();
    }

    btllib::MIBloomFilter<uint64_t> mi_bf(calc_optimal_size(std::max<size_t>(genome_size * 3, 1000000), hash_num, 0.1), hash_num);
    // btllib::MIBloomFilter<uint64_t> mi_bf(calc_optimal_size(1000000000, hash_num, 0.1), hash_num);

    if (verbose_flag)
    {
        std::cerr << "Creating seq_id to ID table" << std::endl;
    }

    
    std::unordered_map<uint32_t, std::pair<std::string, size_t>> miBf_ID_to_seq_ID_and_len;
    std::unordered_map<uint32_t, std::string> miBf_ID_to_seq;
    {
        uint32_t miBf_ID = 1;
        btllib::SeqReader reader(reference_path, btllib::SeqReader::Flag::LONG_MODE);
        for (const auto record : reader)
        {
            // insert record.id into seq_ID_to_miBf_ID with value miBf_ID
            seq_ID_to_miBf_ID[record.id] = miBf_ID;
            // insert miBf_ID into miBf_ID_to_seq_ID_and_len with value record.id and record.seq.size()
            miBf_ID_to_seq_ID_and_len[miBf_ID] = std::make_pair(record.id, record.seq.size());
            miBf_ID_to_seq[miBf_ID] = record.seq;
            /*//std::cerr << "seq name: " << record.id << " miBf_ID: " << miBf_ID << std::endl;
            //std::cerr << "seq size: " << record.seq.size() << std::endl;
            //std::cerr << "seq: " << record.seq << std::endl; */
            ++miBf_ID;
        }
    }
    if (verbose_flag)
    {
        std::cerr << "Making miBF" << std::endl;
    }
    auto sTime = omp_get_wtime();
    for (int stage = 0; stage < 3; stage++)
    {
        if (verbose_flag)
        {
            std::cerr << "stage:" << stage << std::endl;
        }
        btllib::SeqReader reader(reference_path, btllib::SeqReader::Flag::LONG_MODE);
#pragma omp parallel
        for (const auto record : reader)
        {

            btllib::AAHash itr(record.seq, hash_num, kmer_size, 1);
            btllib::AAHash itr2(record.seq, hash_num, kmer_size, 2);
            btllib::AAHash itr3(record.seq, hash_num, kmer_size, 3);
            auto &miBf_ID = seq_ID_to_miBf_ID[record.id];
            

            while (itr.roll() && itr2.roll() && itr3.roll())
            //while (itr.roll())
            {
                if (stage == 0)
                {
                    mi_bf.insert_bv(itr.hashes());
                    mi_bf.insert_bv(itr2.hashes());
                    mi_bf.insert_bv(itr3.hashes());
                }
                else if (stage == 1)
                {
                    uint64_t new_ID = (uint64_t)miBf_ID << 32 | itr.get_pos();
                    mi_bf.insert_id(itr.hashes(), new_ID);
                    mi_bf.insert_id(itr2.hashes(), new_ID);
                    mi_bf.insert_id(itr3.hashes(), new_ID);
                }
                else
                {
                    uint64_t new_ID = (uint64_t)miBf_ID << 32 | itr.get_pos();
                    mi_bf.insert_saturation(itr.hashes(), new_ID);
                    mi_bf.insert_saturation(itr2.hashes(), new_ID);
                    mi_bf.insert_saturation(itr3.hashes(), new_ID);
                }
            }
        }
        if (stage == 0)
        {
            mi_bf.complete_bv_insertion();
        }
    }
    std::cerr << "finished making MiBF" << std::endl;
    std::cerr << "in " << std::setprecision(4) << std::fixed << omp_get_wtime() - sTime
              << "\n";

    mi_bf.save(output_prefix + ".mibf");
    }

    {
    btllib::SeqReader reader(reference_path, btllib::SeqReader::Flag::LONG_MODE);
#pragma omp parallel
        for (const auto record : reader)
        {

            auto mi_bf = make_small_mibf(record.seq, hash_num, kmer_size);
            auto &miBf_ID = seq_ID_to_miBf_ID[record.id];
            

            mi_bf.save(std::to_string(miBf_ID) + ".mibf");
        }

    
    std::cerr << "finished making small MiBF" << std::endl;



    
    }
}