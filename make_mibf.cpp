#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <omp.h>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <argparse/argparse.hpp>

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

// main function
int main(int argc, char* argv[]) {
    // Create argparse parser
    argparse::ArgumentParser program("MiBF_Construction");

    // Add arguments

    program.add_argument("--help")
           .help("Display this help message")
           .default_value(false)
           .implicit_value(true);

    program.add_argument("-i", "--input")
        .help("Input file name")
        .required();
    
    program.add_argument("-o", "--output")
        .help("Output prefix")
        .default_value(std::string("_"));
    
    program.add_argument("-r", "--reference")
        .help("Reference path")
        .required();

    program.add_argument("-t", "--threads")
        .help("Number of threads to use")
        .default_value(size_t(1))
        .scan<'u', size_t>();
    program.add_argument("-h", "--hash")
        .help("Number of hash functions")
        .default_value(uint8_t(9))
        .scan<'u', uint8_t>();
    
    program.add_argument("-k", "--kmer")
        .help("K-mer size")
        .default_value(uint8_t(9))
        .scan<'u', uint8_t>();
    
    program.add_argument("-v", "--verbose")
           .help("Verbose output")
           .default_value(false)
           .implicit_value(true);


    try {
        // Parse arguments
        program.parse_args(argc, argv);
    } catch (const std::exception &e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // Extract parsed arguments
    bool help_flag = program.get<bool>("--help");
    bool verbose_flag = program.get<bool>("--verbose");
    size_t threads = program.get<size_t>("--threads");
    std::string input_file = program.get<std::string>("--input");
    std::string reference_path = program.get<std::string>("--reference");
    std::string output_prefix = program.get<std::string>("--output");
    uint8_t hash_num = program.get<uint8_t>("--hash");
    uint8_t kmer_size = program.get<uint8_t>("--kmer");

    if (help_flag)
    {
        std::cerr << program << std::endl;
        return 0;
    }


    // print log of current parameters if verbose flag is set
    if (verbose_flag)
    {
        std::cerr << "Input file: " << input_file << "\n"
                  << "Output prefix: " << output_prefix << "\n"
                  << "Reference path: " << reference_path << "\n"
                  << "Threads: " << threads << "\n"
                  << "Hash number: " << (uint64_t)hash_num << "\n"
                  << "Kmer size: " << (uint64_t)kmer_size;
    }

    omp_set_num_threads(threads);

    std::unordered_map<std::string, uint32_t> seq_ID_to_miBf_ID;


    if (verbose_flag)
    {
        std::cerr << "Reading reference file: " << reference_path << std::endl;
    }

    // read through reference file which is a fasta file and count the number of characters in the sequences and assign it genome_size
    {
    size_t genome_size = 0;
    btllib::SeqReader ref_reader(reference_path, btllib::SeqReader::Flag::LONG_MODE);
    for (const auto record : ref_reader)
    {
        genome_size += record.seq.size();
    }

    btllib::MIBloomFilter<uint64_t> mi_bf(calc_optimal_size(std::max<size_t>(genome_size * 3, 1000000), hash_num, 0.1), hash_num);

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
            seq_ID_to_miBf_ID[record.id] = miBf_ID;
            miBf_ID_to_seq_ID_and_len[miBf_ID] = std::make_pair(record.id, record.seq.size());
            miBf_ID_to_seq[miBf_ID] = record.seq;
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
