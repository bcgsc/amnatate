#include <fstream>
#include <iostream>
#include <getopt.h>
#include <omp.h>
#include <string>
#include <vector>

#include <btllib/seq.hpp>
#include <btllib/seq_reader.hpp>
#include <Sequence/Translate.hpp>


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
    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"verbose", no_argument, 0, 'v'},
        {0, 0, 0, 0}
    };
    int verbose_flag = 0;
    int help_flag = 0;
    size_t threads = 1;
    std::string input_file = "";
    std::string output_prefix = "_";

    // loop through command line arguments
    while ((opt = getopt_long(argc, argv, "hi:t:o:v", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'h':
                help_flag = 1;
                break;
            case 'i':
                input_file = optarg;
                break;
            case 'o':
                output_prefix = optarg;
                break;
            case 't':
                threads = std::stoul(optarg);
                break;
            case 'v':
                verbose_flag = 1;
                break;
            default:
                // print error message
                break;
        }
    }

    // print help message
    if (help_flag) {
        std::cerr << "Usage: " << argv[0] << " [options] input_file"
                    << "Options:"
                    << "  -h, --help            display this help and exit"
                    << "  -i, --input           input file"
                    << "  -t, --threads         number of threads"
                    << "  -o, --output          output prefix"
                    << "  -v, --verbose         enable verbose output" << std::endl;
        return 0;
    }
    
    // check if input file is provided
    if (input_file.empty()) {
        std::cerr << "Error: input file is not provided" << std::endl;
        exit(1);
    }

    omp_set_num_threads(threads);
    std::ofstream output_file(output_prefix + ".faa");
    if (verbose_flag) {
        std::cerr << "Reading input file: " << input_file << std::endl;
    }
    btllib::SeqReader reader(input_file, btllib::SeqReader::Flag::LONG_MODE);
#pragma omp parallel
  for (const auto record : reader) {
    std::vector<std::string> protein = sixframe_translate(record.seq);
    for (uint8_t i = 0; i < protein.size(); i++) {
#pragma omp critical 
      output_file << ">" << record.id << "_" << i << std::endl << protein[i] << std::endl;
    }
  }




}

