#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Sequence/Translate.hpp>
#include <argparse/argparse.hpp>
#include <btllib/aahash.hpp>
#include <btllib/mi_bloom_filter.hpp>
#include <btllib/seq.hpp>
#include <btllib/seq_reader.hpp>

static constexpr size_t MIN_TARGETTED_MIBF_SIZE = 1000;
static constexpr size_t MIN_MAIN_MIBF_SIZE = 100000;
static constexpr double TARGET_FALSE_POSITIVE_RATE = 0.1;
static constexpr size_t STAGES = 3;
static constexpr size_t SIZE_MULTIPLIER = 3;
static constexpr uint32_t HASH_ID_SHIFT = 32;

size_t calc_optimal_size(size_t entries, unsigned hash_num, double occupancy) {
  size_t approx = static_cast<size_t>(-static_cast<double>(entries) *
                                      static_cast<double>(hash_num) /
                                      std::log(1.0 - occupancy));
  return approx + (64 - approx % 64);
}

btllib::MIBloomFilter<uint64_t>
make_small_mibf(const std::string &seq, size_t hash_num, size_t kmer_size) {
  size_t filter_size =
      std::max(seq.size() * SIZE_MULTIPLIER, MIN_TARGETTED_MIBF_SIZE);
  btllib::MIBloomFilter<uint64_t> mi_bf(
      calc_optimal_size(filter_size, hash_num, TARGET_FALSE_POSITIVE_RATE),
      hash_num);

  for (size_t stage = 0; stage < STAGES; ++stage) {
    btllib::AAHash itr(seq, hash_num, kmer_size, 1);
    btllib::AAHash itr2(seq, hash_num, kmer_size, 2);
    btllib::AAHash itr3(seq, hash_num, kmer_size, 3);
    size_t miBf_ID = 1;

    while (itr.roll() && itr2.roll() && itr3.roll()) {
      uint64_t new_ID =
          (static_cast<uint64_t>(miBf_ID) << HASH_ID_SHIFT) | itr.get_pos();

      if (stage == 0) {
        mi_bf.insert_bv(itr.hashes());
        mi_bf.insert_bv(itr2.hashes());
        mi_bf.insert_bv(itr3.hashes());
      } else if (stage == 1) {
        mi_bf.insert_id(itr.hashes(), new_ID);
        mi_bf.insert_id(itr2.hashes(), new_ID);
        mi_bf.insert_id(itr3.hashes(), new_ID);
      } else {
        mi_bf.insert_saturation(itr.hashes(), new_ID);
        mi_bf.insert_saturation(itr2.hashes(), new_ID);
        mi_bf.insert_saturation(itr3.hashes(), new_ID);
      }
    }

    if (stage == 0) {
      mi_bf.complete_bv_insertion();
    }
  }

  return mi_bf;
}

int main(int argc, char *argv[]) {
  argparse::ArgumentParser program("make_mibf");

  program.add_argument("--help")
      .help("Display this help message")
      .default_value(false)
      .implicit_value(true);

  program.add_argument("-o", "--output")
      .help("Output prefix")
      .default_value(std::string("_"));

  program.add_argument("-r", "--reference").help("Reference path").required();

  program.add_argument("-t", "--threads")
      .help("Number of threads to use")
      .default_value(size_t(1))
      .scan<'u', size_t>();

  program.add_argument("-h", "--hash")
      .help("Number of hash functions")
      .default_value(size_t(9))
      .scan<'u', size_t>();

  program.add_argument("-k", "--kmer")
      .help("K-mer size")
      .default_value(size_t(9))
      .scan<'u', size_t>();

  program.add_argument("-v", "--verbose")
      .help("Verbose output")
      .default_value(false)
      .implicit_value(true);

  program.add_argument("-rks", "--rescue_kmer")
      .help("Rescue k-mer size")
      .default_value(size_t(4))
      .scan<'u', size_t>();

  program.add_argument("--version")
      .help("Display version information")
      .default_value(false)
      .implicit_value(true);

  bool help_flag = std::any_of(argv, argv + argc, [](const char *arg) {
    return std::string(arg) == "--help";
  });

  if (help_flag) {
    std::cerr << program << std::endl;
    return 0;
  }

  try {
    program.parse_args(argc, argv);
  } catch (const std::exception &err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  if (program.get<bool>("--version")) {
    std::cout << "make_mibf version 1.0.0" << std::endl;
    return 0;
  }

  bool verbose_flag = program.get<bool>("--verbose");
  size_t threads = program.get<size_t>("--threads");
  std::string reference_path = program.get<std::string>("--reference");
  std::string output_prefix = program.get<std::string>("--output");
  uint8_t hash_num = static_cast<uint8_t>(program.get<size_t>("--hash"));
  uint8_t kmer_size = static_cast<uint8_t>(program.get<size_t>("--kmer"));
  size_t rescue_kmer_size = program.get<size_t>("--rescue_kmer");

  if (reference_path.empty()) {
    std::cerr
        << "Reference path is required. Use -h or --help for more information."
        << std::endl;
    return 1;
  }

  if (verbose_flag) {
    std::cerr << "Reference path: " << reference_path << "\n"
              << "Output prefix: " << output_prefix << "\n"
              << "Threads: " << threads << "\n"
              << "Hash number: " << static_cast<uint64_t>(hash_num) << "\n"
              << "Kmer size: " << static_cast<uint64_t>(kmer_size) << "\n"
              << "Rescue kmer size: " << static_cast<uint64_t>(rescue_kmer_size)
              << "\n"
              << std::endl;
  }

  omp_set_num_threads(threads);

  std::unordered_map<std::string, uint32_t> seq_ID_to_miBf_ID;

  if (verbose_flag) {
    std::cerr << "Reading reference file: " << reference_path << std::endl;
  }

  size_t genome_size = 0;
  btllib::SeqReader ref_reader(reference_path,
                               btllib::SeqReader::Flag::LONG_MODE);
  for (const auto &record : ref_reader) {
    genome_size += record.seq.size();
  }

  btllib::MIBloomFilter<uint64_t> mi_bf(
      calc_optimal_size(
          std::max<size_t>(genome_size * SIZE_MULTIPLIER, MIN_MAIN_MIBF_SIZE),
          hash_num, TARGET_FALSE_POSITIVE_RATE),
      hash_num);

  if (verbose_flag) {
    std::cerr << "Creating seq_id to ID table" << std::endl;
  }

  std::unordered_map<uint32_t, std::pair<std::string, size_t>>
      miBf_ID_to_seq_ID_and_len;
  std::unordered_map<uint32_t, std::string> miBf_ID_to_seq;
  {
    uint32_t miBf_ID = 1;
    btllib::SeqReader reader(reference_path,
                             btllib::SeqReader::Flag::LONG_MODE);
    for (const auto &record : reader) {
      seq_ID_to_miBf_ID[record.id] = miBf_ID;
      miBf_ID_to_seq_ID_and_len[miBf_ID] = {record.id, record.seq.size()};
      miBf_ID_to_seq[miBf_ID] = record.seq;
      ++miBf_ID;
    }
  }

  if (verbose_flag) {
    std::cerr << "Making miBf" << std::endl;
  }

  auto sTime = omp_get_wtime();
  for (size_t stage = 0; stage < STAGES; ++stage) {
    if (verbose_flag) {
      std::cerr << "stage: " << stage << std::endl;
    }

    btllib::SeqReader reader(reference_path,
                             btllib::SeqReader::Flag::LONG_MODE);
#pragma omp parallel
    for (const auto &record : reader) {
      btllib::AAHash itr(record.seq, hash_num, kmer_size, 1);
      btllib::AAHash itr2(record.seq, hash_num, kmer_size, 2);
      btllib::AAHash itr3(record.seq, hash_num, kmer_size, 3);
      const auto &miBf_ID = seq_ID_to_miBf_ID[record.id];

      while (itr.roll() && itr2.roll() && itr3.roll()) {
        uint64_t new_ID =
            (static_cast<uint64_t>(miBf_ID) << HASH_ID_SHIFT) | itr.get_pos();

        if (stage == 0) {
          mi_bf.insert_bv(itr.hashes());
          mi_bf.insert_bv(itr2.hashes());
          mi_bf.insert_bv(itr3.hashes());
        } else if (stage == 1) {
          mi_bf.insert_id(itr.hashes(), new_ID);
          mi_bf.insert_id(itr2.hashes(), new_ID);
          mi_bf.insert_id(itr3.hashes(), new_ID);
        } else {
          mi_bf.insert_saturation(itr.hashes(), new_ID);
          mi_bf.insert_saturation(itr2.hashes(), new_ID);
          mi_bf.insert_saturation(itr3.hashes(), new_ID);
        }
      }
    }

    if (stage == 0) {
      mi_bf.complete_bv_insertion();
    }
  }

  std::cerr << "finished making miBf" << std::endl;
  std::cerr << "in " << std::fixed << std::setprecision(4)
            << omp_get_wtime() - sTime << std::endl;
  mi_bf.save(output_prefix + ".mibf");

  btllib::SeqReader reader(reference_path, btllib::SeqReader::Flag::LONG_MODE);
#pragma omp parallel
  for (const auto &record : reader) {
    auto mi_bf_small = make_small_mibf(record.seq, hash_num, rescue_kmer_size);
    const auto &miBf_ID = seq_ID_to_miBf_ID[record.id];
    mi_bf_small.save(output_prefix + "." + std::to_string(miBf_ID) + ".mibf");
  }

  std::cerr << "finished making small miBf" << std::endl;
  return 0;
}
