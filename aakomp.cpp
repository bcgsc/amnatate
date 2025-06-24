#include <algorithm>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Third-party libraries
#include <Sequence/Translate.hpp>
#include <argparse/argparse.hpp>
#include <boost/math/distributions/empirical_cumulative_distribution_function.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <btllib/aahash.hpp>
#include <btllib/mi_bloom_filter.hpp>
#include <btllib/seq.hpp>
#include <btllib/seq_reader.hpp>

using boost::math::empirical_cumulative_distribution_function;
using boost::math::quadrature::trapezoidal;

static constexpr uint32_t HASH_ID_SHIFT = 32;
static constexpr uint32_t HASH_POS_MASK = 0xFFFFFFFF;
static constexpr size_t FRAMES = 3;
static constexpr size_t ORIENTATIONS = 2;
static constexpr size_t MINIMUM_CONSECUTIVE_HIT = 5;

struct FrameBlock {
  size_t frame;
  size_t block_id;
  size_t query_start_in_prot_space;

  FrameBlock(size_t frame, size_t block_id, size_t query_start)
      : frame(frame), block_id(block_id),
        query_start_in_prot_space(query_start) {}
};

struct GFFEntry {
  std::string query_name;
  size_t hit_pos_start;
  size_t hit_pos_end;
  double score;
  std::string strand;
  std::string hit_name;

  GFFEntry(const std::string &query, size_t start, size_t end, double score,
           const std::string &strand, const std::string &hit)
      : query_name(query), hit_pos_start(start), hit_pos_end(end), score(score),
        strand(strand), hit_name(hit) {}
};

struct GFFEntryComparator {
  bool operator()(const GFFEntry &lhs, const GFFEntry &rhs) const noexcept {
    if (lhs.query_name == rhs.query_name)
      return lhs.hit_pos_start < rhs.hit_pos_start;
    return lhs.query_name < rhs.query_name;
  }
};

struct FrameBlockComparator {
  bool operator()(const FrameBlock &a, const FrameBlock &b) const noexcept {
    return a.query_start_in_prot_space < b.query_start_in_prot_space;
  }
};

size_t look_ahead(
    const std::vector<std::reference_wrapper<const FrameBlock>> &blocks,
    size_t current_index, size_t query_end,
    const std::unordered_map<
        size_t, std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>>>
        &frame_index,
    size_t offset) {
  if (current_index + 1 >= blocks.size()) {
    return 0;
  }

  const FrameBlock &next = blocks[current_index + 1].get();
  auto next_block_id = next.block_id;
  auto next_frame = next.frame;

  if (query_end < frame_index.at(next_frame).at(next_block_id).first) {
    return 0;
  }

  for (size_t i = 2; i <= offset; ++i) {
    if (current_index + i >= blocks.size()) {
      return 0;
    }

    const FrameBlock &block = blocks[current_index + i].get();
    auto block_id = block.block_id;
    auto frame = block.frame;

    if (query_end < frame_index.at(frame).at(block_id).first) {
      return i - 1;
    }
  }

  return 0;
}

void process_hashes(
    const std::vector<uint64_t> &temp_ID_pos,
    std::unordered_set<uint32_t> &id_set,
    std::unordered_map<uint32_t, std::set<uint32_t>> &id_to_pos_set,
    bool &extend_block, std::vector<uint32_t> &ids_vec,
    std::vector<uint32_t> &temp_pos_vec,
    const btllib::MIBloomFilter<uint64_t> &mi_bf) {
  bool found = false;

  for (const auto &id : id_set) {
    for (const auto &id_pos : temp_ID_pos) {
      auto demasked = id_pos & mi_bf.ANTI_MASK;
      if (id == (demasked >> HASH_ID_SHIFT)) {
        found = true;
        break;
      }
    }
    if (found)
      break;
  }

  if (found) {
    for (const auto &id_pos : temp_ID_pos) {
      auto demasked = id_pos & mi_bf.ANTI_MASK;
      ids_vec.push_back(demasked >> HASH_ID_SHIFT);
      temp_pos_vec.push_back(demasked & HASH_POS_MASK);
    }

    size_t expected_size = ids_vec.size();
    std::vector<uint32_t> new_ids;
    std::vector<uint32_t> new_pos;
    new_ids.reserve(expected_size);
    new_pos.reserve(expected_size);

    for (size_t i = 0; i < expected_size; ++i) {
      uint32_t id = ids_vec[i];
      uint32_t pos = temp_pos_vec[i];

      if (id_set.count(id)) {
        const auto &pos_set = id_to_pos_set[id];
        if (!pos_set.empty() && pos == *pos_set.rbegin() + 1) {
          id_to_pos_set[id].insert(pos);
        }
      } else {
        new_ids.push_back(id);
        new_pos.push_back(pos);
      }
    }

    std::unordered_set<uint32_t> new_id_set(new_ids.begin(), new_ids.end());
    for (const auto &id : new_id_set) {
      std::set<uint32_t> pos_set;
      for (size_t i = 0; i < new_ids.size(); ++i) {
        if (new_ids[i] == id) {
          pos_set.insert(new_pos[i]);
        }
      }
      if (!pos_set.empty()) {
        id_to_pos_set[id].insert(*pos_set.begin());
      }
    }

  } else {
    bool saturated = true;
    for (const auto &id_pos : temp_ID_pos) {
      if (id_pos < mi_bf.MASK) {
        saturated = false;
        break;
      }
    }

    if (!saturated) {
      extend_block = false;
    } else {
      extend_block = true;
      for (auto &[id, pos_set] : id_to_pos_set) {
        if (!pos_set.empty()) {
          pos_set.insert(*pos_set.rbegin() + 1);
        }
      }
    }
  }
}

std::vector<std::string> sixframe_translate(const std::string &dna) {
  std::vector<std::string> protein;
  protein.reserve(6);

  const std::string rev_dna = btllib::get_reverse_complement(dna);

  // Forward frames
  for (size_t frame = 0; frame < FRAMES; ++frame) {
    protein.push_back(Sequence::Translate(dna.begin() + frame, dna.end()));
  }

  // Reverse frames
  for (size_t frame = 0; frame < FRAMES; ++frame) {
    protein.push_back(
        Sequence::Translate(rev_dna.begin() + frame, rev_dna.end()));
  }

  return protein;
}

void fill_in_gaps(
    std::vector<std::tuple<size_t, size_t>> &start_end_pos_vec,
    std::vector<std::tuple<size_t, size_t>> &start_end_pos_in_tar_space_vec,
    size_t &adjusted_kmer_counts, size_t hash_num, size_t rescue_kmer_size,
    const std::vector<std::string> &sixframed_xlated_proteins, size_t ori,
    size_t kmer_size, uint32_t miBf_ID, const std::string &db_path,
    const std::string &mibf_prefix) {
  const size_t MIN_KMER_CHAIN_GAP = rescue_kmer_size - 1;
  const size_t TARGET_GAP_ADJUSTMENT = rescue_kmer_size + 1;

  std::sort(start_end_pos_vec.begin(), start_end_pos_vec.end(),
            [](const auto &a, const auto &b) {
              return std::get<0>(a) < std::get<0>(b);
            });

  std::vector<std::tuple<size_t, size_t>> gap_vec;
  gap_vec.reserve(start_end_pos_vec.size());

  for (size_t i = 0; i + 1 < start_end_pos_vec.size(); ++i) {
    size_t curr_end = std::get<1>(start_end_pos_vec[i]);
    size_t next_start = std::get<0>(start_end_pos_vec[i + 1]);
    if (next_start > curr_end + kmer_size) {
      gap_vec.emplace_back(curr_end, next_start);
    }
  }

  if (gap_vec.empty())
    return;

  std::vector<std::tuple<size_t, size_t>> gap_in_tar_space_vec;
  gap_in_tar_space_vec.reserve(start_end_pos_in_tar_space_vec.size());

  for (size_t i = 0; i + 1 < start_end_pos_in_tar_space_vec.size(); ++i) {
    size_t curr_end =
        std::get<1>(start_end_pos_in_tar_space_vec[i]) + TARGET_GAP_ADJUSTMENT;
    size_t next_start = std::get<0>(start_end_pos_in_tar_space_vec[i + 1]);
    if (next_start > curr_end + kmer_size) {
      gap_in_tar_space_vec.emplace_back(curr_end, next_start);
    } else {
      gap_in_tar_space_vec.emplace_back(0, 0); // dummy
    }
  }

  if (gap_in_tar_space_vec.empty())
    return;

  std::vector<std::unordered_set<size_t>> gap_index_sets;
  gap_index_sets.reserve(gap_in_tar_space_vec.size());

  for (const auto &gap : gap_in_tar_space_vec) {
    std::unordered_set<size_t> index_set;
    for (size_t i = std::get<0>(gap); i < std::get<1>(gap); ++i) {
      index_set.insert(i);
    }
    gap_index_sets.emplace_back(std::move(index_set));
  }

  if (gap_index_sets.empty())
    return;

  std::string small_mibf_path =
      db_path + "/" + mibf_prefix + "." + std::to_string(miBf_ID) + ".mibf";
  btllib::MIBloomFilter<uint64_t> small_mi_bf(small_mibf_path);

  std::vector<std::vector<std::tuple<size_t, size_t>>>
      kmer_pos_per_frame_and_level(FRAMES);

  for (size_t frame = 0; frame < FRAMES; ++frame) {
    const std::string &protein =
        sixframed_xlated_proteins[frame + ori * FRAMES];
    for (size_t level = 1; level <= 3; ++level) {
      for (size_t g = 0; g < gap_vec.size(); ++g) {
        if (gap_in_tar_space_vec[g] == std::make_tuple(0ul, 0ul))
          continue;

        size_t gap_start = std::get<0>(gap_vec[g]);
        size_t gap_end = std::get<1>(gap_vec[g]);
        const auto &index_set = gap_index_sets[g];

        btllib::AAHash aahash(protein, hash_num, rescue_kmer_size, level,
                              gap_start - 1);
        aahash.roll();

        while (aahash.get_pos() <= gap_end + 1) {
          if (small_mi_bf.bv_contains(aahash.hashes())) {
            auto temp_ID_pos = small_mi_bf.get_id(aahash.hashes());
            for (auto id_pos : temp_ID_pos) {
              size_t pos = id_pos & HASH_POS_MASK;
              if (index_set.count(pos)) {
                kmer_pos_per_frame_and_level[frame].emplace_back(
                    pos, aahash.get_pos());
              }
            }
          }
          aahash.roll();
        }
      }
    }
  }

  std::vector<std::tuple<size_t, size_t>> all_kmers;
  std::set<std::tuple<size_t, size_t>> all_kmers_set;
  all_kmers.reserve(100); // 100 is an estimate

  for (size_t frame = 0; frame < FRAMES; ++frame) {
    for (const auto &kmer : kmer_pos_per_frame_and_level[frame]) {
      if (all_kmers_set.insert(kmer).second) {
        all_kmers.emplace_back(kmer);
      }
    }
  }

  std::sort(all_kmers.begin(), all_kmers.end(),
            [](const auto &a, const auto &b) {
              return std::get<1>(a) < std::get<1>(b);
            });

  if (!all_kmers.empty()) {
    size_t prev_pos = std::get<1>(all_kmers[0]);
    for (size_t i = 1; i < all_kmers.size(); ++i) {
      size_t current_pos = std::get<1>(all_kmers[i]);
      size_t gap = current_pos - prev_pos;
      adjusted_kmer_counts +=
          (gap > MIN_KMER_CHAIN_GAP) ? MIN_KMER_CHAIN_GAP : gap;
      prev_pos = current_pos;
    }
    adjusted_kmer_counts += MIN_KMER_CHAIN_GAP;
  }
}

bool explore_frame(btllib::MIBloomFilter<uint64_t> &mi_bf,
                   btllib::AAHash &aahash,
                   std::deque<std::vector<uint32_t>> &miBf_IDs_snapshot,
                   std::deque<std::vector<uint32_t>> &miBf_pos_snapshot,
                   std::unordered_map<uint32_t, size_t> &id_to_count) {
  std::unordered_set<uint32_t> id_set;
  if (miBf_IDs_snapshot.size() >= MINIMUM_CONSECUTIVE_HIT) {
    for (size_t i = 0; i < miBf_IDs_snapshot.front().size(); ++i) {
      id_set.insert(miBf_IDs_snapshot.front()[i]);
    }
    for (auto &ID : id_set) {
      id_to_count[ID]--;
      if (id_to_count[ID] == 0) {
        id_to_count.erase(ID);
      }
    }
    id_set.clear();
    miBf_IDs_snapshot.pop_front();
    miBf_pos_snapshot.pop_front();
  }
  miBf_IDs_snapshot.emplace_back(std::vector<uint32_t>());
  miBf_pos_snapshot.emplace_back(std::vector<uint32_t>());

  if (!mi_bf.bv_contains(aahash.hashes())) {
    miBf_IDs_snapshot.clear();
    miBf_pos_snapshot.clear();
    id_to_count.clear();
    return false;
  }

  auto temp_ID_pos = mi_bf.get_id(aahash.hashes());
  for (auto &ID_pos : temp_ID_pos) {
    auto demasked_ID_pos = ID_pos & mi_bf.ANTI_MASK;
    miBf_IDs_snapshot.back().push_back(demasked_ID_pos >> HASH_ID_SHIFT);
    miBf_pos_snapshot.back().push_back(demasked_ID_pos & HASH_POS_MASK);
  }

  for (size_t j = 0; j < miBf_IDs_snapshot.back().size(); ++j) {
    id_set.insert(miBf_IDs_snapshot.back()[j]);
  }
  for (auto &ID : id_set) {
    if (id_to_count.find(ID) == id_to_count.end()) {
      id_to_count[ID] = 1;
    } else {
      id_to_count[ID]++;
    }
  }
  id_set.clear();

  if (miBf_IDs_snapshot.size() < MINIMUM_CONSECUTIVE_HIT) {
    return false;
  }

  uint32_t temp_mibf_ID = 0;
  size_t temp_max_count = 0;
  for (auto &ID_count : id_to_count) {
    if (ID_count.second > temp_max_count) {
      temp_mibf_ID = ID_count.first;
      temp_max_count = ID_count.second;
    }
  }
  if (temp_mibf_ID == 0 || temp_max_count < MINIMUM_CONSECUTIVE_HIT) {
    return false;
  } else {
    std::vector<uint32_t> ids_to_check;
    for (auto &ID_count : id_to_count) {
      if (ID_count.second == temp_max_count) {
        ids_to_check.push_back(ID_count.first);
      }
    }

    for (auto &ID : ids_to_check) {
      std::set<uint32_t> temp_pos_set;
      for (size_t i = 0; i < miBf_IDs_snapshot.size(); ++i) {
        for (size_t j = 0; j < miBf_IDs_snapshot[i].size(); ++j) {
          if (miBf_IDs_snapshot[i][j] == ID) {
            temp_pos_set.insert(miBf_pos_snapshot[i][j]);
          }
        }
      }

      uint32_t prev_pos = 0;
      size_t counter = 0;
      bool init = false;
      for (auto &pos : temp_pos_set) {
        if (!init) {
          prev_pos = pos;
          init = true;
          continue;
        }
        if (pos - prev_pos == 1) {
          ++counter;
        } else {
          counter = 0;
        }
        prev_pos = pos;

        if (counter >= (MINIMUM_CONSECUTIVE_HIT - 1)) {
          return true;
        }
      }
    }
  }

  return false;
}

int main(int argc, char *argv[]) {
  argparse::ArgumentParser program("aaKomp");

  program.add_argument("--help")
      .help("Print this help message")
      .default_value(false)
      .implicit_value(true);

  program.add_argument("-i", "--input").help("Input file name").required();

  program.add_argument("-o", "--output")
      .help("Output prefix")
      .default_value(std::string("_"));

  program.add_argument("-r", "--reference").help("Reference path").required();

  program.add_argument("-t", "--threads")
      .help("Number of threads to use")
      .default_value(size_t(1))
      .scan<'u', size_t>();

  program.add_argument("-v", "--verbose")
      .help("Verbose output")
      .default_value(false)
      .implicit_value(true);

  program.add_argument("--debug")
      .help("Debug output")
      .default_value(false)
      .implicit_value(true);

  program.add_argument("-m", "--mibf_path").help("miBf file path").required();

  program.add_argument("-h", "--hash")
      .help("Number of hash functions")
      .default_value(size_t(9))
      .scan<'u', size_t>();

  program.add_argument("-k", "--kmer")
      .help("K-mer size")
      .default_value(size_t(9))
      .scan<'u', size_t>();

  program.add_argument("-l", "--lower_bound")
      .help("Lower bound value")
      .default_value(0.7)
      .scan<'g', double>();

  program.add_argument("-rks", "--rescue_kmer")
      .help("Rescue k-mer size")
      .default_value(size_t(4))
      .scan<'u', size_t>();

  program.add_argument("-mo", "--max_offset")
      .help("Maximum offset")
      .default_value(size_t(2))
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
    std::cout << "aaKomp version 1.0.0" << std::endl;
    return 0;
  }

  bool verbose_flag = program.get<bool>("--verbose");
  bool debug_flag = program.get<bool>("--debug");
  size_t threads = program.get<size_t>("--threads");
  std::string input_file = program.get<std::string>("--input");
  std::string reference_path = program.get<std::string>("--reference");
  std::string output_prefix = program.get<std::string>("--output");
  std::string mibf_path = program.get<std::string>("--mibf_path");
  uint8_t hash_num = static_cast<uint8_t>(program.get<size_t>("--hash"));
  uint8_t kmer_size = static_cast<uint8_t>(program.get<size_t>("--kmer"));
  size_t rescue_kmer_size = program.get<size_t>("--rescue_kmer");
  double lower_bound = program.get<double>("--lower_bound");
  size_t max_offset = program.get<size_t>("--max_offset");
  std::string mibf_prefix = std::filesystem::path(mibf_path).stem().string();

  std::string db_path_loc = "./";
  if (!mibf_path.empty()) {
    std::filesystem::path path_obj(mibf_path);
    if (path_obj.has_parent_path()) {
      db_path_loc = path_obj.parent_path().string();
    }
  }

  if (help_flag) {
    std::cerr << program << std::endl;
    return 0;
  }

  if (input_file.empty()) {
    std::cerr
        << "Input file is required. Use -h or --help for more information."
        << std::endl;
    exit(1);
  }

  if (reference_path.empty()) {
    std::cerr
        << "Reference path is required. Use -h or --help for more information."
        << std::endl;
    exit(1);
  }

  if (threads == 0) {
    std::cerr << "Threads is required. Use -h or --help for more information."
              << std::endl;
    exit(1);
  }

  if (input_file.empty()) {
    std::cerr
        << "Input file is required. Use -h or --help for more information."
        << std::endl;
    return 1;
  }

  if (reference_path.empty()) {
    std::cerr
        << "Reference path is required. Use -h or --help for more information."
        << std::endl;
    return 1;
  }

  if (threads == 0) {
    std::cerr << "Threads must be greater than 0. Use -h or --help for more "
                 "information."
              << std::endl;
    return 1;
  }

  if (verbose_flag) {
    std::cerr << "Input file: " << input_file << "\n"
              << "Output prefix: " << output_prefix << "\n"
              << "Reference path: " << reference_path << "\n"
              << "Threads: " << threads << "\n"
              << "Hash number: " << (uint64_t)hash_num << "\n"
              << "Kmer size: " << (uint64_t)kmer_size << "\n"
              << "Rescue kmer size: " << rescue_kmer_size << "\n"
              << "Lower bound: " << lower_bound << "\n"
              << "DB Path Location: " << db_path_loc << "\n"
              << "Max offset: " << max_offset << "\n"
              << std::endl;
  }

  if (debug_flag) {
    std::cerr << "Debugging enabled" << std::endl;
  }

  omp_set_num_threads(threads);

  if (verbose_flag) {
    std::cerr << "Reading reference file: " << reference_path << std::endl;
  }

  if (verbose_flag) {
    std::cerr << "Creating seq_id to ID table" << std::endl;
  }

  std::unordered_map<std::string, uint32_t> seq_ID_to_miBf_ID;
  std::unordered_map<uint32_t, std::pair<std::string, size_t>>
      miBf_ID_to_seq_ID_and_len;
  std::unordered_map<uint32_t, std::string> miBf_ID_to_seq;
  {
    uint32_t miBf_ID = 1;
    btllib::SeqReader reader(reference_path,
                             btllib::SeqReader::Flag::LONG_MODE);
    for (const auto record : reader) {
      if (record.seq.size() < (size_t)kmer_size + 5) {
        continue;
      }
      seq_ID_to_miBf_ID[record.id] = miBf_ID;
      miBf_ID_to_seq_ID_and_len[miBf_ID] =
          std::make_pair(record.id, record.seq.size());
      miBf_ID_to_seq[miBf_ID] = record.seq;
      ++miBf_ID;
    }
  }
  if (verbose_flag) {
    std::cerr << "Reading miBF" << std::endl;
  }

  auto sTime = omp_get_wtime();
  btllib::MIBloomFilter<uint64_t> mi_bf(mibf_path);
  auto sTime2 = omp_get_wtime();
  if (verbose_flag) {
    std::cerr << "finished reading miBf" << std::endl;
    std::cerr << "in " << std::setprecision(4) << std::fixed << sTime2 - sTime
              << "\n";
  }

  std::ofstream output_file;
  std::ofstream gff_file;
  std::ofstream pre_gff_file;

  std::string filename = output_prefix + "_results.tsv";
  output_file.open(filename);
  output_file << "name\tcomplete copies\tincomplete copies\texpected k-mer "
                 "counts\thighest adjusted incomplete k-mer hits"
              << std::endl;
  gff_file.open(output_prefix + ".gff");
  gff_file << "##gff-version 3" << std::endl;
  if (debug_flag) {
    pre_gff_file.open(output_prefix + ".pre.gff");
    pre_gff_file << "##gff-version 3" << std::endl;
  }

  std::set<GFFEntry, GFFEntryComparator> gff_set;
  std::set<GFFEntry, GFFEntryComparator> pre_gff_set;

  btllib::SeqReader reader(input_file, btllib::SeqReader::Flag::LONG_MODE);
  if (verbose_flag) {
    std::cerr << "Reading input file: " << input_file << std::endl;
  }
  struct completeness_struct {
    size_t complete_copies = 0;
    size_t incomplete_copies = 0;
    size_t expected_kmer_counts = 0;
    size_t highest_adjusted_kmer_counts = 0;
  };
  std::unordered_map<std::string, completeness_struct> seq_name_to_completeness;

  for (const auto &seq_ID : seq_ID_to_miBf_ID) {
    seq_name_to_completeness[seq_ID.first] = completeness_struct();
  }

  omp_set_nested(1);

#pragma omp parallel num_threads(threads / 2)
  for (const auto record : reader) {
    std::vector<std::string> sixframed_xlated_proteins =
        sixframe_translate(record.seq);
#pragma omp parallel for num_threads(2)
    for (size_t ori = 0; ori < ORIENTATIONS; ++ori) {
      std::unordered_map<
          size_t, std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>>>
          frame_to_block_id_to_id_and_pos;
      std::map<uint32_t, size_t, std::greater<size_t>>
          id_to_count_across_all_frames;
      std::unordered_map<uint32_t, std::set<FrameBlock, FrameBlockComparator>>
          id_to_FrameBlock_id_and_seq_pos;

      for (size_t frame = 0; frame < FRAMES; ++frame) {
        size_t ori_frame = frame + ori * FRAMES;
        btllib::AAHash aahash(sixframed_xlated_proteins[ori_frame], hash_num,
                              kmer_size, 1);
        btllib::AAHash aahash2(sixframed_xlated_proteins[ori_frame], hash_num,
                               kmer_size, 2);
        btllib::AAHash aahash3(sixframed_xlated_proteins[ori_frame], hash_num,
                               kmer_size, 3);
        aahash.roll();
        aahash2.roll();
        aahash3.roll();
        std::deque<std::vector<uint32_t>> miBf_IDs_snapshot;
        std::deque<std::vector<uint32_t>> miBf_pos_snapshot;
        std::unordered_map<uint32_t, size_t> id_to_count;

        size_t block_id = 0;
        std::unordered_set<uint32_t> id_set;
        while (aahash.get_pos() != std::numeric_limits<size_t>::max()) {
          while (!explore_frame(mi_bf, aahash, miBf_IDs_snapshot,
                                miBf_pos_snapshot, id_to_count) &&
                 aahash.get_pos() != std::numeric_limits<size_t>::max()) {
            aahash.roll();
            aahash2.roll();
            aahash3.roll();
          }
          if (aahash.get_pos() == std::numeric_limits<size_t>::max()) {
            break;
          }
          size_t seq_pos = aahash.get_pos() - (MINIMUM_CONSECUTIVE_HIT - 1);
          size_t temp_max_count = 0;
          for (auto &ID_count : id_to_count) {
            if (ID_count.second > temp_max_count) {
              temp_max_count = ID_count.second;
            }
          }

          for (auto &ID_count : id_to_count) {
            if (ID_count.second == temp_max_count) {
              id_set.insert(ID_count.first);
            }
          }

          std::unordered_map<uint32_t, std::set<uint32_t>> id_to_pos_set;
          for (size_t i = 0; i < miBf_IDs_snapshot.size(); ++i) {
            for (size_t j = 0; j < miBf_IDs_snapshot[i].size(); ++j) {
              if (id_set.find(miBf_IDs_snapshot[i][j]) != id_set.end()) {
                id_to_pos_set[miBf_IDs_snapshot[i][j]].insert(
                    miBf_pos_snapshot[i][j]);
              }
            }
          }

          aahash.roll();
          aahash2.roll();
          aahash3.roll();
          bool extend_block = true;
          while (extend_block &&
                 aahash.get_pos() != std::numeric_limits<size_t>::max()) {
            std::vector<uint32_t> ids_vec;
            std::vector<uint32_t> temp_pos_vec;
            if (mi_bf.bv_contains(aahash.hashes())) {
              auto temp_ID_pos = mi_bf.get_id(aahash.hashes());
              process_hashes(temp_ID_pos, id_set, id_to_pos_set, extend_block,
                             ids_vec, temp_pos_vec, mi_bf);
            } else if (mi_bf.bv_contains(aahash2.hashes())) {
              auto temp_ID_pos = mi_bf.get_id(aahash2.hashes());
              process_hashes(temp_ID_pos, id_set, id_to_pos_set, extend_block,
                             ids_vec, temp_pos_vec, mi_bf);

            } else if (mi_bf.bv_contains(aahash3.hashes())) {
              auto temp_ID_pos = mi_bf.get_id(aahash3.hashes());
              process_hashes(temp_ID_pos, id_set, id_to_pos_set, extend_block,
                             ids_vec, temp_pos_vec, mi_bf);

            } else {
              extend_block = false;
            }
            if (extend_block) {
              aahash.roll();
              aahash2.roll();
              aahash3.roll();
            }
          }

          for (auto &ID_pos_set : id_to_pos_set) {
            if (ID_pos_set.second.size() < 5) {
              continue;
            }
            frame_to_block_id_to_id_and_pos[frame][block_id] = std::make_pair(
                *ID_pos_set.second.begin(), *ID_pos_set.second.rbegin());
            id_to_FrameBlock_id_and_seq_pos[ID_pos_set.first].emplace(
                frame, block_id, seq_pos);

#pragma omp critical
            {

              if (id_to_count_across_all_frames.find(ID_pos_set.first) ==
                  id_to_count_across_all_frames.end()) {
                id_to_count_across_all_frames[ID_pos_set.first] =
                    ID_pos_set.second.size();
              } else {
                id_to_count_across_all_frames[ID_pos_set.first] +=
                    ID_pos_set.second.size();
              }
            }
            ++block_id;
          }

          miBf_IDs_snapshot.clear();
          miBf_pos_snapshot.clear();
          id_to_count.clear();
        }
      }

      if (debug_flag) {
        for (auto &ID_count : id_to_count_across_all_frames) {
          std::cerr << "ID: " << ID_count.first << " count: " << ID_count.second
                    << std::endl;
        }
        for (auto &ID_FrameBlock_id_seq_pos : id_to_FrameBlock_id_and_seq_pos) {
          std::cerr << "ID: " << ID_FrameBlock_id_seq_pos.first << std::endl;
          std::cerr
              << "name: "
              << miBf_ID_to_seq_ID_and_len[ID_FrameBlock_id_seq_pos.first].first
              << std::endl;
          for (auto &FrameBlock_id_seq_pos : ID_FrameBlock_id_seq_pos.second) {
            std::cerr << "frame: " << FrameBlock_id_seq_pos.frame
                      << " block_id: " << FrameBlock_id_seq_pos.block_id
                      << " seq_pos: "
                      << FrameBlock_id_seq_pos.query_start_in_prot_space
                      << std::endl;
          }
        }
        for (auto &FrameBlock_id_to_id_and_pos :
             frame_to_block_id_to_id_and_pos) {
          std::cerr << "frame: " << FrameBlock_id_to_id_and_pos.first
                    << std::endl;
          for (auto &block_id_to_id_and_pos :
               FrameBlock_id_to_id_and_pos.second) {
            std::cerr << "block_id: " << block_id_to_id_and_pos.first
                      << " smallest pos: "
                      << block_id_to_id_and_pos.second.first << " largest pos: "
                      << block_id_to_id_and_pos.second.second << std::endl;
          }
        }
      }

      std::string strand = "+";
      if (ori == 1) {
        strand = "-";
      }

      for (auto &ID_count : id_to_count_across_all_frames) {
        uint32_t miBf_ID = ID_count.first;
        if (miBf_ID_to_seq_ID_and_len.find(miBf_ID) ==
            miBf_ID_to_seq_ID_and_len.end()) {
          continue;
        }

        std::string seq_name = miBf_ID_to_seq_ID_and_len[miBf_ID].first;
        if (debug_flag) {
          std::cerr << "calculating for protein name: " << seq_name
                    << std::endl;
        }
        size_t complete_copies = 0;
        size_t incomplete_copies = 0;
        size_t expected_kmer_counts =
            miBf_ID_to_seq_ID_and_len[miBf_ID].second - kmer_size + 1;
        if (debug_flag) {
          std::cerr << "protein length: "
                    << miBf_ID_to_seq_ID_and_len[miBf_ID].second << std::endl;
        }
        size_t adjusted_kmer_counts = 0;
        size_t end_pos = 0;
        size_t frame = 3;
        size_t seq_start_in_nucleotide = 0;
        size_t seq_end_in_nucleotide = 0;
        size_t block_len = 0;
        size_t prev_block_len = 0;
        size_t block_start = 0;
        size_t prev_block_start = 0;
        std::vector<std::tuple<size_t, size_t>> start_end_pos_vec;
        std::vector<std::tuple<size_t, size_t>> start_end_pos_tar_vec;

        std::vector<std::reference_wrapper<const FrameBlock>> vec;
        for (auto &FrameBlock_id_and_seq_pos :
             id_to_FrameBlock_id_and_seq_pos[miBf_ID]) {
          vec.push_back(FrameBlock_id_and_seq_pos);
        }
        for (size_t ref_idx = 0; ref_idx < vec.size(); ++ref_idx) {
          const FrameBlock &FrameBlock_id_and_seq_pos = vec[ref_idx].get();
          if (frame == 3) {
            frame = FrameBlock_id_and_seq_pos.frame;
            seq_start_in_nucleotide =
                FrameBlock_id_and_seq_pos.query_start_in_prot_space * 3 + frame;
          } else {
            frame = FrameBlock_id_and_seq_pos.frame;
          }

          size_t block_id = FrameBlock_id_and_seq_pos.block_id;

          if (frame_to_block_id_to_id_and_pos[frame].find(block_id) !=
              frame_to_block_id_to_id_and_pos[frame].end()) {

            prev_block_start = block_start;
            block_start = FrameBlock_id_and_seq_pos.query_start_in_prot_space;
            prev_block_len = block_len;
            block_len =
                frame_to_block_id_to_id_and_pos[frame][block_id].second -
                frame_to_block_id_to_id_and_pos[frame][block_id].first + 1;

            if (start_end_pos_vec.empty()) {
              start_end_pos_vec.emplace_back(std::make_tuple(
                  FrameBlock_id_and_seq_pos.query_start_in_prot_space,
                  FrameBlock_id_and_seq_pos.query_start_in_prot_space +
                      block_len - 1));
              start_end_pos_tar_vec.emplace_back(std::make_tuple(
                  frame_to_block_id_to_id_and_pos[frame][block_id].first,
                  frame_to_block_id_to_id_and_pos[frame][block_id].second + 1));
            }

            if (end_pos == 0) {

              end_pos = frame_to_block_id_to_id_and_pos[frame][block_id].second;
              if (debug_flag) {
                std::cerr << "start end_pos: " << end_pos << std::endl;
                std::cerr << std::endl;
              }
              adjusted_kmer_counts = block_len;
            } else {
              if (end_pos <
                  frame_to_block_id_to_id_and_pos[frame][block_id].first) {
                start_end_pos_vec.emplace_back(std::make_tuple(
                    FrameBlock_id_and_seq_pos.query_start_in_prot_space,
                    FrameBlock_id_and_seq_pos.query_start_in_prot_space +
                        block_len - 1));
                start_end_pos_tar_vec.emplace_back(std::make_tuple(
                    frame_to_block_id_to_id_and_pos[frame][block_id].first,
                    frame_to_block_id_to_id_and_pos[frame][block_id].second +
                        1));
                if (frame_to_block_id_to_id_and_pos[frame][block_id].first -
                        end_pos >=
                    kmer_size) {
                  adjusted_kmer_counts +=
                      block_len + kmer_size -
                      1; // TODO adjust for kmer size overlap
                } else {
                  adjusted_kmer_counts +=
                      block_len +
                      frame_to_block_id_to_id_and_pos[frame][block_id].first -
                      end_pos - 1;
                }

                end_pos =
                    frame_to_block_id_to_id_and_pos[frame][block_id].second;
                if (debug_flag) {
                  std::cerr << "update end_pos: " << end_pos << std::endl;
                  std::cerr << std::endl;
                }

                if (ref_idx + 2 < vec.size()) {
                  size_t idx_offset =
                      look_ahead(vec, ref_idx, end_pos,
                                 frame_to_block_id_to_id_and_pos, max_offset);
                  if (idx_offset > 0) {
                    ref_idx += idx_offset;
                  }
                }
              } else {
                seq_end_in_nucleotide =
                    (prev_block_start + prev_block_len + kmer_size - 1) * 3 +
                    frame;
                double prev_score =
                    (double)adjusted_kmer_counts / (double)expected_kmer_counts;
                if (adjusted_kmer_counts > lower_bound * expected_kmer_counts) {
                  if (start_end_pos_vec.size() > 1) {
                    fill_in_gaps(start_end_pos_vec, start_end_pos_tar_vec,
                                 adjusted_kmer_counts, hash_num,
                                 rescue_kmer_size, sixframed_xlated_proteins,
                                 ori, kmer_size, miBf_ID, db_path_loc,
                                 mibf_prefix);
                  }
                  if (adjusted_kmer_counts > 0.95 * expected_kmer_counts) {
                    complete_copies++;
                  } else {
                    incomplete_copies++;
                  }
                }
                double score =
                    (double)adjusted_kmer_counts / (double)expected_kmer_counts;
                if (score > 1) {
                  score = 1;
                }
                {
                  if (strand == "-") {
                    auto temp = seq_start_in_nucleotide;
                    seq_start_in_nucleotide =
                        record.seq.size() - seq_end_in_nucleotide;
                    seq_end_in_nucleotide = record.seq.size() - temp;
                  }
#pragma omp critical
                  {
                    pre_gff_set.emplace(record.id, seq_start_in_nucleotide,
                                        seq_end_in_nucleotide, prev_score,
                                        strand, seq_name);
                    gff_set.emplace(record.id, seq_start_in_nucleotide,
                                    seq_end_in_nucleotide, score, strand,
                                    seq_name);
                  }
                }

                end_pos =
                    frame_to_block_id_to_id_and_pos[frame][block_id].second;
                adjusted_kmer_counts = block_len;
                seq_start_in_nucleotide =
                    FrameBlock_id_and_seq_pos.query_start_in_prot_space * 3 +
                    frame;
                start_end_pos_vec.clear();
                start_end_pos_tar_vec.clear();
                start_end_pos_vec.emplace_back(std::make_tuple(
                    FrameBlock_id_and_seq_pos.query_start_in_prot_space,
                    FrameBlock_id_and_seq_pos.query_start_in_prot_space +
                        block_len - 1));
                start_end_pos_tar_vec.emplace_back(std::make_tuple(
                    frame_to_block_id_to_id_and_pos[frame][block_id].first,
                    frame_to_block_id_to_id_and_pos[frame][block_id].second +
                        1));
                if (debug_flag) {
                  std::cerr << "new end_pos: " << end_pos << std::endl;
                }
              }
            }
          }
        }
        double prev_score =
            (double)adjusted_kmer_counts / (double)expected_kmer_counts;
        if (adjusted_kmer_counts > lower_bound * expected_kmer_counts) {
          if (start_end_pos_vec.size() > 1) {
            fill_in_gaps(start_end_pos_vec, start_end_pos_tar_vec,
                         adjusted_kmer_counts, hash_num, rescue_kmer_size,
                         sixframed_xlated_proteins, ori, kmer_size, miBf_ID,
                         db_path_loc, mibf_prefix);
          }
          if (adjusted_kmer_counts > 0.95 * expected_kmer_counts) {
            complete_copies++;
          } else {
            incomplete_copies++;
          }
        }
#pragma omp atomic
        seq_name_to_completeness[seq_name].complete_copies += complete_copies;

#pragma omp atomic
        seq_name_to_completeness[seq_name].incomplete_copies +=
            incomplete_copies;
        {

          const auto &last_FrameBlock_id_and_seq_pos =
              id_to_FrameBlock_id_and_seq_pos[miBf_ID].rbegin();
          frame = (*last_FrameBlock_id_and_seq_pos).frame;
          seq_end_in_nucleotide =
              ((*last_FrameBlock_id_and_seq_pos).query_start_in_prot_space +
               block_len + kmer_size - 1) *
                  3 +
              frame;
          double score =
              (double)adjusted_kmer_counts / (double)expected_kmer_counts;
          if (score > 1) {
            score = 1;
          }
          if (strand == "-") {
            auto temp = seq_start_in_nucleotide;
            seq_start_in_nucleotide = record.seq.size() - seq_end_in_nucleotide;
            seq_end_in_nucleotide = record.seq.size() - temp;
          }

#pragma omp critical
          {
            pre_gff_set.emplace(record.id, seq_start_in_nucleotide,
                                seq_end_in_nucleotide, prev_score, strand,
                                seq_name);
            gff_set.emplace(record.id, seq_start_in_nucleotide,
                            seq_end_in_nucleotide, score, strand, seq_name);
          }
        }
      }
    }
  }

  std::unordered_map<std::string, double> seq_name_to_score;
  for (auto &seq_name_completeness : seq_name_to_completeness) {
    output_file << seq_name_completeness.first << "\t"
                << seq_name_completeness.second.complete_copies << "\t"
                << seq_name_completeness.second.incomplete_copies << "\t"
                << seq_name_completeness.second.expected_kmer_counts << "\t"
                << seq_name_completeness.second.highest_adjusted_kmer_counts
                << std::endl;
    seq_name_to_score[seq_name_completeness.first] = 0;
  }

  for (auto &gff : gff_set) {
    gff_file << gff.query_name << "\t"
             << "."
             << "\t"
             << "gene"
             << "\t" << gff.hit_pos_start << "\t" << gff.hit_pos_end << "\t"
             << gff.score << "\t" << gff.strand << "\t"
             << "0"
             << "\t"
             << "ID=" << gff.hit_name << std::endl;
    if (seq_name_to_score[gff.hit_name] < gff.score) {
      seq_name_to_score[gff.hit_name] = gff.score;
    }
  }
  std::vector<double> scores_vec;
  for (auto &seq_name_score : seq_name_to_score) {
    scores_vec.push_back(seq_name_score.second);
  }

  double result = 0.0;

  auto ecdf = empirical_cumulative_distribution_function(std::move(scores_vec));
  result = (1 - trapezoidal(ecdf, 0.0, 1.0)) * 100;

  std::cout << "aaKomp score for " << input_file << " is: " << result
            << std::endl;

  // write result to text file
  std::string result_file = output_prefix + "_score.txt";
  std::ofstream result_out(result_file);
  result_out << result << std::endl;

  if (debug_flag) {
    for (auto &gff : pre_gff_set) {
      pre_gff_file << gff.query_name << "\t"
                   << "."
                   << "\t"
                   << "gene"
                   << "\t" << gff.hit_pos_start << "\t" << gff.hit_pos_end
                   << "\t" << gff.score << "\t" << gff.strand << "\t"
                   << "0"
                   << "\t"
                   << "ID=" << gff.hit_name << std::endl;
    }
  }

  return 0;
}
