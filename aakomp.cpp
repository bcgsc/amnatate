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
#include <argparse/argparse.hpp>
#include <btllib/aahash.hpp>
#include <btllib/seq.hpp>
#include <btllib/seq_reader.hpp>
#include <btllib/mi_bloom_filter.hpp>
#include <Sequence/Translate.hpp>

static constexpr uint32_t HASH_ID_SHIFT = 32;
static constexpr uint32_t HASH_POS_MASK = 0xFFFFFFFF;

struct FrameBlock {
    size_t frame;
    size_t block_id;
    size_t query_start_in_prot_space;

    FrameBlock(size_t frame, size_t block_id, size_t query_start)
        : frame(frame), block_id(block_id), query_start_in_prot_space(query_start) {}
};

struct GFFEntry {
    std::string query_name;
    size_t hit_pos_start;
    size_t hit_pos_end;
    double score;
    std::string strand;
    std::string hit_name;

    GFFEntry(const std::string& query, size_t start, size_t end, double score,
             const std::string& strand, const std::string& hit)
        : query_name(query), hit_pos_start(start), hit_pos_end(end),
          score(score), strand(strand), hit_name(hit) {}
};

struct GFFEntryComparator {
    bool operator()(const GFFEntry& lhs, const GFFEntry& rhs) const noexcept {
        if (lhs.query_name == rhs.query_name)
            return lhs.hit_pos_start < rhs.hit_pos_start;
        return lhs.query_name < rhs.query_name;
    }
};

struct FrameBlockComparator {
    bool operator()(const FrameBlock& a, const FrameBlock& b) const noexcept {
        return a.query_start_in_prot_space < b.query_start_in_prot_space;
    }
};

size_t look_ahead(
    const std::vector<std::reference_wrapper<const FrameBlock>>& blocks,
    size_t current_index,
    size_t query_end,
    const std::unordered_map<size_t, std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>>>& frame_index,
    size_t offset)
{
    if (current_index + 1 >= blocks.size()) {
        return 0;
    }

    const FrameBlock& next = blocks[current_index + 1].get();
    auto next_block_id = next.block_id;
    auto next_frame = next.frame;

    if (query_end < frame_index.at(next_frame).at(next_block_id).first) {
        return 0;
    }

    for (size_t i = 2; i <= offset; ++i) {
        if (current_index + i >= blocks.size()) {
            return 0;
        }

        const FrameBlock& block = blocks[current_index + i].get();
        auto block_id = block.block_id;
        auto frame = block.frame;

        if (query_end < frame_index.at(frame).at(block_id).first) {
            return i - 1;
        }
    }

    return 0;
}


void process_hashes(
    const std::vector<uint64_t>& temp_ID_pos,
    std::unordered_set<uint32_t>& id_set,
    std::unordered_map<uint32_t, std::set<uint32_t>>& id_to_pos_set,
    bool& extend_block,
    std::vector<uint32_t>& ids_vec,
    std::vector<uint32_t>& temp_pos_vec,
    const btllib::MIBloomFilter<uint64_t>& mi_bf)
{
    bool found = false;

    for (const auto& id : id_set) {
        for (const auto& id_pos : temp_ID_pos) {
            auto demasked = id_pos & mi_bf.ANTI_MASK;
            if (id == (demasked >> HASH_ID_SHIFT)) {
                found = true;
                break;
            }
        }
        if (found) break;
    }

    if (found) {
        for (const auto& id_pos : temp_ID_pos) {
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
                const auto& pos_set = id_to_pos_set[id];
                if (!pos_set.empty() && pos == *pos_set.rbegin() + 1) {
                    id_to_pos_set[id].insert(pos);
                }
            } else {
                new_ids.push_back(id);
                new_pos.push_back(pos);
            }
        }

        std::unordered_set<uint32_t> new_id_set(new_ids.begin(), new_ids.end());
        for (const auto& id : new_id_set) {
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
        for (const auto& id_pos : temp_ID_pos) {
            if (id_pos < mi_bf.MASK) {
                saturated = false;
                break;
            }
        }

        if (!saturated) {
            extend_block = false;
        } else {
            extend_block = true;
            for (auto& [id, pos_set] : id_to_pos_set) {
                if (!pos_set.empty()) {
                    pos_set.insert(*pos_set.rbegin() + 1);
                }
            }
        }
    }
}

struct CustomComparator {
    bool operator()(const std::tuple<size_t, size_t, size_t>& a, const std::tuple<size_t, size_t, size_t>& b) const {
        return std::get<2>(a) < std::get<2>(b);
    }
};

size_t calc_optimal_size(size_t entries, unsigned hash_num, double occupancy)
{
    size_t non64ApproxVal =
        size_t(-double(entries) * double(hash_num) / log(1.0 - occupancy));
    return non64ApproxVal + (64 - non64ApproxVal % 64);
}

std::vector<std::string> sixframe_translate(const std::string &dna)
{
    std::vector<std::string> protein;
    std::string rev_dna = btllib::get_reverse_complement(dna);
    protein.push_back(Sequence::Translate(dna.begin(), dna.end()));
    protein.push_back(Sequence::Translate(dna.begin() + 1, dna.end()));
    protein.push_back(Sequence::Translate(dna.begin() + 2, dna.end()));
    protein.push_back(Sequence::Translate(rev_dna.begin(), rev_dna.end()));
    protein.push_back(Sequence::Translate(rev_dna.begin() + 1, rev_dna.end()));
    protein.push_back(Sequence::Translate(rev_dna.begin() + 2, rev_dna.end()));
    return protein;
}


void fill_in_gaps(std::vector<std::tuple<size_t, size_t>>& start_end_pos_vec,
                  std::vector<std::tuple<size_t, size_t>>& start_end_pos_in_tar_space_vec,
                  size_t& adjusted_kmer_counts,
                  size_t hash_num, size_t rescue_kmer_size,
                  const std::vector<std::string>& sixframed_xlated_proteins,
                  size_t ori, size_t kmer_size,
                  uint32_t miBf_ID, const std::string& db_path) {

    // Sort start_end_pos_vec by start position
    std::sort(start_end_pos_vec.begin(), start_end_pos_vec.end(), [](const std::tuple<size_t, size_t>& a, const std::tuple<size_t, size_t>& b) {
        return std::get<0>(a) < std::get<0>(b);
    });


    // Find all the gaps between the start and end positions larger than kmer_size
    std::vector<std::tuple<size_t, size_t>> gap_vec;
    for (size_t i = 0; i < start_end_pos_vec.size() - 1; ++i) {
        if ((int)std::get<0>(start_end_pos_vec[i + 1]) - (int)std::get<1>(start_end_pos_vec[i]) > (int)kmer_size) {
            gap_vec.emplace_back(std::make_tuple(std::get<1>(start_end_pos_vec[i]), std::get<0>(start_end_pos_vec[i + 1])));
        }
    }

    if (gap_vec.empty()) {
        return;
    }


    // Find gaps in the target space
    std::vector<std::tuple<size_t, size_t>> gap_in_tar_space_vec;
    for (size_t i = 0; i < start_end_pos_in_tar_space_vec.size() - 1; ++i) {
        if ((int)std::get<0>(start_end_pos_in_tar_space_vec[i + 1]) - (int)(std::get<1>(start_end_pos_in_tar_space_vec[i]) + 5) > (int)kmer_size) {
            gap_in_tar_space_vec.emplace_back(std::make_tuple(std::get<1>(start_end_pos_in_tar_space_vec[i]) + 5, std::get<0>(start_end_pos_in_tar_space_vec[i + 1])));
        } else {
            //insert a dummy gap
            gap_in_tar_space_vec.emplace_back(std::make_tuple(0, 0));
        }
    }

    if (gap_in_tar_space_vec.empty()) {
        return;
    }

    // Create a vector of unordered sets for each gap in the target space
    std::vector<std::unordered_set<size_t>> gap_index_sets(gap_in_tar_space_vec.size());
    // make sure that the index sets are empty
    for (size_t i = 0; i < gap_index_sets.size(); ++i) {
        gap_index_sets[i].clear();
    }

    // Populate the gap index sets
    for (size_t g = 0; g < gap_in_tar_space_vec.size(); ++g) {
        for (size_t i = std::get<0>(gap_in_tar_space_vec[g]); i < std::get<1>(gap_in_tar_space_vec[g]); ++i) {
            gap_index_sets[g].insert(i);
        }
    }

    if (gap_index_sets.empty()) {
        return;
    }

    std::string small_mibf_path = db_path + "/" + std::to_string(miBf_ID) + ".mibf";
    btllib::MIBloomFilter<uint64_t> small_mi_bf(small_mibf_path);

    // Collect k-mers for all frames and levels
    std::vector<std::vector<std::tuple<size_t, size_t>>> kmer_pos_per_frame_and_level(3);

    for (size_t frame = 0; frame < 3; ++frame) {
        for (size_t current_level = 1; current_level <= 3; ++current_level) {
            // Iterate over each gap
            for (size_t g = 0; g < gap_vec.size(); ++g) {
                size_t gap_start = std::get<0>(gap_vec[g]);
                size_t gap_end = std::get<1>(gap_vec[g]);
                // check if gap in target space is 0,0, if so skip
                if (std::get<0>(gap_in_tar_space_vec[g]) == 0 && std::get<1>(gap_in_tar_space_vec[g]) == 0) {
                    continue;
                }

                btllib::AAHash aahash(sixframed_xlated_proteins[frame + ori * 3], hash_num, rescue_kmer_size, current_level, gap_start - 1);
                aahash.roll();

                while (aahash.get_pos() <= gap_end + 1) {
                    if (small_mi_bf.bv_contains(aahash.hashes())) {
                        auto temp_ID_pos = small_mi_bf.get_id(aahash.hashes());
                        for (auto& ID_pos : temp_ID_pos) {
                            auto pos = ID_pos & 0xFFFFFFFF;
                            if (gap_index_sets[g].find(pos) != gap_index_sets[g].end()) {
                                kmer_pos_per_frame_and_level[frame].emplace_back(pos, aahash.get_pos());
                                //gap_index_sets[g].erase(pos);  // Remove the index once it's matched
                                //break;
                            }
                        }
                    }
                    aahash.roll();
                }
            }
        }
    }

    // Perform Dynamic Programming (DP) to find the longest subsequence of valid k-mers with monotonically increasing positions
    std::vector<std::tuple<size_t, size_t>> all_kmers;
    std::set<std::tuple<size_t, size_t>> all_kmers_set;

    for (size_t frame = 0; frame < 3; ++frame) {
        for (const auto& kmer : kmer_pos_per_frame_and_level[frame]) {
            if (all_kmers_set.find(kmer) == all_kmers_set.end()) {
                all_kmers.push_back(kmer);
                all_kmers_set.insert(kmer);
            }

        }
    }

    // Sort k-mers by their sequence positions
    std::sort(all_kmers.begin(), all_kmers.end(), [](const std::tuple<size_t, size_t>& a, const std::tuple<size_t, size_t>& b) {
        return std::get<1>(a) < std::get<1>(b);
    });

    // Dynamic Programming to find the longest subsequence of valid k-mers with monotonically increasing positions
    std::vector<size_t> dp(all_kmers.size(), 1);  // dp[i] = length of longest subsequence ending at i

    for (size_t i = 1; i < all_kmers.size(); ++i) {
        for (size_t j = 0; j < i; ++j) {
            if (std::get<1>(all_kmers[i]) > std::get<1>(all_kmers[j])) {
                dp[i] = std::max(dp[i], dp[j] + 1);
            }
        }
    }

    if (!all_kmers.empty()) {
        size_t prev_pos = std::get<1>(all_kmers[0]);  // safe now
        for (size_t i = 1; i < all_kmers.size(); ++i) {
            size_t current_pos = std::get<1>(all_kmers[i]);
            size_t gap = current_pos - prev_pos;

            if (gap > 3) {
                adjusted_kmer_counts += 3;
            } else {
                adjusted_kmer_counts += gap;
            }

            prev_pos = current_pos;
        }

        adjusted_kmer_counts += 3;  // final adjustment
    }

    return;
}

bool explore_frame(btllib::MIBloomFilter<uint64_t> &mi_bf, btllib::AAHash &aahash, std::deque<std::vector<uint32_t>> &miBf_IDs_snapshot, std::deque<std::vector<uint32_t>> &miBf_pos_snapshot, std::unordered_map<uint32_t, size_t> &id_to_count)
{
    // check size of miBf_IDs_snapshot and miBf_pos_snapshot
    //  if size is more than 10, pop front
    std::unordered_set<uint32_t> id_set;
    if (miBf_IDs_snapshot.size() >= 5)
    {
        // insert id into id_set before removing
        for (size_t i = 0; i < miBf_IDs_snapshot.front().size(); ++i)
        {
            id_set.insert(miBf_IDs_snapshot.front()[i]);
        }
        // remove id from id_to_count
        for (auto &ID : id_set)
        {
            id_to_count[ID]--;
            if (id_to_count[ID] == 0)
            {
                id_to_count.erase(ID);
            }
        }
        id_set.clear();
        miBf_IDs_snapshot.pop_front();
        miBf_pos_snapshot.pop_front();
    }
    // push back new vector
    // consider circular linked list
    miBf_IDs_snapshot.emplace_back(std::vector<uint32_t>());
    miBf_pos_snapshot.emplace_back(std::vector<uint32_t>());


    if (!mi_bf.bv_contains(aahash.hashes()))
    {
        // clear snapshots and id_to_count
        miBf_IDs_snapshot.clear();
        miBf_pos_snapshot.clear();
        id_to_count.clear();
        return false;
    }

    // query mibf and insert into both deques
    auto temp_ID_pos = mi_bf.get_id(aahash.hashes());
    for (auto &ID_pos : temp_ID_pos)
    {
        auto demasked_ID_pos = ID_pos  & mi_bf.ANTI_MASK;
        miBf_IDs_snapshot.back().push_back(demasked_ID_pos >> 32);
        miBf_pos_snapshot.back().push_back(demasked_ID_pos & 0xFFFFFFFF);
    }

    for (size_t j = 0; j < miBf_IDs_snapshot.back().size(); ++j)
    {
        id_set.insert(miBf_IDs_snapshot.back()[j]);
    }
    for (auto &ID : id_set)
    {
        if (id_to_count.find(ID) == id_to_count.end())
        {
            id_to_count[ID] = 1;
        }
        else
        {
            id_to_count[ID]++;
        }
    }
    id_set.clear();

    if (miBf_IDs_snapshot.size() < 5)
    {
        return false;
    }

    uint32_t temp_mibf_ID = 0;
    size_t temp_max_count = 0;
    // iterate id_to_count and find the ID with the highest count
    for (auto &ID_count : id_to_count)
    {
        if (ID_count.second > temp_max_count)
        {
            temp_mibf_ID = ID_count.first;
            temp_max_count = ID_count.second;
        }
    }
    if (temp_mibf_ID == 0 || temp_max_count < 5)
    {
        return false;
    }
    else
    {
        std::vector<uint32_t> ids_to_check;
        for (auto &ID_count : id_to_count)
        {
            if (ID_count.second == temp_max_count)
            {
                ids_to_check.push_back(ID_count.first);
            }
        }
        // check to see if the ID with the highest count has consecutive positions
        // if not, return false
        // if yes, return true
        for (auto &ID : ids_to_check)
        {
            std::set<uint32_t> temp_pos_set;
            for (size_t i = 0; i < miBf_IDs_snapshot.size(); ++i)
            {
                for (size_t j = 0; j < miBf_IDs_snapshot[i].size(); ++j)
                {
                    if (miBf_IDs_snapshot[i][j] == ID)
                    {
                        temp_pos_set.insert(miBf_pos_snapshot[i][j]);
                    }
                }
            }

            // iterate through temp_pos_set and check if it is incrementing by 1
            uint32_t prev_pos = 0;
            size_t counter = 0;
            bool init = false;
            for (auto &pos : temp_pos_set)
            {
                if (!init)
                {
                    prev_pos = pos;
                    init = true;
                    continue;
                }
                if (pos - prev_pos == 1)
                {
                    ++counter;
                }
                else
                {
                    counter = 0;
                }
                prev_pos = pos;

                if (counter >= 4)
                {
                    return true;
                }
            }
        }
    }

    return false;
}


// main function that takes in command line arguments using getopt
    // declare variables
int main(int argc, char* argv[]) {
    argparse::ArgumentParser program("aaKomp");

    program.add_argument("--help")
        .help("Print this help message")
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
    
    program.add_argument("-v", "--verbose")
        .help("Verbose output")
        .default_value(false)
        .implicit_value(true);

    program.add_argument("--debug")
        .help("Debug output")
        .default_value(false)
        .implicit_value(true);
    
    program.add_argument("-m", "--mibf_path")
        .help("MIBF file path")
        .default_value(std::string(""));
    
    program.add_argument("-h", "--hash")
        .help("Number of hash functions")
        .default_value(uint8_t(9))
        .scan<'u', uint8_t>();
    
    program.add_argument("-k", "--kmer")
        .help("K-mer size")
        .default_value(uint8_t(9))
        .scan<'u', uint8_t>();

    
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

    try {
        program.parse_args(argc, argv);
    } catch (const std::exception& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program;
        return 1;
    }

    // Extract values
    bool help_flag = program.get<bool>("--help");
    bool verbose_flag = program.get<bool>("--verbose");
    bool debug_flag = program.get<bool>("--debug");
    size_t threads = program.get<size_t>("--threads");
    std::string input_file = program.get<std::string>("--input");
    std::string reference_path = program.get<std::string>("--reference");
    std::string output_prefix = program.get<std::string>("--output");
    std::string mibf_path = program.get<std::string>("--mibf_path");
    uint8_t hash_num = program.get<uint8_t>("--hash");
    uint8_t kmer_size = program.get<uint8_t>("--kmer");
    size_t rescue_kmer_size = program.get<size_t>("--rescue_kmer");
    double lower_bound = program.get<double>("--lower_bound");
    size_t max_offset = program.get<size_t>("--max_offset");

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

    // Print error messages if required arguments are missing
    if (input_file.empty()) {
        std::cerr << "Input file is required. Use -h or --help for more information." << std::endl;
        return 1;
    }

    if (reference_path.empty()) {
        std::cerr << "Reference path is required. Use -h or --help for more information." << std::endl;
        return 1;
    }

    if (threads == 0) {
        std::cerr << "Threads must be greater than 0. Use -h or --help for more information." << std::endl;
        return 1;
    }
    
    // Print parsed arguments
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


    if (verbose_flag)
    {
        std::cerr << "Reading reference file: " << reference_path << std::endl;
    }


    if (verbose_flag)
    {
        std::cerr << "Creating seq_id to ID table" << std::endl;
    }

    std::unordered_map<std::string, uint32_t> seq_ID_to_miBf_ID;
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
            ++miBf_ID;
        }
    }
    if (verbose_flag)
    {
        std::cerr << "Reading miBF" << std::endl;
    }
    //strip suffix from reference path
    //std::string reference path_no_suffix = reference path.substr(0, reference path.find_last_of("."));
    
  
    auto sTime = omp_get_wtime();
    btllib::MIBloomFilter<uint64_t> mi_bf(mibf_path);
    auto sTime2 = omp_get_wtime();
    if (verbose_flag)
    {
    std::cerr << "finished reading miBf" << std::endl;
    std::cerr << "in " << std::setprecision(4) << std::fixed << sTime2 - sTime
              << "\n";
    }
    
    

    std::vector<std::ofstream> output_files(3);
    std::vector<std::ofstream> gff_files(3);
    std::vector<std::ofstream> pre_gff_files(3);

    // Open each output file stream with a unique filename based on output_prefix
    std::string filename = output_prefix + ".results.tsv";
    output_files[0].open(filename);
    output_files[0] << "name\tcomplete copies\tincomplete copies\texpected k-mer counts\thighest adjusted incomplete k-mer hits" << std::endl;
    gff_files[0].open(output_prefix + ".gff");
    gff_files[0] << "##gff-version 3" << std::endl;
    if (debug_flag) {
        pre_gff_files[0].open(output_prefix + ".pre.gff");
        pre_gff_files[0] << "##gff-version 3" << std::endl;
    }


    // Create and insert three different instances of gff_set into the vector
    std::vector<std::set<GFFEntry, GFFEntryComparator>> gff_set_vector;
    for (int i = 0; i < 3; ++i) {
        gff_set_vector.emplace_back(GFFEntryComparator());
    }

   //std::vector<std::set<std::tuple<std::string, size_t, size_t, double, std::string, std::string>, decltype(gff_comparator)>> pre_gff_set_vector;
    std::vector<std::set<GFFEntry, GFFEntryComparator>> pre_gff_set_vector;
    for (int i = 0; i < 3; ++i) {
        pre_gff_set_vector.emplace_back(GFFEntryComparator());
    }


    btllib::SeqReader reader(input_file, btllib::SeqReader::Flag::LONG_MODE);
    if (verbose_flag)
    {
        std::cerr << "Reading input file: " << input_file << std::endl;
    }
    struct completeness_struct
    {
        size_t complete_copies = 0;
        size_t incomplete_copies = 0;
        size_t expected_kmer_counts = 0;
        size_t highest_adjusted_kmer_counts = 0;
    };
    std::vector<std::unordered_map<std::string, completeness_struct>> seq_name_to_completeness_vec(3);

    // Populate each map in the vector with empty entries
    for (auto &seq_name_to_completeness : seq_name_to_completeness_vec) {
        for (const auto &seq_ID : seq_ID_to_miBf_ID) {
            seq_name_to_completeness[seq_ID.first] = completeness_struct();
        }
    }

#pragma omp parallel num_threads(threads)
    for (const auto record : reader)
    {
        std::vector<std::string> sixframed_xlated_proteins = sixframe_translate(record.seq);
        for (size_t ori = 0; ori < 2; ++ori)
        //for (size_t ori = 0; ori < 2; ++ori) //TODO
        {
            // frame to block id to id and smallest pos and largest pos
            std::vector<std::unordered_map<size_t, std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>>>> frame_to_block_id_to_id_and_pos_vec(3);
            // id to count across all frames sorted by count largest to smallest
            std::vector<std::map<uint32_t, size_t, std::greater<size_t>>> id_to_count_across_all_frames_vec(3);
            // id to set of frame and block id and seq pos, set is sorted by seq pos
            std::vector<std::unordered_map<uint32_t, std::set<FrameBlock, FrameBlockComparator>>> id_to_FrameBlock_id_and_seq_pos_vec(3);
            for (size_t curr_lvl = 1; curr_lvl <= 1; ++curr_lvl) {
                auto& frame_to_block_id_to_id_and_pos = frame_to_block_id_to_id_and_pos_vec[curr_lvl - 1];
                auto& id_to_count_across_all_frames = id_to_count_across_all_frames_vec[curr_lvl - 1];
                auto& id_to_FrameBlock_id_and_seq_pos = id_to_FrameBlock_id_and_seq_pos_vec[curr_lvl - 1];
                auto& gff_set = gff_set_vector[curr_lvl - 1];
                auto& pre_gff_set = pre_gff_set_vector[curr_lvl - 1];
                auto& seq_name_to_completeness = seq_name_to_completeness_vec[curr_lvl - 1];
                for (size_t frame = 0; frame < 3; ++frame)
                {
                    btllib::AAHash aahash(sixframed_xlated_proteins[frame + ori * 3], hash_num, kmer_size, curr_lvl);
                    btllib::AAHash aahash2(sixframed_xlated_proteins[frame + ori * 3], hash_num, kmer_size, 2);
                    btllib::AAHash aahash3(sixframed_xlated_proteins[frame + ori * 3], hash_num, kmer_size, 3);
                    aahash.roll();
                    aahash2.roll();
                    aahash3.roll();
                    std::deque<std::vector<uint32_t>> miBf_IDs_snapshot;
                    std::deque<std::vector<uint32_t>> miBf_pos_snapshot;
                    std::unordered_map<uint32_t, size_t> id_to_count;
                    
                    
                    size_t block_id = 0;
                    std::unordered_set<uint32_t> id_set;
                    while (aahash.get_pos() != std::numeric_limits<size_t>::max())
                    {
                        while (!explore_frame(mi_bf, aahash, miBf_IDs_snapshot, miBf_pos_snapshot, id_to_count) && aahash.get_pos() != std::numeric_limits<size_t>::max())
                        {
                            aahash.roll();
                            aahash2.roll();
                            aahash3.roll();
                        }
                        if (aahash.get_pos() == std::numeric_limits<size_t>::max())
                        {
                            break;
                        }
                        size_t seq_pos = aahash.get_pos() - 4; //TODO
                        // find the largest count in id_to_count
                        size_t temp_max_count = 0;
                        for (auto &ID_count : id_to_count)
                        {
                            if (ID_count.second > temp_max_count)
                            {
                                temp_max_count = ID_count.second;
                            }
                        }
                        // insert id into id_set if count is equal to temp_max_count
                        for (auto &ID_count : id_to_count)
                        {
                            if (ID_count.second == temp_max_count)
                            {
                                id_set.insert(ID_count.first);
                            }
                        }

                        std::unordered_map<uint32_t, std::set<uint32_t>> id_to_pos_set;
                        for (size_t i = 0; i < miBf_IDs_snapshot.size(); ++i)
                        {
                            for (size_t j = 0; j < miBf_IDs_snapshot[i].size(); ++j)
                            {
                                if (id_set.find(miBf_IDs_snapshot[i][j]) != id_set.end())
                                {
                                    id_to_pos_set[miBf_IDs_snapshot[i][j]].insert(miBf_pos_snapshot[i][j]);
                                }
                            }
                        }



                        aahash.roll();
                        aahash2.roll();
                        aahash3.roll();
                        bool extend_block = true;
                        while (extend_block && aahash.get_pos() != std::numeric_limits<size_t>::max())
                        {
                            std::vector<uint32_t> ids_vec;
                            std::vector<uint32_t> temp_pos_vec;
                            if (mi_bf.bv_contains(aahash.hashes()))
                            {
                                auto temp_ID_pos = mi_bf.get_id(aahash.hashes());
                                process_hashes(temp_ID_pos, id_set, id_to_pos_set, extend_block, ids_vec, temp_pos_vec, mi_bf);
                            } else if (mi_bf.bv_contains(aahash2.hashes())) {
                                auto temp_ID_pos = mi_bf.get_id(aahash2.hashes());
                                process_hashes(temp_ID_pos, id_set, id_to_pos_set, extend_block, ids_vec, temp_pos_vec, mi_bf);

                            } else if (mi_bf.bv_contains(aahash3.hashes())) {
                                auto temp_ID_pos = mi_bf.get_id(aahash3.hashes());
                                process_hashes(temp_ID_pos, id_set, id_to_pos_set, extend_block, ids_vec, temp_pos_vec, mi_bf);

                            } else {
                                extend_block = false;
                            }
                            if (extend_block)
                            {
                                aahash.roll();
                                aahash2.roll();
                                aahash3.roll();
                            }
                        }


                        // log the block id, id, and smallest and largest pos
                        for (auto &ID_pos_set : id_to_pos_set)
                        {
                            if (ID_pos_set.second.size() < 5)
                            {
                                continue;
                            }
                            frame_to_block_id_to_id_and_pos[frame][block_id] = std::make_pair(*ID_pos_set.second.begin(), *ID_pos_set.second.rbegin());
                            id_to_FrameBlock_id_and_seq_pos[ID_pos_set.first].emplace(frame, block_id, seq_pos);

    #pragma omp critical
                            {

                                if (id_to_count_across_all_frames.find(ID_pos_set.first) == id_to_count_across_all_frames.end())
                                {
                                    id_to_count_across_all_frames[ID_pos_set.first] = ID_pos_set.second.size();
                                }
                                else
                                {
                                    id_to_count_across_all_frames[ID_pos_set.first] += ID_pos_set.second.size();
                                }
                            }
                            ++block_id;
                        }

                        

                        // clear miBf_IDs_snapshot, miBf_pos_snapshot, and id_to_count
                        miBf_IDs_snapshot.clear();
                        miBf_pos_snapshot.clear();
                        id_to_count.clear();
                    }
                }
                // iterate through id_to_count_across_all_frames and log the completeness
                // print id_to_count_across_all_frames
                if (verbose_flag) {
                    for (auto &ID_count : id_to_count_across_all_frames)
                    {
                        std::cerr << "ID: " << ID_count.first << " count: " << ID_count.second << std::endl;
                    }
                    // print id_to_FrameBlock_id_and_seq_pos
                    for (auto &ID_FrameBlock_id_seq_pos : id_to_FrameBlock_id_and_seq_pos)
                    {
                        std::cerr << "ID: " << ID_FrameBlock_id_seq_pos.first << std::endl;
                        std::cerr << "name: " << miBf_ID_to_seq_ID_and_len[ID_FrameBlock_id_seq_pos.first].first << std::endl;
                        for (auto &FrameBlock_id_seq_pos : ID_FrameBlock_id_seq_pos.second)
                        {
                            std::cerr << "frame: " << FrameBlock_id_seq_pos.frame << " block_id: " << FrameBlock_id_seq_pos.block_id << " seq_pos: " << FrameBlock_id_seq_pos.query_start_in_prot_space << std::endl;
                        }
                    }
                    // print frame_to_block_id_to_id_and_pos
                    for (auto &FrameBlock_id_to_id_and_pos : frame_to_block_id_to_id_and_pos)
                    {
                        std::cerr << "frame: " << FrameBlock_id_to_id_and_pos.first << std::endl;
                        for (auto &block_id_to_id_and_pos : FrameBlock_id_to_id_and_pos.second)
                        {
                            std::cerr << "block_id: " << block_id_to_id_and_pos.first << " smallest pos: " << block_id_to_id_and_pos.second.first << " largest pos: " << block_id_to_id_and_pos.second.second << std::endl;
                        }
                    }
                }
                std::string strand = "+";
                if (ori == 1)
                {
                    strand = "-";
                }

                for (auto &ID_count : id_to_count_across_all_frames)
                {
                    uint32_t miBf_ID = ID_count.first;
                    // if mibf id is not in miBf_ID_to_seq_ID_and_len, continue
                    if (miBf_ID_to_seq_ID_and_len.find(miBf_ID) == miBf_ID_to_seq_ID_and_len.end())
                    {
                        continue;
                    }
                    
                    // iterate through id to frame block id and seq pos and log the completeness
                    std::string seq_name = miBf_ID_to_seq_ID_and_len[miBf_ID].first;
                    if (verbose_flag) {
                        std::cerr << "calculating for protein name: " << seq_name << std::endl;
                    }
                    size_t complete_copies = 0;
                    size_t incomplete_copies = 0;
                    size_t expected_kmer_counts = miBf_ID_to_seq_ID_and_len[miBf_ID].second - kmer_size + 1;
                    if (verbose_flag) {
                        std::cerr << "protein length: " << miBf_ID_to_seq_ID_and_len[miBf_ID].second << std::endl;
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
                    for (auto &FrameBlock_id_and_seq_pos : id_to_FrameBlock_id_and_seq_pos[miBf_ID])
                    {
                        vec.push_back(FrameBlock_id_and_seq_pos);
                    }
                    for (size_t ref_idx = 0; ref_idx < vec.size(); ++ref_idx)
                    {
                        const FrameBlock &FrameBlock_id_and_seq_pos = vec[ref_idx].get();

                        if (verbose_flag) {
                            /*std::cerr << "frame: " << std::get<0>(FrameBlock_id_and_seq_pos) << " block_id: " << std::get<1>(FrameBlock_id_and_seq_pos) << " seq_pos: " << std::get<2>(FrameBlock_id_and_seq_pos) << std::endl;
                            std::cerr << "curr end_pos: " << end_pos << std::endl;
                            std::cerr << "next start_pos: " << frame_to_block_id_to_id_and_pos[std::get<0>(FrameBlock_id_and_seq_pos)][std::get<1>(FrameBlock_id_and_seq_pos)].first << std::endl;
                            std::cerr << std::endl;*/
                        }

                        if (frame == 3) {
                            frame = FrameBlock_id_and_seq_pos.frame;
                            seq_start_in_nucleotide = FrameBlock_id_and_seq_pos.query_start_in_prot_space * 3 + frame;
                        } else {
                            frame = FrameBlock_id_and_seq_pos.frame;
                        }

                        size_t block_id = FrameBlock_id_and_seq_pos.block_id;



                        if (frame_to_block_id_to_id_and_pos[frame].find(block_id) != frame_to_block_id_to_id_and_pos[frame].end())
                        {

                        
                            prev_block_start = block_start;
                            block_start = FrameBlock_id_and_seq_pos.query_start_in_prot_space;
                            prev_block_len = block_len;
                            block_len = frame_to_block_id_to_id_and_pos[frame][block_id].second - frame_to_block_id_to_id_and_pos[frame][block_id].first + 1;
                            


                            if (start_end_pos_vec.empty()) {
                                start_end_pos_vec.emplace_back(std::make_tuple(FrameBlock_id_and_seq_pos.query_start_in_prot_space, FrameBlock_id_and_seq_pos.query_start_in_prot_space + block_len - 1));\
                                start_end_pos_tar_vec.emplace_back(std::make_tuple(frame_to_block_id_to_id_and_pos[frame][block_id].first, frame_to_block_id_to_id_and_pos[frame][block_id].second + 1));
                            }

                            if (end_pos == 0)
                            {
                                
                                end_pos = frame_to_block_id_to_id_and_pos[frame][block_id].second;
                                if (verbose_flag){
                                    std::cerr << "start end_pos: " << end_pos << std::endl;
                                    std::cerr << std::endl;
                                }
                                adjusted_kmer_counts = block_len;
                            }
                            else
                            {
                                if (end_pos < frame_to_block_id_to_id_and_pos[frame][block_id].first)
                                {
                                    start_end_pos_vec.emplace_back(std::make_tuple(FrameBlock_id_and_seq_pos.query_start_in_prot_space, FrameBlock_id_and_seq_pos.query_start_in_prot_space + block_len - 1));
                                    start_end_pos_tar_vec.emplace_back(std::make_tuple(frame_to_block_id_to_id_and_pos[frame][block_id].first, frame_to_block_id_to_id_and_pos[frame][block_id].second + 1));
                                    if (frame_to_block_id_to_id_and_pos[frame][block_id].first - end_pos >= kmer_size) {
                                        adjusted_kmer_counts += block_len + kmer_size - 1; //TODO adjust for kmer size overlap
                                    } else {
                                        adjusted_kmer_counts += block_len + frame_to_block_id_to_id_and_pos[frame][block_id].first - end_pos - 1;
                                    }
                                    

                                    end_pos = frame_to_block_id_to_id_and_pos[frame][block_id].second;
                                    if (verbose_flag){
                                        std::cerr << "update end_pos: " << end_pos << std::endl;
                                        std::cerr << std::endl;
                                    }

                                    // look ahead code
                                    if (ref_idx + 2 < vec.size()) {
                                        size_t idx_offset = look_ahead(vec, ref_idx, end_pos, frame_to_block_id_to_id_and_pos, max_offset);
                                        if (idx_offset > 0) {
                                            ref_idx += idx_offset;
                                        }
                                    }
                                }
                                else
                                {
                                    seq_end_in_nucleotide = (prev_block_start + prev_block_len + kmer_size - 1)* 3 + frame;
                                    // log completeness and reset
                                    double prev_score = (double)adjusted_kmer_counts / (double)expected_kmer_counts;
                                    if (adjusted_kmer_counts > lower_bound * expected_kmer_counts)
                                    {
                                        if (start_end_pos_vec.size() > 1) {
                                            fill_in_gaps(start_end_pos_vec, start_end_pos_tar_vec, adjusted_kmer_counts, hash_num, rescue_kmer_size, sixframed_xlated_proteins, ori, kmer_size, miBf_ID, db_path_loc);
                                        }
                                        if (adjusted_kmer_counts > 0.95 * expected_kmer_counts) {
                                            complete_copies++;
                                        } else {
                                            incomplete_copies++;
                                        }
                                    }
                                    double score = (double)adjusted_kmer_counts / (double)expected_kmer_counts;
                                    if (score > 1) {
                                        score = 1;
                                    }
                                    {
                                        if (strand == "-") {
                                            auto temp = seq_start_in_nucleotide;
                                            seq_start_in_nucleotide = record.seq.size() - seq_end_in_nucleotide;
                                            seq_end_in_nucleotide = record.seq.size() - temp;
                                        }
#pragma omp critical
{
                                            pre_gff_set.emplace(record.id, seq_start_in_nucleotide, seq_end_in_nucleotide, prev_score, strand, seq_name);
                                            gff_set.emplace(record.id, seq_start_in_nucleotide, seq_end_in_nucleotide, score, strand, seq_name);
}
                                    }
                                    
                                    end_pos = frame_to_block_id_to_id_and_pos[frame][block_id].second;
                                    adjusted_kmer_counts = block_len;
                                    seq_start_in_nucleotide = FrameBlock_id_and_seq_pos.query_start_in_prot_space * 3 + frame;
                                    start_end_pos_vec.clear();
                                    start_end_pos_tar_vec.clear();
                                    start_end_pos_vec.emplace_back(std::make_tuple(FrameBlock_id_and_seq_pos.query_start_in_prot_space, FrameBlock_id_and_seq_pos.query_start_in_prot_space + block_len - 1));
                                    start_end_pos_tar_vec.emplace_back(std::make_tuple(frame_to_block_id_to_id_and_pos[frame][block_id].first, frame_to_block_id_to_id_and_pos[frame][block_id].second + 1));
                                    if (verbose_flag){
                                        std::cerr << "new end_pos: " << end_pos << std::endl;
                                    }
                                }
                            }
                        }
                    }
                    double prev_score = (double)adjusted_kmer_counts / (double)expected_kmer_counts;
                     if (adjusted_kmer_counts > lower_bound * expected_kmer_counts)
                    {
                        if (start_end_pos_vec.size() > 1) {
                            fill_in_gaps(start_end_pos_vec, start_end_pos_tar_vec, adjusted_kmer_counts, hash_num, rescue_kmer_size, sixframed_xlated_proteins, ori, kmer_size, miBf_ID, db_path_loc);
                        }
                        if (adjusted_kmer_counts > 0.95 * expected_kmer_counts) {
                            complete_copies++;
                        } else {
                            incomplete_copies++;
                        }
                    }
                    // log completeness
#pragma omp atomic
                    seq_name_to_completeness[seq_name].complete_copies += complete_copies;

#pragma omp atomic
                    seq_name_to_completeness[seq_name].incomplete_copies += incomplete_copies;
                    {

                        const auto &last_FrameBlock_id_and_seq_pos = id_to_FrameBlock_id_and_seq_pos[miBf_ID].rbegin();
                        frame = (*last_FrameBlock_id_and_seq_pos).frame;
                        seq_end_in_nucleotide = ((*last_FrameBlock_id_and_seq_pos).query_start_in_prot_space + block_len + kmer_size - 1) * 3 + frame;
                        // print all the values that make seq_end_in_nucleotide
                        double score = (double)adjusted_kmer_counts / (double)expected_kmer_counts;
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
                            pre_gff_set.emplace(record.id, seq_start_in_nucleotide, seq_end_in_nucleotide, prev_score, strand, seq_name);
                            gff_set.emplace(record.id, seq_start_in_nucleotide, seq_end_in_nucleotide, score, strand, seq_name);
                        }
                    }
                }
            }
        }
    }
    // output the completeness to the output file
    auto& seq_name_to_completeness = seq_name_to_completeness_vec[0];
    for (auto &seq_name_completeness : seq_name_to_completeness)
    {
        output_files[0] << seq_name_completeness.first << "\t" << seq_name_completeness.second.complete_copies << "\t" << seq_name_completeness.second.incomplete_copies << "\t" << seq_name_completeness.second.expected_kmer_counts << "\t" << seq_name_completeness.second.highest_adjusted_kmer_counts << std::endl;
    }


    // output the gff set to the gff file

    auto& gff_set = gff_set_vector[0];
    for (auto &gff : gff_set)
    {
        gff_files[0] << gff.query_name << "\t"
                    << "."
                    << "\t"
                    << "gene"
                    << "\t" << gff.hit_pos_start << "\t" << gff.hit_pos_end << "\t" << gff.score << "\t" << gff.strand << "\t"
                    << "0"
                    << "\t"
                    << "ID=" << gff.hit_name << std::endl;
    }

    if (debug_flag) {
        // output the pre gff set to the pre gff file
        auto& pre_gff_set = pre_gff_set_vector[0];
        for (auto &gff : pre_gff_set)
        {
            pre_gff_files[0] << gff.query_name << "\t"
                        << "."
                        << "\t"
                        << "gene"
                        << "\t" << gff.hit_pos_start << "\t" << gff.hit_pos_end << "\t" << gff.score << "\t" << gff.strand << "\t"
                        << "0"
                        << "\t"
                        << "ID=" << gff.hit_name << std::endl;
        }
    }

    return 0;
}
