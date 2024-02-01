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

bool explore_frame(btllib::MIBloomFilter<uint64_t> &mi_bf, btllib::AAHash &aahash, std::deque<std::vector<uint32_t>> &miBf_IDs_snapshot, std::deque<std::vector<uint32_t>> &miBf_pos_snapshot, std::unordered_map<uint32_t, size_t> &id_to_count)
{
    // check size of miBf_IDs_snapshot and miBf_pos_snapshot
    //  if size is more than 10, pop front
    std::unordered_set<uint32_t> id_set;
    if (miBf_IDs_snapshot.size() >= 10)
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
        miBf_IDs_snapshot.back().push_back(ID_pos >> 32);
        miBf_pos_snapshot.back().push_back(ID_pos & 0xFFFFFFFF);
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

    if (miBf_IDs_snapshot.size() < 10)
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
    if (temp_mibf_ID == 0 || temp_max_count < 10)
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

                if (counter >= 9)
                {
                    return true;
                }
            }
        }
        /*std::set<uint32_t> temp_pos_set;
        for (size_t i = 0; i < miBf_IDs_snapshot.size(); ++i)
        {
            for (size_t j = 0; j < miBf_IDs_snapshot[i].size(); ++j)
            {
                if (miBf_IDs_snapshot[i][j] == temp_mibf_ID)
                {
                    temp_pos_set.insert(miBf_pos_snapshot[i][j]);
                }
            }
        }
        // iterate through temp_pos_set and check if it is incrementing by 1
        uint32_t prev_pos = 0;
        size_t counter = 0;
        for (auto &pos : temp_pos_set)
        {
            if (prev_pos == 0)
            {
                prev_pos = pos;
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
        }
        if (counter < 9)
        {
            return false;
        }*/
        /*std::cerr << "check pos before return true" << std::endl;
        for (auto &pos : temp_pos_set)
        {
            std::cerr << "pos: " << pos << std::endl;
        }*/
    }

    return false;
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
    std::ofstream output_file(output_prefix + ".results.tsv");

    if (verbose_flag)
    {
        std::cerr << "Reading reference file: " << reference_path << std::endl;
    }

    // read through reference file which is a fasta file and count the number of characters in the sequences and assign it genome_size
    btllib::SeqReader ref_reader(reference_path, btllib::SeqReader::Flag::LONG_MODE);
    for (const auto record : ref_reader)
    {
        genome_size += record.seq.size();
    }

    btllib::MIBloomFilter<uint64_t> mi_bf(calc_optimal_size(genome_size / 3 * 6, hash_num, 0.1), hash_num);

    if (verbose_flag)
    {
        std::cerr << "Creating seq_id to ID table" << std::endl;
    }

    std::unordered_map<std::string, uint32_t> seq_ID_to_miBf_ID;
    std::unordered_map<uint32_t, std::pair<std::string, size_t>> miBf_ID_to_seq_ID_and_len;
    {
        uint32_t miBf_ID = 1;
        btllib::SeqReader reader(reference_path, btllib::SeqReader::Flag::LONG_MODE);
        for (const auto record : reader)
        {
            // insert record.id into seq_ID_to_miBf_ID with value miBf_ID
            seq_ID_to_miBf_ID[record.id] = miBf_ID;
            // insert miBf_ID into miBf_ID_to_seq_ID_and_len with value record.id and record.seq.size()
            miBf_ID_to_seq_ID_and_len[miBf_ID] = std::make_pair(record.id, record.seq.size());
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
            auto &miBf_ID = seq_ID_to_miBf_ID[record.id];

            while (itr.roll())
            {
                if (stage == 0)
                {
                    mi_bf.insert_bv(itr.hashes());
                }
                else if (stage == 1)
                {
                    uint64_t new_ID = (uint64_t)miBf_ID << 32 | itr.get_pos();
                    mi_bf.insert_id(itr.hashes(), new_ID);
                }
                else
                {
                    uint64_t new_ID = (uint64_t)miBf_ID << 32 | itr.get_pos();
                    mi_bf.insert_saturation(itr.hashes(), new_ID);
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

    /*output_file << "name\thits\trc_hits\t+1_hits\trc_+1_hits\t+2_hits\trc_+2_hits\texpected_hits\tpct_hits" << std::endl;

    btllib::SeqReader reader(input_file, btllib::SeqReader::Flag::LONG_MODE);
    if (verbose_flag) {
        std::cerr << "Reading input file: " << input_file << std::endl;
    }
pragma omp parallel
  for (const auto record : reader) {
    std::vector<std::string> protein = sixframe_translate(record.seq);
    std::vector<std::map<uint32_t, size_t>> frame_to_id_to_hits(6);
    size_t expected_hits = protein[0].size() - kmer_size + 1;
    for (uint8_t i = 0; i < protein.size(); i++) {
        AAHash itr(protein[i], hash_num, kmer_size);
        auto& id_to_hits = frame_to_id_to_hits[i];
        while (itr != AAHash::end()) {
            auto temp_ID_hits =  mi_bf.get_id(*itr); // change this to avoid reallocating memory
            for (auto& ID_hits : temp_ID_hits) {
                if (id_to_hits.find(ID_hits) == id_to_hits.end()) {
                    id_to_hits[ID_hits] = 1;
                } else {
                    id_to_hits[ID_hits]++;
                }
            }
            ++itr;
        }

    }
    std::vector<uint32_t> max_hits(6, 0);
    for (uint8_t i = 0; i < 6; i++) {
        auto& id_to_hits = frame_to_id_to_hits[i];
        max_hits[i] = std::max_element(id_to_hits.begin(), id_to_hits.end(), [](const auto& a, const auto& b) { return a.second < b.second; })->second;
    }

#pragma omp critical
    {
        output_file << record.id << "\t" << max_hits[0] << "\t" << max_hits[1] << "\t" << max_hits[2] << "\t" << max_hits[3] << "\t" << max_hits[4] << "\t" << max_hits[5] << "\t" << expected_hits << "\t" << (double)max_hits[0] / expected_hits << std::endl;
    }
  }
*/

    output_file << "name\tcomplete copies\tincomplete copies\texpected k-mer counts\thighest adjusted incomplete k-mer hits" << std::endl;

    btllib::SeqReader reader(input_file, btllib::SeqReader::Flag::LONG_MODE);
    if (verbose_flag)
    {
        std::cerr << "Reading input file: " << input_file << std::endl;
    }
    // std::unordered_map<std::string, std::tuple<size_t, size_t, size_t, size_t>> seq_name_to_completeness;
    struct completeness_struct
    {
        size_t complete_copies = 0;
        size_t incomplete_copies = 0;
        size_t expected_kmer_counts = 0;
        size_t highest_adjusted_incomplete_kmer_hits = 0;
    };
    std::unordered_map<std::string, completeness_struct> seq_name_to_completeness;
    // populate with empty entries using seq_ID_to_miBf_ID
    for (auto &seq_ID : seq_ID_to_miBf_ID)
    {
        seq_name_to_completeness[seq_ID.first] = completeness_struct();
    }
#pragma omp parallel num_threads(threads / 2)
    for (const auto record : reader)
    {
        // std::cerr << "seq name: " << record.id << std::endl;
        std::vector<std::string> sixframed_xlated_proteins = sixframe_translate(record.seq);
        // std::cerr << "protein 4: " << sixframed_xlated_proteins[3] << std::endl;
        // std::cerr << "protein 5: " << sixframed_xlated_proteins[4] << std::endl;
        // std::cerr << "protein 6: " << sixframed_xlated_proteins[5] << std::endl;
#pragma omp parallel for num_threads(2)
        for (size_t ori = 0; ori < 2; ++ori)
        {
            /*if (verbose_flag)
            {
                std::cerr << "ori:" << ori << std::endl;
            }*/
            std::vector<btllib::AAHash> aahash_itr_vec;
            std::vector<size_t> pos_vec;
            std::vector<bool> hash_itr_state_vec;
            for (size_t i = 0; i < 3; ++i)
            {
                aahash_itr_vec.emplace_back(btllib::AAHash(sixframed_xlated_proteins[i + 3 * ori], hash_num, kmer_size, 1));
                hash_itr_state_vec.push_back(aahash_itr_vec[i].roll());
                pos_vec.push_back(aahash_itr_vec[i].get_pos());
            }
            size_t min_pos = *std::min_element(pos_vec.begin(), pos_vec.end());
            std::vector<bool> use_frame_vec(3, false);
            for (size_t i = 0; i < 3; ++i)
            {
                if (pos_vec[i] == min_pos)
                {
                    use_frame_vec[i] = true;
                }
            }
            bool while_loop_state = false;
            if (hash_itr_state_vec[0] || hash_itr_state_vec[1] || hash_itr_state_vec[2])
            {
                while_loop_state = true;
            }
            bool exploratory_state = true;
            bool elongation_state = false;
            bool searching_state = false;
            size_t curr_frame = 0;
            // exploratory state helpers
            std::vector<std::deque<std::vector<uint32_t>>> miBf_IDs_snapshot_vec(3);
            std::vector<std::deque<std::vector<uint32_t>>> miBf_pos_snapshot_vec(3);
            std::vector<std::unordered_map<uint32_t, size_t>> id_to_count_vec(3);
            std::vector<bool> frame_to_explore(3, false);

            // elongation state helpers
            std::string candidate_protein = "";
            size_t candidate_protein_len = 0;
            uint32_t candidate_mibf_ID = 0;
            size_t frame_switch_count = 0;
            // pos vec for making sure monontously increasing
            // std::vector<uint32_t> candidate_pos_vec;
            std::set<uint32_t> candidate_pos_set;

            // bool insert = false;
            //  std::string insert_protein = "";
            //  size_t insert_protein_len = 0;
            //  size_t insert_frame_switch_count = 0;
            //  size_t insert_pos_set_size = 0;
            // std::cerr << "pre while loop" << std::endl;
            std::unordered_map<std::string, std::tuple<size_t, size_t, size_t, size_t>> insert_saved_state;
            while (while_loop_state)
            {
                //   advance all iterators by 1
                //   need do for all three frames
                if (exploratory_state)
                {

                    for (size_t a = 0; a < 3; ++a)
                    {
                        if (use_frame_vec[a])
                        {
                            frame_to_explore[a] = explore_frame(mi_bf, aahash_itr_vec[a], miBf_IDs_snapshot_vec[a], miBf_pos_snapshot_vec[a], id_to_count_vec[a]);
                        }
                    }
                    if (frame_to_explore[0] || frame_to_explore[1] || frame_to_explore[2])
                    {
                        // std::cerr << "start elongation" << std::endl;
                        exploratory_state = false;
                        elongation_state = true;

                        // find the frame with to elongate
                        for (size_t a = 0; a < 3; ++a)
                        {
                            if (frame_to_explore[a])
                            {
                                curr_frame = a;
                                break;
                            }
                        }
                        // std::cerr << "curr_frame:" << curr_frame << std::endl;
                        //   find the id with the highest count using id_to_count_vec
                        uint32_t temp_mibf_ID = 0;
                        size_t temp_max_count = 0;
                        for (auto &ID_count : id_to_count_vec[curr_frame])
                        {
                            if (ID_count.second > temp_max_count)
                            {
                                temp_mibf_ID = ID_count.first;
                                temp_max_count = ID_count.second;
                            }
                        }
                        /*if (verbose_flag)
                        {
                            std::cerr << "temp_mibf_ID:" << temp_mibf_ID << std::endl;
                            std::cerr << "temp_max_count:" << temp_max_count << std::endl;
                            std::cerr << "curr_frame: " << curr_frame << std::endl;
                        }*/

                        // get candidate protein from ID
                        candidate_protein = miBf_ID_to_seq_ID_and_len[temp_mibf_ID].first;
                        candidate_protein_len = miBf_ID_to_seq_ID_and_len[temp_mibf_ID].second;
                        candidate_mibf_ID = temp_mibf_ID;
                        // insert the smallest pos larger or equal to the size of the set into candidate_pos_set
                        for (size_t i = 0; i < miBf_IDs_snapshot_vec[curr_frame].size(); ++i)
                        {
                            std::set<uint32_t> temp_pos_set;
                            for (size_t j = 0; j < miBf_IDs_snapshot_vec[curr_frame][i].size(); ++j)
                            {
                                if (miBf_IDs_snapshot_vec[curr_frame][i][j] == candidate_mibf_ID)
                                {
                                    temp_pos_set.insert(miBf_pos_snapshot_vec[curr_frame][i][j]);
                                }
                            }
                            for (auto &pos : temp_pos_set)
                            {
                                if (pos >= candidate_pos_set.size())
                                {
                                    candidate_pos_set.insert(pos);
                                    // std::cerr << "pos: " << pos << std::endl;
                                    break;
                                }
                            }
                        }
                        // clear exploraty state helpers
                        for (size_t a = 0; a < 3; ++a)
                        {
                            miBf_IDs_snapshot_vec[a].clear();
                            miBf_pos_snapshot_vec[a].clear();
                            id_to_count_vec[a].clear();
                            frame_to_explore[a] = false;
                        }
                    }
                }
                else if (elongation_state)
                {
                    // std::cerr << "on elongation" << std::endl;
                    std::vector<uint32_t> ids_vec;
                    std::vector<uint32_t> temp_pos_vec;
                    auto &aaHash = aahash_itr_vec[curr_frame];
                    size_t prev_pos_set_size = candidate_pos_set.size();
                    // std::cerr << "prev_pos_set_size: " << prev_pos_set_size << std::endl;
                    // std::cerr << "1" << std::endl;
                    if (mi_bf.bv_contains(aaHash.hashes()))
                    {
                        auto temp_ID_pos = mi_bf.get_id(aaHash.hashes());
                        // std::cerr << "2" << std::endl;

                        for (auto &ID_pos : temp_ID_pos)
                        {
                            ids_vec.push_back(ID_pos >> 32);
                            temp_pos_vec.push_back(ID_pos & 0xFFFFFFFF);
                        }
                        // std::cerr << "3" << std::endl;
                        std::set<uint32_t> temp_pos_set;
                        for (size_t i = 0; i < ids_vec.size(); ++i)
                        {

                            if (ids_vec[i] == candidate_mibf_ID)
                            {
                                temp_pos_set.insert(temp_pos_vec[i]);
                            }
                        }
                        // std::cerr << "4" << std::endl;
                        for (auto &pos : temp_pos_set)
                        {
                            if (pos == *candidate_pos_set.rbegin() + 1)
                            {
                                // std::cerr << "pos: " << pos << std::endl;
                                candidate_pos_set.insert(pos);
                                break;
                            }
                        }
                        // std::cerr << "5" << std::endl;
                    }
                    // std::cerr << "6" << std::endl;
                    if (candidate_pos_set.size() == prev_pos_set_size)
                    {
                        // std::cerr << "entering searching state" << std::endl;
                        searching_state = true;
                        elongation_state = false;
                    }
                    // std::cerr << "done elongation" << std::endl;
                }
                else if (searching_state)
                {

                    // std::cerr << "on searching" << std::endl;
                    std::vector<bool> frame_to_explore(3, false);
                    for (size_t a = 0; a < 3; ++a)
                    {
                        if (use_frame_vec[a])
                        {
                            frame_to_explore[a] = explore_frame(mi_bf, aahash_itr_vec[a], miBf_IDs_snapshot_vec[a], miBf_pos_snapshot_vec[a], id_to_count_vec[a]);
                        }
                    }
                    if (frame_to_explore[0] || frame_to_explore[1] || frame_to_explore[2])
                    {
                        // find the frame with to elongate
                        for (size_t a = 0; a < 3; ++a)
                        {
                            if (frame_to_explore[a])
                            {
                                curr_frame = a;
                                break;
                            }
                        }
                        // find the id with the highest count using id_to_count_vec
                        uint32_t temp_mibf_ID = 0;
                        size_t temp_max_count = 0;
                        bool saved_state_candidate = false;
                        size_t pos_delta = std::numeric_limits<size_t>::max();
                        // std::cerr << "choosing best ID" << std::endl;
                        size_t largest_count = 0;
                        // search id_to_count_vec for largest count
                        for (auto &ID_count : id_to_count_vec[curr_frame])
                        {
                            if (ID_count.second > largest_count)
                            {
                                largest_count = ID_count.second;
                            }
                        }
                        for (auto &ID_count : id_to_count_vec[curr_frame])
                        {
                            // check if ID_count.second is equal to largest_count, if no, continue to next loop
                            if (ID_count.second != largest_count)
                            {
                                continue;
                            }
                            // check if ID_count.first is in miBf_ID_to_seq_ID_and_len
                            // if no, continue to next loop
                            if (miBf_ID_to_seq_ID_and_len.find(ID_count.first) == miBf_ID_to_seq_ID_and_len.end())
                            {
                                continue;
                            }
                            std::string saved_state_candidate_protein = miBf_ID_to_seq_ID_and_len[ID_count.first].first;
                            // std::cerr << "saved_state_candidate_protein: " << saved_state_candidate_protein << std::endl;
                            //  check if ID_count.first is in insert_saved_state, if yes save result into saved_state_candidate and end_pos else cotinue
                            if (insert_saved_state.find(saved_state_candidate_protein) != insert_saved_state.end())
                            {
                                // std::cerr << "found " << saved_state_candidate_protein << " in insert_saved_state" << std::endl;
                                auto end_pos_on_candidate = std::get<3>(insert_saved_state[saved_state_candidate_protein]);
                                // std::cerr << "end_pos_on_candidate: " << end_pos_on_candidate << std::endl;
                                std::set<uint32_t> temp_pos_set;
                                for (size_t i = 0; i < miBf_IDs_snapshot_vec[curr_frame].size(); ++i)
                                {
                                    for (size_t j = 0; j < miBf_IDs_snapshot_vec[curr_frame][i].size(); ++j)
                                    {
                                        if (miBf_IDs_snapshot_vec[curr_frame][i][j] == ID_count.first)
                                        {
                                            temp_pos_set.insert(miBf_pos_snapshot_vec[curr_frame][i][j]);
                                        }
                                    }
                                }

                                /*std::cerr << "temp_pos_set.begin(): " << *temp_pos_set.begin() << std::endl;
                                std::cerr << "pos_delta: " << pos_delta << std::endl;*/
                                if (pos_delta > *temp_pos_set.begin() - end_pos_on_candidate)
                                {
                                    pos_delta = *temp_pos_set.begin() - end_pos_on_candidate;
                                    temp_mibf_ID = ID_count.first;
                                    temp_max_count = ID_count.second;
                                    saved_state_candidate = true;
                                }
                            }
                        }
                        if (!saved_state_candidate)
                        {
                            for (auto &ID_count : id_to_count_vec[curr_frame])
                            {
                                if (ID_count.second > temp_max_count)
                                {
                                    temp_mibf_ID = ID_count.first;
                                    temp_max_count = ID_count.second;
                                }
                                else if (ID_count.second >= temp_max_count && ID_count.first == candidate_mibf_ID)
                                {
                                    temp_mibf_ID = ID_count.first;
                                    temp_max_count = ID_count.second;
                                }
                            }
                        }

                        // std::cerr << "done choosing best ID" << std::endl;

                        elongation_state = true;
                        searching_state = false;
                        /*std::cerr << "temp_mibf_ID: " << temp_mibf_ID << std::endl;
                        std::cerr << "candidate_mibf_ID: " << candidate_mibf_ID << std::endl;
                        // std::cerr << "curr_frame: " << curr_frame << std::endl;
                        std::cerr << "candidate_pos_set size: " << candidate_pos_set.size() << std::endl;*/
                        // std::cerr << "temp_pos_set size: " << temp_pos_set.size() << std::endl;
                        if (temp_mibf_ID == candidate_mibf_ID)
                        {
                            // found candidate protein
                            // check size of candidate_pos_set matches candidate_protein_len
                            ++frame_switch_count;
                            for (size_t i = 0; i < miBf_IDs_snapshot_vec[curr_frame].size(); ++i)
                            {
                                std::set<uint32_t> temp_pos_set;
                                for (size_t j = 0; j < miBf_IDs_snapshot_vec[curr_frame][i].size(); ++j)
                                {
                                    if (miBf_IDs_snapshot_vec[curr_frame][i][j] == candidate_mibf_ID)
                                    {
                                        temp_pos_set.insert(miBf_pos_snapshot_vec[curr_frame][i][j]);
                                    }
                                }
                                for (auto &pos : temp_pos_set)
                                {
                                    if (pos >= candidate_pos_set.size())
                                    {
                                        // std::cerr << "pos: " << pos << std::endl;
                                        candidate_pos_set.insert(pos);
                                        break;
                                    }
                                }
                            }
                            // std::cerr << "candidate_pos_set size: " << candidate_pos_set.size() << std::endl;
                        }
                        else
                        {
                            // print out candidate id and mibf id
                            /*std::cerr << "switching" << std::endl;
                            std::cerr << "candidate_mibf_ID: " << candidate_mibf_ID << std::endl;
                            std::cerr << "candidate_protein: " << candidate_protein << std::endl;
                            std::cerr << "candidate_protein_len: " << candidate_protein_len << std::endl;
                            std::cerr << "temp_mibf_ID: " << temp_mibf_ID << std::endl;*/
                            // insert = true;
                            /*insert_protein = candidate_protein;
                            insert_protein_len = candidate_protein_len;
                            insert_frame_switch_count = frame_switch_count;
                            insert_pos_set_size = candidate_pos_set.size();*/
                            if (insert_saved_state.find(candidate_protein) != insert_saved_state.end())
                            {
                                // std::cerr << "found " << candidate_protein << " in insert_saved_state" << std::endl;
                                if (std::get<3>(insert_saved_state[candidate_protein]) < *candidate_pos_set.begin())
                                {
                                    size_t prev_frame_switch_count = std::get<1>(insert_saved_state[candidate_protein]);
                                    size_t prev_pos_set_size = std::get<2>(insert_saved_state[candidate_protein]);
                                    insert_saved_state[candidate_protein] = std::make_tuple(candidate_protein_len, prev_frame_switch_count + frame_switch_count + 1, prev_pos_set_size + candidate_pos_set.size(), *candidate_pos_set.rbegin());
                                    /*std::cerr << "updating insert_saved_state" << std::endl;
                                    std::cerr << "candidate_protein_len: " << candidate_protein_len << std::endl;
                                    std::cerr << "prev_frame_switch_count: " << prev_frame_switch_count << std::endl;
                                    std::cerr << "frame_switch_count: " << frame_switch_count << std::endl;
                                    std::cerr << "prev_pos_set_size: " << prev_pos_set_size << std::endl;
                                    std::cerr << "candidate_pos_set.size(): " << candidate_pos_set.size() << std::endl;
                                    std::cerr << "*candidate_pos_set.rbegin(): " << *candidate_pos_set.rbegin() << std::endl;*/
                                }
                                else
                                {
                                    size_t prev_protein_len = std::get<0>(insert_saved_state[candidate_protein]);
                                    size_t prev_frame_switch_count = std::get<1>(insert_saved_state[candidate_protein]);
                                    size_t prev_pos_set_size = std::get<2>(insert_saved_state[candidate_protein]);
                                    size_t adjusted_kmer_hit = prev_pos_set_size + prev_frame_switch_count * (kmer_size - 1);
                                    /*std::cerr << "adjusted_kmer_hit: " << adjusted_kmer_hit << std::endl;
                                    std::cerr << "max_kmer_count: " << prev_protein_len - kmer_size + 1 << std::endl;
                                    std::cerr << "name of protein: " << candidate_protein << std::endl;*/
                                    if (adjusted_kmer_hit >= 0.95 * (prev_protein_len - kmer_size + 1))
                                    {
                                        // update seq_name_to_completeness

                                        seq_name_to_completeness[candidate_protein].complete_copies++;
                                    }
                                    else if (adjusted_kmer_hit >= 0.5 * (prev_protein_len - kmer_size + 1))
                                    {
                                        seq_name_to_completeness[candidate_protein].incomplete_copies++;
                                    }
                                    insert_saved_state[candidate_protein] = std::make_tuple(candidate_protein_len, 0, 0, 0);
                                }
                            }
                            else
                            {
                                insert_saved_state[candidate_protein] = std::make_tuple(candidate_protein_len, frame_switch_count, candidate_pos_set.size(), *candidate_pos_set.rbegin());
                                /*std::cerr << "inserting into insert_saved_state" << std::endl;
                                std::cerr << "candidate_protein_len: " << candidate_protein_len << std::endl;
                                std::cerr << "frame_switch_count: " << frame_switch_count << std::endl;
                                std::cerr << "candidate_pos_set.size(): " << candidate_pos_set.size() << std::endl;
                                std::cerr << "*candidate_pos_set.rbegin(): " << *candidate_pos_set.rbegin() << std::endl;*/
                            }

                            // get candidate protein from ID
                            candidate_pos_set.clear();
                            candidate_protein = miBf_ID_to_seq_ID_and_len[temp_mibf_ID].first;
                            candidate_protein_len = miBf_ID_to_seq_ID_and_len[temp_mibf_ID].second;
                            candidate_mibf_ID = temp_mibf_ID;
                            frame_switch_count = 0;
                            // insert the smallest pos larger or equal to the size of the set into candidate_pos_set
                            for (size_t i = 0; i < miBf_IDs_snapshot_vec[curr_frame].size(); ++i)
                            {
                                std::set<uint32_t> temp_pos_set;
                                for (size_t j = 0; j < miBf_IDs_snapshot_vec[curr_frame][i].size(); ++j)
                                {
                                    if (miBf_IDs_snapshot_vec[curr_frame][i][j] == candidate_mibf_ID)
                                    {
                                        temp_pos_set.insert(miBf_pos_snapshot_vec[curr_frame][i][j]);
                                    }
                                }
                                for (auto &pos : temp_pos_set)
                                {
                                    if (pos >= candidate_pos_set.size())
                                    {
                                        // std::cerr << "pos: " << pos << std::endl;
                                        candidate_pos_set.insert(pos);
                                        break;
                                    }
                                }
                            }
                            // print out content of candidate_pos_set and size
                            // std::cerr << "new candidate_pos_set size: " << candidate_pos_set.size() << std::endl;
                            /*for (auto &pos : candidate_pos_set)
                            {
                                std::cerr << "pos: " << pos << std::endl;
                            }*/
                            // clear exploraty state helpers
                            for (size_t a = 0; a < 3; ++a)
                            {
                                miBf_IDs_snapshot_vec[a].clear();
                                miBf_pos_snapshot_vec[a].clear();
                                id_to_count_vec[a].clear();
                                frame_to_explore[a] = false;
                            }
                        }
                    }
                    // std::cerr << "done searching" << std::endl;
                }
                // std::cerr << "done processing" << std::endl;
                // std::cerr << "curr_frame: " << curr_frame << std::endl;
                bool break_loop = false;
                for (size_t i = 0; i < 3; ++i)
                {
                    if (use_frame_vec[i])
                    {
                        // print seq and seq len
                        size_t prev_hash_pos = aahash_itr_vec[i].get_pos();
                        hash_itr_state_vec[i] = aahash_itr_vec[i].roll();
                        size_t curr_hash_pos = aahash_itr_vec[i].get_pos();
                        if (prev_hash_pos != std::numeric_limits<std::size_t>::max() && prev_hash_pos == curr_hash_pos)
                        {
                            // aahash_itr_vec[i] = btllib::AAHash(sixframed_xlated_proteins[i + 3 * ori].substr(prev_hash_pos + 1), hash_num, kmer_size, 1);
                            break_loop = true;
                        }
                        pos_vec[i] = aahash_itr_vec[i].get_pos();
                    }
                }

                // set min_pos to min value of pos_vec
                /*std::cerr << "size of pos_vec: " << pos_vec.size() << std::endl;
                for (auto &pos : pos_vec)
                {
                    std::cerr << "pos: " << pos << std::endl;
                }*/
                min_pos = *std::min_element(pos_vec.begin(), pos_vec.end());

                for (size_t i = 0; i < 3; ++i)
                {
                    if (pos_vec[i] == min_pos)
                    {
                        use_frame_vec[i] = true;
                    }
                    else
                    {
                        use_frame_vec[i] = false;
                    }
                }
                if (elongation_state)
                {
                    use_frame_vec[curr_frame] = true;
                }

                /*std::string true_or_false = "";
                for (size_t i = 0; i < 3; ++i)
                {
                    if (use_frame_vec[i])
                    {
                        true_or_false = "true";
                    }
                    else
                    {
                        true_or_false = "false";
                    }
                    std::cerr << "frame: " << i << " pos: " << pos_vec[i] << " use_frame: " << true_or_false << std::endl;
                    // print state of hash_itr_state_vec
                    if (hash_itr_state_vec[i])
                    {
                        true_or_false = "true";
                    }
                    else
                    {
                        true_or_false = "false";
                    }
                    std::cerr << "frame: " << i << " hash_itr_state: " << true_or_false << std::endl;
                }*/
                // std::cerr << "done processing" << std::endl;
                while_loop_state = hash_itr_state_vec[0] || hash_itr_state_vec[1] || hash_itr_state_vec[2];
                // create hash table ahead of time with entries

                /*if (insert)
                {
                    // std::cerr << "test1" << std::endl;
                    size_t adjusted_kmer_hit = insert_pos_set_size + insert_frame_switch_count * (kmer_size - 1);
                    std::cerr << "adjusted_kmer_hit: " << adjusted_kmer_hit << std::endl;
                    std::cerr << "max_kmer_count: " << insert_protein_len - kmer_size + 1 << std::endl;
                    std::cerr << "name of protein: " << insert_protein << std::endl;
                    if (adjusted_kmer_hit >= 0.95 * (insert_protein_len - kmer_size + 1))
                    {
                        // update seq_name_to_completeness

                        seq_name_to_completeness[insert_protein].complete_copies++;
                    }
                    else if (adjusted_kmer_hit >= 0.5 * (insert_protein_len - kmer_size + 1))
                    {
                        // update seq_name_to_completeness
                        seq_name_to_completeness[insert_protein].incomplete_copies++;
                        // std::get<1>(seq_name_to_completeness[insert_protein])++;
                        // std::get<3>(seq_name_to_completeness[insert_protein]) += adjusted_kmer_hit;
                    }
                    insert = false;
                }*/

                if (break_loop)
                {
                    break;
                }

                // refactor critical to out of while loop
            }
            /*std::cerr << "candidate_protein: " << candidate_protein << std::endl;
            std::cerr << "*candidate_pos_set.begin(): " << *candidate_pos_set.begin() << std::endl;
            std::cerr << "*candidate_pos_set.begin(): " << *candidate_pos_set.begin() << std::endl;
            std::cerr << "candidate_pos_set.size(): " << candidate_pos_set.size() << std::endl;
            // print out content of candidate_pos_set
            for (auto &pos : candidate_pos_set)
            {
                std::cerr << "pos: " << pos << std::endl;
            }
            std::cerr << "std::get<3>(insert_saved_state[candidate_protein]): " << std::get<3>(insert_saved_state[candidate_protein]) << std::endl;*/
            if (candidate_protein != "")
            {
                // std::cerr << "test2" << std::endl;
                // std::cerr << "candidate protein: " << candidate_protein << std::endl;
                //  update seq_name_to_completeness
                /*size_t adjusted_kmer_hit = candidate_pos_set.size() + frame_switch_count * (kmer_size - 1);
                std::cerr << "adjusted_kmer_hit: " << adjusted_kmer_hit << std::endl;
                std::cerr << "max_kmer_count: " << candidate_protein_len - kmer_size + 1 << std::endl;
                std::cerr << "name of protein: " << candidate_protein << std::endl;
                if (adjusted_kmer_hit >= 0.95 * (candidate_protein_len - kmer_size + 1))
                {
                    // update seq_name_to_completeness

                    seq_name_to_completeness[candidate_protein].complete_copies++;
                }
                else if (adjusted_kmer_hit >= 0.5 * (candidate_protein_len - kmer_size + 1))
                {
                    seq_name_to_completeness[candidate_protein].incomplete_copies++;
                }*/
                if (insert_saved_state.find(candidate_protein) != insert_saved_state.end())
                {
                    // std::cerr << "found " << candidate_protein << " in insert_saved_state" << std::endl;
                    if (std::get<3>(insert_saved_state[candidate_protein]) < *candidate_pos_set.begin())
                    {
                        size_t prev_frame_switch_count = std::get<1>(insert_saved_state[candidate_protein]);
                        size_t prev_pos_set_size = std::get<2>(insert_saved_state[candidate_protein]);
                        insert_saved_state[candidate_protein] = std::make_tuple(candidate_protein_len, prev_frame_switch_count + frame_switch_count + 1, prev_pos_set_size + candidate_pos_set.size(), *candidate_pos_set.rbegin());
                        /*std::cerr << "updating insert_saved_state" << std::endl;
                        std::cerr << "candidate_protein_len: " << candidate_protein_len << std::endl;
                        std::cerr << "prev_frame_switch_count: " << prev_frame_switch_count << std::endl;
                        std::cerr << "frame_switch_count: " << frame_switch_count << std::endl;
                        std::cerr << "prev_pos_set_size: " << prev_pos_set_size << std::endl;
                        std::cerr << "candidate_pos_set.size(): " << candidate_pos_set.size() << std::endl;
                        std::cerr << "*candidate_pos_set.rbegin(): " << *candidate_pos_set.rbegin() << std::endl;*/
                    }
                    else
                    {
                        size_t prev_protein_len = std::get<0>(insert_saved_state[candidate_protein]);
                        size_t prev_frame_switch_count = std::get<1>(insert_saved_state[candidate_protein]);
                        size_t prev_pos_set_size = std::get<2>(insert_saved_state[candidate_protein]);
                        size_t adjusted_kmer_hit = prev_pos_set_size + prev_frame_switch_count * (kmer_size - 1);
                        /*std::cerr << "adjusted_kmer_hit: " << adjusted_kmer_hit << std::endl;
                        std::cerr << "max_kmer_count: " << prev_protein_len - kmer_size + 1 << std::endl;
                        std::cerr << "name of protein: " << candidate_protein << std::endl;*/
                        if (adjusted_kmer_hit >= 0.95 * (prev_protein_len - kmer_size + 1))
                        {
                            // update seq_name_to_completeness

                            seq_name_to_completeness[candidate_protein].complete_copies++;
                        }
                        else if (adjusted_kmer_hit >= 0.5 * (prev_protein_len - kmer_size + 1))
                        {
                            seq_name_to_completeness[candidate_protein].incomplete_copies++;
                        }
                        insert_saved_state[candidate_protein] = std::make_tuple(candidate_protein_len, 0, 0, 0);
                    }
                }
                else
                {
                    insert_saved_state[candidate_protein] = std::make_tuple(candidate_protein_len, frame_switch_count, candidate_pos_set.size(), *candidate_pos_set.rbegin());
                    /*std::cerr << "candidate_protein" << candidate_protein << std::endl;
                    std::cerr << "inserting into insert_saved_state" << std::endl;
                    std::cerr << "candidate_protein_len: " << candidate_protein_len << std::endl;
                    std::cerr << "frame_switch_count: " << frame_switch_count << std::endl;
                    std::cerr << "candidate_pos_set.size(): " << candidate_pos_set.size() << std::endl;
                    std::cerr << "*candidate_pos_set.rbegin(): " << *candidate_pos_set.rbegin() << std::endl;*/
                }
            }

            for (const auto &saved_state : insert_saved_state)
            {
                // std::cerr << "final step" << std::endl;
                size_t adjusted_kmer_hit = std::get<2>(saved_state.second) + std::get<1>(saved_state.second) * (kmer_size - 1);
                /*std::cerr << "adjusted_kmer_hit: " << adjusted_kmer_hit << std::endl;
                std::cerr << "max_kmer_count: " << std::get<0>(saved_state.second) - kmer_size + 1 << std::endl;
                std::cerr << "name of protein: " << saved_state.first << std::endl;*/
                if (adjusted_kmer_hit >= 0.95 * (std::get<0>(saved_state.second) - kmer_size + 1))
                {
                    // update seq_name_to_completeness

                    seq_name_to_completeness[saved_state.first].complete_copies++;
                }
                else if (adjusted_kmer_hit >= 0.5 * (std::get<0>(saved_state.second) - kmer_size + 1))
                {
                    seq_name_to_completeness[saved_state.first].incomplete_copies++;
                }
            }
        }
    }
    std::cerr << "done processing" << std::endl;

    for (auto &seq_name_completeness : seq_name_to_completeness)
    {
        // output_file << seq_name_completeness.first << "\t" << std::get<0>(seq_name_completeness.second) << "\t" << std::get<1>(seq_name_completeness.second) << "\t" << std::get<2>(seq_name_completeness.second) << "\t" << std::get<3>(seq_name_completeness.second) << std::endl;
        output_file << seq_name_completeness.first << "\t" << seq_name_completeness.second.complete_copies << "\t" << seq_name_completeness.second.incomplete_copies << "\t" << seq_name_completeness.second.expected_kmer_counts << "\t" << seq_name_completeness.second.highest_adjusted_incomplete_kmer_hits << std::endl;
    }

    return 0;
}
