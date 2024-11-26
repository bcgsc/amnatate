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

struct frame_block {
    size_t frame;
    size_t block_id;
    size_t query_start_in_prot_space;

    // Constructor for convenience
    frame_block(size_t f, size_t b, size_t q) : frame(f), block_id(b), query_start_in_prot_space(q) {}
};

// Define the comparator for FrameBlock
struct frame_block_comparator {
    bool operator()(const frame_block& a, const frame_block& b) const {
        return a.query_start_in_prot_space < b.query_start_in_prot_space;
    }
};

size_t look_ahead(const std::vector<std::reference_wrapper<const frame_block>> &vec, size_t ref_idx, size_t end_pos, 
                 const std::unordered_map<size_t, std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>>> &frame_to_block_id_to_id_and_pos, size_t offset) {

    if (ref_idx + 1 >= vec.size()) {
        return 0;
    }

    // Special logic for the first position check
    const frame_block &next_frame_block = vec[ref_idx + 1].get();
    auto next_block_id = next_frame_block.block_id;
    auto next_frame = next_frame_block.frame;

    if (end_pos < frame_to_block_id_to_id_and_pos.at(next_frame).at(next_block_id).first) {
        return 0;
    }

    // General logic for subsequent positions
    for (size_t i = 2; i <= offset; ++i) {
        if (ref_idx + i >= vec.size()) {
            return 0;
        }

        const frame_block &frame_block_i = vec[ref_idx + i].get();
        auto block_id = frame_block_i.block_id;
        auto frame = frame_block_i.frame;

        if (end_pos < frame_to_block_id_to_id_and_pos.at(frame).at(block_id).first) {
            return i - 1;
        }
    }

    return offset;
}

void process_hashes(const std::vector<uint64_t>& temp_ID_pos, std::unordered_set<uint32_t>& id_set, 
                    std::unordered_map<uint32_t, std::set<uint32_t>>& id_to_pos_set, bool& extend_block, 
                    std::vector<uint32_t>& ids_vec, std::vector<uint32_t>& temp_pos_vec, const btllib::MIBloomFilter<uint64_t>& mi_bf)
{
    bool found = false;

    for (auto &ID : id_set)
    {
        for (auto &ID_pos : temp_ID_pos)
        {
            auto demasked_ID_pos = ID_pos & mi_bf.ANTI_MASK;
            if (ID == (demasked_ID_pos >> 32))
            {
                found = true;
                break;
            }
        }
        if (found)
        {
            break;
        }
    }

    if (found)
    {
        for (auto &ID_pos : temp_ID_pos)
        {
            auto demasked_ID_pos = ID_pos & mi_bf.ANTI_MASK;
            ids_vec.push_back(demasked_ID_pos >> 32);
            temp_pos_vec.push_back(demasked_ID_pos & 0xFFFFFFFF);
        }
        std::vector<uint32_t> new_pos;
        std::vector<uint32_t> new_ids;

        for (size_t i = 0; i < ids_vec.size(); ++i)
        {
            // check if ID is in id_set
            std::set<uint32_t> temp_pos_set;
            
            if (id_set.find(ids_vec[i]) != id_set.end())
            {
                temp_pos_set.insert(temp_pos_vec[i]);
            } else {                
                new_ids.push_back(ids_vec[i]);
                new_pos.push_back(temp_pos_vec[i]);
            }
            for (auto &pos : temp_pos_set)
            {
                if (id_to_pos_set[ids_vec[i]].size() != 0 && pos == *id_to_pos_set[ids_vec[i]].rbegin() + 1)
                {
                    id_to_pos_set[ids_vec[i]].insert(pos);
                    break;
                }
            }

        }

        // make a set from new_ids
        // insert the smallest pos of the new_pos associated with id from new_ids into id_to_pos_set
        std::unordered_set<uint32_t> new_id_set(new_ids.begin(), new_ids.end()); // TODO IMPROVE
        for (auto &ID : new_id_set)
        {
            std::set<uint32_t> temp_pos_set;
            for (size_t i = 0; i < new_ids.size(); ++i)
            {
                if (new_ids[i] == ID)
                {
                    temp_pos_set.insert(new_pos[i]);
                }
            }
            if (temp_pos_set.size() != 0)
            {
                id_to_pos_set[ID].insert(*temp_pos_set.begin());
            }
        }


        // make a set from id
        // insert the smallest pos of the new_temp_pos_set into id_to_pos_set


        // check if the pos set in id_to_pos_set are all the same size
        /*
        size_t temp_size = 0;
        for (auto &ID_pos_set : id_to_pos_set)
        {
            if (temp_size == 0)
            {
                temp_size = ID_pos_set.second.size();
            }
            else if (temp_size != ID_pos_set.second.size())
            {
                extend_block = false;
                break;
            }
        }
        */
    }
    else
    {
        bool saturated = true;
        for (auto &ID_pos : temp_ID_pos)
        {
            if (ID_pos < mi_bf.MASK)
            {
                saturated = false;
                break;
            }
        }
        if (!saturated)
        {
            extend_block = false;
        }
        else
        {
            extend_block = true;
            for (auto& pos_set : id_to_pos_set) {
                if (pos_set.second.size() == 0) {
                    continue;
                } else {
                    pos_set.second.insert(*pos_set.second.rbegin() + 1);
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

btllib::MIBloomFilter<uint64_t> make_mibf(const std::string& seq, const size_t hash_num, const size_t kmer_size) 
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

void fill_in_gaps(std::vector<std::tuple<size_t, size_t>>& start_end_pos_vec, std::vector<std::tuple<size_t, size_t>>& start_end_pos_in_tar_space_vec, size_t& adjusted_kmer_counts, const std::string& seq, size_t hash_num, size_t rescue_kmer_size, const std::vector<std::string>& sixframed_xlated_proteins, size_t ori, size_t kmer_size, size_t level, size_t threads)
{

    //sort start_end_pos_vec by start_pos
    std::sort(start_end_pos_vec.begin(), start_end_pos_vec.end(), [](const std::tuple<size_t, size_t> &a, const std::tuple<size_t, size_t> &b)
    {
        return std::get<0>(a) < std::get<0>(b);
    });
    // std::cerr << "checkpoint a" << std::endl;

    // make a new vector and find all the gaps between the start and end pos that are larger than kmer_size
    std::vector<std::tuple<size_t, size_t>> gap_vec;
    //output size of start_end_pos_vec
    //std::cerr << "size of start_end_pos_vec: " << start_end_pos_vec.size() << std::endl;
    /*for (const auto& tuple : start_end_pos_vec) {
        std::cout << "(" << std::get<0>(tuple) << ", " << std::get<1>(tuple) << ") ";
    }
    std::cout << std::endl;*/
    for (size_t i = 0; i < start_end_pos_vec.size() - 1; ++i)
    {
        if (std::get<0>(start_end_pos_vec[i + 1]) - std::get<1>(start_end_pos_vec[i]) > kmer_size)
        {
            gap_vec.emplace_back(std::make_tuple(std::get<1>(start_end_pos_vec[i]), std::get<0>(start_end_pos_vec[i + 1])));
        }
    }
    if (gap_vec.size() == 0) {
        return;
    }

    // std::cerr << "checkpoint b" << std::endl;
    /*for (const auto& tuple : gap_vec) {
        std::cout << "(" << std::get<0>(tuple) << ", " << std::get<1>(tuple) << ") ";
    }
    std::cout << std::endl;*/

    // make a new vector and find all the gaps between the start and end pos in tar that are larger than kmer_size
    std::vector<std::tuple<size_t, size_t>> gap_in_tar_space_vec;
    //output size of start_end_pos_in_tar_space_vec
    //std::cerr << "size of start_end_pos_in_tar_space_vec: " << start_end_pos_in_tar_space_vec.size() << std::endl;
    /*for (const auto& tuple : start_end_pos_in_tar_space_vec) {
        std::cout << "(" << std::get<0>(tuple) << ", " << std::get<1>(tuple) << ") ";
    }*/
    for (size_t i = 0; i < start_end_pos_in_tar_space_vec.size() - 1; ++i)
    {
        if (std::get<0>(start_end_pos_in_tar_space_vec[i + 1]) - (std::get<1>(start_end_pos_in_tar_space_vec[i]) + 5 ) > kmer_size) //TODO
        {
            gap_in_tar_space_vec.emplace_back(std::make_tuple(std::get<1>(start_end_pos_in_tar_space_vec[i]) + 5, std::get<0>(start_end_pos_in_tar_space_vec[i + 1])));
        }
    }

    // std::cerr << "checkpoint c" << std::endl;

    // make an unordered_set of all the pos in the gaps in start_end_pos_in_tar_space_vec
    std::unordered_set<size_t> gap_index_set;
    for (auto &gap : gap_in_tar_space_vec)
    {
        //std::cerr << "gaps: " << std::endl;
        for (size_t i = std::get<0>(gap); i < std::get<1>(gap); ++i)
        {
            gap_index_set.insert(i);
            //std::cerr << i << std::endl;
        }
    }

    // std::cerr << "checkpoint d" << std::endl;

    //std::cerr << "done gap index" << std::endl;


    // make small mibf
    //std::cerr << "seq to be inserted into small mibf: " << std::endl;
    //std::cerr << seq << std::endl;
    if (gap_index_set.size() == 0)
    {
        return;
    }
    // std::cerr << "checkpoint d.1" << std::endl;
    // log all of the input to the mibf
    //std::cerr << "seq size: " << seq.size() << std::endl;
    //std::cerr << "hash_num: " << hash_num << std::endl;
    //std::cerr << "rescue_kmer_size: " << rescue_kmer_size << std::endl;
  //omp_set_num_threads(1);
    auto small_mi_bf = make_mibf(seq, hash_num, rescue_kmer_size);
  omp_set_num_threads(threads);
    // std::cerr << "checkpoint d.2" << std::endl;

    //std::cerr << "done mibf" << std::endl;

    // iterate through gap_index_set and check if the pos is in the mibf
    // if it is, add the pos to the start_end_pos_vec
    size_t start_of_first_gap = std::get<0>(gap_vec[0]);
    size_t end_of_last_gap = std::get<1>(gap_vec[gap_vec.size() - 1]);
    std::vector<std::vector<uint32_t>> vec_of_vec_of_pos(3);
    std::vector<std::vector<size_t>> vec_of_vec_of_seq_pos(3);

    // std::cerr << "checkpoint e" << std::endl;

    //std::cerr << "pre frame" << std::endl;

    for (size_t frame = 0; frame < 3; ++frame)
    {
        btllib::AAHash aahash(sixframed_xlated_proteins[frame + ori * 3], hash_num, rescue_kmer_size, level, start_of_first_gap - 1);
        btllib::AAHash aahash2(sixframed_xlated_proteins[frame + ori * 3], hash_num, rescue_kmer_size, 2, start_of_first_gap - 1);
        btllib::AAHash aahash3(sixframed_xlated_proteins[frame + ori * 3], hash_num, rescue_kmer_size, 3, start_of_first_gap - 1);
        aahash.roll();
        aahash2.roll();
        aahash3.roll();
        //bool prev_contains = false;
        while(aahash.get_pos() <= end_of_last_gap + 1)
        {

            if (small_mi_bf.bv_contains(aahash.hashes()))
            {
                // get pos of mibf entry

                auto temp_ID_pos = small_mi_bf.get_id(aahash.hashes());
                for (auto &ID_pos : temp_ID_pos)
                {
                    auto pos = ID_pos & 0xFFFFFFFF;

                    if (gap_index_set.find(pos) != gap_index_set.end())
                    {
                        vec_of_vec_of_pos[frame].push_back(pos);
                        vec_of_vec_of_seq_pos[frame].push_back(aahash.get_pos());
                        // remove entry
                        gap_index_set.erase(pos);
                        break;
                    }
                } 
            } else if (small_mi_bf.bv_contains(aahash2.hashes()))
            {
                // get pos of mibf entry

                auto temp_ID_pos = small_mi_bf.get_id(aahash2.hashes());
                for (auto &ID_pos : temp_ID_pos)
                {
                    auto pos = ID_pos & 0xFFFFFFFF;

                    if (gap_index_set.find(pos) != gap_index_set.end())
                    {
                        vec_of_vec_of_pos[frame].push_back(pos);
                        vec_of_vec_of_seq_pos[frame].push_back(aahash2.get_pos());
                        // remove entry
                        gap_index_set.erase(pos);
                        break;
                    }
                } 
            } else if (small_mi_bf.bv_contains(aahash3.hashes()))
            {
                // get pos of mibf entry

                auto temp_ID_pos = small_mi_bf.get_id(aahash3.hashes());
                for (auto &ID_pos : temp_ID_pos)
                {
                    auto pos = ID_pos & 0xFFFFFFFF;

                    if (gap_index_set.find(pos) != gap_index_set.end())
                    {
                        vec_of_vec_of_pos[frame].push_back(pos);
                        vec_of_vec_of_seq_pos[frame].push_back(aahash3.get_pos());
                        // remove entry
                        gap_index_set.erase(pos);
                        break;
                    }
                } 
            }
            aahash.roll();
            aahash2.roll();
            aahash3.roll();
        }
    }

    // std::cerr << "checkpoint f" << std::endl;

    //std::cerr << "post frame" << std::endl;

    for (size_t i = 0; i < 3; ++i) {
        if (vec_of_vec_of_pos[i].size() == 0) {
            continue;
        }
        size_t start_of_block = 0;
        for (size_t j = 1; j < vec_of_vec_of_pos[i].size(); ++j) {
            if (vec_of_vec_of_pos[i][j] - vec_of_vec_of_pos[i][j - 1] == 1) {
                // check vec_of_vec_of_seq_pos[i][j] - vec_of_vec_of_seq_pos[i][j - 1] == 1
                if (vec_of_vec_of_seq_pos[i][j] - vec_of_vec_of_seq_pos[i][j - 1] == 1) {
                    continue;
                } else {
                    size_t length_of_block = j - 1 - start_of_block + 1;
                    adjusted_kmer_counts += length_of_block + 3;
                    start_of_block = j;
                }
            }
        }
        size_t length_of_block = vec_of_vec_of_pos[i].size() - 1 - start_of_block + 1;
        adjusted_kmer_counts += length_of_block + 3;
    }

    // std::cerr << "checkpoint g" << std::endl;

    //std::cerr << "post processing" << std::endl;
    return;

    // print out vec_of_vec_of_pos with all the pos found
    /*for (size_t frame = 0; frame < 3; ++frame)
    {
        std::cerr << "frame: " << frame << "\t";
        for (auto &pos : vec_of_vec_of_pos[frame])
        {
            std::cerr << pos << "\t";
        }
        std::cerr << std::endl;
    }*/



    /*for (auto &pos : gap_index_set)
    {
        btllib::AAHash aahash(sixframed_xlated_proteins[0 + ori * 3], hash_num, rescue_kmer_size, 1, pos);
        if (small_mi_bf.bv_contains(aahash.hashes()))
        {
            start_end_pos_vec.emplace_back(std::make_tuple(pos, pos + rescue_kmer_size));
        }
    }*/
}

bool explore_frame(btllib::MIBloomFilter<uint64_t> &mi_bf, btllib::AAHash &aahash, std::deque<std::vector<uint32_t>> &miBf_IDs_snapshot, std::deque<std::vector<uint32_t>> &miBf_pos_snapshot, std::unordered_map<uint32_t, size_t> &id_to_count, btllib::AAHash &aahash2, btllib::AAHash &aahash3)
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

    int hash_lvl = 1;
    if (!mi_bf.bv_contains(aahash.hashes()))
    {
        hash_lvl = 2;
    } else if (!mi_bf.bv_contains(aahash2.hashes()))
    {
        hash_lvl = 3;
    } else if (!mi_bf.bv_contains(aahash3.hashes()))
    {
        // clear snapshots and id_to_count
        miBf_IDs_snapshot.clear();
        miBf_pos_snapshot.clear();
        id_to_count.clear();
        return false;
    }

    // query mibf and insert into both deques
    decltype(mi_bf.get_id(aahash.hashes())) temp_ID_pos;
    if (hash_lvl == 1)
    {
        temp_ID_pos = mi_bf.get_id(aahash.hashes());
    } else if (hash_lvl == 2)
    {
        temp_ID_pos = mi_bf.get_id(aahash2.hashes());
    } else if (hash_lvl == 3)
    {
        temp_ID_pos = mi_bf.get_id(aahash3.hashes());
    }
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
        /*////// std::cerr << "check pos before return true" << std::endl;
        for (auto &pos : temp_pos_set)
        {
            //std::cerr << "pos: " << pos << std::endl;
        }*/
    }

    return false;
}
struct gff_entry {
    std::string query_name;
    size_t hit_pos_start;
    size_t hit_pos_end;
    double score;
    std::string strand;
    std::string hit_name;

    gff_entry(const std::string& query_name_, size_t hit_pos_start_, size_t hit_pos_end_, double score_, const std::string& strand_, const std::string& hit_name_)
        : query_name(query_name_), hit_pos_start(hit_pos_start_), hit_pos_end(hit_pos_end_), score(score_), strand(strand_), hit_name(hit_name_) {}
    
    gff_entry(std::string&& query_name_, size_t hit_pos_start_, size_t hit_pos_end_, double score_, std::string&& strand_, std::string&& hit_name_)
        : query_name(std::move(query_name_)), hit_pos_start(hit_pos_start_), hit_pos_end(hit_pos_end_), score(score_), strand(std::move(strand_)), hit_name(std::move(hit_name_)) {}
};

struct gff_entry_comparator {
    bool operator()(const gff_entry& lhs, const gff_entry& rhs) const noexcept {
        if (lhs.query_name == rhs.query_name) {
            return lhs.hit_pos_start < rhs.hit_pos_start;
        }
        return lhs.query_name < rhs.query_name;
    }
};



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
    size_t rescue_kmer_size = 4;
    uint64_t genome_size = 0;
    double lower_bound = 0.8;
    std::string mibf_path = "";
    static struct option long_options[] = {
        {"help", no_argument, &help_flag, 1},
        {"verbose", no_argument, &verbose_flag, 1},
        {"threads", required_argument, 0, 't'},
        {"input", required_argument, 0, 'i'},
        {"reference", required_argument, 0, 'r'},
        {"mibf_path", required_argument, 0, 'm'},
        {"output", required_argument, 0, 'o'},
        {"hash", required_argument, 0, 'h'},
        {"kmer", required_argument, 0, 'k'},
        {"genome", required_argument, 0, 'g'},
        {"lower_bound", required_argument, 0, 'l'},
        {0, 0, 0, 0}};

    // loop through command line arguments
    while ((opt = getopt_long(argc, argv, "g:h:i:t:o:r:k:m:l:", long_options, &option_index)) != -1)
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
        case 'm':
            mibf_path = optarg;
            break;
        case 'g':
            genome_size = std::stoul(optarg);
            break;
        case 'l':
            lower_bound = std::stod(optarg);
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
                  << "Rescue kmer size: " << (uint64_t)rescue_kmer_size << "\n"
                  << "Genome size: " << genome_size << "\n"
                  << "Lower bound: " << lower_bound << "\n" 
                  << std::endl;
    }

    omp_set_num_threads(threads);


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

    //btllib::MIBloomFilter<uint64_t> mi_bf(calc_optimal_size(std::max<size_t>(genome_size * 3, 1000000), hash_num, 0.1), hash_num);
    // btllib::MIBloomFilter<uint64_t> mi_bf(calc_optimal_size(1000000000, hash_num, 0.1), hash_num);

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
            if (record.seq.size() == 0)
            {
                std::cerr << "seq name: " << record.id << std::endl;
                exit(0);
            }
            /*//std::cerr << "seq name: " << record.id << " miBf_ID: " << miBf_ID << std::endl;
            //std::cerr << "seq size: " << record.seq.size() << std::endl;
            //std::cerr << "seq: " << record.seq << std::endl; */
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
    /*for (int stage = 0; stage < 3; stage++)
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
            //while (itr3.roll())
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
    //std::cerr << "Pop sat cnt: " << mi_bf.get_pop_saturated_cnt() << std::endl;
    //exit(0);
    */
    std::cerr << "finished making MiBF" << std::endl;
    std::cerr << "in " << std::setprecision(4) << std::fixed << omp_get_wtime() - sTime
              << "\n";
    
    

    std::vector<std::ofstream> output_files(3);
    std::vector<std::ofstream> gff_files(3);
    std::vector<std::ofstream> pre_gff_files(3);

    // Open each output file stream with a unique filename based on output_prefix
    for (int i = 0; i < 3; ++i) {
        std::string filename = output_prefix + "_lvl" + std::to_string(i + 1) + ".results.tsv";
        output_files[i].open(filename);
        output_files[i] << "name\tcomplete copies\tincomplete copies\texpected k-mer counts\thighest adjusted incomplete k-mer hits" << std::endl;
        gff_files[i].open(output_prefix + "_lvl" + std::to_string(i + 1) + ".gff");
        gff_files[i] << "##gff-version 3" << std::endl;
        pre_gff_files[i].open(output_prefix + "_lvl" + std::to_string(i + 1) + ".pre.gff");
        pre_gff_files[i] << "##gff-version 3" << std::endl;
    }


    // make a gff set sorted by seq name and start pos
    // the columns are seq name, start pos, end pos, score, strand,

    // source, type missing because always the same, no attribute,phase

    // compartor for gff set, sorted by seq name and start pos
    
    /*auto gff_comparator = [](const std::tuple<std::string, size_t, size_t, double, std::string, std::string> &a, const std::tuple<std::string, size_t, size_t, double, std::string, std::string> &b)
    {
        if (std::get<0>(a) == std::get<0>(b))
        {
            return std::get<1>(a) < std::get<1>(b);
        }
        return std::get<0>(a) < std::get<0>(b);
    };

    std::vector<std::set<std::tuple<std::string, size_t, size_t, double, std::string, std::string>, decltype(gff_comparator)>> gff_set_vector;*/

    // Create and insert three different instances of gff_set into the vector
    std::vector<std::set<gff_entry, gff_entry_comparator>> gff_set_vector;
    for (int i = 0; i < 3; ++i) {
        gff_set_vector.emplace_back(gff_entry_comparator());
    }

   //std::vector<std::set<std::tuple<std::string, size_t, size_t, double, std::string, std::string>, decltype(gff_comparator)>> pre_gff_set_vector;
    std::vector<std::set<gff_entry, gff_entry_comparator>> pre_gff_set_vector;
    for (int i = 0; i < 3; ++i) {
        pre_gff_set_vector.emplace_back(gff_entry_comparator());
    }


    //auto& gff_set = gff_set_vector[0];

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
        size_t highest_adjusted_kmer_counts = 0;
    };
    std::vector<std::unordered_map<std::string, completeness_struct>> seq_name_to_completeness_vec(3);

    // Populate each map in the vector with empty entries
    for (auto &seq_name_to_completeness : seq_name_to_completeness_vec) {
        for (const auto &seq_ID : seq_ID_to_miBf_ID) {
            seq_name_to_completeness[seq_ID.first] = completeness_struct();
        }
    }

#pragma omp parallel num_threads(threads / 2)
    for (const auto record : reader)
    {
        // //std::cerr << "seq name: " << record.id << std::endl;
        std::vector<std::string> sixframed_xlated_proteins = sixframe_translate(record.seq);
        /*std::cerr << "length of original sequence: " << record.seq.size() << std::endl;
        std::cerr << "length of protein 1: " << sixframed_xlated_proteins[0].size() << std::endl;
        std::cerr << "length of protein 2: " << sixframed_xlated_proteins[1].size() << std::endl;
        std::cerr << "length of protein 3: " << sixframed_xlated_proteins[2].size() << std::endl;
*/
        //create a protein fa file and write all proteins to it
        /*std::ofstream protein_fa_file(output_prefix + record.id + ".fa");
        protein_fa_file << ">" << record.id << "_1" << std::endl;
        protein_fa_file << sixframed_xlated_proteins[0] << std::endl;
        protein_fa_file << ">" << record.id << "_2" << std::endl;
        protein_fa_file << sixframed_xlated_proteins[1] << std::endl;
        protein_fa_file << ">" << record.id << "_3" << std::endl;
        protein_fa_file << sixframed_xlated_proteins[2] << std::endl;
        protein_fa_file << ">" << record.id << "_4" << std::endl;
        protein_fa_file << sixframed_xlated_proteins[3] << std::endl;
        protein_fa_file << ">" << record.id << "_5" << std::endl;
        protein_fa_file << sixframed_xlated_proteins[4] << std::endl;
        protein_fa_file << ">" << record.id << "_6" << std::endl;
        protein_fa_file << sixframed_xlated_proteins[5] << std::endl;
        protein_fa_file.close();*/

        // //std::cerr << "protein 4: " << sixframed_xlated_proteins[3] << std::endl;
        // //std::cerr << "protein 5: " << sixframed_xlated_proteins[4] << std::endl;
        // //std::cerr << "protein 6: " << sixframed_xlated_proteins[5] << std::endl;
#pragma omp parallel for num_threads(2)
        for (size_t ori = 0; ori < 2; ++ori)
        //for (size_t ori = 0; ori < 2; ++ori) //TODO
        {
            // frame to block id to id and smallest pos and largest pos
            std::vector<std::unordered_map<size_t, std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>>>> frame_to_block_id_to_id_and_pos_vec(3);
            // id to count across all frames sorted by count largest to smallest
            std::vector<std::map<uint32_t, size_t, std::greater<size_t>>> id_to_count_across_all_frames_vec(3);
            // custom comparator for set of pair of size_t and pair of size_t and size_t to sort by seq pos
            /*auto custom_comparator = [](const std::tuple<size_t, size_t, size_t> &a, const std::tuple<size_t, size_t, size_t> &b)
            {
                return std::get<2>(a) < std::get<2>(b);
            };*/
            // id to set of frame and block id and seq pos, set is sorted by seq pos
            //std::vector<std::unordered_map<uint32_t, std::set<std::tuple<size_t, size_t, size_t>, CustomComparator>>> id_to_frame_block_id_and_seq_pos_vec(3);
            std::vector<std::unordered_map<uint32_t, std::set<frame_block, frame_block_comparator>>> id_to_frame_block_id_and_seq_pos_vec(3);
            //for (size_t curr_lvl = 1; curr_lvl <= 3; ++curr_lvl) {
            for (size_t curr_lvl = 1; curr_lvl <= 1; ++curr_lvl) {
                auto& frame_to_block_id_to_id_and_pos = frame_to_block_id_to_id_and_pos_vec[curr_lvl - 1];
                auto& id_to_count_across_all_frames = id_to_count_across_all_frames_vec[curr_lvl - 1];
                auto& id_to_frame_block_id_and_seq_pos = id_to_frame_block_id_and_seq_pos_vec[curr_lvl - 1];
                auto& gff_set = gff_set_vector[curr_lvl - 1];
                auto& pre_gff_set = pre_gff_set_vector[curr_lvl - 1];
                auto& seq_name_to_completeness = seq_name_to_completeness_vec[curr_lvl - 1];
                for (size_t frame = 0; frame < 3; ++frame)
                //for (size_t frame = 0; frame < 1; ++frame)
                {

                    /*if (verbose_flag) {
                        std::cerr << "frame: " << frame << " ori: " << ori << std::endl;
                        btllib::AAHash aahash_test(sixframed_xlated_proteins[frame + ori * 3], hash_num, kmer_size, 1);
                        btllib::AAHash aahash_test2(sixframed_xlated_proteins[frame + ori * 3], hash_num, kmer_size, 2);
                        aahash_test.roll();
                        aahash_test2.roll();
                        for (size_t i = 0; i < 100; ++i) {
                            if (mi_bf.bv_contains(aahash_test.hashes())) {
                                std::cerr << "lvl 1" << std::endl;
                                std::cerr << "bv contains" << i << " th hash" << std::endl;
                                // get id and pos
                                auto temp_ID_pos = mi_bf.get_id(aahash_test.hashes());
                                for (auto &ID_pos : temp_ID_pos) {
                                    std::cerr << "ID: " << (ID_pos >> 32) << " pos: " << (ID_pos & 0xFFFFFFFF) << std::endl;
                                }
                            }
                            if (mi_bf.bv_contains(aahash_test2.hashes())) {
                                std::cerr << "lvl 2" << std::endl;
                                std::cerr << "bv contains" << i << " th hash" << std::endl;
                                // get id and pos
                                auto temp_ID_pos = mi_bf.get_id(aahash_test2.hashes());
                                for (auto &ID_pos : temp_ID_pos) {
                                    std::cerr << "ID: " << (ID_pos >> 32) << " pos: " << (ID_pos & 0xFFFFFFFF) << std::endl;
                                }
                            }
                            aahash_test.roll();
                            aahash_test2.roll();
                        }

                    }*/
                    

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
                        //// std::cerr << "checkpoint 1" << std::endl;
                        while (!explore_frame(mi_bf, aahash, miBf_IDs_snapshot, miBf_pos_snapshot, id_to_count, aahash2, aahash3) && aahash.get_pos() != std::numeric_limits<size_t>::max())
                        {
                            // //std::cerr << "explore_frame returned false" << std::endl;
                            aahash.roll();
                            aahash2.roll();
                            aahash3.roll();
                        }
                        //// std::cerr << "checkpoint 2" << std::endl;
                        if (aahash.get_pos() == std::numeric_limits<size_t>::max())
                        {
                            break;
                        }
                        //// std::cerr << "checkpoint 3" << std::endl;
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
                            id_to_frame_block_id_and_seq_pos[ID_pos_set.first].emplace(frame, block_id, seq_pos);

    #pragma omp critical
                            {

                                if (id_to_count_across_all_frames.find(ID_pos_set.first) == id_to_count_across_all_frames.end())
                                {
                                    //// std::cerr << "checkpoint x.8" << std::endl;
                                    id_to_count_across_all_frames[ID_pos_set.first] = ID_pos_set.second.size();
                                    //// std::cerr << "checkpoint x.9" << std::endl;
                                }
                                else
                                {
                                    //// std::cerr << "checkpoint x.10" << std::endl;
                                    id_to_count_across_all_frames[ID_pos_set.first] += ID_pos_set.second.size();
                                    //// std::cerr << "checkpoint x.11" << std::endl;
                                }
                            }
                            ++block_id;
                        }
                        //// std::cerr << "checkpoint y" << std::endl;
                        

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
                    // print id_to_frame_block_id_and_seq_pos
                    for (auto &ID_frame_block_id_seq_pos : id_to_frame_block_id_and_seq_pos)
                    {
                        std::cerr << "ID: " << ID_frame_block_id_seq_pos.first << std::endl;
                        std::cerr << "name: " << miBf_ID_to_seq_ID_and_len[ID_frame_block_id_seq_pos.first].first << std::endl;
                        for (auto &frame_block_id_seq_pos : ID_frame_block_id_seq_pos.second)
                        {
                            std::cerr << "frame: " << frame_block_id_seq_pos.frame << " block_id: " << frame_block_id_seq_pos.block_id << " seq_pos: " << frame_block_id_seq_pos.query_start_in_prot_space << std::endl;
                        }
                    }
                    // print frame_to_block_id_to_id_and_pos
                    for (auto &frame_block_id_to_id_and_pos : frame_to_block_id_to_id_and_pos)
                    {
                        std::cerr << "frame: " << frame_block_id_to_id_and_pos.first << std::endl;
                        for (auto &block_id_to_id_and_pos : frame_block_id_to_id_and_pos.second)
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
                //std::cerr << "finished processing frame " << curr_lvl << " " << ori << std::endl;

                for (auto &ID_count : id_to_count_across_all_frames)
                {
                    uint32_t miBf_ID = ID_count.first;
                    // if mibf id is not in miBf_ID_to_seq_ID_and_len, continue
                    if (miBf_ID_to_seq_ID_and_len.find(miBf_ID) == miBf_ID_to_seq_ID_and_len.end())
                    {
                        continue;
                    }
                    /*if (verbose_flag){
                        std::cerr << "ID: " << ID_count.first << " count: " << ID_count.second << std::endl;
                    }*/
                    
                    // iterate through id to frame block id and seq pos and log the completeness
                    std::string seq_name = miBf_ID_to_seq_ID_and_len[miBf_ID].first;
                    /*if (verbose_flag) {
                        std::cerr << "calculating for protein name: " << seq_name << std::endl;
                    }*/
                    size_t complete_copies = 0;
                    size_t incomplete_copies = 0;
                    size_t expected_kmer_counts = miBf_ID_to_seq_ID_and_len[miBf_ID].second - kmer_size + 1;
                    /*if (verbose_flag) {
                        std::cerr << "protein length: " << miBf_ID_to_seq_ID_and_len[miBf_ID].second << std::endl;
                    }*/
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
                    //std::cerr << "calculating for protein name: " << seq_name << std::endl;

                    
                    std::vector<std::reference_wrapper<const frame_block>> vec;
                    for (auto &frame_block_id_and_seq_pos : id_to_frame_block_id_and_seq_pos[miBf_ID])
                    {
                        vec.push_back(frame_block_id_and_seq_pos);
                    }
                    for (size_t ref_idx = 0; ref_idx < vec.size(); ++ref_idx)
                    {
                        const frame_block &frame_block_id_and_seq_pos = vec[ref_idx].get();

                        if (verbose_flag) {
                            //std::cerr << "frame: " << std::get<0>(frame_block_id_and_seq_pos) << " block_id: " << std::get<1>(frame_block_id_and_seq_pos) << " seq_pos: " << std::get<2>(frame_block_id_and_seq_pos) << std::endl;
                            //std::cerr << "curr end_pos: " << end_pos << std::endl;
                            //std::cerr << "next start_pos: " << frame_to_block_id_to_id_and_pos[std::get<0>(frame_block_id_and_seq_pos)][std::get<1>(frame_block_id_and_seq_pos)].first << std::endl;
                            //std::cerr << std::endl;
                        }

                        if (frame == 3) {
                            frame = frame_block_id_and_seq_pos.frame;
                            seq_start_in_nucleotide = frame_block_id_and_seq_pos.query_start_in_prot_space * 3 + frame;
                        } else {
                            frame = frame_block_id_and_seq_pos.frame;
                        }

                        size_t block_id = frame_block_id_and_seq_pos.block_id;
                        // std::cerr << "checkpoint 1" << std::endl;



                        if (frame_to_block_id_to_id_and_pos[frame].find(block_id) != frame_to_block_id_to_id_and_pos[frame].end())
                        {

                        
                            prev_block_start = block_start;
                            block_start = frame_block_id_and_seq_pos.query_start_in_prot_space;
                            prev_block_len = block_len;
                            block_len = frame_to_block_id_to_id_and_pos[frame][block_id].second - frame_to_block_id_to_id_and_pos[frame][block_id].first + 1;
                            
                            // std::cerr << "checkpoint 2" << std::endl;

                            if (start_end_pos_vec.empty()) {
                                start_end_pos_vec.emplace_back(std::make_tuple(frame_block_id_and_seq_pos.query_start_in_prot_space, frame_block_id_and_seq_pos.query_start_in_prot_space + block_len - 1));\
                                //std::cerr << "empty" << std::endl;
                                //std::cerr << "start end pos vec start: " << std::get<2>(frame_block_id_and_seq_pos) << std::endl;
                                //std::cerr << "start end pos vec end: " << std::get<2>(frame_block_id_and_seq_pos) + block_len - 1 << std::endl;
                                start_end_pos_tar_vec.emplace_back(std::make_tuple(frame_to_block_id_to_id_and_pos[frame][block_id].first, frame_to_block_id_to_id_and_pos[frame][block_id].second + 1));
                            }

                            // std::cerr << "checkpoint 3" << std::endl;

                            if (end_pos == 0)
                            {
                                
                                end_pos = frame_to_block_id_to_id_and_pos[frame][block_id].second;
                                if (verbose_flag){
                                    std::cerr << "start end_pos: " << end_pos << std::endl;
                                    std::cerr << std::endl;
                                }
                                
                                
                                adjusted_kmer_counts = block_len;
                                // std::cerr << "checkpoint 3.1" << std::endl;
                            }
                            else
                            {
                                if (end_pos < frame_to_block_id_to_id_and_pos[frame][block_id].first)
                                {
                                //std::cerr << "nonempty" << std::endl;
                                //std::cerr << "start end pos vec start: " << std::get<2>(frame_block_id_and_seq_pos) << std::endl;
                                //std::cerr << "start end pos vec end: " << std::get<2>(frame_block_id_and_seq_pos) + block_len - 1 << std::endl;
                                    start_end_pos_vec.emplace_back(std::make_tuple(frame_block_id_and_seq_pos.query_start_in_prot_space, frame_block_id_and_seq_pos.query_start_in_prot_space + block_len - 1));
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
                                    // std::cerr << "checkpoint 3.2.1" << std::endl;

                                    // look ahead code
                                    if (ref_idx + 2 < vec.size()) {
                                        const frame_block &next_frame_block_id_and_seq_pos = vec[ref_idx + 1].get();
                                        auto next_block_id = next_frame_block_id_and_seq_pos.block_id;
                                        auto next_frame = next_frame_block_id_and_seq_pos.frame;
                                        const frame_block &next_next_frame_block_id_and_seq_pos = vec[ref_idx + 2].get();
                                        auto next_next_block_id = next_next_frame_block_id_and_seq_pos.block_id;
                                        auto next_next_frame = next_next_frame_block_id_and_seq_pos.frame;
                                        if (end_pos < frame_to_block_id_to_id_and_pos[next_frame][next_block_id].first) {
                                            continue;
                                        } else {
                                            if (end_pos < frame_to_block_id_to_id_and_pos[next_next_frame][next_next_block_id].first) {
                                                ref_idx += 1;
                                            } else {
                                                continue;
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    seq_end_in_nucleotide = (prev_block_start + prev_block_len + kmer_size - 1)* 3 + frame;
                                    // log completeness and reset
                                    double prev_score = (double)adjusted_kmer_counts / (double)expected_kmer_counts;
                                    // std::cerr << "checkpoint 3.2.2.1" << std::endl;

                                    if (adjusted_kmer_counts > lower_bound * expected_kmer_counts)
                                    {
                                        /*if (adjusted_kmer_counts >= 0.95 * expected_kmer_counts)
                                        {
                                            complete_copies++;
                                        }*/
                                        //std::cerr << "attempting to gap fill1" << std::endl;
                                        //std::cerr << "mibf_id: " << miBf_ID << std::endl;
                                        if (start_end_pos_vec.size() > 1) {
                                            fill_in_gaps(start_end_pos_vec, start_end_pos_tar_vec, adjusted_kmer_counts, miBf_ID_to_seq[miBf_ID], hash_num, rescue_kmer_size, sixframed_xlated_proteins, ori, kmer_size, curr_lvl, threads);
                                        }
                                        //fill_in_gaps(start_end_pos_vec, start_end_pos_tar_vec, adjusted_kmer_counts, miBf_ID_to_seq[miBf_ID], hash_num, rescue_kmer_size, sixframed_xlated_proteins, ori, kmer_size);
                                        //std::cerr << "done gap fill" << std::endl;
                                        if (adjusted_kmer_counts > 0.95 * expected_kmer_counts) {
                                            complete_copies++;
                                        } else {
                                            incomplete_copies++;
                                        }
                                    }
                                    // std::cerr << "checkpoint 3.2.2.2" << std::endl;
                                    double score = (double)adjusted_kmer_counts / (double)expected_kmer_counts;
                                    //if (adjusted_kmer_counts > lower_bound * expected_kmer_counts)
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
                                    // std::cerr << "checkpoint 3.2.2.3" << std::endl;
                                    
                                    end_pos = frame_to_block_id_to_id_and_pos[frame][block_id].second;
                                    adjusted_kmer_counts = block_len;
                                    seq_start_in_nucleotide = frame_block_id_and_seq_pos.query_start_in_prot_space * 3 + frame;
                                    // clear start_end_pos_vec and start_end_pos_tar_vec
                                    start_end_pos_vec.clear();
                                    start_end_pos_tar_vec.clear();
                                    start_end_pos_vec.emplace_back(std::make_tuple(frame_block_id_and_seq_pos.query_start_in_prot_space, frame_block_id_and_seq_pos.query_start_in_prot_space + block_len - 1));
                                    start_end_pos_tar_vec.emplace_back(std::make_tuple(frame_to_block_id_to_id_and_pos[frame][block_id].first, frame_to_block_id_to_id_and_pos[frame][block_id].second + 1));
                                    //frame = 3; //bug fix problem with seq_start in nucleotide, 1 block off
                                    // std::cerr << "checkpoint 3.2.2.4" << std::endl;
                                    if (verbose_flag){
                                        std::cerr << "new end_pos: " << end_pos << std::endl;
                                    }
                                }
                            }
                            // block_len = frame_to_block_id_to_id_and_pos[frame][block_id].second - frame_to_block_id_to_id_and_pos[frame][block_id].first + 1;
                        }
                    }
                    double prev_score = (double)adjusted_kmer_counts / (double)expected_kmer_counts;
                    // std::cerr << "checkpoint 4" << std::endl;
                    /*if (adjusted_kmer_counts >= 0.95 * expected_kmer_counts)
                    {
                        complete_copies++;
                    }
                    else*/ if (adjusted_kmer_counts > lower_bound * expected_kmer_counts)
                    {


                        //std::cerr << "attempting to gap fill2" << std::endl;
                        //std::cerr << "adjusted_kmer_counts: " << adjusted_kmer_counts << std::endl;
                        if (start_end_pos_vec.size() > 1) {
                            fill_in_gaps(start_end_pos_vec, start_end_pos_tar_vec, adjusted_kmer_counts, miBf_ID_to_seq[miBf_ID], hash_num, rescue_kmer_size, sixframed_xlated_proteins, ori, kmer_size, curr_lvl, threads);
                        }
                        //fill_in_gaps(start_end_pos_vec, start_end_pos_tar_vec, adjusted_kmer_counts, miBf_ID_to_seq[miBf_ID], hash_num, rescue_kmer_size, sixframed_xlated_proteins, ori, kmer_size);
                        //std::cerr << "done fill2" << std::endl;
                        //std::cerr << "adjusted_kmer_counts: " << adjusted_kmer_counts << std::endl;
                        if (adjusted_kmer_counts > 0.95 * expected_kmer_counts) {
                            complete_copies++;
                        } else {
                            incomplete_copies++;
                        }

                        // std::cerr << "checkpoint 4.1" << std::endl;
                    }
                    // log completeness
#pragma omp atomic
                    seq_name_to_completeness[seq_name].complete_copies += complete_copies;

#pragma omp atomic
                    seq_name_to_completeness[seq_name].incomplete_copies += incomplete_copies;

                    //if (adjusted_kmer_counts > lower_bound * expected_kmer_counts)
                    {

                        const auto &last_frame_block_id_and_seq_pos = id_to_frame_block_id_and_seq_pos[miBf_ID].rbegin();
                        frame = (*last_frame_block_id_and_seq_pos).frame;
                        seq_end_in_nucleotide = ((*last_frame_block_id_and_seq_pos).query_start_in_prot_space + block_len + kmer_size - 1) * 3 + frame;
                        // print all the values that make seq_end_in_nucleotide
                        double score = (double)adjusted_kmer_counts / (double)expected_kmer_counts;
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
                        // std::cerr << "checkpoint 4.2" << std::endl;
                    }
                }
                // //std::cerr << "done with ori" << std::endl;
            }
        }
    }

    // //std::cerr << "done processing" << std::endl;

    // output the completeness to the output file
    for (int i = 0; i < 3; ++i) {
        auto& seq_name_to_completeness = seq_name_to_completeness_vec[i];
        for (auto &seq_name_completeness : seq_name_to_completeness)
        {
            output_files[i] << seq_name_completeness.first << "\t" << seq_name_completeness.second.complete_copies << "\t" << seq_name_completeness.second.incomplete_copies << "\t" << seq_name_completeness.second.expected_kmer_counts << "\t" << seq_name_completeness.second.highest_adjusted_kmer_counts << std::endl;
        }
    }

    // output the gff set to the gff file
    for (int i = 0; i < 3; ++i) {
        auto& gff_set = gff_set_vector[i];
        for (auto &gff : gff_set)
        {
            gff_files[i] << gff.query_name << "\t"
                        << "."
                        << "\t"
                        << "gene"
                        << "\t" << gff.hit_pos_start << "\t" << gff.hit_pos_end << "\t" << gff.score << "\t" << gff.strand << "\t"
                        << "0"
                        << "\t"
                        << "ID=" << gff.hit_name << std::endl;
        }
        auto& pre_gff_set = pre_gff_set_vector[i];
        for (auto &gff : pre_gff_set)
        {
            pre_gff_files[i] << gff.query_name << "\t"
                        << "."
                        << "\t"
                        << "gene"
                        << "\t" << gff.hit_pos_start << "\t" << gff.hit_pos_end << "\t" << gff.score << "\t" << gff.strand << "\t"
                        << "0"
                        << "\t"
                        << "ID=" << gff.hit_name << std::endl;
        }
    }

    /*for (auto &seq_name_completeness : seq_name_to_completeness)
    {
        // output_file << seq_name_completeness.first << "\t" << std::get<0>(seq_name_completeness.second) << "\t" << std::get<1>(seq_name_completeness.second) << "\t" << std::get<2>(seq_name_completeness.second) << "\t" << std::get<3>(seq_name_completeness.second) << std::endl;
        output_file << seq_name_completeness.first << "\t" << seq_name_completeness.second.complete_copies << "\t" << seq_name_completeness.second.incomplete_copies << "\t" << seq_name_completeness.second.expected_kmer_counts << "\t" << seq_name_completeness.second.highest_adjusted_kmer_counts << std::endl;
    }*/

    /*for (auto &gff : gff_set)
    {
        gff_output_file << std::get<0>(gff) << "\t"
                        << "."
                        << "\t"
                        << "gene"
                        << "\t" << std::get<1>(gff) << "\t" << std::get<2>(gff) << "\t" << std::get<3>(gff) << "\t" << std::get<4>(gff) << "\t"
                        << "0"
                        << "\t"
                        << "ID=" << std::get<5>(gff) << std::endl;
    }*/


    return 0;
}
