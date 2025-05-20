#ifndef HELPERFUNC_H
#define HELPERFUNC_H

#include <fstream>
#include <iostream>
#include <vector>
#include <filesystem>

template <typename T>
bool AppendClustBubblesToFile(const std::string &filename,
                              const std::vector<std::vector<T>> &table,
                              const int l,
                              const int s,
                              const bool newsim,
                              const std::string &head)
{
    std::filesystem::path filepath(filename);
    std::filesystem::create_directories(filepath.parent_path());

    std::ofstream file(filename, newsim ? std::ios::out : std::ios::app);
    if(newsim)
        file << head;


    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return false;
    }
    std::string dashes = "layer " + std::to_string(l) + ", " + "cluster " + std::to_string(s);
    if (!newsim)
    file << "\n"<< dashes << "\n";
    else
    file << dashes << "\n";

    for (const auto &row : table)
    {
        for (const auto &value : row)
        {
            file << value << " ";
        }
        file << "\n";
    }
    file.close();
    return true;
}

inline void save_intrinsic_LLR(const std::vector<std::vector<PoAwN::structures::softdata_t>> &chan_LLR, const std::string &filename)
{
    std::ofstream file(filename, std::ios::app);
    for (const auto &symbol_LLR : chan_LLR)
    {
        for (const auto &llr : symbol_LLR)
            file << llr << " ";
        file << "\n";
    }
}

template <typename T>
bool appendSequenceToFile(const std::string &filename, const std::vector<T> &sequence)
{
    std::ofstream file(filename, std::ios::app);
    if (!file.is_open())
        return false;

    for (const auto &symbol : sequence)
        file << symbol << " ";
    file << "\n";

    return true;
}

#endif // HELPERFUNC_H
