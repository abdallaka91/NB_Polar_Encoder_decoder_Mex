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
                              const bool newsim)
{
    std::filesystem::path filepath(filename);
    std::filesystem::create_directories(filepath.parent_path());

    std::ofstream file(filename, newsim ? std::ios::out : std::ios::app);

    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return false;
    }

    std::string dashes = "layer " + std::to_string(l) + ", " + "cluster " + std::to_string(s);
    if (!newsim)
    //     file << dashes << "\n";
    // else
        file << "\n"
             << dashes << "\n";

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



void writeBtMatrices(const std::string &filename, int n, int nH, int nL,
                     const std::vector<std::vector<std::vector<std::vector<int>>>> &Bt)
{
    std::filesystem::create_directories(std::filesystem::path(filename).parent_path());

    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    for (int l = 0; l < n; l++)
    {
        for (int s = 0; s < (1u << l); s++)
        {
            std::ostringstream lne;
            lne << l << " " << s << ",  ";
            for (int j0 = 0; j0 < nH; j0++)
            {
                for (int j1 = 0; j1 < nL; j1++)
                {
                    if (Bt[l][s][j0][j1] == 1)
                    {
                        lne << j0 << "  " << j1 << std::setw(8);
                    }
                }
            }
            if (!(l == n - 1 && s == (1u << l) - 1))
            {
                lne << "\n";
            }
            file << lne.str();
        }
        file << "\n";
    }

    std::cout << "Bt Matrices written to: " << filename << std::endl;
}

#endif // HELPERFUNC_H
