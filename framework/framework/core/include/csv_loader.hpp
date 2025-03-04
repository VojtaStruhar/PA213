// ################################################################################
// Common Framework for Computer Graphics Courses at FI MUNI.
// 
// Copyright (c) Visitlab (https://visitlab.fi.muni.cz)
// All rights reserved.
// ################################################################################

#pragma once
#include <string>
#include <sstream>
#include <filesystem>
#include <fstream>
#include <vector>

/** The basic CSV loader. */
class CSVLoader {

    // ----------------------------------------------------------------------------
    // Static Methods
    // ----------------------------------------------------------------------------
public:
    /**
     * Loads the raster representations from csv files.
     *
     * @param 	path	The path to the csv file.
     * @return	A loaded raster representation.
     */
    static std::vector<std::vector<float>> from_csv(std::filesystem::path path) {
        std::vector<std::vector<float>> raster_vector;
        std::string word, line;
        std::fstream fin;

        fin.open(path);

        while (std::getline(fin, line)) {
            std::vector<float> tmp;
            std::stringstream s(line);
            while (std::getline(s, word, ',')) {
                tmp.push_back(std::stof(word));
            }
            raster_vector.push_back(tmp);
        }

        fin.close();
        return raster_vector;
    }

    /**
     * Loads the specific row from a csv file.
     * Saves the result as a vector of floats.
     *
     * @param 	path	The path to the csv file.
     * @param   row     Row of the file.
     * @return	A vector of numbers.
     */
    static std::vector<float> read_csv_row(std::filesystem::path path, int row) {
        std::vector<float> result;
        std::string word, line;
        std::fstream fin;

        fin.open(path);
        int count = 0;
        while (std::getline(fin, line)) {
            if (count++ == row) {
                std::stringstream s(line);
                while (std::getline(s, word, ',')) {
                    result.push_back((float)std::stod(word)); // We are using stod to avoid possible problems with long numbers.
                }
            }
        }
        fin.close();
        return result;
    }
};
