/*
 * Copyright 2015 Fadri Furrer, ASL, ETH Zurich, Switzerland
 * Copyright 2015 Michael Burri, ASL, ETH Zurich, Switzerland
 * Copyright 2015 Mina Kamel, ASL, ETH Zurich, Switzerland
 * Copyright 2015 Janosch Nikolic, ASL, ETH Zurich, Switzerland
 * Copyright 2015 Markus Achtelik, ASL, ETH Zurich, Switzerland
 * Copyright 2016 Geoffrey Hunter <gbmhunter@gmail.com>
 * Copyright 2022 He Bai, CoRAL, OSU, USA
 * Copyright 2022 Asma Tabassum , CoRAL, OSU, USA
 * Copyright 2022 Max DeSantis, CoRAL, OSU, USA
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

//TODO Double check licensing before submitting anything anywhere

#ifndef airsim_core_WindModule_hpp
#define airsim_core_WindModule_hpp

#include "common/Common.hpp"
#include "physics/PhysicsEngineBase.hpp"
#include "common/LogFileWriter.hpp"
#include "common/common_utils/Timer.hpp"
#include <chrono>
#include <mutex>

namespace msr
{
namespace airlib
{

    class WindModule
    {
    public:
        // Initiatializes WindModule with details from Airsim settings.json.
        WindModule(const std::string dataPath, const float updateInterval, const Vector3r defaultWind, const bool logEnable)
        {
            WIND_FOLDER_PATH = dataPath;            // Location of wind data files.
            UPDATE_INTERVAL = updateInterval;       // How frequently should new wind be taken from folder. Unit is in seconds.
            DEFAULT_WIND = defaultWind;             // Default wind value if vehicle exits airspace, or no wind is loaded. Vector3r.
            LOG_ENABLED = logEnable;                // Enables/disables the WindModules logging. Primarily for development.


            wind_module_alive = true;
            WindLog("Initializing...");
            // Wind data path was left default "". Only use global wind moving forward.
            if (WIND_FOLDER_PATH.empty()) {
                useSpatialWind = false;
                WindLog(Utils::stringf(" Using global wind value of (%f, %f, %f).", DEFAULT_WIND.x(), DEFAULT_WIND.y(), DEFAULT_WIND.z()));
            }
            // Wind path was specified, attempt to read a file there.
            else { 
                if (!readWindField(WIND_FOLDER_PATH + "wind_001.txt")) {
                    WindLog("Unable to open wind_001.txt at povided wind data path. Reverting to global wind.", Utils::kLogLevelWarn);
                    useSpatialWind = false;
                }
                // Wind path was found and successfully read. Load new wind field information and spin up temporal-updating thread.
                else { 
                    WindLog("First wind data file successfully read.");
                    updateWindField();
                    useSpatialWind = true;

                    // Spin up temporal thread, attempt to read more data files.
                    temporalThread = std::thread(&WindModule::temporalUpdate, this);
                }
            }
        }

        // Safely stop threads, reads, etc. before shutting down.
        ~WindModule()
        {

            WindLog("Shutting down...");

            wind_module_alive = false;          // Shuts down readFile
            temporal_thread_alive = false;      // Shuts down temporal thread
            if (fin.is_open()) fin.close();
            if (temporalThread.joinable()) temporalThread.join(); // Wait for temporal thread to close before continuing

            WindLog("... Successfully shut down.");
        }

        // Sets global wind if not using spatial-temporal wind data, or default (e.g. out of bounds) wind for spatial-temporal. Called by FastphysicsEngine.
        void setDefaultWind(const Vector3r& wind)
        {
            DEFAULT_WIND = wind;
        }

        // Return wind velocity vector at given position in NED
        const Vector3r getLocalWind(const Vector3r& position)
        {

            try {
                if (!useSpatialWind) { // Return default wind if not using spatial wind
                    return DEFAULT_WIND;
                }

                // If temporal thread is updating data, return the previously-calculated wind.
                if(!dataMutex.try_lock()) {
                    return PREVIOUS_WIND;
                }

                Vector3r wind_velocity = Vector3r::Zero();
                Vector3r posTemp;

                // Convert NED to XYZ
                posTemp.x() = position.y();
                posTemp.y() = position.x();
                posTemp.z() = (-1) * position.z();

                // Return default wind if vehicle is outside of wind field
                if (!(posTemp.x() <= max_x_) || !(posTemp.x() >= min_x_) || !(posTemp.y() <= max_y_) || !(posTemp.y() >= min_y_)) {
                    dataMutex.unlock();
                    return DEFAULT_WIND;
                }

                std::size_t x_inf = (size_t)floor((posTemp.x() - min_x_) / res_x_);
                std::size_t y_inf = (size_t)floor((posTemp.y() - min_y_) / res_y_);

                if (x_inf == n_x_ - 1u) x_inf = n_x_ - 2u;
                if (y_inf == n_y_ - 1u) y_inf = n_y_ - 2u;

                std::size_t x_sup = x_inf + 1u;
                std::size_t y_sup = y_inf + 1u;

                constexpr unsigned int n_vertices = 8; // 8 vertices on a cube to do some form of interpolation
                std::size_t idx_x[n_vertices] = { x_inf, x_inf, x_sup, x_sup, x_inf, x_inf, x_sup, x_sup };
                std::size_t idx_y[n_vertices] = { y_inf, y_inf, y_inf, y_inf, y_sup, y_sup, y_sup, y_sup };

                constexpr unsigned int n_columns = 4; // Get vertical factor from four surrounding columns
                float vertical_factors_columns[n_columns];
                for (std::size_t i = 0u; i < n_columns; ++i) {
                    vertical_factors_columns[i] = (posTemp.z() - bottom_z_->at(idx_x[2u * i] + idx_y[2u * i] * n_x_)) /
                                                  (top_z_->at(idx_x[2u * i] + idx_y[2u * i] * n_x_) - bottom_z_->at(idx_x[2u * i] + idx_y[2u * i] + n_x_));
                }

                float vertical_factors_min = std::min(std::min(std::min(
                                                                   vertical_factors_columns[0], vertical_factors_columns[1]),
                                                               vertical_factors_columns[2]),
                                                      vertical_factors_columns[3]);

                float vertical_factors_max = std::max(std::max(std::max(
                                                                   vertical_factors_columns[0], vertical_factors_columns[1]),
                                                               vertical_factors_columns[2]),
                                                      vertical_factors_columns[3]);

                // Check if aircraft is out of wind field or not, and act accordingly.
                if (x_inf >= 0u && y_inf >= 0u && vertical_factors_max >= 0u &&
                    x_sup <= (n_x_ - 1u) && y_sup <= (n_y_ - 1u) && vertical_factors_min <= 1u) {
                    // Find indices in z-direction for each of the vertices. If link is not
                    // within the range of one of the columns, set to lowest or highest two.
                    std::size_t idx_z[n_vertices] = { 0u, static_cast<int>(vertical_spacing_factors_->size()) - 1u, 0u, static_cast<int>(vertical_spacing_factors_->size()) - 1u, 0u, static_cast<int>(vertical_spacing_factors_->size()) - 1u, 0u, static_cast<int>(vertical_spacing_factors_->size()) - 1u };

                    for (std::size_t i = 0u; i < n_columns; ++i) {
                        if (vertical_factors_columns[i] < 0u) {
                            // Link z-position below lowest grid point of that column.
                            idx_z[2u * i + 1u] = 1u;
                        }
                        else if (vertical_factors_columns[i] >= 1u) {
                            // Link z-position above highest grid point of that column.
                            idx_z[2u * i] = vertical_spacing_factors_->size() - 2u;
                        }
                        else {
                            // Link z-position between two grid points in that column.
                            for (std::size_t j = 0u; j < vertical_spacing_factors_->size() - 1u; ++j) {
                                if (vertical_spacing_factors_->at(j) <= vertical_factors_columns[i] &&
                                    vertical_spacing_factors_->at(j + 1u) > vertical_factors_columns[i]) {
                                    idx_z[2u * i] = j;
                                    idx_z[2u * i + 1u] = j + 1u;
                                    break;
                                }
                            }
                        }
                    }

                    //Extract the wind velocities corresponding to each vertex.
                    vector<Vector3r> wind_at_vertices(n_vertices, Vector3r::Zero());

                    for (std::size_t i = 0u; i < n_vertices; ++i) {

                        float xv = (float)(idx_x[i] + idx_y[i] * n_x_ + idx_z[i] * n_x_ * n_y_);
                        float uxv = u_->at((size_t)xv);
                        wind_at_vertices.at(i).x() = uxv;
                        float yv = (float)(idx_x[i] + idx_y[i] * n_x_ + idx_z[i] * n_x_ * n_y_);
                        float vyv = v_->at((size_t)yv);
                        wind_at_vertices.at(i).y() = vyv;
                        float zv = (float)(idx_x[i] + idx_y[i] * n_x_ + idx_z[i] * n_x_ * n_y_);
                        float wzv = w_->at((size_t)zv);
                        wind_at_vertices.at(i).z() = wzv;
                    }

                    // Extract the relevant coordinate of every point needed for trilinear
                    // interpolation (first z-direction, then x-direction, then y-direction).
                    constexpr unsigned int n_points_interp_z = 8;
                    constexpr unsigned int n_points_interp_x = 4;
                    constexpr unsigned int n_points_interp_y = 2;
                    double interpolation_points[n_points_interp_x + n_points_interp_y + n_points_interp_z];
                    for (std::size_t i = 0u; i < n_points_interp_x + n_points_interp_y + n_points_interp_z; ++i) {
                        if (i < n_points_interp_z) {
                            interpolation_points[i] = (top_z_->at(idx_x[i] + idx_y[i] * n_x_) - bottom_z_->at(idx_x[i] + idx_y[i] * n_x_)) * vertical_spacing_factors_->at(idx_z[i]) + bottom_z_->at(idx_x[i] + idx_y[i] * n_x_);
                        }
                        else if (i >= n_points_interp_z && i < n_points_interp_x + n_points_interp_z) {
                            interpolation_points[i] = min_x_ + res_x_ * idx_x[2u * (i - n_points_interp_z)];
                        }
                        else {
                            interpolation_points[i] = min_y_ + res_y_ * idx_y[4u * (i - n_points_interp_z - n_points_interp_x)];
                        }
                    }

                    wind_velocity = TrilinearInterpolation(
                        posTemp, &wind_at_vertices.at(0), interpolation_points);
                }
                else { // Vehicle is outside of wind field. Secondary check.
                    dataMutex.unlock();
                    return DEFAULT_WIND;
                }

                // XYZ back to NED
                Vector3r vecTemp(wind_velocity);
                wind_velocity.x() = vecTemp.y();
                wind_velocity.y() = vecTemp.x();
                wind_velocity.z() = (-1) * vecTemp.z();

                PREVIOUS_WIND = wind_velocity;

                dataMutex.unlock();

                return wind_velocity;
            }
            catch (const std::exception& ex) {
                WindLog((Utils::stringf("Error in getLocalWind(...) - %s", ex.what())));
                dataMutex.unlock();
                return DEFAULT_WIND;
            }

        } //function end

    private:
        // Status flags
        std::mutex dataMutex;                   // Stops temporal thread and wind calculator thread from modifying/accessing variables simultaneously.
        bool temporal_thread_alive = false;     // Whether the temporal-updating thread is active
        bool updating = false;                  // Whether the temporal-updating thread is in the act of updating the variables
        bool useSpatialWind = false;            // Whether to calculate a spatial wind or use the default (global) value
        bool wind_module_alive = false;         // Whether the module is actively running.

        // Temporal-updating support
        std::thread temporalThread; // Thread that updates wind every UPDATE_INTERVAL milliseconds using data from WIND_FOLDER_PATH
        common_utils::Timer temporalThreadTimer{}; // Timer to maintain updating frequency
        Vector3r PREVIOUS_WIND = Vector3r::Zero(); // Used if getLocalWind(...) is called while temporal-updating thread is in the act of updating
        std::ifstream fin; // File stream for data input

        // Basic init values
        std::string WIND_FOLDER_PATH = ""; // File path to folder holding wind data
        Vector3r DEFAULT_WIND = Vector3r::Zero(); // Default wind value if vehicle is out of wind field, or global wind value if not using spatial wind
        float UPDATE_INTERVAL = 5.0; // How long a wind data file is used before a new one is loaded (if using temporal wind)
        bool LOG_ENABLED = false; // Whether to log debug information (errors are always logged)

        // In-use variables
        float min_x_;
        float min_y_;
        float max_x_;
        float max_y_;
        int n_x_;
        int n_y_;
        float res_x_;
        float res_y_;

        // In-use data
        std::unique_ptr<vector<float>> vertical_spacing_factors_ = std::make_unique<vector<float>>();
        std::unique_ptr<vector<float>> bottom_z_ = std::make_unique<vector<float>>();
        std::unique_ptr<vector<float>> top_z_ = std::make_unique<vector<float>>();
        std::unique_ptr<vector<float>> u_ = std::make_unique<vector<float>>();
        std::unique_ptr<vector<float>> v_ = std::make_unique<vector<float>>();
        std::unique_ptr<vector<float>> w_ = std::make_unique<vector<float>>();

        // "Interim" variables (for temporal updates)
        float new_min_x_;
        float new_min_y_;
        float new_max_x_;
        float new_max_y_;
        int new_n_x_;
        int new_n_y_;
        float new_res_x_;
        float new_res_y_;

        // "Interim" data (temporal updates)
        std::unique_ptr<vector<float>> new_vertical_spacing_factors_;
        std::unique_ptr<vector<float>> new_bottom_z_;
        std::unique_ptr<vector<float>> new_top_z_;
        std::unique_ptr<vector<float>> new_u_;
        std::unique_ptr<vector<float>> new_v_;
        std::unique_ptr<vector<float>> new_w_;

        string getFileNumber(int x)
        {
            return (std::to_string(x / 100 % 10) + std::to_string(x / 10 % 10) + std::to_string(x % 10));
        }

        // Manages reading new data files and updating data in parallel thread
        void temporalUpdate()
        {
            int wind_file_count = 2; // start at 2, given 1 is read at beginning of program
            string wind_file = "wind_" + getFileNumber(wind_file_count) + ".txt";
            temporal_thread_alive = true;
            bool readAgain = true;
            int total_data_files = -1;

            WindLog("Spinning up temporal read thread...");

            // Allow exiting when destructor is called or if something goes wrong
            while (temporal_thread_alive) {
                if (readAgain) {
                    if (readWindField(WIND_FOLDER_PATH + wind_file)) { // Successfully read into file
                        WindLog(Utils::stringf("Temporal Update : Queued file (%d).", wind_file_count));
                    }
                    else { // Unsuccessfully read file (e.g. file does not exist or failed to parse)
                        // On unsuccessful read, we now know the maximum file number to be the previous wind_file_count
                        total_data_files = wind_file_count - 1;

                        wind_file_count = 1;
                        wind_file = "wind_" + getFileNumber(wind_file_count) + ".txt";
                        WindLog(Utils::stringf("Temporal Update : Unable to read next data file (%d), set total_data to %d. Reading %s", wind_file_count, total_data_files, wind_file.c_str()));
                        if (total_data_files < 2 || !readWindField(WIND_FOLDER_PATH + wind_file)) {
                            WindLog(Utils::stringf("Unable to read next wind file: %s, exiting temporal thread.", wind_file.c_str()), Utils::kLogLevelError);
                            temporal_thread_alive = false;
                            return;
                        }
                        WindLog(Utils::stringf("Temporal Update : Queued file (%d).", wind_file_count));
                    }

                    wind_file_count += 1;
                    if (total_data_files > 0 && wind_file_count > total_data_files) {
                        // We've gone past maximum, reset to 1
                        WindLog("Temporal Update : Returning to first data file.");
                        wind_file_count = 1;
                    }
                    wind_file = "wind_" + getFileNumber(wind_file_count) + ".txt";
                    WindLog("Temporal Update : Next wind path: " + WIND_FOLDER_PATH + wind_file);
                    readAgain = false;
                }

                // Delay until UPDATE_INTERVAL time has passed. Then update old data with new data held in "interim" variables.
                temporalThreadTimer.stop();
                int delayMS = (int)(UPDATE_INTERVAL * 1000.0) - (int)temporalThreadTimer.milliseconds();
                std::this_thread::sleep_for(std::chrono::milliseconds(delayMS));
                if (temporal_thread_alive) { // Ready to update!
                    updateWindField();
                    readAgain = true;
                    temporalThreadTimer.start();
                }
            }
        }

        // Updates data being used to calculate wind
        void updateWindField()
        {

            updating = false;
            // Spin until can reserve data
            while(!updating && temporal_thread_alive) {
                if(dataMutex.try_lock()) {
                    updating = true;
                }
                else {
                    std::this_thread::sleep_for(std::chrono::microseconds(5));
                }
            }

            // Exit quickly if temporal thread shut down
            if (!temporal_thread_alive) return;

            // Swap primitive data
            min_x_ = new_min_x_;
            min_x_ = new_min_x_;
            min_y_ = new_min_y_;
            max_x_ = new_max_x_;
            max_y_ = new_max_y_;
            n_x_ = new_n_x_;
            n_y_ = new_n_y_;
            res_x_ = new_res_x_;
            res_y_ = new_res_y_;

            // Swap vector-based data
            std::swap(vertical_spacing_factors_, new_vertical_spacing_factors_);
            std::swap(bottom_z_, new_bottom_z_);
            std::swap(top_z_, new_top_z_);
            std::swap(u_, new_u_);
            std::swap(v_, new_v_);
            std::swap(w_, new_w_);

            // Release reservation on data
            dataMutex.unlock();

            // Free up original copy of fresh data
            new_vertical_spacing_factors_.reset(nullptr);
            new_bottom_z_.reset(nullptr);
            new_top_z_.reset(nullptr);
            new_u_.reset(nullptr);
            new_v_.reset(nullptr);
            new_w_.reset(nullptr);
        }

        // Gather wind field data from provided file.
        bool readWindField(std::string wind_field_path)
        {

            fin.open(wind_field_path);
            if (!fin.is_open()) {
                WindLog("Unable to open wind data file!", Utils::kLogLevelError);
                return false;
            }

            new_vertical_spacing_factors_ = std::make_unique<vector<float>>();
            new_bottom_z_ = std::make_unique<vector<float>>();
            new_top_z_ = std::make_unique<vector<float>>();
            new_u_ = std::make_unique<vector<float>>();
            new_v_ = std::make_unique<vector<float>>();
            new_w_ = std::make_unique<vector<float>>();

            std::string data_name;
            float data;

            // Parse file and gather information into "interim" variables
            while (wind_module_alive && fin >> data_name) {
                if (data_name == "min_x:") {
                    fin >> new_min_x_;
                }
                else if (data_name == "min_y:") {
                    fin >> new_min_y_;
                }
                else if (data_name == "n_x:") {
                    fin >> new_n_x_;
                }
                else if (data_name == "n_y:") {
                    fin >> new_n_y_;
                }
                else if (data_name == "res_x:") {
                    fin >> new_res_x_;
                }
                else if (data_name == "res_y:") {
                    fin >> new_res_y_;
                }
                else if (data_name == "vertical_spacing_factors:") {
                    while (fin >> data) {
                        new_vertical_spacing_factors_->push_back(data);
                        if (fin.peek() == '\n') break;
                    }
                }
                else if (data_name == "bottom_z:") {
                    while (fin >> data) {
                        new_bottom_z_->push_back(data);
                        if (fin.peek() == '\n') break;
                    }
                }
                else if (data_name == "top_z:") {
                    while (fin >> data) {
                        new_top_z_->push_back(data);
                        if (fin.peek() == '\n') break;
                    }
                }
                else if (data_name == "u:") {
                    while (fin >> data) {
                        new_u_->push_back(data);
                        if (fin.peek() == '\n') break;
                    }
                }
                else if (data_name == "v:") {
                    while (fin >> data) {
                        new_v_->push_back(data);
                        if (fin.peek() == '\n') break;
                    }
                }
                else if (data_name == "w:") {
                    while (fin >> data) {
                        new_w_->push_back(data);
                        if (fin.peek() == '\n') break;
                    }
                }
                else {
                    std::string restOfLine;
                    getline(fin, restOfLine);
                    WindLog("Unknown content in wind data file: " + restOfLine, Utils::kLogLevelWarn);
                }
            }
            new_max_x_ = new_min_x_ + (float)(new_n_x_ * new_res_x_);
            new_max_y_ = new_min_y_ + (float)(new_n_y_ * new_res_y_);
            fin.close();

            return true;
        }

        Vector3r LinearInterpolation(
            double position, Vector3r* values, double* points)
        {
            Vector3r value = values[0] + (values[1] - values[0]) /
                                             (points[1] - points[0]) * (position - points[0]);
            return value;
        }

        Vector3r BilinearInterpolation(
            double* position, Vector3r* values, double* points)
        {
            Vector3r intermediate_values[2] = { LinearInterpolation(
                                                    position[0], &(values[0]), &(points[0])),
                                                LinearInterpolation(
                                                    position[0], &(values[2]), &(points[2])) };
            Vector3r value = LinearInterpolation(
                position[1], intermediate_values, &(points[4]));
            return value;
        }

        Vector3r TrilinearInterpolation(
            const Vector3r& vehicle_position, Vector3r* values, double* points)
        {
            double position[3] = { vehicle_position.x(), vehicle_position.y(), vehicle_position.z() };
            Vector3r intermediate_values[4] = { LinearInterpolation(
                                                    position[2], &(values[0]), &(points[0])),
                                                LinearInterpolation(
                                                    position[2], &(values[2]), &(points[2])),
                                                LinearInterpolation(
                                                    position[2], &(values[4]), &(points[4])),
                                                LinearInterpolation(
                                                    position[2], &(values[6]), &(points[6])) };
            Vector3r value = BilinearInterpolation(
                &(position[0]), intermediate_values, &(points[8]));
            return value;
        }

        void WindLog(string info, int level = Utils::kLogLevelInfo) {
            string header = "Wind Module: ";
            if (LOG_ENABLED) airlib::Utils::log(Utils::stringf((header.append(info)).c_str()), level);
        }
    };

} //namespace
} //namespace
#endif