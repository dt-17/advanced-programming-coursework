//
//  solution_visualiser.cpp
//  advanced-programming
//
//  Created by Daniel Tompkins on 10/03/2025.
//

#include "solution_visualiser.hpp"    // include header for the visualiser class
#include "gnuplot-iostream.h"         // include gnuplot-iostream to interface with GNU Plot
#include <fstream>                    // include fstream for file I/O operations
#include <sstream>                    // include sstream for string stream parsing
#include <iostream>                   // include iostream for console I/O
#include <string>                     // include string for string handling
#include <cstdio>    // for snprintf and remove
#include <cstdlib>   // for system

// constructor for the Visualiser class
// initialises the filename, grid dimensions and allocates memory for the 3D solution vector
// the vector 'u' is used to store the solution for each time step, organised as u[time][x][y]
Visualiser::Visualiser(const std::string &filename, int nx, int ny, int nt)
    : filename(filename), nx(nx), ny(ny), nt(nt)
{
    // allocate space for the solution: nt time steps of nx by ny grids
    u.resize(nt, std::vector<std::vector<double>>(nx, std::vector<double>(ny, 0.0)));
}

// destructor for the Visualiser class
// there is no dynamic memory to free here as the vector cleans up automatically
Visualiser::~Visualiser() {
    // nothing specific to clean up
}

// load the solution data from the file into the 3D vector 'u'
// the file is expected to contain the solution data organised in blocks corresponding to each time step
bool Visualiser::loadSolution() {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    std::string line;
    
    // read and discard first four lines (the metadata lines)
    for (int k = 0; k < 4; ++k) {
        std::getline(infile, line);
    }

    
    // loop over each time step
    for (int t = 0; t < nt; ++t) {
        // loop over each spatial row in the grid
        for (int i = 0; i < nx; ++i) {
            // read the next non-empty line from the file
            do {
                if (!std::getline(infile, line)) {
                    std::cerr << "Unexpected end of file." << std::endl;
                    return false;
                }
            } while (line.empty());
            
            // use a string stream to extract values from the line
            std::istringstream iss(line);
            for (int j = 0; j < ny; ++j) {
                double val;
                iss >> val;
                u[t][i][j] = val;
            }
        }
        // optionally, skip any extra blank lines that separate time steps
        std::getline(infile, line);
    }
    infile.close();
    return true;
}

// display a 3D plot of the final time step using GNU Plot
// opens a GNU Plot window and plots the data stored in the final grid (u[nt-1])
void Visualiser::showFinalTimeStep() {
    if (nt <= 0) return;

    // retrieve the final time step data
    const auto &finalData = u[nt - 1];

    // define the domain in x and y (e.g. 0 to 10)
    double xMin = 0.0, xMax = 10.0;
    double yMin = 0.0, yMax = 10.0;
    double dx = (xMax - xMin) / (nx - 1);
    double dy = (yMax - yMin) / (ny - 1);

    // create a GNU Plot object
    Gnuplot gp;

    // set up the plot title, labels and other plot properties
    gp << "set title '2D Wave Equation Solution at Final Time'\n";
    gp << "set xlabel 'x'\n";
    gp << "set ylabel 'y'\n";
    gp << "set zlabel 'Amplitude'\n";
    gp << "set hidden3d\n";
    gp << "splot '-' using 1:2:3 with pm3d notitle\n";

    // output the x, y and z (amplitude) data for each grid point
    for (int i = 0; i < nx; i++) {
        double x = xMin + i * dx;
        for (int j = 0; j < ny; j++) {
            double y = yMin + j * dy;
            double z = finalData[i][j];
            gp << x << " " << y << " " << z << "\n";
        }
        gp << "\n"; // separate each row with a blank line
    }
    gp << "e\n";
    gp << std::flush;

    // inform the user and wait for them to press enter before closing the plot
    std::cout << "Press enter to exit visualization..." << std::endl;
    std::cin.get();
}

// save a 3D plot snapshot of the final time step as a PNG file
// this function writes the final time step plot to the specified output path
void Visualiser::saveFinalTimeStep(const std::string &outputPath) {
    if (nt <= 0) return;

    // retrieve the final time step data
    const auto &finalData = u[nt - 1];

    // define the domain in x and y (e.g. 0 to 10)
    double xMin = 0.0, xMax = 10.0;
    double yMin = 0.0, yMax = 10.0;
    double dx = (xMax - xMin) / (nx - 1);
    double dy = (yMax - yMin) / (ny - 1);

    // create a GNU Plot object
    Gnuplot gp;
    
    // set up the PNG output: terminal settings and output file path
    gp << "set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n";
    gp << "set output '" << outputPath << "'\n";

    // set up the plot title, labels and other plot properties
    gp << "set title '2D Wave Equation Solution at Final Time'\n";
    gp << "set xlabel 'x'\n";
    gp << "set ylabel 'y'\n";
    gp << "set zlabel 'Amplitude'\n";
    gp << "set hidden3d\n";
    gp << "splot '-' using 1:2:3 with pm3d notitle\n";

    // output the x, y and z data for each grid point
    for (int i = 0; i < nx; i++) {
        double x = xMin + i * dx;
        for (int j = 0; j < ny; j++) {
            double y = yMin + j * dy;
            double z = finalData[i][j];
            gp << x << " " << y << " " << z << "\n";
        }
        gp << "\n"; // separate each row with a blank line
    }
    gp << "e\n";
    gp << std::flush;

    // inform the user that the snapshot has been saved
    std::cout << "Snapshot saved to " << outputPath << std::endl;
}

// display an animation of the solution evolution over time using GNU Plot
// the function cycles through each time step, plotting the solution on the defined spatial domain
// parameters:
//   xMin, xMax - domain limits for the x-axis
//   yMin, yMax - domain limits for the y-axis
//   pauseDuration - pause between frames in seconds
// Modified function: an extra parameter "videoOutputPath" lets you specify where to save the MP4.
void Visualiser::showAnimation(double xMin, double xMax, double yMin, double yMax, double pauseDuration, const std::string &videoOutputPath) {
    if (nt <= 0) return;

    // Calculate grid spacing based on the provided domain.
    double dx_grid = (xMax - xMin) / (nx - 1);
    double dy_grid = (yMax - yMin) / (ny - 1);

    // Create a Gnuplot instance.
    Gnuplot gp;
    
    // Set common plot properties.
    gp << "set xlabel 'x'\n";
    gp << "set ylabel 'y'\n";
    gp << "set zlabel 'Amplitude'\n";
    gp << "set hidden3d\n";

    // Loop over each time step and save each frame as a PNG file.
    char frameFilename[100];
    for (int t = 0; t < nt; t++) {
        // Create a file name like "frame_0000.png", "frame_0001.png", etc.
        snprintf(frameFilename, sizeof(frameFilename), "frame_%04d.png", t);
        std::string frameFile(frameFilename);
        
        // Set terminal and output file for this frame.
        gp << "set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n";
        gp << "set output '" << frameFile << "'\n";
        gp << "set title '2D Wave Equation Solution Animation (t = " << t << ")'\n";
        
        // Plot the current time step.
        gp << "splot '-' using 1:2:3 with pm3d notitle\n";
        for (int i = 0; i < nx; i++) {
            double x = xMin + i * dx_grid;
            for (int j = 0; j < ny; j++) {
                double y = yMin + j * dy_grid;
                double z = u[t][i][j];
                gp << x << " " << y << " " << z << "\n";
            }
            gp << "\n";  // Separate rows.
        }
        gp << "e\n" << std::flush;
    }

    // Determine the frame rate (assumes pauseDuration is in seconds).
    int frameRate = static_cast<int>(1.0 / pauseDuration);
    
    // Build the ffmpeg command to combine frames into an MP4 video.
    std::stringstream ffmpegCmd;
    ffmpegCmd << "ffmpeg -y -framerate " << frameRate
              << " -i frame_%04d.png -c:v libx264 -pix_fmt yuv420p " << videoOutputPath;
    int ret = system(ffmpegCmd.str().c_str());
    if (ret != 0) {
        std::cerr << "ffmpeg command failed with return code: " << ret << std::endl;
    } else {
        std::cout << "Animation saved to " << videoOutputPath << std::endl;
    }

    // Remove temporary frame images.
    for (int t = 0; t < nt; t++) {
        snprintf(frameFilename, sizeof(frameFilename), "frame_%04d.png", t);
        std::remove(frameFilename);
    }
}

