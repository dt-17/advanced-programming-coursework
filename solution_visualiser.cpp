//
//  solution_visualiser.cpp
//  advanced-programming
//
//  Created by Daniel Tompkins on 10/03/2025.
//

#include "solution_visualiser.hpp"
#include "gnuplot-iostream.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

Visualizer::Visualizer(const std::string &filename, int nx, int ny, int nt)
    : filename(filename), nx(nx), ny(ny), nt(nt)
{
    // Allocate space for the solution: nt time steps of nx x ny grids.
    u.resize(nt, std::vector<std::vector<double>>(nx, std::vector<double>(ny, 0.0)));
}

Visualizer::~Visualizer() {
    // Nothing specific to clean up.
}

bool Visualizer::loadSolution() {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return false;
    }

    std::string line;
    // Loop over time steps.
    for (int t = 0; t < nt; ++t) {
        // Loop over spatial rows.
        for (int i = 0; i < nx; ++i) {
            // Get the next non-empty line.
            do {
                if (!std::getline(infile, line)) {
                    std::cerr << "Unexpected end of file." << std::endl;
                    return false;
                }
            } while (line.empty());
            
            std::istringstream iss(line);
            for (int j = 0; j < ny; ++j) {
                double val;
                iss >> val;
                u[t][i][j] = val;
            }
        }
        // Optionally, skip any extra blank lines separating time steps.
        std::getline(infile, line);
    }
    infile.close();
    return true;
}

void Visualizer::showFinalTimeStep() {
    if (nt <= 0) return;

    // Get the final time step data.
    const auto &finalData = u[nt - 1];

    // Domain in x and y, e.g. 0 to 10
    double xMin = 0.0, xMax = 10.0;
    double yMin = 0.0, yMax = 10.0;
    double dx = (xMax - xMin) / (nx - 1);
    double dy = (yMax - yMin) / (ny - 1);

    Gnuplot gp;

    gp << "set title '2D Wave Equation Solution at Final Time'\n";
    gp << "set xlabel 'x'\n";
    gp << "set ylabel 'y'\n";
    gp << "set zlabel 'Amplitude'\n";
    gp << "set hidden3d\n";
    gp << "splot '-' using 1:2:3 with pm3d notitle\n";

    // Output x, y, z for each grid point
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

    std::cout << "Press enter to exit visualization..." << std::endl;
    std::cin.get();
}

void Visualizer::saveFinalTimeStep(const std::string &outputPath) {
    if (nt <= 0) return;

    // Get the final time step data.
    const auto &finalData = u[nt - 1];

    // Domain in x and y, e.g. 0 to 10
    double xMin = 0.0, xMax = 10.0;
    double yMin = 0.0, yMax = 10.0;
    double dx = (xMax - xMin) / (nx - 1);
    double dy = (yMax - yMin) / (ny - 1);

    Gnuplot gp;
    
    // Set up PNG output. You can adjust size, font, etc.
    gp << "set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n";
    gp << "set output '" << outputPath << "'\n";

    gp << "set title '2D Wave Equation Solution at Final Time'\n";
    gp << "set xlabel 'x'\n";
    gp << "set ylabel 'y'\n";
    gp << "set zlabel 'Amplitude'\n";
    gp << "set hidden3d\n";
    gp << "splot '-' using 1:2:3 with pm3d notitle\n";

    // Output x, y, z for each grid point
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

    std::cout << "Snapshot saved to " << outputPath << std::endl;
}



void Visualizer::showAnimation(double xMin, double xMax, double yMin, double yMax, double pauseDuration) {
    if (nt <= 0) return;

    // Calculate grid spacing based on the domain.
    double dx = (xMax - xMin) / (nx - 1);
    double dy = (yMax - yMin) / (ny - 1);

    // Create a Gnuplot instance.
    Gnuplot gp;
    
    // Set up basic plot properties.
    gp << "set title '2D Wave Equation Solution Animation'\n";
    gp << "set xlabel 'x'\n";
    gp << "set ylabel 'y'\n";
    gp << "set zlabel 'Amplitude'\n";
    gp << "set hidden3d\n";

    // Loop over each time step.
    for (int t = 0; t < nt; t++) {
        gp << "splot '-' using 1:2:3 with pm3d notitle\n";

        // Send the (x, y, z) data for the current time step.
        for (int i = 0; i < nx; i++) {
            double x = xMin + i * dx;
            for (int j = 0; j < ny; j++) {
                double y = yMin + j * dy;
                double z = u[t][i][j];
                gp << x << " " << y << " " << z << "\n";
            }
            gp << "\n";  // Blank line to separate rows.
        }
        gp << "e\n" << std::flush;

        // Pause for the specified duration between frames.
        gp << "pause " << pauseDuration << "\n";
    }
    
    // Final pause to let the user interact before closing.
    gp << "pause mouse\n";
}


