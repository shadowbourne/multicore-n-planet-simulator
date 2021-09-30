// Translate this file with
//
// g++ -O3 assignment-code.cpp -o assignment-code
//
// Run it with
//
// ./demo-code
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2020 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>



#include <cmath>


double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

/**
 * 1e-2/NumberOfBodies (coefficient to determine maximum merging distance)
 */
double C;

/**
 * One coordinate entry along the x direction per molecule/particle.
 */
double* x0;

/**
 * One coordinate entry along the y direction per molecule/particle.
 */
double* x1;

/**
 * One coordinate entry along the z direction per molecule/particle.
 */
double* x2;

/**
 * One velocity entry along the x direction per molecule/particle.
 */
double* v0;

/**
 * One velocity entry along the y direction per molecule/particle.
 */
double* v1;

/**
 * One velocity entry along the z direction per molecule/particle.
 */
double* v2;
/**
 * One mass entry per molecule/particle.
 */
double* mass;

/**
 * One force entry along the x direction per molecule/particle.
 */
double* force0;

/**
 * One force entry along the y direction per molecule/particle.
 */
double* force1;

/**
 * One force entry along the z direction per molecule/particle.
 */
double* force2;

/**
 * Global time step size used.
 */
double timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double maxV;
/**
 * Maximum velocity of all particles.
 */
double maxMass;

/**
 * Minimum distance between two elements.
 */
double minDx;

// /**
//  * Temp variable used for calculating distances within threads (initialised as private).
//  */
// double Dx;

// /**
//  * Temporary indicator/control variable for breaking while loop.
//  */
// bool noMergeFound = true;

/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-4) / 7;
  // WARNING: if you change C you must also change CSqrd within the updateBody function
  C = 1e-2 / NumberOfBodies;

  x0 = new double[NumberOfBodies];
  x1 = new double[NumberOfBodies];
  x2 = new double[NumberOfBodies];
  v0 = new double[NumberOfBodies];
  v1 = new double[NumberOfBodies];
  v2 = new double[NumberOfBodies];
  mass = new double [NumberOfBodies];

  force0 = new double[NumberOfBodies];
  force1 = new double[NumberOfBodies];
  force2 = new double[NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    // x=[i]new double[3];
    // v=[i]new double[3];

    x0[i] = std::stof(argv[readArgument]); readArgument++;
    x1[i] = std::stof(argv[readArgument]); readArgument++;
    x2[i] = std::stof(argv[readArgument]); readArgument++;

    v0[i] = std::stof(argv[readArgument]); readArgument++;
    v1[i] = std::stof(argv[readArgument]); readArgument++;
    v2[i] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;
    maxMass = std::max(maxMass, mass[i]);

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
  
  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;
  }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}





/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x0[i]
        << " "
        << x1[i]
        << " "
        << x2[i]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}



/**
 * This is the main operation you should change in the assignment. You might
 * want to add a few more variables or helper functions, but this is where the
 * magic happens.
 */
void updateBody() {

  #pragma ivdep
  for (int b=0; b < NumberOfBodies; b++) {
    force0[b] = 0.0;
    force1[b] = 0.0;
    force2[b] = 0.0;
  }

  // Calculate forces
  for (int b=0; b < NumberOfBodies-1; b++) {
    #pragma ivdep
    #pragma omp simd reduction(+:force0[b],force1[b],force2[b]) 
    for (int i=b+1; i<NumberOfBodies; i++) {
      const double TDx0 = x0[i]-x0[b];
      const double TDx1 = x1[i]-x1[b];
      const double TDx2 = x2[i]-x2[b];
      const double Dx = sqrt(
        TDx0 * TDx0 +
        TDx1 * TDx1 +
        TDx2 * TDx2
      );
      const double Dx3 = Dx * Dx * Dx;
      const double massDivDist = mass[i] * mass[b] / Dx3;
      // x,y,z forces acting on particle b and therefore i by Newton's Third Law of Motion
      const double Tforce0 = TDx0 * massDivDist;
      force0[b] += Tforce0 ;
      force0[i] -= Tforce0 ;

      const double Tforce1 = TDx1 * massDivDist;
      force1[b] += Tforce1 ;
      force1[i] -= Tforce1 ;

      const double Tforce2 = TDx2 * massDivDist;
      force2[b] += Tforce2 ;
      force2[i] -= Tforce2 ;
    }
  }

  // Update positions with old velocities and then velocities with new forces
  #pragma ivdep
  for (int b=0; b < NumberOfBodies; b++) {
    x0[b] += timeStepSize * v0[b];
    x1[b] += timeStepSize * v1[b];
    x2[b] += timeStepSize * v2[b];

    // const double tDivMass = timeStepSize / mass[b];

    v0[b] += force0[b] * timeStepSize / mass[b];
    v1[b] += force1[b] * timeStepSize / mass[b];
    v2[b] += force2[b] * timeStepSize / mass[b];

  }

  double minDxSqrd;
  // WARNING: If C is changed here it must also be changed in the setUp function
  const double CSqrd = 1e-4 / (NumberOfBodies * NumberOfBodies);
  bool noMergeFound = false;
  while (!noMergeFound) {

    // Section to see if any particles must be merged, by calculating minDx in a while loop

    noMergeFound = true;
    minDxSqrd  = std::numeric_limits<double>::max();
    for (int i=0; i < NumberOfBodies; i++) {
      // if (noMergeFound) {
        #pragma ivdep
        #pragma omp simd reduction(min:minDxSqrd)
        for (int j=i+1; j < NumberOfBodies; j++) {
            const double TDx0 = x0[i]-x0[j];
            const double TDx1 = x1[i]-x1[j];
            const double TDx2 = x2[i]-x2[j];
            const double DxSqrd = (
            TDx0 * TDx0 +
            TDx1 * TDx1 +
            TDx2 * TDx2 
            );
            minDxSqrd = std::min( minDxSqrd, DxSqrd);
          
            // if ( DxSqrd <= CSqrd * (mass[i] + mass[j]) * (mass[i] + mass[j])) {
            //   noMergeFound = false; // This also stops any more inner loops being entered.
            // }
        }
      // }
    }
    // Check to see if objects may need fusing
    if ( minDxSqrd <= CSqrd * (maxMass + maxMass) * (maxMass + maxMass)) {
              noMergeFound = false; 
    }
    if (!noMergeFound) {
      // Section to check if objects definitely need fusing if we discovered above that some may need fusing
      noMergeFound = true;
      minDxSqrd  = std::numeric_limits<double>::max();
      for (int i=0; i < NumberOfBodies; i++) {
        if (noMergeFound) {
          #pragma ivdep
          #pragma omp simd reduction(&&:noMergeFound)
          for (int j=i+1; j < NumberOfBodies; j++) {
              const double TDx0 = x0[i]-x0[j];
              const double TDx1 = x1[i]-x1[j];
              const double TDx2 = x2[i]-x2[j];
              const double DxSqrd = (
              TDx0 * TDx0 +
              TDx1 * TDx1 +
              TDx2 * TDx2 
              );
              // minDxSqrd = std::min( minDxSqrd, DxSqrd);
            
              if ( DxSqrd <= CSqrd * (mass[i] + mass[j]) * (mass[i] + mass[j])) {
                noMergeFound = false; // This also stops any more inner loops being entered.
              }
          }
        }
      }
    }
    
    // Section to fuse objects if we discovered above that some need fusing
    if (!noMergeFound) {
      noMergeFound = true;
      for (int i=0; i < NumberOfBodies; i++) {
        for (int j=i+1; j < NumberOfBodies; j++) {
          const double TDx0 = x0[i]-x0[j];
          const double TDx1 = x1[i]-x1[j];
          const double TDx2 = x2[i]-x2[j];
          // If statement to merge two particles if they are close enough after the end of this timestep
          if ( sqrt(
              TDx0 * TDx0 +
              TDx1 * TDx1 +
              TDx2 * TDx2 
            )
            <= 
            C * (mass[i] + mass[j]) ) {
                      noMergeFound = false;
                      const double massi = mass[i];
                      const double massj = mass[j];

                      // Setting mass of new particle to the sum of both masses
                      mass[i] = massi + massj;

                      // Updating maxMass
                      maxMass = std::max(maxMass, mass[i]);

                      // Would normally do this at the end, however this allows saving time
                      // on doing NumberOfBodies - 1 for each index
                      NumberOfBodies--; 

                      // Applying weighted mean to velocity of the new merged particle
                      v0[i] = (v0[i] * massi + v0[j] * massj) / (massi + massj);
                      v1[i] = (v1[i] * massi + v1[j] * massj) / (massi + massj);
                      v2[i] = (v2[i] * massi + v2[j] * massj) / (massi + massj);

                      // Swapping the particle at the end of the memory array into the now disappeared j'th particle's memory location
                      v0[j] = v0[NumberOfBodies];
                      v1[j] = v1[NumberOfBodies];
                      v2[j] = v2[NumberOfBodies];

                      // Averaging the locations of the two particles
                      x0[i] = (x0[i] * massi + x0[j] * massj) / (massi + massj);
                      x1[i] = (x1[i] * massi + x1[j] * massj) / (massi + massj);
                      x2[i] = (x2[i] * massi + x2[j] * massj) / (massi + massj);


                      // Swapping the particle at the end of the memory array into the now disappeared j'th particle's memory location
                      mass[j] = mass[NumberOfBodies];

                      // Swapping the particle at the end of the memory array into the now disappeared j'th particle's memory location
                      x0[j] = x0[NumberOfBodies];
                      x1[j] = x1[NumberOfBodies];
                      x2[j] = x2[NumberOfBodies];

                      // Move loop on j back one iteration so that it checks the new j
                      j--;
            }
        }
      }
    }
  }
  minDx = std::sqrt(minDxSqrd);
  // Calculating max velocity 
  double maxVSqrd = 0.0;
  #pragma ivdep
  for (int b=0; b < NumberOfBodies; b++) {
    maxVSqrd = std::max(maxVSqrd, ( v0[b]*v0[b] + v1[b]*v1[b] + v2[b]*v2[b] ));
  }
  maxV = std::sqrt(maxVSqrd);
  t += timeStepSize;

}


/**
 * Main routine.
 *
 * No major changes in assignment. You can add a few initialisation
 * or stuff if you feel the need to do so. But keep in mind that you
 * may not alter what the program plots to the terminal.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << "  dt              time step size (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t One spiralling around the other one" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t Three body setup from first lecture" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0 \t Five body setup" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta; 
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      printParaviewSnapshot();
      std::cout << "plot next snapshot"
    		    << ",\t time step=" << timeStepCounter
    		    << ",\t t="         << t
				<< ",\t dt="        << timeStepSize
				<< ",\t v_max="     << maxV
				<< ",\t dx_min="    << minDx
				<< std::endl;

      tPlot += tPlotDelta;
    }
  }

  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: " << x0[0] << ", " << x1[0] << ", " << x2[0] << std::endl;

  closeParaviewVideoFile();

  return 0;
}
