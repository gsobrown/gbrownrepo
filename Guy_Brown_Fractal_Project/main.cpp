#include "EasyBMP.h"
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <random>


//Function to check whether each transformation is valid
double validTrans(double numberOfLines, std::vector<std::vector<double>> w_tot){
    int valid = 0;
    for (int k = 0; k < numberOfLines; k++) {
        double detA = w_tot[k][0] * w_tot[k][3] - w_tot[k][1] * w_tot[k][2];
        if (pow(w_tot[k][0], 2) + pow(w_tot[k][2], 2) >= 1 || pow(w_tot[k][1], 2) + pow(w_tot[k][3], 2) >= 1 ||
            pow(w_tot[k][0], 2) + pow(w_tot[k][2], 2) + pow(w_tot[k][1], 2) + pow(w_tot[k][3], 2) - pow(detA, 2) >= 1) {
            std::cout << "Transformation " << k+1 << " is invalid" << std::endl;
        } else {
            std::cout << "Transformation " << k+1 << " is valid" << std::endl;
            valid++;
        }
    }
    return valid;
}

//Functions used to scale image to fit in x and y directions
double getMin(std::vector<double>&v){
    return *std::min_element(v.begin(), v.end()); //dereferences iterator
}

double getMax(std::vector<double>&v){
    return *std::max_element(v.begin(), v.end());
}

//Function to transform the pixels for each image...
double Trans(std::vector<double>&v, int translation) {
    return *std::transform(v.begin(), v.end(), v.begin(), bind1st(std::plus<double>(), translation - getMin(v)));
}

int main(){
    int sectionInput;
    std::cout << "Enter the number for the section you would like to exlpore..." << std::endl;
    std::cout << "1 - The Chaos Game" << std::endl;
    std::cout << "2 - Generalised Chaos Game" << std::endl;
    std::cout << "3 - The Collage Theorem" << std::endl;
    std::cin >> sectionInput;

    if (sectionInput == 1){
    std::cout << "Chaos Game" << std::endl;
        //Create a 2D array with MxN pixels on the heap
        int M = 1920;
        int N = 1200;
        double **data = new double*[M];

        for(int i = 0; i < M; ++i){

            data[i] = new double[N];

        }

        //Initialise pixel array to zero (i.e. empty image)
        for(int i = 0; i < M; ++i){
            for(int j = 0; j < N; ++j){
                data[i][j] = 0;
            }
        }


        // declare some variables to be determined by the user...
        std::string shapeName;
        double stepSize;
        double desired;
        int iterations;

        // Some sanitised user inputs...
        std::cout << "Which shape would you like to produce? Options: 'triangle', 'rectangle', 'pentagon', 'hexagon'" << std::endl;
        std::cin >> shapeName;
        if(shapeName != "triangle" && shapeName != "rectangle" && shapeName != "pentagon" && shapeName != "hexagon"){
            std::cout << "You have not selected a valid shape" << std::endl;
            return 0;
        }

        std::cout << "Select the step size for each iteration. Step size must be between 0.5 and 0.2 inclusive" << std::endl;
        std::cin >> stepSize;
        if (stepSize > 0.5 || stepSize < 0.2){
            std::cout << "Invalid step size" << std::endl;
            return 0;
        }

        std::cout << "Determine the probability (0-1) of selecting one of the particular vertices"  << std::endl;
        std::cin >> desired;
        if (desired < 0 || desired > 1){
            std::cout << "Probability outside of range" << std::endl;
            return 0;
        }

        std::cout << "Determine the number of iterations to perform: " << std::endl;
        std::cin >> iterations;

        if (shapeName == "triangle") {
            // vertices of triangle
            double vertex_a[2] = {0, 0};
            double vertex_b[2] = {960, 1200};
            double vertex_c[2] = {1920, 0};
            double numSides = 3;
            // choose starting vertex
            double startPoint[2];
            startPoint[0] = vertex_a[0];
            startPoint[1] = vertex_a[1];
            for (int k = 0; k < iterations; k++) {
                double random = 1.0*rand()/ RAND_MAX;
                if (random <= desired) {
                    double xDist = vertex_a[0] + startPoint[0]; //Calculate total horizontal coordinate
                    double yDist = vertex_a[1] + startPoint[1]; //Calculate total vertical coordinate
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)}; //new point is total coordinate scaled by step size
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;  //Fill pixels in data array
                    startPoint[0] = newPoint[0];    //restart iteration from new point
                    startPoint[1] = newPoint[1];
                } else if (desired < random && random <= desired + (1-desired)/(numSides-1)) {
                    double xDist = vertex_b[0] + startPoint[0];
                    double yDist = vertex_b[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                } else if (desired + (1-desired)/(numSides-1) < random) {
                    double xDist = vertex_c[0] + startPoint[0];
                    double yDist = vertex_c[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                }
            }
        } else if (shapeName == "rectangle"){
            // vertices of rectangle
            double vertex_a[2] = {0, 0};
            double vertex_b[2] = {0, 1200};
            double vertex_c[2] = {1920, 1200};
            double vertex_d[2] = {1920,0};
            double numSides = 4;
            // choose starting vertex
            double startPoint[2];
            startPoint[0] = vertex_a[0];
            startPoint[1] = vertex_a[1];
            for (int k = 0; k < iterations; k++) {
                double random = 1.0*rand()/ RAND_MAX;
                if (random <= desired) {
                    double xDist = vertex_a[0] + startPoint[0];
                    double yDist = vertex_a[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                } else if (desired < random && random <= desired + (1-desired)/(numSides-1)) {
                    double xDist = vertex_b[0] + startPoint[0];
                    double yDist = vertex_b[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                } else if (desired + (1-desired)/(numSides-1) < random && random <= desired + 2*(1-desired)/(numSides-1)) {
                    double xDist = vertex_c[0] + startPoint[0];
                    double yDist = vertex_c[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                } else if (desired + 2*(1-desired)/(numSides-1) < random) {
                    double xDist = vertex_d[0] + startPoint[0];
                    double yDist = vertex_d[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                }
            }
        } else if (shapeName == "pentagon") {
            // vertices of pentagon
            double vertex_a[2] = {1045, 855};
            double vertex_b[2] = {864, 986};
            double vertex_c[2] = {932, 1198};
            double vertex_d[2] = {1157, 1198};
            double vertex_e[2] = {1225, 986};
            double numSides = 5;
            // choose starting vertex
            double startPoint[2];
            startPoint[0] = vertex_a[0];
            startPoint[1] = vertex_a[1];
            for (int k = 0; k < iterations; k++) {
                double random = 1.0*rand()/ RAND_MAX;
                if (random <= desired) { //the if / else if statements take the probability weighting of vertices into account
                    double xDist = vertex_a[0] + startPoint[0];
                    double yDist = vertex_a[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                } else if (desired < random && random <= desired + (1-desired)/(numSides-1)) {
                    double xDist = vertex_b[0] + startPoint[0];
                    double yDist = vertex_b[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                } else if (desired + (1-desired)/(numSides-1) < random && random <= desired + 2*(1-desired)/(numSides-1)) {
                    double xDist = vertex_c[0] + startPoint[0];
                    double yDist = vertex_c[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                } else if (desired + 2*(1-desired)/(numSides-1) < random && random <= desired + 3*(1-desired)/(numSides-1)) {
                    double xDist = vertex_d[0] + startPoint[0];
                    double yDist = vertex_d[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                } else if (desired + 3*(1-desired)/(numSides-1) < random) {
                    double xDist = vertex_e[0] + startPoint[0];
                    double yDist = vertex_e[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                }
            }
        } else if (shapeName == "hexagon") {
            // vertices of hexagon
            double vertex_a[2] = {1080,833};
            double vertex_b[2] = {900,833};
            double vertex_c[2] = {810,990};
            double vertex_d[2] = {900,1146};
            double vertex_e[2] = {1080,1146};
            double vertex_f[2] = {1170,990};
            double numSides = 6;
            // choose starting vertex
            double startPoint[2];
            startPoint[0] = vertex_a[0];
            startPoint[1] = vertex_a[1];
            for (int k = 0; k < iterations; k++) {
                double random = 1.0*rand()/ RAND_MAX;
                if (random <= desired) {
                    double xDist = vertex_a[0] + startPoint[0];
                    double yDist = vertex_a[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                } else if (desired < random && random <= desired + (1-desired)/(numSides-1)) {
                    double xDist = vertex_b[0] + startPoint[0];
                    double yDist = vertex_b[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                } else if (desired + (1-desired)/(numSides-1) < random && random <= desired + 2*(1-desired)/(numSides-1)) {
                    double xDist = vertex_c[0] + startPoint[0];
                    double yDist = vertex_c[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                } else if (desired + 2*(1-desired)/(numSides-1) < random && random <= desired + 3*(1-desired)/(numSides-1)) {
                    double xDist = vertex_d[0] + startPoint[0];
                    double yDist = vertex_d[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                } else if (desired + 3*(1-desired)/(numSides-1) < random && random <= desired + 4*(1-desired)/(numSides-1)) {
                    double xDist = vertex_e[0] + startPoint[0];
                    double yDist = vertex_e[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                } else if (desired + 4*(1-desired)/(numSides-1) < random) {
                    double xDist = vertex_f[0] + startPoint[0];
                    double yDist = vertex_f[1] + startPoint[1];
                    double newPoint[2] = {floor(xDist * stepSize), floor(yDist * stepSize)};
                    data[(int) newPoint[0]][(int) newPoint[1]] = 1;
                    startPoint[0] = newPoint[0];
                    startPoint[1] = newPoint[1];
                }
            }
        }


        //Loops over all elements to work out which pixel has largest number (and therefore
        //darkest colour. The entire image is then re-normalised based on the difference between
        //the lightest and darkest pixels to ensure the best contrast
        double min = 9e9;
        double max = -9e9;
        for( int i=0 ; i < M ;i++ )
        {
            for( int j=0; j <N ; j++ )
            {
                if( data[i][j] < min )
                { min = data[i][j]; }
                if( data[i][j] > max )
                { max = data[i][j]; }
            }
        }

        std::cout << "Detected Cell Minimum: " << min << std::endl;
        std::cout << "Detected Cell Maximum: " << max << std::endl;

        //Create output bmp container
        BMP Output;
        Output.SetSize(M,N);
        Output.SetBitDepth(32);

        // plot the pixels
        for( int i=0 ; i < M ; i++ )
        {
            for( int j=0; j < N ; j++ )
            {
                double scaled_value = 1 - ( data[i][j] - min )/( max-min + 1e-16 );
                ebmpBYTE pixel_value = (ebmpBYTE) ( scaled_value * 255.0 );
                Output(i,j)->Red = pixel_value;   //
                Output(i,j)->Green = pixel_value; //for more interesting colours you can play around with these
                Output(i,j)->Blue = pixel_value;  //
            }
        }

        //Save output
        Output.WriteToFile( "ChaosGameOutput.bmp" );

        //Cleanup!

        //Remove data array from heap
        for(int i = 0; i < M; ++i){
            delete[] data[i];
        }

        //Remove array container
        delete[] data;


    } else if (sectionInput == 2){
    std::cout << "Generalised Chaos Game" << std::endl;
        int numdata = 0;

        std::vector<std::vector<double >> affine; //vector of vectors to store input transformations
        std::string line; //read in a data file
        std::ifstream myfile("fernIFS.txt");
        if (myfile.is_open()) {
            while (std::getline(myfile, line)) {
                std::vector<double> storedline;
                std::stringstream rowstream(line);
                double data_input;

                while (rowstream >> data_input) {
                    storedline.push_back(data_input); //populate inner vector
                    numdata++; // total data elements in file
                }
                affine.push_back(storedline); //populate outer vector with inner vectors
            }
            myfile.close();
        } else {  //terminate program to allow new input
            std::cout << "Unable to open file" << std::endl;
            return 0;
        }

        int numlines = numdata / 7;  //number of transformations equals the number of lines in the input file

        //Create a 2D array with MxN pixels on the heap
        int M = 1920;
        int N = 1200;
        double **data = new double *[M];

        for (int i = 0; i < M; ++i) {
            data[i] = new double[N];
        }

        //Initialise pixel array to zero (i.e. empty image)
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                data[i][j] = 0;
            }
        }

        //Check that all defined transformations are valid...
        if (validTrans(numlines,affine) != numlines){
            std::cout << "Not all transformations are valid" << std::endl;
            return 0;
        } else {
            std::cout << "All transformations are valid" << std::endl;
        }

        //User input number of pixels to generate...
        int iterations;
        std::cout << "Select the number of iterations. [For best resolution/computing time choose ~ 1000000]: " << std::endl;
        std::cin >> iterations;

        double start[2] = {};   //Define a start point on the image plane

        //Declare 2 vectors to store pixel information...
        std::vector<double> vecX;
        std::vector<double> vecY;

        //Weighted probabilities for each transformation...
        std::random_device random;
        std::mt19937 gen(random());
        std::uniform_real_distribution<>distr(0,1);
        for (int i = 0; i < iterations; ++i)
        {
            double x;
            double y;
            double probweight = 0;
            double rand = distr(gen);
            for (int j = 0; j < numlines; j++) {
                double prevProb = probweight;
                probweight += affine[j][6];
                if (rand >= prevProb && rand < probweight) { //probability condition to determine which transformation is applied
                    x = affine[j][0] * start[0] + affine[j][1] * start[1] + affine[j][4];
                    y = affine[j][2] * start[0] + affine[j][3] * start[1] + affine[j][5];
                    start[0] = x;
                    start[1] = y;
                    vecX.push_back(start[0]); //populate vectors with transformed pixel coordinates
                    vecY.push_back(start[1]);
                }
            }
        }

        //scaling in x-direction
        double xMax = getMax(vecX);
        double xMin = getMin(vecX);
        double xScale = M/(xMax-xMin);
        Trans(vecX,0); //transform the points

        //scaling in y-direction
        double yMax = getMax(vecY);
        double yMin = getMin(vecY);
        double yScale = N/(yMax-yMin);
        Trans(vecY,0); //transform the points

        //Populate the data array with correctly scaled pixel coordinates...
        for (int k = 0; k < iterations; ++k)
        {
            if((int)(vecX[k]*xScale) >= M){
                continue;
            }
            if((int)(vecY[k]*yScale) >= N){
                continue;
            }
            data[(int)(vecX[k]*xScale)][(int)(vecY[k]*yScale)] = 1;
        }

        //Loops over all elements to work out which pixel has largest number (and therefore
        //darkest colour. The entire image is then re-normalised based on the difference between
        //the lightest and darkest pixels to ensure the best contrast
        double min = 9e9;
        double max = -9e9;
        for( int i=0 ; i < M ;i++ )
        {
            for( int j=0; j <N ; j++ )
            {
                if( data[i][j] < min )
                { min = data[i][j]; }
                if( data[i][j] > max )
                { max = data[i][j]; }
            }
        }

        std::cout << "Detected Cell Minimum: " << min << std::endl;
        std::cout << "Detected Cell Maximum: " << max << std::endl;

        //Create output bmp container
        BMP Output;
        Output.SetSize(M,N);
        Output.SetBitDepth(32);

        // plot the pixels
        for( int i=0 ; i < M ; i++ )
        {
            for( int j=0; j < N ; j++ )
            {
                double scaled_value = 1 - ( data[i][j] - min )/( max-min + 1e-16 );
                ebmpBYTE pixel_value = (ebmpBYTE) ( scaled_value * 255.0 );
                Output(i,j)->Red = pixel_value;   //
                //Output(i,j)->Green = pixel_value; //for more interesting colours you can play around with these
                Output(i,j)->Blue = pixel_value;  //
            }
        }

        //Save output
        Output.WriteToFile( "GeneralisedChaosOutput.bmp" );

        //Cleanup!

        //Remove data array from heap
        for(int i = 0; i < M; ++i){

            delete[] data[i];
        }

        //Remove array container
        delete[] data;


    } else if (sectionInput ==3){
        std::cout << "Collage Theorem" << std::endl;
        int userInput;
        std::cout << "Enter '1' to see transformation modelling with the collage theorem. Enter '2' to reproduce the leaf with IFS code." << std::endl;
        std::cin >> userInput;
        if (userInput == 1) {
            //read in BMP image...
            BMP ImageIn;
            ImageIn.ReadFromFile("maple_leaf.bmp");
            const unsigned int widthPixels = ImageIn.TellWidth();
            const unsigned int heightPixels = ImageIn.TellHeight();

            double **data = new double *[widthPixels];

            for (unsigned int xPixel = 0; xPixel < widthPixels; ++xPixel) {

                data[xPixel] = new double[heightPixels];

                for (unsigned int yPixel = 0; yPixel < heightPixels; ++yPixel) {
                    // convert each pixel to greyscale
                    double Temp = 0.30 * (ImageIn(xPixel, yPixel)->Red) + 0.59 * (ImageIn(xPixel, yPixel)->Green) +
                                  0.11 * (ImageIn(xPixel, yPixel)->Blue);
                    if (Temp != 255) {
                        data[xPixel][yPixel] = (Temp);
                    } else {
                        data[xPixel][yPixel] = 0;
                    }
                }
            }

            double min = 9e9;
            double max = -9e9;
            for (unsigned int i = 0; i < widthPixels; i++) {
                for (unsigned int j = 0; j < heightPixels; j++) {
                    if (data[i][j] < min) { min = data[i][j]; }
                    if (data[i][j] > max) { max = data[i][j]; }
                }
            }

            double **transData = new double *[widthPixels];
            for (unsigned int i = 0; i < widthPixels; ++i) {
                transData[i] = new double[heightPixels];
            }

            for (unsigned int i = 0; i < widthPixels; ++i) {
                for (unsigned int j = 0; j < heightPixels; ++j) {
                    transData[i][j] = data[i][j];
                }
            }

            double transcale = 0.2; //Size scaling of the transformed images

            std::cout << "Detected Cell Minimum: " << min << std::endl;
            std::cout << "Detected Cell Maximum: " << max << std::endl;


            //Declare some vectors to be populated by information from transformed images...
            std::vector<double> dataNew;
            std::vector<double> vecX, vecY, vec2X, vec2Y, vec3X, vec3Y, vec4X, vec4Y, vec5X, vec5Y, vec6X, vec6Y, vec7X, vec7Y, vec8X, vec8Y, vec9X, vec9Y, vec10X, vec10Y;

            for (unsigned int xPixel = 0; xPixel < widthPixels; ++xPixel) {
                for (unsigned int yPixel = 0; yPixel <
                                              heightPixels; ++yPixel) {    //All pixels inside the 2D grid are accounted for in this double loop
                    if (data[xPixel][yPixel] > 0) {     //Only filled pixels are transformed
                        dataNew.push_back(data[xPixel][yPixel]);
                        //10 different transformations (rotation and scaling only)...
                        vecX.push_back((xPixel * transcale * 0.707) + (transcale * 0.707 * yPixel));
                        vecY.push_back((xPixel * transcale * -0.707) + (transcale * 0.707 * yPixel));
                        vec2X.push_back((xPixel * transcale * 0.707) + (transcale * -0.707 * yPixel));
                        vec2Y.push_back((xPixel * transcale * 0.707) + (transcale * 0.707 * yPixel));
                        vec3X.push_back((xPixel * transcale * 0.174) + (transcale * -0.985 * yPixel));
                        vec3Y.push_back((xPixel * transcale * 0.985) + (transcale * 0.174 * yPixel));
                        vec4X.push_back((xPixel * transcale * 0.174) + (transcale * 0.985 * yPixel));
                        vec4Y.push_back((xPixel * transcale * -0.985) + (transcale * 0.174 * yPixel));
                        vec5X.push_back((xPixel * transcale * 1) + (transcale * 0 * yPixel));
                        vec5Y.push_back((xPixel * transcale * 0) + (transcale * 1 * yPixel));
                        vec6X.push_back((xPixel * transcale * 1) + (transcale * 0 * yPixel));
                        vec6Y.push_back((xPixel * transcale * 0) + (transcale * 1 * yPixel));
                        vec7X.push_back((xPixel * transcale * 0.707) + (transcale * 0.707 * yPixel));
                        vec7Y.push_back((xPixel * transcale * -0.707) + (transcale * 0.707 * yPixel));
                        vec8X.push_back((xPixel * transcale * 0.707) + (transcale * -0.707 * yPixel));
                        vec8Y.push_back((xPixel * transcale * 0.707) + (transcale * 0.707 * yPixel));
                        vec9X.push_back((xPixel * transcale * 0.174) + (transcale * -0.985 * yPixel));
                        vec9Y.push_back((xPixel * transcale * 0.985) + (transcale * 0.174 * yPixel));
                        vec10X.push_back((xPixel * transcale * 0.174) + (transcale * 0.985 * yPixel));
                        vec10Y.push_back((xPixel * transcale * -0.985) + (transcale * 0.174 * yPixel));
                    }
                }
            }

            //Perform the transformations and add the translations...
            Trans(vecX, 300);
            Trans(vecY, 300);
            Trans(vec2X, 1350);
            Trans(vec2Y, 300);
            Trans(vec3X, 1350);
            Trans(vec3Y, 825);
            Trans(vec4X, 300);
            Trans(vec4Y, 850);
            Trans(vec5X, 800);
            Trans(vec5Y, 100);
            Trans(vec6X, 800);
            Trans(vec6Y, 800);
            Trans(vec7X, 500);
            Trans(vec7Y, 550);
            Trans(vec8X, 1150);
            Trans(vec8Y, 550);
            Trans(vec9X, 1100);
            Trans(vec9Y, 925);
            Trans(vec10X, 600);
            Trans(vec10Y, 925);

            //Assign new pixels from all the transformations to a new 2D grid placed on top of the original...
            for (unsigned int k = 0; k < vecX.size(); ++k) {
                transData[(int) (vecX[k])][(int) (vecY[k])] = dataNew[k];
                transData[(int) (vec2X[k])][(int) (vec2Y[k])] = dataNew[k];
                transData[(int) (vec3X[k])][(int) (vec3Y[k])] = dataNew[k];
                transData[(int) (vec4X[k])][(int) (vec4Y[k])] = dataNew[k];
                transData[(int) (vec5X[k])][(int) (vec5Y[k])] = dataNew[k];
                transData[(int) (vec6X[k])][(int) (vec6Y[k])] = dataNew[k];
                transData[(int) (vec7X[k])][(int) (vec7Y[k])] = dataNew[k];
                transData[(int) (vec8X[k])][(int) (vec8Y[k])] = dataNew[k];
                transData[(int) (vec9X[k])][(int) (vec9Y[k])] = dataNew[k];
                transData[(int) (vec10X[k])][(int) (vec10Y[k])] = dataNew[k];
            }

            //Create output bmp container
            BMP Output;
            Output.SetSize(widthPixels, heightPixels);
            Output.SetBitDepth(32);


            for (unsigned int i = 0; i < widthPixels; i++) {
                for (unsigned int j = 0; j < heightPixels; j++) {
                    double scaled_value = 1 - (transData[i][j] - min) / (max - min + 1e-16);
                    ebmpBYTE pixel_value = (ebmpBYTE) (scaled_value * 255.0);
                    Output(i, j)->Red = pixel_value;   //
                    Output(i, j)->Green = pixel_value; //for more interesting colours you can play around with these
                    Output(i, j)->Blue = pixel_value;  //
                }
            }

            Output.WriteToFile("CollageOutput.bmp");

            for (unsigned int i = 0; i < widthPixels; ++i) {
                delete[] data[i];
            }
            delete[] data;

            for (unsigned int i = 0; i < widthPixels; ++i) {
                delete[] transData[i];
            }
            delete[] transData;

        } else if (userInput == 2){

            int numdata = 0;

            std::vector<std::vector<double >> affine;
            std::string line;
            std::ifstream myfile("leafifs10trans.txt");
            if (myfile.is_open()) {
                while (std::getline(myfile, line)) {
                    std::vector<double> storedline;
                    std::stringstream rowstream(line);
                    double data_input;
                    while (rowstream >> data_input) {
                        storedline.push_back(data_input);
                        numdata++;
                    }
                    affine.push_back(storedline);
                }
                myfile.close();
            } else {
                std::cout << "Unable to open file" << std::endl;
                return 0;
            }
            int numlines = numdata / 7;  //Number of transformations equals number of lines in file

            //Create a 2D array with MxN pixels on the heap
            int M = 1920;
            int N = 1200;
            double **data = new double *[M];

            for (int i = 0; i < M; ++i) {
                data[i] = new double[N];
            }
            //Initialise pixel array to zero (i.e. empty image)
            for (int i = 0; i < M; ++i) {
                for (int j = 0; j < N; ++j) {
                    data[i][j] = 0;
                }
            }

            //User input number of pixels to generate...
            int iterations;
            std::cout << "Select the number of iterations. [For best resolution/computing time choose ~ 1000000]: " << std::endl;
            std::cin >> iterations;

            double start[2] = {};   //Define a start point on the image plane
            double transcale = 0.3; //Scaling of IFS output. This allows the matrix to be defined purely for rotation.
            // The scaling can be altered here, so as all the smaller images are the same
            // size as eachother.

            //Declare 2 vectors to store pixel information...
            std::vector<double> vecX;
            std::vector<double> vecY;

            //Weighted probabilities for each transformation...
            std::random_device random;
            std::mt19937 gen(random());
            std::uniform_real_distribution<>distr(0,1);
            for (int i = 0; i < iterations; ++i)
            {
                double x;
                double y;
                double probweight = 0;
                double rand = distr(gen);
                for (int j = 0; j < numlines; j++) {
                    double prevProb = probweight;
                    probweight += affine[j][6];
                    if (rand >= prevProb && rand < probweight) { //IFS code for selected transformation...
                        x = transcale*affine[j][0] * start[0] + transcale*affine[j][1] * start[1] + affine[j][4];
                        y = transcale*affine[j][2] * start[0] + transcale*affine[j][3] * start[1] + affine[j][5];
                        start[0] = x;
                        start[1] = y;
                        vecX.push_back(start[0]); //populate the vectors containing pixel information
                        vecY.push_back(start[1]);
                    }
                }
            }

            //scaling in x-direction
            double xMax = getMax(vecX);
            double xMin = getMin(vecX);
            double xScale = M/(xMax-xMin);
            Trans(vecX,0); //Transform the pixels in x-direction

            //scaling in y-direction
            double yMax = getMax(vecY);
            double yMin = getMin(vecY);
            double yScale = N/(yMax-yMin);
            Trans(vecY,0); //Transform the pixels in y-direction

            //Populate the data array with correctly scaled pixel coordinates...
            for (int k = 0; k < iterations; ++k)
            {
                if((int)(vecX[k]*xScale) >= M){
                    continue;
                }
                if((int)(vecY[k]*yScale) >= N){
                    continue;
                }
                data[(int)(vecX[k]*xScale)][(int)(vecY[k]*yScale)] = 1;
            }

            //Loops over all elements to work out which pixel has largest number (and therefore
            //darkest colour. The entire image is then re-normalised based on the difference between
            //the lightest and darkest pixels to ensure the best contrast
            double min = 9e9;
            double max = -9e9;
            for( int i=0 ; i < M ;i++ )
            {
                for( int j=0; j <N ; j++ )
                {
                    if( data[i][j] < min )
                    { min = data[i][j]; }
                    if( data[i][j] > max )
                    { max = data[i][j]; }
                }
            }

            std::cout << "Detected Cell Minimum: " << min << std::endl;
            std::cout << "Detected Cell Maximum: " << max << std::endl;

            //Create output bmp container
            BMP Output;
            Output.SetSize(M,N);
            Output.SetBitDepth(32);

            // plot the pixels
            for( int i=0 ; i < M ; i++ )
            {
                for( int j=0; j < N ; j++ )
                {
                    double scaled_value = 1 - ( data[i][j] - min )/( max-min + 1e-16 );
                    ebmpBYTE pixel_value = (ebmpBYTE) ( scaled_value * 255.0 );
                    Output(i,j)->Red = pixel_value;   //
                    Output(i,j)->Green = pixel_value; //for more interesting colours you can play around with these
                    Output(i,j)->Blue = pixel_value;  //
                }
            }

            //Save output
            Output.WriteToFile( "CollageIFSAttractor.bmp" );

            //Cleanup!

            //Remove data array from heap
            for(int i = 0; i < M; ++i){

                delete[] data[i];
            }

            //Remove array container
            delete[] data;

        } else {
            std::cout << "Invalid input" << std::endl;
            return 0;
        }


    } else {
        std::cout << "You have not selected a valid section. Run again." << std::endl;
        return 0;
    }

}