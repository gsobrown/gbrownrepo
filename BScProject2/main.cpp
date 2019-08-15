
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>

int main() {
    std::string typeOfSim;
    std::cout << "Enter 'MD' or 'SPH' for the type of simulation." << std::endl;
    std::cin >> typeOfSim;
    if (typeOfSim == "MD") {
        std::cout << "This is a molecular dynamics simulation" << std::endl;
    } else if (typeOfSim == "SPH") {
        std::cout << "This is a smooth particle hydrodynamics simulation" << std::endl;
    } else {
        std::cout << "Invalid input" << std::endl;
        return 0;
    }

    bool thermostat;
    std::string thermostatBool;
    std::cout << "Thermostat 'ON' or 'OFF'?" << std::endl;
    std::cin >> thermostatBool;
    if (thermostatBool == "ON") {
        std::cout << "Thermostat is on" << std::endl;
        thermostat = true;
    } else if (thermostatBool == "OFF"){
        std::cout << "Thermostat is off" << std::endl;
        thermostat = false;
    } else {
        std::cout << "Invalid input" << std::endl;
        return 0;
    }

    std::ofstream densityfile;
    densityfile.open("densityfile.txt");

    double i1[100];
    double i2[100];
    double i3[100];
    double i4[100];
    double i5[100];
    double i6[100];
    double grr1[100];
    double grr2[100];
    double grr3[100];
    double grr4[100];
    double grr5[100];
    double grr6[100];

    for (double rho = 0.2; rho <= 1.2; rho+=0.1) {

        double m = 5;
        double N = pow(m, 3);
        //double rho = 0.5;
        double temp = 2.0;
        double vol = N / rho;
        double S = pow(vol, 1. / 3);
        double S5 = S / 2;
        double mx = 0;
        double my = 0;
        double mz = 0;
        double KEbox = 0;
        double KE = 0;
        double cu = 2.5;
        double cu2 = cu * cu;
        double rx[(int) N];
        double ry[(int) N];
        double rz[(int) N];
        double vx[(int) N];
        double vy[(int) N];
        double vz[(int) N];
        double vx1[(int) N];
        double vy1[(int) N];
        double vz1[(int) N];
        double u = 0.0;
        double p = 0.0;
        double DT = 0.005;
        double KEverlet = 0.0;
        double f;
        double r;
        double time = 0;
        double pk;
        double ptotal;
        double tempnow;
        double ke;
        double e;
        double pc;
        double ndist = 100;
        double fixgr = ndist / cu;
        double dr = 1 / fixgr;
        double dvol;
        double gr1;
        double gr[10000] = {};
        std::vector<double> UVECTOR;
        std::vector<double> PVECTOR;
        std::vector<double> KEVECTOR;
        std::vector<double> ETOTVECTOR;
        std::vector<double> TEMPVECTOR;
        double USUM = 0;
        double KESUM = 0;
        double ETOTSUM = 0;
        double TEMPSUM = 0;
        double PSUM = 0;
        double entries = 0;
        double rxDiff[125] = {};
        double ryDiff[125] = {};
        double rzDiff[125] = {};
        double esigmatop = 0;
        double kesigmatop = 0;
        double usigmatop = 0;
        double psigmatop = 0;
        double tsigmatop = 0;

        int ii = 0;
        for (int i = 0; i < m; i++) {

            for (int j = 0; j < m; j++) {

                for (int k = 0; k < m; k++) {

                    rx[ii] = S * (i + 1 - 0.5) / m;
                    ry[ii] = S * (j + 1 - 0.5) / m;
                    rz[ii] = S * (k + 1 - 0.5) / m;
                    vx[ii] = ((double) rand() / (RAND_MAX)) - 0.5;
                    vy[ii] = ((double) rand() / (RAND_MAX)) - 0.5;
                    vz[ii] = ((double) rand() / (RAND_MAX)) - 0.5;
                    mx += vx[ii];
                    my += vy[ii];
                    mz += vz[ii];
                    ii++;
                }
            }
        }

        //**************************************************************************

        //c.......... subtract off the momentum per particle from that of each particle
        for (int i = 0; i < N; i++) {
            vx[i] -= mx / N;
            vy[i] -= my / N;
            vz[i] -= mz / N;
            KEbox += 0.5 * (pow(vx[i], 2) + pow(vy[i], 2) + pow(vz[i], 2));
        }


        for (int i = 0; i < N; i++) {
            vx1[i] = vx[i] * sqrt(1.5 * N * temp / KEbox);
            vy1[i] = vy[i] * sqrt(1.5 * N * temp / KEbox);
            vz1[i] = vz[i] * sqrt(1.5 * N * temp / KEbox);
            KE += 0.5 * (pow(vx1[i], 2) + pow(vy1[i], 2) + pow(vz1[i], 2));
            vx[i] = vx1[i];
            vy[i] = vy1[i];
            vz[i] = vz1[i];
        }

        double temp1 = KE / (1.5 * N);
        std::cout << "************* SETUP **************" << std::endl;
        std::cout << "initial temp: " << temp << std::endl;
        std::cout << "new temp: " << temp1 << std::endl;
        std::cout << "KE of box: " << KEbox << std::endl;
        std::cout << "KE of molecules: " << KE << std::endl;
        std::cout << "**********************************" << std::endl;

        //**************** End of set up section *********************************
        std::ofstream myfile;
        myfile.open("datafile.txt");


        std::ofstream velocityFile;
        velocityFile.open("velocityFile.txt");

        //......BIG LOOP........
        int nt = 5000;

        for (int it = 0; it < nt; it++) {
            time += DT;
            double FX[10000] = {};
            double FY[10000] = {};
            double FZ[10000] = {};
            u = 0;
            p = 0;
            KEverlet = 0;

            //  This is the forces double loop...
            for (int i = 0; i < N; i++) {
                if (it == 0) { //nt-1
                    velocityFile << rx[i] << std::endl;
                    velocityFile << ry[i] << std::endl;
                    velocityFile << rz[i] << std::endl;
                    velocityFile << vx[i] << std::endl;
                    velocityFile << vy[i] << std::endl;
                    velocityFile << vz[i] << std::endl;
                }

                for (int j = i + 1; j < N; j++) {
                    double x = rx[i] - rx[j];
                    double y = ry[i] - ry[j];
                    double z = rz[i] - rz[j];
                    if (x >= S5) {
                        x = x - S;
                    }
                    if (y >= S5) {
                        y = y - S;
                    }
                    if (z >= S5) {
                        z = z - S;
                    }
                    if (x <= -S5) {
                        x = x + S;
                    }
                    if (y <= -S5) {
                        y = y + S;
                    }
                    if (z <= -S5) {
                        z = z + S;
                    }
                    double rr = pow(x, 2) + pow(y, 2) + pow(z, 2);

                    if (rr < cu2) {
                        double rri = 1.0 / rr;
                        double r6i = pow(rri, 3.0);
                        double r12i = pow(r6i, 2.0);
                        if (typeOfSim == "MD") {
                            u += 4.0 * (r12i - r6i);
                            f = 24.0 * (2 * r12i - r6i) * rri; //force divided by r
                        } else if (typeOfSim == "SPH") {
                            double nparam = 12.0;
                            double A = 1 / pow(2, 1. / nparam);
                            u += 1 / pow((A + rr), nparam / 2.0);
                            f = nparam * sqrt(rr) / pow((A + rr), (nparam / 2.0 + 1.0)) * rri;
                        }
                        //.....cartesian components of the force...
                        double FXC = x * f;
                        double FYC = y * f;
                        double FZC = z * f;
                        FX[i] += FXC;
                        FY[i] += FYC;
                        FZ[i] += FZC;
                        FX[j] -= FXC;
                        FY[j] -= FYC;
                        FZ[j] -= FZC;
                        if (typeOfSim == "MD") {
                            p += 24.0 * (2 * r12i - r6i); //pressure
                        } else if (typeOfSim == "SPH") {
                            double nparam = 12.0;
                            double A = 1 / pow(2, 1. / nparam);
                            p += nparam * sqrt(rr) / pow((A + rr), (nparam / 2.0 + 1.0));
                        }
                        //c.....radial distribution function...
                        r = sqrt(rr);
                        double ir = (int) floor(r * fixgr);
                        gr[(int) ir] += 1.0;
                    } else {
                        u += 0.0;
                    }
                    pc = p / (3.0 * pow(S, 3)); //configurational pressure
                }
            }

            double VSCALE;

            if (thermostat == true) {
                KEbox = 0.0;
                for (int i = 0; i < N; i++) {
                    KEbox += 0.5 * (pow(vx[i], 2) + pow(vy[i], 2) + pow(vz[i], 2));
                }
                VSCALE = sqrt(1.5 * N * temp / KEbox);
            } else if (thermostat == false) {
                KEbox = 0.0;
                VSCALE = 1.0;
            }

            for (int i = 0; i < N; i++) {

                //c.....apply Verlet's Leapfrog algorithm
                vx[i] = vx[i] * VSCALE + (FX[i] * DT);
                vy[i] = vy[i] * VSCALE + (FY[i] * DT);
                vz[i] = vz[i] * VSCALE + (FZ[i] * DT);
                rx[i] = rx[i] + (vx[i] * DT);
                ry[i] = ry[i] + (vy[i] * DT);
                rz[i] = rz[i] + (vz[i] * DT);

                //c.....accumulate kinetic energy
                KEverlet += 0.5 * (pow(vx[i], 2) + pow(vy[i], 2) + pow(vz[i], 2));

                //c.....check if a particle has moved out of the box and then apply p.b.c.
                if (rx[i] >= S) {
                    rx[i] = rx[i] - S;
                }
                if (ry[i] >= S) {
                    ry[i] = ry[i] - S;
                }
                if (rz[i] >= S) {
                    rz[i] = rz[i] - S;
                }
                if (rx[i] <= 0.0) {
                    rx[i] = rx[i] + S;
                }
                if (ry[i] <= 0.0) {
                    ry[i] = ry[i] + S;
                }
                if (rz[i] <= 0.0) {
                    rz[i] = rz[i] + S;
                }
                if (rx[i] > S || ry[i] > S || rz[i] > S || rx[i] < 0 || ry[i] < 0 || rz[i] < 0) {
                    std::cout << "OUTSIDE BOX" << std::endl;
                }
                rxDiff[i] += (vx[i] * DT);
                ryDiff[i] += (vy[i] * DT);
                rzDiff[i] += (vz[i] * DT);
            }

            pk = 2.0 * (KEverlet) / (3.0 * pow(S, 3));   //kinetic pressure
            ptotal = pk + pc;
            tempnow = (KEverlet) / (1.5 * N);
            double uN = u / N;    //potential energy per particle
            ke = (KEverlet) / N;  //ke per particle
            e = uN + ke;


            myfile << time << std::endl;
            myfile << tempnow << std::endl;
            myfile << uN << std::endl;
            myfile << ke << std::endl;
            myfile << e << std::endl;
            myfile << ptotal << std::endl;
            UVECTOR.push_back(uN);
            PVECTOR.push_back(ptotal);
            KEVECTOR.push_back(ke);
            ETOTVECTOR.push_back(e);
            TEMPVECTOR.push_back(tempnow);

            // go to BIG LOOP again
        }
        myfile.close();

        velocityFile.close();


        std::ofstream radialfile;
        radialfile.open("radialdistfunc.txt");

        for (int i = 0; i < ndist; i++) {
            r = dr * i;
            dvol = (2.0 * M_PI / 3.0) * ((pow((r + (dr / 2)), 3) - pow((r - (dr / 2)), 3)));
            gr1 = gr[i] / (dvol * rho * nt * N);
            if (rho == 0.2) {
                i1[i] = i;
                grr1[i] = gr1;
            } else if (rho == 0.5) {
                i2[i] = i;
                grr2[i] = gr1;
            } else if (rho == 0.7){ //PROBLEM
                i3[i] = i;
                grr3[i] += gr1;
            } else if (rho == 1.0) { //PROBLEM 0.7
                i4[i] = i;
                grr4[i] = gr1;
            } else if (rho == 1.2) {
                i5[i] = i;
                grr5[i] = gr1;

            }

        }

        for (int i = 0; i < ndist; i++) {
            radialfile << 0.025*i << std::endl;
            radialfile << grr1[i] << std::endl;
        }
        radialfile.close();

        std::ofstream radialfile2;
        radialfile2.open("radial2.txt");

        for (int i = 0; i < ndist; i++){
            radialfile2 << 0.025*i2[i] << std::endl;
            radialfile2 << grr2[i] << std::endl;
        }
        radialfile2.close();

        std::ofstream radialfile3;
        radialfile3.open("radial3.txt");

        for (int i = 0; i < ndist; i++){
            radialfile3 << 0.025*i3[i] << std::endl;
            radialfile3 << grr3[i] << std::endl;
        }
        radialfile3.close();

        std::ofstream radialfile4;
        radialfile4.open("radial4.txt");

        for (int i = 0; i < ndist; i++){
            radialfile4 << 0.025*i4[i] << std::endl;
            radialfile4 << grr4[i] << std::endl;
        }
        radialfile4.close();

        std::ofstream radialfile5;
        radialfile5.open("radial5.txt");

        for (int i = 0; i < ndist; i++){
            radialfile5 << 0.025*i5[i] << std::endl;
            radialfile5 << grr5[i] << std::endl;
        }
        radialfile5.close();

        std::ofstream radialfile6;
        radialfile6.open("radial6.txt");

        for (int i = 0; i < ndist; i++){
            radialfile6 << 0.025*i6[i] << std::endl;
            radialfile6 << grr6[i] << std::endl;
        }
        radialfile6.close();


        for (int n = 1999; n < nt; n++) {
            USUM += UVECTOR[n];
            PSUM += PVECTOR[n];
            KESUM += KEVECTOR[n];
            ETOTSUM += ETOTVECTOR[n];
            TEMPSUM += TEMPVECTOR[n];
            entries++;
        }

        double UAVG = USUM / entries;
        double PAVG = PSUM / entries;
        double KEAVG = KESUM / entries;
        double ETOTAVG = ETOTSUM / entries;
        double TEMPAVG = TEMPSUM / entries;

        for (int i = 0; i < UVECTOR.size(); i++) {
            usigmatop += pow((UVECTOR[i] - UAVG), 2);
            psigmatop += pow((PVECTOR[i] - PAVG), 2);
            kesigmatop += pow((KEVECTOR[i] - KEAVG), 2);
            esigmatop += pow((ETOTVECTOR[i] - ETOTAVG), 2);
            tsigmatop += pow((TEMPVECTOR[i] - TEMPAVG), 2);
        }

        double ustdev = sqrt(usigmatop / UVECTOR.size());
        double userr = ustdev / (N);
        double kestdev = sqrt(kesigmatop / UVECTOR.size());
        double keserr = kestdev / (N);
        double estdev = sqrt(esigmatop / UVECTOR.size());
        double eserr = estdev / (N);
        double tstdev = sqrt(tsigmatop / UVECTOR.size());
        double tserr = tstdev / (N);
        double pstdev = sqrt(psigmatop / UVECTOR.size());
        double pserr = pstdev / (N);

        std::cout << "Density: " << rho << std::endl;
        std::cout << "Average potential energy: " << UAVG << " (std error: " << userr << ")" << std::endl;
        std::cout << "Average kinetic energy: " << KEAVG << " (std error: " << keserr << ")" << std::endl;
        std::cout << "Average total energy: " << ETOTAVG << " (std error: " << eserr << ")" << std::endl;
        std::cout << "Average temperature: " << TEMPAVG << " (std error: " << tserr << ")" << std::endl;
        std::cout << "Average pressure: " << PAVG << " (std error: " << pserr << ")" << std::endl;
        std::cout << "Total time: " << time << std::endl;

        std::ofstream avgfile;
        avgfile.open("avgfile.txt");
        avgfile << UAVG << std::endl;
        avgfile << KEAVG << std::endl;
        avgfile << ETOTAVG << std::endl;
        avgfile << TEMPAVG << std::endl;
        avgfile << PAVG << std::endl;
        avgfile.close();


        densityfile << rho << std::endl;
        densityfile << UAVG << std::endl;
        densityfile << PAVG << std::endl;

        std::vector<double> diffvector;
        double diffsigmatop = 0;
        double diffCoeffSum = 0.0;
        double diffCoeffSum2 = 0.0;

        for (int i = 0; i < N; i++) {
            double distance = sqrt(pow(rxDiff[i], 2) + pow(ryDiff[i], 2) + pow(rzDiff[i], 2));
            double diffCoeff = pow(distance, 2) / (6 * time);
            diffCoeffSum += diffCoeff;
            diffvector.push_back(diffCoeff);
        }

        double diffusionC = diffCoeffSum / N;

        std::cout << "***************************************" << std::endl;
        std::cout << "Self diffusion coefficient: " << diffusionC << std::endl;

        for (int i = 0; i < diffvector.size(); i++) {
            diffsigmatop += pow((diffvector[i] - diffusionC), 2);
        }

        double diffstdev = sqrt(diffsigmatop / diffvector.size());
        double diffserr = diffstdev / sqrt(N);

        std::cout << "Variance: " << pow(diffstdev,2) << std::endl;
        std::cout << "Standard deviation: " << diffstdev << std::endl;
        std::cout << "Standard error: " << diffserr << std::endl;
        densityfile << diffusionC << std::endl;
        densityfile << diffserr << std::endl;

    }
    densityfile.close();

}








