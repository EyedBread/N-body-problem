#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <omp.h>
#include <bits/stdc++.h>
#include "Matrix.h"

using namespace std;



struct point{
    double x=0;
    double y=0;

    point(double x_, double y_) : x(x_), y(y_) {}

    point() : x(0), y(0) {}
};

//initialize the positions, velocities, forces, and masses
double G = 6.67e-11; //gravity

#define DT 0.01
#define SIZE 9e7
#define NUMITER 5


//Generate random double
double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


//calculate total force for every pair of bodies
void calculateForces(vector<point>& p, Matrix<point>& f, vector<double>& m, const int& n){
    

    #pragma omp parallel for schedule(guided) //Stripes evenly distributed
    for (int i=0; i<n-1;i++){
        double distance;
        double magnitude;
        point direction;
        int id = omp_get_thread_num();
        for (int j = i+1; j < n; j++) {
            distance =  sqrt(pow((p[i].x - p[j].x),2) + pow((p[i].y - p[j].y),2));
            magnitude = (G*m[i]*m[j]) / pow(distance,2);
            direction.x = p[j].x-p[i].x; 
            direction.y = p[j].y-p[i].y;

            //False sharing possible
            f(id,i).x = f(id,i).x + magnitude*direction.x/distance;

            f(id, j).x = f(id, j).x - magnitude*direction.x/distance;

            f(id, i).y = f(id, i).y + magnitude*direction.y/distance;
            
            f(id, j).y = f(id, j).y - magnitude*direction.y/distance;
        }
    }
}

//calculate new velocity and position for each body
void moveBodies(vector<point>& p, Matrix<point>& f, vector<double>& m, vector<point>& v, const int& n) {
    point deltav; // dv=f/m * DT
    point deltap; // dp=(v+dv/2) * DT
    point force;

    #pragma omp parallel for private(deltav, deltap) schedule(static)
    for (int i = 0; i < n; i++) {
        point force;
        for (int k = 0; k < omp_get_num_threads(); k++) {
            force.x += f(k,i).x; f(k,i).x = 0;
            force.y += f(k,i).y; f(k,i).y = 0;
        }
        deltav.x = force.x/m[i] * DT;
        deltav.y = force.y/m[i] * DT;
        deltap.x = (v[i].x + deltav.x/2) * DT;
        deltap.y = (v[i].y + deltav.y/2) * DT;
        v[i].x = v[i].x + deltav.x;
        v[i].y = v[i].y + deltav.y;
        p[i].x = p[i].x + deltap.x;
        p[i].y = p[i].y + deltap.y;
        // if (i == 1)
        //     cout << "pos x: " << p[i].x << ", pos y: " << p[i].y << endl;
        force.x = force.y = 0.0; // reset force vector
    }
}

int main(int argc, char** argv)
{

    if (argc != 4) {
        cout << "Incorrect amount of arguments!" << endl;
        return 0;
    }

    srand(time(NULL));

    int gnumBodies = atoi(argv[1]);
    int numSteps = atoi(argv[2]);
    int nThreads = atoi(argv[3]);

    cout << "number of bodies: " << gnumBodies << endl;
    cout << "number of timesteps: " << numSteps << endl;
    cout << "number of threads: " << nThreads << endl; 

    omp_set_num_threads(nThreads);
    vector<float> med(NUMITER); 
    for (int j = 0; j < NUMITER; j++) {
        vector<point> p; //Position
        vector<point> v; //velocity
        Matrix<point> f(nThreads, gnumBodies); //force
        // vector<vector<point>> f(nThreads, vector<point>(gnumBodies)); //force
        vector<double> m; //mass

        for (int i = 0; i < gnumBodies; i++) {
            p.push_back(point(fRand(0, SIZE), fRand(0, SIZE)));
            v.push_back(point(0,0));
            m.push_back(fRand(6e22, 2e30));
        }

        auto start = omp_get_wtime();
        for (int i=0; i<numSteps; i++){
            calculateForces(p, f, m, gnumBodies);
            // cout << "hej" << endl;
            moveBodies(p, f, m, v, gnumBodies);
        }
        med[j] = omp_get_wtime() - start;
    }
    sort(med.begin(), med.end());

    float medianTime;
    if (NUMITER % 2 == 0) {
        medianTime = (med[NUMITER/2-1] + med[NUMITER/2])/2;
    }
    else {
        medianTime = med[NUMITER/2];
    }

    cout << "Time taken was " << medianTime << " s"<< endl;
}