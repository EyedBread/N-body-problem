#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>

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


//Generate random double
double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}


//calculate total force for every pair of bodies
void calculateForces(vector<point>& p, vector<point>& f, vector<double>& m, int n){
    double distance;
    double magnitude;
    point direction;
    for (int i=0; i<n-1;i++){
        for (int j = i+1; j < n; j++) {
            distance =  sqrt(pow((p[i].x - p[j].x),2) + pow((p[i].y - p[j].y),2));
            magnitude = (G*m[i]*m[j]) / pow(distance,2);
            direction.x = p[j].x-p[i].x; 
            direction.y = p[j].y-p[i].y;
            f[i].x = f[i].x + magnitude*direction.x/distance;
            f[j].x = f[j].x - magnitude*direction.x/distance;
            f[i].y = f[i].y + magnitude*direction.y/distance;
            f[j].y = f[j].y - magnitude*direction.y/distance;
        }
    }
}

//calculate new velocity and position for each body
void moveBodies(vector<point>& p, vector<point>& f, vector<double>& m, vector<point>& v, int n) {
    point deltav; // dv=f/m * DT
    point deltap; // dp=(v+dv/2) * DT
    for (int i = 0; i < n; i++) {
        deltav.x = f[i].x/m[i] * DT;
        deltav.y = f[i].y/m[i] * DT;
        deltap.x = (v[i].x + deltav.x/2) * DT;
        deltap.y = (v[i].y + deltav.y/2) * DT;
        v[i].x = v[i].x + deltav.x;
        v[i].y = v[i].y + deltav.y;
        p[i].x = p[i].x + deltap.x;
        p[i].y = p[i].y + deltap.y;
        f[i].x = f[i].y = 0.0; // reset force vector
    }
}

int main(int argc, char** argv)
{

    if (argc != 3) {
        cout << "Incorrect amount of arguments!" << endl;
        return 0;
    }

    srand(time(NULL));

    int gnumBodies = atoi(argv[1]);
    int numSteps = atoi(argv[2]);

    cout << "gnumBodies: " << gnumBodies << endl;
    cout << "numSteps: " << numSteps << endl;

    vector<point> p; //Position
    vector<point> v; //velocity
    vector<point> f; //force
    vector<double> m; //mass

    for (int i = 0; i < gnumBodies; i++) {
        p.push_back(point(fRand(0, SIZE), fRand(0, SIZE)));
        v.push_back(point(0,0));
        f.push_back(point(0,0));
        m.push_back(fRand(6e22, 2e30));
    }
    auto start = std::clock();
    for (int i=0; i<numSteps; i++){
        calculateForces(p, f, m, gnumBodies);
        // cout << "hej" << endl;
        moveBodies(p, f, m, v, gnumBodies);
    }

    cout << "Time taken was " << (std::clock() - start) / (double)(CLOCKS_PER_SEC) << " seconds"<< endl;

}