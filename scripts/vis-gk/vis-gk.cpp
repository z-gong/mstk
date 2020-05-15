#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

void read_xvg(string xvg, vector<double> &time, vector<double> &temp, vector<double> &pxy, vector<double> &pxz,
              vector<double> &pyz) {
    double   n1, n2, n3, n4, n5;
    ifstream fxvg(xvg);
    string   line;
    while (std::getline(fxvg, line)) {
        if (line.substr(0, 1) == "#" or line.substr(0, 1) == "@")
            continue;
        stringstream ss(line);
        ss >> n1 >> n2 >> n3 >> n4 >> n5;
        time.push_back(n1);
        temp.push_back(n2);
        pxy.push_back(n3);
        pxz.push_back(n4);
        pyz.push_back(n5);
    }
    fxvg.close();
}

void calc_acf(vector<double> &pxy, vector<double> &pxz, vector<double> &pyz, vector<double> &acf, int corr_samples) {
    acf.resize(corr_samples, 0.0);
    for (int i = 0; i < corr_samples; i++) {
        double pxyt = 0, pxzt = 0, pyzt = 0;
        int samples = pxy.size() - i;
        #pragma omp parallel for reduction(+:pxyt,pxzt,pyzt)
        for (int j = 0; j < samples; j++) {
            pxyt += pxy[j] * pxy[i + j];
            pxzt += pxz[j] * pxz[i + j];
            pyzt += pyz[j] * pyz[i + j];
        }
        acf[i]   = (pxyt + pxzt + pyzt) / samples / 3;
    }
}

void integrate(vector<double> &acf, vector<double> &vis, double dt, double convert) {
    vis.resize(acf.size(), 0.0);
    double dvis;
    for (int i = 1; i < acf.size(); i++) {
        dvis = (acf[i] + acf[i - 1]) / 2 * dt * convert;
        vis[i] = vis[i-1] + dvis;
    }
}

double mean(vector<double> &vec) {
    double      t = 0;
    for (double item :vec)
        t += item;
    return t / vec.size();
}

void write_data(string file_acf, string file_vis, vector<double> &time, vector<double> &acf, vector<double> &vis,
                int corr_samples) {
    ofstream facf(file_acf);
    ofstream fvis(file_vis);
    for (int i = 0; i < corr_samples; i++) {
        facf << time[i] << "\t" << acf[i] << "\n";
        fvis << time[i] << "\t" << vis[i] << "\n";
    }
    facf.close();
    fvis.close();
}

int main(int argc, char *argv[]) {
    string fxvg(argv[1]); // "energy.xvg"
    string facf(argv[2]); // "acf.txt"
    string fvis(argv[3]); // "vis.txt"
    double corr_time = atof(argv[4]); // 200.0

    vector<double> time, temp, pxy, pxz, pyz, acf, vis;
    read_xvg(fxvg, time, temp, pxy, pxz, pyz);

    double dt           = time[1] - time[0];
    auto   corr_samples = int(corr_time / dt);
    printf("%s\tsamples: %d\tcorr_samples: %d\n", fxvg.c_str(), time.size(), corr_samples);

    calc_acf(pxy, pxz, pyz, acf, corr_samples);
    double V       = 1;
    double convert = 6.022E-3 / 8.314 * V / mean(temp);
    integrate(acf, vis, dt, convert);

    write_data(facf, fvis, time, acf, vis, corr_samples);
}