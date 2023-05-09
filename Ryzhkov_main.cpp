#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <complex>

using namespace std;

double MeanEnergy(const vector<complex<double>> &v){
    double sum = 0.0;
    for(int i = 0; i < v.size(); i++){
        sum += norm(v[i]);
    }
    return sum / v.size();
}

double PAPR(const vector<complex<double>> &v){
    double mean_energy = MeanEnergy(v);
    
    double max_norm = 0.0;
    double temp;
    for(int i = 0; i < v.size(); i++){
        temp = norm(v[i]);
        if (temp > max_norm)
            max_norm = temp;
    }
    return 10. * log10(max_norm / mean_energy);
}

int main()
{
    vector<complex<double>> qam_16(16);
    vector<complex<double>> aqam_16(16);
    
    double a = 1.0;
    complex<double> point;
    double rho, theta;
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            point = {-3. + 2 * i, -3. + 2 * j};
            qam_16[4 * i + j] = a * point;
            
            rho = i ? sqrt(2. + 6 * i) : 2.;
            theta = (i % 2) * M_PI / 4 + j * M_PI / 2;
            aqam_16[4 * i + j] = a * polar(rho, theta);
        }
        
    }
    cout << "QAM-16 constellation: ";
    for (auto p: qam_16)
        cout << p << ' ';
    cout << endl;
    
    cout << "AQAM-16 constellation: ";
    for (auto p: aqam_16)
        cout << p << ' ';
    cout << endl;
        
    cout << "Mean energy for QAM-16: " << MeanEnergy(qam_16) << endl;
    cout << "Mean energy for AQAM-16: " << MeanEnergy(aqam_16) << endl;
    
    cout << "Peak-to-average power ratio for QAM-16: " << PAPR(qam_16) << endl;
    cout << "Peak-to-average power ratio for AQAM-16: " << PAPR(aqam_16) << endl;

    return 0;
}

