#include<iostream>
#include<vector>
#include<stdio.h>
#include <iomanip>
#include<fstream>
using namespace std;

class oneD {
    private:
        long double E= 200000;//MPa
        long double A= 300;//mm2
        long double L=800;//mm
        int Ne;
        long double P0= 15;//N/mm
        vector<long double> F ;
        vector<long double> X;
    public:
        oneD();
        vector<long double> matrixSolver();
        void Display(vector<long double> & arr);
        void formula(long double y,long double x1, long double x2, long double u1, long double u2);
        void UBinary(long double y,vector<long double> & u);
        void SBinary(long double y,vector<long double> & u);
        void plotGraph(vector<long double> & u);
};
oneD::oneD() {
    cout<<"Enter number of elements: ";
    cin>>Ne;
    if(Ne<=0) {
        cout<<"Please enter valid no: of elements";
        exit(EXIT_FAILURE);
    }
    F.resize(Ne+1, 0.0);
    X.resize(Ne+1, 0.0);
    long double step = L / Ne;
    for (int i = 0; i <= Ne; i++) {
        X[i] = i * step; // nodal coordinate vector banaya
    }
    long double F0= (P0/(2*Ne))*L , F1= 2*F0;
    for(int i=0;i<F.size();i++) {
        if(i==0 || (i== F.size()-1)) F[i]=F0;
        else F[i]=F1;
    }
}
vector<long double> oneD::matrixSolver() {
    long double Le=L/Ne,P= (E/Le)*A, Q= 2*P, R= -P;
    int Nn= Ne+1;
    vector<long double> l(Nn - 2, 0.0),u(Nn-1, 0.0),y(Nn-1, 0.0),x(Nn-1, 0.0);

    //initialize karo l and u//
    u[0]=Q;
    for(int i=1;i<Nn-1;i++) {
        l[i - 1] = R / u[i - 1];
        if(i>0 && i<Nn-2) u[i] = Q - l[i - 1] * R;
        else u[i] = P - l[i-1] * R;
    }

    //y nikala//
    y[0] = F[1];
    for (int i = 1; i < Nn-1; i++) {
        y[i] = F[i+1] - l[i - 1] * y[i - 1];
    }

    //solution vector//
    x[Nn - 2] = y[Nn - 2] / u[Nn - 2];
    for (int i = Nn - 3; i >= 0; i--) {
        x[i] = (y[i] - R * x[i + 1]) / u[i];
    }
    return x;
}

void oneD::Display(vector<long double> & arr) {
    cout << "u1"<<": "<<0.00<< endl;
    for (int i = 0; i < arr.size(); i++) {
        cout<<"u"<<i+2<<": "<<arr[i]<<"mm"<<endl;
    }
}

void oneD::formula(long double y,long double x1, long double x2, long double u1, long double u2) {
    cout<<"Displacement at x= "<<y<<"mm: "<<setprecision(15)<<((x2-y)/(x2-x1))*u1 + ((y-x1)/(x2-x1))*u2<<"mm"<<endl;
}

void oneD::UBinary(long double y,vector<long double> & u) {
    int l=0,r=X.size()-1,ans=0;
    while(l<=r) {
        int mid = (l+r)/2;
        if(X[mid]>= y) {
            ans=mid;
            r=mid-1;
        }
        else {
            l=mid+1;
        }
    }
    if(ans==1) return formula(y,X[ans-1],X[ans],0,u[ans-1]);
    else if(ans==0) return formula(0,X[ans+1],X[ans],0,0);
    return formula(y,X[ans-1],X[ans],u[ans-2],u[ans-1]);
}

void oneD::SBinary(long double y,vector<long double> & u) {
    int l=0,r=X.size()-1,ans=0;
    long double strain=0;
    while(l<=r) {
        int mid = (l+r)/2;
        if(X[mid]>= y) {
            ans=mid;
            r=mid-1;
        }
        else {
            l=mid+1;
        }
    }
    if(ans==1 || ans==0)  strain = (u[0]*Ne)/L ;
    else strain = ((u[ans-1]-u[ans-2])*Ne)/L;
    cout<<"Strain at x= "<<y<<"mm: "<<setprecision(15)<<strain<<endl;
    cout<<"Stress at x= "<<y<<"mm: "<<setprecision(15)<<E*strain<<"MPa"<<endl;
}

void oneD::plotGraph(vector<long double> & u) {
    ofstream uFile("uFile.dat");
    ofstream sigmaFile("sigmaFile.dat");

    if (!uFile || !sigmaFile) {
        cout << "Error creating data files!\n";
        return;
    }

    for (int i = 0; i < u.size(); i++) {
        long double strain = (i == 0) ? (u[0] * Ne) / L : ((u[i] - u[i - 1]) * Ne) / L;
        long double stress = E * strain;
        uFile << X[i + 1] << " " << u[i] << endl;
        sigmaFile << X[i + 1] << " " << stress << endl;
    }

    uFile.close();
    sigmaFile.close();

    // Call GNUPlot
    system("gnuplot -p -e \"plot 'uFile.dat' using 1:2 with lines title 'x vs u(x)'\"");
    system("gnuplot -p -e \"plot 'sigmaFile.dat' using 1:2 with lines title 'x vs stress'\"");
}

int main() {
     oneD K;
    vector<long double> r;
    int choice;
    do {
        cout << "\nMenu:\n";
        cout << "1. Solve for nodal displacements\n";
        cout << "2. Display global displacement vector\n";
        cout << "3. Compute displacement,strain and stress at a position\n";
        cout << "4. Plot graphs\n";
        cout << "5. Exit\n";
        cout << "Enter your choice: ";
        cin >> choice;
        switch(choice) {
            case 1:
                r = K.matrixSolver();
                cout << "Solution computed!\n";
                break;
            case 2:
                if(r.empty()) cout << "Please,solve first!!\n";
                else K.Display(r);
                break;
            case 3:
                if(r.empty()) cout << "Solve first!\n";
                else {
                    long double x;
                    cout << "Enter position: ";
                    cin >> x;
                    K.UBinary(x, r);
                    K.SBinary(x, r);
                }
                break;
            case 4:
                if (r.empty()) cout << "Solve first!\n";
                else K.plotGraph(r);
                break;
            case 5:
                cout << "Exiting...\n";
                break;
            default:
                cout << "Please enter any number from the choices available\n";
        }
    } while(choice != 5);
    return 0;
}