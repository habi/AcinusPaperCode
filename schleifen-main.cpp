#include <cstdlib>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    int a;
    int b;
    int erg=1;
    cout<<"Geben sie eine Zahl ein: ";
    cin>>a;
    cout<<"Und jetzt eine andere: ";
    cin>>b;

    int i=0;
    while(i<b)
    {
      i++;
      erg=erg*a;
    }

    cout<<a<<"^"<<b<<"="<<erg<<"\n";
    cin.get();
    system("PAUSE");
    return EXIT_SUCCESS;
}
