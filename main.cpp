#include <cstdlib>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
      cout<<"Bitte nacheinander 3 Zahlen eingeben (jeweils mit 'Enter' bestaetigen):\nAnschliessend gibt das Programm die kleinste und die groesste Zahl aus.\n";
        int a;
        cin>>a;
        cout<<"...noch eine!\n";
        int b;
        cin>>b;
        cout<<"...und noch eine dritte!\n";
        int c;
        cin>>c;

        int max=a;
    if(b>max)
    {
     max=b; 
    }
    if(c>max)
    {
     max=c;
    }
    
    int min=a;
    if(b<min)
    {
     min=b;
    }
    if(c<min)
    {
     min=c;
    }
    cout<<"die groesste zahl ist "<<max<<" und die kleinste zahlt ist "<<min<<"\n";
    system("PAUSE");
    return EXIT_SUCCESS;
}
