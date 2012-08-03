#include <cstdlib>
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    cout<<"bitte eine zahl kleiner 20 eingeben:\nDas Programm gibt anschliessend aus, ob die Zahl eine Primzahl ist\n";
    int zahl;
    cin>>zahl;
    if(zahl<20)
    {
               cout<<zahl<<" ist kleiner als 20";
               switch(zahl)
               {
                 case 2:
                 case 3:
                 case 5:
                 case 7:
                 case 11:
                 case 13:
                 case 17:
                 case 19:
                      cout<<" und eine Primzahl";
                      break;
                 default: // ergibt alle anderen cases (nicht-prim-zahlen)
                      cout<<" und nicht prim.";
                      break;
               }
    }
//    else if(zahl=20)
//        {
//        cout<<zahl<<" ist genau 20, du loeli";
//        }
    else
        {
        cout<<zahl<<" ist groesser als 20, du loeli";
        }                             
     cout<<"\n";
     system("PAUSE");
     return EXIT_SUCCESS;
}
