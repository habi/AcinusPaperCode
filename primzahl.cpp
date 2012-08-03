 #include <iostream>
    using namespace std;

    int main()
    {
    cout<<"Bitte Zahl eingeben: Programm gibt\nalle Primzahlen bis zu dieser Zahl aus: ";
    int max;
    cin>>max;
     for(int i=3; i<max; ++i)
      {
        bool isprim=true;
        for(int j=2; j<i; ++j)
        {
          if(i%j == 0) { isprim=false; break; }
        }
        if(isprim) cout<<i<<" ist eine Primzahl\n";
      }
    }
    
