 #include <iostream>
    using namespace std;

    int main()
    {
    cout<<"Schleife 1 - for: ";
                    for(int i=1; i<11; i++)
                    {
                     cout<<i<<" - ";
                    }
    cout<<"\nSchleife 2 - while: ";
                    int j=0;
                    while(j++<10)
                    {
                     cout<<j<<" - ";
                    }
    cout<<"\nSchleife 3 - do while: ";
                    int k=0;
                    do 
                    {
                     cout<<k+1<<" - ";
                     k++;
                    }                     
                    while(k<10);
                    cin.get();
    }
    
    
