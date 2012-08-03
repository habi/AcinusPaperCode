 #include <iostream>
  using namespace std;

  int main()
  {
    for(int i=0;i<1000;++i)
    {
      cout<<"bla\n";
    }
    //ist aequivalent zu
    //for(;;) //ja, man muss nich ueberall was hinschreiben
        //wenn der mittelteil (die bedingung) leer bleibt
        //dann laeuft die schleife ewig
    //{
      //cout<<"bla\n";
      //++i;
      //if(i>=10) break;
    //}
  }
