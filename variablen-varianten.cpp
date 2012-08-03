    #include<iostream>
    using namespace std;
    int main()
    {
    //Variante 1 - ganz schlecht: variablen mit schlechten namen zu beginn definiert
      int a,b,c;
      /*
        GAAAAAAAAANZ viel Code
      */
      cout<<"gib 2 zahlen ein\n";
      cin>>a;
      cin>>b;
      c=a*b;
      cout<<a<<" mal "<<b<<" ist "<<c<<"\n";

    //Variante 2 - chli besser: wenigstens bessere namen
      int zahl1;      //erste vom user eingegebene Zahl
      int zahl2;      //erste vom user eingegebene Zahl
      int ergebnis;   //produkt von zahl1 und zahl2 (nach der User eingabe)
      /*
        GAAAAAAAAANZ viel Code
      */
      cout<<"gib 2 zahlen ein\n";
      cin>>zahl1;
      cin>>zahl2;
      ergebnis=zahl1*zahl2;
      cout<<zahl1<<" mal "<<zahl2<<" ist "<<ergebnis<<"\n";

    //Variante 3 - am besten - Variablen werden so spät wie möglich definiert und tragen sinnvolle namen
      /*
        GAAAAAAAAANZ viel Code
      */
      cout<<"gib 2 zahlen ein\n";
      int faktor1, faktor2;
      cin>>faktor1;
      cin>>faktor2;
      int produkt=faktor1*faktor2;
      cout<<faktor1<<" mal "<<faktor2<<" ist "<<produkt<<"\n";
    }
