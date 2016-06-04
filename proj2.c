/*
 * Soubor:  proj2.c
 * Datum:   20.11.2011
 * Autor:   Radim Kubis, xkubis03@stud.fit.vutbr.cz
 * Projekt: Iteracni vypocty, projekt c. 2 pro predmet IZP
 * Popis:   Program pocita: arkussinus, logaritmus o zadanem zakladu,
 *                          delku lomene cary a delku lomene cary s chybou.
 *  
 *          Na standardnim vstupu ocekava jedno cislo v pripade arkussinu
 *          a logaritmu, dve cisla (souradnice) v pripade delky lomene cary,
 *          at uz s chybou, nebo bez.
 *          Vystupem je hodnota funkce arkussinu, logaritmu nebo delky lomene
 *          cary bez chyby. V pripade delky lomene cary s chybou jsou vystupni
 *          hodnoty dve.     
 *           
 *          Pro vyber funkce, ktera bude v programu spustena, zadejte
 *          na prikazovem radku: 
 *           
 *          -h
 *          
 *            vypis napovedy programu    
 *  
 * 
 *          --arcsin sigdig
 *           
 *            funkce arkussinus, sigdig je cislo predstavujici presnost vysledku
 *            na pocet platnych cislic
 *              
 * 
 *          --logax sigdig a
 *          
 *            funkce logaritmus, sigdig je cislo predstavujici presnost vysledku
 *            na pocet platnych cislic, cislo a je zaklad logaritmu
 *            
 * 
 *          --lbl
 *          
 *            funkce pro vypocet delky lomene cary bez chyby
 *            
 * 
 *          --lble ERR
 *          
 *            funkce pro vypocet delky lomene cary s chybou, parametr ERR
 *            je velikost chyby
 *            
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>

const double IZP_LN_10 = 2.302585092994045684018;  // Prirozeny logaritmus 10
const double IZP_PI_2  = 1.57079632679489661923;   // pi/2
const double IZP_PI_6  = 0.52359877559829887308;   // pi/6
const double IZP_SR_3  = 1.73205080756887729353;   // sqrt(3)

// Struktura pro uchovani souradnic bodu
typedef struct bod {
  double x;  // Souradnice X
  double y;  // Souradnice Y
} point;

// Struktura pro uchovani max. a min. hodnoty intervalu
typedef struct interval {
  point min;  // Minima
  point max;  // Maxima
} range;

// Struktura pro uchovani hodnot vypoctu delky lomene cary bez chyby
typedef struct longNormal {
  int prvni;     // Byl zadany bod prvni?
  point last;    // Souradnice predchoziho bodu
  double delka;  // Celkova delka lomene cary
} delkaN;

// Struktura pro uchovani hodnot vypoctu delky lomene cary s chybou
typedef struct longError {
  int prvni;        // Byl zadany bod prvni?
  point bod0;       // Souradnice predchoziho bodu
  double error;     // Vzdalenost chyby od bodu
  range aktualni;   // Aktualni rozsah hodnot
  point interval;   // Vysledny interval
  range predchozi;  // Predchozi rozsah hodnot
} delkaE;

// Seznam moznych volanych funkci
enum funkce {
  LBL = 1,  // Delka lomene cary
  LBLE,     // Delka lomene cary bez chyby
  ARCSIN,   // Arkus sinus
  LOGAX     // Logaritmus o zakladu A z X
};

// Seznam moznych druhu vypoctu pro arkus tangens
enum typVypoctu {
  NORMALNI = 0, // Cislo je mensi nez (2.0 - sgrt(3)) = cca 0.268
  PI2_MINUS,    // Cislo je vetsi nebo rovno 1.0
  PI6_PLUS      // Cislo je vetsi nez (2.0 - sqrt(3)) = cca 0.268
};

/**
 * Funkce inicializujici bod
 * @param p - bod pro inicializaci
 */
void initPoint(point *p) {
  p->x = 0.0;
  p->y = 0.0;
}

/**
 * Funkce inicializujici rozsah
 * @param r - rozsah pro inicializaci
 */
void initRange(range *r) {
  initPoint(&(r->min));
  initPoint(&(r->max));
}

/**
 * Funkce inicializujici strukturu pro vypocet delky bez chyby
 * @param r - rozsah pro inicializaci
 */
void initDelkaN(delkaN *l) {
  initPoint(&(l->last));
  l->prvni = 1;
  l->delka = 0.0;
}

/**
 * Funkce inicializujici strukturu pro vypocet delky s chybou
 * @param l - struktura pro inicializaci
 */
void initDelkaE(delkaE *l) {
  initRange(&(l->predchozi));
  initRange(&(l->aktualni));
  initPoint(&(l->interval));
  initPoint(&(l->bod0));
  l->prvni = 1;
  l->error = -1.0;
}

/**
 * Funkce vraci nejvetsi hodnotu ze tri hodnot
 * @param a - vstupni cislo
 * @param b - vstupni cislo
 * @param c - vstupni cislo
 * @return - vraci nejvetsi cislo
 */
double max(double a, double b, double c) {
  return (a < b) ? ((b < c) ? c : b) : ((a < c) ? c : a);
}

/**
 * Funkce vraci nejmensi hodnotu ze tri hodnot
 * @param a - vstupni cislo
 * @param b - vstupni cislo
 * @param c - vstupni cislo
 * @return - vraci nejmensi cislo
 */
double min(double a, double b, double c) {
  return (a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c);
}

/**
 * Funkce pocita absolutni presnost z presnosti relativni (pocet platnych cifer)
 * @param cislo - cislo, ze ktereho se pocita epsilon
 * @param sigdig - pozadovana presnost vysledku na pocet platnych cislic
 * @return - vraci vypocitane epsilon
 */
double dejEpsilon(double cislo, int sigdig) {

  double mocnitel = 0.0;  // Mocnitel cisla deset pro epsilon

  cislo = fabs(cislo);    // Pro vypocet bude lepsi absolutni hodnota 

  // Dokud je cislo mensi jak jedna, nenarazil jsem na platnou cislici
  while (cislo < 1.0) {
    mocnitel -= 1.0;
    cislo    *= 10.0;
  }

  // Vracim vypocitane epsilon
  return pow(10.0, (mocnitel - sigdig));
}

/**
 * Funkce pocita arkus tangens zadaneho cisla
 * na presnost danou druhym parametrem, ktera vyjadruje
 * pocet platnych cislic vysledku. 
 * @param cislo - cislo, ze ktereho se pocita arkus tangens
 * @param sigdig - presnost vysledku na pocet platnych cislic
 * @return - vraci vypocitany arkus tangens
 */
double arctg(double cislo, int sigdig) {

  // Promenne pro vypocet iterace
  double korekce = 1.0;     // Korekce vysledku podle zapornosti/kladnosti cisla
  double mocnina = 0.0;     // Mocnina cisla pouzita v iteraci
  double vysledek = 0.0;    // Aktualni vysledek
  double mocnina2 = 0.0;    // Druha mocnina cisla
  double epsilon = -1.0;    // Absolutni presnost epsilon, prozatim nenastaveno
  double predchozi = 0.0;   // Ulozeny predchozi vysledek iterace
  int vypocet = NORMALNI;   // Typ provadeneho vypoctu
  double znaminko = -1.0;   // V rade se jednou odecita, podruhe pricita
  double jmenovatel = 3.0;  // Jmenovatel v zlomku vypoctu
  
  // Pokud je vstupni cislo nula
  if (cislo == 0.0) {
    return 0.0;  // vracim nulu
  }
  
  // Poznamenam si znamenko vstupniho cisla a pracuji s kladnym
  if (cislo < 0.0) {
    cislo   = -cislo;
    korekce = -korekce;
  }

  // Podle cisla rozhodnu, ktery vypocet provedu
  // a popr. upravim vstupni hodnotu
  // Hodnoty intervalu jsou zvoleny podle algoritmu vypoctu
  if (cislo >= 1.0) {
    cislo = (1.0 / cislo);
    vypocet = PI2_MINUS;
  } else if (cislo > (2.0 - IZP_SR_3)) {  // > cca 0.268
    cislo = ((IZP_SR_3 * cislo - 1.0) / (IZP_SR_3 + cislo));
    vypocet = PI6_PLUS;
  }
  
  // Inicializuji promenne pro vypocet
  vysledek = cislo;
  mocnina2 = cislo * cislo;
  mocnina = mocnina2 * cislo;
  
  do {
    predchozi = vysledek;  // Ulozim si predchozi vysledek
    vysledek += (znaminko * (mocnina / jmenovatel));  // Provedu krok iterace
    znaminko *= (-1.0);  // Zmenim znamenko na opacne pro dalsi iteraci
    jmenovatel += 2.0;  // Zvysim jmenovatele o 2
    mocnina *= mocnina2;  // Vypocitam hodnotu dalsi mocniny zlomku

    // Pokud jsem jeste nevypocital epsilon, provedu tak
    if (epsilon == -1.0) {
      // Je velice mala pravdepodobnost, ze by se zmenil rad vysledku
      epsilon = dejEpsilon(vysledek, sigdig);
    }
  } while (fabs(vysledek - predchozi) >= epsilon);

  // Podle druhu vypoctu vracim vysledek
  if (vypocet == PI2_MINUS) {
    return ((IZP_PI_2 - vysledek) * korekce);
  } else if (vypocet == PI6_PLUS) {
    return ((IZP_PI_6 + vysledek) * korekce);
  }
  return ((vysledek) * korekce);
}  // Konec funkce arctg

/**
 * Funkce pocita arkus sinus zadaneho cisla
 * na presnost danou druhym parametrem, ktera vyjadruje
 * pocet platnych cislic vysledku.
 * Pro urcite intervaly je vyuzit vztah s funkci arkus tangens
 * @param cislo - cislo, ze ktereho se pocita arkus sinus
 * @param sigdig - presnost vysledku na pocet platnych cislic
 * @return - vraci vypocitany arkus sinus
 */
double arcsin(double cislo, int sigdig) {

  double korekce = 1.0;     // Korekce vysledku podle zapornosti/kladnosti cisla
  double mocnina = 0.0;     // Mocnina cisla pouzita v iteraci
  double zlomek = 0.5;      // Pocatecni hodnota zlomku pro nasobeni v iteraci
  double vysledek = 0.0;    // Aktualni vysledek
  double mocnina2 = 0.0;    // Druha mocnina cisla
  double epsilon = -1.0;    // Absolutni presnost epsilon, prozatim nenastaveno
  double predchozi = 0.0;   // Ulozeny predchozi vysledek iterace
  double jmenovatel = 1.0;  // Jmenovatel v zlomku vypoctu
  
  if (cislo == 0.0) {
    return 0.0;  // Pokud je cislo nula, vracim nulu
  } else if (fabs(cislo) == 1.0) {
    return IZP_PI_2 * cislo;  // Pokud je absolutni hodnota cisla 1, vracim pi/2
  } else if (fabs(cislo) > 1.0) {
    return NAN;  // Pokud je absolutni hodnota cisla vetsi jak 1, vracim NaN
  }
  
  // Pokud se jedna o zaporne cislo, poznamenam si znamenko a pracuji s kladnym
  if (cislo < 0.0) {
    korekce = -korekce;
    cislo = fabs(cislo);
  }

  // Pokud je cislo 1/2
  if(cislo == 0.5) {
    return IZP_PI_6*korekce;  // Vracim pi/6
  }

  // Kdyz je cislo vetsi nez 0.71, pouziji vztah s arkus tangens
  // Hodnota 0.71 byla zvolena jako nejvhodnejsi hranice po otestovani
  if (cislo > 0.71) {
    // Prepocet cisla pro arkus tangens
    cislo = (1.0 / (cislo / sqrt(1.0 - cislo * cislo)));
    
    // Pokud bylo zaporne, vracim zaporny vysledek
    if (korekce < 0.0) {
      return ((IZP_PI_2 - arctg(cislo, sigdig)) * korekce);
    }
    // Jinak vracim kladny
    return (IZP_PI_2 - (arctg(cislo, sigdig) * korekce));
  }
  
  // Pokud je cislo mensi nebo rovno 0.71 a nenastal zvlastni pripad,
  // pocitam iteraci pro arkus sinus
  vysledek = cislo;  // Inicializace vysledku na zadane cislo
  mocnina2 = cislo * cislo;  // Druha mocnina cisla, pro iteraci
  mocnina  = vysledek * mocnina2;  // Treti mocnina, pro iteraci

  // Hlavni smycka iterace
  do {
    predchozi = vysledek;  // Ulozim si predchozi vysledek
    
    vysledek += (zlomek * (mocnina / (jmenovatel + 2)));  // Krok iterace
    
    jmenovatel += 2.0;  // Zvysim jmenovatel o 2
    
    // Prepocitam pomocne vypocty iterace
    zlomek  *= (jmenovatel / (jmenovatel + 1));
    mocnina *= mocnina2;
    
    // Pokud jsem jeste nevypocital epsilon, provedu tak
    if (epsilon == -1.0) {
      // Je velice mala pravdepodobnost, ze by se zmenil rad vysledku
      epsilon = dejEpsilon(vysledek, sigdig);
    }
  // Pokud jsem dosahl zadane presnosti, koncim cyklus
  } while (fabs(vysledek - predchozi) >= epsilon);

  // Vracim vysledek arkus sinus
  return (vysledek * korekce);
}  // Konec funkce arcsin

/**
 * Funkce pocita logaritmus cisla o zadanem zakladu
 * na presnost danou tretim parametrem, ktera vyjadruje
 * pocet platnych cislic vysledku.
 * Vysledek je pocitan jako podil dvou prirozenych logaritmu. 
 * @param cislo - cislo, ze ktereho se pocita logaritmus
 * @param zaklad - udava zaklad logaritmu 
 * @param sigdig - presnost vysledku na pocet platnych cislic
 * @return - vraci vypocitany logaritmus
 */
double logax(double cislo, double zaklad, int sigdig)
{
  // Spolecne promenne pro oba logaritmy
  int delitel = 3;
  double epsilon = 0.0;    // Absolutni presnost vysledku vypocitana z relativni
  double vysledek = 0.0;   // Vysledek funkce
  double predchozi = 0.0;  // Prvek i-1 algoritmu

  // Promenne pro prirozeny logaritmus citatele
  double citatel = 0.0;    // Prubezny vysledek
  double pocatekC = 0.0;   // Prvni hodnota iterace
  double nasobekC = 0.0;   // Druha mocnina prvni hodnoty
  double mocninaC = 0.0;   // Mocnina vysledku v iteraci
  double exponentC = 0.0;  // Exponent pro prepocet logaritmu
  
  // Promenne pro prirozeny logaritmus jmenovatele
  double pocatekJ = 0.0;    // Prvni hodnota iterace
  double nasobekJ = 0.0;    // Druha mocnina prvni hodnoty
  double mocninaJ = 0.0;    // Mocnina vysledku v iteraci
  double exponentJ = 0.0;   // Exponent pro prepocet logaritmu
  double jmenovatel = 0.0;  // Prubezny vysledek
  
  // Pokud nastane nektera kombinace ze zvlastnich pripadu,
  // neprovadim vypocet, ale vracim urcitou hodnotu
  if (zaklad == 0.0) {
    if (cislo <= 0.0) {
      return NAN;
    } else if (cislo >= 1.0) {
      return -0.0;
    } else {
      return 0.0;
    }
  } else if (zaklad == 1.0) {
    if (cislo == 1.0 || cislo < 0.0) {
      return NAN;
    } else if (cislo > 1.0) {
      return INFINITY;
    } else {
      return -INFINITY;
    }
  } else if (zaklad < 0.0) {
    return NAN;
  } else if (zaklad < 1.0) {
    if(cislo == 1.0) {
      return -0.0;
    } else if (cislo == 0) {
      return INFINITY;
    }
  }
  
  if (cislo < 0.0) {
    return NAN;
  } else if(cislo == 0.0) {
    return -INFINITY;
  } else if(cislo == INFINITY) {
    return INFINITY;
  }

  // Prepocet cisla do vhodnejsiho intervalu
  // a vypocet exponentu pro vysledek
  // Hodnoty intervalu zvoleny jako nejlepsi mozne po otestovani
  if (cislo > 4.7) {
    while (cislo > 1.0) {
      cislo = cislo * 0.1;
      exponentC = exponentC + 1.0;
    }
  } else if (cislo < 0.25) {
    while ((cislo * 10) < 4.7) {
      cislo = cislo * 10;
      exponentC = exponentC - 1.0;
    }
  }
  
  // Prepocet zakladu do vhodnejsiho intervalu
  // a vypocet exponentu pro vysledek
  // Hodnoty intervalu zvoleny jako nejlepsi mozne po otestovani
  if (zaklad > 4.7) {
    while (zaklad > 1.0) {
      zaklad = zaklad * 0.1;
      exponentJ = exponentJ + 1.0;
    }
  } else if (zaklad < 0.25) {
    while ((zaklad * 10) < 4.7) {
      zaklad = zaklad * 10;
      exponentJ = exponentJ - 1.0;
    }
  }

  // Pokud nenastal nektery zvlastni pripad hodnot,
  // provedu vypocet logaritmu
  
  // Inicializace hodnot pro vypocet logaritmu citatele
  pocatekC = (cislo - 1) / (cislo + 1);  // Prvni clen rady
  citatel = 2 * pocatekC;  // Vysledek se nasobi 2
  nasobekC = pocatekC * pocatekC;  // Druha mocnina
  mocninaC = nasobekC * pocatekC;  // Treti mocnina
  
  // Inicializace hodnot pro vypocet logaritmu jmenovatele
  pocatekJ = (zaklad - 1) / (zaklad + 1);  // Prvni clen rady
  jmenovatel = 2 * pocatekJ;  // Vysledek se nasobi 2
  nasobekJ = pocatekJ * pocatekJ;  // Druha mocnina
  mocninaJ = nasobekJ * pocatekJ;  // Treti mocnina

  // Hlavni smycka vypoctu, ktera konci s dosazenim presnosti
  do
  {
    predchozi = vysledek;  // Ulozim si predchozi vysledek iterace

    // Spocitam novy krok iterace citatele
    citatel += 2 * (mocninaC / delitel);
    mocninaC *= nasobekC;
    
    // Spocitam novy krok iterace jmenovatele
    jmenovatel += 2 * (mocninaJ / delitel);
    mocninaJ *= nasobekJ;
    
    delitel += 2;  // Delitel se zvysuje o 2

    // Vypocet vysledku po aktualni iteraci
    vysledek = (citatel + (exponentC * IZP_LN_10)) / (jmenovatel + (exponentJ * IZP_LN_10));
    
    // Zjistim rad vysledku pro overeni presnosti
    epsilon = dejEpsilon(vysledek, sigdig);

  // Pokud jsem nedosahl pozadovane presnosti, pokracuji
  } while (fabs(predchozi - vysledek) >= epsilon);

  // Vysledek
  return vysledek;
}  // Konec funkce logax

/**
 * Funkce vraci celkovou prubeznou vzdalenost lomene cary
 * @param bod1 - prave nacteny bod
 * @param st - struktura s hodnotami vypoctu delky
 * @return - vraci celkovou prubeznou vzdalenost
 */
double lomenaCara(point *bod1, delkaN *st) {

  // Vypocitam si rozdil dvou po sobe jdoucich souradnic
  st->last.x = bod1->x - st->last.x;
  st->last.y = bod1->y - st->last.y;
  
  // Pokud byl bod prvni
  if (st->prvni == 1) {
    st->prvni = 0;  // Poznacim si a delku nepocitam
  } else {  // Jinak delku pocitam
    st->delka += sqrt(st->last.x * st->last.x + st->last.y * st->last.y);
  }
  
  // Ulozim si aktualne nacteny bod jako posledni
  st->last = *bod1;
  
  // Vracim delku lomene cary
  return st->delka;
}

/**
 * Funkce pocita celkovou prubeznou vzdalenost lomene cary s chybou
 * @param chyba - velikost chyby
 * @param bod - struktura s hodnotami prave nacteneho bodu
 * @param st - struktura s hodnotami pro vypocet
 */
void lomenaCaraChyba(double chyba, point bod, delkaE *st) {

  range vypocet;  // Pomocny rozsah pro vypocet
  double minimum = -1.0;  // Aktualni minimalni vzdalenost
  double vzdalenostBodu = 0.0;  // Vzdalenost poslednich dvou bodu
  double soucin1, soucin2, soucin3;  // Pomocne souciny pro intervalove pocty
  
  // Pokud jeste nemam inicializovanou vzdalenost chyby, inicializuji
  if (st->error == -1.0) {
    st->error = 2 * sqrt((chyba * chyba) + (chyba * chyba));
  }
  
  // Vypocitam vzdalenost poslednich dvou bodu
  vzdalenostBodu = sqrt(
                     (st->bod0.x - bod.x) * (st->bod0.x - bod.x) +
                     (st->bod0.y - bod.y) * (st->bod0.y - bod.y)
                    );
  
  // Vypocitam aktualni max. a min. souradnice chyb
  st->aktualni.min.x = bod.x - chyba;
  st->aktualni.max.x = bod.x + chyba;
  st->aktualni.min.y = bod.y - chyba;
  st->aktualni.max.y = bod.y + chyba;
  
  if(st->prvni == 0) {  // Pokud se nejedna o prvni bod
    if (st->bod0.x == bod.x && st->bod0.y == bod.y) {  // Pokud jsou body stejne
      return;  // Nic nepocitam
    }
    if (fabs(st->bod0.x) == fabs(bod.x) || fabs(st->bod0.y) == fabs(bod.y)) {
      // Pokud se s chybou dostanu na stejne souradnice na nektere z os
      vzdalenostBodu -= 2.0 * chyba;  // Spocitam vzdalenost hranic chyby
      if (vzdalenostBodu > 0.0) {  // Pokud se neprekryvaji
        minimum = vzdalenostBodu;  // Ukladam minimalni vzdalenost
      } else {
        minimum = 0.0;  // Pokud se prekryvaji, ukladam vzdalenost nulovou
      }
    } else if(st->error < vzdalenostBodu) {
      // Pokud je vzdalenost bodu vetsi, nez velikost chyby
      if (
          (st->bod0.x + chyba) == (bod.x - chyba) ||
          (bod.x + chyba) == (st->bod0.x - chyba)
         ) {
         // Pokud se s chybou setkaji souradnice na ose X,
         // bude vzdalenost rozdil souradnic Y
        minimum = fabs(fabs(fabs(st->bod0.y) - chyba) - fabs(fabs(bod.y) - chyba));
      } else if (
                 (st->bod0.y + chyba) == (bod.y - chyba) ||
                 (bod.y + chyba) == (st->bod0.y - chyba)
                ) {
        // Pokud se s chybou setkaji souradnice na ose Y,
         // bude vzdalenost rozdil souradnic X
        minimum = fabs(fabs(fabs(st->bod0.x) - chyba) - fabs(fabs(bod.x) - chyba));
      }
    } else {
      // Pokud se s chybou souradnice nesetkaji na zadne z os
      if (vzdalenostBodu - st->error > 0) {  // Pokud je rozdil nezaporny
        minimum = vzdalenostBodu - st->error;  // Bude to nejmensi hodnota
      } else {  // Jinak je nejmensi hodnota nulova
        minimum = 0.0;
      }
    }
    
    // Vypocitam si rozdil vzdalenosti intervalu
    vypocet.min.x = st->aktualni.min.x - st->predchozi.max.x;
    vypocet.max.x = st->aktualni.max.x - st->predchozi.min.x;
    vypocet.min.y = st->aktualni.min.y - st->predchozi.max.y;
    vypocet.max.y = st->aktualni.max.y - st->predchozi.min.y;
    
    // Provedu souciny pro X souradnice
    soucin1 = vypocet.min.x * vypocet.min.x;
    soucin2 = vypocet.min.x * vypocet.max.x;
    soucin3 = vypocet.max.x * vypocet.max.x;
    
    // Pokud je soucin2 mensi jak 0 - vzdalenost nesmi byt zaporna
    if (soucin2 < 0) {
      soucin2 = 0.0;
    }
    
    // Vyberu nejmensi a nejvetsi vzdalenost 
    vypocet.min.x = min(soucin1, soucin2, soucin3);
    vypocet.max.x = max(soucin1, soucin2, soucin3);
    
    // Vypocitam souciny pro Y souradnice
    soucin1 = vypocet.min.y * vypocet.min.y;
    soucin2 = vypocet.min.y * vypocet.max.y;
    soucin3 = vypocet.max.y * vypocet.max.y;
    
    // Pokud je soucin2 mensi jak 0 - vzdalenost nesmi byt zaporna
    if (soucin2 < 0) {
      soucin2 = 0.0;
    }
    
    // Vyberu nejmensi a nejvetsi vzdalenost
    vypocet.min.y = min(soucin1, soucin2, soucin3);
    vypocet.max.y = max(soucin1, soucin2, soucin3);
    
    // Pokud jsem nevypocital minimalni vzdalenost jiz drive
    if (minimum == -1.0) {
      // Bude nejmensi vzdalenost brana z intervalovych poctu
      st->interval.x += sqrt(fabs(vypocet.min.x + vypocet.min.y));
    } else {
      // Jinak vyberu nejmensi vzdalenost ze specialnich pripadu
      st->interval.x += minimum;
    }
    // Nejvetsi vzdalenost pocitam pomoci intervalovych poctu
    st->interval.y += sqrt(fabs(vypocet.max.x + vypocet.max.y));
  }
  
  // Ulozim si aktualni rozsahy jako predchozi
  st->predchozi.min.x = st->aktualni.min.x;
  st->predchozi.min.y = st->aktualni.min.y;
  st->predchozi.max.x = st->aktualni.max.x;
  st->predchozi.max.y = st->aktualni.max.y;
  
  // Pokud byl bod prvni
  if(st->prvni == 1) {
    st->prvni = 0;  // Nastavim priznak, ze dalsi uz prvni nebude
    st->interval.x = 0.0;  // Vzdalenost je zatim nulova
    st->interval.y = 0.0;  // Vzdalenost je zatim nulova
  }
  
  // Ulozim si souradnice aktualniho bodu
  st->bod0.x = bod.x;
  st->bod0.y = bod.y;
}

int main(int argc, char* argv[]) {
  
  point bod;    // Bod, do ktereho nacitam souradnice ze vstupu
  delkaE stE;   // Struktura pro vypocet delky lomene cary s chybou
  delkaN stN;   // Struktura pro vypocet delky lomene cary bez chyby
  int funkce = 0;  // Promenna urcuje spoustenou funkce podle paramatru
  int result = 0;  // Vysledek operace nacteni cisla
  int sigdig = 0;  // Presnost v poctu platnych cislic
  double chyba = 0.0;  // Uzivatelem zadana chyba pro lomenou caru
  double cislo = 0.0;  // Cislo zadane uzivatelem
  double zaklad = 0.0; // Zaklad logaritmu
  double vysledek = 0.0;  // Vysledek funkce, mimo delky lomene cary s chybou
  
  // Pokud byl spusten program bez parametru
  if(argc < 2) {
    fprintf(stderr, "Nebyl zadan zadny parametr funkce.\n");
    return EXIT_FAILURE;
  } else {
    if(argc == 2) {  // Pokud ma program jen jeden parametr
      if(strcmp(argv[1], "-h") == 0) {  // Jedna se o napovedu?
        fprintf(stdout, "  Iteracni vypocty\n"
                        "  Autor: Radim Kubis, xkubis03@stud.fit.vutbr.cz\n\n"
                        "  Popis: Program pocita: arkussinus, logaritmus\n"
                        "  o zadanem zakladu, delku lomene cary a delku lomene\n"
                        "  cary s chybou.\n\n"
                        "  Na standardnim vstupu ocekava jedno cislo v pripade\n"
                        "  arkussinu a logaritmu, dve cisla (souradnice)\n"
                        "  v pripade delky lomene cary, at uz s chybou,\n"
                        "  nebo bez.\n\n"
                        "  Vystupem je hodnota funkce arkussinu, logaritmu nebo\n"
                        "  delky lomene cary bez chyby. V pripade delky lomene\n"
                        "  cary s chybou jsou vystupni hodnoty dve.\n\n"
                        "  Pro vyber funkce, ktera bude v programu spustena,\n"
                        "  zadejte na prikazovem radku:\n\n"
                        "    -h\n"
                        "      vypis napovedy programu\n\n"
                        "    --arcsin sigdig\n"
                        "      funkce arkussinus, sigdig je cislo predstavujici\n"
                        "      presnost vysledku na pocet platnych cislic\n\n"
                        "    --logax sigdig a\n"
                        "      funkce logaritmus, sigdig je cislo predstavujici\n"
                        "      presnost vysledku na pocet platnych cislic,\n"
                        "      cislo a je zaklad logaritmu\n\n"
                        "    --lbl\n"
                        "      funkce pro vypocet delky lomene cary bez chyby\n\n"
                        "    --lble ERR\n"
                        "      funkce pro vypocet delky lomene cary s chybou,\n"
                        "      parametr ERR je velikost chyby\n\n");
        return EXIT_SUCCESS;
      } else if(strcmp(argv[1], "--lbl") == 0) {  // Jde o delku lomene cary?
        funkce = LBL;
        initDelkaN(&stN);  // Inicializace pro vypocet
      } else {
        fprintf(stderr, "Chybne zadane parametry.\n");
        return EXIT_FAILURE;
      }
    } else if(argc == 3) {  // Pokud ma program dva parametry
      if(strcmp(argv[1], "--arcsin") == 0) {  // Jedna se o sinus?
        funkce = ARCSIN;
        sigdig = atoi(argv[2]);  // Zjistim presnost
        if (sigdig <= 0) {  // Zkontroluji presnost
          fprintf(stderr, "Chybne zadana hodnota sigdig.\n");
          return EXIT_FAILURE;
        } else if(sigdig == INT_MAX || sigdig == INT_MIN || sigdig > DBL_DIG) {
          fprintf(stderr, "Vstupni hodnota sigdig je prilis velka. Max = %d.\n",
                  DBL_DIG);
          return EXIT_FAILURE;
        }
      } else if(strcmp(argv[1], "--lble") == 0) {  // Jde o delku lomene cary?
        funkce = LBLE;
        chyba = atof(argv[2]);  // Zjistim chybu
        if (chyba <= 0.0) {  // Zkontroluji zadanou chybu
          fprintf(stderr, "Chybne zadana hodnota chyby.\n");
          return EXIT_FAILURE;
        } else if(fabs(chyba) == HUGE_VAL) {
          fprintf(stderr, "Vstupni hodnota chyba je prilis velka.\n");
          return EXIT_FAILURE;
        }
        initDelkaE(&stE);  // Inicializuji pro vypocet strukturu
      } else {
        fprintf(stderr, "Chybne zadane parametry.\n");
        return EXIT_FAILURE;
      }
    } else if(argc == 4) {  // Pokud ma program parametry tri
      if(strcmp(argv[1], "--logax") == 0) {  // Jedna se o logaritmus?
        funkce = LOGAX;
        sigdig = atoi(argv[2]);  // Zjistim presnost
        if (sigdig <= 0) {  // Zkontroluji presnostS
          fprintf(stderr, "Chybne zadana hodnota sigdig.\n");
          return EXIT_FAILURE;
        } else if(sigdig == INT_MAX || sigdig == INT_MIN || sigdig > DBL_DIG) {
          fprintf(stderr, "Vstupni hodnota sigdig je prilis velka. Max = %d.\n",
                  DBL_DIG);
          return EXIT_FAILURE;
        }
        zaklad = atof(argv[3]);  // Zjistim zaklad
        if (zaklad <= 0.0) {  // Zkontroluji zaklad
          fprintf(stderr, "Chybne zadana hodnota zakladu.\n");
          return EXIT_FAILURE;
        } else if(fabs(zaklad) == HUGE_VAL) {
          fprintf(stderr, "Vstupni hodnota zaklad je prilis velka.\n");
          return EXIT_FAILURE;
        }
      } else {
        fprintf(stderr, "Chybne zadane parametry.\n");
        return EXIT_FAILURE;
      }
    } else {
      fprintf(stderr, "Prilis mnoho parametru.\n");
      return EXIT_FAILURE;
    }
  }
  
  // Cyklus postupneho nacitani a vyhodnocovani vstupnich hodnot
  while(1) {
    if(funkce > LBLE) {  // Pokud nacitam pouze jednu hodnotu
      result = scanf("%lf", &cislo);
      if(result == EOF) {  // Zkontroluji vstup
        return EXIT_SUCCESS;
      } else if(result == 0) {
        fprintf(stderr, "Spatne zadane cislo.\n");
        return EXIT_FAILURE;
      }
    } else {  // Pokud nacitam hodnoty dve
      result = scanf("%lf", &bod.x);
      if(result == EOF) {  // Zkontroluji prvni
        return EXIT_SUCCESS;
      } else if(result == 0) {
        fprintf(stderr, "Spatne zadana hodnota.\n");
        return EXIT_FAILURE;
      }
      
      result = scanf("%lf", &bod.y);
      if(result == EOF) {  // Zkontroluji druhou
        fprintf(stdout, "%lf\n", NAN);
        return EXIT_SUCCESS;
      } else if(result == 0) {
        fprintf(stderr, "Spatne zadana hodnota.\n");
        return EXIT_FAILURE;
      }
    }
  
    // Podle parametru spustim funkci
    switch(funkce) {
      case LBL: {
    	  vysledek = lomenaCara(&bod, &stN);
    	  break;
      }
      case LBLE: {
  	    lomenaCaraChyba(chyba, bod, &stE);
  		  break;
      }
      case LOGAX: {
        vysledek = logax(cislo, zaklad, sigdig);
        break;
      }
      case ARCSIN: {
        vysledek = arcsin(cislo, sigdig);
        break;
      }
    }
    
    // Rozhodnu, ktere hodnoty budu vypisovat
  	if(funkce == LBLE) {
  		fprintf(stdout, "%.10e\n%.10e\n", stE.interval.x, stE.interval.y);
  	} else {
  		fprintf(stdout, "%.10e\n", vysledek);
  	}
  }  // Konec cyklu while

  return EXIT_SUCCESS;
}
