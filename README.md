# DP_Nezavislost_a_ordinalne_data
Toto obsahuje iba poznamky k jednotlivym scriptom, pre viac informacii treba precitat diplomovu pracu.

Vsetky_Funkcie.R: obsahuje esencialne funkcie ku praci

vzorec_clanok1976.R: obsahuje vzorec pre varianciu tau_b z Agrestiho clanku 1976. Pouzivame na numericke porovnanie hodnot Agrestiho variancie pre tau_b pre jeden multinomicky vyber z tohto clanku, nami odvodenu spravnu varianciu a Agrestiho varianciu z knihy

porovnanie_taub_H0plati.R: porovnanie simulacii pre tau_b - jeden mult. vyber a 3 mult. vybery; pravdepodobnostne tabulky splnaju H0
porovnanie_taub_H1plati.R: porovnanie simulacii pre tau_b - jeden mult. vyber a 3 mult. vybery; pravdepodobnostna tabulka splna H1

chyba_Agresti_H0plati.R: vplyv Agrestiho chyby na simulacie za platnosti H0

Agresti_chyba_H1plati+relativna_chyba.R: obsahuje vysledky ku simuláciám, ak použijeme správne (nami odvodené) phi a nesprávne (Agrestiho) phi.

chyba_Agresti_popul_relat_chyba_popul_taub.R: vypocet populacnej relativnej chyby medzi spravnou (nami odvodenou) varianciou a nespravnou (Agrestiho) varianciou. Vypocet populacneho tau_b.

funkcia_KendallTauB_porovnanie.R: numericke porovnanie hodnot so spravnou (nami odvodenou), nespravnou (Agrestiho) varianciou a varianciou z DescTools balika funkcia KendallTauB.

relat_chyba_hadikove_sigmy_taub.R: vyberova relativna chyba medzi varianciami pre 1 a 3 mult. vybery pre tau_b

zlucovanie_kategorii_taub.R: vplyv zlucovania kategorii na simulacie pre tau_b

gamma_zeta-transform.R: simulacie gamma vs. zeta-transformacia

porovnanie_gamma_H0plati.R: porovnanie simulacii pre gamma - jeden mult. vyber a 3 mult. vybery; pravdepodobnostne tabulky splnaju H0

porovnanie_taub_H1plati.R: porovnanie simulacii pre gamma - jeden mult. vyber a 3 mult. vybery; pravdepodobnostna tabulka splna H1

zlucovanie_kategorii_gamma.R: vplyv zlucovania kategorii na simulacie pre gamma

relat_chyba_popul_gamma.R: relativna chyba medzi varianciami pre jeden vs. 3 mult. vybery; vypocet populacnej gammy

5kap.R: simulacie k 5. kapitole - porovnanie vsetkych testov v simulaciach

rho_rhob_simulacie.R: simulacie k rho a rho_b


## Poznámka k použitiu umelej inteligencie

Niektoré časti kódu (najmä konštrukcia kandidátnych pravdepodobnostných tabuliek) boli vytvorené s podporou nástrojov umelej inteligencie. Tieto nástroje boli využité na preskúmanie možných konfigurácií a urýchlenie procesu hľadania vhodných riešení.

Všetky finálne výsledky, vrátane pravdepodobnostných tabuliek, simulácií a štatistických analýz, boli nezávisle overené, validované a interpretované autorkou práce.
