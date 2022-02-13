# SpooVols
Backs out implied vols for ESmini futures options 

BBG apparenltly stores ES futures options prices...but not the vols. You can see them in ES1 Index OVDV, with the respective strikes but cant just call them via BDH( ).
So building my own to scan sporadically.

Dependencies: 
  Am using FXOption class as the pricer (it's basically same as an equity option if priced in "DomesticPips") and a simple bisection algo to get the vol from the price.
  Some bug on higher strikes ..lower strikes is wihtin 2-3 vols , fix both, also fix the timedelta .  
  
  DownloadData class to pull BBG data for the options prices ..at moment going for ESM2 (2022.06.18) expiries 

