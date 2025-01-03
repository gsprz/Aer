Descrizioni e riferimenti al report:

- weissinger_cessna_da40.zip: in questa cartella zippata è presente il main di weissinger che va compilato assieme a tutte le funzioni e il file .dat presenti all'interno della stessa cartella. A seconda del caso che si vuole analizzare cambiare la variabile "riferimento", per il caso "riferimento = 2" si può cambiare la curvatura del profilo aggiungendo la linea media mettendo "config.lm = [1]", lasciare tale valore a 0 se non si vuole curvatura. Per il caso scelto appariranno il plot della pannellizzazione con la distribuzione di circolazione, la circolazione lungo la span e gli angoli indotti al quarto di corda (SEZIONE 2.1 DEL REPORT).

- analysis_airfoil.zip: in questa zippata sono presenti dei file .dat e un file .m, quest'ultimo va compilato quando i file .dat sono all'interno della sua stessa cartella, il codice all'interno del file calcola alpha di theodorsen e alpha 0 per i profili richiesti e per un profilo sottile aggiunto per confutare i risultati, una volta compilato ci saranno delle celle di alpha_th e alpha_0 (la prima cella = B737D, la seconda = NACA0012 e la terza = NACA4402), in tali celle sono presenti i valori degli alpha variando i numeri di punti interpolati e numero dei punti di integrazione, tale metodo è stato utilizzato per verificare che convergessero allo stesso risultato ad alti numeri di punti di interpolazione; oltre a questo sono presenti gli errori rispetto all'alpha0 e i plot dei tre profili con le corrispondenti linee medie (il file potrebbe metterci molto tempo a compilare e per profili sottili bisogna assicurarsi che il file .dat contenga un numero dispari di pannelli) (SEZIONE 1.2 DEL REPORT).

- weissinger_main_finale.m: lo script deve essere runnato insieme a delle funzioni presenti nella cartella weissinger_cessna_da40.zip. Esso plotta i tre grafici relativi alle polari dei due velivoli con la possibilità di mettere e togliere i piani di coda. Alla fine dello script sono presenti i valori dei parametri geometrici da inserire nelle due sezioni corrispondenti alle code e alle ali dei due velivoli. Le due sezioni in cui modificare i dati sono alle righe 19 e le righe 350 (SEZIONE 2.2 E 2.3 DEL REPORT).

- Grafico_gamma_tandem.m: lo script deve essere runnato insieme a delle funzione presenti nella cartella weissinger_cessna_da40.zip. Esso plotta la distribuzione di circolazione lungo lo span considerando sia l'ala che la coda del velivolo (SEZIONE 2.1 DEL REPORT).

- alaellittica.m: lo script permette di calcolare e fare il grafico della distribuzione ellittica di circolazione associata al velivolo i cui parametri sono stati inseriti nel main. eseguire lo script alaellittica.m dopo aver eseguito il main di weissinger_cessna_da40.zip (SEZIONE 2.1 DEL REPORT).

- hess_smith.zip: il main di questa cartella va compilato con i file .dat e le mat functions tutte in una stessa cartella, esso presenta cp, cl e cm dei profili, ponendo la variabile 'riferimento' pari a 1 si analizza il NACA0012, mentre ponendo la variabile 'riferimento' pari a qualsiasi altro numero si analizza il profilo B737D (SEZIONE 1.1 DEL REPORT).

- alfa_progetto.zip: altra metodologia basata sull'uniformità del bordo di attacco tramite XFOIL per il confronto dell'alpha di Theodorsen calcolato numericamente in analysis_airfoil.zi (SEZIONE 1.2 DEL REPORT).

- alfa_zerolift.zip: altra metodologia basata su XFOIL per il confronto dell'alpha 0 calcolato numericamente in analysis_airfoil.zip (SEZIONE 1.2 DEL REPORT).

- transizione_&_separazione.zip: script che tramite l'utilizzo di XFOIL analizza transizione e separazione dei profili del NACA0012 e del B737D (SEZIONE 1.3 DEL REPORT).
