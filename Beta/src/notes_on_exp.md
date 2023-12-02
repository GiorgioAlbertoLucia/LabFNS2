# Notes on exp

Note sugli esperimenti di Lab2 con proposte di pezzi delle esperienze che possono essere divise.

Topic: caratterizzazione rivelatori al silicio con e senza guadagno interno (LGAD vs PiN)
Setups:
    - Probe station
    - Transient Current Technique (TCT)

## Esperienza 1: Probe station

### Obiettivi

- caratterizzazione current-voltage
- caratterizzazione capacitance-voltage

### Steps

1. curva C-frequency per trovare punto di lavoro (da fare in lab direttamente)
2. curva C-V
    - profilo di drogaggio
    - tensione di svuotamento
3. curva 1/C^2-V
    - tensione di svuotamento substrato e gain layer con fit lineari
4. Curva I-V (?)
    - può essere interessante, si può osservare andamento e fittare per verificare il valore di resistenza e capacità del dispositivo (?) -> eventialmente confrontare con simulazione di circuito: uso parametri di C, R trovati e verifico che l'andamento simulato sia compatibile a quello osservato (?)

## Esperienza 2: TCT

### Obiettivi

- calibrazione sistema laser
- misura del guadagno di un LGAD
- misura del jitter di un LGAD

### Steps

1. curva charge-signal amplitude (calibrazione laser) -> carica generata dal laser regolata con intensità modificabile
2. curva gain-Vbb
    - osservare andamento moltiplicazione, geiger mode etc
    - dati da oscilloscopio 
3. curva resolution-gain (andamento dominato dal gain)

## Esperienza 3: beta source

### Obiettivi

- misura risoluzione temporale di un ufsd con sorgente 90Sr

### Steps

1. a partire dal .root creare un dataset con ampiezza, baseline, rms noise, ToA (una per ciascun detector)
2. ottenere la risoluzione dalla distribuzione ToA1, ToA2 (singolarmente) per i tempi di arrivo delle MIP (i.e. elettroni da 2.2 MeV):
la risoluzione effettiva è la somma in quadratura delle due!
    - senza tagli 
    - con tagli su distribuzioni di ampiezze
    - eventualmente facendo threshold scan (verifichiamo per tagli diversi cosa succeda alla distribuzione)
3. misura rise time del sensore
    - facendo fit
4. riduzione noise con fft per sopprimerlo (?)
5. studio sistematici cercando eventuale shift nella distribuzione delle energie nel tempo (?)