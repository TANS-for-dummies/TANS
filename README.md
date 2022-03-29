# TANS
Ricostruzione del vertice di impatto per un collider.

Per eseguire le 2 macro di simulazione e ricostruzione prima compilare CompilaLibreria.C
Le simulazioni sono ottenibili da ALICE_Junior.C, la cui macro è MonteCarlo()

  MonteCarlo(int N_esp = 1000000, const char* output_file = "MonteCarlo.root", int gen = 1, bool scat = 1, bool smear = 1, int N_false_hit = 0, unsigned int seed = 125, int verboseEvent = 3154)
    N_esp è il numero degli esperimenti
    output_file è il nome del file di output generato (solo.root)
    gen definisce il metodo di generazione delle molteplicità, 1 per la distribuzione fornita, 2 da una distribuzione uniforme, o 3 per molteplicità fissa
    scat se messo a 0 spegne il multiscattering
    smear se messo a 0 spegne lo smearing
    N_false_hit è il numero di false hit
    seed è il seed per la generazione di numeri casuali
    
    
Per la ricostruzione si fa riferimento a Ricostruzione_Vertice.C, la cui macro è Ricostruzione_Vertice()

  Ricostruzione_Vertice(const char* input = "MonteCarlo.root", double window_size = 0.5, double window_step = 0.25, int n_sigma = 3)
    input è il nome del file in input (solo .root)
    window_size è la dimensione della running window che scorre sul vector dei tracklet (è in CM)
    window_step è il passo con cui avanza la window
    n_sigma serve poichè escludiamo le z dei tracklet lasciando solo quelle per cui z < n_sigma * 5.3
    

Il metodo di ricostruzione è la running window che scorre su tutto il vettore ordinato delle z dei tracklet.
Molte parti commentate sono per il debugging del programma, sono state lasciate nel caso servissero.
