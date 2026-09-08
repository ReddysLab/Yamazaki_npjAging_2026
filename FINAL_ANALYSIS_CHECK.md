# Final analysis check

The clean package was tested on 12 August 2026 using the bundled public inputs and the locally downloaded Schaum et al. GSE132040 object.

## Participant-aware sleep-restriction model

- 19,541 genes, 427 samples, 26 participants.
- `duplicateCorrelation` consensus correlation: 0.520 for the phase-factor model and 0.513 for the cosinor model.
- Final dual-FDR SR-down set: 2,163 genes.
- Heat-stress response module: median delta MESOR = -0.0339; empirical P = 0.0008.
- HSF1 heat-shock regulation module: median delta MESOR = -0.0336; empirical P = 0.0026.

## Downstream checks

- Fig. 3e excitatory neurons: universe 13,930; SR-down 1,691; aging-down 1,095; overlap 201; OR 1.7118; P = 3.49e-10.
- Fig. 3e inhibitory neurons: universe 13,930; SR-down 1,691; aging-down 877; overlap 127; OR 1.2439; P = 0.0178.
- Fig. 4c brain: universe 5,086; SR-down 741; RNA/protein aging-decline 1,634; overlap 202; OR 0.7624; one-sided P = 0.9992. This overlap is not enriched.
- Supplementary Fig. 4e: universe 15,442; SR-down 1,947; JenAge aging-down 7,916; overlap 1,363; OR 2.4723; P = 4.36e-72.

All final Reactome exports include the query set, background universe, raw P value, BH FDR, overlap count, and overlap genes. Displayed network terms use nominal Reactome P < 0.05 and at least two overlapping genes; no term is forced.
