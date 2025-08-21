## Pipeline for curating intracranial EEG data

Our research aims to uncover mechanistic explanations of the neural basis of human behavior, that is, move from where to how. Our goals are multifaceted: (1) advance fundamental science by discovering new knowledge using rigorous, reproducible methods; and (2) advance translational applications in neurotechnology, precision medicine, and product development that are grounded in rigorous science. Intracranial EEG (iEEG) is an invaluable tool to study the human brain because it provides data with the high spatiotemporal resolution and signal quality traditionally limited to animal neurophysiology. Wrangling iEEG data is complex because the data come from neurosurgical patients, typically with medication-resistant epilepsy, and every dataset is different. With careful cleaning to remove seizure-related data, iEEG data can represent healthy brain activity. 

These scripts read and process raw iEEG data to standardize and prepare the data for analysis. The scripts also read and process MRI and CT scans to reconstruct the iEEG electrode positions in native and MNI space, and merge the MNI coordinates with the ready-to-analyze iEEG data. The resulting data are in a standardized MATLAB structure in FieldTrip format, enabling use in MATLAB (with or without FieldTrip) and straightforward conversion to other formats such as Python.

The pipeline is described in:
- Johnson, EL, Lin, JJ, King-Stephens, D, Weber, PB, Laxer, KD, Saez, I, Girgis, F, D’Esposito, M, Knight, RT, Badre, D. A rapid theta network mechanism for flexible information encoding. _Nature Communications_ 14 (2023). [DOI](https://doi.org/10.1038/s41467-023-38574-7)
- Cross, ZR, Gray, SM, Dede, AJO, Rivera, YM, Yin, Q, Vahidi, P, Rau, EMB, Cyr, C, Holubecki, AM, Asano, E, Lin, JJ, Kim McManus, O, Sattar, S, Saez, I, Girgis, F, King-Stephens, D, Weber, PB, Laxer, KD, Schuele, SU, Rosenow, JM, Wu, JY, Lam, SK, Raskin, JS, Chang, EF, Shaikhouni, A, Brunner, P, Roland, JL, Braga, RM, Knight, RT, Ofen, N, Johnson, EL. The development of aperiodic neural activity in the human brain. _Nature Human Behaviour_ (2025). [DOI](https://doi.org/10.1038/s41562-025-02270-x)
- Rau, EMB, Fellner, M-C, Heinen, R, Zhang, H, Yin, Q, Vahidi, P, Kobelt, M, Asano, E, Kim McManus, O, Sattar, S, Lin, JJ, Auguste, KI, Chang, EF, King-Stephens, D, Weber, PB, Laxer, KD, Knight, RT, Johnson, EL, Ofen, N, Axmacher, N. Reinstatement and transformation of memory traces for recognition. _Science Advances_ 11 (2025). [DOI](https://doi.org/10.1126/sciadv.adp9336)
- Dede, AJO, Cross, ZR, Gray, SM, Kelly, JP, Yin, Q, Vahidi, P, Asano, E, Schuele, SU, Rosenow, JM, Wu, JY, Lam, SK, Raskin, JS, Lin, JJ, Kim McManus, O, Sattar, S, Shaikhouni, A, King-Stephens, D, Weber, PB, Laxer, KD, Brunner, P, Roland, JL, Saez, I, Girgis, F, Knight, RT, Ofen, N, Johnson, EL. Declarative memory through the lens of single-trial peaks in high frequency power. _bioRxiv_ (2025). [DOI](https://doi.org/10.1101/2025.01.02.631123)
- Shi, L, Chattopadhyay, K, Gray, SM, Yarbrough, JB, King-Stephens, D, Saez, I, Girgis, F, Shaikhouni, A, Schuele, SU, Rosenow, JM, Asano, E, Knight, RT, Johnson, EL. Distributed theta networks support the control of working memory: Evidence from scalp and intracranial EEG. _bioRxiv_ (2025). [DOI](https://doi.org/10.1101/2025.08.14.670214)

Software:
- MATLAB 8.6 (R2021b; last tested with R2024b)
- FieldTrip (last tested with fieldtrip-20250114) - [download](https://www.fieldtriptoolbox.org/download/)
- FreeSurfer 7.1 - [download](https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall)

Notes:
- Functions must be run in order from #0-6.
- The two functions #0 can be run in any order. Function #1 reads the outputs of both functions #0.
- Function #1 reconstructs the iEEG electrode positions, analysis described in: Stolk, A, Griffin, S, van der Meij, R et al. Integrated analysis of anatomical and electrophysiological human intracranial data. _Nature Protocols_ 13 (2018). [DOI](https://doi.org/10.1038/s41596-018-0009-6)
- Function #1 requires FieldTrip templates of the standard MNI brain (which are too big for GitHub) - [download](https://drive.google.com/file/d/1qWP-v-ytWxXUuKWydSksg1X1hvZjBi8P/view?usp=sharing)
- Function #4 timestamps the data and is necessarily task-specific. It corresponds to the [working memory delayed match to sample experimental paradigm](https://github.com/elizljohnson-projects/paradigm-working-memory-dms.git). It may be used as a template for other tasks, permitting read in from the outputs of function #3 and read out to function #5.
- All other functions are task-agnostic.
- These functions assume that events were marked using a photodiode sensor. For a step-by-step guide to build one, see photodiode_how_to.
