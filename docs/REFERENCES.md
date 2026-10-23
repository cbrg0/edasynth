# References (APA 7) and what each source contributes

| Source | Used for |
|---|---|
| Birjandtalab et al. (2016, 2017) | UTD data: protocol, device, signals |
| Zhang et al. (2017) | expert labels, 5-s windows, F0 and FA features, LOSO, the two explanations that define the map's axes, prevalence 3.8 % and 15 % |
| Greco et al. (2016); lciti (n.d.) | Bateman response, τ1 = 0.7 s, τ0 ∈ [2, 4] s, α, γ, Poisson driver, tonic form; decomposition used to estimate parameters |
| Hossain et al. (2022) | high-frequency walking artifacts vs in-band electrode-shift artifacts; P_n detector |
| Kleckner et al. (2018) | rules 1, 2, 4 as a detector; 5-s transition zone; 0.05 uS minimum response; range limit for class C5 |
| Chen et al. (2015) | SWT detector; pressure-induced steep rise (class C1); additive artifact model; accelerometer stillness criterion (comparison) |
| Taylor et al. (2015) | 5-s epochs; artifact criteria (classes C1, C2, C5; quantisation); Kocielnik relative-rate rule (secondary citation); the argument against using the accelerometer that goal 4 tests |
| Posada-Quintero & Chon (2020) | response latency 1-5 s and 0.05 uS (Boucsein, 2012, secondary citation); 0.045-0.15 Hz band (G3); movement can evoke real responses; disagreement between decomposition methods |
| Donoho & Johnstone (1994) | noise estimate from level-1 wavelet details |
| Amin et al. (2024) | accelerometer as a noise reference for artifact reduction (context for F_B) |
| Pedregosa et al. (2011) | scikit-learn |
| Pollard et al. (2026) | standard PhysioNet citation requested by the data set page |

Amin, M. R., Alam, S., Khazaei, S., Azgomi, H. F., & Faghih, R. T. (2024). Skin conductance response artifact reduction: Leveraging accelerometer noise reference and deep breath detection. *IEEE Access, 12*, 68208–68231. https://doi.org/10.1109/ACCESS.2024.3393515

Birjandtalab, J., Cogan, D., Baran Pouyan, M., & Nourani, M. (2016). A non-EEG biosignals dataset for assessment and visualization of neurological status. In *2016 IEEE International Workshop on Signal Processing Systems (SiPS)* (pp. 110–114). IEEE. https://doi.org/10.1109/SiPS.2016.27

Birjandtalab, J., Cogan, D., Baran Pouyan, M., & Nourani, M. (2017). *Non-EEG dataset for assessment of neurological status* (Version 1.0.0) [Data set]. PhysioNet. https://doi.org/10.13026/C26Q2Z

Chen, W., Jaques, N., Taylor, S., Sano, A., Fedor, S., & Picard, R. W. (2015). Wavelet-based motion artifact removal for electrodermal activity. In *2015 37th Annual International Conference of the IEEE Engineering in Medicine and Biology Society (EMBC)* (pp. 6223–6226). IEEE. https://doi.org/10.1109/EMBC.2015.7319814

Donoho, D. L., & Johnstone, I. M. (1994). Ideal spatial adaptation by wavelet shrinkage. *Biometrika, 81*(3), 425–455. https://doi.org/10.1093/biomet/81.3.425

Greco, A., Valenza, G., Lanata, A., Scilingo, E. P., & Citi, L. (2016). cvxEDA: A convex optimization approach to electrodermal activity processing. *IEEE Transactions on Biomedical Engineering, 63*(4), 797–804. https://doi.org/10.1109/TBME.2015.2474131

Hossain, M.-B., Kong, Y., Posada-Quintero, H. F., & Chon, K. H. (2022). Comparison of electrodermal activity from multiple body locations based on standard EDA indices' quality and robustness against motion artifact. *Sensors, 22*(9), Article 3177. https://doi.org/10.3390/s22093177

IdeasLabUT. (n.d.). *EDA-Artifact-Detection* [Computer software]. GitHub. https://github.com/IdeasLabUT/EDA-Artifact-Detection

Kleckner, I. R., Jones, R. M., Wilder-Smith, O., Wormwood, J. B., Akcakaya, M., Quigley, K. S., Lord, C., & Goodwin, M. S. (2018). Simple, transparent, and flexible automated quality assessment procedures for ambulatory electrodermal activity data. *IEEE Transactions on Biomedical Engineering, 65*(7), 1460–1467. https://doi.org/10.1109/TBME.2017.2758643

lciti. (n.d.). *cvxEDA* [Computer software]. GitHub. https://github.com/lciti/cvxEDA

Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., Blondel, M., Prettenhofer, P., Weiss, R., Dubourg, V., Vanderplas, J., Passos, A., Cournapeau, D., Brucher, M., Perrot, M., & Duchesnay, É. (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research, 12*, 2825–2830. https://jmlr.org/papers/v12/pedregosa11a.html

Pollard, T., Moody, B. E., Lehman, L., Gow, B., Fernandes, C., Xie, C., Johnson, A., Mark, R. G., & Heldt, T. (2026). PhysioNet as a global platform for biomedical research. *Nature Health, 1*(8), 792–795. https://doi.org/10.1038/s44360-026-00096-z

Posada-Quintero, H. F., & Chon, K. H. (2020). Innovations in electrodermal activity data collection and signal processing: A systematic review. *Sensors, 20*(2), Article 479. https://doi.org/10.3390/s20020479

Taylor, S., Jaques, N., Chen, W., Fedor, S., Sano, A., & Picard, R. (2015). Automatic identification of artifacts in electrodermal activity data. In *2015 37th Annual International Conference of the IEEE Engineering in Medicine and Biology Society (EMBC)* (pp. 1934–1937). IEEE. https://doi.org/10.1109/EMBC.2015.7318762

Zhang, Y., Haghdan, M., & Xu, K. S. (2017). Unsupervised motion artifact detection in wrist-measured electrodermal activity data. In *Proceedings of the 2017 ACM International Symposium on Wearable Computers* (pp. 54–57). Association for Computing Machinery. https://doi.org/10.1145/3123021.3123054

Secondary citations (read through the citing work only): Boucsein (2012) via Posada-Quintero & Chon (2020); Hammerla & Plötz (2015) via Zhang et al. (2017); Jain et al. via Posada-Quintero & Chon (2020); Kocielnik et al. via Taylor et al. (2015).
