# References

This page lists every source used by the project, what it is used for, and whether
the specific claim was checked against the source text. References follow APA 7.

## What each source contributes

"Checked" means that the statement attributed to the source in this repository was
compared with the text of the source on 23 October 2026. "Not re-checked" means the
source is cited for a general point that has not been compared with the text again.

| Source | Used for | Checked |
|---|---|---|
| Birjandtalab et al. (2016, 2017) | UTD data: 20 participants, wrist sensors, protocol stages | checked (data set page) |
| Zhang et al. (2017) | expert labels for UTD and AWW; labelling criteria; 5-s windows; F0 and FA features; artifact shares 3.8 % (UTD) and 15 % (AWW rest); the gain from acceleration and its two explanations; device | checked |
| IdeasLabUT (n.d.) | the label files and feature matrices of Zhang et al., commit `7c4bd46` | used directly |
| Greco et al. (2016) | Bateman response; τ<sub>1</sub> = 0.7 s and τ<sub>0</sub> ∈ [2, 4] s; α = 8·10⁻⁴, γ = 10⁻², 10-s knots; τ<sub>0</sub> chosen by the residual norm (from which this project departs); synthetic tonic = linear trend + sinusoid with period 45–90 s | checked |
| lciti (n.d.) | reference implementation of cvxEDA, used as an optional cross-check | not used yet |
| Kleckner et al. (2018) | rules 1, 2 and 4 (0.05–60 µS, ±10 µS/s, 5 s around invalid data); 0.35 Hz low-pass FIR filter of 2 057 taps at 32 Hz; 0.05 µS as accepted minimum response amplitude | checked |
| Taylor et al. (2015) | 5-s epochs; expert criteria (peak without exponential decay, quantisation error, sudden change correlated with motion, level ≤ 0); acceleration not given to their classifier and the seizure example; relative-rate rule of Kocielnik et al. (+20 %/s, −10 %/s) | checked |
| Posada-Quintero & Chon (2020) | response onset 1–5 s after a stimulus and 0.05 µS threshold (both citing Boucsein et al., 2012); disagreement between decomposition methods; general background on EDA | checked |
| Hossain et al. (2022) | P<sub>n</sub>, the share of power above 0.4 Hz; walking artifacts in foot EDA being mostly high-frequency | checked |
| Chen et al. (2015) | wavelet-based artifact removal, basis of the SWT reference detector (Haar, 8 levels, Gaussian-mixture thresholds, δ = 0.01, 400-s windows) | not re-checked |
| Donoho & Johnstone (1994) | noise estimate from the median absolute level-1 wavelet detail | not re-checked (standard estimator) |
| Amin et al. (2024) | accelerometer as a noise reference for artifact reduction (background) | not re-checked beyond the title |
| U.S. Food and Drug Administration (2018) | Embrace cleared as an adjunct to seizure monitoring "during periods of rest", using EDA and motion data | checked |
| Pedregosa et al. (2011) | scikit-learn | — |
| Pollard et al. (2026) | citation requested by PhysioNet for its data | checked (data set page) |

Two sources are cited only through other works and have not been read directly:
Boucsein et al. (2012), through Posada-Quintero & Chon (2020), and Kocielnik et al.,
through Taylor et al. (2015).

The parameters of the SWT detector (Chen et al., 2015) and the role of Amin et al.
(2024) should be compared with the full texts before results that depend on them are
reported.

## Reference list

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

Pollard, T., Moody, B. E., Lehman, L., Gow, B., Fernandes, C., Xie, C., Johnson, A., Mark, R. G., & Heldt, T. (2026). PhysioNet as a global platform for biomedical research. *Nature Health*. https://doi.org/10.1038/s44360-026-00096-z

Posada-Quintero, H. F., & Chon, K. H. (2020). Innovations in electrodermal activity data collection and signal processing: A systematic review. *Sensors, 20*(2), Article 479. https://doi.org/10.3390/s20020479

Taylor, S., Jaques, N., Chen, W., Fedor, S., Sano, A., & Picard, R. (2015). Automatic identification of artifacts in electrodermal activity data. In *2015 37th Annual International Conference of the IEEE Engineering in Medicine and Biology Society (EMBC)* (pp. 1934–1937). IEEE. https://doi.org/10.1109/EMBC.2015.7318762

U.S. Food and Drug Administration. (2018). *510(k) premarket notification K181861: Embrace*. https://www.accessdata.fda.gov/cdrh_docs/pdf18/K181861.pdf

Zhang, Y., Haghdan, M., & Xu, K. S. (2017). Unsupervised motion artifact detection in wrist-measured electrodermal activity data. In *Proceedings of the 2017 ACM International Symposium on Wearable Computers* (pp. 54–57). Association for Computing Machinery. https://doi.org/10.1145/3123021.3123054
