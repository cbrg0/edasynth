# References

This page lists every source the project uses and says what each source is used for.
The last column says whether the statement that this repository attributes to the
source was compared with the text of the source. "Checked" means that it was.
"Background" means that the source is cited for a general point that was not compared
with its text again. The references follow the APA style, 7th edition.

## Source usage

| Source | What it is used for | Compared with the text |
|---|---|---|
| Birjandtalab et al. (2016, 2017) | The UTD recordings: 20 participants, the sensors on the wrist and the stages of the protocol. | checked on the page of the data set |
| Zhang et al. (2017) | The expert labels for UTD and the criteria the experts used, the five-second windows, the feature sets F0 and FA, the share of artifacts of 3.8 % in UTD and 15 % in the resting recordings of their second data set, the small gain from acceleration with its two possible explanations, and the device. | checked |
| IdeasLabUT (n.d.) | The files with the labels and the feature tables of Zhang et al., used at commit `7c4bd46`. | used directly |
| Greco et al. (2016) | The shape of a response, the constants τ₁ = 0.7 s, τ₀ between 2 and 4 s, α = 8·10⁻⁴ and γ = 10⁻², the curve pieces joined every 10 seconds, the choice of τ₀ by the smallest difference between model and signal, from which this project departs, and the synthetic tonic level made of a linear trend and a sine wave with a period of 45 to 90 seconds. | checked |
| lciti (n.d.) | The authors' own code of cvxEDA, which can be used instead of the built-in solver for comparison. | not used so far |
| Kleckner et al. (2018) | The valid range of EDA from 0.05 to 60 µS (their rule 1), the margin around invalid data (their rule 4) and 0.05 µS as the smallest size of a response. | checked |
| Taylor et al. (2015) | Five-second pieces of the signal, the criteria experts use to mark artifacts, and the fact that their classifier did not receive the acceleration. | checked |
| Posada-Quintero & Chon (2020) | The start of a response one to five seconds after a stimulus and the threshold of 0.05 µS, both of which they cite from Boucsein et al. (2012); the disagreement between decomposition methods; general background on EDA. | checked |
| Hossain et al. (2022) | The measure *P<sub>n</sub>*, the share of power above 0.4 Hz, and the finding that walking artifacts in EDA measured on the foot are mostly fast. | checked |
| Chen et al. (2015) | Removal of artifacts with wavelets. | background |
| Donoho & Johnstone (1994) | The estimate of noise from the median of the finest wavelet details. | background, a standard estimator |
| Amin et al. (2024) | The accelerometer as a reference signal for reducing artifacts. | background, checked only against the title |
| U.S. Food and Drug Administration (2018) | The Embrace wristband was cleared as an aid to seizure monitoring "during periods of rest", using EDA and motion data. | checked |
| Llanes-Jurado et al. (2023) | The EDABE data set: 43 recordings from a Shimmer3 GSR+ device with expert artifact labels. | checked on the page of the data set |
| Schmidt et al. (2018) | The WESAD data set: EDA and acceleration from an Empatica E4 wristband, the stages of the protocol, and acceleration stored in units of 1/64 g. | checked in its description file |
| Pedregosa et al. (2011) | The scikit-learn library, which trains the detectors. | not applicable |
| Pollard et al. (2026) | The citation that PhysioNet asks users of its data to include. | checked on the page of the data set |

Boucsein et al. (2012) is cited only as Posada-Quintero & Chon (2020) summarise it.

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

Llanes-Jurado, J., Carrasco-Ribelles, L. A., Alcañiz, M., Soria-Olivas, E., & Marín-Morales, J. (2023). Automatic artifact recognition and correction for electrodermal activity based on LSTM-CNN models. *Expert Systems with Applications, 230*, Article 120581. https://doi.org/10.1016/j.eswa.2023.120581

Pedregosa, F., Varoquaux, G., Gramfort, A., Michel, V., Thirion, B., Grisel, O., Blondel, M., Prettenhofer, P., Weiss, R., Dubourg, V., Vanderplas, J., Passos, A., Cournapeau, D., Brucher, M., Perrot, M., & Duchesnay, É. (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research, 12*, 2825–2830. https://jmlr.org/papers/v12/pedregosa11a.html

Pollard, T., Moody, B. E., Lehman, L., Gow, B., Fernandes, C., Xie, C., Johnson, A., Mark, R. G., & Heldt, T. (2026). PhysioNet as a global platform for biomedical research. *Nature Health*. https://doi.org/10.1038/s44360-026-00096-z

Posada-Quintero, H. F., & Chon, K. H. (2020). Innovations in electrodermal activity data collection and signal processing: A systematic review. *Sensors, 20*(2), Article 479. https://doi.org/10.3390/s20020479

Schmidt, P., Reiss, A., Duerichen, R., Marberger, C., & Van Laerhoven, K. (2018). Introducing WESAD, a multimodal dataset for wearable stress and affect detection. In *Proceedings of the 20th ACM International Conference on Multimodal Interaction* (pp. 400–408). https://doi.org/10.1145/3242969.3242985

Taylor, S., Jaques, N., Chen, W., Fedor, S., Sano, A., & Picard, R. (2015). Automatic identification of artifacts in electrodermal activity data. In *2015 37th Annual International Conference of the IEEE Engineering in Medicine and Biology Society (EMBC)* (pp. 1934–1937). IEEE. https://doi.org/10.1109/EMBC.2015.7318762

U.S. Food and Drug Administration. (2018). *510(k) premarket notification K181861: Embrace*. https://www.accessdata.fda.gov/cdrh_docs/pdf18/K181861.pdf

Zhang, Y., Haghdan, M., & Xu, K. S. (2017). Unsupervised motion artifact detection in wrist-measured electrodermal activity data. In *Proceedings of the 2017 ACM International Symposium on Wearable Computers* (pp. 54–57). Association for Computing Machinery. https://doi.org/10.1145/3123021.3123054
