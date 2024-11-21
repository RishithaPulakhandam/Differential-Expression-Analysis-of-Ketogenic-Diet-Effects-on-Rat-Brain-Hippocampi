# README: Differential Expression Analysis of Rat Ketogenic Diet Dataset

## Project Overview
This project conducts a differential expression analysis to identify genes or probes that are significantly different between rats on a ketogenic diet (KD) and those on a control diet. The dataset is sourced from the GEO database and processed using various statistical and visualization techniques to derive insights.

### Objectives:
1. Identify genes/probes with significant expression changes between the two groups.
2. Understand the biological implications of these genes, focusing on the differences induced by the ketogenic diet in the hippocampal region of rat brains.
3. Visualize the differential expression using histograms and volcano plots for exploratory and explanatory purposes.

---

## Data Source:
- Dataset: Rat Ketogenic Diet Data (`rat_KD.txt`)
- Format: Affymetrix Array Expression Data
- Groups: 
  - KD (Ketogenic Diet): 6 samples
  - Control Diet: 6 samples

---

## Steps Performed:

### 1. Data Preprocessing
- Log2 transformation of the dataset to normalize the expression values.

```r
rat_KD_log2 <- log2(rat_KD)
```

### 2. Differential Expression Analysis
- A Student’s t-test was applied to each gene/probe to compare the two groups.

```r
t.test.all.genes <- function(x,s1,s2) {
  x1 <- as.numeric(x[s1])
  x2 <- as.numeric(x[s2])
  t.out <- t.test(x1, x2, alternative="two.sided", var.equal=TRUE)
  return(as.numeric(t.out$p.value))
}
pv <- apply(rat_KD_log2, 1, t.test.all.genes, s1=1:6, s2=7:12)
```

### 3. Visualization of p-values
- A histogram of p-values and their `-log10` transformations was plotted to observe the distribution of significance levels.

```r
hist(pv, col="lightblue", xlab="p-values", main="P-value Distribution")
hist(-log10(pv), col="lightblue", xlab="log10(p-values)", main="-log10(p-values)")
```

### 4. Thresholding and Bonferroni Correction
- Identified the number of probes below p-value thresholds (e.g., p < 0.05, p < Bonferroni-corrected threshold).

---

## Results:
1. **Significantly Differential Probesets**:
   - 5,160 probes had a p-value < 0.05.
   - 2,414 probes had a p-value < 0.01.
   - Using Bonferroni correction, 12 probes were identified as significant.

2. **Biological Insights**:
   - Top significant probes belonged to the hemoglobin gene family, including alpha and beta chains.
   - Biological Function: Hemoglobins act as respiratory carriers, transporting oxygen and carbon dioxide.

---

## Visualizations:

### 1. Histograms:
- **p-value Distribution**:
  - A majority of probes had low p-values, indicating strong differences between the two groups.<br>
  <img width="320" alt="image" src="https://github.com/user-attachments/assets/f3690dd6-c140-4687-a242-f77620d3f9d1">
 
- **-log10(p-values)**:
  - Highlighted the significance of small p-values.<br>
  <img width="307" alt="image" src="https://github.com/user-attachments/assets/2fbe2c1f-3348-4c41-90b5-2de5e7b28888">


### 2. Volcano Plot:
- **Highlights**:
  - X-axis: Transformed p-values (`-log10`).
  - Y-axis: Fold changes in log2 scale.
  - Vertical line: Threshold for significance (`p = 0.05`).
  - Horizontal lines: Fold change thresholds (log2 |fold| > 1).

```r
plot(range(p_value), range(fold), type = 'n', xlab = '-1*log10(p values)', ylab = 'Fold change', main = 'Volcano Plot')
points(p_value, fold, col='black')
abline(v =-log10(0.05), col ='blue', lty=2)
abline(h =1, col="red", lty=2)
abline(h=-1, col="red", lty=2)
```
  <img width="541" alt="image" src="https://github.com/user-attachments/assets/b8568cb4-0d56-4778-955a-ec1c5aa3b165">


---

## Key Findings:
- The ketogenic diet significantly alters the expression of genes related to respiratory functions and possibly energy metabolism.
- Hemoglobin-related genes, such as **Hba** and **Hbb**, showed the most significant changes in expression levels.
- These genes' differential expression could contribute to adaptive responses to altered energy states in the brain.

---

## Files Included:
1. **rat_KD.txt**: Original expression data.
2. **Scripts**: R code for analysis and visualization.
3. **Visualizations**:
   - Histograms of p-values and `-log10(p-values)`
   - Volcano Plot

---

## Acknowledgments:
This project was conducted as part of Lab 5 coursework. The dataset is sourced from publicly available GEO data, and analysis techniques follow standard bioinformatics protocols.
