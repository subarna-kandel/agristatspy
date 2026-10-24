# Statistical Methods

This document describes every statistical method implemented in AgriStatsPy, including the model equations, assumptions, limitations, and references.

---

## 1. Experimental designs

### 1.1 Completely Randomized Design (CRD)

**Definition.** Experimental units are assigned to treatments using simple random allocation with no blocking factor.

**Model.**  
$$y_{ij} = \mu + \tau_i + \varepsilon_{ij}, \quad \varepsilon_{ij} \sim N(0, \sigma^2)$$

where $\mu$ is the grand mean, $\tau_i$ is the fixed effect of treatment $i$, and $\varepsilon_{ij}$ is the random error for the $j$-th unit in treatment $i$.

**Invariants.**

- Total units = $k \times r$ (treatments × replications).
- Each treatment appears exactly $r$ times.
- Plot IDs are unique.
- Treatment assignments are randomized using a reproducible seed.

**Limitations.** CRD is most efficient when experimental material is uniform. If the environment is heterogeneous and blocking is not used, the error variance will be inflated.

---

### 1.2 Randomized Complete Block Design (RCBD)

**Definition.** Experimental units are grouped into blocks of homogeneous material. Each block contains every treatment exactly once. Treatments are randomized independently within each block.

**Model (additive).**  
$$y_{ij} = \mu + \tau_i + \beta_j + \varepsilon_{ij}, \quad \varepsilon_{ij} \sim N(0, \sigma^2)$$

where $\beta_j$ is the fixed effect of block $j$.

**Invariants.**

- Total plots = $k \times b$ (treatments × blocks).
- Every treatment appears exactly once per block.
- Each block contains exactly $k$ plots.
- Plot IDs are unique.

**Limitations.** The standard additive model does not separately estimate a treatment-by-block interaction. If interaction is present, the error term is inflated and F-tests are conservative. Check treatment-by-block interaction plots before accepting the additive model.

---

## 2. Descriptive statistics

For each treatment group:

| Statistic | Symbol | Formula |
|---|---|---|
| Count | $n_i$ | Number of non-missing observations |
| Mean | $\bar{y}_i$ | $\frac{1}{n_i}\sum_j y_{ij}$ |
| Standard deviation | $s_i$ | $\sqrt{\frac{\sum_j (y_{ij} - \bar{y}_i)^2}{n_i - 1}}$ |
| Standard error | $\text{SE}_i$ | $s_i / \sqrt{n_i}$ |
| Minimum | — | $\min_j y_{ij}$ |
| Median | — | Empirical median |
| Maximum | — | $\max_j y_{ij}$ |

**Coefficient of Variation (CV%).** Agricultural ANOVA reports use the residual-based CV:

$$\text{CV\%} = \frac{\sqrt{MS_{\text{error}}}}{\bar{y}_{\cdot\cdot}} \times 100$$

where $\bar{y}_{\cdot\cdot}$ is the grand mean and $MS_{\text{error}}$ is the error mean square from the ANOVA. When CV is reported outside the ANOVA context (e.g., from `treatment_summary`), the pooled within-group standard deviation is used instead. **CV is undefined when the grand mean is zero**; in that case it is returned as `NaN` with an explanatory note.

---

## 3. ANOVA

### 3.1 One-way ANOVA (CRD)

**Hypotheses.**
$$H_0: \tau_1 = \tau_2 = \cdots = \tau_k = 0$$
$$H_1: \text{at least one } \tau_i \neq 0$$

**ANOVA table.**

| Source | df | SS | MS | F | p |
|---|---|---|---|---|---|
| Treatment | $k-1$ | $SS_T$ | $MS_T = SS_T/(k-1)$ | $MS_T/MS_E$ | from $F_{k-1,\,n-k}$ |
| Error | $n-k$ | $SS_E$ | $MS_E = SS_E/(n-k)$ | | |
| Total | $n-1$ | $SS_{\text{Total}}$ | | | |

**Assumptions.**

1. Observations are independent.
2. Residuals are approximately normally distributed.
3. Variances are equal across treatment groups (homoscedasticity).

A significant F-test indicates that not all treatment means are equal. It does not identify which treatments differ — use Tukey HSD for pairwise comparisons.

**Statistical significance does not establish agronomic importance.**

---

### 3.2 Two-way additive ANOVA (RCBD)

Implemented via `statsmodels` OLS with the formula:

```
response ~ C(treatment) + C(block)
```

**ANOVA table.**

| Source | df | SS | MS | F | p |
|---|---|---|---|---|---|
| Treatment | $k-1$ | $SS_T$ | $MS_T$ | $MS_T/MS_E$ | |
| Block | $b-1$ | $SS_B$ | $MS_B$ | $MS_B/MS_E$ | |
| Error | $(k-1)(b-1)$ | $SS_E$ | $MS_E$ | | |
| Total | $n-1$ | $SS_{\text{Total}}$ | | | |

**Assumptions.** Same as CRD, plus: additive treatment and block effects (no interaction).

---

## 4. Mean comparisons

### 4.1 Tukey Honest Significant Difference (HSD)

**When to use.** After a significant omnibus F-test, when all pairwise comparisons are of interest.

**Family-wise error rate.** Tukey HSD controls the probability that any comparison is falsely declared significant at the specified $\alpha$ level.

**Implementation.** `statsmodels.stats.multicomp.pairwise_tukeyhsd`, which uses the Tukey–Kramer extension for unequal group sizes.

**Output per comparison.**

- Group labels.
- Mean difference ($\bar{y}_i - \bar{y}_j$).
- Simultaneous 95% confidence interval.
- Adjusted p-value.
- Reject / do not reject decision.

**Assumptions.**

1. Independent observations.
2. Approximately normal residuals.
3. Homogeneous variances.

**Limitations.** Letter groupings (compact letter displays) are a convenience summary and can have multiple valid arrangements. Pairwise comparisons in the table are the authoritative output.

**Only LSD and Duncan's tests will be added in later phases**, after their method specifications and limitations have been independently reviewed and validated.

---

## 5. Statistical safety notices

The following notices apply to all analysis functions in AgriStatsPy:

1. **Non-significance ≠ no effect.** A non-significant test does not prove the null hypothesis.
2. **Assumptions are not proven by diagnostic tests.** A non-significant Levene or Shapiro–Wilk test does not prove variance homogeneity or normality.
3. **Significance ≠ importance.** Statistical significance depends on sample size and does not measure agronomic relevance.
4. **Subsamples ≠ replications.** Multiple measurements from the same plot are subsamples; treating them as independent replications is pseudoreplication.
5. **Missing data.** All functions report missing observations and do not silently remove them.
6. **Ordinary ANOVA is not always appropriate.** Count, proportion, binary, and repeated-measures outcomes may require generalized or mixed models.

---

## 6. References

- Cochran, W. G., & Cox, G. M. (1957). *Experimental Designs* (2nd ed.). Wiley.
- Montgomery, D. C. (2017). *Design and Analysis of Experiments* (9th ed.). Wiley.
- Steel, R. G. D., Torrie, J. H., & Dickey, D. A. (1997). *Principles and Procedures of Statistics: A Biometrical Approach* (3rd ed.). McGraw-Hill.
- Tukey, J. W. (1949). Comparing individual means in the analysis of variance. *Biometrics*, 5(2), 99–114.
