# Mercaldo Calculations Module

## Overview
The Mercaldo Calculations Module provides a robust Python library for statistical analysis of binary classification tests. It includes functionality to compute key metrics such as Sensitivity (Se), Specificity (Sp), Positive Predictive Value (PPV), and Negative Predictive Value (NPV), as well as their respective confidence intervals using advanced statistical methods. This module is highly suitable for medical statistics, diagnostic test evaluation, and fields where accurate diagnostic performance metrics are crucial.

## Features
- **Sensitivity (Se) and Specificity (Sp) Calculation:** Calculate the true positive rate (sensitivity) and true negative rate (specificity) of tests.
- **Predictive Values:** Compute both PPV and NPV, which indicate the performance of a diagnostic test depending on the disease prevalence.
- **Confidence Intervals:** Generate confidence intervals for Se, Sp, PPV, and NPV using both standard and logit transformed methods to handle the skewness of binomial distributions effectively, especially near boundary values.
- **Logit Transformations:** Employ logit transformations for variance stabilization in PPV and NPV estimations, providing more reliable confidence interval estimations under various conditions including small sample sizes and near-boundary probability estimates.

## Installation
To use the Mercaldo Calculations Module, Python along with the `numpy` and `scipy` packages are required. Install them using pip if not already installed:

```bash
pip install numpy scipy
```

You can then download this module directly from our GitHub repository or clone it:

```bash
git clone https://github.com/yourgithubrepo/mercaldo.git
```

## Usage
Import the `calculate_estimates` function from the module:

```python
from mercaldo.calculations import calculate_estimates
```

Example of how to call the `calculate_estimates` function:

```python
# Example data
x11 = 90  # True Positives
x10 = 10  # False Positives
x01 = 5   # False Negatives
x00 = 95  # True Negatives
prevalence = 0.1
confidence = 0.95

# Calculate metrics
results = calculate_estimates(x11, x10, x01, x00, prevalence, confidence)

# Output results
for key, value in results.items():
    print(f"{key}: {value}")
```

## Detailed Methodology
This module calculates PPV and NPV by considering both the intrinsic test accuracy (Se and Sp) and the disease prevalence. Confidence intervals are derived using exact binomial methods (Clopper-Pearson) for Se and Sp, and adjusted logit transformations for PPV and NPV to manage the skewness and discretization issues associated with standard binomial estimates, particularly useful in scenarios with small sample sizes or extreme probabilities.

## Contributing
Contributions are welcome. Please fork the repository, make changes, and submit pull requests. You can also open an issue for bugs or suggestions.

## License
This software is provided "as is", without warranty of any kind, express or implied. See the LICENSE file for more details.
